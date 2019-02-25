//! The parser reads a scene description more or less line by line and
//! stores the read information by calling API functions starting with
//! *pbrt_*.

// std
use std;
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;
// pbrt
use accelerators::bvh::{BVHAccel, SplitMethod};
use cameras::environment::EnvironmentCamera;
use cameras::orthographic::OrthographicCamera;
use cameras::perspective::PerspectiveCamera;
use cameras::realistic::RealisticCamera;
use core::camera::Camera;
use core::film::Film;
use core::filter::Filter;
use core::geometry::{vec3_coordinate_system, vec3_cross_vec3};
use core::geometry::{Bounds2f, Bounds2i, Normal3f, Point2f, Point2i, Point3f, Vector3f};
use core::integrator::SamplerIntegrator;
use core::light::Light;
use core::material::Material;
use core::medium::get_medium_scattering_properties;
use core::medium::{Medium, MediumInterface};
use core::mipmap::ImageWrap;
use core::paramset::{ParamSet, TextureParams};
use core::pbrt::{clamp_t, lerp};
use core::pbrt::{Float, Spectrum};
use core::primitive::{GeometricPrimitive, Primitive, TransformedPrimitive};
use core::reflection::FourierBSDFTable;
use core::sampler::Sampler;
use core::scene::Scene;
use core::shape::Shape;
use core::texture::{PlanarMapping2D, Texture, TextureMapping2D, UVMapping2D};
use core::transform::{AnimatedTransform, Matrix4x4, Transform};
use filters::boxfilter::BoxFilter;
use filters::gaussian::GaussianFilter;
use filters::mitchell::MitchellNetravali;
use filters::triangle::TriangleFilter;
use integrators::ao::AOIntegrator;
use integrators::bdpt::render_bdpt;
use integrators::bdpt::BDPTIntegrator;
use integrators::directlighting::{DirectLightingIntegrator, LightStrategy};
use integrators::mlt::render_mlt;
use integrators::mlt::MLTIntegrator;
use integrators::path::PathIntegrator;
use integrators::render;
use integrators::sppm::render_sppm;
use integrators::sppm::SPPMIntegrator;
use lights::diffuse::DiffuseAreaLight;
use lights::distant::DistantLight;
use lights::infinite::InfiniteAreaLight;
use lights::point::PointLight;
use lights::spot::SpotLight;
use materials::disney::DisneyMaterial;
use materials::fourier::FourierMaterial;
use materials::glass::GlassMaterial;
use materials::hair::HairMaterial;
use materials::matte::MatteMaterial;
use materials::metal::MetalMaterial;
use materials::mirror::MirrorMaterial;
use materials::mixmat::MixMaterial;
use materials::plastic::PlasticMaterial;
use materials::substrate::SubstrateMaterial;
use materials::subsurface::SubsurfaceMaterial;
use materials::uber::UberMaterial;
use media::homogeneous::HomogeneousMedium;
use samplers::halton::HaltonSampler;
use samplers::random::RandomSampler;
use samplers::sobol::SobolSampler;
use samplers::zerotwosequence::ZeroTwoSequenceSampler;
use shapes::curve::create_curve_shape;
use shapes::cylinder::Cylinder;
use shapes::disk::Disk;
use shapes::loopsubdiv::loop_subdivide;
use shapes::nurbs::nurbs_evaluate_surface;
use shapes::nurbs::Homogeneous3;
use shapes::plymesh::create_ply_mesh;
use shapes::sphere::Sphere;
use shapes::triangle::{Triangle, TriangleMesh};
use textures::checkerboard::Checkerboard2DTexture;
use textures::constant::ConstantTexture;
use textures::imagemap::ImageTexture;
use textures::imagemap::{convert_to_float, convert_to_spectrum};
use textures::scale::ScaleTexture;

// see api.cpp

pub struct BsdfState {
    pub loaded_bsdfs: HashMap<String, Arc<FourierBSDFTable>>,
}

impl Default for BsdfState {
    fn default() -> Self {
        BsdfState {
            loaded_bsdfs: HashMap::new(),
        }
    }
}

pub struct ApiState {
    number_of_threads: u8,
    pub search_directory: Option<Box<PathBuf>>,
    cur_transform: TransformSet,
    active_transform_bits: u8,
    named_coordinate_systems: HashMap<&'static str, TransformSet>,
    render_options: RenderOptions,
    graphics_state: GraphicsState,
    pushed_graphics_states: Vec<GraphicsState>,
    pushed_transforms: Vec<TransformSet>,
    pushed_active_transform_bits: Vec<u8>,
    param_set: ParamSet,
}

impl Default for ApiState {
    fn default() -> Self {
        ApiState {
            number_of_threads: 0_u8,
            search_directory: None,
            cur_transform: TransformSet {
                t: [Transform {
                    m: Matrix4x4 {
                        m: [
                            [1.0, 0.0, 0.0, 0.0],
                            [0.0, 1.0, 0.0, 0.0],
                            [0.0, 0.0, 1.0, 0.0],
                            [0.0, 0.0, 0.0, 1.0],
                        ],
                    },
                    m_inv: Matrix4x4 {
                        m: [
                            [1.0, 0.0, 0.0, 0.0],
                            [0.0, 1.0, 0.0, 0.0],
                            [0.0, 0.0, 1.0, 0.0],
                            [0.0, 0.0, 0.0, 1.0],
                        ],
                    },
                }; 2],
            },
            named_coordinate_systems: HashMap::new(),
            active_transform_bits: 3_u8, // 0x11 for MaxTransforms = 2
            render_options: RenderOptions::default(),
            graphics_state: GraphicsState::new(),
            pushed_graphics_states: Vec::new(),
            pushed_transforms: Vec::new(),
            pushed_active_transform_bits: Vec::new(),
            param_set: ParamSet::default(),
        }
    }
}

#[derive(Debug, Default, Copy, Clone)]
pub struct TransformSet {
    pub t: [Transform; 2],
}

impl TransformSet {
    pub fn is_animated(&self) -> bool {
        // for (int i = 0; i < MaxTransforms - 1; ++i)
        //     if (t[i] != t[i + 1]) return true;
        // return false;

        // we have only 2 transforms
        if self.t[0] != self.t[1] {
            true
        } else {
            false
        }
    }
}

pub struct RenderOptions {
    pub transform_start_time: Float,
    pub transform_end_time: Float,
    pub filter_name: String, // "box"
    pub filter_params: ParamSet,
    pub film_name: String, // "image"
    pub film_params: ParamSet,
    pub sampler_name: String, // "halton";
    pub sampler_params: ParamSet,
    pub accelerator_name: String, // "bvh";
    pub accelerator_params: ParamSet,
    pub integrator_name: String, // "path";
    pub integrator_params: ParamSet,
    pub camera_name: String, // "perspective";
    pub camera_params: ParamSet,
    pub camera_to_world: TransformSet,
    pub named_media: HashMap<String, Arc<Medium + Sync + Send>>,
    pub lights: Vec<Arc<Light + Sync + Send>>,
    pub primitives: Vec<Arc<Primitive + Sync + Send>>,
    pub instances: HashMap<String, Vec<Arc<Primitive + Sync + Send>>>,
    pub current_instance: String,
    pub have_scattering_media: bool, // false
}

impl Default for RenderOptions {
    fn default() -> RenderOptions {
        RenderOptions {
            transform_start_time: 0.0 as Float,
            transform_end_time: 1.0 as Float,
            filter_name: String::from("box"),
            filter_params: ParamSet::default(),
            film_name: String::from("image"),
            film_params: ParamSet::default(),
            sampler_name: String::from("halton"),
            sampler_params: ParamSet::default(),
            accelerator_name: String::from("bvh"),
            accelerator_params: ParamSet::default(),
            integrator_name: String::from("path"),
            integrator_params: ParamSet::default(),
            camera_name: String::from("perspective"),
            camera_params: ParamSet::default(),
            camera_to_world: TransformSet {
                t: [Transform {
                    m: Matrix4x4 {
                        m: [
                            [1.0, 0.0, 0.0, 0.0],
                            [0.0, 1.0, 0.0, 0.0],
                            [0.0, 0.0, 1.0, 0.0],
                            [0.0, 0.0, 0.0, 1.0],
                        ],
                    },
                    m_inv: Matrix4x4 {
                        m: [
                            [1.0, 0.0, 0.0, 0.0],
                            [0.0, 1.0, 0.0, 0.0],
                            [0.0, 0.0, 1.0, 0.0],
                            [0.0, 0.0, 0.0, 1.0],
                        ],
                    },
                }; 2],
            },
            named_media: HashMap::new(),
            lights: Vec::new(),
            primitives: Vec::new(),
            instances: HashMap::new(),
            current_instance: String::from(""),
            have_scattering_media: false,
        }
    }
}

#[derive(Default)]
pub struct GraphicsState {
    pub current_inside_medium: String,
    pub current_outside_medium: String,
    pub float_textures: HashMap<String, Arc<Texture<Float> + Send + Sync>>,
    pub spectrum_textures: HashMap<String, Arc<Texture<Spectrum> + Send + Sync>>,
    pub material_params: ParamSet,
    pub material: String,
    pub named_materials: HashMap<String, Option<Arc<Material + Send + Sync>>>,
    pub current_material: String,
    pub area_light_params: ParamSet,
    pub area_light: String,
    pub reverse_orientation: bool,
}

impl GraphicsState {
    pub fn new() -> Self {
        let float_textures: HashMap<String, Arc<Texture<Float> + Send + Sync>> = HashMap::new();
        let spectrum_textures: HashMap<String, Arc<Texture<Spectrum> + Send + Sync>> =
            HashMap::new();
        let mut tp: TextureParams = TextureParams::new(
            ParamSet::default(),
            ParamSet::default(),
            float_textures.clone(),
            spectrum_textures.clone(),
        );
        let mtl: Arc<Material + Send + Sync> = MatteMaterial::create(&mut tp);
        let mut named_materials: HashMap<String, Option<Arc<Material + Send + Sync>>> =
            HashMap::new();
        named_materials.insert(String::from("matte"), Some(mtl));
        let current_material: String = String::from("matte");
        GraphicsState {
            current_inside_medium: String::from(""),
            current_outside_medium: String::from(""),
            float_textures: float_textures.clone(),
            spectrum_textures: spectrum_textures.clone(),
            material_params: ParamSet::default(),
            material: String::from(""),
            named_materials: named_materials,
            current_material: current_material,
            area_light_params: ParamSet::default(),
            area_light: String::from(""),
            reverse_orientation: false,
        }
    }
}

fn create_material(
    api_state: &ApiState,
    bsdf_state: &mut BsdfState,
) -> Option<Arc<Material + Send + Sync>> {
    // CreateMaterial
    let mut material_params = ParamSet::default();
    material_params.copy_from(&api_state.graphics_state.material_params);
    let mut mp: TextureParams = TextureParams {
        float_textures: api_state.graphics_state.float_textures.clone(),
        spectrum_textures: api_state.graphics_state.spectrum_textures.clone(),
        geom_params: ParamSet::default(),
        material_params: material_params,
    };
    if api_state.graphics_state.current_material != String::new() {
        match api_state
            .graphics_state
            .named_materials
            .get(api_state.graphics_state.current_material.as_str())
        {
            Some(named_material) => {
                return named_material.clone();
            }
            None => {
                println!(
                    "WARNING: Named material \"{}\" not defined. Using \"matte\".",
                    api_state.graphics_state.current_material
                );
            }
        }
    } else {
        // MakeMaterial
        if api_state.graphics_state.material == "" || api_state.graphics_state.material == "none" {
            return None;
        } else if api_state.graphics_state.material == "matte" {
            return Some(MatteMaterial::create(&mut mp));
        } else if api_state.graphics_state.material == "plastic" {
            return Some(PlasticMaterial::create(&mut mp));
        } else if api_state.graphics_state.material == "translucent" {
            println!("TODO: CreateTranslucentMaterial");
        } else if api_state.graphics_state.material == "glass" {
            return Some(GlassMaterial::create(&mut mp));
        } else if api_state.graphics_state.material == "mirror" {
            return Some(MirrorMaterial::create(&mut mp));
        } else if api_state.graphics_state.material == "hair" {
            return Some(HairMaterial::create(&mut mp));
        } else if api_state.graphics_state.material == "mix" {
            let m1: String = mp.find_string("namedmaterial1", String::from(""));
            let m2: String = mp.find_string("namedmaterial2", String::from(""));
            let mat1 = match api_state.graphics_state.named_materials.get(&m1) {
                Some(named_material) => named_material,
                None => {
                    panic!("Material \"{}\" unknown.", m1);
                }
            };
            let mat2 = match api_state.graphics_state.named_materials.get(&m2) {
                Some(named_material) => named_material,
                None => {
                    panic!("Material \"{}\" unknown.", m2);
                }
            };
            let scale: Arc<Texture<Spectrum> + Send + Sync> =
                mp.get_spectrum_texture("amount", Spectrum::new(0.5));
            if let Some(m1) = mat1 {
                if let Some(m2) = mat2 {
                    let mix = Arc::new(MixMaterial::new(m1.clone(), m2.clone(), scale));
                    return Some(mix);
                }
            }
            return None;
        } else if api_state.graphics_state.material == "metal" {
            return Some(MetalMaterial::create(&mut mp));
        } else if api_state.graphics_state.material == "substrate" {
            return Some(SubstrateMaterial::create(&mut mp));
        } else if api_state.graphics_state.material == "uber" {
            return Some(UberMaterial::create(&mut mp));
        } else if api_state.graphics_state.material == "subsurface" {
            return Some(SubsurfaceMaterial::create(&mut mp));
        } else if api_state.graphics_state.material == "kdsubsurface" {
            println!("TODO: CreateKdsubsurfaceMaterial");
        } else if api_state.graphics_state.material == "fourier" {
            return Some(FourierMaterial::create(&mut mp, bsdf_state));
        } else if api_state.graphics_state.material == "disney" {
            return Some(DisneyMaterial::create(&mut mp));
        } else {
            panic!(
                "Material \"{}\" unknown.",
                api_state.graphics_state.material
            );
        }
    }
    let kd = Arc::new(ConstantTexture::new(Spectrum::new(0.5)));
    let sigma = Arc::new(ConstantTexture::new(0.0 as Float));
    Some(Arc::new(MatteMaterial::new(kd, sigma, None)))
}

fn create_medium_interface(api_state: &ApiState) -> MediumInterface {
    let mut m: MediumInterface = MediumInterface::default();
    if api_state.graphics_state.current_inside_medium != String::from("") {
        match api_state
            .render_options
            .named_media
            .get(&api_state.graphics_state.current_inside_medium)
        {
            Some(inside_medium_arc) => m.inside = Some(inside_medium_arc.clone()),
            None => {
                panic!(
                    "ERROR: Named medium \"{:?}\" undefined.",
                    api_state.graphics_state.current_inside_medium
                );
            }
        }
    }
    if api_state.graphics_state.current_outside_medium != String::from("") {
        match api_state
            .render_options
            .named_media
            .get(&api_state.graphics_state.current_outside_medium)
        {
            Some(outside_medium_arc) => m.outside = Some(outside_medium_arc.clone()),
            None => {
                panic!(
                    "ERROR: Named medium \"{:?}\" undefined.",
                    api_state.graphics_state.current_outside_medium
                );
            }
        }
    }
    m
}

fn make_light(api_state: &mut ApiState, medium_interface: &MediumInterface) {
    // MakeLight (api.cpp:591)
    if api_state.param_set.name == "point" {
        let i: Spectrum = api_state
            .param_set
            .find_one_spectrum("I", Spectrum::new(1.0 as Float));
        let sc: Spectrum = api_state
            .param_set
            .find_one_spectrum("scale", Spectrum::new(1.0 as Float));
        let p: Point3f = api_state
            .param_set
            .find_one_point3f("from", Point3f::default());
        let l2w: Transform = Transform::translate(&Vector3f {
            x: p.x,
            y: p.y,
            z: p.z,
        }) * api_state.cur_transform.t[0];
        let point_light = Arc::new(PointLight::new(&l2w, medium_interface, &(i * sc)));
        api_state.render_options.lights.push(point_light);
    } else if api_state.param_set.name == "spot" {
        // CreateSpotLight
        let i: Spectrum = api_state
            .param_set
            .find_one_spectrum("I", Spectrum::new(1.0 as Float));
        let sc: Spectrum = api_state
            .param_set
            .find_one_spectrum("scale", Spectrum::new(1.0 as Float));
        let coneangle: Float = api_state
            .param_set
            .find_one_float("coneangle", 30.0 as Float);
        let conedelta: Float = api_state
            .param_set
            .find_one_float("conedeltaangle", 5.0 as Float);
        // compute spotlight world to light transformation
        let from: Point3f = api_state.param_set.find_one_point3f(
            "from",
            Point3f {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
        );
        let to: Point3f = api_state.param_set.find_one_point3f(
            "to",
            Point3f {
                x: 0.0,
                y: 0.0,
                z: 1.0,
            },
        );
        let dir: Vector3f = (to - from).normalize();
        let mut du: Vector3f = Vector3f::default();
        let mut dv: Vector3f = Vector3f::default();
        vec3_coordinate_system(&dir, &mut du, &mut dv);
        let dir_to_z: Transform = Transform::new(
            du.x, du.y, du.z, 0.0, dv.x, dv.y, dv.z, 0.0, dir.x, dir.y, dir.z, 0.0, 0.0, 0.0, 0.0,
            1.0,
        );
        let light2world: Transform = api_state.cur_transform.t[0]
            * Transform::translate(&Vector3f {
                x: from.x,
                y: from.y,
                z: from.z,
            })
            * Transform::inverse(&dir_to_z);
        // return std::make_shared<SpotLight>(light2world, medium, I * sc, coneangle, coneangle - conedelta);
        let spot_light = Arc::new(SpotLight::new(
            &light2world,
            medium_interface,
            &(i * sc),
            coneangle,
            coneangle - conedelta,
        ));
        api_state.render_options.lights.push(spot_light);
    } else if api_state.param_set.name == "goniometric" {
        println!("TODO: CreateGoniometricLight");
    } else if api_state.param_set.name == "projection" {
        println!("TODO: CreateProjectionLight");
    } else if api_state.param_set.name == "distant" {
        // CreateDistantLight
        let l: Spectrum = api_state
            .param_set
            .find_one_spectrum("L", Spectrum::new(1.0 as Float));
        let sc: Spectrum = api_state
            .param_set
            .find_one_spectrum("scale", Spectrum::new(1.0 as Float));
        let from: Point3f = api_state.param_set.find_one_point3f(
            "from",
            Point3f {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
        );
        let to: Point3f = api_state.param_set.find_one_point3f(
            "to",
            Point3f {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
        );
        let dir: Vector3f = from - to;
        // return std::make_shared<DistantLight>(light2world, L * sc, dir);
        let distant_light = Arc::new(DistantLight::new(
            &api_state.cur_transform.t[0],
            &(l * sc),
            &dir,
        ));
        api_state.render_options.lights.push(distant_light);
    } else if api_state.param_set.name == "infinite" || api_state.param_set.name == "exinfinite" {
        let l: Spectrum = api_state
            .param_set
            .find_one_spectrum("L", Spectrum::new(1.0 as Float));
        let sc: Spectrum = api_state
            .param_set
            .find_one_spectrum("scale", Spectrum::new(1.0 as Float));
        let mut texmap: String = api_state
            .param_set
            .find_one_filename("mapname", String::from(""));
        if texmap != String::from("") {
            if let Some(ref search_directory) = api_state.search_directory {
                // texmap = AbsolutePath(ResolveFilename(texmap));
                let mut path_buf: PathBuf = PathBuf::from("/");
                path_buf.push(search_directory.as_ref());
                path_buf.push(texmap);
                texmap = String::from(path_buf.to_str().unwrap());
            }
        }
        let n_samples: i32 = api_state.param_set.find_one_int("nsamples", 1 as i32);
        // TODO: if (PbrtOptions.quickRender) nSamples = std::max(1, nSamples / 4);

        // return std::make_shared<InfiniteAreaLight>(light2world, L * sc, nSamples, texmap);
        let infinte_light = Arc::new(InfiniteAreaLight::new(
            &api_state.cur_transform.t[0],
            &(l * sc),
            n_samples,
            texmap,
        ));
        api_state.render_options.lights.push(infinte_light);
    } else {
        panic!("MakeLight: unknown name {}", api_state.param_set.name);
    }
}

fn make_medium(api_state: &mut ApiState) {
    let medium_type: String = api_state.param_set.find_one_string("type", String::new());
    if medium_type == "" {
        panic!("ERROR: No parameter string \"type\" found in MakeNamedMedium");
    }
    // MakeMedium (api.cpp:685)
    let sig_a_rgb: [Float; 3] = [0.0011, 0.0024, 0.014];
    let sig_s_rgb: [Float; 3] = [2.55, 3.21, 3.77];
    let mut sig_a: Spectrum = Spectrum::from_rgb(&sig_a_rgb);
    let mut sig_s: Spectrum = Spectrum::from_rgb(&sig_s_rgb);
    let preset: String = api_state.param_set.find_one_string("preset", String::new());
    let found: bool = get_medium_scattering_properties(&preset, &mut sig_a, &mut sig_s);
    if preset != String::from("") && !found {
        println!(
            "WARNING: Material preset \"{:?}\" not found.  Using defaults.",
            preset
        );
    }
    let scale: Float = api_state.param_set.find_one_float("scale", 1.0 as Float);
    let g: Float = api_state.param_set.find_one_float("g", 0.0 as Float);
    sig_a = api_state.param_set.find_one_spectrum("sigma_a", sig_a) * scale;
    sig_s = api_state.param_set.find_one_spectrum("sigma_s", sig_s) * scale;
    let some_medium: Option<Arc<Medium + Sync + Send>>;
    if medium_type == "homogeneous" {
        some_medium = Some(Arc::new(HomogeneousMedium::new(&sig_a, &sig_s, g)));
    } else if medium_type == "heterogeneous" {
        panic!("TODO: make_medium(\"heterogeneous\")");
    } else {
        panic!("MakeMedium: unknown name {}", medium_type);
    }
    if let Some(medium) = some_medium {
        api_state
            .render_options
            .named_media
            .insert(api_state.param_set.name.clone(), medium);
    }
}

fn make_texture(api_state: &mut ApiState) {
    // pbrtTexture (api.cpp:1049)
    let mut geom_params: ParamSet = ParamSet::default();
    let mut material_params: ParamSet = ParamSet::default();
    geom_params.copy_from(&api_state.param_set);
    material_params.copy_from(&api_state.param_set);
    let mut tp: TextureParams = TextureParams {
        float_textures: api_state.graphics_state.float_textures.clone(),
        spectrum_textures: api_state.graphics_state.spectrum_textures.clone(),
        geom_params: geom_params,
        material_params: material_params,
    };
    if api_state.param_set.tex_type == "float" {
        match api_state
            .graphics_state
            .spectrum_textures
            .get(api_state.param_set.name.as_str())
        {
            Some(_spectrum_texture) => {
                println!("Texture \"{}\" being redefined", api_state.param_set.name);
            }
            None => {}
        }
        // TODO: WARN_IF_ANIMATED_TRANSFORM("Texture");
        // MakeFloatTexture(texname, curTransform[0], tp);
        if api_state.param_set.tex_name == "constant" {
            println!("TODO: CreateConstantFloatTexture");
        } else if api_state.param_set.tex_name == "scale" {
            let ft = Arc::new(ScaleTexture::<Float>::new(
                tp.get_float_texture("tex1", 1.0 as Float),
                tp.get_float_texture("tex2", 1.0 as Float),
            ));
            api_state
                .graphics_state
                .float_textures
                .insert(api_state.param_set.name.clone(), ft);
        } else if api_state.param_set.tex_name == "mix" {
            println!("TODO: CreateMixFloatTexture");
        } else if api_state.param_set.tex_name == "bilerp" {
            println!("TODO: CreateBilerpFloatTexture");
        } else if api_state.param_set.tex_name == "imagemap" {
            // CreateImageFloatTexture
            let mut map: Option<Box<TextureMapping2D + Send + Sync>> = None;
            let mapping: String = tp.find_string("mapping", String::from("uv"));
            if mapping == "uv" {
                let su: Float = tp.find_float("uscale", 1.0);
                let sv: Float = tp.find_float("vscale", 1.0);
                let du: Float = tp.find_float("udelta", 0.0);
                let dv: Float = tp.find_float("vdelta", 0.0);
                map = Some(Box::new(UVMapping2D {
                    su: su,
                    sv: sv,
                    du: du,
                    dv: dv,
                }));
            } else if mapping == "spherical" {
                println!("TODO: SphericalMapping2D");
            } else if mapping == "cylindrical" {
                println!("TODO: CylindricalMapping2D");
            } else if mapping == "planar" {
                map = Some(Box::new(PlanarMapping2D {
                    vs: tp.find_vector3f(
                        "v1",
                        Vector3f {
                            x: 1.0,
                            y: 0.0,
                            z: 0.0,
                        },
                    ),
                    vt: tp.find_vector3f(
                        "v2",
                        Vector3f {
                            x: 0.0,
                            y: 1.0,
                            z: 0.0,
                        },
                    ),
                    ds: tp.find_float("udelta", 0.0),
                    dt: tp.find_float("vdelta", 0.0),
                }));
            } else {
                panic!("2D texture mapping \"{}\" unknown", mapping);
            }
            // initialize _ImageTexture_ parameters
            let max_aniso: Float = tp.find_float("maxanisotropy", 8.0);
            let do_trilinear: bool = tp.find_bool("trilinear", false);
            let wrap: String = tp.find_string("wrap", String::from("repeat"));
            let mut wrap_mode: ImageWrap = ImageWrap::Repeat;
            if wrap == "black" {
                wrap_mode = ImageWrap::Black;
            } else if wrap == "clamp" {
                wrap_mode = ImageWrap::Clamp;
            }
            let scale: Float = tp.find_float("scale", 1.0);
            let mut filename: String = tp.find_filename("filename", String::new());
            if let Some(ref search_directory) = api_state.search_directory {
                // filename = AbsolutePath(ResolveFilename(filename));
                let mut path_buf: PathBuf = PathBuf::from("/");
                path_buf.push(search_directory.as_ref());
                path_buf.push(filename);
                filename = String::from(path_buf.to_str().unwrap());
            }
            // TODO: default depends on:
            // HasExtension(filename,
            // ".tga") ||
            // HasExtension(filename,
            // ".png"));
            let gamma: bool = tp.find_bool("gamma", true);

            if let Some(mapping) = map {
                let ft = Arc::new(ImageTexture::new(
                    mapping,
                    filename,
                    do_trilinear,
                    max_aniso,
                    wrap_mode,
                    scale,
                    gamma,
                    convert_to_float,
                ));
                api_state
                    .graphics_state
                    .float_textures
                    .insert(api_state.param_set.name.clone(), ft);
            }
        } else if api_state.param_set.tex_name == "uv" {
            println!("TODO: CreateUVFloatTexture");
        } else if api_state.param_set.tex_name == "checkerboard" {
            println!("TODO: CreateCheckerboardFloatTexture");
        } else if api_state.param_set.tex_name == "dots" {
            println!("TODO: CreateDotsFloatTexture");
        } else if api_state.param_set.tex_name == "fbm" {
            println!("TODO: CreateFBmFloatTexture");
        } else if api_state.param_set.tex_name == "wrinkled" {
            println!("TODO: CreateWrinkledFloatTexture");
        } else if api_state.param_set.tex_name == "marble" {
            println!("TODO: CreateMarbleFloatTexture");
        } else if api_state.param_set.tex_name == "windy" {
            println!("TODO: CreateWindyFloatTexture");
        } else if api_state.param_set.tex_name == "ptex" {
            println!("TODO: CreatePtexFloatTexture");
        } else {
            println!(
                "Float texture \"{}\" unknown.",
                api_state.param_set.tex_name
            );
        }
    } else if api_state.param_set.tex_type == "color" || api_state.param_set.tex_type == "spectrum"
    {
        match api_state
            .graphics_state
            .spectrum_textures
            .get(api_state.param_set.name.as_str())
        {
            Some(_spectrum_texture) => {
                println!("Texture \"{}\" being redefined", api_state.param_set.name);
            }
            None => {}
        }
        // TODO: WARN_IF_ANIMATED_TRANSFORM("Texture");
        // MakeSpectrumTexture(texname, curTransform[0], tp);
        if api_state.param_set.tex_name == "constant" {
            let ct = Arc::new(ConstantTexture::new(
                tp.find_spectrum("value", Spectrum::new(1.0)),
            ));
            api_state
                .graphics_state
                .spectrum_textures
                .insert(api_state.param_set.name.clone(), ct);
        } else if api_state.param_set.tex_name == "scale" {
            let tex1: Arc<Texture<Spectrum> + Send + Sync> =
                tp.get_spectrum_texture("tex1", Spectrum::new(1.0));
            let tex2: Arc<Texture<Spectrum> + Send + Sync> =
                tp.get_spectrum_texture("tex2", Spectrum::new(0.0));
            let st = Arc::new(ScaleTexture::<Spectrum>::new(tex1, tex2));
            api_state
                .graphics_state
                .spectrum_textures
                .insert(api_state.param_set.name.clone(), st);
        } else if api_state.param_set.tex_name == "mix" {
            println!("TODO: CreateMixSpectrumTexture");
        } else if api_state.param_set.tex_name == "bilerp" {
            println!("TODO: CreateBilerpSpectrumTexture");
        } else if api_state.param_set.tex_name == "imagemap" {
            // CreateImageSpectrumTexture
            let mut map: Option<Box<TextureMapping2D + Send + Sync>> = None;
            let mapping: String = tp.find_string("mapping", String::from("uv"));
            if mapping == "uv" {
                let su: Float = tp.find_float("uscale", 1.0);
                let sv: Float = tp.find_float("vscale", 1.0);
                let du: Float = tp.find_float("udelta", 0.0);
                let dv: Float = tp.find_float("vdelta", 0.0);
                map = Some(Box::new(UVMapping2D {
                    su: su,
                    sv: sv,
                    du: du,
                    dv: dv,
                }));
            } else if mapping == "spherical" {
                println!("TODO: SphericalMapping2D");
            } else if mapping == "cylindrical" {
                println!("TODO: CylindricalMapping2D");
            } else if mapping == "planar" {
                map = Some(Box::new(PlanarMapping2D {
                    vs: tp.find_vector3f(
                        "v1",
                        Vector3f {
                            x: 1.0,
                            y: 0.0,
                            z: 0.0,
                        },
                    ),
                    vt: tp.find_vector3f(
                        "v2",
                        Vector3f {
                            x: 0.0,
                            y: 1.0,
                            z: 0.0,
                        },
                    ),
                    ds: tp.find_float("udelta", 0.0),
                    dt: tp.find_float("vdelta", 0.0),
                }));
            } else {
                panic!("2D texture mapping \"{}\" unknown", mapping);
            }
            // initialize _ImageTexture_ parameters
            let max_aniso: Float = tp.find_float("maxanisotropy", 8.0);
            let do_trilinear: bool = tp.find_bool("trilinear", false);
            let wrap: String = tp.find_string("wrap", String::from("repeat"));
            let mut wrap_mode: ImageWrap = ImageWrap::Repeat;
            if wrap == "black" {
                wrap_mode = ImageWrap::Black;
            } else if wrap == "clamp" {
                wrap_mode = ImageWrap::Clamp;
            }
            let scale: Float = tp.find_float("scale", 1.0);
            let mut filename: String = tp.find_filename("filename", String::new());
            if let Some(ref search_directory) = api_state.search_directory {
                // filename = AbsolutePath(ResolveFilename(filename));
                let mut path_buf: PathBuf = PathBuf::from("/");
                path_buf.push(search_directory.as_ref());
                path_buf.push(filename);
                filename = String::from(path_buf.to_str().unwrap());
            }
            // TODO: default depends on:
            // HasExtension(filename,
            // ".tga") ||
            // HasExtension(filename,
            // ".png"));
            let gamma: bool = tp.find_bool("gamma", true);

            if let Some(mapping) = map {
                let st = Arc::new(ImageTexture::new(
                    mapping,
                    filename,
                    do_trilinear,
                    max_aniso,
                    wrap_mode,
                    scale,
                    gamma,
                    convert_to_spectrum,
                ));
                api_state
                    .graphics_state
                    .spectrum_textures
                    .insert(api_state.param_set.name.clone(), st);
            }
        } else if api_state.param_set.tex_name == "uv" {
            println!("TODO: CreateUVSpectrumTexture");
        } else if api_state.param_set.tex_name == "checkerboard" {
            // CreateCheckerboardSpectrumTexture
            let dim: i32 = tp.find_int("dimension", 2);
            if dim != 2 && dim != 3 {
                panic!("{} dimensional checkerboard texture not supported", dim);
            }
            let tex1: Arc<Texture<Spectrum> + Send + Sync> =
                tp.get_spectrum_texture("tex1", Spectrum::new(1.0));
            let tex2: Arc<Texture<Spectrum> + Send + Sync> =
                tp.get_spectrum_texture("tex2", Spectrum::new(0.0));
            if dim == 2 {
                let mut map: Option<Box<TextureMapping2D + Send + Sync>> = None;
                let mapping: String = tp.find_string("mapping", String::from("uv"));
                if mapping == "uv" {
                    let su: Float = tp.find_float("uscale", 1.0);
                    let sv: Float = tp.find_float("vscale", 1.0);
                    let du: Float = tp.find_float("udelta", 0.0);
                    let dv: Float = tp.find_float("vdelta", 0.0);
                    map = Some(Box::new(UVMapping2D {
                        su: su,
                        sv: sv,
                        du: du,
                        dv: dv,
                    }));
                } else if mapping == "spherical" {
                    println!("TODO: SphericalMapping2D");
                } else if mapping == "cylindrical" {
                    println!("TODO: CylindricalMapping2D");
                } else if mapping == "planar" {
                    map = Some(Box::new(PlanarMapping2D {
                        vs: tp.find_vector3f(
                            "v1",
                            Vector3f {
                                x: 1.0,
                                y: 0.0,
                                z: 0.0,
                            },
                        ),
                        vt: tp.find_vector3f(
                            "v2",
                            Vector3f {
                                x: 0.0,
                                y: 1.0,
                                z: 0.0,
                            },
                        ),
                        ds: tp.find_float("udelta", 0.0),
                        dt: tp.find_float("vdelta", 0.0),
                    }));
                } else {
                    panic!("2D texture mapping \"{}\" unknown", mapping);
                }
                // TODO: aamode
                if let Some(mapping) = map {
                    let st = Arc::new(Checkerboard2DTexture::new(mapping, tex1, tex2));
                    api_state
                        .graphics_state
                        .spectrum_textures
                        .insert(api_state.param_set.name.clone(), st);
                }
            } else {
                // dim == 3
                println!("TODO: TextureMapping3D");
            }
        } else if api_state.param_set.tex_name == "dots" {
            println!("TODO: CreateDotsSpectrumTexture");
        } else if api_state.param_set.tex_name == "fbm" {
            println!("TODO: CreateFBmSpectrumTexture");
        } else if api_state.param_set.tex_name == "wrinkled" {
            println!("TODO: CreateWrinkledSpectrumTexture");
        } else if api_state.param_set.tex_name == "marble" {
            println!("TODO: CreateMarbleSpectrumTexture");
        } else if api_state.param_set.tex_name == "windy" {
            println!("TODO: CreateWindySpectrumTexture");
        } else {
            println!(
                "Spectrum texture \"{}\" unknown.",
                api_state.param_set.tex_name
            );
        }
    } else {
        panic!("Texture type \"{}\" unknown.", api_state.param_set.tex_type);
    }
    // MakeFloatTexture(texname, curTransform[0], tp);
    // or
    // MakeSpectrumTexture(texname, curTransform[0], tp);
}

fn get_shapes_and_materials(
    api_state: &ApiState,
    bsdf_state: &mut BsdfState,
) -> (
    Vec<Arc<Shape + Send + Sync>>,
    Vec<Option<Arc<Material + Send + Sync>>>,
) {
    let mut shapes: Vec<Arc<Shape + Send + Sync>> = Vec::new();
    let mut materials: Vec<Option<Arc<Material + Send + Sync>>> = Vec::new();
    // pbrtShape (api.cpp:1153)
    // TODO: if (!curTransform.IsAnimated()) { ... }
    // TODO: transformCache.Lookup(curTransform[0], &ObjToWorld, &WorldToObj);
    let mut obj_to_world: Transform = Transform {
        m: api_state.cur_transform.t[0].m,
        m_inv: api_state.cur_transform.t[0].m_inv,
    };
    let mut world_to_obj: Transform = Transform {
        m: api_state.cur_transform.t[0].m_inv,
        m_inv: api_state.cur_transform.t[0].m,
    };
    if api_state.cur_transform.is_animated() {
        if api_state.graphics_state.area_light != String::from("") {
            println!("WARNING: Ignoring currently set area light when creating animated shape",);
        }
        // WORK
        // set both transforms to identity
        obj_to_world = Transform::default();
        world_to_obj = Transform::default();
    }
    // MakeShapes (api.cpp:296)
    if api_state.param_set.name == "sphere" {
        // CreateSphereShape
        let radius: Float = api_state.param_set.find_one_float("radius", 1.0 as Float);
        let z_min: Float = api_state.param_set.find_one_float("zmin", -radius);
        let z_max: Float = api_state.param_set.find_one_float("zmax", radius);
        let phi_max: Float = api_state.param_set.find_one_float("phimax", 360.0 as Float);
        let sphere = Arc::new(Sphere::new(
            obj_to_world,
            world_to_obj,
            false,
            false,
            radius,
            z_min,
            z_max,
            phi_max,
        ));
        let mtl: Option<Arc<Material + Send + Sync>> = create_material(&api_state, bsdf_state);
        shapes.push(sphere.clone());
        materials.push(mtl);
    } else if api_state.param_set.name == "cylinder" {
        let radius: Float = api_state.param_set.find_one_float("radius", 1.0);
        let z_min: Float = api_state.param_set.find_one_float("zmin", -radius);
        let z_max: Float = api_state.param_set.find_one_float("zmax", radius);
        let phi_max: Float = api_state.param_set.find_one_float("phimax", 360.0 as Float);
        let cylinder = Arc::new(Cylinder::new(
            obj_to_world,
            world_to_obj,
            false,
            radius,
            z_min,
            z_max,
            phi_max,
        ));
        let mtl: Option<Arc<Material + Send + Sync>> = create_material(&api_state, bsdf_state);
        shapes.push(cylinder.clone());
        materials.push(mtl.clone());
    } else if api_state.param_set.name == "disk" {
        let height: Float = api_state.param_set.find_one_float("height", 0.0);
        let radius: Float = api_state.param_set.find_one_float("radius", 1.0);
        let inner_radius: Float = api_state.param_set.find_one_float("innerradius", 0.0);
        let phi_max: Float = api_state.param_set.find_one_float("phimax", 360.0);
        let disk = Arc::new(Disk::new(
            obj_to_world,
            world_to_obj,
            false,
            false,
            height,
            radius,
            inner_radius,
            phi_max,
        ));
        let mtl: Option<Arc<Material + Send + Sync>> = create_material(&api_state, bsdf_state);
        shapes.push(disk.clone());
        materials.push(mtl.clone());
    } else if api_state.param_set.name == "cone" {
        println!("TODO: CreateConeShape");
    } else if api_state.param_set.name == "paraboloid" {
        println!("TODO: CreateParaboloidShape");
    } else if api_state.param_set.name == "hyperboloid" {
        println!("TODO: CreateHyperboloidShape");
    } else if api_state.param_set.name == "curve" {
        let mtl: Option<Arc<Material + Send + Sync>> = create_material(&api_state, bsdf_state);
        let curve_shapes: Vec<Arc<Shape + Send + Sync>> = create_curve_shape(
            &obj_to_world,
            &world_to_obj,
            false, // reverse_orientation
            &api_state.param_set,
        );
        for shape in curve_shapes {
            shapes.push(shape.clone());
            materials.push(mtl.clone());
        }
    } else if api_state.param_set.name == "trianglemesh" {
        let vi = api_state.param_set.find_int("indices");
        let p = api_state.param_set.find_point3f("P");
        // try "uv" with Point2f
        let mut uvs = api_state.param_set.find_point2f("uv");
        if uvs.is_empty() {
            // try "st" with Point2f
            uvs = api_state.param_set.find_point2f("st");
        }
        if uvs.is_empty() {
            // try "uv" with float
            let mut fuv = api_state.param_set.find_float("uv");
            if fuv.is_empty() {
                // try "st" with float
                fuv = api_state.param_set.find_float("st");
            }
            if !fuv.is_empty() {
                // found some float UVs
                for i in 0..(fuv.len() / 2) {
                    uvs.push(Point2f {
                        x: fuv[2 * i],
                        y: fuv[2 * i + 1],
                    });
                }
            }
        }
        if !uvs.is_empty() {
            // TODO: if (nuvi < npi) {...} else if (nuvi > npi) ...
            assert!(uvs.len() == p.len());
        }
        assert!(vi.len() > 0_usize);
        assert!(p.len() > 0_usize);
        let s = api_state.param_set.find_vector3f("S");
        let mut s_ws: Vec<Vector3f> = Vec::new();
        if !s.is_empty() {
            assert!(s.len() == p.len());
            // transform tangents to world space
            let n_tangents: usize = s.len();
            for i in 0..n_tangents {
                s_ws.push(obj_to_world.transform_vector(&s[i]));
            }
        }
        let n = api_state.param_set.find_normal3f("N");
        let mut n_ws: Vec<Normal3f> = Vec::new();
        if !n.is_empty() {
            assert!(n.len() == p.len());
            // transform normals to world space
            let n_normals: usize = n.len();
            for i in 0..n_normals {
                n_ws.push(obj_to_world.transform_normal(&n[i]));
            }
        }
        for i in 0..vi.len() {
            if vi[i] as usize >= p.len() {
                panic!(
                    "trianglemesh has out of-bounds vertex index {} ({} \"P\" values were given)",
                    vi[i],
                    p.len()
                );
            }
        }
        // TODO: alpha
        // CreateTriangleMesh
        // transform mesh vertices to world space
        let mut p_ws: Vec<Point3f> = Vec::new();
        let n_vertices: usize = p.len();
        for i in 0..n_vertices {
            p_ws.push(obj_to_world.transform_point(&p[i]));
        }
        // vertex indices are expected as usize, not i32
        let mut vertex_indices: Vec<usize> = Vec::new();
        for i in 0..vi.len() {
            vertex_indices.push(vi[i] as usize);
        }
        let mesh = Arc::new(TriangleMesh::new(
            obj_to_world,
            world_to_obj,
            api_state.graphics_state.reverse_orientation,
            false,        // transform_swaps_handedness
            vi.len() / 3, // n_triangles
            vertex_indices,
            n_vertices,
            p_ws, // in world space
            s_ws, // in world space
            n_ws, // in world space
            uvs,
        ));
        let mtl: Option<Arc<Material + Send + Sync>> = create_material(&api_state, bsdf_state);
        for id in 0..mesh.n_triangles {
            let triangle = Arc::new(Triangle::new(
                mesh.object_to_world,
                mesh.world_to_object,
                mesh.reverse_orientation,
                mesh.clone(),
                id,
            ));
            shapes.push(triangle.clone());
            materials.push(mtl.clone());
        }
    } else if api_state.param_set.name == "plymesh" {
        if let Some(ref search_directory) = api_state.search_directory {
            let mtl: Option<Arc<Material + Send + Sync>> = create_material(&api_state, bsdf_state);
            let ply_shapes: Vec<Arc<Shape + Send + Sync>> = create_ply_mesh(
                &obj_to_world,
                &world_to_obj,
                false, // reverse_orientation
                &api_state.param_set,
                api_state.graphics_state.float_textures.clone(),
                // additional parameters:
                Some(search_directory),
            );
            for shape in ply_shapes {
                shapes.push(shape.clone());
                materials.push(mtl.clone());
            }
        } else {
            panic!("No search directory for plymesh.");
        }
    } else if api_state.param_set.name == "heightfield" {
        println!("TODO: CreateHeightfield");
    } else if api_state.param_set.name == "loopsubdiv" {
        // CreateLoopSubdiv
        let n_levels: i32 = api_state
            .param_set
            .find_one_int("levels", api_state.param_set.find_one_int("nlevels", 3));
        // int nps, nIndices;
        let vertex_indices: Vec<i32> = api_state.param_set.find_int("indices");
        let p = api_state.param_set.find_point3f("P");
        if vertex_indices.is_empty() {
            panic!("Vertex indices \"indices\" not provided for LoopSubdiv shape.");
        }
        if p.is_empty() {
            panic!("Vertex positions \"P\" not provided for LoopSubdiv shape.");
        }
        // don't actually use this for now...
        let _scheme: String = api_state
            .param_set
            .find_one_string("scheme", String::from("loop"));
        let mesh = loop_subdivide(
            &obj_to_world,
            &world_to_obj,
            api_state.graphics_state.reverse_orientation,
            n_levels,
            &vertex_indices,
            &p,
        );
        let mtl: Option<Arc<Material + Send + Sync>> = create_material(&api_state, bsdf_state);
        for id in 0..mesh.n_triangles {
            let triangle = Arc::new(Triangle::new(
                mesh.object_to_world,
                mesh.world_to_object,
                mesh.reverse_orientation,
                mesh.clone(),
                id,
            ));
            shapes.push(triangle.clone());
            materials.push(mtl.clone());
        }
    } else if api_state.param_set.name == "nurbs" {
        // CreateNURBS
        let nu: i32 = api_state.param_set.find_one_int("nu", -1);
        if nu == -1_i32 {
            panic!("Must provide number of control points \"nu\" with NURBS shape.");
        }
        let uorder: i32 = api_state.param_set.find_one_int("uorder", -1);
        if uorder == -1_i32 {
            panic!("Must provide u order \"uorder\" with NURBS shape.");
        }
        let uknots: Vec<Float> = api_state.param_set.find_float("uknots");
        if uknots.is_empty() {
            panic!("Must provide u knot vector \"uknots\" with NURBS shape.");
        }
        if uknots.len() != (nu + uorder) as usize {
            panic!("Number of knots in u knot vector {} doesn't match sum of number of u control points {} and u order {}.",
                   uknots.len(), nu, uorder);
        }
        let u0: Float = api_state
            .param_set
            .find_one_float("u0", uknots[(uorder - 1) as usize]);
        let u1: Float = api_state
            .param_set
            .find_one_float("u1", uknots[nu as usize]);
        let nv: i32 = api_state.param_set.find_one_int("nv", -1);
        if nv == -1_i32 {
            panic!("Must provide number of control points \"nv\" with NURBS shape.");
        }
        let vorder: i32 = api_state.param_set.find_one_int("vorder", -1);
        if vorder == -1_i32 {
            panic!("Must provide u order \"vorder\" with NURBS shape.");
        }
        let vknots: Vec<Float> = api_state.param_set.find_float("vknots");
        if vknots.is_empty() {
            panic!("Must provide u knot vector \"vknots\" with NURBS shape.");
        }
        if vknots.len() != (nv + vorder) as usize {
            panic!("Number of knots in v knot vector {} doesn't match sum of number of v control points {} and v order {}.",
                   vknots.len(), nv, vorder);
        }
        let v0: Float = api_state
            .param_set
            .find_one_float("v0", vknots[(vorder - 1) as usize]);
        let v1: Float = api_state
            .param_set
            .find_one_float("v1", vknots[nv as usize]);
        let mut is_homogeneous: bool = false;
        let p: Vec<Point3f> = api_state.param_set.find_point3f("P");
        let mut pw: Vec<Float> = Vec::new();
        let mut npts: usize = p.len();
        if p.is_empty() {
            pw = api_state.param_set.find_float("Pw");
            if pw.is_empty() {
                panic!("Must provide control points via \"P\" or \"Pw\" parameter to NURBS shape.");
            }
            if pw.len() % 4 != 0 {
                panic!("Number of \"Pw\" control points provided to NURBS shape must be multiple of four");
            }
            npts = pw.len() / 4_usize;
            is_homogeneous = true;
        }
        if npts != (nu * nv) as usize {
            panic!(
                "NURBS shape was expecting {}x{}={} control points, was given {}",
                nu,
                nv,
                nu * nv,
                npts
            );
        }
        // compute NURBS dicing rates
        let diceu: usize = 30;
        let dicev: usize = 30;
        let mut ueval: Vec<Float> = Vec::with_capacity(diceu);
        let mut veval: Vec<Float> = Vec::with_capacity(dicev);
        let mut eval_ps: Vec<Point3f> = Vec::with_capacity(diceu * dicev);
        let mut eval_ns: Vec<Normal3f> = Vec::with_capacity(diceu * dicev);
        for i in 0..diceu {
            ueval.push(lerp(i as Float / (diceu - 1) as Float, u0, u1));
        }
        for i in 0..dicev {
            veval.push(lerp(i as Float / (dicev - 1) as Float, v0, v1));
        }
        // evaluate NURBS over grid of points
        let mut uvs: Vec<Point2f> = Vec::with_capacity(diceu * dicev);
        // turn NURBS into triangles
        let mut hom3: Vec<Homogeneous3> = Vec::with_capacity((nu * nv) as usize);
        if is_homogeneous {
            for i in 0..(nu * nv) as usize {
                hom3.push(Homogeneous3 {
                    x: pw[4 * i],
                    y: pw[4 * i + 1],
                    z: pw[4 * i + 2],
                    w: pw[4 * i + 3],
                });
            }
        } else {
            for i in 0..(nu * nv) as usize {
                hom3.push(Homogeneous3 {
                    x: p[i].x,
                    y: p[i].y,
                    z: p[i].z,
                    w: 1.0 as Float,
                });
            }
        }
        for v in 0..dicev {
            for u in 0..diceu {
                uvs.push(Point2f {
                    x: ueval[u],
                    y: veval[v],
                });
                let mut dpdu: Vector3f = Vector3f::default();
                let mut dpdv: Vector3f = Vector3f::default();
                let pt: Point3f = nurbs_evaluate_surface(
                    uorder,
                    &uknots,
                    nu,
                    ueval[u],
                    vorder,
                    &vknots,
                    nv,
                    veval[v],
                    &hom3,
                    Some(&mut dpdu),
                    Some(&mut dpdv),
                );
                eval_ps.push(Point3f {
                    x: pt.x,
                    y: pt.y,
                    z: pt.z,
                });
                eval_ns.push(Normal3f::from(vec3_cross_vec3(&dpdu, &dpdv).normalize()));
            }
        }
        // generate points-polygons mesh
        let n_tris: usize = 2 * (diceu - 1) * (dicev - 1);
        let mut vertices: Vec<usize> = Vec::with_capacity(3 * n_tris);
        // compute the vertex offset numbers for the triangles
        for v in 0_usize..(dicev - 1) as usize {
            for u in 0_usize..(diceu - 1) as usize {
                vertices.push(v * diceu + u);
                vertices.push(v * diceu + u + 1);
                vertices.push((v + 1) * diceu + u + 1);
                vertices.push(v * diceu + u);
                vertices.push((v + 1) * diceu + u + 1);
                vertices.push((v + 1) * diceu + u);
            }
        }
        // transform mesh vertices to world space
        let mut p_ws: Vec<Point3f> = Vec::new();
        let n_vertices: usize = eval_ps.len();
        for i in 0..n_vertices {
            p_ws.push(obj_to_world.transform_point(&eval_ps[i]));
        }
        // transform normals to world space
        let mut n_ws: Vec<Normal3f> = Vec::new();
        let n_normals: usize = eval_ns.len();
        for i in 0..n_normals {
            n_ws.push(obj_to_world.transform_normal(&eval_ns[i]));
        }
        let mesh = Arc::new(TriangleMesh::new(
            obj_to_world,
            world_to_obj,
            api_state.graphics_state.reverse_orientation,
            false,  // transform_swaps_handedness
            n_tris, // n_triangles
            vertices,
            n_vertices,
            p_ws,       // in world space
            Vec::new(), // in world space
            n_ws,       // in world space
            uvs,
        ));
        let mtl: Option<Arc<Material + Send + Sync>> = create_material(&api_state, bsdf_state);
        for id in 0..mesh.n_triangles {
            let triangle = Arc::new(Triangle::new(
                mesh.object_to_world,
                mesh.world_to_object,
                mesh.reverse_orientation,
                mesh.clone(),
                id,
            ));
            shapes.push(triangle.clone());
            materials.push(mtl.clone());
        }
    } else {
        panic!("Shape \"{}\" unknown.", api_state.param_set.name);
    }
    (shapes, materials)
}

fn print_params(params: &ParamSet) {
    for p in &params.strings {
        if p.n_values == 1_usize {
            println!("  \"string {}\" [\"{}\"]", p.name, p.values[0]);
        }
    }
    for p in &params.bools {
        if p.n_values == 1_usize {
            println!("  \"bool {}\" [{}]", p.name, p.values[0]);
        } else {
            print!("  \"bool {}\" [ ", p.name);
            for i in 0..p.n_values {
                print!("{} ", p.values[i]);
            }
            println!("]");
        }
    }
    for p in &params.ints {
        if p.n_values == 1_usize {
            println!("  \"integer {}\" [{}]", p.name, p.values[0]);
        } else {
            print!("  \"integer {}\" [ ", p.name);
            for i in 0..p.n_values {
                print!("{} ", p.values[i]);
            }
            println!("]");
        }
    }
    for p in &params.floats {
        if p.n_values == 1_usize {
            println!("  \"float {}\" [{}]", p.name, p.values[0]);
        } else {
            print!("  \"float {}\" [ ", p.name);
            for i in 0..p.n_values {
                print!("{} ", p.values[i]);
            }
            println!("]");
        }
    }
    for p in &params.point3fs {
        if p.n_values == 1_usize {
            println!(
                "  \"point {}\" [{} {} {}]",
                p.name, p.values[0].x, p.values[0].y, p.values[0].z
            );
        } else {
            println!("  \"point {}\" [", p.name);
            for i in 0..p.n_values {
                println!("    {} {} {} ", p.values[i].x, p.values[i].y, p.values[i].z);
            }
            println!("  ]");
        }
    }
    for p in &params.vector3fs {
        if p.n_values == 1_usize {
            println!(
                "  \"vector {}\" [{} {} {}]",
                p.name, p.values[0].x, p.values[0].y, p.values[0].z
            );
        }
    }
    for p in &params.normals {
        if p.n_values == 1_usize {
            println!(
                "  \"normal {}\" [{} {} {}]",
                p.name, p.values[0].x, p.values[0].y, p.values[0].z
            );
        } else {
            println!("  \"normal {}\" [", p.name);
            for i in 0..p.n_values {
                println!("    {} {} {} ", p.values[i].x, p.values[i].y, p.values[i].z);
            }
            println!("  ]");
        }
    }
    for p in &params.spectra {
        if p.n_values == 1_usize {
            println!(
                "  \"rgb {}\" [{} {} {}]",
                p.name, p.values[0].c[0], p.values[0].c[1], p.values[0].c[2]
            );
        }
    }
    for p in &params.textures {
        if p.n_values == 1_usize {
            println!("  \"texture {}\" \"{}\"", p.name, p.values[0]);
        }
    }
}

pub fn pbrt_init(number_of_threads: u8) -> (ApiState, BsdfState) {
    let mut api_state: ApiState = ApiState::default();
    let bsdf_state: BsdfState = BsdfState::default();
    api_state.number_of_threads = number_of_threads;
    (api_state, bsdf_state)
}

pub fn pbrt_cleanup(api_state: &ApiState) {
    // println!("WorldEnd");
    assert!(
        api_state.pushed_graphics_states.len() == 0_usize,
        "Missing end to pbrtAttributeBegin()"
    );
    assert!(
        api_state.pushed_transforms.len() == 0_usize,
        "Missing end to pbrtTransformBegin()"
    );
    // MakeFilter
    let mut some_filter: Option<Arc<Filter + Sync + Send>> = None;
    if api_state.render_options.filter_name == "box" {
        some_filter = Some(BoxFilter::create(&api_state.render_options.filter_params));
    } else if api_state.render_options.filter_name == "gaussian" {
        some_filter = Some(GaussianFilter::create(
            &api_state.render_options.filter_params,
        ));
    } else if api_state.render_options.filter_name == "mitchell" {
        some_filter = Some(MitchellNetravali::create(
            &api_state.render_options.filter_params,
        ));
    } else if api_state.render_options.filter_name == "sinc" {
        println!("TODO: CreateSincFilter");
    } else if api_state.render_options.filter_name == "triangle" {
        some_filter = Some(TriangleFilter::create(
            &api_state.render_options.filter_params,
        ));
    } else {
        panic!(
            "Filter \"{}\" unknown.",
            api_state.render_options.filter_name
        );
    }
    // MakeFilm
    if api_state.render_options.film_name == "image" {
        let filename: String = api_state
            .render_options
            .film_params
            .find_one_string("filename", String::new());
        let xres: i32 = api_state
            .render_options
            .film_params
            .find_one_int("xresolution", 1280);
        let yres: i32 = api_state
            .render_options
            .film_params
            .find_one_int("yresolution", 720);
        // TODO: if (PbrtOptions.quickRender) xres = std::max(1, xres / 4);
        // TODO: if (PbrtOptions.quickRender) yres = std::max(1, yres / 4);
        let mut crop: Bounds2f = Bounds2f {
            p_min: Point2f { x: 0.0, y: 0.0 },
            p_max: Point2f { x: 1.0, y: 1.0 },
        };
        // TODO: const Float *cr = params.FindFloat("cropwindow", &cwi);
        let cr: Vec<Float> = api_state
            .render_options
            .film_params
            .find_float("cropwindow");
        if cr.len() == 4 {
            crop.p_min.x = clamp_t(cr[0].min(cr[1]), 0.0, 1.0);
            crop.p_max.x = clamp_t(cr[0].max(cr[1]), 0.0, 1.0);
            crop.p_min.y = clamp_t(cr[2].min(cr[3]), 0.0, 1.0);
            crop.p_max.y = clamp_t(cr[2].max(cr[3]), 0.0, 1.0);
        } else if cr.len() != 0 {
            panic!(
                "{:?} values supplied for \"cropwindow\". Expected 4.",
                cr.len()
            );
        }
        let scale: Float = api_state
            .render_options
            .film_params
            .find_one_float("scale", 1.0);
        let diagonal: Float = api_state
            .render_options
            .film_params
            .find_one_float("diagonal", 35.0);
        let max_sample_luminance: Float = api_state
            .render_options
            .film_params
            .find_one_float("maxsampleluminance", std::f32::INFINITY);
        if let Some(filter) = some_filter {
            let film: Arc<Film> = Arc::new(Film::new(
                Point2i { x: xres, y: yres },
                crop,
                filter,
                diagonal,
                filename,
                scale,
                max_sample_luminance,
            ));
            // MakeCamera
            // TODO: let mut some_camera: Option<Arc<Camera + Sync + Send>> = None;
            let mut some_camera: Option<Arc<Camera + Sync + Send>>;
            let medium_interface: MediumInterface = create_medium_interface(&api_state);
            let animated_cam_to_world: AnimatedTransform = AnimatedTransform::new(
                &api_state.render_options.camera_to_world.t[0],
                api_state.render_options.transform_start_time,
                &api_state.render_options.camera_to_world.t[1],
                api_state.render_options.transform_end_time,
            );
            if api_state.render_options.camera_name == "perspective" {
                let camera: Arc<Camera + Send + Sync> = PerspectiveCamera::create(
                    &api_state.render_options.camera_params,
                    animated_cam_to_world,
                    film,
                    medium_interface.outside,
                );
                some_camera = Some(camera);
            } else if api_state.render_options.camera_name == "orthographic" {
                let camera: Arc<Camera + Send + Sync> = OrthographicCamera::create(
                    &api_state.render_options.camera_params,
                    animated_cam_to_world,
                    film,
                    medium_interface.outside,
                );
                some_camera = Some(camera);
            } else if api_state.render_options.camera_name == "realistic" {
                if let Some(ref search_directory) = api_state.search_directory {
                    let camera: Arc<Camera + Send + Sync> = RealisticCamera::create(
                        &api_state.render_options.camera_params,
                        animated_cam_to_world,
                        film,
                        medium_interface.outside,
                        // additional parameters:
                        Some(search_directory),
                    );
                    some_camera = Some(camera);
                } else {
                    let camera: Arc<Camera + Send + Sync> = RealisticCamera::create(
                        &api_state.render_options.camera_params,
                        animated_cam_to_world,
                        film,
                        medium_interface.outside,
                        // additional parameters:
                        None,
                    );
                    some_camera = Some(camera);
                }
            } else if api_state.render_options.camera_name == "environment" {
                let camera: Arc<Camera + Send + Sync> = EnvironmentCamera::create(
                    &api_state.render_options.camera_params,
                    animated_cam_to_world,
                    film,
                    medium_interface.outside,
                );
                some_camera = Some(camera);
            } else {
                panic!(
                    "Camera \"{}\" unknown.",
                    api_state.render_options.camera_name
                );
            }
            if let Some(camera) = some_camera {
                // MakeSampler
                let mut some_sampler: Option<Box<Sampler + Sync + Send>> = None;
                if api_state.render_options.sampler_name == "lowdiscrepancy"
                    || api_state.render_options.sampler_name == "02sequence"
                {
                    let nsamp: i32 = api_state
                        .render_options
                        .sampler_params
                        .find_one_int("pixelsamples", 16);
                    let sd: i32 = api_state
                        .render_options
                        .sampler_params
                        .find_one_int("dimensions", 4);
                    // TODO: if (PbrtOptions.quickRender) nsamp = 1;
                    let sampler = Box::new(ZeroTwoSequenceSampler::new(nsamp as i64, sd as i64));
                    some_sampler = Some(sampler);
                } else if api_state.render_options.sampler_name == "maxmindist" {
                    println!("TODO: CreateMaxMinDistSampler");
                } else if api_state.render_options.sampler_name == "halton" {
                    let nsamp: i32 = api_state
                        .render_options
                        .sampler_params
                        .find_one_int("pixelsamples", 16);
                    // TODO: if (PbrtOptions.quickRender) nsamp = 1;
                    let sample_at_center: bool = api_state
                        .render_options
                        .integrator_params
                        .find_one_bool("samplepixelcenter", false);
                    let sample_bounds: Bounds2i = camera.get_film().get_sample_bounds();
                    let sampler = Box::new(HaltonSampler::new(
                        nsamp as i64,
                        sample_bounds,
                        sample_at_center,
                    ));
                    some_sampler = Some(sampler);
                } else if api_state.render_options.sampler_name == "sobol" {
                    let nsamp: i32 = api_state
                        .render_options
                        .sampler_params
                        .find_one_int("pixelsamples", 16);
                    let sample_bounds: Bounds2i = camera.get_film().get_sample_bounds();
                    let sampler = Box::new(SobolSampler::new(nsamp as i64, sample_bounds));
                    some_sampler = Some(sampler);
                } else if api_state.render_options.sampler_name == "random" {
                    let nsamp: i32 = api_state
                        .render_options
                        .sampler_params
                        .find_one_int("pixelsamples", 4);
                    let sampler = Box::new(RandomSampler::new(nsamp as i64));
                    some_sampler = Some(sampler);
                } else if api_state.render_options.sampler_name == "stratified" {
                    println!("TODO: CreateStratifiedSampler");
                } else {
                    panic!(
                        "Sampler \"{}\" unknown.",
                        api_state.render_options.sampler_name
                    );
                }
                if let Some(mut sampler) = some_sampler {
                    // MakeIntegrator
                    // if let Some(mut sampler) = some_sampler {
                    let mut some_integrator: Option<Box<SamplerIntegrator + Sync + Send>> = None;
                    let mut some_bdpt_integrator: Option<Box<BDPTIntegrator>> = None;
                    let mut some_mlt_integrator: Option<Box<MLTIntegrator>> = None;
                    let mut some_sppm_integrator: Option<Box<SPPMIntegrator>> = None;
                    if api_state.render_options.integrator_name == "whitted" {
                        println!("TODO: CreateWhittedIntegrator");
                    } else if api_state.render_options.integrator_name == "directlighting" {
                        // CreateDirectLightingIntegrator
                        let max_depth: i32 = api_state
                            .render_options
                            .integrator_params
                            .find_one_int("maxdepth", 5);
                        let st: String = api_state
                            .render_options
                            .integrator_params
                            .find_one_string("strategy", String::from("all"));
                        let strategy: LightStrategy;
                        if st == "one" {
                            strategy = LightStrategy::UniformSampleOne;
                        } else if st == "all" {
                            strategy = LightStrategy::UniformSampleAll;
                        } else {
                            panic!("Strategy \"{}\" for direct lighting unknown.", st);
                        }
                        // TODO: const int *pb = params.FindInt("pixelbounds", &np);
                        let pixel_bounds: Bounds2i = Bounds2i {
                            p_min: Point2i { x: 0, y: 0 },
                            p_max: Point2i { x: xres, y: yres },
                        };
                        let integrator = Box::new(DirectLightingIntegrator::new(
                            strategy,
                            max_depth as i64,
                            pixel_bounds,
                        ));
                        some_integrator = Some(integrator);
                    } else if api_state.render_options.integrator_name == "path" {
                        // CreatePathIntegrator
                        let max_depth: i32 = api_state
                            .render_options
                            .integrator_params
                            .find_one_int("maxdepth", 5);
                        let pb: Vec<i32> = api_state
                            .render_options
                            .integrator_params
                            .find_int("pixelbounds");
                        let np: usize = pb.len();
                        let pixel_bounds: Bounds2i = camera.get_film().get_sample_bounds();
                        if np > 0 as usize {
                            if np != 4 as usize {
                                panic!(
                                    "Expected four values for \"pixelbounds\" parameter. Got {}.",
                                    np
                                );
                            } else {
                                println!("TODO: pixelBounds = Intersect(...)");
                                // pixelBounds = Intersect(pixelBounds,
                                //                         Bounds2i{{pb[0], pb[2]}, {pb[1], pb[3]}});
                                // if (pixelBounds.Area() == 0)
                                //     Error("Degenerate \"pixelbounds\" specified.");
                            }
                        }
                        let rr_threshold: Float = api_state
                            .render_options
                            .integrator_params
                            .find_one_float("rrthreshold", 1.0 as Float);
                        // std::string lightStrategy =
                        //     params.FindOneString("lightsamplestrategy", "spatial");
                        let light_strategy: String = api_state
                            .render_options
                            .integrator_params
                            .find_one_string("lightsamplestrategy", String::from("spatial"));
                        let integrator = Box::new(PathIntegrator::new(
                            max_depth as u32,
                            pixel_bounds,
                            rr_threshold,
                            light_strategy,
                        ));
                        some_integrator = Some(integrator);
                    } else if api_state.render_options.integrator_name == "volpath" {
                        println!("TODO: CreateVolPathIntegrator");
                    } else if api_state.render_options.integrator_name == "bdpt" {
                        // CreateBDPTIntegrator
                        let mut max_depth: i32 = api_state
                            .render_options
                            .integrator_params
                            .find_one_int("maxdepth", 5);
                        let visualize_strategies: bool = api_state
                            .render_options
                            .integrator_params
                            .find_one_bool("visualizestrategies", false);
                        let visualize_weights: bool = api_state
                            .render_options
                            .integrator_params
                            .find_one_bool("visualizeweights", false);
                        if (visualize_strategies || visualize_weights) && max_depth > 5_i32 {
                            println!("WARNING: visualizestrategies/visualizeweights was enabled, limiting maxdepth to 5");
                            max_depth = 5;
                        }
                        let pixel_bounds: Bounds2i = camera.get_film().get_sample_bounds();
                        let light_strategy: String = api_state
                            .render_options
                            .integrator_params
                            .find_one_string("lightsamplestrategy", String::from("power"));
                        let mut integrator = Box::new(BDPTIntegrator::new(
                            max_depth as u32,
                            // visualize_strategies,
                            // visualize_weights,
                            pixel_bounds,
                            light_strategy,
                        ));
                        some_bdpt_integrator = Some(integrator);
                    } else if api_state.render_options.integrator_name == "mlt" {
                        // CreateMLTIntegrator
                        let mut max_depth: i32 = api_state
                            .render_options
                            .integrator_params
                            .find_one_int("maxdepth", 5);
                        let mut n_bootstrap: i32 = api_state
                            .render_options
                            .integrator_params
                            .find_one_int("bootstrapsamples", 100000);
                        let mut n_chains: i32 = api_state
                            .render_options
                            .integrator_params
                            .find_one_int("chains", 1000);
                        let mut mutations_per_pixel: i32 = api_state
                            .render_options
                            .integrator_params
                            .find_one_int("mutationsperpixel", 100);
                        let large_step_probability: Float = api_state
                            .render_options
                            .integrator_params
                            .find_one_float("largestepprobability", 0.3 as Float);
                        let sigma: Float = api_state
                            .render_options
                            .integrator_params
                            .find_one_float("sigma", 0.01 as Float);
                        let mut integrator = Box::new(MLTIntegrator::new(
                            camera.clone(),
                            max_depth as u32,
                            n_bootstrap as u32,
                            n_chains as u32,
                            mutations_per_pixel as u32,
                            sigma,
                            large_step_probability,
                        ));
                        some_mlt_integrator = Some(integrator);
                    } else if api_state.render_options.integrator_name == "ambientocclusion" {
                        // CreateAOIntegrator
                        let pb: Vec<i32> = api_state
                            .render_options
                            .integrator_params
                            .find_int("pixelbounds");
                        let np: usize = pb.len();
                        let pixel_bounds: Bounds2i = camera.get_film().get_sample_bounds();
                        if np > 0 as usize {
                            if np != 4 as usize {
                                panic!(
                                    "Expected four values for \"pixelbounds\" parameter. Got {}.",
                                    np
                                );
                            } else {
                                println!("TODO: pixelBounds = Intersect(...)");
                                // pixelBounds = Intersect(pixelBounds,
                                //                         Bounds2i{{pb[0], pb[2]}, {pb[1], pb[3]}});
                                // if (pixelBounds.Area() == 0)
                                //     Error("Degenerate \"pixelbounds\" specified.");
                            }
                        }
                        let cos_sample: bool = api_state
                            .render_options
                            .integrator_params
                            .find_one_bool("cossample", true);
                        let n_samples: i32 = api_state
                            .render_options
                            .integrator_params
                            .find_one_int("nsamples", 64 as i32);
                        let integrator =
                            Box::new(AOIntegrator::new(cos_sample, n_samples, pixel_bounds));
                        some_integrator = Some(integrator);
                    } else if api_state.render_options.integrator_name == "sppm" {
                        // CreateSPPMIntegrator
                        let mut n_iterations: i32 = api_state
                            .render_options
                            .integrator_params
                            .find_one_int("numiterations", 64);
                        n_iterations = api_state
                            .render_options
                            .integrator_params
                            .find_one_int("iterations", n_iterations);
                        let max_depth: i32 = api_state
                            .render_options
                            .integrator_params
                            .find_one_int("maxdepth", 5);
                        let photons_per_iter: i32 = api_state
                            .render_options
                            .integrator_params
                            .find_one_int("photonsperiteration", -1);
                        let write_freq: i32 = api_state
                            .render_options
                            .integrator_params
                            .find_one_int("imagewritefrequency", 1 << 31);
                        let radius: Float = api_state
                            .render_options
                            .integrator_params
                            .find_one_float("radius", 1.0 as Float);
                        // TODO: if (PbrtOptions.quickRender) nIterations = std::max(1, nIterations / 16);
                        let mut integrator = Box::new(SPPMIntegrator::new(
                            camera.clone(),
                            n_iterations,
                            photons_per_iter,
                            max_depth as u32,
                            radius,
                            write_freq,
                        ));
                        some_sppm_integrator = Some(integrator);
                    } else {
                        panic!(
                            "Integrator \"{}\" unknown.",
                            api_state.render_options.integrator_name
                        );
                    }
                    if api_state.render_options.have_scattering_media
                        && api_state.render_options.integrator_name != String::from("volpath")
                        && api_state.render_options.integrator_name != String::from("bdpt")
                        && api_state.render_options.integrator_name != String::from("mlt")
                    {
                        print!("WARNING: Scene has scattering media but \"{}\" integrator doesn't support ",
                               api_state.render_options.integrator_name);
                        print!("volume scattering. Consider using \"volpath\", \"bdpt\", or ");
                        println!("\"mlt\".");
                    }
                    if let Some(mut integrator) = some_integrator {
                        // MakeIntegrator
                        // TODO: if (renderOptions->haveScatteringMedia && ...)
                        if api_state.render_options.lights.is_empty() {
                            // warn if no light sources are defined
                            println!("WARNING: No light sources defined in scene; rendering a black image.",);
                        }
                        // MakeAccelerator
                        if api_state.render_options.accelerator_name == "bvh" {
                            //  CreateBVHAccelerator
                            let split_method_name: String = api_state
                                .render_options
                                .accelerator_params
                                .find_one_string("splitmethod", String::from("sah"));
                            let split_method;
                            if split_method_name == "sah" {
                                split_method = SplitMethod::SAH;
                            } else if split_method_name == "hlbvh" {
                                split_method = SplitMethod::HLBVH;
                            } else if split_method_name == "middle" {
                                split_method = SplitMethod::Middle;
                            } else if split_method_name == "equal" {
                                split_method = SplitMethod::EqualCounts;
                            } else {
                                println!(
                                    "WARNING: BVH split method \"{}\" unknown.  Using \"sah\".",
                                    split_method_name
                                );
                                split_method = SplitMethod::SAH;
                            }
                            let max_prims_in_node: i32 = api_state
                                .render_options
                                .accelerator_params
                                .find_one_int("maxnodeprims", 4);
                            let accelerator = Arc::new(BVHAccel::new(
                                api_state.render_options.primitives.clone(),
                                max_prims_in_node as usize,
                                split_method,
                            ));
                            // MakeScene
                            let scene: Scene = Scene::new(
                                accelerator.clone(),
                                api_state.render_options.lights.clone(),
                            );
                            // TODO: primitives.erase(primitives.begin(), primitives.end());
                            // TODO: lights.erase(lights.begin(), lights.end());
                            let num_threads: u8 = api_state.number_of_threads;
                            render(
                                &scene,
                                &camera.clone(),
                                &mut sampler,
                                &mut integrator,
                                num_threads,
                            );
                        } else if api_state.render_options.accelerator_name == "kdtree" {
                            // println!("TODO: CreateKdTreeAccelerator");
                            // WARNING: Use BVHAccel for now !!!
                            let accelerator = Arc::new(BVHAccel::new(
                                api_state.render_options.primitives.clone(),
                                4,
                                SplitMethod::SAH,
                            ));
                            // MakeScene
                            let scene: Scene = Scene::new(
                                accelerator.clone(),
                                api_state.render_options.lights.clone(),
                            );
                            // TODO: primitives.erase(primitives.begin(), primitives.end());
                            // TODO: lights.erase(lights.begin(), lights.end());
                            let num_threads: u8 = api_state.number_of_threads;
                            render(&scene, &camera, &mut sampler, &mut integrator, num_threads);
                        } else {
                            panic!(
                                "Accelerator \"{}\" unknown.",
                                api_state.render_options.accelerator_name
                            );
                        }
                    } else if let Some(mut integrator) = some_bdpt_integrator {
                        // because we can't call
                        // integrator.render() yet,
                        // let us repeat some code and
                        // call render_bdpt(...)
                        // instead:

                        // MakeIntegrator
                        // TODO: if (renderOptions->haveScatteringMedia && ...)
                        if api_state.render_options.lights.is_empty() {
                            // warn if no light sources are defined
                            println!("WARNING: No light sources defined in scene; rendering a black image.",);
                        }
                        // MakeAccelerator
                        if api_state.render_options.accelerator_name == "bvh" {
                            //  CreateBVHAccelerator
                            let split_method_name: String = api_state
                                .render_options
                                .accelerator_params
                                .find_one_string("splitmethod", String::from("sah"));
                            let split_method;
                            if split_method_name == "sah" {
                                split_method = SplitMethod::SAH;
                            } else if split_method_name == "hlbvh" {
                                split_method = SplitMethod::HLBVH;
                            } else if split_method_name == "middle" {
                                split_method = SplitMethod::Middle;
                            } else if split_method_name == "equal" {
                                split_method = SplitMethod::EqualCounts;
                            } else {
                                println!(
                                    "WARNING: BVH split method \"{}\" unknown.  Using \"sah\".",
                                    split_method_name
                                );
                                split_method = SplitMethod::SAH;
                            }
                            let max_prims_in_node: i32 = api_state
                                .render_options
                                .accelerator_params
                                .find_one_int("maxnodeprims", 4);
                            let accelerator = Arc::new(BVHAccel::new(
                                api_state.render_options.primitives.clone(),
                                max_prims_in_node as usize,
                                split_method,
                            ));
                            // MakeScene
                            let scene: Scene = Scene::new(
                                accelerator.clone(),
                                api_state.render_options.lights.clone(),
                            );
                            // TODO: primitives.erase(primitives.begin(), primitives.end());
                            // TODO: lights.erase(lights.begin(), lights.end());
                            let num_threads: u8 = api_state.number_of_threads;
                            render_bdpt(
                                &scene,
                                &camera,
                                &mut sampler,
                                &mut integrator,
                                num_threads,
                            );
                        } else if api_state.render_options.accelerator_name == "kdtree" {
                            // println!("TODO: CreateKdTreeAccelerator");
                            // WARNING: Use BVHAccel for now !!!
                            let accelerator = Arc::new(BVHAccel::new(
                                api_state.render_options.primitives.clone(),
                                4,
                                SplitMethod::SAH,
                            ));
                            // MakeScene
                            let scene: Scene = Scene::new(
                                accelerator.clone(),
                                api_state.render_options.lights.clone(),
                            );
                            // TODO: primitives.erase(primitives.begin(), primitives.end());
                            // TODO: lights.erase(lights.begin(), lights.end());
                            let num_threads: u8 = api_state.number_of_threads;
                            render_bdpt(
                                &scene,
                                &camera,
                                &mut sampler,
                                &mut integrator,
                                num_threads,
                            );
                        } else {
                            panic!(
                                "Accelerator \"{}\" unknown.",
                                api_state.render_options.accelerator_name
                            );
                        }
                    } else if let Some(mut integrator) = some_mlt_integrator {
                        // because we can't call
                        // integrator.render() yet,
                        // let us repeat some code and
                        // call render_mlt(...)
                        // instead:

                        // MakeIntegrator
                        // TODO: if (renderOptions->haveScatteringMedia && ...)
                        if api_state.render_options.lights.is_empty() {
                            // warn if no light sources are defined
                            println!("WARNING: No light sources defined in scene; rendering a black image.",);
                        }
                        // MakeAccelerator
                        if api_state.render_options.accelerator_name == "bvh" {
                            //  CreateBVHAccelerator
                            let split_method_name: String = api_state
                                .render_options
                                .accelerator_params
                                .find_one_string("splitmethod", String::from("sah"));
                            let split_method;
                            if split_method_name == "sah" {
                                split_method = SplitMethod::SAH;
                            } else if split_method_name == "hlbvh" {
                                split_method = SplitMethod::HLBVH;
                            } else if split_method_name == "middle" {
                                split_method = SplitMethod::Middle;
                            } else if split_method_name == "equal" {
                                split_method = SplitMethod::EqualCounts;
                            } else {
                                println!(
                                    "WARNING: BVH split method \"{}\" unknown.  Using \"sah\".",
                                    split_method_name
                                );
                                split_method = SplitMethod::SAH;
                            }
                            let max_prims_in_node: i32 = api_state
                                .render_options
                                .accelerator_params
                                .find_one_int("maxnodeprims", 4);
                            let accelerator = Arc::new(BVHAccel::new(
                                api_state.render_options.primitives.clone(),
                                max_prims_in_node as usize,
                                split_method,
                            ));
                            // MakeScene
                            let scene: Scene = Scene::new(
                                accelerator.clone(),
                                api_state.render_options.lights.clone(),
                            );
                            // TODO: primitives.erase(primitives.begin(), primitives.end());
                            // TODO: lights.erase(lights.begin(), lights.end());
                            let num_threads: u8 = api_state.number_of_threads;
                            render_mlt(&scene, &camera, &mut sampler, &mut integrator, num_threads);
                        } else if api_state.render_options.accelerator_name == "kdtree" {
                            // println!("TODO: CreateKdTreeAccelerator");
                            // WARNING: Use BVHAccel for now !!!
                            let accelerator = Arc::new(BVHAccel::new(
                                api_state.render_options.primitives.clone(),
                                4,
                                SplitMethod::SAH,
                            ));
                            // MakeScene
                            let scene: Scene = Scene::new(
                                accelerator.clone(),
                                api_state.render_options.lights.clone(),
                            );
                            // TODO: primitives.erase(primitives.begin(), primitives.end());
                            // TODO: lights.erase(lights.begin(), lights.end());
                            let num_threads: u8 = api_state.number_of_threads;
                            render_mlt(&scene, &camera, &mut sampler, &mut integrator, num_threads);
                        } else {
                            panic!(
                                "Accelerator \"{}\" unknown.",
                                api_state.render_options.accelerator_name
                            );
                        }
                    } else if let Some(mut integrator) = some_sppm_integrator {
                        // because we can't call
                        // integrator.render() yet,
                        // let us repeat some code and
                        // call render_sppm(...)
                        // instead:

                        // MakeIntegrator
                        // TODO: if (renderOptions->haveScatteringMedia && ...)
                        if api_state.render_options.lights.is_empty() {
                            // warn if no light sources are defined
                            println!("WARNING: No light sources defined in scene; rendering a black image.",);
                        }
                        // MakeAccelerator
                        if api_state.render_options.accelerator_name == "bvh" {
                            //  CreateBVHAccelerator
                            let split_method_name: String = api_state
                                .render_options
                                .accelerator_params
                                .find_one_string("splitmethod", String::from("sah"));
                            let split_method;
                            if split_method_name == "sah" {
                                split_method = SplitMethod::SAH;
                            } else if split_method_name == "hlbvh" {
                                split_method = SplitMethod::HLBVH;
                            } else if split_method_name == "middle" {
                                split_method = SplitMethod::Middle;
                            } else if split_method_name == "equal" {
                                split_method = SplitMethod::EqualCounts;
                            } else {
                                println!(
                                    "WARNING: BVH split method \"{}\" unknown.  Using \"sah\".",
                                    split_method_name
                                );
                                split_method = SplitMethod::SAH;
                            }
                            let max_prims_in_node: i32 = api_state
                                .render_options
                                .accelerator_params
                                .find_one_int("maxnodeprims", 4);
                            let accelerator = Arc::new(BVHAccel::new(
                                api_state.render_options.primitives.clone(),
                                max_prims_in_node as usize,
                                split_method,
                            ));
                            // MakeScene
                            let scene: Scene = Scene::new(
                                accelerator.clone(),
                                api_state.render_options.lights.clone(),
                            );
                            // TODO: primitives.erase(primitives.begin(), primitives.end());
                            // TODO: lights.erase(lights.begin(), lights.end());
                            let num_threads: u8 = api_state.number_of_threads;
                            render_sppm(
                                &scene,
                                &camera,
                                &mut sampler,
                                &mut integrator,
                                num_threads,
                            );
                        } else if api_state.render_options.accelerator_name == "kdtree" {
                            // println!("TODO: CreateKdTreeAccelerator");
                            // WARNING: Use BVHAccel for now !!!
                            let accelerator = Arc::new(BVHAccel::new(
                                api_state.render_options.primitives.clone(),
                                4,
                                SplitMethod::SAH,
                            ));
                            // MakeScene
                            let scene: Scene = Scene::new(
                                accelerator.clone(),
                                api_state.render_options.lights.clone(),
                            );
                            // TODO: primitives.erase(primitives.begin(), primitives.end());
                            // TODO: lights.erase(lights.begin(), lights.end());
                            let num_threads: u8 = api_state.number_of_threads;
                            render_sppm(
                                &scene,
                                &camera,
                                &mut sampler,
                                &mut integrator,
                                num_threads,
                            );
                        } else {
                            panic!(
                                "Accelerator \"{}\" unknown.",
                                api_state.render_options.accelerator_name
                            );
                        }
                    } else {
                        panic!("Unable to create integrator.");
                    }
                } else {
                    panic!("Unable to create sampler.");
                }
            } else {
                panic!("Unable to create camera.");
            }
        } else {
            panic!("Unable to create film.");
        }
    } else {
        panic!("Film \"{}\" unknown.", api_state.render_options.film_name);
    }
}

pub fn pbrt_translate(api_state: &mut ApiState, dx: Float, dy: Float, dz: Float) {
    // println!("Translate {} {} {}", dx, dy, dz);
    let translate: Transform = Transform::translate(&Vector3f {
        x: dx,
        y: dy,
        z: dz,
    });
    if api_state.active_transform_bits & 1_u8 > 0_u8 {
        // 0x?1
        api_state.cur_transform.t[0] = api_state.cur_transform.t[0] * translate;
    }
    if api_state.active_transform_bits & 2_u8 > 0_u8 {
        // 0x1?
        api_state.cur_transform.t[1] = api_state.cur_transform.t[1] * translate;
    }
}

pub fn pbrt_transform(api_state: &mut ApiState, tr: &Transform) {
    // println!("{:?}", tr);
    if api_state.active_transform_bits & 1_u8 > 0_u8 {
        // 0x?1
        api_state.cur_transform.t[0] = *tr;
    }
    if api_state.active_transform_bits & 2_u8 > 0_u8 {
        // 0x1?
        api_state.cur_transform.t[1] = *tr;
    }
}

pub fn pbrt_concat_transform(api_state: &mut ApiState, tr: &Transform) {
    // println!("Concat{:?}", tr);
    if api_state.active_transform_bits & 1_u8 > 0_u8 {
        // 0x?1
        api_state.cur_transform.t[0] = api_state.cur_transform.t[0] * *tr;
    }
    if api_state.active_transform_bits & 2_u8 > 0_u8 {
        // 0x1?
        api_state.cur_transform.t[1] = api_state.cur_transform.t[1] * *tr;
    }
}

pub fn pbrt_rotate(api_state: &mut ApiState, angle: Float, dx: Float, dy: Float, dz: Float) {
    // println!("Rotate {} {} {} {}", angle, dx, dy, dz);
    let rotate: Transform = Transform::rotate(
        angle,
        &Vector3f {
            x: dx,
            y: dy,
            z: dz,
        },
    );
    if api_state.active_transform_bits & 1_u8 > 0_u8 {
        // 0x?1
        api_state.cur_transform.t[0] = api_state.cur_transform.t[0] * rotate;
    }
    if api_state.active_transform_bits & 2_u8 > 0_u8 {
        // 0x1?
        api_state.cur_transform.t[1] = api_state.cur_transform.t[1] * rotate;
    }
}

pub fn pbrt_scale(api_state: &mut ApiState, sx: Float, sy: Float, sz: Float) {
    // println!("Scale {} {} {}", sx, sy, sz);
    let scale: Transform = Transform::scale(sx, sy, sz);
    if api_state.active_transform_bits & 1_u8 > 0_u8 {
        // 0x?1
        api_state.cur_transform.t[0] = api_state.cur_transform.t[0] * scale;
    }
    if api_state.active_transform_bits & 2_u8 > 0_u8 {
        // 0x1?
        api_state.cur_transform.t[1] = api_state.cur_transform.t[1] * scale;
    }
}

pub fn pbrt_look_at(
    api_state: &mut ApiState,
    ex: Float,
    ey: Float,
    ez: Float,
    lx: Float,
    ly: Float,
    lz: Float,
    ux: Float,
    uy: Float,
    uz: Float,
) {
    // println!(
    //     "LookAt {} {} {} {} {} {} {} {} {}",
    //     ex, ey, ez, lx, ly, lz, ux, uy, uz
    // );
    let pos: Point3f = Point3f {
        x: ex,
        y: ey,
        z: ez,
    };
    let look: Point3f = Point3f {
        x: lx,
        y: ly,
        z: lz,
    };
    let up: Vector3f = Vector3f {
        x: ux,
        y: uy,
        z: uz,
    };
    let look_at: Transform = Transform::look_at(&pos, &look, &up);
    if api_state.active_transform_bits & 1_u8 > 0_u8 {
        // 0x?1
        api_state.cur_transform.t[0] = api_state.cur_transform.t[0] * look_at;
    }
    if api_state.active_transform_bits & 2_u8 > 0_u8 {
        // 0x1?
        api_state.cur_transform.t[1] = api_state.cur_transform.t[1] * look_at;
    }
}

pub fn pbrt_coord_sys_transform(api_state: &mut ApiState, params: ParamSet) {
    // println!("CoordSysTransform \"{}\"", params.name);
    api_state.param_set = params;
    match api_state
        .named_coordinate_systems
        .get(api_state.param_set.name.as_str())
    {
        Some(transform_set) => {
            api_state.cur_transform.t[0] = transform_set.t[0];
            api_state.cur_transform.t[1] = transform_set.t[1];
        }
        None => {
            println!(
                "Couldn't find named coordinate system \"{}\"",
                api_state.param_set.name
            );
        }
    };
}

pub fn pbrt_active_transform_all(api_state: &mut ApiState) {
    // println!("ActiveTransform All");
    api_state.active_transform_bits = 3_u8 // 0x11
}

pub fn pbrt_active_transform_end_time(api_state: &mut ApiState) {
    // println!("ActiveTransform EndTime");
    api_state.active_transform_bits = 2_u8 // 0x10
}

pub fn pbrt_active_transform_start_time(api_state: &mut ApiState) {
    // println!("ActiveTransform StartTime");
    api_state.active_transform_bits = 1_u8 // 0x01
}

pub fn pbrt_transform_times(api_state: &mut ApiState, start: Float, end: Float) {
    println!("TransformTimes {} {}", start, end);
    api_state.render_options.transform_start_time = start;
    api_state.render_options.transform_end_time = end;
}

pub fn pbrt_pixel_filter(api_state: &mut ApiState, params: ParamSet) {
    // println!("PixelFilter \"{}\"", params.name);
    print_params(&params);
    api_state.render_options.filter_name = params.name.clone();
    api_state.param_set = params;
    api_state
        .render_options
        .filter_params
        .copy_from(&api_state.param_set);
}

pub fn pbrt_film(api_state: &mut ApiState, params: ParamSet) {
    println!("Film \"{}\"", params.name);
    print_params(&params);
    api_state.render_options.film_name = params.name.clone();
    api_state.param_set = params;
    api_state
        .render_options
        .film_params
        .copy_from(&api_state.param_set);
}

pub fn pbrt_sampler(api_state: &mut ApiState, params: ParamSet) {
    println!("Sampler \"{}\"", params.name);
    print_params(&params);
    api_state.render_options.sampler_name = params.name.clone();
    api_state.param_set = params;
    api_state
        .render_options
        .sampler_params
        .copy_from(&api_state.param_set);
}

pub fn pbrt_integrator(api_state: &mut ApiState, params: ParamSet) {
    println!("Integrator \"{}\"", params.name);
    print_params(&params);
    api_state.render_options.integrator_name = params.name.clone();
    api_state.param_set = params;
    api_state
        .render_options
        .integrator_params
        .copy_from(&api_state.param_set);
}

pub fn pbrt_camera(api_state: &mut ApiState, params: ParamSet) {
    // println!("Camera \"{}\"", params.name);
    // print_params(&params);
    api_state.render_options.camera_name = params.name.clone();
    api_state.param_set = params;
    api_state.render_options.camera_to_world.t[0] =
        Transform::inverse(&api_state.cur_transform.t[0]);
    api_state.render_options.camera_to_world.t[1] =
        Transform::inverse(&api_state.cur_transform.t[1]);
    api_state.named_coordinate_systems.insert(
        "camera",
        TransformSet {
            t: [
                api_state.render_options.camera_to_world.t[0],
                api_state.render_options.camera_to_world.t[1],
            ],
        },
    );
    api_state
        .render_options
        .camera_params
        .copy_from(&api_state.param_set);
}

pub fn pbrt_make_named_medium(api_state: &mut ApiState, params: ParamSet) {
    // println!("MakeNamedMedium \"{}\"", params.name);
    print_params(&api_state.param_set);
    api_state.param_set = params;
    make_medium(api_state);
}

pub fn pbrt_medium_interface(
    api_state: &mut ApiState,
    inside_name: &String,
    outside_name: &String,
) {
    // println!("MediumInterface \"{}\" \"{}\"", inside_name, outside_name);
    api_state.graphics_state.current_inside_medium = inside_name.clone();
    api_state.graphics_state.current_outside_medium = outside_name.clone();
    api_state.render_options.have_scattering_media = true;
}

pub fn pbrt_world_begin(api_state: &mut ApiState) {
    // println!("WorldBegin");
    api_state.cur_transform.t[0] = Transform::default();
    api_state.cur_transform.t[1] = Transform::default();
    api_state.active_transform_bits = 3_u8; // 0x11
    api_state.named_coordinate_systems.insert(
        "world",
        TransformSet {
            t: [Transform::default(); 2],
        },
    );
}

pub fn pbrt_attribute_begin(api_state: &mut ApiState) {
    // println!("AttributeBegin");
    let mut material_param_set: ParamSet = ParamSet::default();
    material_param_set.copy_from(&api_state.graphics_state.material_params);
    let mut area_light_param_set: ParamSet = ParamSet::default();
    area_light_param_set.copy_from(&api_state.graphics_state.area_light_params);
    api_state.pushed_graphics_states.push(GraphicsState {
        current_inside_medium: api_state.graphics_state.current_inside_medium.clone(),
        current_outside_medium: api_state.graphics_state.current_outside_medium.clone(),
        float_textures: api_state.graphics_state.float_textures.clone(),
        spectrum_textures: api_state.graphics_state.spectrum_textures.clone(),
        material_params: material_param_set,
        material: String::from(api_state.graphics_state.material.as_ref()),
        named_materials: api_state.graphics_state.named_materials.clone(),
        current_material: String::from(api_state.graphics_state.current_material.as_ref()),
        area_light_params: area_light_param_set,
        area_light: String::from(api_state.graphics_state.area_light.as_ref()),
        reverse_orientation: api_state.graphics_state.reverse_orientation,
    });
    api_state.pushed_transforms.push(TransformSet {
        t: [
            Transform {
                m: api_state.cur_transform.t[0].m,
                m_inv: api_state.cur_transform.t[0].m_inv,
            },
            Transform {
                m: api_state.cur_transform.t[1].m,
                m_inv: api_state.cur_transform.t[1].m_inv,
            },
        ],
    });
    api_state
        .pushed_active_transform_bits
        .push(api_state.active_transform_bits);
}

pub fn pbrt_attribute_end(api_state: &mut ApiState) {
    // println!("AttributeEnd");
    if !(api_state.pushed_graphics_states.len() >= 1_usize) {
        panic!("Unmatched pbrtAttributeEnd() encountered.")
    }
    let pgs: GraphicsState = api_state.pushed_graphics_states.pop().unwrap();
    // current_inside_medium
    api_state.graphics_state.current_inside_medium =
        String::from(pgs.current_inside_medium.as_ref());
    // current_outside_medium
    api_state.graphics_state.current_outside_medium =
        String::from(pgs.current_outside_medium.as_ref());
    // material_params
    api_state.graphics_state.material_params.reset(
        String::new(),
        String::from(""),
        String::from(""),
        String::new(),
    );
    api_state
        .graphics_state
        .material_params
        .copy_from(&pgs.material_params);
    // material
    api_state.graphics_state.material = String::from(pgs.material.as_ref());
    // area_light_params
    api_state.graphics_state.area_light_params.reset(
        String::new(),
        String::from(""),
        String::from(""),
        String::new(),
    );
    api_state
        .graphics_state
        .area_light_params
        .copy_from(&pgs.area_light_params);
    // area_light
    api_state.graphics_state.area_light = String::from(pgs.area_light.as_ref());
    // reverse_orientation
    api_state.graphics_state.reverse_orientation = pgs.reverse_orientation;
    let popped_transform_set: TransformSet = api_state.pushed_transforms.pop().unwrap();
    api_state.cur_transform.t[0] = popped_transform_set.t[0];
    api_state.cur_transform.t[1] = popped_transform_set.t[1];
    let active_transform_bits: u8 = api_state.pushed_active_transform_bits.pop().unwrap();
    api_state.active_transform_bits = active_transform_bits;
}

pub fn pbrt_transform_begin(api_state: &mut ApiState) {
    // println!("TransformBegin");
    api_state.pushed_transforms.push(TransformSet {
        t: [
            Transform {
                m: api_state.cur_transform.t[0].m,
                m_inv: api_state.cur_transform.t[0].m_inv,
            },
            Transform {
                m: api_state.cur_transform.t[1].m,
                m_inv: api_state.cur_transform.t[1].m_inv,
            },
        ],
    });
    api_state
        .pushed_active_transform_bits
        .push(api_state.active_transform_bits);
}

pub fn pbrt_transform_end(api_state: &mut ApiState) {
    // println!("TransformEnd");
    let popped_transform_set: TransformSet = api_state.pushed_transforms.pop().unwrap();
    api_state.cur_transform.t[0] = popped_transform_set.t[0];
    api_state.cur_transform.t[1] = popped_transform_set.t[1];
    let active_transform_bits: u8 = api_state.pushed_active_transform_bits.pop().unwrap();
    api_state.active_transform_bits = active_transform_bits;
}

pub fn pbrt_texture(api_state: &mut ApiState, params: ParamSet) {
    // println!(
    //     "Texture \"{}\" \"{}\" \"{}\"",
    //     params.name, params.tex_type, params.tex_name
    // );
    // print_params(&params);
    api_state.param_set = params;
    make_texture(api_state);
}

pub fn pbrt_material(api_state: &mut ApiState, params: ParamSet) {
    // println!("MakeMaterial \"{}\"", params.name);
    // print_params(&params);
    api_state.param_set = params;
    api_state.graphics_state.material = api_state.param_set.name.clone();
    api_state
        .graphics_state
        .material_params
        .copy_from(&api_state.param_set);
    api_state.graphics_state.current_material = String::new();
}

pub fn pbrt_make_named_material(
    api_state: &mut ApiState,
    bsdf_state: &mut BsdfState,
    params: ParamSet,
) {
    // println!("MakeNamedMaterial \"{}\"", params.name);
    // print_params(&params);
    api_state.param_set = params;
    let mat_type: String = api_state.param_set.find_one_string("type", String::new());
    if mat_type == "" {
        panic!("No parameter string \"type\" found in MakeNamedMaterial");
    }
    api_state.graphics_state.material = mat_type.clone();
    api_state
        .graphics_state
        .material_params
        .copy_from(&api_state.param_set);
    api_state.graphics_state.current_material = String::new();
    let mtl: Option<Arc<Material + Send + Sync>> = create_material(&api_state, bsdf_state);
    match api_state
        .graphics_state
        .named_materials
        .get(api_state.param_set.name.as_str())
    {
        Some(_named_material) => {
            println!("Named material \"{}\" redefined", mat_type);
        }
        None => {}
    }
    api_state
        .graphics_state
        .named_materials
        .insert(api_state.param_set.name.clone(), mtl);
}

pub fn pbrt_named_material(api_state: &mut ApiState, params: ParamSet) {
    // println!("NamedMaterial \"{}\"", params.name);
    api_state.param_set = params;
    api_state.graphics_state.current_material = api_state.param_set.name.clone();
}

pub fn pbrt_light_source(api_state: &mut ApiState, params: ParamSet) {
    // println!("LightSource \"{}\"", params.name);
    // print_params(&params);
    api_state.param_set = params;
    let mi: MediumInterface = create_medium_interface(&api_state);
    make_light(api_state, &mi);
}

pub fn pbrt_area_light_source(api_state: &mut ApiState, params: ParamSet) {
    // println!("AreaLightSource \"{}\"", params.name);
    // print_params(&params);
    api_state.param_set = params;
    api_state.graphics_state.area_light = api_state.param_set.name.clone();
    api_state
        .graphics_state
        .area_light_params
        .copy_from(&api_state.param_set);
}

pub fn pbrt_shape(api_state: &mut ApiState, bsdf_state: &mut BsdfState, params: ParamSet) {
    // println!("Shape \"{}\"", params.name);
    // print_params(&params);
    api_state.param_set = params;
    // collect area lights
    let mut prims: Vec<Arc<Primitive + Send + Sync>> = Vec::new();
    let mut area_lights: Vec<Arc<Light + Send + Sync>> = Vec::new();
    // possibly create area light for shape (see pbrtShape())
    if api_state.graphics_state.area_light != String::new() {
        // MakeAreaLight
        if api_state.graphics_state.area_light == "area"
            || api_state.graphics_state.area_light == "diffuse"
        {
            // first create the shape
            let (shapes, materials) = get_shapes_and_materials(&api_state, bsdf_state);
            assert_eq!(shapes.len(), materials.len());
            // MediumInterface
            let mi: MediumInterface = create_medium_interface(&api_state);
            for i in 0..shapes.len() {
                let shape = &shapes[i];
                let material = &materials[i];
                // CreateDiffuseAreaLight
                let light_to_world: Transform = api_state.cur_transform.t[0];
                let l: Spectrum = api_state
                    .graphics_state
                    .area_light_params
                    .find_one_spectrum("L", Spectrum::new(1.0));
                let sc: Spectrum = api_state
                    .graphics_state
                    .area_light_params
                    .find_one_spectrum("scale", Spectrum::new(1.0));
                let n_samples: i32 = // try "nsamples" first
                    api_state.graphics_state.area_light_params.find_one_int("nsamples",
                                                                  1);
                let n_samples: i32 = // try "samples"next
                    api_state.graphics_state.area_light_params.find_one_int("samples",
                                                                  n_samples);
                let two_sided: bool = api_state
                    .graphics_state
                    .area_light_params
                    .find_one_bool("twosided", false);
                // TODO: if (PbrtOptions.quickRender) nSamples = std::max(1, nSamples / 4);
                let l_emit: Spectrum = l * sc;
                let area_light: Arc<DiffuseAreaLight> = Arc::new(DiffuseAreaLight::new(
                    &light_to_world,
                    &mi,
                    &l_emit,
                    n_samples,
                    shape.clone(),
                    two_sided,
                ));
                area_lights.push(area_light.clone());
                let geo_prim = Arc::new(GeometricPrimitive::new(
                    shape.clone(),
                    material.clone(),
                    Some(area_light.clone()),
                    Some(Arc::new(mi.clone())),
                ));
                prims.push(geo_prim.clone());
            }
        }
    } else {
        // continue with shape itself
        let (shapes, materials) = get_shapes_and_materials(&api_state, bsdf_state);
        assert_eq!(shapes.len(), materials.len());
        // MediumInterface
        let mi: MediumInterface = create_medium_interface(&api_state);
        for i in 0..shapes.len() {
            let shape = &shapes[i];
            let material = &materials[i];
            let geo_prim = Arc::new(GeometricPrimitive::new(
                shape.clone(),
                material.clone(),
                None,
                Some(Arc::new(mi.clone())),
            ));
            prims.push(geo_prim.clone());
        }
        // animated?
        if api_state.cur_transform.is_animated() {
            let animated_object_to_world: AnimatedTransform = AnimatedTransform::new(
                &api_state.cur_transform.t[0],
                api_state.render_options.transform_start_time,
                &api_state.cur_transform.t[1],
                api_state.render_options.transform_end_time,
            );
            if prims.len() > 1 {
                let bvh: Arc<Primitive + Send + Sync> =
                    Arc::new(BVHAccel::new(prims.clone(), 4, SplitMethod::SAH));
                prims.clear();
                prims.push(bvh.clone());
            }
            if let Some(primitive) = prims.pop() {
                let geo_prim = Arc::new(TransformedPrimitive::new(
                    primitive,
                    animated_object_to_world,
                ));
                prims.push(geo_prim.clone());
            }
        }
    }
    // add _prims_ and _areaLights_ to scene or current instance
    if api_state.render_options.current_instance != String::from("") {
        if area_lights.len() > 0 {
            println!("WARNING: Area lights not supported with object instancing");
        }
        if let Some(instance_vec) = api_state
            .render_options
            .instances
            .get_mut(&api_state.render_options.current_instance.clone())
        {
            for prim in prims {
                instance_vec.push(prim.clone());
            }
        }
    } else {
        for prim in prims {
            api_state.render_options.primitives.push(prim.clone());
        }
        if area_lights.len() > 0 {
            for area_light in area_lights {
                api_state.render_options.lights.push(area_light);
            }
        }
    }
}

pub fn pbrt_reverse_orientation(api_state: &mut ApiState) {
    // println!("ReverseOrientation");
    api_state.graphics_state.reverse_orientation = !api_state.graphics_state.reverse_orientation;
}

pub fn pbrt_object_begin(api_state: &mut ApiState, params: ParamSet) {
    // println!("ObjectBegin \"{}\"", params.name);
    api_state.param_set = params;
    pbrt_attribute_begin(api_state);
    if api_state.render_options.current_instance != String::from("") {
        println!("ERROR: ObjectBegin called inside of instance definition");
    }
    api_state
        .render_options
        .instances
        .insert(api_state.param_set.name.clone(), Vec::new());
    api_state.render_options.current_instance = api_state.param_set.name.clone();
}

pub fn pbrt_object_end(api_state: &mut ApiState) {
    // println!("ObjectEnd");
    if api_state.render_options.current_instance == "" {
        println!("ERROR: ObjectEnd called outside of instance definition");
    }
    api_state.render_options.current_instance = String::from("");
    pbrt_attribute_end(api_state);
}

pub fn pbrt_object_instance(api_state: &mut ApiState, params: ParamSet) {
    // println!("ObjectInstance \"{}\"", params.name);
    api_state.param_set = params;
    // perform object instance error checking
    if api_state.render_options.current_instance != String::from("") {
        println!("ERROR: ObjectInstance can't be called inside instance definition");
        return;
    }
    if let Some(instance_vec) = api_state
        .render_options
        .instances
        .get_mut(&api_state.param_set.name.clone())
    {
        if instance_vec.is_empty() {
            return;
        }
        // TODO: ++nObjectInstancesUsed;
        if instance_vec.len() > 1_usize {
            // create aggregate for instance _Primitive_s
            if api_state.render_options.accelerator_name == "bvh" {
                //  CreateBVHAccelerator
                let split_method_name: String = api_state
                    .render_options
                    .accelerator_params
                    .find_one_string("splitmethod", String::from("sah"));
                let split_method;
                if split_method_name == "sah" {
                    split_method = SplitMethod::SAH;
                } else if split_method_name == "hlbvh" {
                    split_method = SplitMethod::HLBVH;
                } else if split_method_name == "middle" {
                    split_method = SplitMethod::Middle;
                } else if split_method_name == "equal" {
                    split_method = SplitMethod::EqualCounts;
                } else {
                    println!(
                        "WARNING: BVH split method \"{}\" unknown.  Using \"sah\".",
                        split_method_name
                    );
                    split_method = SplitMethod::SAH;
                }
                let max_prims_in_node: i32 = api_state
                    .render_options
                    .accelerator_params
                    .find_one_int("maxnodeprims", 4);
                let accelerator: Arc<Primitive + Sync + Send> = Arc::new(BVHAccel::new(
                    instance_vec.clone(),
                    max_prims_in_node as usize,
                    split_method,
                ));
                instance_vec.clear();
                instance_vec.push(accelerator);
            } else if api_state.render_options.accelerator_name == "kdtree" {
                // println!("TODO: CreateKdTreeAccelerator");
                // WARNING: Use BVHAccel for now !!!
                let accelerator: Arc<Primitive + Sync + Send> =
                    Arc::new(BVHAccel::new(instance_vec.clone(), 4, SplitMethod::SAH));
                instance_vec.clear();
                instance_vec.push(accelerator);
            } else {
                panic!(
                    "Accelerator \"{}\" unknown.",
                    api_state.render_options.accelerator_name
                );
            }
        }
        // create _animatedInstanceToWorld_ transform for instance
        let animated_instance_to_world: AnimatedTransform = AnimatedTransform::new(
            &api_state.cur_transform.t[0],
            api_state.render_options.transform_start_time,
            &api_state.cur_transform.t[1],
            api_state.render_options.transform_end_time,
        );
        let prim: Arc<Primitive + Send + Sync> = Arc::new(TransformedPrimitive::new(
            instance_vec[0].clone(),
            animated_instance_to_world,
        ));
        api_state.render_options.primitives.push(prim.clone());
    } else {
        println!(
            "ERROR: Unable to find instance named {:?}",
            api_state.param_set.name.clone()
        );
        return;
    }
}
