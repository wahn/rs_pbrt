// std
use std::collections::HashMap;
use std::sync::Arc;
// pbrt
use core::geometry::{Point3f, Vector3f};
use core::light::Light;
use core::material::Material;
use core::medium::Medium;
use core::paramset::{ParamSet, TextureParams};
use core::pbrt::{Float, Spectrum};
use core::primitive::Primitive;
use core::texture::Texture;
use core::transform::{Matrix4x4, Transform};
use materials::glass::GlassMaterial;
use materials::hair::HairMaterial;
use materials::matte::MatteMaterial;
use materials::metal::MetalMaterial;
use materials::mirror::MirrorMaterial;
use materials::mixmat::MixMaterial;
use materials::plastic::PlasticMaterial;
use materials::substrate::SubstrateMaterial;
use materials::uber::UberMaterial;
use textures::constant::ConstantTexture;

// see api.cpp

pub struct ApiState {
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
    // TODO: std::map<std::string, std::vector<std::shared_ptr<Primitive>>> instances;
    // TODO: std::vector<std::shared_ptr<Primitive>> *currentInstance = nullptr;
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
            integrator_name: String::from("image"),
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

fn create_material(api_state: &ApiState) -> Option<Arc<Material + Send + Sync>> {
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
        if api_state.graphics_state.material == String::from("")
            || api_state.graphics_state.material == String::from("none")
        {
            return None;
        } else if api_state.graphics_state.material == String::from("matte") {
            return Some(MatteMaterial::create(&mut mp));
        } else if api_state.graphics_state.material == String::from("plastic") {
            let kd = mp.get_spectrum_texture(String::from("Kd"), Spectrum::new(0.25 as Float));
            let ks = mp.get_spectrum_texture(String::from("Ks"), Spectrum::new(0.25 as Float));
            let roughness = mp.get_float_texture(String::from("roughness"), 0.1 as Float);
            // TODO: std::shared_ptr<Texture<Float>> bumpMap = mp.GetFloatTextureOrNull("bumpmap");
            let remap_roughness: bool = mp.find_bool(String::from("remaproughness"), true);
            let plastic = Arc::new(PlasticMaterial::new(kd, ks, roughness, remap_roughness));
            return Some(plastic);
        } else if api_state.graphics_state.material == String::from("translucent") {
            println!("TODO: CreateTranslucentMaterial");
        } else if api_state.graphics_state.material == String::from("glass") {
            let kr = mp.get_spectrum_texture(String::from("Kr"), Spectrum::new(1.0 as Float));
            let kt = mp.get_spectrum_texture(String::from("Kt"), Spectrum::new(1.0 as Float));
            // let some_eta = mp.get_float_texture_or_null(String::from("eta"));
            // if let Some(eta) = some_eta {
            //     println!("some eta");
            // } else {
            let eta = mp.get_float_texture(String::from("index"), 1.5);
            // }
            // std::shared_ptr<Texture<Float>> roughu =
            //     mp.GetFloatTexture("uroughness", 0.f);
            let roughu = mp.get_float_texture(String::from("uroughness"), 0.0 as Float);
            // std::shared_ptr<Texture<Float>> roughv =
            //     mp.GetFloatTexture("vroughness", 0.f);
            let roughv = mp.get_float_texture(String::from("vroughness"), 0.0 as Float);
            // std::shared_ptr<Texture<Float>> bumpMap =
            //     mp.GetFloatTextureOrNull("bumpmap");
            let remap_roughness: bool = mp.find_bool(String::from("remaproughness"), true);
            let glass = Arc::new(GlassMaterial {
                kr: kr,
                kt: kt,
                u_roughness: roughu,
                v_roughness: roughv,
                index: eta,
                remap_roughness: remap_roughness,
            });
            return Some(glass);
        } else if api_state.graphics_state.material == String::from("mirror") {
            let kr = mp.get_spectrum_texture(String::from("Kr"), Spectrum::new(0.9 as Float));
            // TODO: std::shared_ptr<Texture<Float>> bumpMap = mp.GetFloatTextureOrNull("bumpmap");
            let mirror = Arc::new(MirrorMaterial { kr: kr });
            return Some(mirror);
        } else if api_state.graphics_state.material == String::from("hair") {
            return Some(HairMaterial::create(&mut mp));
        } else if api_state.graphics_state.material == String::from("mix") {
            let m1: String = mp.find_string(String::from("namedmaterial1"), String::from(""));
            let m2: String = mp.find_string(String::from("namedmaterial2"), String::from(""));
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
                mp.get_spectrum_texture(String::from("amount"), Spectrum::new(0.5));
            if let Some(m1) = mat1 {
                if let Some(m2) = mat2 {
                    let mix = Arc::new(MixMaterial::new(m1.clone(), m2.clone(), scale));
                    return Some(mix);
                }
            }
            return None;
        } else if api_state.graphics_state.material == String::from("metal") {
            return Some(MetalMaterial::create(&mut mp));
        } else if api_state.graphics_state.material == String::from("substrate") {
            return Some(SubstrateMaterial::create(&mut mp));
        } else if api_state.graphics_state.material == String::from("uber") {
            return Some(UberMaterial::create(&mut mp));
        } else if api_state.graphics_state.material == String::from("subsurface") {
            println!("TODO: CreateSubsurfaceMaterial");
        } else if api_state.graphics_state.material == String::from("kdsubsurface") {
            println!("TODO: CreateKdsubsurfaceMaterial");
        } else if api_state.graphics_state.material == String::from("fourier") {
            println!("TODO: CreateFourierMaterial");
        } else {
            panic!(
                "Material \"{}\" unknown.",
                api_state.graphics_state.material
            );
        }
    }
    let kd = Arc::new(ConstantTexture::new(Spectrum::new(0.5)));
    let sigma = Arc::new(ConstantTexture::new(0.0 as Float));
    Some(Arc::new(MatteMaterial::new(kd, sigma)))
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

pub fn pbrt_init() -> ApiState {
    ApiState::default()
}

pub fn pbrt_cleanup() {
    println!("WorldEnd");
}

pub fn pbrt_transform(api_state: &mut ApiState, tr: &Transform) {
    println!("{:?}", tr);
    if api_state.active_transform_bits & 1_u8 > 0_u8 {
        // 0x?1
        api_state.cur_transform.t[0] = api_state.cur_transform.t[0] * *tr;
    }
    if api_state.active_transform_bits & 2_u8 > 0_u8 {
        // 0x1?
        api_state.cur_transform.t[1] = api_state.cur_transform.t[1] * *tr;
    }
}

pub fn pbrt_scale(api_state: &mut ApiState, sx: Float, sy: Float, sz: Float) {
    println!("Scale {} {} {}", sx, sy, sz);
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
    println!(
        "LookAt {} {} {} {} {} {} {} {} {}",
        ex, ey, ez, lx, ly, lz, ux, uy, uz
    );
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

pub fn pbrt_film(api_state: &mut ApiState, params: ParamSet) {
    println!("Film \"{}\" ", params.name);
    print_params(&params);
    api_state.render_options.film_name = params.name.clone();
    api_state.param_set = params;
}

pub fn pbrt_sampler(api_state: &mut ApiState, params: ParamSet) {
    println!("Sampler \"{}\" ", params.name);
    print_params(&params);
    api_state.render_options.sampler_name = params.name.clone();
    api_state.param_set = params;
}

pub fn pbrt_integrator(api_state: &mut ApiState, params: ParamSet) {
    println!("Integrator \"{}\" ", params.name);
    print_params(&params);
    api_state.render_options.integrator_name = params.name.clone();
    api_state.param_set = params;
}

pub fn pbrt_camera(api_state: &mut ApiState, params: ParamSet) {
    println!("Camera \"{}\" ", params.name);
    print_params(&params);
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
}

pub fn pbrt_world_begin(api_state: &mut ApiState) {
    println!("WorldBegin");
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
    println!("AttributeBegin");
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
    println!("AttributeEnd");
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

pub fn pbrt_make_named_material(api_state: &mut ApiState, params: ParamSet) {
    println!("MakeNamedMaterial \"{}\" ", params.name);
    print_params(&params);
    api_state.param_set = params;
    let mat_type: String = api_state
        .param_set
        .find_one_string(String::from("type"), String::new());
    if mat_type == String::new() {
        panic!("No parameter string \"type\" found in MakeNamedMaterial");
    }
    api_state.graphics_state.material = mat_type.clone();
    api_state
        .graphics_state
        .material_params
        .copy_from(&api_state.param_set);
    api_state.graphics_state.current_material = String::new();
    let mtl: Option<Arc<Material + Send + Sync>> = create_material(&api_state);
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
    println!("NamedMaterial \"{}\" ", params.name);
    api_state.param_set = params;
}

pub fn pbrt_area_light_source(api_state: &mut ApiState, params: ParamSet) {
    println!("AreaLightSource \"{}\" ", params.name);
    print_params(&params);
    api_state.param_set = params;
}

pub fn pbrt_shape(api_state: &mut ApiState, params: ParamSet) {
    println!("Shape \"{}\" ", params.name);
    print_params(&params);
    api_state.param_set = params;
}
