// std
use std::collections::HashMap;
use std::sync::Arc;
// pbrt
use core::light::Light;
use core::material::Material;
use core::medium::Medium;
use core::paramset::{ParamSet, TextureParams};
use core::pbrt::{Float, Spectrum};
use core::primitive::Primitive;
use core::texture::Texture;
use core::transform::{Matrix4x4, Transform};
use materials::matte::MatteMaterial;

// see api.cpp

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

pub fn pbrt_transform(tr: Transform) {
    println!("DONE: {:?}", tr);
    // TODO
    // VERIFY_INITIALIZED("Transform");
    // FOR_ACTIVE_TRANSFORMS(
    //     curTransform[i] = Transform(Matrix4x4(
    //         tr[0], tr[4], tr[8], tr[12], tr[1], tr[5], tr[9], tr[13], tr[2],
    //         tr[6], tr[10], tr[14], tr[3], tr[7], tr[11], tr[15]));)
    // if (PbrtOptions.cat || PbrtOptions.toPly) {
    //     printf("%*sTransform [ ", catIndentCount, "");
    //     for (int i = 0; i < 16; ++i) printf("%.9g ", tr[i]);
    //     printf(" ]\n");
    // }
}

pub fn pbrt_scale(sx: Float, sy: Float, sz: Float) {
    println!("DONE: Scale {:?} {:?} {:?}", sx, sy, sz);
    // TODO
    //     VERIFY_INITIALIZED("Scale");
    //     FOR_ACTIVE_TRANSFORMS(curTransform[i] =
    //                               curTransform[i] * Scale(sx, sy, sz);)
    //     if (PbrtOptions.cat || PbrtOptions.toPly)
    //         printf("%*sScale %.9g %.9g %.9g\n", catIndentCount, "", sx, sy, sz);
}

pub fn pbrt_look_at(
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
        "DONE: LookAt {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?}",
        ex, ey, ez, lx, ly, lz, ux, uy, uz
    );
    // TODO
    // VERIFY_INITIALIZED("LookAt");
    // Transform lookAt =
    //     LookAt(Point3f(ex, ey, ez), Point3f(lx, ly, lz), Vector3f(ux, uy, uz));
    // FOR_ACTIVE_TRANSFORMS(curTransform[i] = curTransform[i] * lookAt;);
    // if (PbrtOptions.cat || PbrtOptions.toPly)
    //     printf(
    //         "%*sLookAt %.9g %.9g %.9g\n%*s%.9g %.9g %.9g\n"
    //         "%*s%.9g %.9g %.9g\n",
    //         catIndentCount, "", ex, ey, ez, catIndentCount + 8, "", lx, ly, lz,
    //         catIndentCount + 8, "", ux, uy, uz);
}

pub fn pbrt_world_begin() {
    println!("DONE: WorldBegin");
    // TODO
    // VERIFY_OPTIONS("WorldBegin");
    // currentApiState = APIState::WorldBlock;
    // for (int i = 0; i < MaxTransforms; ++i) curTransform[i] = Transform();
    // activeTransformBits = AllTransformsBits;
    // namedCoordinateSystems["world"] = curTransform;
    // if (PbrtOptions.cat || PbrtOptions.toPly)
    //     printf("\n\nWorldBegin\n\n");
}

pub fn pbrt_attribute_begin() {
    println!("DONE: AttributeBegin");
    // TODO
    // VERIFY_WORLD("AttributeBegin");
    // pushedGraphicsStates.push_back(graphicsState);
    // graphicsState.floatTexturesShared = graphicsState.spectrumTexturesShared =
    //     graphicsState.namedMaterialsShared = true;
    // pushedTransforms.push_back(curTransform);
    // pushedActiveTransformBits.push_back(activeTransformBits);
    // if (PbrtOptions.cat || PbrtOptions.toPly) {
    //     printf("\n%*sAttributeBegin\n", catIndentCount, "");
    //     catIndentCount += 4;
    // }
}

pub fn pbrt_attribute_end() {
    println!("DONE: AttributeEnd");
    // TODO
    // VERIFY_WORLD("AttributeEnd");
    // if (!pushedGraphicsStates.size()) {
    //     Error(
    //         "Unmatched pbrtAttributeEnd() encountered. "
    //         "Ignoring it.");
    //     return;
    // }
    // graphicsState = std::move(pushedGraphicsStates.back());
    // pushedGraphicsStates.pop_back();
    // curTransform = pushedTransforms.back();
    // pushedTransforms.pop_back();
    // activeTransformBits = pushedActiveTransformBits.back();
    // pushedActiveTransformBits.pop_back();
    // if (PbrtOptions.cat || PbrtOptions.toPly) {
    //     catIndentCount -= 4;
    //     printf("%*sAttributeEnd\n", catIndentCount, "");
    // }
}
