#![recursion_limit="2000"]

extern crate pest;
#[macro_use]
extern crate pest_derive;
extern crate getopts;
extern crate pbrt;

use pbrt::accelerators::bvh::{BVHAccel, SplitMethod};
use pbrt::cameras::perspective::PerspectiveCamera;
use pbrt::core::api::{GraphicsState, RenderOptions, TransformSet};
use pbrt::core::camera::Camera;
use pbrt::core::filter::Filter;
use pbrt::core::geometry::{Bounds2f, Bounds2i, Normal3f, Point2f, Point2i, Point3f, Vector3f};
use pbrt::core::integrator::SamplerIntegrator;
use pbrt::core::light::Light;
use pbrt::core::material::Material;
use pbrt::core::mipmap::ImageWrap;
use pbrt::core::paramset::{ParamSet, TextureParams};
use pbrt::core::primitive::{GeometricPrimitive, Primitive, TransformedPrimitive};
use pbrt::core::pbrt::{Float, Spectrum};
use pbrt::core::pbrt::clamp_t;
use pbrt::core::transform::{AnimatedTransform, Matrix4x4, Transform};
use pbrt::core::film::Film;
use pbrt::core::sampler::Sampler;
use pbrt::core::scene::Scene;
use pbrt::core::shape::Shape;
use pbrt::core::texture::{PlanarMapping2D, Texture, TextureMapping2D, UVMapping2D};
use pbrt::filters::boxfilter::BoxFilter;
use pbrt::filters::gaussian::GaussianFilter;
use pbrt::filters::triangle::TriangleFilter;
use pbrt::integrators::ao::AOIntegrator;
use pbrt::integrators::directlighting::{DirectLightingIntegrator, LightStrategy};
use pbrt::integrators::path::PathIntegrator;
use pbrt::materials::glass::GlassMaterial;
use pbrt::materials::hair::HairMaterial;
use pbrt::materials::matte::MatteMaterial;
use pbrt::materials::metal::MetalMaterial;
use pbrt::materials::mirror::MirrorMaterial;
use pbrt::materials::mixmat::MixMaterial;
use pbrt::materials::plastic::PlasticMaterial;
use pbrt::materials::substrate::SubstrateMaterial;
use pbrt::materials::uber::UberMaterial;
use pbrt::lights::diffuse::DiffuseAreaLight;
use pbrt::lights::distant::DistantLight;
use pbrt::lights::infinite::InfiniteAreaLight;
use pbrt::lights::point::PointLight;
use pbrt::samplers::halton::HaltonSampler;
use pbrt::samplers::random::RandomSampler;
use pbrt::samplers::sobol::SobolSampler;
use pbrt::samplers::zerotwosequence::ZeroTwoSequenceSampler;
use pbrt::shapes::curve::create_curve_shape;
use pbrt::shapes::cylinder::Cylinder;
use pbrt::shapes::disk::Disk;
use pbrt::shapes::plymesh::create_ply_mesh;
use pbrt::shapes::sphere::Sphere;
use pbrt::shapes::triangle::{Triangle, TriangleMesh};
use pbrt::textures::checkerboard::Checkerboard2DTexture;
use pbrt::textures::constant::ConstantTexture;
use pbrt::textures::imagemap::ImageTexture;

// parser
use pest::Parser;
// getopts
use getopts::Options;
// std
use std::str::FromStr;
use std::collections::{HashMap, LinkedList};
use std::env;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::sync::Arc;

pub const VERSION: &'static str = env!("CARGO_PKG_VERSION");

static mut NUMBER_OF_THREADS: u8 = 0_u8;
static mut SEARCH_DIRECTORY: Option<Box<PathBuf>> = None;
static mut CUR_TRANSFORM: TransformSet = TransformSet {
    t: [Transform {
        m: Matrix4x4 {
            m: [[1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0]],
        },
        m_inv: Matrix4x4 {
            m: [[1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0]],
        },
    }; 2],
};
static mut ACTIVE_TRANSFORM_BITS: u8 = 3_u8; // 0x11 for MaxTransforms = 2
static mut NAMED_COORDINATE_SYSTEMS: Option<Box<HashMap<&str, TransformSet>>> = None;
static mut RENDER_OPTIONS: Option<Box<RenderOptions>> = None;
static mut GRAPHICS_STATE: Option<Box<GraphicsState>> = None;
static mut PUSHED_GRAPHICS_STATES: Option<Box<Vec<GraphicsState>>> = None;
static mut PUSHED_TRANSFORMS: Option<Box<Vec<TransformSet>>> = None;
static mut PUSHED_ACTIVE_TRANSFORM_BITS: Option<Box<Vec<u8>>> = None;
// not used in original C++ code:
static mut PARAM_SET: Option<Box<ParamSet>> = None;

#[derive(Debug, PartialEq)]
pub enum IntNode {
    Ints(LinkedList<IntNode>),
    OneInt(i32),
}

#[derive(Debug, PartialEq)]
pub enum FloatNode {
    Floats(LinkedList<FloatNode>),
    OneFloat(Float),
}

#[derive(Parser)]
#[grammar = "../examples/pbrt.pest"]
struct PbrtParser;

fn print_usage(program: &str, opts: Options) {
    let brief = format!("Usage: {} [options]", program);
    print!("{}", opts.usage(&brief));
}

fn print_version(program: &str) {
    println!("{} {}", program, VERSION);
}

fn print_params(params: &ParamSet) {
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
            println!("  \"point {}\" [{} {} {}]",
                     p.name,
                     p.values[0].x,
                     p.values[0].y,
                     p.values[0].z);
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
            println!("  \"vector {}\" [{} {} {}]",
                     p.name,
                     p.values[0].x,
                     p.values[0].y,
                     p.values[0].z);
        }
    }
    for p in &params.normals {
        if p.n_values == 1_usize {
            println!("  \"normal {}\" [{} {} {}]",
                     p.name,
                     p.values[0].x,
                     p.values[0].y,
                     p.values[0].z);
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
            println!("  \"rgb {}\" [{} {} {}]",
                     p.name,
                     p.values[0].c[0],
                     p.values[0].c[1],
                     p.values[0].c[2]);
        }
    }
    for p in &params.strings {
        if p.n_values == 1_usize {
            println!("  \"string {}\" [\"{}\"]", p.name, p.values[0]);
        }
    }
    for p in &params.textures {
        if p.n_values == 1_usize {
            println!("  \"texture {}\" \"{}\"", p.name, p.values[0]);
        }
    }
}

fn create_material() -> Arc<Material + Send + Sync> {
    unsafe {
        if let Some(ref mut graphics_state) = GRAPHICS_STATE {
            // CreateMaterial
            let mut material_params = ParamSet::default();
            material_params.copy_from(&graphics_state.material_params);
            let mut mp: TextureParams = TextureParams {
                float_textures: graphics_state.float_textures.clone(),
                spectrum_textures: graphics_state.spectrum_textures.clone(),
                geom_params: ParamSet::default(),
                material_params: material_params,
            };
            if graphics_state.current_material != String::new() {
                match graphics_state
                          .named_materials
                          .get(graphics_state.current_material.as_str()) {
                    Some(named_material) => {
                        return named_material.clone();
                    }
                    None => {
                        println!("WARNING: Named material \"{}\" not defined. Using \"matte\".",
                                 graphics_state.current_material);
                    }
                }
            } else {
                // MakeMaterial
                assert_ne!(graphics_state.material, String::new());
                assert_ne!(graphics_state.material, String::from("none"));
                if graphics_state.material == String::from("matte") {
                    return MatteMaterial::create(&mut mp);
                } else if graphics_state.material == String::from("plastic") {
                    let kd = mp.get_spectrum_texture(String::from("Kd"),
                                                     Spectrum::new(0.25 as Float));
                    let ks = mp.get_spectrum_texture(String::from("Ks"),
                                                     Spectrum::new(0.25 as Float));
                    let roughness = mp.get_float_texture(String::from("roughness"), 0.1 as Float);
                    // TODO: std::shared_ptr<Texture<Float>> bumpMap = mp.GetFloatTextureOrNull("bumpmap");
                    let remap_roughness: bool = mp.find_bool(String::from("remaproughness"), true);
                    let plastic =
                        Arc::new(PlasticMaterial::new(kd, ks, roughness, remap_roughness));
                    return plastic;
                } else if graphics_state.material == String::from("translucent") {
                    println!("TODO: CreateTranslucentMaterial");
                } else if graphics_state.material == String::from("glass") {
                    let kr = mp.get_spectrum_texture(String::from("Kr"),
                                                     Spectrum::new(1.0 as Float));
                    let kt = mp.get_spectrum_texture(String::from("Kt"),
                                                     Spectrum::new(1.0 as Float));
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
                    return glass;
                } else if graphics_state.material == String::from("mirror") {
                    let kr = mp.get_spectrum_texture(String::from("Kr"),
                                                     Spectrum::new(0.9 as Float));
                    // TODO: std::shared_ptr<Texture<Float>> bumpMap = mp.GetFloatTextureOrNull("bumpmap");
                    let mirror = Arc::new(MirrorMaterial { kr: kr });
                    return mirror;
                } else if graphics_state.material == String::from("hair") {
                    return HairMaterial::create(&mut mp);
                } else if graphics_state.material == String::from("mix") {
                    let m1: String = mp.find_string(String::from("namedmaterial1"),
                                                    String::from(""));
                    let m2: String = mp.find_string(String::from("namedmaterial2"),
                                                    String::from(""));
                    let mat1 = match graphics_state.named_materials.get(&m1) {
                        Some(named_material) => {
                            named_material
                        },
                        None => {
                            panic!("Material \"{}\" unknown.", m1);
                        },
                    };
                    let mat2 = match graphics_state.named_materials.get(&m2) {
                        Some(named_material) => {
                            named_material
                        },
                        None => {
                            panic!("Material \"{}\" unknown.", m2);
                        },
                    };
                    let scale: Arc<Texture<Spectrum> + Send + Sync> =
                        mp.get_spectrum_texture(String::from("amount"), Spectrum::new(0.5));
                    let mix = Arc::new(MixMaterial::new(mat1.clone(), mat2.clone(), scale));
                    return mix;
                } else if graphics_state.material == String::from("metal") {
                    return MetalMaterial::create(&mut mp);
                } else if graphics_state.material == String::from("substrate") {
                    return SubstrateMaterial::create(&mut mp);
                } else if graphics_state.material == String::from("uber") {
                    return UberMaterial::create(&mut mp);
                } else if graphics_state.material == String::from("subsurface") {
                    println!("TODO: CreateSubsurfaceMaterial");
                } else if graphics_state.material == String::from("kdsubsurface") {
                    println!("TODO: CreateKdsubsurfaceMaterial");
                } else if graphics_state.material == String::from("fourier") {
                    println!("TODO: CreateFourierMaterial");
                } else {
                    panic!("Material \"{}\" unknown.", graphics_state.material);
                }
            }
        }
    }
    let kd = Arc::new(ConstantTexture::new(Spectrum::new(0.5)));
    let sigma = Arc::new(ConstantTexture::new(0.0 as Float));
    Arc::new(MatteMaterial::new(kd, sigma))
}

fn make_light(param_set: &ParamSet, ro: &mut Box<RenderOptions>) {
    // MakeLight (api.cpp:591)
    if param_set.name == String::from("point") {
        let i: Spectrum = param_set
            .find_one_spectrum(String::from("I"), Spectrum::new(1.0 as Float));
        // Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
        // Point3f P = paramSet.FindOnePoint3f("from", Point3f(0, 0, 0));
        // Transform l2w = Translate(Vector3f(P.x, P.y, P.z)) * light2world;
        // return std::make_shared<PointLight>(l2w, medium, I * sc);
        unsafe {
            let point_light = Arc::new(PointLight::new(&CUR_TRANSFORM.t[0], &i));
            ro.lights.push(point_light);
        }
    } else if param_set.name == String::from("spot") {
        println!("TODO: CreateSpotLight");
    } else if param_set.name == String::from("goniometric") {
        println!("TODO: CreateGoniometricLight");
    } else if param_set.name == String::from("projection") {
        println!("TODO: CreateProjectionLight");
    } else if param_set.name == String::from("distant") {
        // CreateDistantLight
        let l: Spectrum = param_set
            .find_one_spectrum(String::from("L"), Spectrum::new(1.0 as Float));
        let sc: Spectrum =
            param_set.find_one_spectrum(String::from("scale"), Spectrum::new(1.0 as Float));
        let from: Point3f = param_set.find_one_point3f(String::from("from"),
                                                       Point3f {
                                                           x: 0.0,
                                                           y: 0.0,
                                                           z: 0.0,
                                                       });
        let to: Point3f = param_set.find_one_point3f(String::from("to"),
                                                     Point3f {
                                                         x: 0.0,
                                                         y: 0.0,
                                                         z: 0.0,
                                                     });
        let dir: Vector3f = from - to;
        // return std::make_shared<DistantLight>(light2world, L * sc, dir);
        unsafe {
            let distant_light = Arc::new(DistantLight::new(&CUR_TRANSFORM.t[0], &(l * sc), &dir));
            ro.lights.push(distant_light);
        }
    } else if param_set.name == String::from("infinite") ||
        param_set.name == String::from("exinfinite")
    {
        let l: Spectrum = param_set
            .find_one_spectrum(String::from("L"), Spectrum::new(1.0 as Float));
        let sc: Spectrum =
            param_set.find_one_spectrum(String::from("scale"), Spectrum::new(1.0 as Float));
        let mut texmap: String = param_set
            .find_one_filename(String::from("mapname"), String::from(""));
        if texmap != String::from("") {
            unsafe {
                if let Some(ref search_directory) = SEARCH_DIRECTORY {
                    // texmap = AbsolutePath(ResolveFilename(texmap));
                    let mut path_buf: PathBuf = PathBuf::from("/");
                    path_buf.push(search_directory.as_ref());
                path_buf.push(texmap);
                    texmap = String::from(path_buf.to_str().unwrap());
                }
            }
        }
        let n_samples: i32 = param_set.find_one_int(String::from("nsamples"), 1 as i32);
        // TODO: if (PbrtOptions.quickRender) nSamples = std::max(1, nSamples / 4);

        // return std::make_shared<InfiniteAreaLight>(light2world, L * sc, nSamples, texmap);
        unsafe {
            let infinte_light =
                Arc::new(InfiniteAreaLight::new(&CUR_TRANSFORM.t[0], &(l * sc), n_samples, texmap));
            ro.lights.push(infinte_light);
        }
    } else {
        panic!("MakeLight: unknown name {}", param_set.name);
    }
}

fn make_texture(param_set: &ParamSet) {
    unsafe {
        // pbrtTexture (api.cpp:1049)
        if let Some(ref mut graphics_state) = GRAPHICS_STATE {
            let mut geom_params: ParamSet = ParamSet::default();
            let mut material_params: ParamSet = ParamSet::default();
            geom_params.copy_from(param_set);
            material_params.copy_from(param_set);
            let mut tp: TextureParams = TextureParams {
                float_textures: graphics_state.float_textures.clone(),
                spectrum_textures: graphics_state.spectrum_textures.clone(),
                geom_params: geom_params,
                material_params: material_params,
            };
            if param_set.tex_type == String::from("float") {
                println!("TODO: MakeFloatTexture");
            } else if param_set.tex_type == String::from("color") ||
                      param_set.tex_type == String::from("spectrum") {
                match graphics_state
                          .spectrum_textures
                          .get(param_set.name.as_str()) {
                    Some(_spectrum_texture) => {
                        println!("Texture \"{}\" being redefined", param_set.name);
                    }
                    None => {}
                }
                // MakeSpectrumTexture(texname, curTransform[0], tp);
                if param_set.tex_name == String::from("constant") {
                    println!("TODO: CreateConstantSpectrumTexture");
                } else if param_set.tex_name == String::from("scale") {
                    println!("TODO: CreateScaleSpectrumTexture");
                } else if param_set.tex_name == String::from("mix") {
                    println!("TODO: CreateMixSpectrumTexture");
                } else if param_set.tex_name == String::from("bilerp") {
                    println!("TODO: CreateBilerpSpectrumTexture");
                } else if param_set.tex_name == String::from("imagemap") {
                    // CreateImageSpectrumTexture
                    let mut map: Option<Box<TextureMapping2D + Send + Sync>> = None;
                    let mapping: String = tp.find_string(String::from("mapping"),
                                                         String::from("uv"));
                    if mapping == String::from("uv") {
                        let su: Float = tp.find_float(String::from("uscale"), 1.0);
                        let sv: Float = tp.find_float(String::from("vscale"), 1.0);
                        let du: Float = tp.find_float(String::from("udelta"), 0.0);
                        let dv: Float = tp.find_float(String::from("vdelta"), 0.0);
                        map = Some(Box::new(UVMapping2D {
                                                su: su,
                                                sv: sv,
                                                du: du,
                                                dv: dv,
                                            }));
                    } else if mapping == String::from("spherical") {
                        println!("TODO: SphericalMapping2D");
                    } else if mapping == String::from("cylindrical") {
                        println!("TODO: CylindricalMapping2D");
                    } else if mapping == String::from("planar") {
                        map = Some(Box::new(PlanarMapping2D {
                                                vs: tp.find_vector3f(String::from("v1"),
                                                                     Vector3f {
                                                                         x: 1.0,
                                                                         y: 0.0,
                                                                         z: 0.0,
                                                                     }),
                                                vt: tp.find_vector3f(String::from("v2"),
                                                                     Vector3f {
                                                                         x: 0.0,
                                                                         y: 1.0,
                                                                         z: 0.0,
                                                                     }),
                                                ds: tp.find_float(String::from("udelta"), 0.0),
                                                dt: tp.find_float(String::from("vdelta"), 0.0),
                                            }));
                    } else {
                        panic!("2D texture mapping \"{}\" unknown", mapping);
                    }
                    // initialize _ImageTexture_ parameters
                    let max_aniso: Float = tp.find_float(String::from("maxanisotropy"), 8.0);
                    let do_trilinear: bool = tp.find_bool(String::from("trilinear"), false);
                    let wrap: String = tp.find_string(String::from("wrap"), String::from("repeat"));
                    let mut wrap_mode: ImageWrap = ImageWrap::Repeat;
                    if wrap == String::from("black") {
                        wrap_mode = ImageWrap::Black;
                    } else if wrap == String::from("clamp") {
                        wrap_mode = ImageWrap::Clamp;
                    }
                    let scale: Float = tp.find_float(String::from("scale"), 1.0);
                    let mut filename: String =
                        tp.find_filename(String::from("filename"), String::new());
                    if let Some(ref search_directory) = SEARCH_DIRECTORY {
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
                    let gamma: bool = tp.find_bool(String::from("gamma"), true);

                    if let Some(mapping) = map {
                        let st = Arc::new(ImageTexture::new(mapping,
                                                            filename,
                                                            do_trilinear,
                                                            max_aniso,
                                                            wrap_mode,
                                                            scale,
                                                            gamma));
                        graphics_state
                            .spectrum_textures
                            .insert(param_set.name.clone(), st);
                    }
                } else if param_set.tex_name == String::from("uv") {
                    println!("TODO: CreateUVSpectrumTexture");
                } else if param_set.tex_name == String::from("checkerboard") {
                    // CreateCheckerboardSpectrumTexture
                    let dim: i32 = tp.find_int(String::from("dimension"), 2);
                    if dim != 2 && dim != 3 {
                        panic!("{} dimensional checkerboard texture not supported", dim);
                    }
                    let tex1: Arc<Texture<Spectrum> + Send + Sync> =
                        tp.get_spectrum_texture(String::from("tex1"), Spectrum::new(1.0));
                    let tex2: Arc<Texture<Spectrum> + Send + Sync> =
                        tp.get_spectrum_texture(String::from("tex2"), Spectrum::new(0.0));
                    if dim == 2 {
                        let mut map: Option<Box<TextureMapping2D + Send + Sync>> = None;
                        let mapping: String =
                            tp.find_string(String::from("mapping"), String::from("uv"));
                        if mapping == String::from("uv") {
                            let su: Float = tp.find_float(String::from("uscale"), 1.0);
                            let sv: Float = tp.find_float(String::from("vscale"), 1.0);
                            let du: Float = tp.find_float(String::from("udelta"), 0.0);
                            let dv: Float = tp.find_float(String::from("vdelta"), 0.0);
                            map = Some(Box::new(UVMapping2D {
                                su: su,
                                sv: sv,
                                du: du,
                                dv: dv,
                            }));
                        } else if mapping == String::from("spherical") {
                            println!("TODO: SphericalMapping2D");
                        } else if mapping == String::from("cylindrical") {
                            println!("TODO: CylindricalMapping2D");
                        } else if mapping == String::from("planar") {
                            map = Some(Box::new(PlanarMapping2D {
                                                    vs: tp.find_vector3f(String::from("v1"),
                                                                         Vector3f {
                                                                             x: 1.0,
                                                                             y: 0.0,
                                                                             z: 0.0,
                                                                         }),
                                                    vt: tp.find_vector3f(String::from("v2"),
                                                                         Vector3f {
                                                                             x: 0.0,
                                                                             y: 1.0,
                                                                             z: 0.0,
                                                                         }),
                                                    ds: tp.find_float(String::from("udelta"), 0.0),
                                                    dt: tp.find_float(String::from("vdelta"), 0.0),
                                                }));
                        } else {
                            panic!("2D texture mapping \"{}\" unknown", mapping);
                        }
                        // TODO: aamode
                        if let Some(mapping) = map {
                            let st = Arc::new(Checkerboard2DTexture::new(mapping, tex1, tex2));
                            graphics_state
                                .spectrum_textures
                                .insert(param_set.name.clone(), st);
                        }
                    } else {
                        // dim == 3
                        println!("TODO: TextureMapping3D");
                    }
                } else if param_set.tex_name == String::from("dots") {
                    println!("TODO: CreateDotsSpectrumTexture");
                } else if param_set.tex_name == String::from("fbm") {
                    println!("TODO: CreateFBmSpectrumTexture");
                } else if param_set.tex_name == String::from("wrinkled") {
                    println!("TODO: CreateWrinkledSpectrumTexture");
                } else if param_set.tex_name == String::from("marble") {
                    println!("TODO: CreateMarbleSpectrumTexture");
                } else if param_set.tex_name == String::from("windy") {
                    println!("TODO: CreateWindySpectrumTexture");
                } else {
                    println!("Spectrum texture \"{}\" unknown.", param_set.tex_name);
                }
            } else {
                panic!("Texture type \"{}\" unknown.", param_set.tex_type);
            }
        }
        // MakeFloatTexture(texname, curTransform[0], tp);
        // or
        // MakeSpectrumTexture(texname, curTransform[0], tp);
    }
}

fn pbrt_bool_parameter<R, I>(pairs: &mut pest::iterators::Pairs<R, I>) -> (String, bool)
    where I: pest::inputs::Input,
          R: pest::RuleType
{
    // single string with or without brackets
    let ident = pairs.next();
    let string: String = String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
    let option = pairs.next();
    let lbrack = option.clone().unwrap();
    let string2: String;
    if lbrack.as_str() == String::from("[") {
        // check for brackets
        let string = pairs.next();
        let pair = string.unwrap().clone();
        let ident = pair.into_inner().next();
        string2 = String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
    } else {
        // no brackets
        let string = option.clone();
        let pair = string.unwrap().clone();
        let ident = pair.into_inner().next();
        string2 = String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
    }
    // return boolean (instead of string)
    let b: bool;
    if string2 == String::from("true") {
        b = true;
    } else if string2 == String::from("false") {
        b = false
    } else {
        println!("WARNING: parameter {:?} not well defined, defaulting to false",
                 string);
        b = false
    }
    (string, b)
}

fn pbrt_float_parameter<R, I>(pairs: &mut pest::iterators::Pairs<R, I>) -> (String, Vec<Float>)
    where I: pest::inputs::Input,
          R: pest::RuleType
{
    let mut floats: Vec<Float> = Vec::new();
    // single float or several floats using brackets
    let ident = pairs.next();
    let string: String = String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
    let option = pairs.next();
    let lbrack = option.clone().unwrap();
    if lbrack.as_str() == String::from("[") {
        // check for brackets
        let mut number = pairs.next();
        while number.is_some() {
            let pair = number.unwrap().clone();
            if pair.as_str() == String::from("]") {
                // closing bracket found
                break;
            } else {
                let float: Float = f32::from_str(pair.into_span().as_str()).unwrap();
                floats.push(float);
            }
            number = pairs.next();
        }
    } else {
        // no brackets
        let mut number = option.clone();
        while number.is_some() {
            let pair = number.unwrap().clone();
            let float: Float = f32::from_str(pair.into_span().as_str()).unwrap();
            floats.push(float);
            number = pairs.next();
        }
    }
    (string, floats)
}

fn pbrt_integer_parameter<R, I>(pairs: &mut pest::iterators::Pairs<R, I>) -> (String, Vec<i32>)
    where I: pest::inputs::Input,
          R: pest::RuleType
{
    let mut integers: Vec<i32> = Vec::new();
    // single integer or several integers using brackets
    let ident = pairs.next();
    let string: String = String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
    let option = pairs.next();
    let lbrack = option.clone().unwrap();
    if lbrack.as_str() == String::from("[") {
        // check for brackets
        let mut number = pairs.next();
        while number.is_some() {
            let pair = number.unwrap().clone();
            if pair.as_str() == String::from("]") {
                // closing bracket found
                break;
            } else {
                let integer: i32 = i32::from_str(pair.into_span().as_str()).unwrap();
                integers.push(integer);
            }
            number = pairs.next();
        }
    } else {
        // no brackets
        let mut number = option.clone();
        while number.is_some() {
            let pair = number.unwrap().clone();
            let integer: i32 = i32::from_str(pair.into_span().as_str()).unwrap();
            integers.push(integer);
            number = pairs.next();
        }
    }
    (string, integers)
}

fn pbrt_string_parameter<R, I>(pairs: &mut pest::iterators::Pairs<R, I>) -> (String, String)
    where I: pest::inputs::Input,
          R: pest::RuleType
{
    // single string with or without brackets
    let ident = pairs.next();
    let string1: String = String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
    let option = pairs.next();
    let lbrack = option.clone().unwrap();
    let string2: String;
    if lbrack.as_str() == String::from("[") {
        // check for brackets
        let string = pairs.next();
        let pair = string.unwrap().clone();
        let ident = pair.into_inner().next();
        string2 = String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
    } else {
        // no brackets
        let string = option.clone();
        let pair = string.unwrap().clone();
        let ident = pair.into_inner().next();
        string2 = String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
    }
    (string1, string2)
}

fn pbrt_texture_parameter<R, I>(pairs: &mut pest::iterators::Pairs<R, I>) -> (String, String)
    where I: pest::inputs::Input,
          R: pest::RuleType
{
    // single string with or without brackets
    let ident = pairs.next();
    let string1: String = String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
    let option = pairs.next();
    let lbrack = option.clone().unwrap();
    let string2: String;
    if lbrack.as_str() == String::from("[") {
        // check for brackets
        let string = pairs.next();
        let pair = string.unwrap().clone();
        let ident = pair.into_inner().next();
        string2 = String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
    } else {
        // no brackets
        let string = option.clone();
        let pair = string.unwrap().clone();
        let ident = pair.into_inner().next();
        string2 = String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
    }
    (string1, string2)
}

fn pbrt_shape(param_set: &ParamSet)
              -> (Vec<Arc<Shape + Send + Sync>>, Vec<Arc<Material + Send + Sync>>) {
    let mut shapes: Vec<Arc<Shape + Send + Sync>> = Vec::new();
    let mut materials: Vec<Arc<Material + Send + Sync>> = Vec::new();
    // pbrtShape (api.cpp:1153)
    // TODO: if (!curTransform.IsAnimated()) { ... }
    // TODO: transformCache.Lookup(curTransform[0], &ObjToWorld, &WorldToObj);
    unsafe {
        let mut obj_to_world: Transform = Transform {
            m: CUR_TRANSFORM.t[0].m,
            m_inv: CUR_TRANSFORM.t[0].m_inv,
        };
        let mut world_to_obj: Transform = Transform {
            m: CUR_TRANSFORM.t[0].m_inv,
            m_inv: CUR_TRANSFORM.t[0].m,
        };
        if CUR_TRANSFORM.is_animated() {
            if let Some(ref mut graphics_state) = GRAPHICS_STATE {
                if graphics_state.area_light != String::from("") {
                    println!("WARNING: Ignoring currently set area light when creating animated shape",);
                }
            }
            // WORK
            // set both transforms to identity
            obj_to_world = Transform::default();
            world_to_obj = Transform::default();
        }
        // MakeShapes (api.cpp:296)
        if param_set.name == String::from("sphere") {
            // CreateSphereShape
            let radius: Float = param_set.find_one_float(String::from("radius"), 1.0 as Float);
            let z_min: Float = param_set.find_one_float(String::from("zmin"), -radius);
            let z_max: Float = param_set.find_one_float(String::from("zmax"), radius);
            let phi_max: Float = param_set.find_one_float(String::from("phimax"), 360.0 as Float);
            let sphere = Arc::new(Sphere::new(obj_to_world,
                                              world_to_obj,
                                              false,
                                              false,
                                              radius,
                                              z_min,
                                              z_max,
                                              phi_max));
            let mtl: Arc<Material + Send + Sync> = create_material();
            shapes.push(sphere.clone());
            materials.push(mtl.clone());
        } else if param_set.name == String::from("cylinder") {
            let radius: Float = param_set.find_one_float(String::from("radius"), 1.0);
            let z_min: Float = param_set.find_one_float(String::from("zmin"), -radius);
            let z_max: Float = param_set.find_one_float(String::from("zmax"), radius);
            let phi_max: Float = param_set.find_one_float(String::from("phimax"), 360.0 as Float);
            let cylinder = Arc::new(Cylinder::new(obj_to_world,
                                                  world_to_obj,
                                                  false,
                                                  radius,
                                                  z_min,
                                                  z_max,
                                                  phi_max));
            let mtl: Arc<Material + Send + Sync> = create_material();
            shapes.push(cylinder.clone());
            materials.push(mtl.clone());
        } else if param_set.name == String::from("disk") {
            let height: Float = param_set.find_one_float(String::from("height"), 0.0);
            let radius: Float = param_set.find_one_float(String::from("radius"), 1.0);
            let inner_radius: Float = param_set.find_one_float(String::from("innerradius"), 0.0);
            let phi_max: Float = param_set.find_one_float(String::from("phimax"), 360.0);
            let disk = Arc::new(Disk::new(obj_to_world,
                                          world_to_obj,
                                          false,
                                          false,
                                          height,
                                          radius,
                                          inner_radius,
                                          phi_max));
            let mtl: Arc<Material + Send + Sync> = create_material();
            shapes.push(disk.clone());
            materials.push(mtl.clone());
        } else if param_set.name == String::from("cone") {
            println!("TODO: CreateConeShape");
        } else if param_set.name == String::from("paraboloid") {
            println!("TODO: CreateParaboloidShape");
        } else if param_set.name == String::from("hyperboloid") {
            println!("TODO: CreateHyperboloidShape");
        } else if param_set.name == String::from("curve") {
            let mtl: Arc<Material + Send + Sync> = create_material();
            let ply_shapes: Vec<Arc<Shape + Send + Sync>> =
                create_curve_shape(obj_to_world,
                                   world_to_obj,
                                   false, // reverse_orientation
                                   param_set);
            for shape in ply_shapes {
                shapes.push(shape.clone());
                materials.push(mtl.clone());
            }
        } else if param_set.name == String::from("trianglemesh") {
            let vi = param_set.find_int(String::from("indices"));
            let p = param_set.find_point3f(String::from("P"));
            // try "uv" with Point2f
            let mut uvs = param_set.find_point2f(String::from("uv"));
            if uvs.is_empty() {
                // try "st" with Point2f
                uvs = param_set.find_point2f(String::from("st"));
            }
            if uvs.is_empty() {
                // try "uv" with float
                let mut fuv = param_set.find_float(String::from("uv"));
                if fuv.is_empty() {
                    // try "st" with float
                    fuv = param_set.find_float(String::from("st"));
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
            let s = param_set.find_vector3f(String::from("S"));
            let mut s_ws: Vec<Vector3f> = Vec::new();
            if !s.is_empty() {
                assert!(s.len() == p.len());
                // transform tangents to world space
                let n_tangents: usize = s.len();
                for i in 0..n_tangents {
                    s_ws.push(obj_to_world.transform_vector(s[i]));
                }
            }
            let n = param_set.find_normal3f(String::from("N"));
            let mut n_ws: Vec<Normal3f> = Vec::new();
            if !n.is_empty() {
                assert!(n.len() == p.len());
                // transform normals to world space
                let n_normals: usize = n.len();
                for i in 0..n_normals {
                    n_ws.push(obj_to_world.transform_normal(n[i]));
                }
            }
            for i in 0..vi.len() {
                if vi[i] as usize >= p.len() {
                    panic!("trianglemesh has out of-bounds vertex index {} ({} \"P\" values were given)",
                           vi[i],
                           p.len());
                }
            }
            // TODO: alpha
            // CreateTriangleMesh
            // transform mesh vertices to world space
            let mut p_ws: Vec<Point3f> = Vec::new();
            let n_vertices: usize = p.len();
            for i in 0..n_vertices {
                p_ws.push(obj_to_world.transform_point(p[i]));
            }
            // vertex indices are expected as usize, not i32
            let mut vertex_indices: Vec<usize> = Vec::new();
            for i in 0..vi.len() {
                vertex_indices.push(vi[i] as usize);
            }
            if let Some(ref mut graphics_state) = GRAPHICS_STATE {
                let mesh = Arc::new(TriangleMesh::new(obj_to_world,
                                                      world_to_obj,
                                                      graphics_state.reverse_orientation,
                                                      false, // transform_swaps_handedness
                                                      vi.len() / 3, // n_triangles
                                                      vertex_indices,
                                                      n_vertices,
                                                      p_ws, // in world space
                                                      s_ws, // in world space
                                                      n_ws, // in world space
                                                      uvs));
                let mtl: Arc<Material + Send + Sync> = create_material();
                for id in 0..mesh.n_triangles {
                    let triangle = Arc::new(Triangle::new(mesh.object_to_world,
                                                          mesh.world_to_object,
                                                          mesh.reverse_orientation,
                                                          mesh.clone(),
                                                          id));
                    shapes.push(triangle.clone());
                    materials.push(mtl.clone());
                }
            }
        } else if param_set.name == String::from("plymesh") {
            if let Some(ref mut graphics_state) = GRAPHICS_STATE {
                if let Some(ref search_directory) = SEARCH_DIRECTORY {
                    let mtl: Arc<Material + Send + Sync> = create_material();
                    let ply_shapes: Vec<Arc<Shape + Send + Sync>> =
                        create_ply_mesh(obj_to_world,
                                        world_to_obj,
                                        false, // reverse_orientation
                                        param_set,
                                        graphics_state.float_textures.clone(),
                                        // additional parameters:
                                        Some(search_directory));
                    for shape in ply_shapes {
                        shapes.push(shape.clone());
                        materials.push(mtl.clone());
                    }
                }
            }
        } else if param_set.name == String::from("heightfield") {
            println!("TODO: CreateHeightfield");
        } else if param_set.name == String::from("loopsubdiv") {
            println!("TODO: CreateLoopSubdiv");
        } else if param_set.name == String::from("nurbs") {
            println!("TODO: CreateNURBS");
        } else {
            panic!("Shape \"{}\" unknown.", param_set.name);
        }
    }
    (shapes, materials)
}

fn pbrt_world_end() {
    // println!("WorldEnd");
    unsafe {
        if let Some(ref mut pushed_graphics_states) = PUSHED_GRAPHICS_STATES {
            assert!(pushed_graphics_states.len() == 0_usize,
                    "Missing end to pbrtAttributeBegin()");
            if let Some(ref mut pt) = PUSHED_TRANSFORMS {
                assert!(pt.len() == 0_usize, "Missing end to pbrtTransformBegin()");
                if let Some(ref mut ro) = RENDER_OPTIONS {
                    // MakeFilter
                    let mut some_filter: Option<Arc<Filter + Sync + Send>> = None;
                    if ro.filter_name == String::from("box") {
                        some_filter = Some(BoxFilter::create(&mut ro.filter_params));
                    } else if ro.filter_name == String::from("gaussian") {
                        some_filter = Some(GaussianFilter::create(&mut ro.filter_params));
                    } else if ro.filter_name == String::from("mitchell") {
                        println!("TODO: CreateMitchellFilter");
                    } else if ro.filter_name == String::from("sinc") {
                        println!("TODO: CreateSincFilter");
                    } else if ro.filter_name == String::from("triangle") {
                        some_filter = Some(TriangleFilter::create(&mut ro.filter_params));
                    } else {
                        panic!("Filter \"{}\" unknown.", ro.filter_name);
                    }
                    // MakeFilm
                    if ro.film_name == String::from("image") {
                        let filename: String =
                            ro.film_params
                                .find_one_string(String::from("filename"), String::new());
                        let xres: i32 = ro.film_params
                            .find_one_int(String::from("xresolution"), 1280);
                        let yres: i32 = ro.film_params
                            .find_one_int(String::from("yresolution"), 720);
                        // TODO: if (PbrtOptions.quickRender) xres = std::max(1, xres / 4);
                        // TODO: if (PbrtOptions.quickRender) yres = std::max(1, yres / 4);
                        let mut crop: Bounds2f = Bounds2f {
                            p_min: Point2f { x: 0.0, y: 0.0 },
                            p_max: Point2f { x: 1.0, y: 1.0 },
                        };
                        // TODO: const Float *cr = params.FindFloat("cropwindow", &cwi);
                        let cr: Vec<Float> = ro.film_params.find_float(String::from("cropwindow"));
                        if cr.len() == 4 {
                            crop.p_min.x = clamp_t(cr[0].min(cr[1]), 0.0, 1.0);
                            crop.p_max.x = clamp_t(cr[0].max(cr[1]), 0.0, 1.0);
                            crop.p_min.y = clamp_t(cr[2].min(cr[3]), 0.0, 1.0);
                            crop.p_max.y = clamp_t(cr[2].max(cr[3]), 0.0, 1.0);
                        } else {
                            panic!("{:?} values supplied for \"cropwindow\". Expected 4.", cr.len());
                        }
                        let scale: Float = ro.film_params.find_one_float(String::from("scale"),
                                                                         1.0);
                        let diagonal: Float = ro.film_params
                            .find_one_float(String::from("diagonal"), 35.0);
                        let max_sample_luminance: Float =
                            ro.film_params
                                .find_one_float(String::from("maxsampleluminance"),
                                                std::f32::INFINITY);
                        if let Some(filter) = some_filter {
                            let film: Arc<Film> = Arc::new(Film::new(Point2i { x: xres, y: yres },
                                                                     crop,
                                                                     filter,
                                                                     diagonal,
                                                                     filename,
                                                                     scale,
                                                                     max_sample_luminance));
                            // MakeCamera
                            // TODO: let mut some_camera: Option<Arc<Camera + Sync + Send>> = None;
                            let mut some_camera: Option<Box<Camera + Sync + Send>> = None;
                            // TODO: MediumInterface mediumInterface = graphicsState.CreateMediumInterface();
                            let animated_cam_to_world: AnimatedTransform =
                                AnimatedTransform::new(&ro.camera_to_world.t[0],
                                                       ro.transform_start_time,
                                                       &ro.camera_to_world.t[1],
                                                       ro.transform_end_time);
                            if ro.camera_name == String::from("perspective") {
                                let camera: Box<Camera + Send + Sync> =
                                    PerspectiveCamera::create(&ro.camera_params,
                                                              animated_cam_to_world,
                                    film);
                                some_camera = Some(camera);
                            } else if ro.camera_name == String::from("orthographic") {
                                println!("TODO: CreateOrthographicCamera");
                            } else if ro.camera_name == String::from("realistic") {
                                println!("TODO: CreateRealisticCamera");
                            } else if ro.camera_name == String::from("environment") {
                                println!("TODO: CreateEnvironmentCamera");
                            } else {
                                panic!("Camera \"{}\" unknown.", ro.camera_name);
                            }
                            if let Some(camera) = some_camera {
                                // MakeSampler
                                let mut some_sampler: Option<Box<Sampler + Sync + Send>> = None;
                                if ro.sampler_name == String::from("lowdiscrepancy") ||
                                   ro.sampler_name == String::from("02sequence") {
                                    let nsamp: i32 =
                                        ro.sampler_params
                                            .find_one_int(String::from("pixelsamples"), 16);
                                    let sd: i32 = ro.sampler_params
                                        .find_one_int(String::from("dimensions"), 4);
                                    // TODO: if (PbrtOptions.quickRender) nsamp = 1;
                                    let sampler = Box::new(ZeroTwoSequenceSampler::new(nsamp as
                                                                                       i64,
                                                                                       sd as i64));
                                    some_sampler = Some(sampler);
                                } else if ro.sampler_name == String::from("maxmindist") {
                                    println!("TODO: CreateMaxMinDistSampler");
                                } else if ro.sampler_name == String::from("halton") {
                                    let nsamp: i32 =
                                        ro.sampler_params
                                            .find_one_int(String::from("pixelsamples"), 16);
                                    // TODO: if (PbrtOptions.quickRender) nsamp = 1;
                                    let sample_at_center: bool =
                                        ro.integrator_params
                                            .find_one_bool(String::from("samplepixelcenter"),
                                                           false);
                                    let sample_bounds: Bounds2i =
                                        camera.get_film().get_sample_bounds();
                                    let sampler = Box::new(HaltonSampler::new(nsamp as i64,
                                                                              sample_bounds,
                                                                              sample_at_center));
                                    some_sampler = Some(sampler);
                                } else if ro.sampler_name == String::from("sobol") {
                                    let nsamp: i32 =
                                        ro.sampler_params
                                            .find_one_int(String::from("pixelsamples"), 16);
                                    let sample_bounds: Bounds2i =
                                        camera.get_film().get_sample_bounds();
                                    let sampler = Box::new(SobolSampler::new(nsamp as i64,
                                                                             sample_bounds));
                                    some_sampler = Some(sampler);
                                } else if ro.sampler_name == String::from("random") {
                                    let nsamp: i32 =
                                        ro.sampler_params
                                            .find_one_int(String::from("pixelsamples"), 4);
                                    let sampler = Box::new(RandomSampler::new(nsamp as i64));
                                    some_sampler = Some(sampler);
                                } else if ro.sampler_name == String::from("stratified") {
                                    println!("TODO: CreateStratifiedSampler");
                                } else {
                                    panic!("Sampler \"{}\" unknown.", ro.sampler_name);
                                }
                                if let Some(mut sampler) = some_sampler {
                                    // MakeIntegrator
                                    // if let Some(mut sampler) = some_sampler {
                                    let mut some_integrator: Option<Box<SamplerIntegrator + Sync + Send>> = None;
                                    if ro.integrator_name == String::from("whitted") {
                                        println!("TODO: CreateWhittedIntegrator");
                                    } else if ro.integrator_name == String::from("directlighting") {
                                        // CreateDirectLightingIntegrator
                                        let max_depth: i32 =
                                            ro.integrator_params
                                                .find_one_int(String::from("maxdepth"), 5);
                                        let st: String = ro.integrator_params
                                            .find_one_string(String::from("strategy"),
                                                             String::from("all"));
                                        let strategy: LightStrategy;
                                        if st == String::from("one") {
                                            strategy = LightStrategy::UniformSampleOne;
                                        } else if st == String::from("all") {
                                            strategy = LightStrategy::UniformSampleAll;
                                        } else {
                                            panic!("Strategy \"{}\" for direct lighting unknown.",
                                                   st);
                                        }
                                        // TODO: const int *pb = params.FindInt("pixelbounds", &np);
                                        let pixel_bounds: Bounds2i = Bounds2i {
                                            p_min: Point2i { x: 0, y: 0 },
                                            p_max: Point2i { x: xres, y: yres },
                                        };
                                        let integrator =
                                            Box::new(DirectLightingIntegrator::new(strategy,
                                                                                   max_depth as
                                                                                   i64,
                                                                                   pixel_bounds));
                                        some_integrator = Some(integrator);
                                    } else if ro.integrator_name == String::from("path") {
                                        // CreatePathIntegrator
                                        let max_depth: i32 =
                                            ro.integrator_params
                                                .find_one_int(String::from("maxdepth"), 5);
                                        let pb: Vec<i32> = ro.integrator_params
                                            .find_int(String::from("pixelbounds"));
                                        let np: usize = pb.len();
                                        let pixel_bounds: Bounds2i = camera.get_film().get_sample_bounds();
                                        if np > 0 as usize {
                                            if np != 4 as usize {
                                                panic!("Expected four values for \"pixelbounds\" parameter. Got {}.",
                                                       np);
                                            } else {
                                                println!("TODO: pixelBounds = Intersect(...)");
                                                // pixelBounds = Intersect(pixelBounds,
                                                //                         Bounds2i{{pb[0], pb[2]}, {pb[1], pb[3]}});
                                                // if (pixelBounds.Area() == 0)
                                                //     Error("Degenerate \"pixelbounds\" specified.");
                                            }
                                        }
                                        let rr_threshold: Float =
                                            ro.integrator_params
                                            .find_one_float(String::from("rrthreshold"),
                                                            1.0 as Float);
                                        println!("DEBUG: rr_threshold = {:?}", rr_threshold);
                                        // std::string lightStrategy =
                                        //     params.FindOneString("lightsamplestrategy", "spatial");
                                        let light_strategy: String =
                                            ro.integrator_params
                                            .find_one_string(String::from("lightsamplestrategy"),
                                                             String::from("spatial"));
                                        let integrator =
                                            Box::new(PathIntegrator::new(max_depth as u32,
                                                                         pixel_bounds,
                                                                         rr_threshold,
                                                                         light_strategy));
                                        some_integrator = Some(integrator);
                                    } else if ro.integrator_name == String::from("volpath") {
                                        println!("TODO: CreateVolPathIntegrator");
                                    } else if ro.integrator_name == String::from("bdpt") {
                                        println!("TODO: CreateBDPTIntegrator");
                                    } else if ro.integrator_name == String::from("mlt") {
                                        println!("TODO: CreateMLTIntegrator");
                                    } else if ro.integrator_name ==
                                              String::from("ambientocclusion") {
                                        // CreateAOIntegrator
                                        let pb: Vec<i32> = ro.integrator_params
                                            .find_int(String::from("pixelbounds"));
                                        let np: usize = pb.len();
                                        let pixel_bounds: Bounds2i = camera.get_film().get_sample_bounds();
                                        if np > 0 as usize {
                                            if np != 4 as usize {
                                                panic!("Expected four values for \"pixelbounds\" parameter. Got {}.",
                                                       np);
                                            } else {
                                                println!("TODO: pixelBounds = Intersect(...)");
                                                // pixelBounds = Intersect(pixelBounds,
                                                //                         Bounds2i{{pb[0], pb[2]}, {pb[1], pb[3]}});
                                                // if (pixelBounds.Area() == 0)
                                                //     Error("Degenerate \"pixelbounds\" specified.");
                                            }
                                        }
                                        let rr_threshold: Float =
                                            ro.integrator_params
                                            .find_one_float(String::from("rrthreshold"),
                                                            1.0 as Float);
                                        println!("DEBUG: rr_threshold = {:?}", rr_threshold);
                                        let cos_sample: bool =
                                            ro.integrator_params
                                            .find_one_bool(String::from("cossample"), true);
                                        println!("DEBUG: cos_sample = {:?}", cos_sample);
                                        // int nSamples = params.Find_One_Int("nsamples", 64);
                                        let n_samples: i32 =
                                            ro.integrator_params
                                                .find_one_int(String::from("nsamples"), 64 as i32);
                                        // return new AOIntegrator(cosSample, nSamples, camera, sampler, pixelBounds);

                                        let integrator = Box::new(AOIntegrator::new(cos_sample,
                                                                                    n_samples,
                                                                                    pixel_bounds));
                                        some_integrator = Some(integrator);
                                    } else if ro.integrator_name == String::from("sppm") {
                                        println!("TODO: CreateSPPMIntegrator");
                                    } else {
                                        panic!("Integrator \"{}\" unknown.", ro.integrator_name);
                                    }
                                    if let Some(mut integrator) = some_integrator {
                                        // MakeIntegrator
                                        // TODO: if (renderOptions->haveScatteringMedia && ...)
                                        if ro.lights.is_empty() {
                                            // warn if no light sources are defined
                                            println!("WARNING: No light sources defined in scene; rendering a black image.",);
                                        }
                                        // MakeAccelerator
                                        if ro.accelerator_name == String::from("bvh") {
                                            //  CreateBVHAccelerator
                                            let split_method_name: String =
                                                ro.accelerator_params.find_one_string(String::from("splitmethod"),
                                                                                      String::from("sah"));
                                            let split_method;
                                            if split_method_name == String::from("sah") {
                                                split_method = SplitMethod::SAH;
                                            } else if split_method_name == String::from("hlbvh") {
                                                split_method = SplitMethod::HLBVH;
                                            } else if split_method_name == String::from("middle") {
                                                split_method = SplitMethod::Middle;
                                            } else if split_method_name == String::from("equal") {
                                                split_method = SplitMethod::EqualCounts;
                                            } else {
                                                println!("WARNING: BVH split method \"{}\" unknown.  Using \"sah\".",
                                                         split_method_name);
                                                split_method = SplitMethod::SAH;
                                            }
                                            let max_prims_in_node: i32 =
                                                ro.accelerator_params.find_one_int(String::from("maxnodeprims"), 4);
                                            let accelerator =
                                                Arc::new(BVHAccel::new(ro.primitives.clone(),
                                                                       max_prims_in_node as usize,
                                                                       split_method));
                                            // MakeScene
                                            let scene: Scene = Scene::new(accelerator.clone(),
                                                                          ro.lights.clone());
                                            // TODO: primitives.erase(primitives.begin(), primitives.end());
                                            // TODO: lights.erase(lights.begin(), lights.end());
                                            let num_threads: u8 = NUMBER_OF_THREADS;
                                            pbrt::render(&scene,
                                                         camera,
                                                         &mut sampler,
                                                         &mut integrator,
                                                         num_threads);
                                        } else if ro.accelerator_name == String::from("kdtree") {
                                            // println!("TODO: CreateKdTreeAccelerator");
                                            // WARNING: Use BVHAccel for now !!!
                                            let accelerator =
                                                Arc::new(BVHAccel::new(ro.primitives.clone(),
                                                                       4,
                                                                       SplitMethod::SAH));
                                            // MakeScene
                                            let scene: Scene = Scene::new(accelerator.clone(),
                                                                          ro.lights.clone());
                                            // TODO: primitives.erase(primitives.begin(), primitives.end());
                                            // TODO: lights.erase(lights.begin(), lights.end());
                                            let num_threads: u8 = NUMBER_OF_THREADS;
                                            pbrt::render(&scene,
                                                         camera,
                                                         &mut sampler,
                                                         &mut integrator,
                                                         num_threads);
                                        } else {
                                            panic!("Accelerator \"{}\" unknown.",
                                                   ro.accelerator_name);
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
                        panic!("Film \"{}\" unknown.", ro.film_name);
                    }
                }
            }
        }
    }
}

fn main() {
    // handle command line options
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();
    let mut opts = Options::new();
    opts.optflag("h", "help", "print this help menu");
    opts.optopt("i", "", "parse an input file", "FILE");
    opts.optopt("t",
                "nthreads",
                "use specified number of threads for rendering",
                "NUM");
    opts.optflag("v", "version", "print version number");
    let matches = match opts.parse(&args[1..]) {
        Ok(m) => m,
        Err(f) => panic!(f.to_string()),
    };
    if matches.opt_present("h") {
        print_usage(&program, opts);
        return;
    } else if matches.opt_present("i") {
        if matches.opt_present("t") {
            let nthreads = matches.opt_str("t");
            match nthreads {
                Some(x) => {
                    let number_result = x.parse::<u8>();
                    assert!(!number_result.is_err(),
                            "ERROR: 8 bit unsigned integer expected");
                    let num_threads: u8 = number_result.unwrap();
                    println!("nthreads = {:?}", num_threads);
                    unsafe {
                        NUMBER_OF_THREADS = num_threads;
                    }
                }
                None => panic!("No argument for number of threads given."),
            }
        }
        let infile = matches.opt_str("i");
        match infile {
            Some(x) => {
                println!("FILE = {}", x);
                let f = File::open(x.clone()).unwrap();
                let ip: &Path = Path::new(x.as_str());
                if ip.is_relative() {
                    let cp: PathBuf = env::current_dir().unwrap();
                    let pb: PathBuf = cp.join(ip);
                    let search_directory: &Path = pb.as_path().parent().unwrap();
                    println!("search_directory is {}", search_directory.display());
                    unsafe {
                        SEARCH_DIRECTORY = Some(Box::new(PathBuf::from(search_directory)));
                    }
                }
                let mut reader = BufReader::new(f);
                let mut str_buf: String = String::default();
                let num_bytes = reader.read_to_string(&mut str_buf);
                if num_bytes.is_ok() {
                    let n_bytes = num_bytes.unwrap();
                    println!("{} bytes read", n_bytes);
                }
                unsafe {
                    // render options
                    NAMED_COORDINATE_SYSTEMS = Some(Box::new(HashMap::new()));
                    RENDER_OPTIONS = Some(Box::new(RenderOptions::default()));
                    GRAPHICS_STATE = Some(Box::new(GraphicsState::new()));
                    PUSHED_GRAPHICS_STATES = Some(Box::new(Vec::new()));
                    PUSHED_TRANSFORMS = Some(Box::new(Vec::new()));
                    PUSHED_ACTIVE_TRANSFORM_BITS = Some(Box::new(Vec::new()));
                    PARAM_SET = Some(Box::new(ParamSet::default()));
                    // parser
                    let pairs = PbrtParser::parse_str(Rule::pbrt, &str_buf)
                        .unwrap_or_else(|e| panic!("{}", e));
                    // assert!(parser.pbrt());
                    // assert!(parser.end());
                    // println!("{:?}", parser.queue());
                    println!("do something with created tokens ...");
                    for pair in pairs {
                        // a pair is a combination of the rule which matched and a span of input
                        match pair.as_rule() {
                            Rule::statement => {
                                for statement_pair in pair.into_inner() {
                                    match statement_pair.as_rule() {
                                        Rule::active_transform => {
                                            for active_transform_pair in
                                                statement_pair.into_inner() {
                                                match active_transform_pair.as_rule() {
                                                    Rule::all => {
                                                        // println!("ActiveTransform All");
                                                        ACTIVE_TRANSFORM_BITS = 3_u8 // 0x11
                                                    }
                                                    Rule::start_time => {
                                                        // println!("ActiveTransform StartTime");
                                                        ACTIVE_TRANSFORM_BITS = 1_u8 // 0x01
                                                    }
                                                    Rule::end_time => {
                                                        // println!("ActiveTransform EndTime");
                                                        ACTIVE_TRANSFORM_BITS = 2_u8 // 0x10
                                                    }
                                                    _ => unreachable!(),
                                                }
                                            }
                                        }
                                        Rule::concat_transform => {
                                            let mut numbers: Vec<Float> = Vec::new();
                                            for concat_transform_pair in
                                                statement_pair.into_inner() {
                                                // ignore brackets
                                                let not_opening: bool = concat_transform_pair.as_str() != String::from("[");
                                                let not_closing: bool = concat_transform_pair.as_str() != String::from("]");
                                                if not_opening && not_closing {
                                                    let number: Float =
                                                        f32::from_str(concat_transform_pair.clone().into_span().as_str()).unwrap();
                                                    numbers.push(number);
                                                }
                                            }
                                            let m00: Float = numbers[0];
                                            let m01: Float = numbers[1];
                                            let m02: Float = numbers[2];
                                            let m03: Float = numbers[3];
                                            let m10: Float = numbers[4];
                                            let m11: Float = numbers[5];
                                            let m12: Float = numbers[6];
                                            let m13: Float = numbers[7];
                                            let m20: Float = numbers[8];
                                            let m21: Float = numbers[9];
                                            let m22: Float = numbers[10];
                                            let m23: Float = numbers[11];
                                            let m30: Float = numbers[12];
                                            let m31: Float = numbers[13];
                                            let m32: Float = numbers[14];
                                            let m33: Float = numbers[15];
                                            let transform: Transform = Transform::new(m00, m10, m20, m30,
                                                                                      m01, m11, m21, m31,
                                                                                      m02, m12, m22, m32,
                                                                                      m03, m13, m23, m33);
                                            if ACTIVE_TRANSFORM_BITS & 1_u8 > 0_u8 {
                                                // 0x?1
                                                CUR_TRANSFORM.t[0] = CUR_TRANSFORM.t[0] * transform;
                                            }
                                            if ACTIVE_TRANSFORM_BITS & 2_u8 > 0_u8 {
                                                // 0x1?
                                                CUR_TRANSFORM.t[1] = CUR_TRANSFORM.t[1] * transform;
                                            }
                                        }
                                        Rule::keyword => {
                                            for keyword_pair in statement_pair.into_inner() {
                                                match keyword_pair.as_rule() {
                                                    Rule::attribute_begin => {
                                                        if let Some(ref mut graphics_state) =
                                                            GRAPHICS_STATE {
                                                            if let Some(ref mut pushed_graphics_states) = PUSHED_GRAPHICS_STATES {
                                                                let mut material_param_set: ParamSet = ParamSet::default();
                                                                material_param_set.copy_from(&graphics_state.material_params);
                                                                let mut area_light_param_set: ParamSet = ParamSet::default();
                                                                area_light_param_set.copy_from(&graphics_state.area_light_params);
                                                                pushed_graphics_states.push(GraphicsState {
                                                                    float_textures: graphics_state.float_textures.clone(),
                                                                    spectrum_textures: graphics_state.spectrum_textures.clone(),
                                                                    material_params: material_param_set,
                                                                    material: String::from(graphics_state.material.as_ref()),
                                                                    named_materials: graphics_state.named_materials.clone(),
                                                                    current_material: String::from(graphics_state.current_material.as_ref()),
                                                                    area_light_params: area_light_param_set,
                                                                    area_light: String::from(graphics_state.area_light.as_ref()),
                                                                    reverse_orientation: graphics_state.reverse_orientation,
                                                                });
                                                            }
                                                            if let Some(ref mut pt) =
                                                                PUSHED_TRANSFORMS {
                                                                pt.push(TransformSet {
                                                                    t: [
                                                                        Transform {
                                                                            m: CUR_TRANSFORM.t[0].m,
                                                                            m_inv: CUR_TRANSFORM.t[0].m_inv,},
                                                                        Transform {
                                                                            m: CUR_TRANSFORM.t[1].m,
                                                                            m_inv: CUR_TRANSFORM.t[1].m_inv,},
                                                                    ]
                                                                });
                                                            }
                                                            if let Some(ref mut patb) =
                                                                PUSHED_ACTIVE_TRANSFORM_BITS {
                                                                patb.push(ACTIVE_TRANSFORM_BITS);
                                                            }
                                                        }
                                                    }
                                                    Rule::attribute_end => {
                                                        if let Some(ref mut graphics_state) =
                                                            GRAPHICS_STATE {
                                                            if let Some(ref mut pushed_graphics_states) = PUSHED_GRAPHICS_STATES {
                                                                if !(pushed_graphics_states.len() >= 1_usize) {
                                                                    panic!("Unmatched pbrtAttributeEnd() encountered.")
                                                                }
                                                                let pgs: GraphicsState = pushed_graphics_states.pop().unwrap();
                                                                // material_params
                                                                graphics_state.material_params.reset(String::new(),
                                                                                                     String::from(""),
                                                                                                     String::from(""),
                                                                                                     String::new());
                                                                graphics_state.material_params.copy_from(&pgs.material_params);
                                                                // material
                                                                graphics_state.material = String::from(pgs.material.as_ref());
                                                                // area_light_params
                                                                graphics_state.area_light_params.reset(String::new(),
                                                                                                       String::from(""),
                                                                                                       String::from(""),
                                                                                                       String::new());
                                                                graphics_state.area_light_params.copy_from(&pgs.area_light_params);
                                                                // area_light
                                                                graphics_state.area_light = String::from(pgs.area_light.as_ref());
                                                                // reverse_orientation
                                                                graphics_state.reverse_orientation = pgs.reverse_orientation;
                                                            }
                                                            if let Some(ref mut pt) =
                                                                PUSHED_TRANSFORMS {
                                                                let popped_transform_set: TransformSet = pt.pop().unwrap();
                                                                CUR_TRANSFORM.t[0] =
                                                                    popped_transform_set.t[0];
                                                                CUR_TRANSFORM.t[1] =
                                                                    popped_transform_set.t[1];
                                                            }
                                                            if let Some(ref mut patb) =
                                                                PUSHED_ACTIVE_TRANSFORM_BITS {
                                                                let active_transform_bits: u8 = patb.pop().unwrap();
                                                                ACTIVE_TRANSFORM_BITS =
                                                                    active_transform_bits;
                                                            }
                                                        }
                                                    }
                                                    Rule::transform_begin => {
                                                        if let Some(ref mut pt) =
                                                            PUSHED_TRANSFORMS {
                                                                pt.push(TransformSet {
                                                                    t: [
                                                                        Transform {
                                                                            m: CUR_TRANSFORM.t[0].m,
                                                                            m_inv: CUR_TRANSFORM.t[0].m_inv,},
                                                                        Transform {
                                                                            m: CUR_TRANSFORM.t[1].m,
                                                                            m_inv: CUR_TRANSFORM.t[1].m_inv,},
                                                                    ]
                                                                });
                                                            }
                                                        if let Some(ref mut patb) =
                                                            PUSHED_ACTIVE_TRANSFORM_BITS {
                                                                patb.push(ACTIVE_TRANSFORM_BITS);
                                                            }
                                                    }
                                                    Rule::transform_end => {
                                                        if let Some(ref mut pt) =
                                                            PUSHED_TRANSFORMS {
                                                                let popped_transform_set: TransformSet = pt.pop().unwrap();
                                                                CUR_TRANSFORM.t[0] =
                                                                    popped_transform_set.t[0];
                                                                CUR_TRANSFORM.t[1] =
                                                                    popped_transform_set.t[1];
                                                            }
                                                        if let Some(ref mut patb) =
                                                            PUSHED_ACTIVE_TRANSFORM_BITS {
                                                                let active_transform_bits: u8 = patb.pop().unwrap();
                                                                ACTIVE_TRANSFORM_BITS =
                                                                    active_transform_bits;
                                                            }
                                                    }
                                                    Rule::reverse_orientation => {
                                                        println!("ReverseOrientation");
                                                        if let Some(ref mut graphics_state) = GRAPHICS_STATE {
                                                            graphics_state.reverse_orientation = !graphics_state.reverse_orientation;
                                                        }
                                                    }
                                                    Rule::world_begin => {
                                                        println!("WorldBegin");
                                                        CUR_TRANSFORM.t[0] = Transform::default();
                                                        CUR_TRANSFORM.t[1] = Transform::default();
                                                        ACTIVE_TRANSFORM_BITS = 3_u8; // 0x11
                                                        if let Some(ref mut named_coordinate_systems) = NAMED_COORDINATE_SYSTEMS {
                                                            named_coordinate_systems.insert("world",
                                                                                            TransformSet {
                                                                                                t: [Transform::default(); 2]
                                                                                            });
                                                        }
                                                    }
                                                    _ => unreachable!(),
                                                }
                                            }
                                        }
                                        Rule::look_at => {
                                            let mut numbers: Vec<Float> = Vec::new();
                                            for look_at_pair in statement_pair.into_inner() {
                                                let number: Float =
                                                    f32::from_str(look_at_pair.clone().into_span().as_str()).unwrap();
                                                numbers.push(number);
                                            }
                                            let pos: Point3f = Point3f { x: numbers[0], y: numbers[1], z: numbers[2], };
                                            let look: Point3f = Point3f { x: numbers[3], y: numbers[4], z: numbers[5], };
                                            let up: Vector3f = Vector3f { x: numbers[6], y: numbers[7], z: numbers[8], };
                                            let look_at: Transform = Transform::look_at(pos, look, up);
                                            if ACTIVE_TRANSFORM_BITS & 1_u8 > 0_u8 {
                                                // 0x?1
                                                CUR_TRANSFORM.t[0] = CUR_TRANSFORM.t[0] * look_at;
                                            }
                                            if ACTIVE_TRANSFORM_BITS & 2_u8 > 0_u8 {
                                                // 0x1?
                                                CUR_TRANSFORM.t[1] = CUR_TRANSFORM.t[1] * look_at;
                                            }
                                        }
                                        Rule::named_statement => {
                                            for named_statement_pair in
                                                statement_pair.into_inner() {
                                                match named_statement_pair.as_rule() {
                                                    Rule::accelerator => {
                                                        println!("TODO: Rule::accelerator")
                                                    }
                                                    Rule::area_light_source => {
                                                        for area_light_source_pair in
                                                            named_statement_pair.into_inner() {
                                                            match area_light_source_pair.as_rule() {
                                                                Rule::string => {
                                                                    let mut string_pairs =
                                                                        area_light_source_pair
                                                                            .into_inner();
                                                                    let ident = string_pairs.next();
                                                                    let name: String =
                                                                        String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
                                                                    if let Some(ref mut param_set) = PARAM_SET {
                                                                        param_set.reset(String::from("AreaLightSource"),
                                                                                        String::from(name),
                                                                                        String::from(""),
                                                                                        String::from(""));
                                                                    } else {
                                                                        panic!("Can't get parameter set.");
                                                                    }
                                                                }
                                                                Rule::parameter => {
                                                                    for parameter_pair in
                                                                        area_light_source_pair
                                                                            .into_inner() {
                                                                        match parameter_pair
                                                                                  .as_rule() {
                                                                            Rule::bool_param => {
                                                                                let tuple: (String, bool) = pbrt_bool_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let b: bool = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_bool(string, b);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::blackbody_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_blackbody_spectrum(string, floats);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::float_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 1 {
                                                                                        param_set.add_float(string, floats[0]);
                                                                                    } else {
                                                                                        param_set.add_floats(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::integer_param => {
                                                                                let tuple: (String, Vec<i32>) = pbrt_integer_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let integers: Vec<i32> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if integers.len() == 1 {
                                                                                        param_set.add_int(string, integers[0]);
                                                                                    } else {
                                                                                        param_set.add_ints(string, integers);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::point_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_point3f(string,
                                                                                                              Point3f {
                                                                                                                  x: floats[0],
                                                                                                                  y: floats[1],
                                                                                                                  z: floats[2],
                                                                                                              });
                                                                                    } else {
                                                                                        param_set.add_point3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::normal_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_normal3f(string,
                                                                                                               Normal3f {
                                                                                                                   x: floats[0],
                                                                                                                   y: floats[1],
                                                                                                                   z: floats[2],
                                                                                                               });
                                                                                    } else {
                                                                                        param_set.add_normal3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::rgb_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_rgb_spectrum(string,
                                                                                                               Spectrum {
                                                                                                                   c: [floats[0], floats[1], floats[2]],
                                                                                                               });
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::spectrum_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_rgb_spectrum(string,
                                                                                                               Spectrum {
                                                                                                                   c: [floats[0], floats[1], floats[2]],
                                                                                                               });
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::string_param => {
                                                                                let tuple: (String, String) = pbrt_string_parameter(&mut parameter_pair.into_inner());
                                                                                let string1: String = tuple.0;
                                                                                let string2: String = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_string(string1, string2);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::texture_param => {
                                                                                let tuple: (String, String) = pbrt_texture_parameter(&mut parameter_pair.into_inner());
                                                                                let string1: String = tuple.0;
                                                                                let string2: String = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_texture(string1, string2);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::vector_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_vector3f(string,
                                                                                                               Vector3f {
                                                                                                                   x: floats[0],
                                                                                                                   y: floats[1],
                                                                                                                   z: floats[2],
                                                                                                               });
                                                                                    } else {
                                                                                        param_set.add_point3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            _ => unreachable!(),
                                                                        };
                                                                    }
                                                                }
                                                                _ => unreachable!(),
                                                            }
                                                        }
                                                        // we should have the area_light_source parameters by now
                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                            println!("AreaLightSource \"{}\" ",
                                                                     param_set.name);
                                                            print_params(&param_set);
                                                            if let Some(ref mut graphics_state) =
                                                                GRAPHICS_STATE {
                                                                graphics_state.area_light =
                                                                    param_set.name.clone();
                                                                graphics_state
                                                                    .area_light_params
                                                                    .copy_from(&param_set);
                                                            }
                                                        } else {
                                                            panic!("Can't get parameter set.");
                                                        }
                                                    }
                                                    Rule::camera => {
                                                        for camera_pair in named_statement_pair
                                                                .into_inner() {
                                                            match camera_pair.as_rule() {
                                                                Rule::string => {
                                                                    let mut string_pairs =
                                                                        camera_pair.into_inner();
                                                                    let ident = string_pairs.next();
                                                                    let name: String =
                                                                        String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
                                                                    if let Some(ref mut ro) =
                                                                        RENDER_OPTIONS {
                                                                        ro.camera_name = name;
                                                                        ro.camera_to_world.t[0] =
                                                                            Transform::inverse(CUR_TRANSFORM.t[0]);
                                                                        ro.camera_to_world.t[1] =
                                                                            Transform::inverse(CUR_TRANSFORM.t[1]);
                                                                        if let Some(ref mut named_coordinate_systems) = NAMED_COORDINATE_SYSTEMS {
                                                                            named_coordinate_systems.insert("camera",
                                                                                                            TransformSet {
                                                                                                                t: [ro.camera_to_world.t[0],
                                                                                                                    ro.camera_to_world.t[1]]
                                                                                                            });
                                                                        }
                                                                    } else {
                                                                        panic!("Can't get render options.");
                                                                    }
                                                                    if let Some(ref mut param_set) = PARAM_SET {
                                                                        param_set.reset(String::from("Camera"),
                                                                                        String::from(""),
                                                                                        String::from(""),
                                                                                        String::from(""));
                                                                    } else {
                                                                        panic!("Can't get parameter set.");
                                                                    }
                                                                }
                                                                Rule::parameter => {
                                                                    for parameter_pair in
                                                                        camera_pair.into_inner() {
                                                                        match parameter_pair
                                                                                  .as_rule() {
                                                                            Rule::bool_param => {
                                                                                let tuple: (String, bool) = pbrt_bool_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let b: bool = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_bool(string, b);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::blackbody_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_blackbody_spectrum(string, floats);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::float_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 1 {
                                                                                        param_set.add_float(string, floats[0]);
                                                                                    } else {
                                                                                        param_set.add_floats(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::integer_param => {
                                                                                let tuple: (String, Vec<i32>) = pbrt_integer_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let integers: Vec<i32> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if integers.len() == 1 {
                                                                                        param_set.add_int(string, integers[0]);
                                                                                    } else {
                                                                                        param_set.add_ints(string, integers);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::point_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_point3f(string,
                                                                                                              Point3f {
                                                                                                                  x: floats[0],
                                                                                                                  y: floats[1],
                                                                                                                  z: floats[2],
                                                                                                              });
                                                                                    } else {
                                                                                        param_set.add_point3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::normal_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_normal3f(string,
                                                                                                               Normal3f {
                                                                                                                   x: floats[0],
                                                                                                                   y: floats[1],
                                                                                                                   z: floats[2],
                                                                                                               });
                                                                                    } else {
                                                                                        param_set.add_normal3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::rgb_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_rgb_spectrum(string,
                                                                                                               Spectrum {
                                                                                                                   c: [floats[0], floats[1], floats[2]],
                                                                                                               });
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::spectrum_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_rgb_spectrum(string,
                                                                                                               Spectrum {
                                                                                                                   c: [floats[0], floats[1], floats[2]],
                                                                                                               });
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::string_param => {
                                                                                let tuple: (String, String) = pbrt_string_parameter(&mut parameter_pair.into_inner());
                                                                                let string1: String = tuple.0;
                                                                                let string2: String = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_string(string1, string2);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::texture_param => {
                                                                                let tuple: (String, String) = pbrt_texture_parameter(&mut parameter_pair.into_inner());
                                                                                let string1: String = tuple.0;
                                                                                let string2: String = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_texture(string1, string2);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::vector_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_vector3f(string,
                                                                                                               Vector3f {
                                                                                                                   x: floats[0],
                                                                                                                   y: floats[1],
                                                                                                                   z: floats[2],
                                                                                                               });
                                                                                    } else {
                                                                                        param_set.add_point3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            _ => unreachable!(),
                                                                        };
                                                                    }
                                                                }
                                                                _ => unreachable!(),
                                                            }
                                                        }
                                                        // we should have the camera parameters by now
                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                            if let Some(ref mut ro) =
                                                                RENDER_OPTIONS {
                                                                println!("Camera \"{}\" ",
                                                                         ro.camera_name);
                                                                ro.camera_params
                                                                    .copy_from(param_set);
                                                                print_params(&ro.camera_params);
                                                            } else {
                                                                panic!("Can't get render options.");
                                                            }
                                                        } else {
                                                            panic!("Can't get parameter set.");
                                                        }
                                                    }
                                                    Rule::coord_sys_transform => {
                                                        for coord_sys_transform_pair in
                                                            named_statement_pair.into_inner() {
                                                            match coord_sys_transform_pair
                                                                      .as_rule() {
                                                                Rule::string => {
                                                                    let mut string_pairs =
                                                                        coord_sys_transform_pair
                                                                            .into_inner();
                                                                    let ident = string_pairs.next();
                                                                    let name: String =
                                                                        String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
                                                                    if let Some(ref mut named_coordinate_systems) = NAMED_COORDINATE_SYSTEMS {
                                                                        match named_coordinate_systems.get(name.as_str()) {
                                                                            Some(transform_set) => {
                                                                                CUR_TRANSFORM.t[0] = transform_set.t[0];
                                                                                CUR_TRANSFORM.t[1] = transform_set.t[1];
                                                                            },
                                                                            None => {
                                                                                println!("Couldn't find named coordinate system \"{}\"", name);
                                                                            },
                                                                        };
                                                                    }
                                                                }
                                                                _ => unreachable!(),
                                                            }
                                                        }
                                                    }
                                                    Rule::film => {
                                                        for film_pair in named_statement_pair
                                                                .into_inner() {
                                                            match film_pair.as_rule() {
                                                                Rule::string => {
                                                                    let mut string_pairs =
                                                                        film_pair.into_inner();
                                                                    let ident = string_pairs.next();
                                                                    let name: String =
                                                                        String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
                                                                    if let Some(ref mut ro) =
                                                                        RENDER_OPTIONS {
                                                                        ro.film_name = name;
                                                                    } else {
                                                                        panic!("Can't get render options.");
                                                                    }
                                                                    if let Some(ref mut param_set) = PARAM_SET {
                                                                        param_set.reset(String::from("Film"),
                                                                                        String::from(""),
                                                                                        String::from(""),
                                                                                        String::from(""));
                                                                    } else {
                                                                        panic!("Can't get parameter set.");
                                                                    }
                                                                }
                                                                Rule::parameter => {
                                                                    for parameter_pair in
                                                                        film_pair.into_inner() {
                                                                        match parameter_pair
                                                                                  .as_rule() {
                                                                            Rule::bool_param => {
                                                                                let tuple: (String, bool) = pbrt_bool_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let b: bool = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_bool(string, b);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::blackbody_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_blackbody_spectrum(string, floats);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::float_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 1 {
                                                                                        param_set.add_float(string, floats[0]);
                                                                                    } else {
                                                                                        param_set.add_floats(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::integer_param => {
                                                                                let tuple: (String, Vec<i32>) = pbrt_integer_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let integers: Vec<i32> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if integers.len() == 1 {
                                                                                        param_set.add_int(string, integers[0]);
                                                                                    } else {
                                                                                        param_set.add_ints(string, integers);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::point_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_point3f(string,
                                                                                                              Point3f {
                                                                                                                  x: floats[0],
                                                                                                                  y: floats[1],
                                                                                                                  z: floats[2],
                                                                                                              });
                                                                                    } else {
                                                                                        param_set.add_point3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::normal_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_normal3f(string,
                                                                                                               Normal3f {
                                                                                                                   x: floats[0],
                                                                                                                   y: floats[1],
                                                                                                                   z: floats[2],
                                                                                                               });
                                                                                    } else {
                                                                                        param_set.add_normal3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::rgb_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_rgb_spectrum(string,
                                                                                                               Spectrum {
                                                                                                                   c: [floats[0], floats[1], floats[2]],
                                                                                                               });
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::spectrum_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_rgb_spectrum(string,
                                                                                                               Spectrum {
                                                                                                                   c: [floats[0], floats[1], floats[2]],
                                                                                                               });
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::string_param => {
                                                                                let tuple: (String, String) = pbrt_string_parameter(&mut parameter_pair.into_inner());
                                                                                let string1: String = tuple.0;
                                                                                let string2: String = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_string(string1, string2);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::texture_param => {
                                                                                let tuple: (String, String) = pbrt_texture_parameter(&mut parameter_pair.into_inner());
                                                                                let string1: String = tuple.0;
                                                                                let string2: String = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_texture(string1, string2);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::vector_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_vector3f(string,
                                                                                                               Vector3f {
                                                                                                                   x: floats[0],
                                                                                                                   y: floats[1],
                                                                                                                   z: floats[2],
                                                                                                               });
                                                                                    } else {
                                                                                        param_set.add_point3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            _ => unreachable!(),
                                                                        };
                                                                    }
                                                                }
                                                                _ => unreachable!(),
                                                            }
                                                        }
                                                        // we should have the film parameters by now
                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                            if let Some(ref mut ro) =
                                                                RENDER_OPTIONS {
                                                                println!("Film \"{}\" ",
                                                                         ro.film_name);
                                                                ro.film_params.copy_from(param_set);
                                                                print_params(&ro.film_params);
                                                            } else {
                                                                panic!("Can't get render options.");
                                                            }
                                                        } else {
                                                            panic!("Can't get parameter set.");
                                                        }
                                                    }
                                                    Rule::integrator => {
                                                        for integrator_pair in
                                                            named_statement_pair.into_inner() {
                                                            match integrator_pair.as_rule() {
                                                                Rule::string => {
                                                                    let mut string_pairs =
                                                                        integrator_pair
                                                                            .into_inner();
                                                                    let ident = string_pairs.next();
                                                                    let name: String =
                                                                        String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
                                                                    if let Some(ref mut ro) =
                                                                        RENDER_OPTIONS {
                                                                        ro.integrator_name = name;
                                                                    } else {
                                                                        panic!("Can't get render options.");
                                                                    }
                                                                    if let Some(ref mut param_set) = PARAM_SET {
                                                                        param_set.reset(String::from("Integrator"),
                                                                                        String::from(""),
                                                                                        String::from(""),
                                                                                        String::from(""));
                                                                    } else {
                                                                        panic!("Can't get parameter set.");
                                                                    }
                                                                }
                                                                Rule::parameter => {
                                                                    for parameter_pair in
                                                                        integrator_pair.into_inner() {
                                                                        match parameter_pair
                                                                                  .as_rule() {
                                                                            Rule::bool_param => {
                                                                                let tuple: (String, bool) = pbrt_bool_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let b: bool = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_bool(string, b);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::blackbody_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_blackbody_spectrum(string, floats);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::float_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 1 {
                                                                                        param_set.add_float(string, floats[0]);
                                                                                    } else {
                                                                                        param_set.add_floats(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::integer_param => {
                                                                                let tuple: (String, Vec<i32>) = pbrt_integer_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let integers: Vec<i32> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if integers.len() == 1 {
                                                                                        param_set.add_int(string, integers[0]);
                                                                                    } else {
                                                                                        param_set.add_ints(string, integers);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::point_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_point3f(string,
                                                                                                              Point3f {
                                                                                                                  x: floats[0],
                                                                                                                  y: floats[1],
                                                                                                                  z: floats[2],
                                                                                                              });
                                                                                    } else {
                                                                                        param_set.add_point3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::normal_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_normal3f(string,
                                                                                                               Normal3f {
                                                                                                                   x: floats[0],
                                                                                                                   y: floats[1],
                                                                                                                   z: floats[2],
                                                                                                               });
                                                                                    } else {
                                                                                        param_set.add_normal3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::rgb_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_rgb_spectrum(string,
                                                                                                               Spectrum {
                                                                                                                   c: [floats[0], floats[1], floats[2]],
                                                                                                               });
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::spectrum_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_rgb_spectrum(string,
                                                                                                               Spectrum {
                                                                                                                   c: [floats[0], floats[1], floats[2]],
                                                                                                               });
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::string_param => {
                                                                                let tuple: (String, String) = pbrt_string_parameter(&mut parameter_pair.into_inner());
                                                                                let string1: String = tuple.0;
                                                                                let string2: String = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_string(string1, string2);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::texture_param => {
                                                                                let tuple: (String, String) = pbrt_texture_parameter(&mut parameter_pair.into_inner());
                                                                                let string1: String = tuple.0;
                                                                                let string2: String = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_texture(string1, string2);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::vector_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_vector3f(string,
                                                                                                               Vector3f {
                                                                                                                   x: floats[0],
                                                                                                                   y: floats[1],
                                                                                                                   z: floats[2],
                                                                                                               });
                                                                                    } else {
                                                                                        param_set.add_point3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            _ => unreachable!(),
                                                                        };
                                                                    }
                                                                }
                                                                _ => unreachable!(),
                                                            }
                                                        }
                                                        // we should have the integrator parameters by now
                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                            if let Some(ref mut ro) =
                                                                RENDER_OPTIONS {
                                                                println!("Integrator \"{}\" ",
                                                                         ro.integrator_name);
                                                                ro.integrator_params
                                                                    .copy_from(param_set);
                                                                print_params(&ro.integrator_params);
                                                            } else {
                                                                panic!("Can't get render options.");
                                                            }
                                                        } else {
                                                            panic!("Can't get parameter set.");
                                                        }
                                                    }
                                                    Rule::light_source => {
                                                        for light_source_pair in
                                                            named_statement_pair.into_inner() {
                                                            match light_source_pair.as_rule() {
                                                                Rule::string => {
                                                                    let mut string_pairs =
                                                                        light_source_pair
                                                                            .into_inner();
                                                                    let ident = string_pairs.next();
                                                                    let name: String =
                                                                        String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
                                                                    if let Some(ref mut param_set) = PARAM_SET {
                                                                        param_set.reset(String::from("Light_Source"),
                                                                                        String::from(name),
                                                                                        String::from(""),
                                                                                        String::from(""));
                                                                    } else {
                                                                        panic!("Can't get parameter set.");
                                                                    }
                                                                }
                                                                Rule::parameter => {
                                                                    for parameter_pair in
                                                                        light_source_pair
                                                                            .into_inner() {
                                                                        match parameter_pair
                                                                                  .as_rule() {
                                                                            Rule::bool_param => {
                                                                                let tuple: (String, bool) = pbrt_bool_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let b: bool = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_bool(string, b);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::blackbody_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_blackbody_spectrum(string, floats);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::float_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 1 {
                                                                                        param_set.add_float(string, floats[0]);
                                                                                    } else {
                                                                                        param_set.add_floats(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::integer_param => {
                                                                                let tuple: (String, Vec<i32>) = pbrt_integer_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let integers: Vec<i32> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if integers.len() == 1 {
                                                                                        param_set.add_int(string, integers[0]);
                                                                                    } else {
                                                                                        param_set.add_ints(string, integers);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::point_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_point3f(string,
                                                                                                              Point3f {
                                                                                                                  x: floats[0],
                                                                                                                  y: floats[1],
                                                                                                                  z: floats[2],
                                                                                                              });
                                                                                    } else {
                                                                                        param_set.add_point3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::normal_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_normal3f(string,
                                                                                                               Normal3f {
                                                                                                                   x: floats[0],
                                                                                                                   y: floats[1],
                                                                                                                   z: floats[2],
                                                                                                               });
                                                                                    } else {
                                                                                        param_set.add_normal3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::rgb_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_rgb_spectrum(string,
                                                                                                               Spectrum {
                                                                                                                   c: [floats[0], floats[1], floats[2]],
                                                                                                               });
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::spectrum_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_rgb_spectrum(string,
                                                                                                               Spectrum {
                                                                                                                   c: [floats[0], floats[1], floats[2]],
                                                                                                               });
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::string_param => {
                                                                                let tuple: (String, String) = pbrt_string_parameter(&mut parameter_pair.into_inner());
                                                                                let string1: String = tuple.0;
                                                                                let string2: String = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_string(string1, string2);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::texture_param => {
                                                                                let tuple: (String, String) = pbrt_texture_parameter(&mut parameter_pair.into_inner());
                                                                                let string1: String = tuple.0;
                                                                                let string2: String = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_texture(string1, string2);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::vector_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_vector3f(string,
                                                                                                               Vector3f {
                                                                                                                   x: floats[0],
                                                                                                                   y: floats[1],
                                                                                                                   z: floats[2],
                                                                                                               });
                                                                                    } else {
                                                                                        param_set.add_point3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            _ => unreachable!(),
                                                                        };
                                                                    }
                                                                }
                                                                _ => unreachable!(),
                                                            }
                                                        }
                                                        // we should have the light_source parameters by now
                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                            if let Some(ref mut ro) =
                                                                RENDER_OPTIONS {
                                                                println!("LightSource \"{}\" ",
                                                                         param_set.name);
                                                                print_params(&param_set);
                                                                make_light(&param_set, ro);
                                                            } else {
                                                                panic!("Can't get render options.");
                                                            }
                                                        } else {
                                                            panic!("Can't get parameter set.");
                                                        }
                                                    }
                                                    Rule::make_named_material => {
                                                        for make_named_material_pair in
                                                            named_statement_pair.into_inner() {
                                                            match make_named_material_pair
                                                                      .as_rule() {
                                                                Rule::string => {
                                                                    let mut string_pairs =
                                                                        make_named_material_pair
                                                                            .into_inner();
                                                                    let ident = string_pairs.next();
                                                                    let name: String =
                                                                        String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
                                                                    if let Some(ref mut param_set) = PARAM_SET {
                                                                        param_set.reset(String::from("MakeNamedMaterial"),
                                                                                        String::from(name.clone()),
                                                                                        String::from(""),
                                                                                        String::from(""));
                                                                    } else {
                                                                        panic!("Can't get parameter set.");
                                                                    }
                                                                }
                                                                Rule::parameter => {
                                                                    for parameter_pair in
                                                                        make_named_material_pair
                                                                            .into_inner() {
                                                                        match parameter_pair
                                                                                  .as_rule() {
                                                                            Rule::bool_param => {
                                                                                let tuple: (String, bool) = pbrt_bool_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let b: bool = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_bool(string, b);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::blackbody_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_blackbody_spectrum(string, floats);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::float_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 1 {
                                                                                        param_set.add_float(string, floats[0]);
                                                                                    } else {
                                                                                        param_set.add_floats(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::integer_param => {
                                                                                let tuple: (String, Vec<i32>) = pbrt_integer_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let integers: Vec<i32> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if integers.len() == 1 {
                                                                                        param_set.add_int(string, integers[0]);
                                                                                    } else {
                                                                                        param_set.add_ints(string, integers);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::point_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_point3f(string,
                                                                                                              Point3f {
                                                                                                                  x: floats[0],
                                                                                                                  y: floats[1],
                                                                                                                  z: floats[2],
                                                                                                              });
                                                                                    } else {
                                                                                        param_set.add_point3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::normal_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_normal3f(string,
                                                                                                               Normal3f {
                                                                                                                   x: floats[0],
                                                                                                                   y: floats[1],
                                                                                                                   z: floats[2],
                                                                                                               });
                                                                                    } else {
                                                                                        param_set.add_normal3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::rgb_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_rgb_spectrum(string,
                                                                                                               Spectrum {
                                                                                                                   c: [floats[0], floats[1], floats[2]],
                                                                                                               });
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::spectrum_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_rgb_spectrum(string,
                                                                                                               Spectrum {
                                                                                                                   c: [floats[0], floats[1], floats[2]],
                                                                                                               });
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::string_param => {
                                                                                let tuple: (String, String) = pbrt_string_parameter(&mut parameter_pair.into_inner());
                                                                                let string1: String = tuple.0;
                                                                                let string2: String = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_string(string1, string2);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::texture_param => {
                                                                                let tuple: (String, String) = pbrt_texture_parameter(&mut parameter_pair.into_inner());
                                                                                let string1: String = tuple.0;
                                                                                let string2: String = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_texture(string1, string2);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::vector_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_vector3f(string,
                                                                                                               Vector3f {
                                                                                                                   x: floats[0],
                                                                                                                   y: floats[1],
                                                                                                                   z: floats[2],
                                                                                                               });
                                                                                    } else {
                                                                                        param_set.add_point3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            _ => unreachable!(),
                                                                        };
                                                                    }
                                                                }
                                                                _ => unreachable!(),
                                                            }
                                                        }
                                                        // we should have the make_named_material parameters by now
                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                            let mat_type: String = param_set.find_one_string(String::from("type"),
                                                                                                             String::new());
                                                            if mat_type == String::new() {
                                                                panic!("No parameter string \"type\" found in MakeNamedMaterial");
                                                            }
                                                            if let Some(ref mut graphics_state) =
                                                                GRAPHICS_STATE {
                                                                graphics_state.material =
                                                                    mat_type.clone();
                                                                graphics_state
                                                                    .material_params
                                                                    .copy_from(&param_set);
                                                                graphics_state.current_material = String::new();
                                                                let mtl: Arc<Material + Send + Sync> = create_material();
                                                                match graphics_state.named_materials.get(param_set.name.as_str()) {
                                                                    Some(_named_material) => {
                                                                        println!("Named material \"{}\" redefined",
                                                                                 mat_type);
                                                                    },
                                                                    None => {},
                                                                }
                                                                graphics_state
                                                                    .named_materials
                                                                    .insert(param_set.name.clone(),
                                                                            mtl);
                                                            }
                                                        } else {
                                                            panic!("Can't get parameter set.");
                                                        }
                                                    }
                                                    Rule::material => {
                                                        for material_pair in named_statement_pair
                                                                .into_inner() {
                                                            match material_pair.as_rule() {
                                                                Rule::string => {
                                                                    let mut string_pairs =
                                                                        material_pair.into_inner();
                                                                    let ident = string_pairs.next();
                                                                    let name: String =
                                                                        String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
                                                                    if let Some(ref mut param_set) = PARAM_SET {
                                                                        param_set.reset(String::from("Material"),
                                                                                        String::from(name.clone()),
                                                                                        String::from(""),
                                                                                        String::from(""));
                                                                    } else {
                                                                        panic!("Can't get parameter set.");
                                                                    }
                                                                }
                                                                Rule::parameter => {
                                                                    for parameter_pair in
                                                                        material_pair.into_inner() {
                                                                        match parameter_pair
                                                                                  .as_rule() {
                                                                            Rule::bool_param => {
                                                                                let tuple: (String, bool) = pbrt_bool_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let b: bool = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_bool(string, b);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::blackbody_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_blackbody_spectrum(string, floats);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::float_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 1 {
                                                                                        param_set.add_float(string, floats[0]);
                                                                                    } else {
                                                                                        param_set.add_floats(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::integer_param => {
                                                                                let tuple: (String, Vec<i32>) = pbrt_integer_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let integers: Vec<i32> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if integers.len() == 1 {
                                                                                        param_set.add_int(string, integers[0]);
                                                                                    } else {
                                                                                        param_set.add_ints(string, integers);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::point_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_point3f(string,
                                                                                                              Point3f {
                                                                                                                  x: floats[0],
                                                                                                                  y: floats[1],
                                                                                                                  z: floats[2],
                                                                                                              });
                                                                                    } else {
                                                                                        param_set.add_point3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::normal_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_normal3f(string,
                                                                                                               Normal3f {
                                                                                                                   x: floats[0],
                                                                                                                   y: floats[1],
                                                                                                                   z: floats[2],
                                                                                                               });
                                                                                    } else {
                                                                                        param_set.add_normal3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::rgb_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_rgb_spectrum(string,
                                                                                                               Spectrum {
                                                                                                                   c: [floats[0], floats[1], floats[2]],
                                                                                                               });
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::spectrum_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_rgb_spectrum(string,
                                                                                                               Spectrum {
                                                                                                                   c: [floats[0], floats[1], floats[2]],
                                                                                                               });
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::string_param => {
                                                                                let tuple: (String, String) = pbrt_string_parameter(&mut parameter_pair.into_inner());
                                                                                let string1: String = tuple.0;
                                                                                let string2: String = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_string(string1, string2);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::texture_param => {
                                                                                let tuple: (String, String) = pbrt_texture_parameter(&mut parameter_pair.into_inner());
                                                                                let string1: String = tuple.0;
                                                                                let string2: String = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_texture(string1, string2);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::vector_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_vector3f(string,
                                                                                                               Vector3f {
                                                                                                                   x: floats[0],
                                                                                                                   y: floats[1],
                                                                                                                   z: floats[2],
                                                                                                               });
                                                                                    } else {
                                                                                        param_set.add_point3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            _ => unreachable!(),
                                                                        };
                                                                    }
                                                                }
                                                                _ => unreachable!(),
                                                            }
                                                        }
                                                        // we should have the material parameters by now
                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                            if let Some(ref mut graphics_state) =
                                                                GRAPHICS_STATE {
                                                                graphics_state.material =
                                                                    param_set.name.clone();
                                                                graphics_state
                                                                    .material_params
                                                                    .copy_from(&param_set);
                                                                graphics_state.current_material = String::new();
                                                            }
                                                        } else {
                                                            panic!("Can't get parameter set.");
                                                        }
                                                    }
                                                    Rule::named_material => {
                                                        for named_material_pair in
                                                            named_statement_pair.into_inner() {
                                                            match named_material_pair.as_rule() {
                                                                Rule::string => {
                                                                    let mut string_pairs =
                                                                        named_material_pair
                                                                            .into_inner();
                                                                    let ident = string_pairs.next();
                                                                    let name: String =
                                                                        String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
                                                                    if let Some(ref mut param_set) = PARAM_SET {
                                                                        param_set.reset(String::from("MakeNamedMaterial"),
                                                                                        String::from(name.clone()),
                                                                                        String::from(""),
                                                                                        String::from(""));
                                                                    } else {
                                                                        panic!("Can't get parameter set.");
                                                                    }
                                                                }
                                                                Rule::parameter => {
                                                                    for parameter_pair in
                                                                        named_material_pair
                                                                            .into_inner() {
                                                                        match parameter_pair
                                                                                  .as_rule() {
                                                                            Rule::bool_param => {
                                                                                let tuple: (String, bool) = pbrt_bool_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let b: bool = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_bool(string, b);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::blackbody_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_blackbody_spectrum(string, floats);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::float_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 1 {
                                                                                        param_set.add_float(string, floats[0]);
                                                                                    } else {
                                                                                        param_set.add_floats(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::integer_param => {
                                                                                let tuple: (String, Vec<i32>) = pbrt_integer_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let integers: Vec<i32> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if integers.len() == 1 {
                                                                                        param_set.add_int(string, integers[0]);
                                                                                    } else {
                                                                                        param_set.add_ints(string, integers);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::point_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_point3f(string,
                                                                                                              Point3f {
                                                                                                                  x: floats[0],
                                                                                                                  y: floats[1],
                                                                                                                  z: floats[2],
                                                                                                              });
                                                                                    } else {
                                                                                        param_set.add_point3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::normal_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_normal3f(string,
                                                                                                               Normal3f {
                                                                                                                   x: floats[0],
                                                                                                                   y: floats[1],
                                                                                                                   z: floats[2],
                                                                                                               });
                                                                                    } else {
                                                                                        param_set.add_normal3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::rgb_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_rgb_spectrum(string,
                                                                                                               Spectrum {
                                                                                                                   c: [floats[0], floats[1], floats[2]],
                                                                                                               });
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::spectrum_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_rgb_spectrum(string,
                                                                                                               Spectrum {
                                                                                                                   c: [floats[0], floats[1], floats[2]],
                                                                                                               });
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::string_param => {
                                                                                let tuple: (String, String) = pbrt_string_parameter(&mut parameter_pair.into_inner());
                                                                                let string1: String = tuple.0;
                                                                                let string2: String = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_string(string1, string2);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::texture_param => {
                                                                                let tuple: (String, String) = pbrt_texture_parameter(&mut parameter_pair.into_inner());
                                                                                let string1: String = tuple.0;
                                                                                let string2: String = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_texture(string1, string2);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::vector_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_vector3f(string,
                                                                                                               Vector3f {
                                                                                                                   x: floats[0],
                                                                                                                   y: floats[1],
                                                                                                                   z: floats[2],
                                                                                                               });
                                                                                    } else {
                                                                                        param_set.add_point3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            _ => unreachable!(),
                                                                        };
                                                                    }
                                                                }
                                                                _ => unreachable!(),
                                                            }
                                                        }
                                                        // we should have the named_material parameters by now
                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                            if let Some(ref mut graphics_state) =
                                                                GRAPHICS_STATE {
                                                                graphics_state.current_material = param_set.name.clone();
                                                            }
                                                        } else {
                                                            panic!("Can't get parameter set.");
                                                        }
                                                    }
                                                    Rule::pixel_filter => {
                                                        for pixel_filter_pair in
                                                            named_statement_pair.into_inner() {
                                                            match pixel_filter_pair.as_rule() {
                                                                Rule::string => {
                                                                    let mut string_pairs =
                                                                        pixel_filter_pair
                                                                            .into_inner();
                                                                    let ident = string_pairs.next();
                                                                    let name: String =
                                                                        String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
                                                                    if let Some(ref mut ro) =
                                                                        RENDER_OPTIONS {
                                                                        ro.filter_name = name;
                                                                    } else {
                                                                        panic!("Can't get render options.");
                                                                    }
                                                                    if let Some(ref mut param_set) = PARAM_SET {
                                                                        param_set.reset(String::from("PixelFilter"),
                                                                                        String::from(""),
                                                                                        String::from(""),
                                                                                        String::from(""));
                                                                    } else {
                                                                        panic!("Can't get parameter set.");
                                                                    }
                                                                }
                                                                Rule::parameter => {
                                                                    for parameter_pair in
                                                                        pixel_filter_pair
                                                                            .into_inner() {
                                                                        match parameter_pair
                                                                                  .as_rule() {
                                                                            Rule::bool_param => {
                                                                                let tuple: (String, bool) = pbrt_bool_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let b: bool = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_bool(string, b);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::blackbody_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_blackbody_spectrum(string, floats);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::float_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 1 {
                                                                                        param_set.add_float(string, floats[0]);
                                                                                    } else {
                                                                                        param_set.add_floats(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::integer_param => {
                                                                                let tuple: (String, Vec<i32>) = pbrt_integer_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let integers: Vec<i32> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if integers.len() == 1 {
                                                                                        param_set.add_int(string, integers[0]);
                                                                                    } else {
                                                                                        param_set.add_ints(string, integers);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::point_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_point3f(string,
                                                                                                              Point3f {
                                                                                                                  x: floats[0],
                                                                                                                  y: floats[1],
                                                                                                                  z: floats[2],
                                                                                                              });
                                                                                    } else {
                                                                                        param_set.add_point3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::normal_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_normal3f(string,
                                                                                                               Normal3f {
                                                                                                                   x: floats[0],
                                                                                                                   y: floats[1],
                                                                                                                   z: floats[2],
                                                                                                               });
                                                                                    } else {
                                                                                        param_set.add_normal3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::rgb_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_rgb_spectrum(string,
                                                                                                               Spectrum {
                                                                                                                   c: [floats[0], floats[1], floats[2]],
                                                                                                               });
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::spectrum_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_rgb_spectrum(string,
                                                                                                               Spectrum {
                                                                                                                   c: [floats[0], floats[1], floats[2]],
                                                                                                               });
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::string_param => {
                                                                                let tuple: (String, String) = pbrt_string_parameter(&mut parameter_pair.into_inner());
                                                                                let string1: String = tuple.0;
                                                                                let string2: String = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_string(string1, string2);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::texture_param => {
                                                                                let tuple: (String, String) = pbrt_texture_parameter(&mut parameter_pair.into_inner());
                                                                                let string1: String = tuple.0;
                                                                                let string2: String = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_texture(string1, string2);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::vector_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_vector3f(string,
                                                                                                               Vector3f {
                                                                                                                   x: floats[0],
                                                                                                                   y: floats[1],
                                                                                                                   z: floats[2],
                                                                                                               });
                                                                                    } else {
                                                                                        param_set.add_point3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            _ => unreachable!(),
                                                                        };
                                                                    }
                                                                }
                                                                _ => unreachable!(),
                                                            }
                                                        }
                                                        // we should have the pixel_filter parameters by now
                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                            if let Some(ref mut ro) =
                                                                RENDER_OPTIONS {
                                                                println!("PixelFilter \"{}\" ",
                                                                         ro.filter_name);
                                                                ro.filter_params
                                                                    .copy_from(param_set);
                                                                print_params(&ro.filter_params);
                                                            } else {
                                                                panic!("Can't get render options.");
                                                            }
                                                        } else {
                                                            panic!("Can't get parameter set.");
                                                        }
                                                    }
                                                    Rule::sampler => {
                                                        for sampler_pair in named_statement_pair
                                                                .into_inner() {
                                                            match sampler_pair.as_rule() {
                                                                Rule::string => {
                                                                    let mut string_pairs =
                                                                        sampler_pair.into_inner();
                                                                    let ident = string_pairs.next();
                                                                    let name: String =
                                                                        String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
                                                                    if let Some(ref mut ro) =
                                                                        RENDER_OPTIONS {
                                                                        ro.sampler_name = name;
                                                                    } else {
                                                                        panic!("Can't get render options.");
                                                                    }
                                                                    if let Some(ref mut param_set) = PARAM_SET {
                                                                        param_set.reset(String::from("Sampler"),
                                                                                        String::from(""),
                                                                                        String::from(""),
                                                                                        String::from(""));
                                                                    } else {
                                                                        panic!("Can't get parameter set.");
                                                                    }
                                                                }
                                                                Rule::parameter => {
                                                                    for parameter_pair in
                                                                        sampler_pair.into_inner() {
                                                                        match parameter_pair
                                                                                  .as_rule() {
                                                                            Rule::bool_param => {
                                                                                let tuple: (String, bool) = pbrt_bool_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let b: bool = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_bool(string, b);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::blackbody_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_blackbody_spectrum(string, floats);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::float_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 1 {
                                                                                        param_set.add_float(string, floats[0]);
                                                                                    } else {
                                                                                        param_set.add_floats(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::integer_param => {
                                                                                let tuple: (String, Vec<i32>) = pbrt_integer_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let integers: Vec<i32> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if integers.len() == 1 {
                                                                                        param_set.add_int(string, integers[0]);
                                                                                    } else {
                                                                                        param_set.add_ints(string, integers);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::point_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_point3f(string,
                                                                                                              Point3f {
                                                                                                                  x: floats[0],
                                                                                                                  y: floats[1],
                                                                                                                  z: floats[2],
                                                                                                              });
                                                                                    } else {
                                                                                        param_set.add_point3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::normal_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_normal3f(string,
                                                                                                               Normal3f {
                                                                                                                   x: floats[0],
                                                                                                                   y: floats[1],
                                                                                                                   z: floats[2],
                                                                                                               });
                                                                                    } else {
                                                                                        param_set.add_normal3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::rgb_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_rgb_spectrum(string,
                                                                                                               Spectrum {
                                                                                                                   c: [floats[0], floats[1], floats[2]],
                                                                                                               });
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::spectrum_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_rgb_spectrum(string,
                                                                                                               Spectrum {
                                                                                                                   c: [floats[0], floats[1], floats[2]],
                                                                                                               });
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::string_param => {
                                                                                let tuple: (String, String) = pbrt_string_parameter(&mut parameter_pair.into_inner());
                                                                                let string1: String = tuple.0;
                                                                                let string2: String = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_string(string1, string2);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::texture_param => {
                                                                                let tuple: (String, String) = pbrt_texture_parameter(&mut parameter_pair.into_inner());
                                                                                let string1: String = tuple.0;
                                                                                let string2: String = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_texture(string1, string2);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::vector_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_vector3f(string,
                                                                                                               Vector3f {
                                                                                                                   x: floats[0],
                                                                                                                   y: floats[1],
                                                                                                                   z: floats[2],
                                                                                                               });
                                                                                    } else {
                                                                                        param_set.add_point3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            _ => unreachable!(),
                                                                        };
                                                                    }
                                                                }
                                                                _ => unreachable!(),
                                                            }
                                                        }
                                                        // we should have the sampler parameters by now
                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                            if let Some(ref mut ro) =
                                                                RENDER_OPTIONS {
                                                                println!("Sampler \"{}\" ",
                                                                         ro.sampler_name);
                                                                ro.sampler_params
                                                                    .copy_from(param_set);
                                                                print_params(&ro.sampler_params);
                                                            } else {
                                                                panic!("Can't get render options.");
                                                            }
                                                        } else {
                                                            panic!("Can't get parameter set.");
                                                        }
                                                    }
                                                    Rule::shape => {
                                                        for shape_pair in named_statement_pair
                                                                .into_inner() {
                                                            match shape_pair.as_rule() {
                                                                Rule::string => {
                                                                    let mut string_pairs =
                                                                        shape_pair.into_inner();
                                                                    let ident = string_pairs.next();
                                                                    let name: String =
                                                                        String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
                                                                    if let Some(ref mut param_set) = PARAM_SET {
                                                                        // WARNING: Reset BEFORE calling pbrt_shape() !
                                                                        param_set.reset(String::from("Shape"),
                                                                                        String::from(name.clone()),
                                                                                        String::from(""),
                                                                                        String::from(""));
                                                                    } else {
                                                                        panic!("Can't get parameter set.");
                                                                    }
                                                                }
                                                                Rule::parameter => {
                                                                    for parameter_pair in
                                                                        shape_pair.into_inner() {
                                                                        match parameter_pair
                                                                                  .as_rule() {
                                                                            Rule::bool_param => {
                                                                                let tuple: (String, bool) = pbrt_bool_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let b: bool = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_bool(string, b);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::blackbody_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_blackbody_spectrum(string, floats);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::float_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 1 {
                                                                                        param_set.add_float(string, floats[0]);
                                                                                    } else {
                                                                                        param_set.add_floats(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::integer_param => {
                                                                                let tuple: (String, Vec<i32>) = pbrt_integer_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let integers: Vec<i32> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if integers.len() == 1 {
                                                                                        param_set.add_int(string, integers[0]);
                                                                                    } else {
                                                                                        param_set.add_ints(string, integers);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::point_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_point3f(string,
                                                                                                              Point3f {
                                                                                                                  x: floats[0],
                                                                                                                  y: floats[1],
                                                                                                                  z: floats[2],
                                                                                                              });
                                                                                    } else {
                                                                                        param_set.add_point3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::normal_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_normal3f(string,
                                                                                                               Normal3f {
                                                                                                                   x: floats[0],
                                                                                                                   y: floats[1],
                                                                                                                   z: floats[2],
                                                                                                               });
                                                                                    } else {
                                                                                        param_set.add_normal3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::rgb_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_rgb_spectrum(string,
                                                                                                               Spectrum {
                                                                                                                   c: [floats[0], floats[1], floats[2]],
                                                                                                               });
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::spectrum_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_rgb_spectrum(string,
                                                                                                               Spectrum {
                                                                                                                   c: [floats[0], floats[1], floats[2]],
                                                                                                               });
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::string_param => {
                                                                                let tuple: (String, String) = pbrt_string_parameter(&mut parameter_pair.into_inner());
                                                                                let string1: String = tuple.0;
                                                                                let string2: String = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_string(string1, string2);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::texture_param => {
                                                                                let tuple: (String, String) = pbrt_texture_parameter(&mut parameter_pair.into_inner());
                                                                                let string1: String = tuple.0;
                                                                                let string2: String = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_texture(string1, string2);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::vector_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_vector3f(string,
                                                                                                               Vector3f {
                                                                                                                   x: floats[0],
                                                                                                                   y: floats[1],
                                                                                                                   z: floats[2],
                                                                                                               });
                                                                                    } else {
                                                                                        param_set.add_point3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            _ => unreachable!(),
                                                                        };
                                                                    }
                                                                }
                                                                _ => unreachable!(),
                                                            }
                                                        }
                                                        // we should have the shape parameters by now
                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                            // println!("Shape \"{}\" ", param_set.name);
                                                            // print_params(&param_set);
                                                            // collect area lights
                                                            let mut prims: Vec<Arc<Primitive + Send + Sync>> = Vec::new();
                                                            let mut area_lights: Vec<Arc<Light + Send + Sync>> = Vec::new();
                                                            // possibly create area light for shape (see pbrtShape())
                                                            if let Some(ref mut graphics_state) =
                                                                GRAPHICS_STATE {
                                                                if graphics_state.area_light !=
                                                                   String::new() {
                                                                    // MakeAreaLight
                                                                    if graphics_state.area_light ==
                                                                       String::from("area") ||
                                                                       graphics_state.area_light ==
                                                                       String::from("diffuse") {
                                                                        // first create the shape
                                                                        let (shapes, materials) =
                                                                            pbrt_shape(&param_set);
                                                                        assert_eq!(shapes.len(),
                                                                                   materials.len());
                                                                        for i in 0..shapes.len() {
                                                                            let shape = &shapes[i];
                                                                            let material =
                                                                                &materials[i];
                                                                            // CreateDiffuseAreaLight
                                                                            let light_to_world: Transform = CUR_TRANSFORM.t[0];
                                                                            let l: Spectrum =
                                                                                graphics_state.area_light_params.find_one_spectrum(String::from("L"),
                                                                                                                                   Spectrum::new(1.0));
                                                                            let sc: Spectrum =
                                                                                graphics_state.area_light_params.find_one_spectrum(String::from("scale"),
                                                                                                                                   Spectrum::new(1.0));
                                                                            let n_samples: i32 = // try "nsamples" first
                                                                                graphics_state.area_light_params.find_one_int(String::from("nsamples"),
                                                                                                                              1);
                                                                            let n_samples: i32 = // try "samples"next
                                                                                graphics_state.area_light_params.find_one_int(String::from("samples"),
                                                                                                                              n_samples);
                                                                            let two_sided: bool =
                                                                                graphics_state.area_light_params.find_one_bool(String::from("twosided"),
                                                                                                                               false);
                                                                            // TODO: if (PbrtOptions.quickRender) nSamples = std::max(1, nSamples / 4);
                                                                            let l_emit: Spectrum = l * sc;
                                                                            let area_light: Arc<DiffuseAreaLight> =
                                                                                Arc::new(DiffuseAreaLight::new(
                                                                                    &light_to_world,
                                                                                    &l_emit,
                                                                                    n_samples,
                                                                                    shape.clone(),
                                                                                    two_sided
                                                                                ));
                                                                            area_lights.push(area_light.clone());
                                                                            let geo_prim = Arc::new(GeometricPrimitive::new(shape.clone(),
                                                                                                                            material.clone(),
                                                                                                                            Some(area_light.clone())));
                                                                            prims.push(geo_prim.clone());
                                                                        }
                                                                    }
                                                                } else {
                                                                    // continue with shape itself
                                                                    let (shapes, materials) =
                                                                        pbrt_shape(&param_set);
                                                                    assert_eq!(shapes.len(),
                                                                               materials.len());
                                                                    for i in 0..shapes.len() {
                                                                        let shape = &shapes[i];
                                                                        let material = &materials
                                                                                            [i];
                                                                        let geo_prim = Arc::new(GeometricPrimitive::new(shape.clone(),
                                                                                                                        material.clone(),
                                                                                                                        None));
                                                                        prims
                                                                            .push(geo_prim.clone());
                                                                    }
                                                                    // animated?
                                                                    if CUR_TRANSFORM.is_animated() {
                                                                        if let Some(ref mut ro) =
                                                                            RENDER_OPTIONS {
                                                                            let animated_object_to_world: AnimatedTransform =
                                                                                AnimatedTransform::new(&CUR_TRANSFORM.t[0],
                                                                                                       ro.transform_start_time,
                                                                                                       &CUR_TRANSFORM.t[1],
                                                                                                       ro.transform_end_time);
                                                                            if prims.len() > 1 {
                                                                                println!("TODO: prims.len() > 1");
                                                                                let bvh: Arc<Primitive + Send + Sync> = Arc::new(BVHAccel::new(prims.clone(), 4, SplitMethod::SAH));
                                                                                prims.clear();
                                                                                prims.push(bvh.clone());
                                                                            } else {
                                                                                if let Some(primitive) = prims.pop() {
                                                                                    let geo_prim = Arc::new(TransformedPrimitive::new(primitive,
                                                                                                                                      animated_object_to_world));
                                                                                    prims.push(geo_prim.clone());
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            } else {
                                                                // continue with shape itself
                                                                let (shapes, materials) =
                                                                    pbrt_shape(&param_set);
                                                                assert_eq!(shapes.len(),
                                                                           materials.len());
                                                                for i in 0..shapes.len() {
                                                                    let shape = &shapes[i];
                                                                    let material = &materials[i];
                                                                    let geo_prim = Arc::new(GeometricPrimitive::new(shape.clone(),
                                                                                                                    material.clone(),
                                                                                                                    None));
                                                                    prims.push(geo_prim.clone());
                                                                }
                                                            }
                                                            // add _prims_ and _areaLights_ to scene or current instance
                                                            // if (renderOptions->currentInstance) {
                                                            //     if (areaLights.size())
                                                            //         Warning("Area lights not supported with object instancing");
                                                            //     renderOptions->currentInstance->insert(
                                                            //         renderOptions->currentInstance->end(), prims.begin(), prims.end());
                                                            // } else {
                                                            if let Some(ref mut ro) =
                                                                RENDER_OPTIONS {
                                                                // renderOptions->primitives.insert(renderOptions->primitives.end(),
                                                                //     prims.begin(), prims.end());
                                                                for prim in prims {
                                                                    ro.primitives
                                                                        .push(prim.clone());
                                                                }
                                                                // ro.primitives.insert(ro.primitives.end(),
                                                                //                      prims.begin(), prims.end());
                                                                if area_lights.len() > 0 {
                                                                    for area_light in area_lights {
                                                                        ro.lights.push(area_light);
                                                                    }
                                                                }
                                                            }
                                                        } else {
                                                            panic!("Can't get parameter set.");
                                                        }
                                                    }
                                                    Rule::texture => {
                                                        let mut counter: u8 = 0_u8;
                                                        let mut name: String = String::from("undefined");
                                                        let mut tex_type: String = String::from("undefined");
                                                        let mut tex_name;
                                                        for texture_pair in named_statement_pair
                                                                .into_inner() {
                                                            match texture_pair.as_rule() {
                                                                Rule::string => {
                                                                    match counter {
                                                                        0 => {
                                                                            // name
                                                                            let mut string_pairs =
                                                                                texture_pair
                                                                                    .into_inner();
                                                                            let ident =
                                                                                string_pairs.next();
                                                                            name = String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
                                                                        }
                                                                        1 => {
                                                                            // tex_type
                                                                            let mut string_pairs =
                                                                                texture_pair
                                                                                    .into_inner();
                                                                            let ident =
                                                                                string_pairs.next();
                                                                            tex_type = String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
                                                                        }
                                                                        2 => {
                                                                            // tex_name
                                                                            let mut string_pairs =
                                                                                texture_pair
                                                                                    .into_inner();
                                                                            let ident =
                                                                                string_pairs.next();
                                                                            tex_name = String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
                                                                            if let Some(ref mut param_set) = PARAM_SET {
                                                                                param_set.reset(String::from("Texture"),
                                                                                                String::from(name.clone()),
                                                                                                String::from(tex_type.clone()),
                                                                                                String::from(tex_name.clone()));
                                                                            } else {
                                                                                panic!("Can't get parameter set.");
                                                                            }
                                                                        }
                                                                        _ => unreachable!(),
                                                                    };
                                                                    counter += 1_u8;
                                                                }
                                                                Rule::parameter => {
                                                                    for parameter_pair in
                                                                        texture_pair.into_inner() {
                                                                        match parameter_pair
                                                                                  .as_rule() {
                                                                            Rule::bool_param => {
                                                                                let tuple: (String, bool) = pbrt_bool_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let b: bool = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_bool(string, b);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::blackbody_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_blackbody_spectrum(string, floats);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::float_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 1 {
                                                                                        param_set.add_float(string, floats[0]);
                                                                                    } else {
                                                                                        param_set.add_floats(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::integer_param => {
                                                                                let tuple: (String, Vec<i32>) = pbrt_integer_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let integers: Vec<i32> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if integers.len() == 1 {
                                                                                        param_set.add_int(string, integers[0]);
                                                                                    } else {
                                                                                        param_set.add_ints(string, integers);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::point_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_point3f(string,
                                                                                                              Point3f {
                                                                                                                  x: floats[0],
                                                                                                                  y: floats[1],
                                                                                                                  z: floats[2],
                                                                                                              });
                                                                                    } else {
                                                                                        param_set.add_point3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::normal_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_normal3f(string,
                                                                                                               Normal3f {
                                                                                                                   x: floats[0],
                                                                                                                   y: floats[1],
                                                                                                                   z: floats[2],
                                                                                                               });
                                                                                    } else {
                                                                                        param_set.add_normal3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::rgb_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_rgb_spectrum(string,
                                                                                                               Spectrum {
                                                                                                                   c: [floats[0], floats[1], floats[2]],
                                                                                                               });
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::spectrum_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_rgb_spectrum(string,
                                                                                                               Spectrum {
                                                                                                                   c: [floats[0], floats[1], floats[2]],
                                                                                                               });
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::string_param => {
                                                                                let tuple: (String, String) = pbrt_string_parameter(&mut parameter_pair.into_inner());
                                                                                let string1: String = tuple.0;
                                                                                let string2: String = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_string(string1, string2);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::texture_param => {
                                                                                let tuple: (String, String) = pbrt_texture_parameter(&mut parameter_pair.into_inner());
                                                                                let string1: String = tuple.0;
                                                                                let string2: String = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    param_set.add_texture(string1, string2);
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            Rule::vector_param => {
                                                                                let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                let string: String = tuple.0;
                                                                                let floats: Vec<Float> = tuple.1;
                                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                                    if floats.len() == 3 {
                                                                                        param_set.add_vector3f(string,
                                                                                                               Vector3f {
                                                                                                                   x: floats[0],
                                                                                                                   y: floats[1],
                                                                                                                   z: floats[2],
                                                                                                               });
                                                                                    } else {
                                                                                        param_set.add_point3fs(string, floats);
                                                                                    }
                                                                                } else {
                                                                                    panic!("Can't get parameter set.");
                                                                                }
                                                                            }
                                                                            _ => unreachable!(),
                                                                        };
                                                                    }
                                                                }
                                                                _ => unreachable!(),
                                                            }
                                                        }
                                                        // we should have the texture parameters by now
                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                            // println!("Texture \"{}\" \"{}\" \"{}\" ",
                                                            //          param_set.name,
                                                            //          param_set.tex_type,
                                                            //          param_set.tex_name);
                                                            // print_params(&param_set);
                                                            make_texture(&param_set);
                                                        } else {
                                                            panic!("Can't get parameter set.");
                                                        }
                                                    }
                                                    _ => unreachable!(),
                                                }
                                            }
                                        }
                                        Rule::rotate => {
                                            let mut numbers: Vec<Float> = Vec::new();
                                            for rotate_pair in statement_pair.into_inner() {
                                                let number: Float =
                                                    f32::from_str(rotate_pair.clone().into_span().as_str()).unwrap();
                                                numbers.push(number);
                                            }
                                            let angle: Float = numbers[0];
                                            let x: Float = numbers[1];
                                            let y: Float = numbers[2];
                                            let z: Float = numbers[3];
                                            let rotate: Transform = Transform::rotate(angle, Vector3f { x: x, y: y, z: z, });
                                            if ACTIVE_TRANSFORM_BITS & 1_u8 > 0_u8 {
                                                // 0x?1
                                                CUR_TRANSFORM.t[0] = CUR_TRANSFORM.t[0] * rotate;
                                            }
                                            if ACTIVE_TRANSFORM_BITS & 2_u8 > 0_u8 {
                                                // 0x1?
                                                CUR_TRANSFORM.t[1] = CUR_TRANSFORM.t[1] * rotate;
                                            }
                                        }
                                        Rule::scale => {
                                            let mut numbers: Vec<Float> = Vec::new();
                                            for scale_pair in statement_pair.into_inner() {
                                                let number: Float =
                                                    f32::from_str(scale_pair.clone().into_span().as_str()).unwrap();
                                                numbers.push(number);
                                            }
                                            let x: Float = numbers[0];
                                            let y: Float = numbers[1];
                                            let z: Float = numbers[2];
                                            let scale: Transform = Transform::scale(x, y, z);
                                            if ACTIVE_TRANSFORM_BITS & 1_u8 > 0_u8 {
                                                // 0x?1
                                                CUR_TRANSFORM.t[0] = CUR_TRANSFORM.t[0] * scale;
                                            }
                                            if ACTIVE_TRANSFORM_BITS & 2_u8 > 0_u8 {
                                                // 0x1?
                                                CUR_TRANSFORM.t[1] = CUR_TRANSFORM.t[1] * scale;
                                            }
                                        }
                                        Rule::transform => {
                                            let mut numbers: Vec<Float> = Vec::new();
                                            for transform_pair in statement_pair.into_inner() {
                                                // ignore brackets
                                                let not_opening: bool = transform_pair.as_str() != String::from("[");
                                                let not_closing: bool = transform_pair.as_str() != String::from("]");
                                                if not_opening && not_closing {
                                                    let number: Float =
                                                        f32::from_str(transform_pair.clone().into_span().as_str()).unwrap();
                                                    numbers.push(number);
                                                }
                                            }
                                            let m00: Float = numbers[0];
                                            let m01: Float = numbers[1];
                                            let m02: Float = numbers[2];
                                            let m03: Float = numbers[3];
                                            let m10: Float = numbers[4];
                                            let m11: Float = numbers[5];
                                            let m12: Float = numbers[6];
                                            let m13: Float = numbers[7];
                                            let m20: Float = numbers[8];
                                            let m21: Float = numbers[9];
                                            let m22: Float = numbers[10];
                                            let m23: Float = numbers[11];
                                            let m30: Float = numbers[12];
                                            let m31: Float = numbers[13];
                                            let m32: Float = numbers[14];
                                            let m33: Float = numbers[15];
                                            let transform: Transform = Transform::new(m00, m10, m20, m30,
                                                                                      m01, m11, m21, m31,
                                                                                      m02, m12, m22, m32,
                                                                                      m03, m13, m23, m33);
                                            if ACTIVE_TRANSFORM_BITS & 1_u8 > 0_u8 {
                                                // 0x?1
                                                CUR_TRANSFORM.t[0] = CUR_TRANSFORM.t[0] * transform;
                                            }
                                            if ACTIVE_TRANSFORM_BITS & 2_u8 > 0_u8 {
                                                // 0x1?
                                                CUR_TRANSFORM.t[1] = CUR_TRANSFORM.t[1] * transform;
                                            }
                                        }
                                        Rule::transform_times => {
                                            let mut numbers: Vec<Float> = Vec::new();
                                            for transform_times_pair in
                                                statement_pair.into_inner() {
                                                let number: Float =
                                                    f32::from_str(transform_times_pair.clone().into_span().as_str()).unwrap();
                                                numbers.push(number);
                                            }
                                            let start: Float = numbers[0];
                                            let end: Float = numbers[1];
                                            if let Some(ref mut ro) = RENDER_OPTIONS {
                                                // TODO: VERIFY_OPTIONS("TransformTimes");
                                                ro.transform_start_time = start;
                                                ro.transform_end_time = end;
                                                println!("TransformTimes {} {}",
                                                         ro.transform_start_time,
                                                         ro.transform_end_time);
                                            }
                                        }
                                        Rule::translate => {
                                            let mut numbers: Vec<Float> = Vec::new();
                                            for translate_pair in statement_pair.into_inner() {
                                                let number: Float =
                                                    f32::from_str(translate_pair.clone().into_span().as_str()).unwrap();
                                                numbers.push(number);
                                            }
                                            let x: Float = numbers[0];
                                            let y: Float = numbers[1];
                                            let z: Float = numbers[2];
                                            let translate: Transform = Transform::translate(Vector3f { x: x, y: y, z: z, });
                                            if ACTIVE_TRANSFORM_BITS & 1_u8 > 0_u8 {
                                                // 0x?1
                                                CUR_TRANSFORM.t[0] = CUR_TRANSFORM.t[0] * translate;
                                            }
                                            if ACTIVE_TRANSFORM_BITS & 2_u8 > 0_u8 {
                                                // 0x1?
                                                CUR_TRANSFORM.t[1] = CUR_TRANSFORM.t[1] * translate;
                                            }
                                        }
                                        _ => unreachable!(),
                                    };
                                }
                            }
                            Rule::last_statement => {
                                println!("WorldEnd");
                                pbrt_world_end();
                            }
                            _ => unreachable!(),
                        }
                    }
                    // parser.main();
                    println!("done.");
                }
            }
            None => panic!("No input file name."),
        }
        return;
    } else if matches.opt_present("v") {
        print_version(&program);
        return;
    } else {
        print_usage(&program, opts);
        return;
    }
}
