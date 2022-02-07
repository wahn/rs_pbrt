// Assumptions:
// 1. objects, data, and materials have the same base name:
//    e.g. OBcornellbox, MEcornellbox and MAcornellbox
//    TODO: Find a better solution (following the Object->data pointer?)
// 2. Right now we search for "Camera" for Transform::look_at(...)
//    TODO: get the render camera/transform from the scene (self.scene.camera)
// 3. Smoothness is valid for the whole mesh
//    TODO: store smoothness per polygon, split mesh into smooth and non-smooth parts
// 4. There might be one IM block with a HDR image, which is used for IBL
//    TODO: check if World actually uses "Environment Lighting"
// 5. Lights with type LA_SUN have a position but point always to the origin
//    TODO: use an empty (compass) object as target

// command line options
use structopt::StructOpt;
// std
use std::collections::HashMap;
use std::convert::TryInto;
use std::ffi::OsString;
use std::path::Path;
use std::sync::Arc;
// others
use blend_info::{
    calc_mem_tlen, get_char, get_float, get_float2, get_float3, get_id_name, get_int, get_matrix,
    get_pointer, get_short, get_short3, read_dna, use_dna, DnaStrC,
};
// pbrt
use rs_pbrt::core::api::{make_accelerator, make_camera, make_film, make_filter, make_sampler};
use rs_pbrt::core::camera::Camera;
use rs_pbrt::core::film::Film;
use rs_pbrt::core::geometry::{Bounds2f, Bounds2i, Normal3f, Point2f, Point2i, Point3f, Vector3f};
use rs_pbrt::core::integrator::{Integrator, SamplerIntegrator};
use rs_pbrt::core::light::Light;
use rs_pbrt::core::material::Material;
use rs_pbrt::core::medium::MediumInterface;
use rs_pbrt::core::mipmap::ImageWrap;
use rs_pbrt::core::paramset::ParamSet;
use rs_pbrt::core::pbrt::degrees;
use rs_pbrt::core::pbrt::{Float, Spectrum};
use rs_pbrt::core::primitive::{GeometricPrimitive, Primitive};
use rs_pbrt::core::sampler::Sampler;
use rs_pbrt::core::scene::Scene;
use rs_pbrt::core::shape::Shape;
use rs_pbrt::core::texture::{Texture, TextureMapping2D, UVMapping2D};
use rs_pbrt::core::transform::{AnimatedTransform, Transform};
use rs_pbrt::integrators::ao::AOIntegrator;
use rs_pbrt::integrators::bdpt::BDPTIntegrator;
use rs_pbrt::integrators::directlighting::{DirectLightingIntegrator, LightStrategy};
use rs_pbrt::integrators::mlt::MLTIntegrator;
use rs_pbrt::integrators::path::PathIntegrator;
use rs_pbrt::integrators::sppm::SPPMIntegrator;
use rs_pbrt::integrators::volpath::VolPathIntegrator;
use rs_pbrt::integrators::whitted::WhittedIntegrator;
use rs_pbrt::lights::diffuse::DiffuseAreaLight;
use rs_pbrt::lights::distant::DistantLight;
use rs_pbrt::lights::infinite::InfiniteAreaLight;
use rs_pbrt::lights::point::PointLight;
use rs_pbrt::materials::glass::GlassMaterial;
use rs_pbrt::materials::matte::MatteMaterial;
use rs_pbrt::materials::metal::MetalMaterial;
use rs_pbrt::materials::metal::{COPPER_K, COPPER_N, COPPER_SAMPLES, COPPER_WAVELENGTHS};
use rs_pbrt::materials::mirror::MirrorMaterial;
// use rs_pbrt::shapes::cylinder::Cylinder;
// use rs_pbrt::shapes::disk::Disk;
// use rs_pbrt::shapes::sphere::Sphere;
use rs_pbrt::shapes::triangle::{Triangle, TriangleMesh};
use rs_pbrt::textures::constant::ConstantTexture;
use rs_pbrt::textures::imagemap::convert_to_spectrum;
use rs_pbrt::textures::imagemap::ImageTexture;

pub const VERSION: &'static str = env!("CARGO_PKG_VERSION");

/// Parse a Blender scene file and render it.
#[derive(StructOpt)]
struct Cli {
    /// camera name
    #[structopt(short = "c", long = "camera_name")]
    camera_name: Option<String>,
    /// global light scaling
    #[structopt(short = "l", long = "light_scale", default_value = "1.0")]
    light_scale: f32,
    /// pixel samples
    #[structopt(short = "s", long = "samples", default_value = "1")]
    samples: u32,
    /// ao, directlighting, whitted, path, bdpt, mlt, sppm, volpath
    #[structopt(short = "i", long = "integrator")]
    integrator: Option<String>,
    /// max length of a light-carrying path
    #[structopt(short = "m", long = "max_depth", default_value = "5")]
    max_depth: u32,
    /// bootstrap samples [MLT]
    #[structopt(long = "bootstrap_samples", default_value = "100000")]
    bootstrap_samples: u32,
    /// number of Markov chains [MLT]
    #[structopt(long = "chains", default_value = "1000")]
    chains: u32,
    /// number of path mutations [MLT]
    #[structopt(long = "mutations_per_pixel", default_value = "100")]
    mutations_per_pixel: u32,
    /// prob of discarding path [MLT]
    #[structopt(long = "step_probability", default_value = "0.3")]
    step_probability: f32,
    /// perturbation deviation [MLT]
    #[structopt(long = "sigma", default_value = "0.01")]
    sigma: f32,
    /// frequency to write image [SPPM]
    #[structopt(long = "write_frequency", default_value = "1")]
    write_frequency: i32,
    /// The path to the file to read
    #[structopt(parse(from_os_str))]
    path: std::path::PathBuf,
}

// PBRT

// #[derive(Debug, Default, Copy, Clone)]
// struct PbrtSphere {
//     pub radius: f32,
//     pub zmin: f32,
//     pub zmax: f32,
//     pub phimax: f32,
// }

// impl PbrtSphere {
//     fn new(radius: f32, zmin: f32, zmax: f32, phimax: f32) -> Self {
//         PbrtSphere {
//             radius,
//             zmin,
//             zmax,
//             phimax,
//         }
//     }
// }

// #[derive(Debug, Default, Copy, Clone)]
// struct PbrtCylinder {
//     pub radius: f32,
//     pub zmin: f32,
//     pub zmax: f32,
//     pub phimax: f32,
// }

// impl PbrtCylinder {
//     fn new(radius: f32, zmin: f32, zmax: f32, phimax: f32) -> Self {
//         PbrtCylinder {
//             radius,
//             zmin,
//             zmax,
//             phimax,
//         }
//     }
// }

// #[derive(Debug, Default, Copy, Clone)]
// struct PbrtDisk {
//     pub height: f32,
//     pub radius: f32,
//     pub innerradius: f32,
//     pub phimax: f32,
// }

// impl PbrtDisk {
//     fn new(height: f32, radius: f32, innerradius: f32, phimax: f32) -> Self {
//         PbrtDisk {
//             height,
//             radius,
//             innerradius,
//             phimax,
//         }
//     }
// }

// Blender

#[derive(Debug, Default, Copy, Clone)]
struct BlendCamera {
    pub lens: f32,
    // pub angle_x: f32,
    // pub angle_y: f32,
    pub clipsta: f32,
}

#[derive(Debug, Default, Copy, Clone)]
struct Blend279Material {
    pub r: f32,
    pub g: f32,
    pub b: f32,
    // pub a: f32,
    pub specr: f32,
    pub specg: f32,
    pub specb: f32,
    pub mirr: f32,
    pub mirg: f32,
    pub mirb: f32,
    pub emit: f32,
    pub ang: f32, // IOR
    pub ray_mirror: f32,
    pub roughness: f32,
}

fn focallength_to_fov(focal_length: f32, sensor: f32) -> f32 {
    2.0_f32 * ((sensor / 2.0_f32) / focal_length).atan()
}

// TMP (see pbrt_spheres_differentials_texfilt.rs)

struct SceneDescription {
    // mesh
    mesh_names: Vec<String>,
    meshes: Vec<Arc<TriangleMesh>>,
    triangle_colors: Vec<Vec<Spectrum>>,
    // cylinder
    cylinder_names: Vec<String>,
    cylinders: Vec<Arc<Shape>>,
    // disk
    disk_names: Vec<String>,
    disks: Vec<Arc<Shape>>,
    // sphere
    sphere_names: Vec<String>,
    spheres: Vec<Arc<Shape>>,
    // lights
    lights: Vec<Arc<Light>>,
}

struct SceneDescriptionBuilder {
    // mesh
    mesh_names: Vec<String>,
    meshes: Vec<Arc<TriangleMesh>>,
    triangle_colors: Vec<Vec<Spectrum>>,
    // cylinder
    cylinder_names: Vec<String>,
    cylinders: Vec<Arc<Shape>>,
    // disk
    disk_names: Vec<String>,
    disks: Vec<Arc<Shape>>,
    // sphere
    sphere_names: Vec<String>,
    spheres: Vec<Arc<Shape>>,
    // lights
    lights: Vec<Arc<Light>>,
}

impl SceneDescriptionBuilder {
    fn new() -> SceneDescriptionBuilder {
        SceneDescriptionBuilder {
            mesh_names: Vec::with_capacity(256),
            meshes: Vec::with_capacity(256),
            triangle_colors: Vec::with_capacity(256),
            cylinder_names: Vec::new(),
            cylinders: Vec::new(),
            disk_names: Vec::new(),
            disks: Vec::new(),
            sphere_names: Vec::new(),
            spheres: Vec::new(),
            lights: Vec::new(),
        }
    }
    fn add_mesh(
        &mut self,
        base_name: String,
        object_to_world: Transform,
        world_to_object: Transform,
        n_triangles: u32,
        vertex_indices: Vec<u32>,
        n_vertices: u32,
        p_ws: Vec<Point3f>,
        s: Vec<Vector3f>,
        n_ws: Vec<Normal3f>,
        uv: Vec<Point2f>,
        triangle_colors: Vec<Spectrum>,
    ) -> &mut SceneDescriptionBuilder {
        self.mesh_names.push(base_name);
        let triangle_mesh = Arc::new(TriangleMesh::new(
            object_to_world,
            world_to_object,
            false,
            n_triangles.try_into().unwrap(),
            vertex_indices.try_into().unwrap(),
            n_vertices,
            p_ws, // in world space
            s,    // empty
            n_ws, // in world space
            uv,
            None,
            None,
        ));
        self.meshes.push(triangle_mesh);
        self.triangle_colors.push(triangle_colors);
        self
    }
    // fn add_cylinder(
    //     &mut self,
    //     base_name: String,
    //     object_to_world: Transform,
    //     world_to_object: Transform,
    //     radius: Float,
    //     z_min: Float,
    //     z_max: Float,
    //     phi_max: Float,
    // ) -> &mut SceneDescriptionBuilder {
    //     self.cylinder_names.push(base_name);
    //     let cylinder = Arc::new(Shape::Clndr(Cylinder::new(
    //         object_to_world,
    //         world_to_object,
    //         false,
    //         radius,
    //         z_min,
    //         z_max,
    //         phi_max,
    //     )));
    //     self.cylinders.push(cylinder);
    //     self
    // }
    // fn add_disk(
    //     &mut self,
    //     base_name: String,
    //     object_to_world: Transform,
    //     world_to_object: Transform,
    //     height: Float,
    //     radius: Float,
    //     inner_radius: Float,
    //     phi_max: Float,
    // ) -> &mut SceneDescriptionBuilder {
    //     self.disk_names.push(base_name);
    //     let disk = Arc::new(Shape::Dsk(Disk::new(
    //         object_to_world,
    //         world_to_object,
    //         false,
    //         height,
    //         radius,
    //         inner_radius,
    //         phi_max,
    //     )));
    //     self.disks.push(disk);
    //     self
    // }
    // fn add_sphere(
    //     &mut self,
    //     base_name: String,
    //     object_to_world: Transform,
    //     world_to_object: Transform,
    //     radius: Float,
    //     z_min: Float,
    //     z_max: Float,
    //     phi_max: Float,
    // ) -> &mut SceneDescriptionBuilder {
    //     self.sphere_names.push(base_name);
    //     let sphere = Arc::new(Shape::Sphr(Sphere::new(
    //         object_to_world,
    //         world_to_object,
    //         false,
    //         radius,
    //         z_min,
    //         z_max,
    //         phi_max,
    //     )));
    //     self.spheres.push(sphere);
    //     self
    // }
    fn add_hdr_light(
        &mut self,
        light_to_world: Transform,
        texmap: String,
        light_scale: f32,
    ) -> &mut SceneDescriptionBuilder {
        let l: Spectrum = Spectrum::new(1.0 as Float);
        let sc: Spectrum = Spectrum::new(light_scale as Float);
        let n_samples: i32 = 1;
        let infinte_light = Arc::new(Light::InfiniteArea(Box::new(InfiniteAreaLight::new(
            &light_to_world,
            &(l * sc),
            n_samples,
            texmap,
        ))));
        self.lights.push(infinte_light);
        self
    }
    fn add_distant_light(
        &mut self,
        light_to_world: Transform,
        l: Spectrum,
        light_scale: f32,
    ) -> &mut SceneDescriptionBuilder {
        let sc: Spectrum = Spectrum::new(light_scale as Float);
        let mut from: Point3f = Point3f {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };
        let to: Point3f = Point3f {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };
        from = light_to_world.transform_point(&from);
        let dir: Vector3f = from - to;
        let object_to_world: Transform = Transform::new(
            1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
        );
        let distant_light = Arc::new(Light::Distant(Box::new(DistantLight::new(
            &object_to_world,
            &(l * sc),
            &dir,
        ))));
        self.lights.push(distant_light);
        self
    }
    fn add_point_light(
        &mut self,
        light_to_world: Transform,
        l: Spectrum,
        light_scale: f32,
    ) -> &mut SceneDescriptionBuilder {
        let sc: Spectrum = Spectrum::new(light_scale as Float);
        let medium_interface: MediumInterface = MediumInterface::default();
        let point_light = Arc::new(Light::Point(Box::new(PointLight::new(
            &light_to_world,
            &medium_interface,
            &(l * sc),
        ))));
        self.lights.push(point_light);
        self
    }
    fn finalize(self) -> SceneDescription {
        SceneDescription {
            mesh_names: self.mesh_names,
            meshes: self.meshes,
            triangle_colors: self.triangle_colors,
            cylinder_names: self.cylinder_names,
            cylinders: self.cylinders,
            disk_names: self.disk_names,
            disks: self.disks,
            sphere_names: self.sphere_names,
            spheres: self.spheres,
            lights: self.lights,
        }
    }
}

struct RenderOptions {
    has_emitters: bool,
    primitives: Vec<Arc<Primitive>>,
    shapes: Vec<Arc<Shape>>,
    shape_materials: Vec<Arc<Material>>,
    shape_lights: Vec<Option<Arc<Light>>>,
    lights: Vec<Arc<Light>>,
}

impl RenderOptions {
    fn new(
        scene: SceneDescription,
        material_hm: &HashMap<String, Blend279Material>,
        texture_hm: &HashMap<String, OsString>,
        light_scale: Float,
    ) -> RenderOptions {
        let mut has_emitters: bool = false;
        let primitives: Vec<Arc<Primitive>> = Vec::new();
        let mut shapes: Vec<Arc<Shape>> = Vec::new();
        let mut shape_materials: Vec<Arc<Material>> = Vec::new();
        let mut shape_lights: Vec<Option<Arc<Light>>> = Vec::new();
        let mut lights: Vec<Arc<Light>> = Vec::new();
        // default material
        let kd = Arc::new(ConstantTexture::new(Spectrum::new(1.0)));
        let sigma = Arc::new(ConstantTexture::new(0.0 as Float));
        let default_material = Arc::new(Material::Matte(Box::new(MatteMaterial::new(
            kd, sigma, None,
        ))));
        // lights
        for light in &scene.lights {
            lights.push(light.clone());
        }
        // cylinders
        for cylinder_idx in 0..scene.cylinders.len() {
            let cylinder = &scene.cylinders[cylinder_idx];
            let cylinder_name = &scene.cylinder_names[cylinder_idx];
            if let Some(mat) = get_material(cylinder_name, material_hm) {
                // println!("{:?}: {:?}", cylinder_name, mat);
                if mat.emit > 0.0 {
                    has_emitters = true;
                    let mi: MediumInterface = MediumInterface::default();
                    let l_emit: Spectrum = Spectrum::rgb(
                        mat.r * mat.emit * light_scale,
                        mat.g * mat.emit * light_scale,
                        mat.b * mat.emit * light_scale,
                    );
                    let n_samples: i32 = 1;
                    let two_sided: bool = false;
                    let area_light: Arc<Light> =
                        Arc::new(Light::DiffuseArea(Box::new(DiffuseAreaLight::new(
                            &cylinder.get_object_to_world(),
                            &mi,
                            &l_emit,
                            n_samples,
                            cylinder.clone(),
                            two_sided,
                        ))));
                    lights.push(area_light.clone());
                    shapes.push(cylinder.clone());
                    shape_materials.push(default_material.clone());
                    let shape_light: Option<Arc<Light>> = Some(area_light.clone());
                    shape_lights.push(shape_light);
                } else {
                    if mat.ang != 1.0 {
                        // GlassMaterial
                        let kr = Arc::new(ConstantTexture::new(Spectrum::new(1.0)));
                        let kt = Arc::new(ConstantTexture::new(Spectrum::rgb(
                            mat.specr, mat.specg, mat.specb,
                        )));
                        let u_roughness = Arc::new(ConstantTexture::new(0.0 as Float));
                        let v_roughness = Arc::new(ConstantTexture::new(0.0 as Float));
                        let index = Arc::new(ConstantTexture::new(mat.ang as Float));
                        let glass = Arc::new(Material::Glass(Box::new(GlassMaterial {
                            kr: kr,
                            kt: kt,
                            u_roughness: u_roughness,
                            v_roughness: v_roughness,
                            index: index,
                            bump_map: None,
                            remap_roughness: true,
                        })));
                        shapes.push(cylinder.clone());
                        shape_materials.push(glass.clone());
                        shape_lights.push(None);
                    } else if mat.ray_mirror > 0.0 {
                        if mat.roughness > 0.0 {
                            // MetalMaterial
                            let copper_n: Spectrum = Spectrum::from_sampled(
                                &COPPER_WAVELENGTHS,
                                &COPPER_N,
                                COPPER_SAMPLES as i32,
                            );
                            let eta: Arc<dyn Texture<Spectrum> + Send + Sync> =
                                Arc::new(ConstantTexture::new(copper_n));
                            let copper_k: Spectrum = Spectrum::from_sampled(
                                &COPPER_WAVELENGTHS,
                                &COPPER_K,
                                COPPER_SAMPLES as i32,
                            );
                            let k: Arc<dyn Texture<Spectrum> + Send + Sync> =
                                Arc::new(ConstantTexture::new(copper_k));
                            let remap_roughness: bool = true;
                            let metal = Arc::new(Material::Metal(Box::new(MetalMaterial::new(
                                eta,
                                k,
                                Arc::new(ConstantTexture::new(mat.roughness as Float)),
                                None,
                                None,
                                None,
                                remap_roughness,
                            ))));
                            shapes.push(cylinder.clone());
                            shape_materials.push(metal.clone());
                            shape_lights.push(None);
                        } else {
                            // MirrorMaterial
                            let kr = Arc::new(ConstantTexture::new(Spectrum::rgb(
                                mat.mirr * mat.ray_mirror,
                                mat.mirg * mat.ray_mirror,
                                mat.mirb * mat.ray_mirror,
                            )));
                            let mirror =
                                Arc::new(Material::Mirror(Box::new(MirrorMaterial::new(kr, None))));
                            shapes.push(cylinder.clone());
                            shape_materials.push(mirror.clone());
                            shape_lights.push(None);
                        }
                    } else {
                        // MatteMaterial
                        let mut kd: Arc<dyn Texture<Spectrum> + Send + Sync> =
                            Arc::new(ConstantTexture::new(Spectrum::rgb(mat.r, mat.g, mat.b)));
                        if let Some(tex) = texture_hm.get(cylinder_name) {
                            // first try texture with exactly the same name as the mesh
                            let su: Float = 1.0;
                            let sv: Float = 1.0;
                            let du: Float = 0.0;
                            let dv: Float = 0.0;
                            let mapping: Box<TextureMapping2D> =
                                Box::new(TextureMapping2D::UV(UVMapping2D {
                                    su: su,
                                    sv: sv,
                                    du: du,
                                    dv: dv,
                                }));
                            let filename: String = String::from(tex.to_str().unwrap());
                            let do_trilinear: bool = false;
                            let max_aniso: Float = 8.0;
                            let wrap_mode: ImageWrap = ImageWrap::Repeat;
                            let scale: Float = 1.0;
                            let gamma: bool = false;
                            kd = Arc::new(ImageTexture::new(
                                mapping,
                                filename,
                                do_trilinear,
                                max_aniso,
                                wrap_mode,
                                scale,
                                gamma,
                                convert_to_spectrum,
                            ));
                        } else {
                            // then remove trailing digits from mesh name
                            let mut ntd: String = String::new();
                            let mut chars = cylinder_name.chars();
                            let mut digits: String = String::new(); // many digits
                            while let Some(c) = chars.next() {
                                if c.is_digit(10_u32) {
                                    // collect digits
                                    digits.push(c);
                                } else {
                                    // push collected digits (if any)
                                    ntd += &digits;
                                    // and reset
                                    digits = String::new();
                                    // push non-digit
                                    ntd.push(c);
                                }
                            }
                            // try no trailing digits (ntd)
                            if let Some(tex) = texture_hm.get(&ntd) {
                                let su: Float = 1.0;
                                let sv: Float = 1.0;
                                let du: Float = 0.0;
                                let dv: Float = 0.0;
                                let mapping: Box<TextureMapping2D> =
                                    Box::new(TextureMapping2D::UV(UVMapping2D {
                                        su: su,
                                        sv: sv,
                                        du: du,
                                        dv: dv,
                                    }));
                                let filename: String = String::from(tex.to_str().unwrap());
                                let do_trilinear: bool = false;
                                let max_aniso: Float = 8.0;
                                let wrap_mode: ImageWrap = ImageWrap::Repeat;
                                let scale: Float = 1.0;
                                let gamma: bool = false;
                                kd = Arc::new(ImageTexture::new(
                                    mapping,
                                    filename,
                                    do_trilinear,
                                    max_aniso,
                                    wrap_mode,
                                    scale,
                                    gamma,
                                    convert_to_spectrum,
                                ));
                            }
                        }
                        let sigma = Arc::new(ConstantTexture::new(0.0 as Float));
                        let matte = Arc::new(Material::Matte(Box::new(MatteMaterial::new(
                            kd,
                            sigma.clone(),
                            None,
                        ))));
                        shapes.push(cylinder.clone());
                        shape_materials.push(matte.clone());
                        shape_lights.push(None);
                    }
                }
            }
        }
        // disks
        for disk_idx in 0..scene.disks.len() {
            let disk = &scene.disks[disk_idx];
            let disk_name = &scene.disk_names[disk_idx];
            if let Some(mat) = get_material(disk_name, material_hm) {
                // println!("{:?}: {:?}", disk_name, mat);
                if mat.emit > 0.0 {
                    has_emitters = true;
                    let mi: MediumInterface = MediumInterface::default();
                    let l_emit: Spectrum = Spectrum::rgb(
                        mat.r * mat.emit * light_scale,
                        mat.g * mat.emit * light_scale,
                        mat.b * mat.emit * light_scale,
                    );
                    let n_samples: i32 = 1;
                    let two_sided: bool = false;
                    let area_light: Arc<Light> =
                        Arc::new(Light::DiffuseArea(Box::new(DiffuseAreaLight::new(
                            &disk.get_object_to_world(),
                            &mi,
                            &l_emit,
                            n_samples,
                            disk.clone(),
                            two_sided,
                        ))));
                    lights.push(area_light.clone());
                    shapes.push(disk.clone());
                    shape_materials.push(default_material.clone());
                    let shape_light: Option<Arc<Light>> = Some(area_light.clone());
                    shape_lights.push(shape_light);
                } else {
                    if mat.ang != 1.0 {
                        // GlassMaterial
                        let kr = Arc::new(ConstantTexture::new(Spectrum::new(1.0)));
                        let kt = Arc::new(ConstantTexture::new(Spectrum::rgb(
                            mat.specr, mat.specg, mat.specb,
                        )));
                        let u_roughness = Arc::new(ConstantTexture::new(0.0 as Float));
                        let v_roughness = Arc::new(ConstantTexture::new(0.0 as Float));
                        let index = Arc::new(ConstantTexture::new(mat.ang as Float));
                        let glass = Arc::new(Material::Glass(Box::new(GlassMaterial {
                            kr: kr,
                            kt: kt,
                            u_roughness: u_roughness,
                            v_roughness: v_roughness,
                            index: index,
                            bump_map: None,
                            remap_roughness: true,
                        })));
                        shapes.push(disk.clone());
                        shape_materials.push(glass.clone());
                        shape_lights.push(None);
                    } else if mat.ray_mirror > 0.0 {
                        if mat.roughness > 0.0 {
                            // MetalMaterial
                            let copper_n: Spectrum = Spectrum::from_sampled(
                                &COPPER_WAVELENGTHS,
                                &COPPER_N,
                                COPPER_SAMPLES as i32,
                            );
                            let eta: Arc<dyn Texture<Spectrum> + Send + Sync> =
                                Arc::new(ConstantTexture::new(copper_n));
                            let copper_k: Spectrum = Spectrum::from_sampled(
                                &COPPER_WAVELENGTHS,
                                &COPPER_K,
                                COPPER_SAMPLES as i32,
                            );
                            let k: Arc<dyn Texture<Spectrum> + Send + Sync> =
                                Arc::new(ConstantTexture::new(copper_k));
                            let remap_roughness: bool = true;
                            let metal = Arc::new(Material::Metal(Box::new(MetalMaterial::new(
                                eta,
                                k,
                                Arc::new(ConstantTexture::new(mat.roughness as Float)),
                                None,
                                None,
                                None,
                                remap_roughness,
                            ))));
                            shapes.push(disk.clone());
                            shape_materials.push(metal.clone());
                            shape_lights.push(None);
                        } else {
                            // MirrorMaterial
                            let kr = Arc::new(ConstantTexture::new(Spectrum::rgb(
                                mat.mirr * mat.ray_mirror,
                                mat.mirg * mat.ray_mirror,
                                mat.mirb * mat.ray_mirror,
                            )));
                            let mirror =
                                Arc::new(Material::Mirror(Box::new(MirrorMaterial::new(kr, None))));
                            shapes.push(disk.clone());
                            shape_materials.push(mirror.clone());
                            shape_lights.push(None);
                        }
                    } else {
                        // MatteMaterial
                        let mut kd: Arc<dyn Texture<Spectrum> + Send + Sync> =
                            Arc::new(ConstantTexture::new(Spectrum::rgb(mat.r, mat.g, mat.b)));
                        if let Some(tex) = texture_hm.get(disk_name) {
                            // first try texture with exactly the same name as the mesh
                            let su: Float = 1.0;
                            let sv: Float = 1.0;
                            let du: Float = 0.0;
                            let dv: Float = 0.0;
                            let mapping: Box<TextureMapping2D> =
                                Box::new(TextureMapping2D::UV(UVMapping2D {
                                    su: su,
                                    sv: sv,
                                    du: du,
                                    dv: dv,
                                }));
                            let filename: String = String::from(tex.to_str().unwrap());
                            let do_trilinear: bool = false;
                            let max_aniso: Float = 8.0;
                            let wrap_mode: ImageWrap = ImageWrap::Repeat;
                            let scale: Float = 1.0;
                            let gamma: bool = false;
                            kd = Arc::new(ImageTexture::new(
                                mapping,
                                filename,
                                do_trilinear,
                                max_aniso,
                                wrap_mode,
                                scale,
                                gamma,
                                convert_to_spectrum,
                            ));
                        } else {
                            // then remove trailing digits from mesh name
                            let mut ntd: String = String::new();
                            let mut chars = disk_name.chars();
                            let mut digits: String = String::new(); // many digits
                            while let Some(c) = chars.next() {
                                if c.is_digit(10_u32) {
                                    // collect digits
                                    digits.push(c);
                                } else {
                                    // push collected digits (if any)
                                    ntd += &digits;
                                    // and reset
                                    digits = String::new();
                                    // push non-digit
                                    ntd.push(c);
                                }
                            }
                            // try no trailing digits (ntd)
                            if let Some(tex) = texture_hm.get(&ntd) {
                                let su: Float = 1.0;
                                let sv: Float = 1.0;
                                let du: Float = 0.0;
                                let dv: Float = 0.0;
                                let mapping: Box<TextureMapping2D> =
                                    Box::new(TextureMapping2D::UV(UVMapping2D {
                                        su: su,
                                        sv: sv,
                                        du: du,
                                        dv: dv,
                                    }));
                                let filename: String = String::from(tex.to_str().unwrap());
                                let do_trilinear: bool = false;
                                let max_aniso: Float = 8.0;
                                let wrap_mode: ImageWrap = ImageWrap::Repeat;
                                let scale: Float = 1.0;
                                let gamma: bool = false;
                                kd = Arc::new(ImageTexture::new(
                                    mapping,
                                    filename,
                                    do_trilinear,
                                    max_aniso,
                                    wrap_mode,
                                    scale,
                                    gamma,
                                    convert_to_spectrum,
                                ));
                            }
                        }
                        let sigma = Arc::new(ConstantTexture::new(0.0 as Float));
                        let matte = Arc::new(Material::Matte(Box::new(MatteMaterial::new(
                            kd,
                            sigma.clone(),
                            None,
                        ))));
                        shapes.push(disk.clone());
                        shape_materials.push(matte.clone());
                        shape_lights.push(None);
                    }
                }
            }
        }
        // spheres
        for sphere_idx in 0..scene.spheres.len() {
            let sphere = &scene.spheres[sphere_idx];
            let sphere_name = &scene.sphere_names[sphere_idx];
            if let Some(mat) = get_material(sphere_name, material_hm) {
                // println!("{:?}: {:?}", sphere_name, mat);
                if mat.emit > 0.0 {
                    has_emitters = true;
                    let mi: MediumInterface = MediumInterface::default();
                    let l_emit: Spectrum = Spectrum::rgb(
                        mat.r * mat.emit * light_scale,
                        mat.g * mat.emit * light_scale,
                        mat.b * mat.emit * light_scale,
                    );
                    let n_samples: i32 = 1;
                    let two_sided: bool = false;
                    let area_light: Arc<Light> =
                        Arc::new(Light::DiffuseArea(Box::new(DiffuseAreaLight::new(
                            &sphere.get_object_to_world(),
                            &mi,
                            &l_emit,
                            n_samples,
                            sphere.clone(),
                            two_sided,
                        ))));
                    lights.push(area_light.clone());
                    shapes.push(sphere.clone());
                    shape_materials.push(default_material.clone());
                    let shape_light: Option<Arc<Light>> = Some(area_light.clone());
                    shape_lights.push(shape_light);
                } else {
                    if mat.ang != 1.0 {
                        // GlassMaterial
                        let kr = Arc::new(ConstantTexture::new(Spectrum::new(1.0)));
                        let kt = Arc::new(ConstantTexture::new(Spectrum::rgb(
                            mat.specr, mat.specg, mat.specb,
                        )));
                        let u_roughness = Arc::new(ConstantTexture::new(0.0 as Float));
                        let v_roughness = Arc::new(ConstantTexture::new(0.0 as Float));
                        let index = Arc::new(ConstantTexture::new(mat.ang as Float));
                        let glass = Arc::new(Material::Glass(Box::new(GlassMaterial {
                            kr: kr,
                            kt: kt,
                            u_roughness: u_roughness,
                            v_roughness: v_roughness,
                            index: index,
                            bump_map: None,
                            remap_roughness: true,
                        })));
                        shapes.push(sphere.clone());
                        shape_materials.push(glass.clone());
                        shape_lights.push(None);
                    } else if mat.ray_mirror > 0.0 {
                        if mat.roughness > 0.0 {
                            // MetalMaterial
                            let copper_n: Spectrum = Spectrum::from_sampled(
                                &COPPER_WAVELENGTHS,
                                &COPPER_N,
                                COPPER_SAMPLES as i32,
                            );
                            let eta: Arc<dyn Texture<Spectrum> + Send + Sync> =
                                Arc::new(ConstantTexture::new(copper_n));
                            let copper_k: Spectrum = Spectrum::from_sampled(
                                &COPPER_WAVELENGTHS,
                                &COPPER_K,
                                COPPER_SAMPLES as i32,
                            );
                            let k: Arc<dyn Texture<Spectrum> + Send + Sync> =
                                Arc::new(ConstantTexture::new(copper_k));
                            let remap_roughness: bool = true;
                            let metal = Arc::new(Material::Metal(Box::new(MetalMaterial::new(
                                eta,
                                k,
                                Arc::new(ConstantTexture::new(mat.roughness as Float)),
                                None,
                                None,
                                None,
                                remap_roughness,
                            ))));
                            shapes.push(sphere.clone());
                            shape_materials.push(metal.clone());
                            shape_lights.push(None);
                        } else {
                            // MirrorMaterial
                            let kr = Arc::new(ConstantTexture::new(Spectrum::rgb(
                                mat.mirr * mat.ray_mirror,
                                mat.mirg * mat.ray_mirror,
                                mat.mirb * mat.ray_mirror,
                            )));
                            let mirror =
                                Arc::new(Material::Mirror(Box::new(MirrorMaterial::new(kr, None))));
                            shapes.push(sphere.clone());
                            shape_materials.push(mirror.clone());
                            shape_lights.push(None);
                        }
                    } else {
                        // MatteMaterial
                        let mut kd: Arc<dyn Texture<Spectrum> + Send + Sync> =
                            Arc::new(ConstantTexture::new(Spectrum::rgb(mat.r, mat.g, mat.b)));
                        if let Some(tex) = texture_hm.get(sphere_name) {
                            // first try texture with exactly the same name as the mesh
                            let su: Float = 1.0;
                            let sv: Float = 1.0;
                            let du: Float = 0.0;
                            let dv: Float = 0.0;
                            let mapping: Box<TextureMapping2D> =
                                Box::new(TextureMapping2D::UV(UVMapping2D {
                                    su: su,
                                    sv: sv,
                                    du: du,
                                    dv: dv,
                                }));
                            let filename: String = String::from(tex.to_str().unwrap());
                            let do_trilinear: bool = false;
                            let max_aniso: Float = 8.0;
                            let wrap_mode: ImageWrap = ImageWrap::Repeat;
                            let scale: Float = 1.0;
                            let gamma: bool = false;
                            kd = Arc::new(ImageTexture::new(
                                mapping,
                                filename,
                                do_trilinear,
                                max_aniso,
                                wrap_mode,
                                scale,
                                gamma,
                                convert_to_spectrum,
                            ));
                        } else {
                            // then remove trailing digits from mesh name
                            let mut ntd: String = String::new();
                            let mut chars = sphere_name.chars();
                            let mut digits: String = String::new(); // many digits
                            while let Some(c) = chars.next() {
                                if c.is_digit(10_u32) {
                                    // collect digits
                                    digits.push(c);
                                } else {
                                    // push collected digits (if any)
                                    ntd += &digits;
                                    // and reset
                                    digits = String::new();
                                    // push non-digit
                                    ntd.push(c);
                                }
                            }
                            // try no trailing digits (ntd)
                            if let Some(tex) = texture_hm.get(&ntd) {
                                let su: Float = 1.0;
                                let sv: Float = 1.0;
                                let du: Float = 0.0;
                                let dv: Float = 0.0;
                                let mapping: Box<TextureMapping2D> =
                                    Box::new(TextureMapping2D::UV(UVMapping2D {
                                        su: su,
                                        sv: sv,
                                        du: du,
                                        dv: dv,
                                    }));
                                let filename: String = String::from(tex.to_str().unwrap());
                                let do_trilinear: bool = false;
                                let max_aniso: Float = 8.0;
                                let wrap_mode: ImageWrap = ImageWrap::Repeat;
                                let scale: Float = 1.0;
                                let gamma: bool = false;
                                kd = Arc::new(ImageTexture::new(
                                    mapping,
                                    filename,
                                    do_trilinear,
                                    max_aniso,
                                    wrap_mode,
                                    scale,
                                    gamma,
                                    convert_to_spectrum,
                                ));
                            }
                        }
                        let sigma = Arc::new(ConstantTexture::new(0.0 as Float));
                        let matte = Arc::new(Material::Matte(Box::new(MatteMaterial::new(
                            kd,
                            sigma.clone(),
                            None,
                        ))));
                        shapes.push(sphere.clone());
                        shape_materials.push(matte.clone());
                        shape_lights.push(None);
                    }
                }
            }
        }
        // meshes
        for mesh_idx in 0..scene.meshes.len() {
            let mesh = &scene.meshes[mesh_idx];
            let mesh_name = &scene.mesh_names[mesh_idx];
            let triangle_colors = &scene.triangle_colors[mesh_idx];
            // create individual triangles
            let mut triangles: Vec<Arc<Shape>> = Vec::new();
            for id in 0..mesh.n_triangles {
                let triangle = Arc::new(Shape::Trngl(Triangle::new(mesh.clone(), id)));
                triangles.push(triangle.clone());
                shapes.push(triangle.clone());
            }
            if let Some(mat) = get_material(mesh_name, material_hm) {
                // println!("{:?}: {:?}", mesh_name, mat);
                if mat.emit > 0.0 {
                    has_emitters = true;
                    for i in 0..triangles.len() {
                        let triangle = &triangles[i];
                        let mi: MediumInterface = MediumInterface::default();
                        let l_emit: Spectrum = Spectrum::rgb(
                            mat.r * mat.emit * light_scale,
                            mat.g * mat.emit * light_scale,
                            mat.b * mat.emit * light_scale,
                        );
                        let n_samples: i32 = 1;
                        let two_sided: bool = false;
                        let area_light: Arc<Light> =
                            Arc::new(Light::DiffuseArea(Box::new(DiffuseAreaLight::new(
                                &mesh.object_to_world,
                                &mi,
                                &l_emit,
                                n_samples,
                                triangle.clone(),
                                two_sided,
                            ))));
                        lights.push(area_light.clone());
                        shape_materials.push(default_material.clone());
                        let triangle_light: Option<Arc<Light>> = Some(area_light.clone());
                        shape_lights.push(triangle_light);
                    }
                } else {
                    if mat.ang != 1.0 {
                        // GlassMaterial
                        let kr = Arc::new(ConstantTexture::new(Spectrum::new(1.0)));
                        let kt = Arc::new(ConstantTexture::new(Spectrum::rgb(
                            mat.specr, mat.specg, mat.specb,
                        )));
                        let u_roughness = Arc::new(ConstantTexture::new(0.0 as Float));
                        let v_roughness = Arc::new(ConstantTexture::new(0.0 as Float));
                        let index = Arc::new(ConstantTexture::new(mat.ang as Float));
                        let glass = Arc::new(Material::Glass(Box::new(GlassMaterial {
                            kr: kr,
                            kt: kt,
                            u_roughness: u_roughness,
                            v_roughness: v_roughness,
                            index: index,
                            bump_map: None,
                            remap_roughness: true,
                        })));
                        for _i in 0..triangles.len() {
                            shape_materials.push(glass.clone());
                            shape_lights.push(None);
                        }
                    } else if mat.ray_mirror > 0.0 {
                        if mat.roughness > 0.0 {
                            // MetalMaterial
                            let copper_n: Spectrum = Spectrum::from_sampled(
                                &COPPER_WAVELENGTHS,
                                &COPPER_N,
                                COPPER_SAMPLES as i32,
                            );
                            let eta: Arc<dyn Texture<Spectrum> + Send + Sync> =
                                Arc::new(ConstantTexture::new(copper_n));
                            let copper_k: Spectrum = Spectrum::from_sampled(
                                &COPPER_WAVELENGTHS,
                                &COPPER_K,
                                COPPER_SAMPLES as i32,
                            );
                            let k: Arc<dyn Texture<Spectrum> + Send + Sync> =
                                Arc::new(ConstantTexture::new(copper_k));
                            let remap_roughness: bool = true;
                            let metal = Arc::new(Material::Metal(Box::new(MetalMaterial::new(
                                eta,
                                k,
                                Arc::new(ConstantTexture::new(mat.roughness as Float)),
                                None,
                                None,
                                None,
                                remap_roughness,
                            ))));
                            for _i in 0..triangles.len() {
                                shape_materials.push(metal.clone());
                                shape_lights.push(None);
                            }
                        } else {
                            // MirrorMaterial
                            let kr = Arc::new(ConstantTexture::new(Spectrum::rgb(
                                mat.mirr * mat.ray_mirror,
                                mat.mirg * mat.ray_mirror,
                                mat.mirb * mat.ray_mirror,
                            )));
                            let mirror =
                                Arc::new(Material::Mirror(Box::new(MirrorMaterial::new(kr, None))));
                            for _i in 0..triangles.len() {
                                shape_materials.push(mirror.clone());
                                shape_lights.push(None);
                            }
                        }
                    } else {
                        // MatteMaterial
                        let mut kd: Arc<dyn Texture<Spectrum> + Send + Sync> =
                            Arc::new(ConstantTexture::new(Spectrum::rgb(mat.r, mat.g, mat.b)));
                        if let Some(tex) = texture_hm.get(mesh_name) {
                            // first try texture with exactly the same name as the mesh
                            let su: Float = 1.0;
                            let sv: Float = 1.0;
                            let du: Float = 0.0;
                            let dv: Float = 0.0;
                            let mapping: Box<TextureMapping2D> =
                                Box::new(TextureMapping2D::UV(UVMapping2D {
                                    su: su,
                                    sv: sv,
                                    du: du,
                                    dv: dv,
                                }));
                            let filename: String = String::from(tex.to_str().unwrap());
                            let do_trilinear: bool = false;
                            let max_aniso: Float = 8.0;
                            let wrap_mode: ImageWrap = ImageWrap::Repeat;
                            let scale: Float = 1.0;
                            let gamma: bool = false;
                            kd = Arc::new(ImageTexture::new(
                                mapping,
                                filename,
                                do_trilinear,
                                max_aniso,
                                wrap_mode,
                                scale,
                                gamma,
                                convert_to_spectrum,
                            ));
                        } else {
                            // then remove trailing digits from mesh name
                            let mut ntd: String = String::new();
                            let mut chars = mesh_name.chars();
                            let mut digits: String = String::new(); // many digits
                            while let Some(c) = chars.next() {
                                if c.is_digit(10_u32) {
                                    // collect digits
                                    digits.push(c);
                                } else {
                                    // push collected digits (if any)
                                    ntd += &digits;
                                    // and reset
                                    digits = String::new();
                                    // push non-digit
                                    ntd.push(c);
                                }
                            }
                            // try no trailing digits (ntd)
                            if let Some(tex) = texture_hm.get(&ntd) {
                                let su: Float = 1.0;
                                let sv: Float = 1.0;
                                let du: Float = 0.0;
                                let dv: Float = 0.0;
                                let mapping: Box<TextureMapping2D> =
                                    Box::new(TextureMapping2D::UV(UVMapping2D {
                                        su: su,
                                        sv: sv,
                                        du: du,
                                        dv: dv,
                                    }));
                                let filename: String = String::from(tex.to_str().unwrap());
                                let do_trilinear: bool = false;
                                let max_aniso: Float = 8.0;
                                let wrap_mode: ImageWrap = ImageWrap::Repeat;
                                let scale: Float = 1.0;
                                let gamma: bool = false;
                                kd = Arc::new(ImageTexture::new(
                                    mapping,
                                    filename,
                                    do_trilinear,
                                    max_aniso,
                                    wrap_mode,
                                    scale,
                                    gamma,
                                    convert_to_spectrum,
                                ));
                            }
                        }
                        let sigma = Arc::new(ConstantTexture::new(0.0 as Float));
                        let mut matte = Arc::new(Material::Matte(Box::new(MatteMaterial::new(
                            kd,
                            sigma.clone(),
                            None,
                        ))));
                        if triangle_colors.len() != 0_usize {
                            assert!(triangle_colors.len() == triangles.len());
                            // ignore textures, use triangle colors
                            for i in 0..triangles.len() {
                                // overwrite kd
                                kd = Arc::new(ConstantTexture::new(triangle_colors[i]));
                                matte = Arc::new(Material::Matte(Box::new(MatteMaterial::new(
                                    kd,
                                    sigma.clone(),
                                    None,
                                ))));
                                shape_materials.push(matte.clone());
                                shape_lights.push(None);
                            }
                        } else {
                            for _i in 0..triangles.len() {
                                shape_materials.push(matte.clone());
                                shape_lights.push(None);
                            }
                        }
                    }
                }
            } else {
                println!("{:?}: no mat", mesh_name);
                for _i in 0..triangles.len() {
                    shape_materials.push(default_material.clone());
                    shape_lights.push(None);
                }
            }
        }
        RenderOptions {
            has_emitters: has_emitters,
            primitives: primitives,
            shapes: shapes,
            shape_materials: shape_materials,
            shape_lights: shape_lights,
            lights: lights,
        }
    }
}

// TMP

fn get_material<'s, 'h>(
    mesh_name: &'s String,
    material_hm: &'h HashMap<String, Blend279Material>,
) -> Option<&'h Blend279Material> {
    // first try material with exactly the same name as the mesh
    if let Some(mat) = material_hm.get(mesh_name) {
        return Some(mat);
    } else {
        // then remove trailing digits from mesh name
        let mut ntd: String = String::new();
        let mut chars = mesh_name.chars();
        let mut digits: String = String::new(); // many digits
        while let Some(c) = chars.next() {
            if c.is_digit(10_u32) {
                // collect digits
                digits.push(c);
            } else {
                // push collected digits (if any)
                ntd += &digits;
                // and reset
                digits = String::new();
                // push non-digit
                ntd.push(c);
            }
        }
        // try no trailing digits (ntd)
        if let Some(mat) = material_hm.get(&ntd) {
            return Some(mat);
        } else {
            // finally try adding a '1' at the end
            ntd.push('1');
            if let Some(mat) = material_hm.get(&ntd) {
                return Some(mat);
            } else {
                println!("WARNING: No material found for {:?}", mesh_name);
                return None;
            }
        }
    }
    // material_hm.get(mesh_name)
}

fn read_mesh(
    base_name: &String,
    object_to_world_hm: &HashMap<String, Transform>,
    object_to_world: &mut Transform,
    p: &Vec<Point3f>,
    n: &Vec<Normal3f>,
    uvs: &mut Vec<Point2f>,
    loops: &Vec<u8>,
    vertex_indices: Vec<u32>,
    vertex_colors: Vec<u8>,
    is_smooth: bool,
    builder: &mut SceneDescriptionBuilder,
) {
    let mut instances: Vec<Transform> = Vec::with_capacity(u8::MAX as usize);
    let mut triangle_colors: Vec<Spectrum> = Vec::with_capacity(loops.len() * 2);
    if vertex_colors.len() != 0_usize {
        // let n_vertex_colors: usize = vertex_colors.len() / 4;
        // println!("{:?}: {} vertex colors found", base_name, n_vertex_colors);
        // println!("{:?}: {} loops found", base_name, loops.len());
        // println!("{:?}: {} vertices", base_name, vertex_indices.len());
        let mut rgba: usize = 0;
        for l in loops {
            for i in 0..*l {
                let r: Float = vertex_colors[rgba * 4 + 0] as Float / 255.0 as Float;
                let g: Float = vertex_colors[rgba * 4 + 1] as Float / 255.0 as Float;
                let b: Float = vertex_colors[rgba * 4 + 2] as Float / 255.0 as Float;
                // let a: Float = vertex_colors[rgba*4 + 3] as Float / 255.0 as Float;
                let s: Spectrum = Spectrum::rgb(r, g, b);
                // println!("{}: {:?}", rgba, s);
                if i == 0 {
                    // only store first color in loop (triangle or quad)
                    triangle_colors.push(s);
                    if *l == 4 {
                        // store it twice for quads (two triangles)
                        triangle_colors.push(s);
                    }
                }
                rgba += 1;
            }
        }
    }
    if let Some(o2w) = object_to_world_hm.get(base_name) {
        // one instance
        instances.push(*o2w);
    } else {
        // potentially many instances
        let mut counter: u8 = 0;
        let mut search_name: String;
        loop {
            counter += 1;
            search_name = base_name.to_string() + &counter.to_string();
            if let Some(o2w) = object_to_world_hm.get(&search_name) {
                instances.push(*o2w);
            } else {
                break;
            }
        }
    }
    for o2w in instances {
        *object_to_world = o2w;
        let world_to_object: Transform = Transform::inverse(&object_to_world);
        let n_triangles: usize = vertex_indices.len() / 3;
        // transform mesh vertices to world space
        let mut p_ws: Vec<Point3f> = Vec::with_capacity(p.len());
        let mut n_vertices: usize = p.len();
        for i in 0..n_vertices {
            p_ws.push(object_to_world.transform_point(&p[i]));
        }
        let mut n_ws: Vec<Normal3f> = Vec::with_capacity(n.len());
        if is_smooth {
            // println!("  is_smooth = {}", is_smooth);
            assert!(n.len() == p.len());
            if !n.is_empty() {
                for i in 0..n_vertices {
                    n_ws.push(object_to_world.transform_normal(&n[i]));
                }
            }
        }
        let s: Vec<Vector3f> = Vec::new();
        let mut uv: Vec<Point2f> = Vec::with_capacity(uvs.len());
        let mut new_vertex_indices: Vec<u32> = Vec::with_capacity(vertex_indices.len());
        if !uvs.is_empty() {
            let mut p_ws_vi: Vec<Point3f> = Vec::with_capacity(vertex_indices.len());
            let mut vertex_counter: u32 = 0;
            for vi in &vertex_indices {
                p_ws_vi.push(p_ws[*vi as usize]);
                new_vertex_indices.push(vertex_counter);
                vertex_counter += 1;
            }
            if is_smooth {
                let mut n_ws_vi: Vec<Normal3f> = Vec::with_capacity(vertex_indices.len());
                for vi in &vertex_indices {
                    n_ws_vi.push(n_ws[*vi as usize]);
                }
                n_ws = n_ws_vi;
            }
            let mut new_uvs: Vec<Point2f> = Vec::with_capacity(uvs.len());
            let mut loop_idx: usize = 0;
            for poly in loops {
                // triangle
                if *poly == 3_u8 {
                    for _i in 0..3 {
                        new_uvs.push(uvs[loop_idx]);
                        loop_idx += 1;
                    }
                }
                // quad
                else if *poly == 4_u8 {
                    new_uvs.push(uvs[loop_idx + 0]);
                    new_uvs.push(uvs[loop_idx + 1]);
                    new_uvs.push(uvs[loop_idx + 2]);
                    new_uvs.push(uvs[loop_idx + 0]);
                    new_uvs.push(uvs[loop_idx + 2]);
                    new_uvs.push(uvs[loop_idx + 3]);
                    loop_idx += 4;
                } else {
                    println!("WARNING: quads or triangles expected (poly = {})", poly)
                }
            }
            *uvs = new_uvs;
            assert!(
                uvs.len() == p_ws_vi.len(),
                "{} != {}",
                uvs.len(),
                p_ws_vi.len()
            );
            p_ws = p_ws_vi;
            n_vertices = p_ws.len();
            for i in 0..n_vertices {
                uv.push(uvs[i]);
            }
        }
        if new_vertex_indices.len() != 0 {
            builder.add_mesh(
                base_name.clone(),
                *object_to_world,
                world_to_object,
                n_triangles.try_into().unwrap(),
                new_vertex_indices.clone(),
                n_vertices.try_into().unwrap(),
                p_ws, // in world space
                s,    // empty
                n_ws, // in world space
                uv,
                triangle_colors.clone(),
            );
        } else {
            builder.add_mesh(
                base_name.clone(),
                *object_to_world,
                world_to_object,
                n_triangles.try_into().unwrap(),
                vertex_indices.clone(),
                n_vertices.try_into().unwrap(),
                p_ws, // in world space
                s,    // empty
                n_ws, // in world space
                uv,
                triangle_colors.clone(),
            );
        }
    }
}

pub fn make_perspective_camera(
    filter_width: Float,
    xres: i32,
    yres: i32,
    fov: Float,
    animated_cam_to_world: AnimatedTransform,
    clipsta: Float,
) -> Option<Arc<Camera>> {
    let mut some_camera: Option<Arc<Camera>> = None;
    let mut filter_params: ParamSet = ParamSet::default();
    filter_params.add_float(String::from("xwidth"), filter_width);
    filter_params.add_float(String::from("ywidth"), filter_width);
    let some_filter = make_filter(&String::from("gaussian"), &filter_params);
    if let Some(filter) = some_filter {
        let film_name: String = String::from("image");
        let mut film_params: ParamSet = ParamSet::default();
        film_params.add_int(String::from("xresolution"), xres);
        film_params.add_int(String::from("yresolution"), yres);
        let crop_window: Bounds2f = Bounds2f {
            p_min: Point2f { x: 0.0, y: 0.0 },
            p_max: Point2f { x: 1.0, y: 1.0 },
        };
        let some_film: Option<Arc<Film>> =
            make_film(&film_name, &film_params, filter, &crop_window);
        if let Some(film) = some_film {
            let camera_name: String = String::from("perspective");
            let mut camera_params: ParamSet = ParamSet::default();
            camera_params.add_float(String::from("fov"), fov);
            some_camera = make_camera(
                &camera_name,
                &camera_params,
                animated_cam_to_world,
                film,
                clipsta,
            );
        }
    }
    some_camera
}

fn make_integrator(
    integrator_name: &String,
    filter_width: Float,
    xres: i32,
    yres: i32,
    fov: Float,
    clipsta: Float,
    animated_cam_to_world: AnimatedTransform,
    pixelsamples: i32,
    integrator_params: ParamSet,
) -> Option<Box<Integrator>> {
    let mut some_integrator: Option<Box<Integrator>> = None;
    let some_camera: Option<Arc<Camera>> = make_perspective_camera(
        filter_width,
        xres,
        yres,
        fov,
        animated_cam_to_world,
        clipsta,
    );
    if let Some(camera) = some_camera {
        let sampler_name: String = String::from("halton");
        let mut sampler_params: ParamSet = ParamSet::default();
        sampler_params.add_int(String::from("pixelsamples"), pixelsamples);
        let some_sampler: Option<Box<Sampler>> =
            make_sampler(&sampler_name, &sampler_params, camera.get_film());
        if let Some(sampler) = some_sampler {
            print!("integrator = {:?} [", integrator_name);
            if integrator_name == "whitted" {
                println!("Whitted’s Ray-Tracing]");
                println!("  pixelsamples = {}", pixelsamples);
                // CreateWhittedIntegrator
                let max_depth: i32 = integrator_params.find_one_int("maxdepth", 5);
                let pixel_bounds: Bounds2i = camera.get_film().get_sample_bounds();
                let integrator = Box::new(Integrator::Sampler(SamplerIntegrator::Whitted(
                    WhittedIntegrator::new(max_depth as u32, camera, sampler, pixel_bounds),
                )));
                some_integrator = Some(integrator);
            } else if integrator_name == "directlighting" {
                println!("Direct Lighting]");
                println!("  pixelsamples = {}", pixelsamples);
                // CreateDirectLightingIntegrator
                let max_depth: i32 = integrator_params.find_one_int("maxdepth", 5);
                println!("  max_depth = {}", max_depth);
                let st: String = integrator_params.find_one_string("strategy", String::from("all"));
                let strategy: LightStrategy;
                if st == "one" {
                    strategy = LightStrategy::UniformSampleOne;
                } else if st == "all" {
                    strategy = LightStrategy::UniformSampleAll;
                } else {
                    panic!("Strategy \"{}\" for direct lighting unknown.", st);
                }
                let pixel_bounds: Bounds2i = Bounds2i {
                    p_min: Point2i { x: 0, y: 0 },
                    p_max: Point2i { x: xres, y: yres },
                };
                let integrator = Box::new(Integrator::Sampler(SamplerIntegrator::DirectLighting(
                    DirectLightingIntegrator::new(
                        strategy,
                        max_depth as u32,
                        camera,
                        sampler,
                        pixel_bounds,
                    ),
                )));
                some_integrator = Some(integrator);
            } else if integrator_name == "path" {
                println!("(Unidirectional) Path Tracing]");
                println!("  pixelsamples = {}", pixelsamples);
                // CreatePathIntegrator
                let max_depth: i32 = integrator_params.find_one_int("maxdepth", 5);
                println!("  max_depth = {}", max_depth);
                let pixel_bounds: Bounds2i = camera.get_film().get_sample_bounds();
                let rr_threshold: Float = 1.0;
                let light_strategy: String = String::from("spatial");
                let integrator = Box::new(Integrator::Sampler(SamplerIntegrator::Path(
                    PathIntegrator::new(
                        max_depth as u32,
                        camera,
                        sampler,
                        pixel_bounds,
                        rr_threshold,
                        light_strategy,
                    ),
                )));
                some_integrator = Some(integrator);
            } else if integrator_name == "volpath" {
                println!("Path Tracing (Participating Media)]");
                println!("  pixelsamples = {}", pixelsamples);
                // CreateVolPathIntegrator
                let max_depth: i32 = integrator_params.find_one_int("maxdepth", 5);
                let pixel_bounds: Bounds2i = camera.get_film().get_sample_bounds();
                let rr_threshold: Float = 1.0;
                let light_strategy: String = String::from("spatial");
                let integrator = Box::new(Integrator::Sampler(SamplerIntegrator::VolPath(
                    VolPathIntegrator::new(
                        max_depth as u32,
                        camera,
                        sampler,
                        pixel_bounds,
                        rr_threshold,
                        light_strategy,
                    ),
                )));
                some_integrator = Some(integrator);
            } else if integrator_name == "bdpt" {
                println!("Bidirectional Path Tracing (BDPT)]");
                println!("  pixelsamples = {}", pixelsamples);
                // CreateBDPTIntegrator
                let max_depth: i32 = integrator_params.find_one_int("maxdepth", 5);
                println!("  max_depth = {}", max_depth);
                let pixel_bounds: Bounds2i = camera.get_film().get_sample_bounds();
                let light_strategy: String = String::from("power");
                let integrator = Box::new(Integrator::BDPT(BDPTIntegrator::new(
                    camera,
                    sampler,
                    pixel_bounds,
                    max_depth as u32,
                    light_strategy,
                )));
                some_integrator = Some(integrator);
            } else if integrator_name == "mlt" {
                println!("Metropolis Light Transport (MLT)]");
                println!("  pixelsamples = {}", pixelsamples);
                // CreateMLTIntegrator
                let max_depth: i32 = integrator_params.find_one_int("maxdepth", 5);
                println!("  max_depth = {}", max_depth);
                let n_bootstrap: i32 = integrator_params.find_one_int("bootstrapsamples", 100000);
                println!("  bootstrap_samples = {}", n_bootstrap);
                let n_chains: i32 = integrator_params.find_one_int("chains", 1000);
                println!("  chains = {}", n_chains);
                let mutations_per_pixel: i32 =
                    integrator_params.find_one_int("mutationsperpixel", 100);
                println!("  mutations_per_pixel = {}", mutations_per_pixel);
                let large_step_probability: Float =
                    integrator_params.find_one_float("largestepprobability", 0.3 as Float);
                println!("  step_probability = {}", large_step_probability);
                let sigma: Float = integrator_params.find_one_float("sigma", 0.01 as Float);
                println!("  sigma = {}", sigma);
                let integrator = Box::new(Integrator::MLT(MLTIntegrator::new(
                    camera.clone(),
                    max_depth as u32,
                    n_bootstrap as u32,
                    n_chains as u32,
                    mutations_per_pixel as u32,
                    sigma,
                    large_step_probability,
                )));
                some_integrator = Some(integrator);
            } else if integrator_name == "ao" || integrator_name == "ambientocclusion" {
                println!("Ambient Occlusion (AO)]");
                println!("  pixelsamples = {}", pixelsamples);
                // CreateAOIntegrator
                let pixel_bounds: Bounds2i = camera.get_film().get_sample_bounds();
                let cos_sample: bool = integrator_params.find_one_bool("cossample", true);
                let n_samples: i32 = integrator_params.find_one_int("nsamples", 64 as i32);
                let integrator = Box::new(Integrator::Sampler(SamplerIntegrator::AO(
                    AOIntegrator::new(cos_sample, n_samples, camera, sampler, pixel_bounds),
                )));
                some_integrator = Some(integrator);
            } else if integrator_name == "sppm" {
                println!("Stochastic Progressive Photon Mapping (SPPM)]");
                println!("  pixelsamples = {}", pixelsamples);
                // CreateSPPMIntegrator
                let mut n_iterations: i32 = integrator_params.find_one_int("numiterations", 64);
                n_iterations = integrator_params.find_one_int("iterations", n_iterations);
                let max_depth: i32 = integrator_params.find_one_int("maxdepth", 5);
                let photons_per_iter: i32 =
                    integrator_params.find_one_int("photonsperiteration", -1);
                let write_freq: i32 =
                    integrator_params.find_one_int("imagewritefrequency", 1 << 31);
                println!("  imagewritefrequency = {}", write_freq);
                let radius: Float = integrator_params.find_one_float("radius", 1.0 as Float);
                // TODO: if (PbrtOptions.quickRender) nIterations = std::max(1, nIterations / 16);
                let integrator = Box::new(Integrator::SPPM(SPPMIntegrator::new(
                    camera.clone(),
                    n_iterations,
                    photons_per_iter,
                    max_depth as u32,
                    radius,
                    write_freq,
                )));
                some_integrator = Some(integrator);
            } else {
                println!("Integrator \"{}\" unknown.", integrator_name);
            }
        } else {
            panic!("Unable to create sampler.");
        }
    } else {
        panic!("Unable to create camera.");
    }
    some_integrator
}

fn make_scene(primitives: &Vec<Arc<Primitive>>, lights: Vec<Arc<Light>>) -> Scene {
    let accelerator_name: String = String::from("bvh");
    let some_accelerator = make_accelerator(&accelerator_name, &primitives, &ParamSet::default());
    if let Some(accelerator) = some_accelerator {
        return Scene::new(accelerator, lights);
    } else {
        panic!("Unable to create accelerator.");
    }
}

fn main() -> std::io::Result<()> {
    let args = Cli::from_args();
    let git_describe = option_env!("GIT_DESCRIBE").unwrap_or("unknown");
    let num_threads: u8 = num_cpus::get() as u8;
    println!(
        "parse_blend_file version {} ({}) [Detected {} cores]",
        VERSION, git_describe, num_threads
    );
    // PBRT
    let mut scale_length: f32 = 1.0;
    let mut resolution_x: u32 = 640;
    let mut resolution_y: u32 = 480;
    let mut resolution_percentage: i16 = 100;
    let mut angle_x: f32 = 45.0;
    let mut angle_y: f32 = 45.0;
    let mut base_name = String::new();
    let mut camera_hm: HashMap<String, BlendCamera> = HashMap::new();
    let mut texture_hm: HashMap<String, OsString> = HashMap::new();
    // let mut spheres_hm: HashMap<String, PbrtSphere> = HashMap::new();
    // let mut cylinders_hm: HashMap<String, PbrtCylinder> = HashMap::new();
    // let mut disks_hm: HashMap<String, PbrtDisk> = HashMap::new();
    let mut ob_count: usize = 0;
    let mut object_to_world: Transform = Transform::new(
        1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    );
    let mut p: Vec<Point3f> = Vec::with_capacity(1048576);
    let mut n: Vec<Normal3f> = Vec::with_capacity(1048576);
    let mut uvs: Vec<Point2f> = Vec::with_capacity(1048576);
    let mut loops: Vec<u8> = Vec::with_capacity(1048576);
    let mut vertex_indices: Vec<u32> = Vec::with_capacity(1048576);
    let mut vertex_colors: Vec<u8> = Vec::with_capacity(1048576);
    let mut loop_indices: Vec<i32> = Vec::with_capacity(4194304);
    // let mut prop_height: f64 = 0.0;
    // let mut prop_radius: f64 = 1.0;
    // let mut prop_innerradius: f64 = 0.0;
    // let mut prop_zmin: f64 = -1.0;
    // let mut prop_zmax: f64 = 1.0;
    // let mut prop_phimax: f64 = 360.0;
    let mut hdr_path: OsString = OsString::new();
    // read DNA
    let mut dna_types_hm: HashMap<String, u16> = HashMap::new();
    let mut dna_structs_hm: HashMap<String, DnaStrC> = HashMap::new();
    let mut dna_pointers_hm: HashMap<usize, usize> = HashMap::new();
    let mut dna_2_type_id: Vec<u16> = Vec::new();
    let mut types: Vec<String> = Vec::new();
    let mut num_bytes_read: usize = 0;
    let print_dna: bool = false;
    let print_pointers: bool = false;
    let verbose: bool = false;
    read_dna(
        print_dna,
        print_pointers,
        &args.path,
        &mut dna_types_hm,
        &mut dna_structs_hm,
        &mut dna_pointers_hm,
        &mut dna_2_type_id,
        &mut types,
        &mut num_bytes_read,
    )?;
    if verbose {
        println!("{} {:?}", num_bytes_read, &args.path);
    }
    let names: Vec<String> = vec![
        "Scene".to_string(),
        "Object".to_string(),
        "Camera".to_string(),
        "Lamp".to_string(),
        "Material".to_string(),
        "Image".to_string(),
        "Mesh".to_string(),
        "MPoly".to_string(),
        "MVert".to_string(),
        "MLoop".to_string(),
        "MLoopUV".to_string(),
        "MLoopCol".to_string(),
        "bNodeTree".to_string(),
        "bNode".to_string(),
        "bNodeSocket".to_string(),
        "bNodeSocketValueFloat".to_string(),
    ];
    // then use the DNA
    let mut bytes_read: Vec<u8> = Vec::with_capacity(num_bytes_read);
    let mut structs_read: Vec<String> = Vec::with_capacity(names.len());
    let mut data_read: Vec<u32> = Vec::with_capacity(names.len());
    let mut pointers_read: Vec<(usize, u32)> = Vec::with_capacity(names.len());
    use_dna(
        print_dna,
        &args.path,
        &dna_types_hm,
        &dna_structs_hm,
        &names,
        &dna_2_type_id,
        &types,
        &mut bytes_read,
        &mut structs_read,
        &mut data_read,
        &mut pointers_read,
    )?;
    if verbose {
        println!(
            "Do something with the {} bytes returned by use_dna(...).",
            bytes_read.len()
        );
    }
    // store data for pbrt here
    let mut material_hm: HashMap<String, Blend279Material> = HashMap::with_capacity(ob_count);
    let mut object_to_world_hm: HashMap<String, Transform> = HashMap::with_capacity(ob_count);
    let mut builder: SceneDescriptionBuilder = SceneDescriptionBuilder::new();
    let mut data_following_mesh: bool = false;
    let mut data_following_material: bool = false;
    let mut is_smooth: bool = false;
    let parent = args.path.parent().unwrap();
    // emit (use nodes or old material settings?)
    let mut emit: f32 = 0.0;
    let mut search_for_emit: bool = false;
    let mut emit_default_value: usize = 0;
    // structs_read
    let mut byte_index: usize = 0;
    let mut struct_index: usize = 0;
    for struct_read in structs_read {
        if verbose {
            println!(
                "{} ({} - {})",
                struct_read, byte_index, data_read[struct_index]
            );
        }
        if let Some(tlen) = dna_types_hm.get(&struct_read) {
            if data_read[struct_index] == *tlen as u32 {
                // single structs
                if let Some(struct_found) = dna_structs_hm.get(&struct_read) {
                    match struct_read.as_str() {
                        "Scene" => {
                            for member in &struct_found.members {
                                match member.mem_name.as_str() {
                                    "id" => {
                                        let id: String = get_id_name(
                                            member,
                                            &bytes_read,
                                            byte_index,
                                            &dna_structs_hm,
                                            &dna_types_hm,
                                        );
                                        if verbose {
                                            println!("  ID.name = {:?}", id);
                                        }
                                        base_name = id.clone()[2..].to_string();
                                    }
                                    "r" => {
                                        if let Some(struct_found2) =
                                            dna_structs_hm.get(member.mem_type.as_str())
                                        {
                                            let mut byte_index2: usize = 0;
                                            for member2 in &struct_found2.members {
                                                if let Some(type_found2) =
                                                    dna_types_hm.get(&member2.mem_type)
                                                {
                                                    let mem_tlen2: u16 =
                                                        calc_mem_tlen(member2, *type_found2);
                                                    if member2.mem_name.as_str() == "size" {
                                                        resolution_percentage = get_short(
                                                            member2,
                                                            &bytes_read,
                                                            byte_index + byte_index2,
                                                        );
                                                        byte_index2 += mem_tlen2 as usize;
                                                    } else if member2.mem_name.as_str() == "xsch" {
                                                        let xsch = get_int(
                                                            member2,
                                                            &bytes_read,
                                                            byte_index + byte_index2,
                                                        );
                                                        resolution_x = xsch as u32;
                                                        byte_index2 += mem_tlen2 as usize;
                                                    } else if member2.mem_name.as_str() == "ysch" {
                                                        let ysch = get_int(
                                                            member2,
                                                            &bytes_read,
                                                            byte_index + byte_index2,
                                                        );
                                                        resolution_y = ysch as u32;
                                                        byte_index2 += mem_tlen2 as usize;
                                                    } else {
                                                        byte_index2 += mem_tlen2 as usize;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    "unit" => {
                                        if let Some(struct_found2) =
                                            dna_structs_hm.get(member.mem_type.as_str())
                                        {
                                            let mut byte_index2: usize = 0;
                                            for member2 in &struct_found2.members {
                                                if let Some(type_found2) =
                                                    dna_types_hm.get(&member2.mem_type)
                                                {
                                                    let mem_tlen2: u16 =
                                                        calc_mem_tlen(member2, *type_found2);
                                                    if member2.mem_name.contains("scale_length") {
                                                        scale_length = get_float(
                                                            member2,
                                                            &bytes_read,
                                                            byte_index + byte_index2,
                                                        );
                                                        byte_index2 += mem_tlen2 as usize;
                                                    } else {
                                                        byte_index2 += mem_tlen2 as usize;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    _ => {}
                                }
                                // find mem_type in dna_types.names
                                if let Some(type_found) = dna_types_hm.get(&member.mem_type) {
                                    let mem_tlen: u16 = calc_mem_tlen(member, *type_found);
                                    byte_index += mem_tlen as usize;
                                }
                            }
                            if verbose {
                                println!("  scale_length = {}", scale_length);
                            }
                            // reset booleans
                            data_following_mesh = false;
                            data_following_material = false;
                            is_smooth = false;
                        }
                        "Object" => {
                            ob_count += 1;
                            for member in &struct_found.members {
                                match member.mem_name.as_str() {
                                    "id" => {
                                        let id: String = get_id_name(
                                            member,
                                            &bytes_read,
                                            byte_index,
                                            &dna_structs_hm,
                                            &dna_types_hm,
                                        );
                                        if verbose {
                                            println!("  ID.name = {:?}", id);
                                        }
                                        base_name = id.clone()[2..].to_string();
                                    }
                                    "obmat[4][4]" => {
                                        let obmat: [f32; 16] =
                                            get_matrix(member, &bytes_read, byte_index);
                                        if verbose {
                                            println!("  obmat[4][4] = {:?}", obmat);
                                            println!("  scale_length = {}", scale_length);
                                        }
                                        object_to_world = Transform::new(
                                            obmat[0],
                                            obmat[4],
                                            obmat[8],
                                            obmat[12] * scale_length,
                                            obmat[1],
                                            obmat[5],
                                            obmat[9],
                                            obmat[13] * scale_length,
                                            obmat[2],
                                            obmat[6],
                                            obmat[10],
                                            obmat[14] * scale_length,
                                            obmat[3],
                                            obmat[7],
                                            obmat[11],
                                            obmat[15],
                                        );
                                        if verbose {
                                            println!("  object_to_world = {:?}", object_to_world);
                                        }
                                        object_to_world_hm
                                            .insert(base_name.clone(), object_to_world);
                                    }
                                    _ => {}
                                }
                                // find mem_type in dna_types.names
                                if let Some(type_found) = dna_types_hm.get(&member.mem_type) {
                                    let mem_tlen: u16 = calc_mem_tlen(member, *type_found);
                                    byte_index += mem_tlen as usize;
                                }
                            }
                            // reset booleans
                            data_following_mesh = false;
                            data_following_material = false;
                            is_smooth = false;
                        }
                        "Camera" => {
                            let mut lens: f32 = 0.0;
                            let mut clipsta: f32 = 0.0;
                            let mut sensor_x: f32 = 0.0;
                            let mut sensor_y: f32 = 0.0;
                            for member in &struct_found.members {
                                match member.mem_name.as_str() {
                                    "id" => {
                                        let id: String = get_id_name(
                                            member,
                                            &bytes_read,
                                            byte_index,
                                            &dna_structs_hm,
                                            &dna_types_hm,
                                        );
                                        if verbose {
                                            println!("  ID.name = {:?}", id);
                                        }
                                        base_name = id.clone()[2..].to_string();
                                    }
                                    "lens" => {
                                        lens = get_float(member, &bytes_read, byte_index);
                                    }
                                    "clipsta" => {
                                        clipsta = get_float(member, &bytes_read, byte_index);
                                    }
                                    "sensor_x" => {
                                        sensor_x = get_float(member, &bytes_read, byte_index);
                                    }
                                    "sensor_y" => {
                                        sensor_y = get_float(member, &bytes_read, byte_index);
                                    }
                                    _ => {}
                                }
                                // find mem_type in dna_types.names
                                if let Some(type_found) = dna_types_hm.get(&member.mem_type) {
                                    let mem_tlen: u16 = calc_mem_tlen(member, *type_found);
                                    byte_index += mem_tlen as usize;
                                }
                            }
                            // calculate angle_x and angle_y
                            angle_x = degrees(focallength_to_fov(lens, sensor_x) as Float);
                            angle_y = degrees(focallength_to_fov(lens, sensor_y) as Float);
                            let cam: BlendCamera = BlendCamera {
                                lens,
                                // angle_x,
                                // angle_y,
                                clipsta,
                            };
                            if verbose {
                                println!("  {:?}", cam);
                            }
                            camera_hm.insert(base_name.clone(), cam);
                            // reset booleans
                            data_following_mesh = false;
                            data_following_material = false;
                            is_smooth = false;
                        }
                        "Lamp" => {
                            let mut la_type: i16 = 0;
                            let mut r: f32 = 0.0;
                            let mut g: f32 = 0.0;
                            let mut b: f32 = 0.0;
                            let mut energy: f32 = 0.0;
                            for member in &struct_found.members {
                                match member.mem_name.as_str() {
                                    "id" => {
                                        let id: String = get_id_name(
                                            member,
                                            &bytes_read,
                                            byte_index,
                                            &dna_structs_hm,
                                            &dna_types_hm,
                                        );
                                        if verbose {
                                            println!("  ID.name = {:?}", id);
                                        }
                                        base_name = id.clone()[2..].to_string();
                                    }
                                    "type" => {
                                        la_type = get_short(member, &bytes_read, byte_index);
                                        if verbose {
                                            println!("  type = {}", la_type);
                                        }
                                    }
                                    "r" => {
                                        r = get_float(member, &bytes_read, byte_index);
                                    }
                                    "g" => {
                                        g = get_float(member, &bytes_read, byte_index);
                                    }
                                    "b" => {
                                        b = get_float(member, &bytes_read, byte_index);
                                    }
                                    "energy" => {
                                        energy = get_float(member, &bytes_read, byte_index);
                                        if verbose {
                                            println!("  energy = {}", energy);
                                        }
                                    }
                                    _ => {}
                                }
                                // find mem_type in dna_types.names
                                if let Some(type_found) = dna_types_hm.get(&member.mem_type) {
                                    let mem_tlen: u16 = calc_mem_tlen(member, *type_found);
                                    byte_index += mem_tlen as usize;
                                }
                            }
                            if la_type == 0 {
                                // LA_LOCAL
                                if let Some(o2w) = object_to_world_hm.get(&base_name) {
                                    object_to_world = *o2w;
                                } else {
                                    println!(
                                        "WARNING: looking up object_to_world by name ({:?}) failed",
                                        base_name
                                    );
                                }
                                let l: Spectrum = Spectrum::rgb(r, g, b);
                                if verbose {
                                    println!("  l = {:?}", l);
                                }
                                // point light
                                builder.add_point_light(
                                    object_to_world,
                                    l,
                                    args.light_scale * energy,
                                );
                            } else if la_type == 1 {
                                // LA_SUN
                                if let Some(o2w) = object_to_world_hm.get(&base_name) {
                                    object_to_world = *o2w;
                                } else {
                                    println!(
                                        "WARNING: looking up object_to_world by name ({:?}) failed",
                                        base_name
                                    );
                                }
                                let l: Spectrum = Spectrum::rgb(r, g, b);
                                if verbose {
                                    println!("  l = {:?}", l);
                                }
                                // distant light
                                builder.add_distant_light(
                                    object_to_world,
                                    l,
                                    args.light_scale * energy,
                                );
                            } else {
                                println!("WARNING: la_type = {} not supported (yet)", la_type);
                            }
                            // reset booleans
                            data_following_mesh = false;
                            data_following_material = false;
                            is_smooth = false;
                        }
                        "Material" => {
			    emit = 0.0; // reset
                            if data_following_mesh {
                                // time to use the gathered data to create a mesh
                                read_mesh(
                                    &base_name,
                                    &object_to_world_hm,
                                    &mut object_to_world,
                                    &p,
                                    &n,
                                    &mut uvs,
                                    &loops,
                                    vertex_indices.clone(),
                                    vertex_colors.clone(),
                                    is_smooth,
                                    &mut builder,
                                );
                                // clear all Vecs
                                p.clear();
                                n.clear();
                                uvs.clear();
                                loops.clear();
                                vertex_indices.clear();
                                vertex_colors.clear();
                                loop_indices.clear();
                            }
                            let mut r: f32 = 0.0;
                            let mut g: f32 = 0.0;
                            let mut b: f32 = 0.0;
                            let mut specr: f32 = 0.0;
                            let mut specg: f32 = 0.0;
                            let mut specb: f32 = 0.0;
                            let mut mirr: f32 = 0.0;
                            let mut mirg: f32 = 0.0;
                            let mut mirb: f32 = 0.0;
                            let mut ang: f32 = 0.0;
                            let mut ray_mirror: f32 = 0.0;
                            let mut roughness: f32 = 0.0;
                            let mut use_nodes: u8 = 0;
                            let mut nodetree: usize = 0;
                            for member in &struct_found.members {
                                match member.mem_name.as_str() {
                                    "id" => {
                                        let id: String = get_id_name(
                                            member,
                                            &bytes_read,
                                            byte_index,
                                            &dna_structs_hm,
                                            &dna_types_hm,
                                        );
                                        if verbose {
                                            println!("  ID.name = {:?}", id);
                                        }
                                        base_name = id.clone()[2..].to_string();
                                    }
                                    "r" => {
                                        r = get_float(member, &bytes_read, byte_index);
                                    }
                                    "g" => {
                                        g = get_float(member, &bytes_read, byte_index);
                                    }
                                    "b" => {
                                        b = get_float(member, &bytes_read, byte_index);
                                    }
                                    "specr" => {
                                        specr = get_float(member, &bytes_read, byte_index);
                                    }
                                    "specg" => {
                                        specg = get_float(member, &bytes_read, byte_index);
                                    }
                                    "specb" => {
                                        specb = get_float(member, &bytes_read, byte_index);
                                    }
                                    "mirr" => {
                                        mirr = get_float(member, &bytes_read, byte_index);
                                    }
                                    "mirg" => {
                                        mirg = get_float(member, &bytes_read, byte_index);
                                    }
                                    "mirb" => {
                                        mirb = get_float(member, &bytes_read, byte_index);
                                    }
                                    "emit" => {
                                        emit = get_float(member, &bytes_read, byte_index);
                                    }
                                    "ang" => {
                                        ang = get_float(member, &bytes_read, byte_index);
                                    }
                                    "ray_mirror" => {
                                        ray_mirror = get_float(member, &bytes_read, byte_index);
                                    }
                                    "roughness" => {
                                        roughness = get_float(member, &bytes_read, byte_index);
                                    }
                                    "use_nodes" => {
                                        use_nodes = get_char(member, &bytes_read, byte_index);
                                    }
                                    "*nodetree" => {
                                        nodetree = get_pointer(member, &bytes_read, byte_index);
                                    }
                                    _ => {}
                                }
                                // find mem_type in dna_types.names
                                if let Some(type_found) = dna_types_hm.get(&member.mem_type) {
                                    let mem_tlen: u16 = calc_mem_tlen(member, *type_found);
                                    byte_index += mem_tlen as usize;
                                }
                            }
                            if emit == 0.0 && use_nodes == 1 {
                                search_for_emit = true;
                                println!(
                                    "{} (SDNAnr = {}) ({:#018x})",
                                    struct_read,
                                    pointers_read[struct_index].1,
                                    pointers_read[struct_index].0
                                );
                                println!("use_nodes = {}", use_nodes);
                                if let Some(pointer_found) = dna_pointers_hm.get(&nodetree) {
                                    println!(
                                        "nodetree = {:#010x} ({:#018x})",
                                        pointer_found, nodetree
                                    );
                                }
                            } else {
                                search_for_emit = false;
                            }
                            // Blend279Material
                            let mat: Blend279Material = Blend279Material {
                                r: r,
                                g: g,
                                b: b,
                                // a: 1.0,
                                specr: specr,
                                specg: specg,
                                specb: specb,
                                mirr: mirr,
                                mirg: mirg,
                                mirb: mirb,
                                emit: emit,
                                ang: ang,
                                ray_mirror: ray_mirror,
                                roughness: roughness,
                            };
                            if verbose {
                                println!("  mat[{:?}] = {:?}", base_name, mat);
                            }
                            material_hm.insert(base_name.clone(), mat);
                            // reset booleans
                            data_following_mesh = false;
                            // data_following_material = false;
                            is_smooth = false;
                            // data
                            data_following_material = true;
                        }
                        "Image" => {
                            for member in &struct_found.members {
                                match member.mem_name.as_str() {
                                    "id" => {
                                        let id: String = get_id_name(
                                            member,
                                            &bytes_read,
                                            byte_index,
                                            &dna_structs_hm,
                                            &dna_types_hm,
                                        );
                                        if verbose {
                                            println!("  ID.name = {:?}", id);
                                        }
                                        base_name = id.clone()[2..].to_string();
                                    }
                                    "name[1024]" => {
                                        if let Some(type_found) = dna_types_hm.get(&member.mem_type)
                                        {
                                            let mem_tlen: u16 = calc_mem_tlen(member, *type_found);
                                            let mut name = String::with_capacity(mem_tlen as usize);
                                            for i in 0..mem_tlen as usize {
                                                if bytes_read[byte_index + i] == 0 {
                                                    break;
                                                }
                                                name.push(bytes_read[byte_index + i] as char);
                                            }
                                            if name.len() > 2 {
                                                // println!("  name = {}", name);
                                                let image_path: &Path = Path::new(&name);
                                                // println!("  image_path = {:?}", image_path);
                                                if let Some(img_ext) = image_path.extension() {
                                                    // println!("  img_ext = {:?}", img_ext);
                                                    if img_ext == "hdr" {
                                                        if image_path.starts_with("//") {
                                                            if let Ok(relative) =
                                                                image_path.strip_prefix("//")
                                                            {
                                                                let canonicalized = parent
                                                                    .join(relative.clone())
                                                                    .canonicalize()
                                                                    .unwrap();
                                                                println!("{:?}", canonicalized);
                                                                hdr_path =
                                                                    canonicalized.into_os_string();
                                                            }
                                                        }
                                                    } else {
                                                        if image_path.starts_with("//") {
                                                            if let Ok(relative) =
                                                                image_path.strip_prefix("//")
                                                            {
                                                                let canonicalized = parent
                                                                    .join(relative.clone())
                                                                    .canonicalize()
                                                                    .unwrap();
                                                                println!("{:?}", canonicalized);
                                                                texture_hm.insert(
                                                                    base_name.clone(),
                                                                    canonicalized.into_os_string(),
                                                                );
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    _ => {}
                                }
                                // find mem_type in dna_types.names
                                if let Some(type_found) = dna_types_hm.get(&member.mem_type) {
                                    let mem_tlen: u16 = calc_mem_tlen(member, *type_found);
                                    byte_index += mem_tlen as usize;
                                }
                            }
                            // reset booleans
                            data_following_mesh = false;
                            data_following_material = false;
                            is_smooth = false;
                        }
                        "Mesh" => {
                            if data_following_mesh {
                                // time to use the gathered data to create a mesh
                                read_mesh(
                                    &base_name,
                                    &object_to_world_hm,
                                    &mut object_to_world,
                                    &p,
                                    &n,
                                    &mut uvs,
                                    &loops,
                                    vertex_indices.clone(),
                                    vertex_colors.clone(),
                                    is_smooth,
                                    &mut builder,
                                );
                                // clear all Vecs
                                p.clear();
                                n.clear();
                                uvs.clear();
                                loops.clear();
                                vertex_indices.clear();
                                vertex_colors.clear();
                                loop_indices.clear();
                            }
                            let mut totvert: i32;
                            let mut totedge: i32;
                            let mut totface: i32;
                            let mut totselect: i32;
                            let mut totpoly: i32;
                            let mut totloop: i32;
                            for member in &struct_found.members {
                                match member.mem_name.as_str() {
                                    "id" => {
                                        let id: String = get_id_name(
                                            member,
                                            &bytes_read,
                                            byte_index,
                                            &dna_structs_hm,
                                            &dna_types_hm,
                                        );
                                        if verbose {
                                            println!("  ID.name = {:?}", id);
                                        }
                                        base_name = id.clone()[2..].to_string();
                                    }
                                    "totvert" => {
                                        totvert = get_int(member, &bytes_read, byte_index);
                                        if verbose {
                                            println!("  totvert = {:?}", totvert);
                                        }
                                    }
                                    "totedge" => {
                                        totedge = get_int(member, &bytes_read, byte_index);
                                        if verbose {
                                            println!("  totedge = {:?}", totedge);
                                        }
                                    }
                                    "totface" => {
                                        totface = get_int(member, &bytes_read, byte_index);
                                        if verbose {
                                            println!("  totface = {:?}", totface);
                                        }
                                    }
                                    "totselect" => {
                                        totselect = get_int(member, &bytes_read, byte_index);
                                        if verbose {
                                            println!("  totselect = {:?}", totselect);
                                        }
                                    }
                                    "totpoly" => {
                                        totpoly = get_int(member, &bytes_read, byte_index);
                                        if verbose {
                                            println!("  totpoly = {:?}", totpoly);
                                        }
                                    }
                                    "totloop" => {
                                        totloop = get_int(member, &bytes_read, byte_index);
                                        if verbose {
                                            println!("  totloop = {:?}", totloop);
                                        }
                                    }
                                    _ => {}
                                }
                                // find mem_type in dna_types.names
                                if let Some(type_found) = dna_types_hm.get(&member.mem_type) {
                                    let mem_tlen: u16 = calc_mem_tlen(member, *type_found);
                                    byte_index += mem_tlen as usize;
                                }
                            }
                            // reset booleans
                            // data_following_mesh = false;
                            data_following_material = false;
                            is_smooth = false;
                            // data
                            data_following_mesh = true;
                        }
                        "MPoly" => {
                            let mut loopstart: i32 = 0;
                            let mut totloop: i32 = 0;
                            for member in &struct_found.members {
                                match member.mem_name.as_str() {
                                    "loopstart" => {
                                        loopstart = get_int(member, &bytes_read, byte_index);
                                        if verbose {
                                            println!("  loopstart = {:?}", loopstart);
                                        }
                                    }
                                    "totloop" => {
                                        totloop = get_int(member, &bytes_read, byte_index);
                                        if verbose {
                                            println!("  totloop = {:?}", totloop);
                                        }
                                    }
                                    "flag" => {
                                        let flag: u8 = bytes_read[byte_index];
                                        if verbose {
                                            println!("  flag = {}", flag);
                                        }
                                        is_smooth = flag % 2 == 1;
                                        if verbose {
                                            println!("  is_smooth = {}", is_smooth);
                                        }
                                    }
                                    _ => {}
                                }
                                // find mem_type in dna_types.names
                                if let Some(type_found) = dna_types_hm.get(&member.mem_type) {
                                    let mem_tlen: u16 = calc_mem_tlen(member, *type_found);
                                    byte_index += mem_tlen as usize;
                                }
                            }
                            if totloop == 3_i32 {
                                loops.push(totloop as u8);
                                // triangle
                                for i in 0..3 {
                                    vertex_indices
                                        .push(loop_indices[(loopstart + i) as usize] as u32);
                                }
                            } else if totloop == 4_i32 {
                                loops.push(totloop as u8);
                                // quads
                                vertex_indices.push(loop_indices[(loopstart + 0) as usize] as u32);
                                vertex_indices.push(loop_indices[(loopstart + 1) as usize] as u32);
                                vertex_indices.push(loop_indices[(loopstart + 2) as usize] as u32);
                                vertex_indices.push(loop_indices[(loopstart + 0) as usize] as u32);
                                vertex_indices.push(loop_indices[(loopstart + 2) as usize] as u32);
                                vertex_indices.push(loop_indices[(loopstart + 3) as usize] as u32);
                            } else {
                                println!(
                                    "WARNING: quads or triangles expected (totloop = {}): {:?}",
                                    totloop, base_name
                                )
                            }
                        }
                        "bNode" => {
                            let mut next: usize;
                            let mut prev: usize;
                            let mut idname: String;
                            let mut inputs_first: usize;
                            let mut inputs_last: usize;
                            println!(
                                "{} (SDNAnr = {}) ({:#018x})",
                                struct_read,
                                pointers_read[struct_index].1,
                                pointers_read[struct_index].0
                            );
                            for member in &struct_found.members {
                                match member.mem_name.as_str() {
                                    "*next" => {
                                        next = get_pointer(member, &bytes_read, byte_index);
                                        if let Some(pointer_found) = dna_pointers_hm.get(&next) {
                                            println!(
                                                "next = {:#010x} ({:#018x})",
                                                pointer_found, next
                                            );
                                        }
                                    }
                                    "*prev" => {
                                        prev = get_pointer(member, &bytes_read, byte_index);
                                        if let Some(pointer_found) = dna_pointers_hm.get(&prev) {
                                            println!(
                                                "prev = {:#010x} ({:#018x})",
                                                pointer_found, prev
                                            );
                                        }
                                    }
                                    "idname[64]" => match member.mem_type.as_str() {
                                        "char" => {
                                            if let Some(type_found) =
                                                dna_types_hm.get(&member.mem_type)
                                            {
                                                let mem_tlen: u16 =
                                                    calc_mem_tlen(member, *type_found);
                                                idname = String::with_capacity(mem_tlen as usize);
                                                for i in 0..mem_tlen as usize {
                                                    if bytes_read[byte_index + i] == 0 {
                                                        break;
                                                    }
                                                    idname.push(bytes_read[byte_index + i] as char);
                                                }
                                                println!("idname[64] = {:?}", idname);
                                            }
                                        }
                                        _ => {}
                                    },
                                    "inputs" => {
                                        if let Some(struct_found2) =
                                            dna_structs_hm.get(&member.mem_type)
                                        {
                                            for member2 in &struct_found2.members {
                                                match member2.mem_name.as_str() {
                                                    "*first" => {
                                                        inputs_first = get_pointer(
                                                            member,
                                                            &bytes_read,
                                                            byte_index,
                                                        );
                                                        if let Some(pointer_found) =
                                                            dna_pointers_hm.get(&inputs_first)
                                                        {
                                                            println!(
                                                                "inputs_first = {:#010x} ({:#018x})",
                                                                pointer_found, inputs_first
                                                            );
                                                        }
                                                    }
                                                    "*last" => {
                                                        inputs_last = get_pointer(
                                                            member,
                                                            &bytes_read,
                                                            byte_index,
                                                        );
                                                        if let Some(pointer_found) =
                                                            dna_pointers_hm.get(&inputs_last)
                                                        {
                                                            println!(
                                                                "inputs_last = {:#010x} ({:#018x})",
                                                                pointer_found, inputs_last
                                                            );
                                                        }
                                                    }
                                                    _ => {}
                                                }
                                                // find mem_type in dna_types.names
                                                if let Some(type_found) =
                                                    dna_types_hm.get(&member2.mem_type)
                                                {
                                                    let mem_tlen: u16 =
                                                        calc_mem_tlen(member2, *type_found);
                                                    byte_index += mem_tlen as usize;
                                                }
                                            }
                                            // find mem_type in dna_types.names
                                            if let Some(type_found) =
                                                dna_types_hm.get(&member.mem_type)
                                            {
                                                let mem_tlen: u16 =
                                                    calc_mem_tlen(member, *type_found);
                                                // subtract (because it gets added below)
                                                byte_index -= mem_tlen as usize;
                                            }
                                        }
                                    }
                                    _ => {}
                                }
                                // find mem_type in dna_types.names
                                if let Some(type_found) = dna_types_hm.get(&member.mem_type) {
                                    let mem_tlen: u16 = calc_mem_tlen(member, *type_found);
                                    byte_index += mem_tlen as usize;
                                }
                            }
                        }
                        "bNodeTree" => {
                            let mut first: usize;
                            let mut last: usize;
                            println!(
                                "{} (SDNAnr = {}) ({:#018x})",
                                struct_read,
                                pointers_read[struct_index].1,
                                pointers_read[struct_index].0
                            );
                            for member in &struct_found.members {
                                match member.mem_name.as_str() {
                                    "nodes" => {
                                        if let Some(struct_found2) =
                                            dna_structs_hm.get(&member.mem_type)
                                        {
                                            for member2 in &struct_found2.members {
                                                match member2.mem_name.as_str() {
                                                    "*first" => {
                                                        first = get_pointer(
                                                            member,
                                                            &bytes_read,
                                                            byte_index,
                                                        );
                                                        if let Some(pointer_found) =
                                                            dna_pointers_hm.get(&first)
                                                        {
                                                            println!(
                                                                "first = {:#010x} ({:#018x})",
                                                                pointer_found, first
                                                            );
                                                        }
                                                    }
                                                    "*last" => {
                                                        last = get_pointer(
                                                            member,
                                                            &bytes_read,
                                                            byte_index,
                                                        );
                                                        if let Some(pointer_found) =
                                                            dna_pointers_hm.get(&last)
                                                        {
                                                            println!(
                                                                "last = {:#010x} ({:#018x})",
                                                                pointer_found, last
                                                            );
                                                        }
                                                    }
                                                    _ => {}
                                                }
                                                // find mem_type in dna_types.names
                                                if let Some(type_found) =
                                                    dna_types_hm.get(&member2.mem_type)
                                                {
                                                    let mem_tlen: u16 =
                                                        calc_mem_tlen(member2, *type_found);
                                                    byte_index += mem_tlen as usize;
                                                }
                                            }
                                            // find mem_type in dna_types.names
                                            if let Some(type_found) =
                                                dna_types_hm.get(&member.mem_type)
                                            {
                                                let mem_tlen: u16 =
                                                    calc_mem_tlen(member, *type_found);
                                                // subtract (because it gets added below)
                                                byte_index -= mem_tlen as usize;
                                            }
                                        }
                                    }
                                    _ => {}
                                }
                                // find mem_type in dna_types.names
                                if let Some(type_found) = dna_types_hm.get(&member.mem_type) {
                                    let mem_tlen: u16 = calc_mem_tlen(member, *type_found);
                                    byte_index += mem_tlen as usize;
                                }
                            }
                        }
                        "bNodeSocket" => {
                            let mut idname: String;
                            let mut default_value: usize;
                            println!(
                                "{} (SDNAnr = {}) ({:#018x})",
                                struct_read,
                                pointers_read[struct_index].1,
                                pointers_read[struct_index].0
                            );
                            for member in &struct_found.members {
                                match member.mem_name.as_str() {
                                    "idname[64]" => match member.mem_type.as_str() {
                                        "char" => {
                                            if let Some(type_found) =
                                                dna_types_hm.get(&member.mem_type)
                                            {
                                                let mem_tlen: u16 =
                                                    calc_mem_tlen(member, *type_found);
                                                idname = String::with_capacity(mem_tlen as usize);
                                                for i in 0..mem_tlen as usize {
                                                    if bytes_read[byte_index + i] == 0 {
                                                        break;
                                                    }
                                                    idname.push(bytes_read[byte_index + i] as char);
                                                }
                                                println!("idname[64] = {:?}", idname);
                                            }
                                        }
                                        _ => {}
                                    },
                                    "*default_value" => {
                                        default_value =
                                            get_pointer(member, &bytes_read, byte_index);
                                        if let Some(pointer_found) =
                                            dna_pointers_hm.get(&default_value)
                                        {
                                            println!(
                                                "default_value = {:#010x} ({:#018x})",
                                                pointer_found, default_value
                                            );
                                            if data_following_material
                                                && search_for_emit
                                                && default_value != 0
                                            {
                                                emit_default_value = default_value;
                                                println!(
                                                    "emit_default_value = {:#018x}",
                                                    emit_default_value
                                                );
                                            }
                                        }
                                    }
                                    _ => {}
                                }
                                // find mem_type in dna_types.names
                                if let Some(type_found) = dna_types_hm.get(&member.mem_type) {
                                    let mem_tlen: u16 = calc_mem_tlen(member, *type_found);
                                    byte_index += mem_tlen as usize;
                                }
                            }
                        }
                        "bNodeSocketValueFloat" => {
                            let mut value: f32;
                            for member in &struct_found.members {
                                match member.mem_name.as_str() {
                                    "value" => {
                                        value = get_float(member, &bytes_read, byte_index);
                                        println!(
                                            "{} (SDNAnr = {}) ({:#018x})",
                                            struct_read,
                                            pointers_read[struct_index].1,
                                            pointers_read[struct_index].0
                                        );
                                        println!("{} = {:?}", member.mem_name, value);
                                        if data_following_material
                                            && search_for_emit
                                            && emit_default_value == pointers_read[struct_index].0
                                        {
                                            emit = value;
                                            println!("emit = {:?}", emit);
                                        }
                                    }
                                    _ => {}
                                }
                                // find mem_type in dna_types.names
                                if let Some(type_found) = dna_types_hm.get(&member.mem_type) {
                                    let mem_tlen: u16 = calc_mem_tlen(member, *type_found);
                                    byte_index += mem_tlen as usize;
                                }
                            }
                        }
                        _ => {
                            byte_index += data_read[struct_index] as usize;
                        }
                    }
                }
            } else {
                // several structs (from DATA chunk)
                let num_structs: u32 = data_read[struct_index] / (*tlen as u32);
                if let Some(struct_found) = dna_structs_hm.get(&struct_read) {
                    match struct_read.as_str() {
                        "MVert" => {
                            for s in 0..num_structs {
                                for member in &struct_found.members {
                                    match member.mem_name.as_str() {
                                        "co[3]" => {
                                            let co: [f32; 3] =
                                                get_float3(member, &bytes_read, byte_index);
                                            if verbose {
                                                println!("  co[{}] = {:?}", s, co);
                                            }
                                            p.push(Point3f {
                                                x: (co[0] * scale_length) as Float,
                                                y: (co[1] * scale_length) as Float,
                                                z: (co[2] * scale_length) as Float,
                                            });
                                        }
                                        "no[3]" => {
                                            let no: [i16; 3] =
                                                get_short3(member, &bytes_read, byte_index);
                                            // convert short values to floats
                                            let factor: f32 = 1.0 / 32767.0;
                                            let mut nof: [f32; 3] = [0.0; 3];
                                            for i in 0..3 {
                                                nof[i] = no[i] as f32 * factor;
                                            }
                                            if verbose {
                                                println!("  no[{}] = {:?}", s, nof);
                                            }
                                            n.push(Normal3f {
                                                x: nof[0] as Float,
                                                y: nof[1] as Float,
                                                z: nof[2] as Float,
                                            });
                                        }
                                        _ => {}
                                    }
                                    // find mem_type in dna_types.names
                                    if let Some(type_found) = dna_types_hm.get(&member.mem_type) {
                                        let mem_tlen: u16 = calc_mem_tlen(member, *type_found);
                                        byte_index += mem_tlen as usize;
                                    }
                                }
                            }
                        }
                        "MLoop" => {
                            for s in 0..num_structs {
                                for member in &struct_found.members {
                                    match member.mem_name.as_str() {
                                        "v" => {
                                            let v: i32 = get_int(member, &bytes_read, byte_index);
                                            if verbose {
                                                println!("  v[{}] = {:?}", s, v);
                                            }
                                            loop_indices.push(v);
                                        }
                                        "e" => {
                                            let e: i32 = get_int(member, &bytes_read, byte_index);
                                            if verbose {
                                                println!("  e[{}] = {:?}", s, e);
                                            }
                                        }
                                        _ => {}
                                    }
                                    // find mem_type in dna_types.names
                                    if let Some(type_found) = dna_types_hm.get(&member.mem_type) {
                                        let mem_tlen: u16 = calc_mem_tlen(member, *type_found);
                                        byte_index += mem_tlen as usize;
                                    }
                                }
                            }
                        }
                        "MLoopUV" => {
                            for s in 0..num_structs {
                                for member in &struct_found.members {
                                    match member.mem_name.as_str() {
                                        "uv[2]" => {
                                            let uv: [f32; 2] =
                                                get_float2(member, &bytes_read, byte_index);
                                            if verbose {
                                                println!("  uv[{}] = {:?}", s, uv);
                                            }
                                            uvs.push(Point2f {
                                                x: uv[0] as Float,
                                                y: uv[1] as Float,
                                            });
                                        }
                                        _ => {}
                                    }
                                    // find mem_type in dna_types.names
                                    if let Some(type_found) = dna_types_hm.get(&member.mem_type) {
                                        let mem_tlen: u16 = calc_mem_tlen(member, *type_found);
                                        byte_index += mem_tlen as usize;
                                    }
                                }
                            }
                        }
                        "MPoly" => {
                            for s in 0..num_structs {
                                let mut loopstart: i32 = 0;
                                let mut totloop: i32 = 0;
                                for member in &struct_found.members {
                                    match member.mem_name.as_str() {
                                        "loopstart" => {
                                            loopstart = get_int(member, &bytes_read, byte_index);
                                            if verbose {
                                                println!("  loopstart[{}] = {:?}", s, loopstart);
                                            }
                                        }
                                        "totloop" => {
                                            totloop = get_int(member, &bytes_read, byte_index);
                                            if verbose {
                                                println!("  totloop[{}] = {:?}", s, totloop);
                                            }
                                        }
                                        "flag" => {
                                            let flag: u8 = bytes_read[byte_index];
                                            if verbose {
                                                println!("  flag[{}] = {}", s, flag);
                                            }
                                            is_smooth = flag % 2 == 1;
                                            if verbose {
                                                println!("  is_smooth[{}] = {}", s, is_smooth);
                                            }
                                        }
                                        _ => {}
                                    }
                                    // find mem_type in dna_types.names
                                    if let Some(type_found) = dna_types_hm.get(&member.mem_type) {
                                        let mem_tlen: u16 = calc_mem_tlen(member, *type_found);
                                        byte_index += mem_tlen as usize;
                                    }
                                }
                                if totloop == 3_i32 {
                                    loops.push(totloop as u8);
                                    // triangle
                                    for i in 0..3 {
                                        vertex_indices
                                            .push(loop_indices[(loopstart + i) as usize] as u32);
                                    }
                                } else if totloop == 4_i32 {
                                    loops.push(totloop as u8);
                                    // quads
                                    vertex_indices
                                        .push(loop_indices[(loopstart + 0) as usize] as u32);
                                    vertex_indices
                                        .push(loop_indices[(loopstart + 1) as usize] as u32);
                                    vertex_indices
                                        .push(loop_indices[(loopstart + 2) as usize] as u32);
                                    vertex_indices
                                        .push(loop_indices[(loopstart + 0) as usize] as u32);
                                    vertex_indices
                                        .push(loop_indices[(loopstart + 2) as usize] as u32);
                                    vertex_indices
                                        .push(loop_indices[(loopstart + 3) as usize] as u32);
                                } else {
                                    println!(
                                        "WARNING: quads or triangles expected (totloop = {}): {:?}",
                                        totloop, base_name
                                    )
                                }
                            }
                        }
                        _ => {
                            println!("{} * {}", num_structs, struct_read);
                            byte_index += data_read[struct_index] as usize;
                        }
                    }
                } else {
                }
            }
        } else {
            byte_index += data_read[struct_index] as usize;
        }
        struct_index += 1;
    }
    if verbose {
        println!("byte_index = {}", byte_index);
    }
    // use HDR image if one was found
    if !hdr_path.is_empty() {
        let axis: Vector3f = Vector3f {
            x: 0.0 as Float,
            y: 0.0 as Float,
            z: 1.0 as Float,
        };
        let light_to_world: Transform = Transform::rotate(180.0 as Float, &axis);
        builder.add_hdr_light(
            light_to_world,
            String::from(hdr_path.to_str().unwrap()),
            args.light_scale,
        );
    }
    let scene_description: SceneDescription = builder.finalize();
    let mut render_options: RenderOptions = RenderOptions::new(
        scene_description,
        &material_hm,
        &texture_hm,
        args.light_scale as Float,
    );
    assert!(render_options.shapes.len() == render_options.shape_lights.len());
    for shape_idx in 0..render_options.shapes.len() {
        let shape = &render_options.shapes[shape_idx];
        let shape_material = &render_options.shape_materials[shape_idx];
        let shape_light = &render_options.shape_lights[shape_idx];
        let geo_prim = Arc::new(Primitive::Geometric(Box::new(GeometricPrimitive::new(
            shape.clone(),
            Some(shape_material.clone()),
            shape_light.clone(),
            None,
        ))));
        render_options.primitives.push(geo_prim.clone());
    }
    println!("number of lights = {:?}", render_options.lights.len());
    println!(
        "number of primitives = {:?}",
        render_options.primitives.len()
    );
    let mut pos = Point3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };
    let mut look = Point3f {
        x: 1.0,
        y: 0.0,
        z: 0.0,
    };
    let mut up = Vector3f {
        x: 0.0,
        y: 0.0,
        z: 1.0,
    };
    base_name = String::from("Camera");
    if let Some(camera_name) = args.camera_name {
        base_name = camera_name.clone();
    }
    if let Some(o2w) = object_to_world_hm.get(&base_name) {
        pos = Point3f {
            x: o2w.m.m[0][3],
            y: o2w.m.m[1][3],
            z: o2w.m.m[2][3],
        };
        let forwards: Point3f = Point3f {
            x: -o2w.m.m[0][2] * scale_length,
            y: -o2w.m.m[1][2] * scale_length,
            z: -o2w.m.m[2][2] * scale_length,
        };
        look = pos + forwards;
        up = Vector3f {
            x: o2w.m.m[0][1],
            y: o2w.m.m[1][1],
            z: o2w.m.m[2][1],
        };
    } else {
        println!(
            "WARNING: looking up object_to_world by name ({:?}) failed",
            base_name
        );
    }
    let t: Transform = Transform::scale(-1.0, 1.0, 1.0) * Transform::look_at(&pos, &look, &up);
    let it: Transform = Transform {
        m: t.m_inv.clone(),
        m_inv: t.m.clone(),
    };
    let animated_cam_to_world: AnimatedTransform = AnimatedTransform::new(&it, 0.0, &it, 1.0);
    let aspect: Float = resolution_x as Float / resolution_y as Float;
    let mut fov: Float;
    let mut clipsta: Float = 0.0;
    if aspect > 1.0 {
        fov = angle_y;
    } else {
        fov = angle_x;
    }
    if let Some(cam) = camera_hm.get(&base_name) {
        // overwrite fov
        if aspect > 1.0 {
            // fov = angle_x / 2.0;
            fov = 2.0 as Float * degrees((16.0 as Float / (aspect * cam.lens)).atan());
        } else {
            fov = 2.0 as Float * degrees(((aspect * 16.0 as Float) / cam.lens).atan());
        }
        clipsta = cam.clipsta;
        // println!("fov[{}] overwritten", fov);
        // println!("clipsta[{}] overwritten", clipsta);
    }
    let frame: Float = resolution_x as Float / resolution_y as Float;
    let mut screen: Bounds2f = Bounds2f::default();
    if frame > 1.0 {
        screen.p_min.x = -frame;
        screen.p_max.x = frame;
        screen.p_min.y = -1.0;
        screen.p_max.y = 1.0;
    } else {
        screen.p_min.x = -1.0;
        screen.p_max.x = 1.0;
        screen.p_min.y = -1.0 / frame;
        screen.p_max.y = 1.0 / frame;
    }
    let render_x: u32 = resolution_x * resolution_percentage as u32 / 100_u32;
    let render_y: u32 = resolution_y * resolution_percentage as u32 / 100_u32;
    println!(
        "{}x{} [{}%] = {}x{}",
        resolution_x, resolution_y, resolution_percentage, render_x, render_y
    );
    if let Some(integrator_name) = args.integrator {
        let mut integrator_params: ParamSet = ParamSet::default();
        if integrator_name == "mlt" {
            integrator_params.add_int(String::from("maxdepth"), args.max_depth as i32);
            // MLT
            integrator_params.add_int(
                String::from("bootstrapsamples"),
                args.bootstrap_samples as i32,
            );
            integrator_params.add_int(String::from("chains"), args.chains as i32);
            integrator_params.add_int(
                String::from("mutationsperpixel"),
                args.mutations_per_pixel as i32,
            );
            integrator_params.add_float(
                String::from("largestepprobability"),
                args.step_probability as Float,
            );
            integrator_params.add_float(String::from("sigma"), args.sigma as Float);
        } else if integrator_name == "sppm" {
            integrator_params.add_int(String::from("maxdepth"), args.max_depth as i32);
            // SPPM
            integrator_params.add_int(String::from("imagewritefrequency"), args.write_frequency);
        } else {
            integrator_params.add_int(String::from("maxdepth"), args.max_depth as i32);
        }
        let some_integrator: Option<Box<Integrator>> = make_integrator(
            &integrator_name,
            2.0 as Float,
            render_x as i32,
            render_y as i32,
            fov,
            clipsta,
            animated_cam_to_world,
            args.samples as i32,
            integrator_params,
        );
        if let Some(mut integrator) = some_integrator {
            let scene = make_scene(&render_options.primitives, render_options.lights);
            let num_threads: u8 = num_cpus::get() as u8;
            integrator.render(&scene, num_threads);
        } else {
            panic!("Unable to create integrator.");
        }
    } else {
        let mut integrator_params: ParamSet = ParamSet::default();
        integrator_params.add_int(String::from("maxdepth"), args.max_depth as i32);
        let integrator_name: String;
        if render_options.has_emitters || render_options.lights.len() > 0 {
            integrator_name = String::from("path");
        } else {
            integrator_name = String::from("ao");
        }
        let some_integrator: Option<Box<Integrator>> = make_integrator(
            &integrator_name,
            2.0 as Float,
            render_x as i32,
            render_y as i32,
            fov,
            clipsta,
            animated_cam_to_world,
            args.samples as i32,
            integrator_params,
        );
        if let Some(mut integrator) = some_integrator {
            let scene = make_scene(&render_options.primitives, render_options.lights);
            let num_threads: u8 = num_cpus::get() as u8;
            integrator.render(&scene, num_threads);
        } else {
            panic!("Unable to create integrator.");
        }
    }
    Ok(())
}
