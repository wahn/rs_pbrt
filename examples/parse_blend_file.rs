// Assumptions:
// 1. objects, data, and materials have the same base name:
//    e.g. OBcornellbox, MEcornellbox and MAcornellbox
//    TODO: Find a better solution (following the Object->data pointer?)
// 2. Right now we search for "Camera" for Transform::look_at(...)
//    TODO: get the render camera/transform from the scene (self.scene.camera)
// 3. Smoothness is valid for the whole mesh
//    TODO: store smoothness per polygon, split mesh into smooth and non-smooth parts

extern crate num_cpus;
extern crate pbrt;
extern crate structopt;

// std
use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::mem;
use std::sync::Arc;
use structopt::StructOpt;
// pbrt
use pbrt::accelerators::bvh::{BVHAccel, SplitMethod};
use pbrt::cameras::perspective::PerspectiveCamera;
use pbrt::core::camera::Camera;
use pbrt::core::film::Film;
use pbrt::core::filter::Filter;
use pbrt::core::geometry::{
    Bounds2f, Bounds2i, Normal3f, Point2f, Point2i, Point3f, Vector2f, Vector3f,
};
use pbrt::core::integrator::SamplerIntegrator;
use pbrt::core::light::{AreaLight, Light};
use pbrt::core::material::Material;
use pbrt::core::medium::MediumInterface;
use pbrt::core::pbrt::degrees;
use pbrt::core::pbrt::{Float, Spectrum};
use pbrt::core::primitive::{GeometricPrimitive, Primitive};
use pbrt::core::sampler::Sampler;
use pbrt::core::scene::Scene;
use pbrt::core::shape::Shape;
use pbrt::core::transform::{AnimatedTransform, Transform};
use pbrt::filters::boxfilter::BoxFilter;
use pbrt::integrators::ao::AOIntegrator;
use pbrt::integrators::bdpt::render_bdpt;
use pbrt::integrators::bdpt::BDPTIntegrator;
use pbrt::integrators::directlighting::{DirectLightingIntegrator, LightStrategy};
use pbrt::integrators::mlt::render_mlt;
use pbrt::integrators::mlt::MLTIntegrator;
use pbrt::integrators::path::PathIntegrator;
use pbrt::integrators::render;
use pbrt::lights::diffuse::DiffuseAreaLight;
use pbrt::materials::matte::MatteMaterial;
use pbrt::samplers::zerotwosequence::ZeroTwoSequenceSampler;
use pbrt::shapes::triangle::{Triangle, TriangleMesh};
use pbrt::textures::constant::ConstantTexture;

pub const VERSION: &'static str = env!("CARGO_PKG_VERSION");

/// Parse a Blender scene file and render it.
#[derive(StructOpt)]
struct Cli {
    /// pixel samples
    #[structopt(short = "s", long = "samples", default_value = "1")]
    samples: u32,
    /// ao, directlighting, path, bdpt, mlt
    #[structopt(short = "i", long = "integrator")]
    integrator: Option<String>,
    /// The path to the file to read
    #[structopt(parse(from_os_str))]
    path: std::path::PathBuf,
}

// Blender

#[derive(Debug, Default, Copy, Clone)]
struct Blend279Material {
    pub r: f32,
    pub g: f32,
    pub b: f32,
    pub a: f32,
    pub emit: f32,
}

fn focallength_to_fov(focal_length: f32, sensor: f32) -> f32 {
    2.0_f32 * ((sensor / 2.0_f32) / focal_length).atan()
}

// TMP (see pbrt_spheres_differentials_texfilt.rs)

struct SceneDescription {
    mesh_names: Vec<String>,
    meshes: Vec<Arc<TriangleMesh>>,
    lights: Vec<Arc<Light + Sync + Send>>,
}

struct SceneDescriptionBuilder {
    mesh_names: Vec<String>,
    meshes: Vec<Arc<TriangleMesh>>,
    lights: Vec<Arc<Light + Sync + Send>>,
}

impl SceneDescriptionBuilder {
    fn new() -> SceneDescriptionBuilder {
        SceneDescriptionBuilder {
            mesh_names: Vec::new(),
            meshes: Vec::new(),
            lights: Vec::new(),
        }
    }
    fn add_mesh(
        &mut self,
        base_name: String,
        object_to_world: Transform,
        world_to_object: Transform,
        n_triangles: usize,
        vertex_indices: Vec<usize>,
        n_vertices: usize,
        p_ws: Vec<Point3f>,
        s: Vec<Vector3f>,
        n_ws: Vec<Normal3f>,
        uv: Vec<Point2f>,
    ) -> &mut SceneDescriptionBuilder {
        self.mesh_names.push(base_name);
        let triangle_mesh = Arc::new(TriangleMesh::new(
            object_to_world,
            world_to_object,
            false,
            n_triangles,
            vertex_indices,
            n_vertices,
            p_ws, // in world space
            s,    // empty
            n_ws, // in world space
            uv,   // empty
        ));
        self.meshes.push(triangle_mesh);
        self
    }
    fn finalize(self) -> SceneDescription {
        SceneDescription {
            mesh_names: self.mesh_names,
            meshes: self.meshes,
            lights: self.lights,
        }
    }
}

struct RenderOptions {
    has_emitters: bool,
    primitives: Vec<Arc<Primitive + Sync + Send>>,
    triangles: Vec<Arc<Triangle>>,
    triangle_materials: Vec<Arc<Material + Sync + Send>>,
    triangle_lights: Vec<Option<Arc<AreaLight + Sync + Send>>>,
    lights: Vec<Arc<Light + Sync + Send>>,
}

impl RenderOptions {
    fn new(
        scene: SceneDescription,
        material_hm: &HashMap<String, Blend279Material>,
    ) -> RenderOptions {
        let mut has_emitters: bool = false;
        let primitives: Vec<Arc<Primitive + Sync + Send>> = Vec::new();
        let mut triangles: Vec<Arc<Triangle>> = Vec::new();
        let mut triangle_materials: Vec<Arc<Material + Sync + Send>> = Vec::new();
        let mut triangle_lights: Vec<Option<Arc<AreaLight + Sync + Send>>> = Vec::new();
        let mut lights: Vec<Arc<Light + Sync + Send>> = Vec::new();
        // default material
        let kd = Arc::new(ConstantTexture::new(Spectrum::new(1.0)));
        let sigma = Arc::new(ConstantTexture::new(0.0 as Float));
        let default_material = Arc::new(MatteMaterial::new(kd, sigma, None));
        // lights
        for light in &scene.lights {
            lights.push(light.clone());
        }
        // meshes
        for mesh_idx in 0..scene.meshes.len() {
            let mesh = &scene.meshes[mesh_idx];
            let mesh_name = &scene.mesh_names[mesh_idx];
            // create individual triangles
            let mut shapes: Vec<Arc<Shape + Send + Sync>> = Vec::new();
            for id in 0..mesh.n_triangles {
                let triangle = Arc::new(Triangle::new(
                    mesh.object_to_world,
                    mesh.world_to_object,
                    mesh.transform_swaps_handedness,
                    mesh.clone(),
                    id,
                ));
                triangles.push(triangle.clone());
                shapes.push(triangle.clone());
            }
            if let Some(mat) = material_hm.get(mesh_name) {
                // println!("{:?}: {:?}", mesh_name, mat);
                if mat.emit > 0.0 {
                    has_emitters = true;
                    for i in 0..shapes.len() {
                        let shape = &shapes[i];
                        let mi: MediumInterface = MediumInterface::default();
                        let l_emit: Spectrum =
                            Spectrum::rgb(mat.r * mat.emit, mat.g * mat.emit, mat.b * mat.emit);
                        let n_samples: i32 = 1;
                        let two_sided: bool = false;
                        let area_light: Arc<DiffuseAreaLight> = Arc::new(DiffuseAreaLight::new(
                            &mesh.object_to_world,
                            &mi,
                            &l_emit,
                            n_samples,
                            shape.clone(),
                            two_sided,
                        ));
                        lights.push(area_light.clone());
                        triangle_materials.push(default_material.clone());
                        let triangle_light: Option<Arc<AreaLight + Sync + Send>> =
                            Some(area_light.clone());
                        triangle_lights.push(triangle_light);
                    }
                } else {
                    for _i in 0..shapes.len() {
                        let kd = Arc::new(ConstantTexture::new(Spectrum::rgb(mat.r, mat.g, mat.b)));
                        let sigma = Arc::new(ConstantTexture::new(0.0 as Float));
                        let matte = Arc::new(MatteMaterial::new(kd, sigma, None));
                        triangle_materials.push(matte);
                        triangle_lights.push(None);
                    }
                }
            } else {
                // println!("{:?}: no mat", mesh_name);
                for _i in 0..shapes.len() {
                    triangle_materials.push(default_material.clone());
                    triangle_lights.push(None);
                }
            }
        }
        RenderOptions {
            has_emitters: has_emitters,
            primitives: primitives,
            triangles: triangles,
            triangle_materials: triangle_materials,
            triangle_lights: triangle_lights,
            lights: lights,
        }
    }
}

// TMP

fn decode_blender_header(header: &[u8], version: &mut u32, print_it: bool) -> bool {
    // BLENDER
    match header[0] as char {
        'B' => {
            if print_it {
                print!("B");
            }
        }
        _ => return false,
    }
    match header[1] as char {
        'L' => {
            if print_it {
                print!("L");
            }
        }
        _ => return false,
    }
    match header[2] as char {
        'E' => {
            if print_it {
                print!("E");
            }
        }
        _ => return false,
    }
    match header[3] as char {
        'N' => {
            if print_it {
                print!("N");
            }
        }
        _ => return false,
    }
    match header[4] as char {
        'D' => {
            if print_it {
                print!("D");
            }
        }
        _ => return false,
    }
    match header[5] as char {
        'E' => {
            if print_it {
                print!("E");
            }
        }
        _ => return false,
    }
    match header[6] as char {
        'R' => {
            if print_it {
                print!("R");
            }
        }
        _ => return false,
    }
    // [_|-]
    match header[7] as char {
        '_' => {
            if print_it {
                print!("_");
            }
        }
        '-' => {
            if print_it {
                print!("-");
            }
        }
        _ => return false,
    }
    // [v|V]
    match header[8] as char {
        'v' => {
            if print_it {
                print!("v");
            }
        }
        'V' => {
            if print_it {
                print!("V");
            }
        }
        _ => return false,
    }
    for i in 9..12 {
        if header[i].is_ascii_digit() {
            if print_it {
                print!("{:?}", (header[i] as char).to_digit(10).unwrap());
            }
        } else {
            return false;
        }
    }
    if print_it {
        print!("\n");
    }
    // get the version number
    let last3c = vec![header[9], header[10], header[11]];
    let version_str = String::from_utf8(last3c).unwrap();
    // convert to u32 and return
    *version = version_str.parse::<u32>().unwrap();
    true
}

fn make_id(code: &[u8]) -> String {
    let mut id = String::with_capacity(4);
    for i in 0..4 {
        if (code[i] as char).is_ascii_alphanumeric() {
            id.push(code[i] as char);
        }
    }
    id
}

fn read_names(
    f: &mut File,
    nr_names: usize,
    names: &mut Vec<String>,
    byte_counter: &mut usize,
) -> std::io::Result<()> {
    // let mut name_counter: usize = 0;
    let mut buffer = [0; 1];
    loop {
        if names.len() == nr_names {
            break;
        } else {
            let mut name = String::new();
            loop {
                // read only one char/byte
                f.read(&mut buffer)?;
                *byte_counter += 1;
                if buffer[0] == 0 {
                    break;
                } else {
                    name.push(buffer[0] as char);
                }
            }
            // println!("  {:?}", name);
            names.push(name);
            // name_counter += 1;
        }
    }
    // println!("  {} names found in {} bytes", name_counter, byte_counter);
    Ok(())
}

fn read_type_names(
    f: &mut File,
    nr_types: usize,
    type_names: &mut Vec<String>,
    byte_counter: &mut usize,
) -> std::io::Result<()> {
    // let mut name_counter: usize = 0;
    let mut buffer = [0; 1];
    loop {
        if type_names.len() == nr_types {
            break;
        } else {
            let mut name = String::new();
            loop {
                // read only one char/byte
                f.read(&mut buffer)?;
                *byte_counter += 1;
                if buffer[0] == 0 {
                    break;
                } else {
                    name.push(buffer[0] as char);
                }
            }
            // println!("  {:?}", name);
            type_names.push(name);
            // name_counter += 1;
        }
    }
    // println!(
    //     "  {} type names found in {} bytes",
    //     name_counter, byte_counter
    // );
    Ok(())
}

fn main() -> std::io::Result<()> {
    let args = Cli::from_args();
    let num_threads: u8 = num_cpus::get() as u8;
    println!(
        "parse_blend_file version {} [Detected {} cores]",
        VERSION, num_threads
    );
    // PBRT
    let mut scale_length: f32 = 1.0;
    let mut resolution_x: u32 = 640;
    let mut resolution_y: u32 = 480;
    let mut angle_x: f32 = 45.0;
    let mut angle_y: f32 = 45.0;
    let mut base_name = String::new();
    let mut material_hm: HashMap<String, Blend279Material> = HashMap::new();
    let mut object_to_world_hm: HashMap<String, Transform> = HashMap::new();
    let mut object_to_world: Transform = Transform::new(
        1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    );
    let mut p: Vec<Point3f> = Vec::new();
    let mut n: Vec<Normal3f> = Vec::new();
    let mut vertex_indices: Vec<usize> = Vec::new();
    // TODO: let mut primitives: Vec<Arc<Primitive + Sync + Send>> = Vec::new();
    // TODO: let mut lights: Vec<Arc<Light + Sync + Send>> = Vec::new();
    // first get the DNA
    let mut names: Vec<String> = Vec::new();
    let mut names_len: usize = 0;
    let mut types: Vec<String> = Vec::new();
    let mut dna_2_type_id: Vec<u16> = Vec::new();
    let mut types_len: usize = 0;
    let mut tlen: Vec<u16> = Vec::new();
    {
        let mut f = File::open(&args.path)?;
        // read exactly 12 bytes
        // let mut counter: usize = 0;
        let mut buffer = [0; 12];
        f.read(&mut buffer)?;
        // counter += 12;
        let mut blender_version: u32 = 0;
        if !decode_blender_header(&buffer, &mut blender_version, true) {
            println!("ERROR: Not a .blend file");
            println!("First 12 bytes:");
            println!("{:?}", buffer);
        } else {
            loop {
                // code
                let mut buffer = [0; 4];
                f.read(&mut buffer)?;
                // counter += 4;
                let code = make_id(&buffer);
                // len
                let mut buffer = [0; 4];
                f.read(&mut buffer)?;
                // counter += 4;
                let mut len: u32 = 0;
                len += (buffer[0] as u32) << 0;
                len += (buffer[1] as u32) << 8;
                len += (buffer[2] as u32) << 16;
                len += (buffer[3] as u32) << 24;
                // for now ignore the old entry
                let mut buffer = [0; 8];
                f.read(&mut buffer)?;
                // counter += 8;
                // get SDNAnr
                let mut buffer = [0; 4];
                f.read(&mut buffer)?;
                // counter += 4;
                // for now ignore the nr entry
                let mut buffer = [0; 4];
                f.read(&mut buffer)?;
                // counter += 4;
                // are we done?
                if code == String::from("ENDB") {
                    break;
                }
                if code == String::from("DNA1") {
                    // println!("{} ({})", code, len);
                    // "SDNANAME" in first 8 bytes
                    let mut buffer = [0; 8];
                    f.read(&mut buffer)?;
                    // counter += 8;
                    let mut sdna_name = String::with_capacity(8);
                    for i in 0..8 {
                        if (buffer[i] as char).is_ascii_alphabetic() {
                            sdna_name.push(buffer[i] as char);
                        }
                    }
                    if sdna_name != String::from("SDNANAME") {
                        // read remaining bytes
                        let mut buffer = vec![0; (len - 8) as usize];
                        f.read(&mut buffer)?;
                    // counter += (len - 8) as usize;
                    } else {
                        let mut buffer = [0; 4];
                        f.read(&mut buffer)?;
                        // counter += 4;
                        let mut nr_names: u32 = 0;
                        nr_names += (buffer[0] as u32) << 0;
                        nr_names += (buffer[1] as u32) << 8;
                        nr_names += (buffer[2] as u32) << 16;
                        nr_names += (buffer[3] as u32) << 24;
                        read_names(&mut f, nr_names as usize, &mut names, &mut names_len)?;
                        // counter += names_len;
                        let mut remaining_bytes: usize = (len - 12) as usize - names_len;
                        // skip pad bytes, read "TYPE" and nr_types
                        let mut buffer = [0; 1];
                        loop {
                            f.read(&mut buffer)?;
                            // counter += 1;
                            if buffer[0] == 0 {
                                // skip pad byte
                                remaining_bytes -= 1;
                            } else if buffer[0] as char == 'T' {
                                remaining_bytes -= 1;
                                break;
                            }
                        }
                        // match 'YPE' ('T' was matched above)
                        let mut buffer = [0; 3];
                        f.read(&mut buffer)?;
                        // counter += 3;
                        remaining_bytes -= 3;
                        if buffer[0] as char == 'Y'
                            && buffer[1] as char == 'P'
                            && buffer[2] as char == 'E'
                        {
                            // nr_types
                            let mut buffer = [0; 4];
                            f.read(&mut buffer)?;
                            // counter += 4;
                            remaining_bytes -= 4;
                            let mut nr_types: u32 = 0;
                            nr_types += (buffer[0] as u32) << 0;
                            nr_types += (buffer[1] as u32) << 8;
                            nr_types += (buffer[2] as u32) << 16;
                            nr_types += (buffer[3] as u32) << 24;
                            read_type_names(&mut f, nr_types as usize, &mut types, &mut types_len)?;
                            // counter += types_len;
                            remaining_bytes -= types_len;
                            // skip pad bytes, read "TLEN"
                            let mut buffer = [0; 1];
                            loop {
                                f.read(&mut buffer)?;
                                // counter += 1;
                                if buffer[0] == 0 {
                                    // skip pad byte
                                    remaining_bytes -= 1;
                                } else if buffer[0] as char == 'T' {
                                    remaining_bytes -= 1;
                                    break;
                                }
                            }
                            // match 'LEN' ('T' was matched above)
                            let mut buffer = [0; 3];
                            f.read(&mut buffer)?;
                            // counter += 3;
                            remaining_bytes -= 3;
                            if buffer[0] as char == 'L'
                                && buffer[1] as char == 'E'
                                && buffer[2] as char == 'N'
                            {
                                // read short (16 bits = 2 bytes) for each type
                                for _i in 0..nr_types as usize {
                                    let mut buffer = [0; 2];
                                    f.read(&mut buffer)?;
                                    // counter += 2;
                                    remaining_bytes -= 2;
                                    let mut type_size: u16 = 0;
                                    type_size += (buffer[0] as u16) << 0;
                                    type_size += (buffer[1] as u16) << 8;
                                    tlen.push(type_size);
                                }
                                // skip pad bytes, read "STRC"
                                let mut buffer = [0; 1];
                                loop {
                                    f.read(&mut buffer)?;
                                    // counter += 1;
                                    if buffer[0] == 0 {
                                        // skip pad byte
                                        remaining_bytes -= 1;
                                    } else if buffer[0] as char == 'S' {
                                        remaining_bytes -= 1;
                                        break;
                                    }
                                }
                                // match 'TRC' ('S' was matched above)
                                let mut buffer = [0; 3];
                                f.read(&mut buffer)?;
                                // counter += 3;
                                remaining_bytes -= 3;
                                if buffer[0] as char == 'T'
                                    && buffer[1] as char == 'R'
                                    && buffer[2] as char == 'C'
                                {
                                    // nr_structs
                                    let mut buffer = [0; 4];
                                    f.read(&mut buffer)?;
                                    // counter += 4;
                                    remaining_bytes -= 4;
                                    let mut nr_structs: u32 = 0;
                                    nr_structs += (buffer[0] as u32) << 0;
                                    nr_structs += (buffer[1] as u32) << 8;
                                    nr_structs += (buffer[2] as u32) << 16;
                                    nr_structs += (buffer[3] as u32) << 24;
                                    for _s in 0..nr_structs as usize {
                                        // read two short values
                                        let mut buffer = [0; 2];
                                        f.read(&mut buffer)?;
                                        // counter += 2;
                                        remaining_bytes -= 2;
                                        let mut type_idx: u16 = 0;
                                        type_idx += (buffer[0] as u16) << 0;
                                        type_idx += (buffer[1] as u16) << 8;
                                        f.read(&mut buffer)?;
                                        // counter += 2;
                                        remaining_bytes -= 2;
                                        let mut short2: u16 = 0;
                                        short2 += (buffer[0] as u16) << 0;
                                        short2 += (buffer[1] as u16) << 8;
                                        dna_2_type_id.push(type_idx);
                                        let tuple_counter: usize = short2 as usize;
                                        for _t in 0..tuple_counter {
                                            // read two short values
                                            let mut buffer = [0; 2];
                                            f.read(&mut buffer)?;
                                            // counter += 2;
                                            remaining_bytes -= 2;
                                            f.read(&mut buffer)?;
                                            // counter += 2;
                                            remaining_bytes -= 2;
                                        }
                                    }
                                } else {
                                    println!("ERROR: \"STRC\" expected, \"S\"{:?} found", buffer)
                                }
                            } else {
                                println!("ERROR: \"TLEN\" expected, \"T\"{:?} found", buffer)
                            }
                        } else {
                            println!("ERROR: \"TYPE\" expected, \"T\"{:?} found", buffer)
                        }
                        // read remaining bytes
                        let mut buffer = vec![0; remaining_bytes];
                        f.read(&mut buffer)?;
                        // counter += remaining_bytes;
                    }
                } else {
                    // read len bytes
                    let mut buffer = vec![0; len as usize];
                    f.read(&mut buffer)?;
                    // counter += len as usize;
                }
            }
            // println!("{} bytes read", counter);
        }
    }
    // then use the DNA
    let mut builder: SceneDescriptionBuilder = SceneDescriptionBuilder::new();
    {
        let mut f = File::open(&args.path)?;
        // read exactly 12 bytes
        let mut counter: usize = 0;
        let mut buffer = [0; 12];
        f.read(&mut buffer)?;
        counter += 12;
        let mut blender_version: u32 = 0;
        if !decode_blender_header(&buffer, &mut blender_version, false) {
            println!("ERROR: Not a .blend file");
            println!("First 12 bytes:");
            println!("{:?}", buffer);
        } else {
            let mut data_following_mesh: bool = false;
            let mut is_smooth: bool = false;
            // Blender
            let mut loop_indices: Vec<u32> = Vec::new();
            loop {
                // code
                let mut buffer = [0; 4];
                f.read(&mut buffer)?;
                counter += 4;
                let code = make_id(&buffer);
                // len
                let mut buffer = [0; 4];
                f.read(&mut buffer)?;
                counter += 4;
                let mut len: u32 = 0;
                len += (buffer[0] as u32) << 0;
                len += (buffer[1] as u32) << 8;
                len += (buffer[2] as u32) << 16;
                len += (buffer[3] as u32) << 24;
                // for now ignore the old entry
                let mut buffer = [0; 8];
                f.read(&mut buffer)?;
                counter += 8;
                // get SDNAnr
                let mut buffer = [0; 4];
                f.read(&mut buffer)?;
                counter += 4;
                let mut sdna_nr: u32 = 0;
                sdna_nr += (buffer[0] as u32) << 0;
                sdna_nr += (buffer[1] as u32) << 8;
                sdna_nr += (buffer[2] as u32) << 16;
                sdna_nr += (buffer[3] as u32) << 24;
                // get data len
                let mut buffer = [0; 4];
                f.read(&mut buffer)?;
                counter += 4;
                let mut data_len: u32 = 0;
                data_len += (buffer[0] as u32) << 0;
                data_len += (buffer[1] as u32) << 8;
                data_len += (buffer[2] as u32) << 16;
                data_len += (buffer[3] as u32) << 24;
                // are we done?
                if code == String::from("ENDB") {
                    break;
                }
                if code == String::from("DNA1") {
                    // read len bytes
                    let mut buffer = vec![0; len as usize];
                    f.read(&mut buffer)?;
                    counter += len as usize;
                    if data_following_mesh {
                        if let Some(o2w) = object_to_world_hm.get(&base_name) {
                            object_to_world = *o2w;
                        } else {
                            println!(
                                "WARNING: looking up object_to_world by name ({:?}) failed",
                                base_name
                            );
                        }
                        let world_to_object: Transform = Transform::inverse(&object_to_world);
                        let n_triangles: usize = vertex_indices.len() / 3;
                        // transform mesh vertices to world space
                        let mut p_ws: Vec<Point3f> = Vec::new();
                        let n_vertices: usize = p.len();
                        for i in 0..n_vertices {
                            p_ws.push(object_to_world.transform_point(&p[i]));
                        }
                        let mut n_ws: Vec<Normal3f> = Vec::new();
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
                        let uv: Vec<Point2f> = Vec::new();
                        builder.add_mesh(
                            base_name.clone(),
                            object_to_world,
                            world_to_object,
                            n_triangles,
                            vertex_indices.clone(),
                            n_vertices,
                            p_ws, // in world space
                            s,    // empty
                            n_ws, // in world space
                            uv,   // empty
                        );
                    }
                    // reset booleans
                    data_following_mesh = false;
                    is_smooth = false;
                } else {
                    // read len bytes
                    let mut buffer = vec![0; len as usize];
                    f.read(&mut buffer)?;
                    counter += len as usize;
                    if code == String::from("OB") {
                        // OB
                        // println!("{} ({})", code, len);
                        // println!("  SDNAnr = {}", sdna_nr);
                        // v279: Object (len=1440) { ... }
                        // v280: Object (len=1408) { ... }
                        let mut skip_bytes: usize = 0;
                        // id
                        let mut id_name = String::new();
                        base_name = String::new();
                        for i in 32..(32 + 66) {
                            if buffer[i] == 0 {
                                break;
                            }
                            if (buffer[i] as char).is_ascii_alphanumeric() {
                                id_name.push(buffer[i] as char);
                                if i != 32 && i != 33 {
                                    base_name.push(buffer[i] as char);
                                }
                            }
                        }
                        // println!("  id_name = {}", id_name);
                        // println!("  base_name = {}", base_name);
                        if blender_version < 280 {
                            skip_bytes += 120;
                        } else {
                            skip_bytes += 152;
                        }
                        // adt
                        skip_bytes += 8;
                        if blender_version < 280 {
                            // nothing there
                        } else {
                            // DrawDataList (len=16)
                            skip_bytes += 16;
                        }
                        // sculpt
                        skip_bytes += 8;
                        // type
                        // let mut ob_type: u16 = 0;
                        // ob_type += (buffer[skip_bytes] as u16) << 0;
                        // ob_type += (buffer[skip_bytes + 1] as u16) << 8;
                        skip_bytes += 2;
                        // match ob_type {
                        //     0 => println!("  ob_type = {}", "OB_EMPTY"),
                        //     1 => println!("  ob_type = {}", "OB_MESH"),
                        //     11 => println!("  ob_type = {}", "OB_CAMERA"),
                        //     _ => println!("  ob_type = {}", ob_type),
                        // }
                        // partype
                        skip_bytes += 2;
                        // par1, par2, par3
                        skip_bytes += 4 * 3;
                        // parsubstr[64]
                        skip_bytes += 64;
                        // parent, track, proxy, proxy_group, proxy_from
                        skip_bytes += 8 * 5;
                        // ipo
                        skip_bytes += 8;
                        if blender_version < 280 {
                            // bb
                            skip_bytes += 8;
                        } else {
                            // nothing there
                        }
                        // action, poselib, pose, data, gpd
                        skip_bytes += 8 * 5;
                        // v279: bAnimVizSettings (len=48)
                        // v280: bAnimVizSettings (len=32)
                        if blender_version < 280 {
                            skip_bytes += 48;
                        } else {
                            skip_bytes += 32;
                        }
                        // mpath
                        skip_bytes += 8;
                        if blender_version < 280 {
                            // ListBase * 4
                            skip_bytes += 16 * 4;
                        } else {
                            // _pad0
                            skip_bytes += 8;
                            // ListBase * 7
                            skip_bytes += 16 * 7;
                        }
                        // mode, restore_mode
                        skip_bytes += 4 * 2;
                        // mat, matbits
                        skip_bytes += 8 * 2;
                        // totcol, actcol
                        skip_bytes += 4 * 2;
                        // loc
                        skip_bytes += 4 * 3;
                        // dloc
                        skip_bytes += 4 * 3;
                        if blender_version < 280 {
                            // orig
                            skip_bytes += 4 * 3;
                        } else {
                            // nothing there
                        }
                        // size
                        skip_bytes += 4 * 3;
                        // dsize
                        skip_bytes += 4 * 3;
                        // dscale
                        skip_bytes += 4 * 3;
                        // rot
                        for _i in 0..3 {
                            let mut rot_buf: [u8; 4] = [0_u8; 4];
                            for i in 0..4 as usize {
                                rot_buf[i] = buffer[skip_bytes + i];
                            }
                            let _rot: f32 = unsafe { mem::transmute(rot_buf) };
                            // println!("  rot[{}] = {}", i, rot);
                            skip_bytes += 4;
                        }
                        //skip_bytes += 4 * 3;
                        // drot
                        skip_bytes += 4 * 3;
                        // quat
                        skip_bytes += 4 * 4;
                        // dquat
                        skip_bytes += 4 * 4;
                        // rotAxis
                        skip_bytes += 4 * 3;
                        // drotAxis
                        skip_bytes += 4 * 3;
                        // rotAngle
                        let mut rot_angle_buf: [u8; 4] = [0_u8; 4];
                        for i in 0..4 as usize {
                            rot_angle_buf[i] = buffer[skip_bytes + i];
                        }
                        let _rot_angle: f32 = unsafe { mem::transmute(rot_angle_buf) };
                        // println!("  rot_angle = {}", rot_angle);
                        skip_bytes += 4;
                        // drotAngle
                        skip_bytes += 4;
                        // obmat
                        let mut mat_values: [f32; 16] = [0.0_f32; 16];
                        for i in 0..4 {
                            for j in 0..4 {
                                let mut obmat_buf: [u8; 4] = [0_u8; 4];
                                for i in 0..4 as usize {
                                    obmat_buf[i] = buffer[skip_bytes + i];
                                }
                                let obmat: f32 = unsafe { mem::transmute(obmat_buf) };
                                // println!("  obmat[{}][{}] = {}", i, j, obmat);
                                mat_values[i * 4 + j] = obmat;
                                skip_bytes += 4;
                            }
                        }
                        object_to_world = Transform::new(
                            mat_values[0],
                            mat_values[4],
                            mat_values[8],
                            mat_values[12] * scale_length,
                            mat_values[1],
                            mat_values[5],
                            mat_values[9],
                            mat_values[13] * scale_length,
                            mat_values[2],
                            mat_values[6],
                            mat_values[10],
                            mat_values[14] * scale_length,
                            mat_values[3],
                            mat_values[7],
                            mat_values[11],
                            mat_values[15],
                        );
                        object_to_world_hm.insert(base_name.clone(), object_to_world);
                        // parentinv
                        for _i in 0..4 {
                            for _j in 0..4 {
                                let mut parentinv_buf: [u8; 4] = [0_u8; 4];
                                for i in 0..4 as usize {
                                    parentinv_buf[i] = buffer[skip_bytes + i];
                                }
                                let _parentinv: f32 = unsafe { mem::transmute(parentinv_buf) };
                                // println!("  parentinv[{}][{}] = {}", i, j, parentinv);
                                skip_bytes += 4;
                            }
                        }
                        // constinv
                        for _i in 0..4 {
                            for _j in 0..4 {
                                let mut constinv_buf: [u8; 4] = [0_u8; 4];
                                for i in 0..4 as usize {
                                    constinv_buf[i] = buffer[skip_bytes + i];
                                }
                                let _constinv: f32 = unsafe { mem::transmute(constinv_buf) };
                                // println!("  constinv[{}][{}] = {}", i, j, constinv);
                                skip_bytes += 4;
                            }
                        }
                        // imat
                        for _i in 0..4 {
                            for _j in 0..4 {
                                let mut imat_buf: [u8; 4] = [0_u8; 4];
                                for i in 0..4 as usize {
                                    imat_buf[i] = buffer[skip_bytes + i];
                                }
                                let _imat: f32 = unsafe { mem::transmute(imat_buf) };
                                // println!("  imat[{}][{}] = {}", i, j, imat);
                                skip_bytes += 4;
                            }
                        }
                        // imat_ren
                        for _i in 0..4 {
                            for _j in 0..4 {
                                let mut imat_ren_buf: [u8; 4] = [0_u8; 4];
                                for i in 0..4 as usize {
                                    imat_ren_buf[i] = buffer[skip_bytes + i];
                                }
                                let _imat_ren: f32 = unsafe { mem::transmute(imat_ren_buf) };
                                // println!("  imat_ren[{}][{}] = {}", i, j, imat_ren);
                                skip_bytes += 4;
                            }
                        }
                        // reset booleans
                        data_following_mesh = false;
                        is_smooth = false;
                    } else if code == String::from("ME") {
                        if data_following_mesh {
                            if let Some(o2w) = object_to_world_hm.get(&base_name) {
                                object_to_world = *o2w;
                            } else {
                                println!(
                                    "WARNING: looking up object_to_world by name ({:?}) failed",
                                    base_name
                                );
                            }
                            let world_to_object: Transform = Transform::inverse(&object_to_world);
                            let n_triangles: usize = vertex_indices.len() / 3;
                            // transform mesh vertices to world space
                            let mut p_ws: Vec<Point3f> = Vec::new();
                            let n_vertices: usize = p.len();
                            for i in 0..n_vertices {
                                p_ws.push(object_to_world.transform_point(&p[i]));
                            }
                            let mut n_ws: Vec<Normal3f> = Vec::new();
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
                            let uv: Vec<Point2f> = Vec::new();
                            builder.add_mesh(
                                base_name.clone(),
                                object_to_world,
                                world_to_object,
                                n_triangles,
                                vertex_indices.clone(),
                                n_vertices,
                                p_ws, // in world space
                                s,    // empty
                                n_ws, // in world space
                                uv,   // empty
                            );
                        }
                        // ME
                        // println!("{} ({})", code, len);
                        // println!("  SDNAnr = {}", sdna_nr);
                        // Mesh (len=1416) { ... }
                        // let mut skip_bytes: usize = 0;
                        // id
                        let mut id_name = String::new();
                        base_name = String::new();
                        for i in 32..(32 + 66) {
                            if buffer[i] == 0 {
                                break;
                            }
                            if (buffer[i] as char).is_ascii_alphanumeric() {
                                id_name.push(buffer[i] as char);
                                if i != 32 && i != 33 {
                                    base_name.push(buffer[i] as char);
                                }
                            }
                        }
                        // println!("  id_name = {}", id_name);
                        // println!("  base_name = {}", base_name);
                        // if blender_version < 280 {
                        //     skip_bytes += 120;
                        // } else {
                        //     skip_bytes += 152;
                        // }
                        // adt
                        // skip_bytes += 8;
                        // bb, ipo, key, mat, mselect, mpoly, mtpoly, mloop, mloopuv, mloopcol
                        // skip_bytes += 8 * 10;
                        // mface, mtface, tface, mvert, medge, dvert, mcol, texcomesh, edit_btmesh
                        // skip_bytes += 8 * 9;
                        // CustomData * 5
                        // skip_bytes += 208 * 5;
                        // totvert
                        // let mut totvert: u32 = 0;
                        // totvert += (buffer[skip_bytes] as u32) << 0;
                        // totvert += (buffer[skip_bytes + 1] as u32) << 8;
                        // totvert += (buffer[skip_bytes + 2] as u32) << 16;
                        // totvert += (buffer[skip_bytes + 3] as u32) << 24;
                        // println!("  totvert = {}", totvert);
                        // skip_bytes += 4;
                        // totedge
                        // let mut totedge: u32 = 0;
                        // totedge += (buffer[skip_bytes] as u32) << 0;
                        // totedge += (buffer[skip_bytes + 1] as u32) << 8;
                        // totedge += (buffer[skip_bytes + 2] as u32) << 16;
                        // totedge += (buffer[skip_bytes + 3] as u32) << 24;
                        // println!("  totedge = {}", totedge);
                        // skip_bytes += 4;
                        // totface
                        // let mut totface: u32 = 0;
                        // totface += (buffer[skip_bytes] as u32) << 0;
                        // totface += (buffer[skip_bytes + 1] as u32) << 8;
                        // totface += (buffer[skip_bytes + 2] as u32) << 16;
                        // totface += (buffer[skip_bytes + 3] as u32) << 24;
                        // println!("  totface = {}", totface);
                        // skip_bytes += 4;
                        // totselect
                        // let mut totselect: u32 = 0;
                        // totselect += (buffer[skip_bytes] as u32) << 0;
                        // totselect += (buffer[skip_bytes + 1] as u32) << 8;
                        // totselect += (buffer[skip_bytes + 2] as u32) << 16;
                        // totselect += (buffer[skip_bytes + 3] as u32) << 24;
                        // println!("  totselect = {}", totselect);
                        // skip_bytes += 4;
                        // totpoly
                        // let mut totpoly: u32 = 0;
                        // totpoly += (buffer[skip_bytes] as u32) << 0;
                        // totpoly += (buffer[skip_bytes + 1] as u32) << 8;
                        // totpoly += (buffer[skip_bytes + 2] as u32) << 16;
                        // totpoly += (buffer[skip_bytes + 3] as u32) << 24;
                        // println!("  totpoly = {}", totpoly);
                        // // skip_bytes += 4;
                        // totloop
                        // let mut totloop: u32 = 0;
                        // totloop += (buffer[skip_bytes] as u32) << 0;
                        // totloop += (buffer[skip_bytes + 1] as u32) << 8;
                        // totloop += (buffer[skip_bytes + 2] as u32) << 16;
                        // totloop += (buffer[skip_bytes + 3] as u32) << 24;
                        // println!("  totloop = {}", totloop);
                        // skip_bytes += 4;
                        // check tlen
                        // for n in 0..types.len() {
                        //     if types[n] == "CustomData" {
                        //         println!("  {:?} = types[{}] needs {} bytes", types[n], n, tlen[n]);
                        //     }
                        // }
                        // data_following_mesh
                        data_following_mesh = true;
                        // clear all Vecs
                        p.clear();
                        n.clear();
                        vertex_indices.clear();
                        loop_indices.clear();
                    } else if code == String::from("SC") {
                        // SC
                        // println!("{} ({})", code, len);
                        // println!("  SDNAnr = {}", sdna_nr);
                        // v279: Scene (len=5472) { ... }
                        // v280: Scene (len=6392) { ... }
                        let mut skip_bytes: usize = 0;
                        // v279: ID (len=120)
                        // v280: ID (len=152)
                        let mut id_name = String::new();
                        base_name = String::new();
                        for i in 32..(32 + 66) {
                            if buffer[i] == 0 {
                                break;
                            }
                            if (buffer[i] as char).is_ascii_alphanumeric() {
                                id_name.push(buffer[i] as char);
                                if i != 32 && i != 33 {
                                    base_name.push(buffer[i] as char);
                                }
                            }
                        }
                        // println!("  id_name = {}", id_name);
                        // println!("  base_name = {}", base_name);
                        if blender_version < 280 {
                            skip_bytes += 120;
                        } else {
                            skip_bytes += 152;
                        }
                        // adt
                        skip_bytes += 8;
                        // camera, world, set
                        skip_bytes += 8 * 3;
                        // ListBase * 1
                        skip_bytes += 16 * 1;
                        // basact
                        skip_bytes += 8;
                        if blender_version < 280 {
                            // obedit
                            skip_bytes += 8;
                            // cursor
                            skip_bytes += 4 * 3;
                            // twcent
                            skip_bytes += 4 * 3;
                            // twmin
                            skip_bytes += 4 * 3;
                            // twmax
                            skip_bytes += 4 * 3;
                        } else {
                            // _pad1
                            skip_bytes += 8;
                            // View3DCursor (len=64)
                            skip_bytes += 64;
                        }
                        // lay, layact, lay_updated/_pad2[4]
                        skip_bytes += 4 * 3;
                        // flag
                        skip_bytes += 2;
                        // use_nodes
                        skip_bytes += 1;
                        // pad/_pad3
                        skip_bytes += 1;
                        // nodetree, ed, toolsettings, stats/_pad4
                        skip_bytes += 8 * 4;
                        // DisplaySafeAreas (len=32)
                        skip_bytes += 32;
                        // v279: RenderData (len=4432)
                        // v280: RenderData (len=4192)
                        let mut render_data_bytes: usize = 0;
                        // ImageFormatData (len=256)
                        skip_bytes += 256;
                        render_data_bytes += 256;
                        // avicodecdata
                        skip_bytes += 8;
                        render_data_bytes += 8;
                        if blender_version < 280 {
                            // qtcodecdata
                            skip_bytes += 8;
                            render_data_bytes += 8;
                            // QuicktimeCodecSettings (len=64)
                            skip_bytes += 64;
                            render_data_bytes += 64;
                        } else {
                            // nothing there
                        }
                        // FFMpegCodecData (len=88)
                        skip_bytes += 88;
                        render_data_bytes += 88;
                        // cfra
                        // let mut cfra: u32 = 0;
                        // cfra += (buffer[skip_bytes] as u32) << 0;
                        // cfra += (buffer[skip_bytes + 1] as u32) << 8;
                        // cfra += (buffer[skip_bytes + 2] as u32) << 16;
                        // cfra += (buffer[skip_bytes + 3] as u32) << 24;
                        // println!("    cfra = {}", cfra);
                        skip_bytes += 4;
                        render_data_bytes += 4;
                        // sfra
                        // let mut sfra: u32 = 0;
                        // sfra += (buffer[skip_bytes] as u32) << 0;
                        // sfra += (buffer[skip_bytes + 1] as u32) << 8;
                        // sfra += (buffer[skip_bytes + 2] as u32) << 16;
                        // sfra += (buffer[skip_bytes + 3] as u32) << 24;
                        // println!("    sfra = {}", sfra);
                        skip_bytes += 4;
                        render_data_bytes += 4;
                        // efra
                        // let mut efra: u32 = 0;
                        // efra += (buffer[skip_bytes] as u32) << 0;
                        // efra += (buffer[skip_bytes + 1] as u32) << 8;
                        // efra += (buffer[skip_bytes + 2] as u32) << 16;
                        // efra += (buffer[skip_bytes + 3] as u32) << 24;
                        // println!("    efra = {}", efra);
                        skip_bytes += 4;
                        render_data_bytes += 4;
                        // subframe
                        skip_bytes += 4;
                        render_data_bytes += 4;
                        // psfra, pefra, images, framapto
                        skip_bytes += 4 * 4;
                        render_data_bytes += 4 * 4;
                        // flag, threads
                        skip_bytes += 2 * 2;
                        render_data_bytes += 2 * 2;
                        // framelen, blurfac,
                        skip_bytes += 4 * 2;
                        render_data_bytes += 4 * 2;
                        if blender_version < 280 {
                            // edgeR, edgeG, edgeB
                            skip_bytes += 4 * 3;
                            render_data_bytes += 4 * 3;
                            // fullscreen, xplay, yplay, freqplay, depth, attrib
                            skip_bytes += 2 * 6;
                            render_data_bytes += 2 * 6;
                        } else {
                            // nothing there
                        }
                        // frame_step
                        skip_bytes += 4;
                        render_data_bytes += 4;
                        // stereomode, dimensionspreset
                        skip_bytes += 2 * 2;
                        render_data_bytes += 2 * 2;
                        if blender_version < 280 {
                            // filtertype
                            skip_bytes += 2;
                            render_data_bytes += 2;
                        } else {
                            // nothing there
                        }
                        // size
                        skip_bytes += 2;
                        render_data_bytes += 2;
                        if blender_version < 280 {
                            // maximsize, pad6
                            skip_bytes += 2;
                            render_data_bytes += 2;
                        } else {
                            // nothing there
                        }
                        // pad6
                        skip_bytes += 2;
                        render_data_bytes += 2;
                        // xsch
                        let mut xsch: u32 = 0;
                        xsch += (buffer[skip_bytes] as u32) << 0;
                        xsch += (buffer[skip_bytes + 1] as u32) << 8;
                        xsch += (buffer[skip_bytes + 2] as u32) << 16;
                        xsch += (buffer[skip_bytes + 3] as u32) << 24;
                        // println!("    xsch = {}", xsch);
                        skip_bytes += 4;
                        render_data_bytes += 4;
                        resolution_x = xsch;
                        // ysch
                        let mut ysch: u32 = 0;
                        ysch += (buffer[skip_bytes] as u32) << 0;
                        ysch += (buffer[skip_bytes + 1] as u32) << 8;
                        ysch += (buffer[skip_bytes + 2] as u32) << 16;
                        ysch += (buffer[skip_bytes + 3] as u32) << 24;
                        // println!("    ysch = {}", ysch);
                        skip_bytes += 4;
                        render_data_bytes += 4;
                        resolution_y = ysch;
                        // skip remaining RenderData
                        if blender_version < 280 {
                            skip_bytes += 4432 - render_data_bytes;
                        } else {
                            skip_bytes += 4192 - render_data_bytes;
                        }
                        // AudioData (len=32)
                        skip_bytes += 32;
                        // ListBase * 2
                        skip_bytes += 16 * 2;
                        if blender_version < 280 {
                            // nothing there
                        } else {
                            // TransformOrientationSlot (16) * 4
                            skip_bytes += 16 * 4;
                        }
                        // sound_scene, playback_handle, sound_scrub_handle, speaker_handles, fps_info
                        skip_bytes += 8 * 5;
                        // depsgraph/depsgraph_hash
                        skip_bytes += 8;
                        if blender_version < 280 {
                            // pad1, theDag
                            skip_bytes += 8 * 2;
                            // dagflags, pad3
                            skip_bytes += 2 * 2;
                        } else {
                            // _pad7
                            skip_bytes += 4;
                        }
                        // active_keyingset
                        skip_bytes += 4;
                        // ListBase * 1
                        skip_bytes += 16 * 1;
                        if blender_version < 280 {
                            // GameFraming (len=16)
                            skip_bytes += 16;
                            // GameData (len=192)
                            skip_bytes += 192;
                        } else {
                            // nothing there
                        }
                        // v279: UnitSettings (len=8)
                        // v280: UnitSettings (len=16)
                        // scale_length
                        let mut scale_length_buf: [u8; 4] = [0_u8; 4];
                        for i in 0..4 as usize {
                            scale_length_buf[i] = buffer[skip_bytes + i];
                        }
                        scale_length = unsafe { mem::transmute(scale_length_buf) };
                        // println!("    scale_length = {}", scale_length);
                        // skip_bytes += 4;
                        // reset booleans
                        data_following_mesh = false;
                        is_smooth = false;
                    } else if code == String::from("CA") {
                        // CA
                        // println!("{} ({})", code, len);
                        // println!("  SDNAnr = {}", sdna_nr);
                        // v279: Camera (len=248) { ... }
                        // v280: Camera (len=520) { ... }
                        let mut skip_bytes: usize = 0;
                        // v279: ID (len=120)
                        // v280: ID (len=152)
                        let mut id_name = String::new();
                        base_name = String::new();
                        for i in 32..(32 + 66) {
                            if buffer[i] == 0 {
                                break;
                            }
                            if (buffer[i] as char).is_ascii_alphanumeric() {
                                id_name.push(buffer[i] as char);
                                if i != 32 && i != 33 {
                                    base_name.push(buffer[i] as char);
                                }
                            }
                        }
                        // println!("  id_name = {}", id_name);
                        // println!("  base_name = {}", base_name);
                        if blender_version < 280 {
                            skip_bytes += 120;
                        } else {
                            skip_bytes += 152;
                        }
                        // adt
                        skip_bytes += 8;
                        // type, dtx
                        skip_bytes += 2;
                        // flag
                        skip_bytes += 2;
                        // passepartalpha
                        skip_bytes += 4;
                        // clipsta
                        // let mut clipsta_buf: [u8; 4] = [0_u8; 4];
                        // for i in 0..4 as usize {
                        //     clipsta_buf[i] = buffer[skip_bytes + i];
                        // }
                        // let clipsta: f32 = unsafe { mem::transmute(clipsta_buf) };
                        // println!("  clipsta = {}", clipsta);
                        skip_bytes += 4;
                        // clipend
                        // let mut clipend_buf: [u8; 4] = [0_u8; 4];
                        // for i in 0..4 as usize {
                        //     clipend_buf[i] = buffer[skip_bytes + i];
                        // }
                        // let clipend: f32 = unsafe { mem::transmute(clipend_buf) };
                        // println!("  clipend = {}", clipend);
                        skip_bytes += 4;
                        // lens
                        let mut lens_buf: [u8; 4] = [0_u8; 4];
                        for i in 0..4 as usize {
                            lens_buf[i] = buffer[skip_bytes + i];
                        }
                        let lens: f32 = unsafe { mem::transmute(lens_buf) };
                        // println!("  lens = {}", lens);
                        skip_bytes += 4;
                        // ortho_scale
                        // let mut ortho_scale_buf: [u8; 4] = [0_u8; 4];
                        // for i in 0..4 as usize {
                        //     ortho_scale_buf[i] = buffer[skip_bytes + i];
                        // }
                        // let ortho_scale: f32 = unsafe { mem::transmute(ortho_scale_buf) };
                        // println!("  ortho_scale = {}", ortho_scale);
                        skip_bytes += 4;
                        // drawsize
                        // let mut drawsize_buf: [u8; 4] = [0_u8; 4];
                        // for i in 0..4 as usize {
                        //     drawsize_buf[i] = buffer[skip_bytes + i];
                        // }
                        // let drawsize: f32 = unsafe { mem::transmute(drawsize_buf) };
                        // println!("  drawsize = {}", drawsize);
                        skip_bytes += 4;
                        // sensor_x
                        let mut sensor_x_buf: [u8; 4] = [0_u8; 4];
                        for i in 0..4 as usize {
                            sensor_x_buf[i] = buffer[skip_bytes + i];
                        }
                        let sensor_x: f32 = unsafe { mem::transmute(sensor_x_buf) };
                        // println!("  sensor_x = {}", sensor_x);
                        skip_bytes += 4;
                        // sensor_y
                        let mut sensor_y_buf: [u8; 4] = [0_u8; 4];
                        for i in 0..4 as usize {
                            sensor_y_buf[i] = buffer[skip_bytes + i];
                        }
                        let sensor_y: f32 = unsafe { mem::transmute(sensor_y_buf) };
                        // println!("  sensor_y = {}", sensor_y);
                        // skip_bytes += 4;
                        // calculate angle_x and angle_y
                        angle_x = degrees(focallength_to_fov(lens, sensor_x) as Float);
                        angle_y = degrees(focallength_to_fov(lens, sensor_y) as Float);
                        // println!("  angle_x = {}", angle_x);
                        // println!("  angle_y = {}", angle_y);
                        // shiftx
                        // let mut shiftx_buf: [u8; 4] = [0_u8; 4];
                        // for i in 0..4 as usize {
                        //     shiftx_buf[i] = buffer[skip_bytes + i];
                        // }
                        // let shiftx: f32 = unsafe { mem::transmute(shiftx_buf) };
                        // println!("  shiftx = {}", shiftx);
                        // skip_bytes += 4;
                        // shifty
                        // let mut shifty_buf: [u8; 4] = [0_u8; 4];
                        // for i in 0..4 as usize {
                        //     shifty_buf[i] = buffer[skip_bytes + i];
                        // }
                        // let shifty: f32 = unsafe { mem::transmute(shifty_buf) };
                        // println!("  shifty = {}", shifty);
                        // skip_bytes += 4;
                        // reset booleans
                        data_following_mesh = false;
                        is_smooth = false;
                    } else if code == String::from("MA") {
                        if data_following_mesh {
                            if let Some(o2w) = object_to_world_hm.get(&base_name) {
                                object_to_world = *o2w;
                            } else {
                                println!(
                                    "WARNING: looking up object_to_world by name ({:?}) failed",
                                    base_name
                                );
                            }
                            let world_to_object: Transform = Transform::inverse(&object_to_world);
                            let n_triangles: usize = vertex_indices.len() / 3;
                            // transform mesh vertices to world space
                            let mut p_ws: Vec<Point3f> = Vec::new();
                            let n_vertices: usize = p.len();
                            for i in 0..n_vertices {
                                p_ws.push(object_to_world.transform_point(&p[i]));
                            }
                            let mut n_ws: Vec<Normal3f> = Vec::new();
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
                            let uv: Vec<Point2f> = Vec::new();
                            builder.add_mesh(
                                base_name.clone(),
                                object_to_world,
                                world_to_object,
                                n_triangles,
                                vertex_indices.clone(),
                                n_vertices,
                                p_ws, // in world space
                                s,    // empty
                                n_ws, // in world space
                                uv,   // empty
                            );
                        }
                        // MA
                        // println!("{} ({})", code, len);
                        // println!("  SDNAnr = {}", sdna_nr);
                        // v279: Material (len=1528) { ... }
                        // v280: Material (len=320) { ... }
                        let mut skip_bytes: usize = 0;
                        // v279: ID (len=120)
                        // v280: ID (len=152)
                        let mut id_name = String::new();
                        base_name = String::new();
                        for i in 32..(32 + 66) {
                            if buffer[i] == 0 {
                                break;
                            }
                            if (buffer[i] as char).is_ascii_alphanumeric() {
                                id_name.push(buffer[i] as char);
                                if i != 32 && i != 33 {
                                    base_name.push(buffer[i] as char);
                                }
                            }
                        }
                        // println!("  id_name = {}", id_name);
                        // println!("  base_name = {}", base_name);
                        if blender_version < 280 {
                            skip_bytes += 120;
                        } else {
                            skip_bytes += 152;
                        }
                        // adt
                        skip_bytes += 8;
                        if blender_version < 280 {
                            // material_type, flag
                            skip_bytes += 2 * 2;
                            // r
                            let mut r_buf: [u8; 4] = [0_u8; 4];
                            for i in 0..4 as usize {
                                r_buf[i] = buffer[skip_bytes + i];
                            }
                            let r: f32 = unsafe { mem::transmute(r_buf) };
                            // println!("  r = {}", r);
                            skip_bytes += 4;
                            // g
                            let mut g_buf: [u8; 4] = [0_u8; 4];
                            for i in 0..4 as usize {
                                g_buf[i] = buffer[skip_bytes + i];
                            }
                            let g: f32 = unsafe { mem::transmute(g_buf) };
                            // println!("  g = {}", g);
                            skip_bytes += 4;
                            // b
                            let mut b_buf: [u8; 4] = [0_u8; 4];
                            for i in 0..4 as usize {
                                b_buf[i] = buffer[skip_bytes + i];
                            }
                            let b: f32 = unsafe { mem::transmute(b_buf) };
                            // println!("  b = {}", b);
                            skip_bytes += 4;
                            // specr, specg, specb
                            skip_bytes += 4 * 3;
                            // mirr, mirg, mirb
                            skip_bytes += 4 * 3;
                            // ambr, ambg, ambb
                            skip_bytes += 4 * 3;
                            // amb
                            skip_bytes += 4;
                            // emit
                            let mut emit_buf: [u8; 4] = [0_u8; 4];
                            for i in 0..4 as usize {
                                emit_buf[i] = buffer[skip_bytes + i];
                            }
                            let emit: f32 = unsafe { mem::transmute(emit_buf) };
                            // println!("  emit = {}", emit);
                            let mat: Blend279Material = Blend279Material {
                                r: r,
                                g: g,
                                b: b,
                                a: 1.0,
                                emit: emit,
                            };
                            // println!("  mat[{:?}] = {:?}", base_name, mat);
                            material_hm.insert(base_name.clone(), mat);
                        // skip_bytes += 4;
                        } else {
                            // flag
                            skip_bytes += 2;
                            // _pad1
                            skip_bytes += 2;
                            // r
                            let mut r_buf: [u8; 4] = [0_u8; 4];
                            for i in 0..4 as usize {
                                r_buf[i] = buffer[skip_bytes + i];
                            }
                            let r: f32 = unsafe { mem::transmute(r_buf) };
                            // println!("  r = {}", r);
                            skip_bytes += 4;
                            // g
                            let mut g_buf: [u8; 4] = [0_u8; 4];
                            for i in 0..4 as usize {
                                g_buf[i] = buffer[skip_bytes + i];
                            }
                            let g: f32 = unsafe { mem::transmute(g_buf) };
                            // println!("  g = {}", g);
                            skip_bytes += 4;
                            // b
                            let mut b_buf: [u8; 4] = [0_u8; 4];
                            for i in 0..4 as usize {
                                b_buf[i] = buffer[skip_bytes + i];
                            }
                            let b: f32 = unsafe { mem::transmute(b_buf) };
                            // println!("  b = {}", b);
                            skip_bytes += 4;
                            // a
                            let mut a_buf: [u8; 4] = [0_u8; 4];
                            for i in 0..4 as usize {
                                a_buf[i] = buffer[skip_bytes + i];
                            }
                            let a: f32 = unsafe { mem::transmute(a_buf) };
                            // println!("  a = {}", a);
                            skip_bytes += 4;
                            // specr, specg, specb, alpha
                            skip_bytes += 4 * 4;
                            // ray_mirror, spec, gloss_mir, roughness, metallic
                            skip_bytes += 4 * 5;
                            // use_nodes
                            let _use_nodes: u8 = buffer[skip_bytes] as u8;
                            // println!("  use_nodes = {}", use_nodes);
                            // skip_bytes += 1;
                            let mat: Blend279Material = Blend279Material {
                                r: r,
                                g: g,
                                b: b,
                                a: a,
                                emit: 0.0,
                            };
                            println!("  mat[{:?}] = {:?}", base_name, mat);
                            material_hm.insert(base_name.clone(), mat);
                        }
                        // reset booleans
                        data_following_mesh = false;
                        is_smooth = false;
                    } else if code == String::from("DATA") {
                        // DATA
                        if data_following_mesh {
                            // type_id
                            let type_id: usize = dna_2_type_id[sdna_nr as usize] as usize;
                            if types[type_id] == "MPoly" {
                                // println!("{}[{}] ({})", code, data_len, len);
                                // println!("  SDNAnr = {}", sdna_nr);
                                // println!("  {} ({})", types[type_id], tlen[type_id]);
                                let mut skip_bytes: usize = 0;
                                for _p in 0..data_len {
                                    // println!("  {}:", p + 1);
                                    // loopstart
                                    let mut loopstart: u32 = 0;
                                    loopstart += (buffer[skip_bytes] as u32) << 0;
                                    loopstart += (buffer[skip_bytes + 1] as u32) << 8;
                                    loopstart += (buffer[skip_bytes + 2] as u32) << 16;
                                    loopstart += (buffer[skip_bytes + 3] as u32) << 24;
                                    // println!("    loopstart = {}", loopstart);
                                    skip_bytes += 4;
                                    // totloop
                                    let mut totloop: u32 = 0;
                                    totloop += (buffer[skip_bytes] as u32) << 0;
                                    totloop += (buffer[skip_bytes + 1] as u32) << 8;
                                    totloop += (buffer[skip_bytes + 2] as u32) << 16;
                                    totloop += (buffer[skip_bytes + 3] as u32) << 24;
                                    // println!("    totloop = {}", totloop);
                                    skip_bytes += 4;
                                    // mat_nr
                                    // let mut mat_nr: u16 = 0;
                                    // mat_nr += (buffer[skip_bytes] as u16) << 0;
                                    // mat_nr += (buffer[skip_bytes + 1] as u16) << 8;
                                    // println!("    mat_nr = {}", mat_nr);
                                    skip_bytes += 2;
                                    // flag
                                    let flag: u8 = buffer[skip_bytes];
                                    // println!("    flag = {}", flag);
                                    is_smooth = flag % 2 == 1;
                                    // println!("    is_smooth = {}", is_smooth);
                                    skip_bytes += 1;
                                    // pad
                                    // println!("    pad = {}", buffer[skip_bytes]);
                                    skip_bytes += 1;
                                    // PBRT
                                    if totloop == 3_u32 {
                                        // triangle
                                        for i in 0..3 {
                                            vertex_indices.push(
                                                loop_indices[(loopstart + i) as usize] as usize,
                                            );
                                        }
                                    } else if totloop == 4_u32 {
                                        // quads
                                        vertex_indices
                                            .push(loop_indices[(loopstart + 0) as usize] as usize);
                                        vertex_indices
                                            .push(loop_indices[(loopstart + 1) as usize] as usize);
                                        vertex_indices
                                            .push(loop_indices[(loopstart + 2) as usize] as usize);
                                        vertex_indices
                                            .push(loop_indices[(loopstart + 0) as usize] as usize);
                                        vertex_indices
                                            .push(loop_indices[(loopstart + 2) as usize] as usize);
                                        vertex_indices
                                            .push(loop_indices[(loopstart + 3) as usize] as usize);
                                    } else {
                                        println!(
                                            "WARNING: quads or triangles expected (totloop = {})",
                                            totloop
                                        )
                                    }
                                }
                            // println!("  vertex_indices = {:?}", vertex_indices);
                            } else if types[type_id] == "MVert" {
                                // println!("{}[{}] ({})", code, data_len, len);
                                // println!("  SDNAnr = {}", sdna_nr);
                                // println!("  {} ({})", types[type_id], tlen[type_id]);
                                let mut skip_bytes: usize = 0;
                                let factor: f32 = 1.0 / 32767.0;
                                let mut coords: [f32; 3] = [0.0_f32; 3];
                                for _v in 0..data_len {
                                    // println!("  {}:", v + 1);
                                    // co
                                    for i in 0..3 {
                                        let mut co_buf: [u8; 4] = [0_u8; 4];
                                        for b in 0..4 as usize {
                                            co_buf[b] = buffer[skip_bytes + b];
                                        }
                                        let co: f32 = unsafe { mem::transmute(co_buf) };
                                        // println!("    co[{}] = {}", i, co);
                                        coords[i] = co;
                                        skip_bytes += 4;
                                    }
                                    p.push(Point3f {
                                        x: (coords[0] * scale_length) as Float,
                                        y: (coords[1] * scale_length) as Float,
                                        z: (coords[2] * scale_length) as Float,
                                    });
                                    // no
                                    for i in 0..3 {
                                        let mut no: i16 = 0;
                                        no += (buffer[skip_bytes] as i16) << 0;
                                        no += (buffer[skip_bytes + 1] as i16) << 8;
                                        let nof: f32 = no as f32 * factor;
                                        // println!("    no[{}] = {}", i, nof);
                                        coords[i] = nof;
                                        skip_bytes += 2;
                                    }
                                    n.push(Normal3f {
                                        x: coords[0] as Float,
                                        y: coords[1] as Float,
                                        z: coords[2] as Float,
                                    });
                                    // flag
                                    // println!("    flag = {}", buffer[skip_bytes]);
                                    skip_bytes += 1;
                                    // bweight
                                    // println!("    bweight = {}", buffer[skip_bytes]);
                                    skip_bytes += 1;
                                }
                            // for v in 0..data_len as usize {
                            //     println!("  {}:", v + 1);
                            //     println!("    co: {:?}", p[v]);
                            //     println!("    no: {:?}", n[v]);
                            // }
                            } else if types[type_id] == "MLoop" {
                                // println!("{}[{}] ({})", code, data_len, len);
                                // println!("  SDNAnr = {}", sdna_nr);
                                // println!("  {} ({})", types[type_id], tlen[type_id]);
                                let mut skip_bytes: usize = 0;
                                for _l in 0..data_len {
                                    // println!("  {}:", l + 1);
                                    // v
                                    let mut v: u32 = 0;
                                    v += (buffer[skip_bytes] as u32) << 0;
                                    v += (buffer[skip_bytes + 1] as u32) << 8;
                                    v += (buffer[skip_bytes + 2] as u32) << 16;
                                    v += (buffer[skip_bytes + 3] as u32) << 24;
                                    // println!("    v = {}", v);
                                    loop_indices.push(v);
                                    skip_bytes += 4;
                                    // e
                                    // let mut e: u32 = 0;
                                    // e += (buffer[skip_bytes] as u32) << 0;
                                    // e += (buffer[skip_bytes + 1] as u32) << 8;
                                    // e += (buffer[skip_bytes + 2] as u32) << 16;
                                    // e += (buffer[skip_bytes + 3] as u32) << 24;
                                    // println!("    e = {}", e);
                                    skip_bytes += 4;
                                }
                                // println!("    loop_indices = {:?}", loop_indices);
                            }
                        }
                    } else {
                        if data_following_mesh {
                            if let Some(o2w) = object_to_world_hm.get(&base_name) {
                                object_to_world = *o2w;
                            } else {
                                println!(
                                    "WARNING: looking up object_to_world by name ({:?}) failed",
                                    base_name
                                );
                            }
                            let world_to_object: Transform = Transform::inverse(&object_to_world);
                            let n_triangles: usize = vertex_indices.len() / 3;
                            // transform mesh vertices to world space
                            let mut p_ws: Vec<Point3f> = Vec::new();
                            let n_vertices: usize = p.len();
                            for i in 0..n_vertices {
                                p_ws.push(object_to_world.transform_point(&p[i]));
                            }
                            let mut n_ws: Vec<Normal3f> = Vec::new();
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
                            let uv: Vec<Point2f> = Vec::new();
                            builder.add_mesh(
                                base_name.clone(),
                                object_to_world,
                                world_to_object,
                                n_triangles,
                                vertex_indices.clone(),
                                n_vertices,
                                p_ws, // in world space
                                s,    // empty
                                n_ws, // in world space
                                uv,   // empty
                            );
                        }
                        // reset booleans
                        data_following_mesh = false;
                        is_smooth = false;
                    }
                    if code != String::from("DATA")
                        && code != String::from("REND")
                        && code != String::from("TEST")
                    {
                        let type_id: usize = dna_2_type_id[sdna_nr as usize] as usize;
                        if len != tlen[type_id] as u32 {
                            println!("WARNING: {} ({} != {})", code, len, tlen[type_id]);
                        }
                    }
                }
            }
            println!("{} bytes read", counter);
        }
    }
    let scene_description: SceneDescription = builder.finalize();
    let mut render_options: RenderOptions = RenderOptions::new(scene_description, &material_hm);
    assert!(render_options.triangles.len() == render_options.triangle_lights.len());
    for triangle_idx in 0..render_options.triangles.len() {
        let triangle = &render_options.triangles[triangle_idx];
        let triangle_material = &render_options.triangle_materials[triangle_idx];
        let triangle_light = &render_options.triangle_lights[triangle_idx];
        let geo_prim = Arc::new(GeometricPrimitive::new(
            triangle.clone(),
            Some(triangle_material.clone()),
            triangle_light.clone(),
            None,
        ));
        render_options.primitives.push(geo_prim.clone());
    }
    println!("number of lights = {:?}", render_options.lights.len());
    println!(
        "number of primitives = {:?}",
        render_options.primitives.len()
    );
    let accelerator = Arc::new(BVHAccel::new(
        render_options.primitives,
        4,
        SplitMethod::SAH,
    ));
    let scene: Scene = Scene::new(accelerator.clone(), render_options.lights);
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
    let fov: Float;
    if aspect > 1.0 {
        fov = angle_y;
    } else {
        fov = angle_x;
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
    let shutteropen: Float = 0.0;
    let shutterclose: Float = 1.0;
    let lensradius: Float = 0.0;
    let focaldistance: Float = 1e6;
    let crop: Bounds2f = Bounds2f {
        p_min: Point2f { x: 0.0, y: 0.0 },
        p_max: Point2f { x: 1.0, y: 1.0 },
    };
    let xw: Float = 0.5;
    let yw: Float = 0.5;
    let filter: Arc<Filter + Sync + Send> = Arc::new(BoxFilter {
        radius: Vector2f { x: xw, y: yw },
        inv_radius: Vector2f {
            x: 1.0 / xw,
            y: 1.0 / yw,
        },
    });
    let filename: String = String::from("spheres-differentials-texfilt.exr");
    let film: Arc<Film> = Arc::new(Film::new(
        Point2i {
            x: resolution_x as i32,
            y: resolution_y as i32,
        },
        crop,
        filter,
        35.0,
        filename,
        1.0,
        std::f32::INFINITY,
    ));
    let camera: Arc<Camera + Send + Sync> = Arc::new(PerspectiveCamera::new(
        animated_cam_to_world,
        screen,
        shutteropen,
        shutterclose,
        lensradius,
        focaldistance,
        fov,
        film.clone(),
        None,
    ));
    let mut sampler: Box<Sampler + Sync + Send> =
        Box::new(ZeroTwoSequenceSampler::new(args.samples as i64, 4));
    let sample_bounds: Bounds2i = film.get_sample_bounds();
    let mut integrator: Box<SamplerIntegrator + Send + Sync>;
    if let Some(integrator_str) = args.integrator {
        print!("integrator = {:?} [", integrator_str);
        if integrator_str == String::from("ao") {
            println!("Ambient Occlusion (AO)]");
            // AOIntegrator
            integrator = Box::new(AOIntegrator::new(true, 64_i32, sample_bounds));
            // in the end we want to call render()
            render(
                &scene,
                &camera.clone(),
                &mut sampler,
                &mut integrator,
                num_threads,
            );
        } else if integrator_str == String::from("directlighting") {
            println!("Direct Lighting]");
            // DirectLightingIntegrator
            let max_depth: i32 = 5;
            let strategy: LightStrategy = LightStrategy::UniformSampleAll;
            let pixel_bounds: Bounds2i = camera.get_film().get_sample_bounds();
            integrator = Box::new(DirectLightingIntegrator::new(
                strategy,
                max_depth as i64,
                pixel_bounds,
            ));
            // in the end we want to call render()
            render(
                &scene,
                &camera.clone(),
                &mut sampler,
                &mut integrator,
                num_threads,
            );
        } else if integrator_str == String::from("path") {
            println!("(Unidirectional) Path Tracing]");
            // PathIntegrator
            let max_depth: i32 = 5;
            let pixel_bounds: Bounds2i = camera.get_film().get_sample_bounds();
            let rr_threshold: Float = 1.0;
            let light_strategy: String = String::from("spatial");
            integrator = Box::new(PathIntegrator::new(
                max_depth as u32,
                pixel_bounds,
                rr_threshold,
                light_strategy,
            ));
            // in the end we want to call render()
            render(
                &scene,
                &camera.clone(),
                &mut sampler,
                &mut integrator,
                num_threads,
            );
        } else if integrator_str == String::from("bdpt") {
            println!("Bidirectional Path Tracing (BDPT)]");
            // BDPTIntegrator
            let max_depth: i32 = 5;
            let pixel_bounds: Bounds2i = camera.get_film().get_sample_bounds();
            let light_strategy: String = String::from("power");
            let mut integrator: Box<BDPTIntegrator> = Box::new(BDPTIntegrator::new(
                max_depth as u32,
                pixel_bounds,
                light_strategy,
            ));
            // in the end we want to call render()
            render_bdpt(
                &scene,
                &camera.clone(),
                &mut sampler,
                &mut integrator,
                num_threads,
            );
        } else if integrator_str == String::from("mlt") {
            println!("Metropolis Light Transport (MLT)]");
            // CreateMLTIntegrator
            let max_depth: i32 = 5;
            let mut n_bootstrap: i32 = 100000;
            let mut n_chains: i32 = 1000;
            let mut mutations_per_pixel: i32 = 100;
            let sigma: Float = 0.01;
            let large_step_probability: Float = 0.3;
            let mut integrator: Box<MLTIntegrator> = Box::new(MLTIntegrator::new(
                camera.clone(),
                max_depth as u32,
                n_bootstrap as u32,
                n_chains as u32,
                mutations_per_pixel as u32,
                sigma,
                large_step_probability,
            ));
            // in the end we want to call render()
            render_mlt(
                &scene,
                &camera.clone(),
                &mut sampler,
                &mut integrator,
                num_threads,
            );
        } else {
            println!("unknown]");
        }
    } else {
        if render_options.has_emitters {
            // PathIntegrator
            let pixel_bounds: Bounds2i = camera.get_film().get_sample_bounds();
            let rr_threshold: Float = 1.0;
            let light_strategy: String = String::from("spatial");
            let max_depth: i32 = 5;
            integrator = Box::new(PathIntegrator::new(
                max_depth as u32,
                pixel_bounds,
                rr_threshold,
                light_strategy,
            ));
        } else {
            // AOIntegrator
            integrator = Box::new(AOIntegrator::new(true, 64_i32, sample_bounds));
        }
        // in the end we want to call render()
        render(
            &scene,
            &camera.clone(),
            &mut sampler,
            &mut integrator,
            num_threads,
        );
    }
    Ok(())
}
