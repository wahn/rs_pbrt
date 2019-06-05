// Assumptions:
// 1. objects and data have the same base name:
//    e.g. OBcornellbox and MEcornellbox
//    TODO: Find a better solution (following the Object->data pointer?)
// 2. bpy.context.scene.unit_settings.scale_length = 0.001
//    TODO: get the scale_length from Blender file
// 3. Right now we search for "Camera" for Transform::look_at(...)
//    TODO: get the render camera/transform from the scene (self.scene.camera)

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
use pbrt::core::light::Light;
use pbrt::core::pbrt::{Float, Spectrum};
use pbrt::core::primitive::{GeometricPrimitive, Primitive};
use pbrt::core::sampler::Sampler;
use pbrt::core::scene::Scene;
use pbrt::core::transform::{AnimatedTransform, Transform};
use pbrt::filters::boxfilter::BoxFilter;
use pbrt::integrators::ao::AOIntegrator;
use pbrt::integrators::render;
use pbrt::materials::matte::MatteMaterial;
use pbrt::samplers::zerotwosequence::ZeroTwoSequenceSampler;
use pbrt::shapes::triangle::{Triangle, TriangleMesh};
use pbrt::textures::constant::ConstantTexture;

pub const VERSION: &'static str = env!("CARGO_PKG_VERSION");

/// Parse a Blender scene file and render it.
#[derive(StructOpt)]
struct Cli {
    /// The path to the file to read
    #[structopt(parse(from_os_str))]
    path: std::path::PathBuf,
}

// TMP (see pbrt_spheres_differentials_texfilt.rs)

struct SceneDescription {
    meshes: Vec<Arc<TriangleMesh>>,
    lights: Vec<Arc<Light + Sync + Send>>,
}

struct SceneDescriptionBuilder {
    meshes: Vec<Arc<TriangleMesh>>,
    lights: Vec<Arc<Light + Sync + Send>>,
}

impl SceneDescriptionBuilder {
    fn new() -> SceneDescriptionBuilder {
        SceneDescriptionBuilder {
            meshes: Vec::new(),
            lights: Vec::new(),
        }
    }
    fn add_mesh(
        &mut self,
        object_to_world: Transform,
        world_to_object: Transform,
        n_triangles: usize,
        vertex_indices: Vec<usize>,
        n_vertices: usize,
        p_ws: Vec<Point3f>,
        s: Vec<Vector3f>,
        n: Vec<Normal3f>,
        uv: Vec<Point2f>,
    ) -> &mut SceneDescriptionBuilder {
        let triangle_mesh = Arc::new(TriangleMesh::new(
            object_to_world,
            world_to_object,
            false,
            n_triangles,
            vertex_indices,
            n_vertices,
            p_ws, // in world space
            s,    // empty
            n,    // empty
            uv,   // empty
        ));
        self.meshes.push(triangle_mesh);
        self
    }
    fn finalize(self) -> SceneDescription {
        SceneDescription {
            meshes: self.meshes,
            lights: self.lights,
        }
    }
}

struct RenderOptions {
    primitives: Vec<Arc<Primitive + Sync + Send>>,
    triangles: Vec<Arc<Triangle>>,
    lights: Vec<Arc<Light + Sync + Send>>,
}

impl RenderOptions {
    fn new(scene: SceneDescription) -> RenderOptions {
        let primitives: Vec<Arc<Primitive + Sync + Send>> = Vec::new();
        let mut triangles: Vec<Arc<Triangle>> = Vec::new();
        let mut lights: Vec<Arc<Light + Sync + Send>> = Vec::new();
        // lights
        for light in &scene.lights {
            lights.push(light.clone());
        }
        // meshes
        for mesh in scene.meshes {
            // create individual triangles
            for id in 0..mesh.n_triangles {
                let triangle = Arc::new(Triangle::new(
                    mesh.object_to_world,
                    mesh.world_to_object,
                    mesh.transform_swaps_handedness,
                    mesh.clone(),
                    id,
                ));
                triangles.push(triangle);
            }
        }
        RenderOptions {
            primitives: primitives,
            triangles: triangles,
            lights: lights,
        }
    }
}

// TMP

fn decode_blender_header(header: &[u8], version: &mut u32) -> bool {
    // BLENDER
    match header[0] as char {
        'B' => print!("B"),
        _ => return false,
    }
    match header[1] as char {
        'L' => print!("L"),
        _ => return false,
    }
    match header[2] as char {
        'E' => print!("E"),
        _ => return false,
    }
    match header[3] as char {
        'N' => print!("N"),
        _ => return false,
    }
    match header[4] as char {
        'D' => print!("D"),
        _ => return false,
    }
    match header[5] as char {
        'E' => print!("E"),
        _ => return false,
    }
    match header[6] as char {
        'R' => print!("R"),
        _ => return false,
    }
    // [_|-]
    match header[7] as char {
        '_' => print!("_"),
        '-' => print!("-"),
        _ => return false,
    }
    // [v|V]
    match header[8] as char {
        'v' => print!("v"),
        'V' => print!("V"),
        _ => return false,
    }
    for i in 9..12 {
        if header[i].is_ascii_digit() {
            print!("{:?}", (header[i] as char).to_digit(10).unwrap());
        } else {
            return false;
        }
    }
    print!("\n");
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
    let mut name_counter: usize = 0;
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
            name_counter += 1;
        }
    }
    println!("  {} names found in {} bytes", name_counter, byte_counter);
    Ok(())
}

fn read_type_names(
    f: &mut File,
    nr_types: usize,
    type_names: &mut Vec<String>,
    byte_counter: &mut usize,
) -> std::io::Result<()> {
    let mut name_counter: usize = 0;
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
            name_counter += 1;
        }
    }
    println!(
        "  {} type names found in {} bytes",
        name_counter, byte_counter
    );
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
    let scale_length: f32 = 1.0; // 0.001;
    let mut base_name = String::new();
    let mut object_to_world_hm: HashMap<String, Transform> = HashMap::new();
    let mut object_to_world: Transform = Transform::new(
        1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    );
    let mut p: Vec<Point3f> = Vec::new();
    let mut n: Vec<Normal3f> = Vec::new();
    let mut vertex_indices: Vec<usize> = Vec::new();
    // TODO: let mut primitives: Vec<Arc<Primitive + Sync + Send>> = Vec::new();
    let mut lights: Vec<Arc<Light + Sync + Send>> = Vec::new();
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
        let mut counter: usize = 0;
        let mut buffer = [0; 12];
        f.read(&mut buffer)?;
        counter += 12;
        let mut blender_version: u32 = 0;
        if !decode_blender_header(&buffer, &mut blender_version) {
            println!("ERROR: Not a .blend file");
            println!("First 12 bytes:");
            println!("{:?}", buffer);
        } else {
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
                // for now ignore the nr entry
                let mut buffer = [0; 4];
                f.read(&mut buffer)?;
                counter += 4;
                // are we done?
                if code == String::from("ENDB") {
                    break;
                }
                if code == String::from("DNA1") {
                    println!("{} ({})", code, len);
                    // "SDNANAME" in first 8 bytes
                    let mut buffer = [0; 8];
                    f.read(&mut buffer)?;
                    counter += 8;
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
                        counter += (len - 8) as usize;
                    } else {
                        let mut buffer = [0; 4];
                        f.read(&mut buffer)?;
                        counter += 4;
                        let mut nr_names: u32 = 0;
                        nr_names += (buffer[0] as u32) << 0;
                        nr_names += (buffer[1] as u32) << 8;
                        nr_names += (buffer[2] as u32) << 16;
                        nr_names += (buffer[3] as u32) << 24;
                        read_names(&mut f, nr_names as usize, &mut names, &mut names_len)?;
                        counter += names_len;
                        let mut remaining_bytes: usize = (len - 12) as usize - names_len;
                        // skip pad bytes, read "TYPE" and nr_types
                        let mut buffer = [0; 1];
                        loop {
                            f.read(&mut buffer)?;
                            counter += 1;
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
                        counter += 3;
                        remaining_bytes -= 3;
                        if buffer[0] as char == 'Y'
                            && buffer[1] as char == 'P'
                            && buffer[2] as char == 'E'
                        {
                            // nr_types
                            let mut buffer = [0; 4];
                            f.read(&mut buffer)?;
                            counter += 4;
                            remaining_bytes -= 4;
                            let mut nr_types: u32 = 0;
                            nr_types += (buffer[0] as u32) << 0;
                            nr_types += (buffer[1] as u32) << 8;
                            nr_types += (buffer[2] as u32) << 16;
                            nr_types += (buffer[3] as u32) << 24;
                            read_type_names(&mut f, nr_types as usize, &mut types, &mut types_len)?;
                            counter += types_len;
                            remaining_bytes -= types_len;
                            // skip pad bytes, read "TLEN"
                            let mut buffer = [0; 1];
                            loop {
                                f.read(&mut buffer)?;
                                counter += 1;
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
                            counter += 3;
                            remaining_bytes -= 3;
                            if buffer[0] as char == 'L'
                                && buffer[1] as char == 'E'
                                && buffer[2] as char == 'N'
                            {
                                // read short (16 bits = 2 bytes) for each type
                                for _i in 0..nr_types as usize {
                                    let mut buffer = [0; 2];
                                    f.read(&mut buffer)?;
                                    counter += 2;
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
                                    counter += 1;
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
                                counter += 3;
                                remaining_bytes -= 3;
                                if buffer[0] as char == 'T'
                                    && buffer[1] as char == 'R'
                                    && buffer[2] as char == 'C'
                                {
                                    // nr_structs
                                    let mut buffer = [0; 4];
                                    f.read(&mut buffer)?;
                                    counter += 4;
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
                                        counter += 2;
                                        remaining_bytes -= 2;
                                        let mut type_idx: u16 = 0;
                                        type_idx += (buffer[0] as u16) << 0;
                                        type_idx += (buffer[1] as u16) << 8;
                                        f.read(&mut buffer)?;
                                        counter += 2;
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
                                            counter += 2;
                                            remaining_bytes -= 2;
                                            f.read(&mut buffer)?;
                                            counter += 2;
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
                        counter += remaining_bytes;
                    }
                } else {
                    // read len bytes
                    let mut buffer = vec![0; len as usize];
                    f.read(&mut buffer)?;
                    counter += len as usize;
                }
            }
            println!("{} bytes read", counter);
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
        if !decode_blender_header(&buffer, &mut blender_version) {
            println!("ERROR: Not a .blend file");
            println!("First 12 bytes:");
            println!("{:?}", buffer);
        } else {
            let mut data_following_mesh: bool = false;
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
                        let s: Vec<Vector3f> = Vec::new();
                        let n: Vec<Normal3f> = Vec::new();
                        let uv: Vec<Point2f> = Vec::new();
                        builder.add_mesh(
                            object_to_world,
                            world_to_object,
                            n_triangles,
                            vertex_indices.clone(),
                            n_vertices,
                            p_ws, // in world space
                            s,    // empty
                            n,    // empty
                            uv,   // empty
                        );
                    }
                    data_following_mesh = false;
                } else {
                    // read len bytes
                    let mut buffer = vec![0; len as usize];
                    f.read(&mut buffer)?;
                    counter += len as usize;
                    if code == String::from("OB") {
                        // OB
                        println!("{} ({})", code, len);
                        println!("  SDNAnr = {}", sdna_nr);
                        // Object (len=1440) { ... }
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
                        println!("  id_name = {}", id_name);
                        println!("  base_name = {}", base_name);
                        skip_bytes += 120;
                        // adt
                        skip_bytes += 8;
                        // sculpt
                        skip_bytes += 8;
                        // type
                        let mut ob_type: u16 = 0;
                        ob_type += (buffer[skip_bytes] as u16) << 0;
                        ob_type += (buffer[skip_bytes + 1] as u16) << 8;
                        skip_bytes += 2;
                        match ob_type {
                            0 => println!("  ob_type = {}", "OB_EMPTY"),
                            1 => println!("  ob_type = {}", "OB_MESH"),
                            11 => println!("  ob_type = {}", "OB_CAMERA"),
                            _ => println!("  ob_type = {}", ob_type),
                        }
                        // partype
                        skip_bytes += 2;
                        // par1, par2, par3
                        skip_bytes += 4 * 3;
                        // parsubstr[64]
                        skip_bytes += 64;
                        // parent, track, proxy, proxy_group, proxy_from
                        skip_bytes += 8 * 5;
                        // ipo, bb, action, poselib, pose, data, gpd
                        skip_bytes += 8 * 7;
                        // bAnimVizSettings
                        skip_bytes += 48;
                        // mpath
                        skip_bytes += 8;
                        // ListBase * 4
                        skip_bytes += 16 * 4;
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
                        // orig
                        skip_bytes += 4 * 3;
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
                        data_following_mesh = false;
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
                            let s: Vec<Vector3f> = Vec::new();
                            let n: Vec<Normal3f> = Vec::new();
                            let uv: Vec<Point2f> = Vec::new();
                            builder.add_mesh(
                                object_to_world,
                                world_to_object,
                                n_triangles,
                                vertex_indices.clone(),
                                n_vertices,
                                p_ws, // in world space
                                s,    // empty
                                n,    // empty
                                uv,   // empty
                            );
                        }
                        // ME
                        println!("{} ({})", code, len);
                        println!("  SDNAnr = {}", sdna_nr);
                        // Mesh (len=1416) { ... }
                        let mut skip_bytes: usize = 0;
                        // id
                        let mut id_name = String::with_capacity(4);
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
                        println!("  id_name = {}", id_name);
                        println!("  base_name = {}", base_name);
                        skip_bytes += 120;
                        // adt
                        skip_bytes += 8;
                        // bb, ipo, key, mat, mselect, mpoly, mtpoly, mloop, mloopuv, mloopcol
                        skip_bytes += 8 * 10;
                        // mface, mtface, tface, mvert, medge, dvert, mcol, texcomesh, edit_btmesh
                        skip_bytes += 8 * 9;
                        // CustomData * 5
                        skip_bytes += 208 * 5;
                        // totvert
                        let mut totvert: u32 = 0;
                        totvert += (buffer[skip_bytes] as u32) << 0;
                        totvert += (buffer[skip_bytes + 1] as u32) << 8;
                        totvert += (buffer[skip_bytes + 2] as u32) << 16;
                        totvert += (buffer[skip_bytes + 3] as u32) << 24;
                        println!("  totvert = {}", totvert);
                        skip_bytes += 4;
                        // totedge
                        // let mut totedge: u32 = 0;
                        // totedge += (buffer[skip_bytes] as u32) << 0;
                        // totedge += (buffer[skip_bytes + 1] as u32) << 8;
                        // totedge += (buffer[skip_bytes + 2] as u32) << 16;
                        // totedge += (buffer[skip_bytes + 3] as u32) << 24;
                        // println!("  totedge = {}", totedge);
                        skip_bytes += 4;
                        // totface
                        // let mut totface: u32 = 0;
                        // totface += (buffer[skip_bytes] as u32) << 0;
                        // totface += (buffer[skip_bytes + 1] as u32) << 8;
                        // totface += (buffer[skip_bytes + 2] as u32) << 16;
                        // totface += (buffer[skip_bytes + 3] as u32) << 24;
                        // println!("  totface = {}", totface);
                        skip_bytes += 4;
                        // totselect
                        // let mut totselect: u32 = 0;
                        // totselect += (buffer[skip_bytes] as u32) << 0;
                        // totselect += (buffer[skip_bytes + 1] as u32) << 8;
                        // totselect += (buffer[skip_bytes + 2] as u32) << 16;
                        // totselect += (buffer[skip_bytes + 3] as u32) << 24;
                        // println!("  totselect = {}", totselect);
                        skip_bytes += 4;
                        // totpoly
                        let mut totpoly: u32 = 0;
                        totpoly += (buffer[skip_bytes] as u32) << 0;
                        totpoly += (buffer[skip_bytes + 1] as u32) << 8;
                        totpoly += (buffer[skip_bytes + 2] as u32) << 16;
                        totpoly += (buffer[skip_bytes + 3] as u32) << 24;
                        println!("  totpoly = {}", totpoly);
                        // skip_bytes += 4;
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
                        data_following_mesh = true;
                        // clear all Vecs
                        p.clear();
                        n.clear();
                        vertex_indices.clear();
                        loop_indices.clear();
                    } else if code == String::from("DATA") {
                        // DATA
                        if data_following_mesh {
                            // type_id
                            let type_id: usize = dna_2_type_id[sdna_nr as usize] as usize;
                            if types[type_id] == "MPoly" {
                                println!("{}[{}] ({})", code, data_len, len);
                                println!("  SDNAnr = {}", sdna_nr);
                                println!("  {} ({})", types[type_id], tlen[type_id]);
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
                                    // println!("    flag = {}", buffer[skip_bytes]);
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
                                println!("{}[{}] ({})", code, data_len, len);
                                println!("  SDNAnr = {}", sdna_nr);
                                println!("  {} ({})", types[type_id], tlen[type_id]);
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
                                println!("{}[{}] ({})", code, data_len, len);
                                println!("  SDNAnr = {}", sdna_nr);
                                println!("  {} ({})", types[type_id], tlen[type_id]);
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
                            let s: Vec<Vector3f> = Vec::new();
                            let n: Vec<Normal3f> = Vec::new();
                            let uv: Vec<Point2f> = Vec::new();
                            builder.add_mesh(
                                object_to_world,
                                world_to_object,
                                n_triangles,
                                vertex_indices.clone(),
                                n_vertices,
                                p_ws, // in world space
                                s,    // empty
                                n,    // empty
                                uv,   // empty
                            );
                        }
                        data_following_mesh = false;
                    }
                    if code != String::from("DATA")
                        && code != String::from("REND")
                        && code != String::from("TEST")
                    {
                        let type_id: usize = dna_2_type_id[sdna_nr as usize] as usize;
                        if len != tlen[type_id] as u32 {
                            println!("{} ({} != {})", code, len, tlen[type_id]);
                        }
                    }
                }
            }
            println!("{} bytes read", counter);
        }
    }
    let scene_description: SceneDescription = builder.finalize();
    let mut render_options: RenderOptions = RenderOptions::new(scene_description);
    let kd = Arc::new(ConstantTexture::new(Spectrum::new(0.5)));
    let sigma = Arc::new(ConstantTexture::new(0.0 as Float));
    let matte = Arc::new(MatteMaterial::new(kd, sigma, None));
    for triangle in render_options.triangles {
        let geo_prim = Arc::new(GeometricPrimitive::new(
            triangle,
            Some(matte.clone()),
            None,
            None,
        ));
        render_options.primitives.push(geo_prim.clone());
    }
    // WORK
    println!("number of lights = {:?}", lights.len());
    println!(
        "number of primitives = {:?}",
        render_options.primitives.len()
    );
    // TMP
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
    let fov: Float = 39.14625166082039;
    let xres = 500;
    let yres = 500;
    let frame: Float = xres as Float / yres as Float;
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
        Point2i { x: xres, y: yres },
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
    let mut sampler: Box<Sampler + Sync + Send> = Box::new(ZeroTwoSequenceSampler::default());
    let sample_bounds: Bounds2i = film.get_sample_bounds();
    let mut integrator: Box<SamplerIntegrator + Send + Sync> =
        Box::new(AOIntegrator::new(true, 64_i32, sample_bounds));
    // TMP
    // in the end we want to call render()
    render(
        &scene,
        &camera.clone(),
        &mut sampler,
        &mut integrator,
        num_threads,
    );
    Ok(())
}
