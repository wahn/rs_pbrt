extern crate num_cpus;
extern crate pbrt;
extern crate structopt;

// std
use std::fs::File;
use std::io::Read;
use std::sync::Arc;
use structopt::StructOpt;
// pbrt
use pbrt::core::light::Light;
use pbrt::core::primitive::Primitive;

pub const VERSION: &'static str = env!("CARGO_PKG_VERSION");

/// Parse a Blender scene file and render it.
#[derive(StructOpt)]
struct Cli {
    /// The path to the file to read
    #[structopt(parse(from_os_str))]
    path: std::path::PathBuf,
}

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

fn main() -> std::io::Result<()> {
    let args = Cli::from_args();
    let num_threads: u8 = num_cpus::get() as u8;
    println!(
        "parse_blend_file version {} [Detected {} cores]",
        VERSION, num_threads
    );
    let mut primitives: Vec<Arc<Primitive + Sync + Send>> = Vec::new();
    let mut lights: Vec<Arc<Light + Sync + Send>> = Vec::new();
    // first get the DNA
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
                // println!("{} ({})", code, len);
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
                // println!("  SDNAnr = {}", sdna_nr);
                // for now ignore the nr entry
                let mut buffer = [0; 4];
                f.read(&mut buffer)?;
                counter += 4;
                // are we done?
                if code == String::from("ENDB") {
                    break;
                }
                // read len bytes
                let mut buffer = vec![0; len as usize];
                f.read(&mut buffer)?;
                counter += len as usize;
            }
            println!("{} bytes read", counter);
        }
    }
    // then use the DNA
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
                // println!("{} ({})", code, len);
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
                // println!("  SDNAnr = {}", sdna_nr);
                // for now ignore the nr entry
                let mut buffer = [0; 4];
                f.read(&mut buffer)?;
                counter += 4;
                // are we done?
                if code == String::from("ENDB") {
                    break;
                }
                // read len bytes
                let mut buffer = vec![0; len as usize];
                f.read(&mut buffer)?;
                counter += len as usize;
            }
            println!("{} bytes read", counter);
        }
    }
    // WORK
    println!("number of lights = {:?}", lights.len());
    println!("number of primitives = {:?}", primitives.len());
    // let accelerator = Arc::new(BVHAccel::new(
    //     primitives.clone(),
    //     max_prims_in_node as usize,
    //     split_method,
    // ));
    // let scene: Scene = Scene::new(accelerator.clone(), lights.clone());
    // in the end we want to call render()
    // render(
    //     &scene,
    //     &camera.clone(),
    //     &mut sampler,
    //     &mut integrator,
    //     num_threads,
    // );
    Ok(())
}
