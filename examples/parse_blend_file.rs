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
    let mut primitives: Vec<Arc<Primitive + Sync + Send>> = Vec::new();
    let mut lights: Vec<Arc<Light + Sync + Send>> = Vec::new();
    // first get the DNA
    let mut names: Vec<String> = Vec::new();
    let mut names_len: usize = 0;
    let mut types: Vec<String> = Vec::new();
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
                let mut sdna_nr: u32 = 0;
                sdna_nr += (buffer[0] as u32) << 0;
                sdna_nr += (buffer[1] as u32) << 8;
                sdna_nr += (buffer[2] as u32) << 16;
                sdna_nr += (buffer[3] as u32) << 24;
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
                                    for s in 0..nr_structs as usize {
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
                                        let tuple_counter: usize = short2 as usize;
                                        for _t in 0..tuple_counter {
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
                                            let mut name_idx: u16 = 0;
                                            name_idx += (buffer[0] as u16) << 0;
                                            name_idx += (buffer[1] as u16) << 8;
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
