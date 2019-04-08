extern crate getopts;
extern crate pbrt;
// pest
extern crate pest;
#[macro_use]
extern crate pest_derive;

pub const VERSION: &'static str = env!("CARGO_PKG_VERSION");

#[derive(Parser)]
#[grammar = "../examples/parse_pbrt.pest"]
struct PbrtParser;

// parser
use pest::Parser;

// getopts
use getopts::Options;
// pbrt
use pbrt::core::api::{pbrt_cleanup, pbrt_init, pbrt_world_begin};
// std
use std::env;
use std::fs::File;
use std::io::BufReader;
use std::io::Read;
use std::path::{Path, PathBuf};

fn print_usage(program: &str, opts: Options) {
    let brief = format!("Usage: {} [options]", program);
    print!("{}", opts.usage(&brief));
}

fn print_version(program: &str) {
    println!("{} {}", program, VERSION);
}

// AttributeBegin
// AttributeEnd
// ActiveTransform
// AreaLightSource
// Accelerator
// ConcatTransform
// CoordinateSystem
// CoordSysTransform
// Camera
// Integrator
// Include
// Identity
// LightSource
// LookAt
// MakeNamedMaterial
// MakeNamedMedium
// Material
// MediumInterface
// NamedMaterial
// ObjectBegin
// ObjectEnd
// ObjectInstance
// PixelFilter
// ReverseOrientation
// Rotate
// Shape
// Sampler
// Scale
// TransformBegin
// TransformEnd
// Transform
// Translate
// TransformTimes
// Texture

fn main() {
    // handle command line options
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();
    let mut opts = Options::new();
    opts.optflag("h", "help", "print this help menu");
    opts.optopt("i", "", "parse an input file", "FILE");
    opts.optopt(
        "t",
        "nthreads",
        "use specified number of threads for rendering",
        "NUM",
    );
    opts.optflag("v", "version", "print version number");
    let matches = match opts.parse(&args[1..]) {
        Ok(m) => m,
        Err(f) => panic!(f.to_string()),
    };
    if matches.opt_present("h") {
        print_usage(&program, opts);
        return;
    } else if matches.opt_present("i") {
        let mut number_of_threads: u8 = 0_u8;
        if matches.opt_present("t") {
            let nthreads = matches.opt_str("t");
            match nthreads {
                Some(x) => {
                    let number_result = x.parse::<u8>();
                    assert!(
                        !number_result.is_err(),
                        "ERROR: 8 bit unsigned integer expected"
                    );
                    let num_threads: u8 = number_result.unwrap();
                    println!("nthreads = {:?}", num_threads);
                    number_of_threads = num_threads;
                }
                None => panic!("No argument for number of threads given."),
            }
        }
        let infile = matches.opt_str("i");
        match infile {
            Some(x) => {
                let num_cores = num_cpus::get();
                println!("pbrt version {} [Detected {} cores]", VERSION, num_cores);
                println!("Copyright (c) 2016-2019 Jan Douglas Bert Walter.");
                println!(
                    "Rust code based on C++ code by Matt Pharr, Greg Humphreys, and Wenzel Jakob."
                );
                println!("FILE = {}", x);
                let f = File::open(x.clone()).unwrap();
                let ip: &Path = Path::new(x.as_str());
                let (mut api_state, mut bsdf_state) = pbrt_init(number_of_threads);
                if ip.is_relative() {
                    let cp: PathBuf = env::current_dir().unwrap();
                    let pb: PathBuf = cp.join(ip);
                    let search_directory: &Path = pb.as_path().parent().unwrap();
                    // println!("search_directory is {}", search_directory.display());
                    api_state.search_directory = Some(Box::new(PathBuf::from(search_directory)));
                }
                let mut reader = BufReader::new(f);
                let mut str_buf: String = String::default();
                let num_bytes = reader.read_to_string(&mut str_buf);
                if num_bytes.is_ok() {
                    let n_bytes = num_bytes.unwrap();
                    println!("{} bytes read", n_bytes);
                }
                let pairs = PbrtParser::parse(Rule::pbrt, &str_buf)
                    .expect("unsuccessful parse")
                    .next()
                    .unwrap();
                let mut identifier: &str = "";
                let mut comment_count: u64 = 0;
                let mut empty_count: u64 = 0;
                let mut todo_count: u64 = 0;
                let mut parse_again: String = String::default();
                // first parse file line by line
                for inner_pair in pairs.into_inner() {
                    match inner_pair.as_rule() {
                        // comment lines (starting with '#')
                        Rule::comment_line => {
                            comment_count += 1;
                        }
                        Rule::statement_line => {
                            for statement_pair in inner_pair.into_inner() {
                                match statement_pair.as_rule() {
                                    Rule::identifier => {
                                        if identifier != "" {
                                            println!("{} {}", identifier, parse_again);
                                        }
                                        identifier = statement_pair.as_str();
                                        parse_again = String::default();
                                    }
                                    Rule::remaining_line => {
                                        if parse_again != "" {
                                            parse_again =
                                                parse_again + " " + statement_pair.as_str();
                                        } else {
                                            parse_again += statement_pair.as_str();
                                        }
                                    }
                                    _ => println!("TODO: {:?}", statement_pair.as_rule()),
                                }
                            }
                        }
                        Rule::empty_line => {
                            empty_count += 1;
                        }
                        Rule::todo_line => {
                            todo_count += 1;
                            for params_pair in inner_pair.into_inner() {
                                match params_pair.as_rule() {
                                    Rule::remaining_params => {
                                        if parse_again != "" {
                                            parse_again = parse_again + " " + params_pair.as_str();
                                        } else {
                                            parse_again += params_pair.as_str();
                                        }
                                    }
                                    _ => println!("TODO: {:?}", params_pair.as_rule()),
                                }
                            }
                        }
                        Rule::EOI => (),
                        _ => unreachable!(),
                    }
                }
                println!("Number of comment line(s):   {}", comment_count);
                println!("Number of line(s) left TODO: {}", todo_count);
                println!("Number of empty line(s):     {}", empty_count);
                if todo_count == 0 {
                    // usually triggered by WorldEnd
                    pbrt_cleanup(&mut api_state);
                    // TODO: rename pbrt_cleanup() to pbrt_world_end()?
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
