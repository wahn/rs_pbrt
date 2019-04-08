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
use pbrt::core::api::ApiState;
use pbrt::core::api::{
    pbrt_attribute_begin, pbrt_attribute_end, pbrt_cleanup, pbrt_init, pbrt_look_at, pbrt_rotate,
    pbrt_translate, pbrt_world_begin,
};
use pbrt::core::pbrt::Float;
// std
use std::env;
use std::fs::File;
use std::io::BufReader;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::str::FromStr;

fn print_usage(program: &str, opts: Options) {
    let brief = format!("Usage: {} [options]", program);
    print!("{}", opts.usage(&brief));
}

fn print_version(program: &str) {
    println!("{} {}", program, VERSION);
}

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

fn parse_line(api_state: &mut ApiState, identifier: &str, str_buf: String) {
    if str_buf == "" {
        // no additional arguments
        match identifier {
            "AttributeBegin" => {
                // AttributeBegin
                println!("{} {}", identifier, str_buf);
                pbrt_attribute_begin(api_state);
            }
            "AttributeEnd" => {
                // AttributeEnd
                println!("{} {}", identifier, str_buf);
                pbrt_attribute_end(api_state);
            }
            "WorldBegin" => {
                // WorldBegin
                println!("{} {}", identifier, str_buf);
                pbrt_world_begin(api_state);
            }
            "WorldEnd" => {
                // WorldEnd
                println!("{} {}", identifier, str_buf);
                // TODO: pbrt_cleanup(api_state);
            }
            _ => println!("{} {}", identifier, str_buf),
        }
    } else {
        let pairs = PbrtParser::parse(Rule::name_and_or_params, &str_buf)
            .expect("unsuccessful parse")
            .next()
            .unwrap();
        for inner_pair in pairs.into_inner() {
            match inner_pair.as_rule() {
                Rule::type_params => {
                    // identifier "type" parameter-list
                    println!("> {} {}", identifier, inner_pair.as_str());
                }
                Rule::look_at => {
                    // LookAt eye_x eye_y eye_z look_x look_y look_z up_x up_y up_z
                    let mut v: Vec<Float> = Vec::new();
                    for rule_pair in inner_pair.into_inner() {
                        let number: Float =
                            f32::from_str(rule_pair.clone().into_span().as_str()).unwrap();
                        v.push(number);
                    }
                    println!(
                        "LookAt {} {} {} {} {} {} {} {} {}",
                        v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8],
                    );
                    pbrt_look_at(
                        api_state, v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8],
                    );
                }
                Rule::rotate => {
                    // Rotate angle x y z
                    let mut v: Vec<Float> = Vec::new();
                    for rule_pair in inner_pair.into_inner() {
                        let number: Float =
                            f32::from_str(rule_pair.clone().into_span().as_str()).unwrap();
                        v.push(number);
                    }
                    println!("Rotate {} {} {} {}", v[0], v[1], v[2], v[3]);
                    pbrt_rotate(api_state, v[0], v[1], v[2], v[3]);
                }
                Rule::translate => {
                    // Translate x y z
                    let mut v: Vec<Float> = Vec::new();
                    for rule_pair in inner_pair.into_inner() {
                        let number: Float =
                            f32::from_str(rule_pair.clone().into_span().as_str()).unwrap();
                        v.push(number);
                    }
                    println!("Translate {} {} {}", v[0], v[1], v[2]);
                    pbrt_translate(api_state, v[0], v[1], v[2]);
                }
                Rule::remaining_line => {
                    // predetermined number of arguments of predetermined type
                    println!("< {} {}", identifier, inner_pair.as_str());
                }
                // _ => unreachable!(),
                _ => println!("TODO: {:?}", inner_pair.as_rule()),
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
                let (mut api_state, mut _bsdf_state) = pbrt_init(number_of_threads);
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
                                            parse_line(
                                                &mut api_state,
                                                identifier,
                                                parse_again.clone(),
                                            );
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
                                    Rule::trailing_comment => {
                                        // ignore
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
                                    Rule::trailing_comment => {
                                        // ignore
                                    }
                                    _ => println!("TODO: {:?}", params_pair.as_rule()),
                                }
                            }
                        }
                        Rule::EOI => parse_line(&mut api_state, identifier, parse_again.clone()),
                        _ => unreachable!(),
                    }
                }
                println!("Number of comment line(s):   {}", comment_count);
                println!("Number of parameter line(s): {}", todo_count);
                println!("Number of empty line(s):     {}", empty_count);
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
