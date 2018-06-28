extern crate getopts;
extern crate num_cpus;
extern crate pbrt;
// pest
extern crate pest;
#[macro_use]
extern crate pest_derive;

// parser
use pest::Parser;

// getopts
use getopts::Options;
// pbrt
use pbrt::core::api::ApiState;
use pbrt::core::api::{
    pbrt_area_light_source, pbrt_attribute_begin, pbrt_attribute_end, pbrt_camera, pbrt_cleanup,
    pbrt_film, pbrt_init, pbrt_integrator, pbrt_look_at, pbrt_make_named_material,
    pbrt_named_material, pbrt_sampler, pbrt_scale, pbrt_shape, pbrt_transform, pbrt_world_begin,
};
use pbrt::core::geometry::Point3f;
use pbrt::core::paramset::ParamSet;
use pbrt::core::pbrt::{Float, Spectrum};
use pbrt::core::transform::Transform;
// std
use std::env;
use std::fs::File;
use std::io::BufReader;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::str::FromStr;

pub const VERSION: &'static str = env!("CARGO_PKG_VERSION");

static mut NUMBER_OF_THREADS: u8 = 0_u8;
static mut SEARCH_DIRECTORY: Option<Box<PathBuf>> = None;

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

fn pbrt_float_parameter(pairs: &mut pest::iterators::Pairs<Rule>) -> (String, Vec<Float>) {
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

fn pbrt_integer_parameter(pairs: &mut pest::iterators::Pairs<Rule>) -> (String, Vec<i32>) {
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

fn pbrt_string_parameter(pairs: &mut pest::iterators::Pairs<Rule>) -> (String, String) {
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

fn extract_params(key_word: String, pairs: pest::iterators::Pair<Rule>) -> ParamSet {
    let mut name: String = String::from("");
    let mut params: ParamSet = ParamSet::default();
    params.key_word = key_word;
    for pair in pairs.into_inner() {
        // let span = pair.clone().into_span();
        // println!("Rule:    {:?}", pair.as_rule());
        // println!("Span:    {:?}", span);
        // println!("Text:    {}", span.as_str());
        match pair.as_rule() {
            Rule::string => {
                let mut string_pairs = pair.into_inner();
                let ident = string_pairs.next();
                name = String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
                params.name = name;
            }
            Rule::parameter => {
                for parameter_pair in pair.into_inner() {
                    match parameter_pair.as_rule() {
                        Rule::float_param => {
                            let tuple: (String, Vec<Float>) =
                                pbrt_float_parameter(&mut parameter_pair.into_inner());
                            let string: String = tuple.0;
                            let floats: Vec<Float> = tuple.1;
                            if floats.len() == 1 {
                                params.add_float(string, floats[0]);
                            } else {
                                params.add_floats(string, floats);
                            }
                        }
                        Rule::integer_param => {
                            let tuple: (String, Vec<i32>) =
                                pbrt_integer_parameter(&mut parameter_pair.into_inner());
                            let string: String = tuple.0;
                            let integers: Vec<i32> = tuple.1;
                            if integers.len() == 1 {
                                params.add_int(string, integers[0]);
                            } else {
                                params.add_ints(string, integers);
                            }
                        }
                        Rule::point_param => {
                            let tuple: (String, Vec<Float>) =
                                pbrt_float_parameter(&mut parameter_pair.into_inner());
                            let string: String = tuple.0;
                            let floats: Vec<Float> = tuple.1;
                            if floats.len() == 3 {
                                params.add_point3f(
                                    string,
                                    Point3f {
                                        x: floats[0],
                                        y: floats[1],
                                        z: floats[2],
                                    },
                                );
                            } else {
                                params.add_point3fs(string, floats);
                            }
                        }
                        Rule::rgb_param => {
                            let tuple: (String, Vec<Float>) =
                                pbrt_float_parameter(&mut parameter_pair.into_inner());
                            let string: String = tuple.0;
                            let floats: Vec<Float> = tuple.1;
                            params.add_rgb_spectrum(
                                string,
                                Spectrum {
                                    c: [floats[0], floats[1], floats[2]],
                                },
                            );
                        }
                        Rule::string_param => {
                            let tuple: (String, String) =
                                pbrt_string_parameter(&mut parameter_pair.into_inner());
                            let string1: String = tuple.0;
                            let string2: String = tuple.1;
                            params.add_string(string1, string2);
                        }
                        // TODO: more rules
                        _ => println!("TODO: {:?}", parameter_pair.as_rule()),
                    }
                }
            }
            _ => println!("TODO: {:?}", pair.as_rule()),
        }
    }
    params
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
                let num_cores = num_cpus::get();
                println!("pbrt version {} [Detected {} cores]", VERSION, num_cores);
                println!("Copyright (c)2016-2018 Jan Douglas Bert Walter.");
                println!(
                    "Rust code based on C++ code by Matt Pharr, Greg Humphreys, and Wenzel Jakob."
                );
                // println!("FILE = {}", x);
                let f = File::open(x.clone()).unwrap();
                let ip: &Path = Path::new(x.as_str());
                if ip.is_relative() {
                    let cp: PathBuf = env::current_dir().unwrap();
                    let pb: PathBuf = cp.join(ip);
                    let search_directory: &Path = pb.as_path().parent().unwrap();
                    // println!("search_directory is {}", search_directory.display());
                    unsafe {
                        SEARCH_DIRECTORY = Some(Box::new(PathBuf::from(search_directory)));
                    }
                }
                let mut reader = BufReader::new(f);
                let mut str_buf: String = String::default();
                let _num_bytes = reader.read_to_string(&mut str_buf);
                // if num_bytes.is_ok() {
                //     let n_bytes = num_bytes.unwrap();
                //     println!("{} bytes read", n_bytes);
                // }
                // parser
                let pairs =
                    PbrtParser::parse(Rule::pbrt, &str_buf).unwrap_or_else(|e| panic!("{}", e));
                // println!("do something with created tokens ...");
                let mut api_state: ApiState = pbrt_init();
                for pair in pairs {
                    // let span = pair.clone().into_span();
                    // println!("Rule:    {:?}", pair.as_rule());
                    // println!("Span:    {:?}", span);
                    // println!("Text:    {}", span.as_str());
                    for inner_pair in pair.into_inner() {
                        match inner_pair.as_rule() {
                            Rule::keyword => {
                                for rule_pair in inner_pair.into_inner() {
                                    match rule_pair.as_rule() {
                                        Rule::attribute_begin => {
                                            pbrt_attribute_begin(&mut api_state);
                                        }
                                        Rule::attribute_end => {
                                            pbrt_attribute_end(&mut api_state);
                                        }
                                        Rule::world_begin => {
                                            pbrt_world_begin(&mut api_state);
                                        }
                                        _ => println!("TODO: {:?}", rule_pair.as_rule()),
                                    }
                                }
                            }
                            Rule::look_at => {
                                // LookAt eye_x eye_y eye_z look_x look_y look_z up_x up_y up_z
                                let mut v: Vec<Float> = Vec::new();
                                for rule_pair in inner_pair.into_inner() {
                                    let number: Float = f32::from_str(
                                        rule_pair.clone().into_span().as_str(),
                                    ).unwrap();
                                    v.push(number);
                                }
                                pbrt_look_at(&mut api_state, v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]);
                            }
                            Rule::named_statement => {
                                for rule_pair in inner_pair.into_inner() {
                                    match rule_pair.as_rule() {
                                        Rule::area_light_source => {
                                            let params = extract_params(String::from("AreaLightSource"), rule_pair);
                                            pbrt_area_light_source(&mut api_state, params);
                                        }
                                        Rule::camera => {
                                            let params = extract_params(String::from("Camera"), rule_pair);
                                            pbrt_camera(&mut api_state, params);
                                        }
                                        Rule::film => {
                                            let params = extract_params(String::from("Film"), rule_pair);
                                            pbrt_film(&mut api_state, params);
                                        }
                                        Rule::integrator => {
                                            let params = extract_params(String::from("Integrator"), rule_pair);
                                            pbrt_integrator(&mut api_state, params);
                                        }
                                        Rule::make_named_material => {
                                            let params = extract_params(String::from("MakeNamedMaterial"), rule_pair);
                                            pbrt_make_named_material(&mut api_state, params);
                                        }
                                        Rule::named_material => {
                                            let params = extract_params(String::from("NamedMaterial"), rule_pair);
                                            pbrt_named_material(&mut api_state, params);
                                        }
                                        Rule::sampler => {
                                            let params = extract_params(String::from("Sampler"), rule_pair);
                                            pbrt_sampler(&mut api_state, params);
                                        }
                                        Rule::shape => {
                                            let params = extract_params(String::from("Shape"), rule_pair);
                                            pbrt_shape(&mut api_state, params);
                                        }
                                        _ => println!("TODO: {:?}", rule_pair.as_rule()),
                                    }
                                }
                            }
                            Rule::scale => {
                                // Scale x y z
                                let mut v: Vec<Float> = Vec::new();
                                for rule_pair in inner_pair.into_inner() {
                                    let number: Float = f32::from_str(
                                        rule_pair.clone().into_span().as_str(),
                                    ).unwrap();
                                    v.push(number);
                                }
                                pbrt_scale(&mut api_state, v[0], v[1], v[2]);
                            }
                            Rule::transform => {
                                let mut m: Vec<Float> = Vec::new();
                                for rule_pair in inner_pair.into_inner() {
                                    // ignore brackets
                                    let not_opening: bool = rule_pair.as_str() != String::from("[");
                                    let not_closing: bool = rule_pair.as_str() != String::from("]");
                                    if not_opening && not_closing {
                                        let number: Float =
                                            f32::from_str(rule_pair.clone().into_span().as_str())
                                                .unwrap();
                                        m.push(number);
                                    }
                                }
                                let m00: Float = m[0];
                                let m01: Float = m[1];
                                let m02: Float = m[2];
                                let m03: Float = m[3];
                                let m10: Float = m[4];
                                let m11: Float = m[5];
                                let m12: Float = m[6];
                                let m13: Float = m[7];
                                let m20: Float = m[8];
                                let m21: Float = m[9];
                                let m22: Float = m[10];
                                let m23: Float = m[11];
                                let m30: Float = m[12];
                                let m31: Float = m[13];
                                let m32: Float = m[14];
                                let m33: Float = m[15];
                                let tr: Transform = Transform::new(
                                    m00, m10, m20, m30, m01, m11, m21, m31, m02, m12, m22, m32,
                                    m03, m13, m23, m33,
                                );
                                pbrt_transform(&mut api_state, &tr);
                            }
                            // WORK
                            _ => println!("TODO: {:?}", inner_pair.as_rule()),
                        };
                    }
                }
                pbrt_cleanup();
                // println!("done.");
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
