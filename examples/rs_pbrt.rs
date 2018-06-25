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
use pbrt::core::api::{
    pbrt_attribute_begin, pbrt_attribute_end, pbrt_look_at, pbrt_scale, pbrt_transform,
    pbrt_world_begin,
};
use pbrt::core::pbrt::Float;
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
                // parser
                let pairs =
                    PbrtParser::parse(Rule::pbrt, &str_buf).unwrap_or_else(|e| panic!("{}", e));
                println!("do something with created tokens ...");
                for pair in pairs {
                    let span = pair.clone().into_span();
                    println!("Rule:    {:?}", pair.as_rule());
                    println!("Span:    {:?}", span);
                    println!("Text:    {}", span.as_str());
                    for inner_pair in pair.into_inner() {
                        match inner_pair.as_rule() {
                            Rule::keyword => {
                                for rule_pair in inner_pair.into_inner() {
                                    match rule_pair.as_rule() {
                                        Rule::attribute_begin => {
                                            pbrt_attribute_begin();
                                        }
                                        Rule::attribute_end => {
                                            pbrt_attribute_end();
                                        }
                                        Rule::world_begin => {
                                            pbrt_world_begin();
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
                                pbrt_look_at(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]);
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
                                pbrt_scale(v[0], v[1], v[2]);
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
                                pbrt_transform(tr);
                            }
                            // WORK
                            _ => println!("TODO: {:?}", inner_pair.as_rule()),
                        };
                    }
                }
                println!("done.");
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
