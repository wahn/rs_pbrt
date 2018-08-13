extern crate getopts;
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
use pbrt::core::geometry::Point2i;
// std
use std::env;
use std::fs::File;
use std::io::BufReader;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::str::FromStr;

#[derive(Parser)]
#[grammar = "../examples/ass.pest"]
struct AssParser;

fn print_usage(program: &str, opts: Options) {
    let brief = format!("Usage: {} [options]", program);
    print!("{}", opts.usage(&brief));
}

pub const VERSION: &'static str = env!("CARGO_PKG_VERSION");

fn print_version(program: &str) {
    println!("{} {}", program, VERSION);
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();
    let mut opts = Options::new();
    opts.optflag("h", "help", "print this help menu");
    opts.optopt("i", "", "parse an input file", "FILE");
    opts.optflag("v", "version", "print version number");
    let matches = match opts.parse(&args[1..]) {
        Ok(m) => m,
        Err(f) => panic!(f.to_string()),
    };
    if matches.opt_present("h") {
        print_usage(&program, opts);
        return;
    } else if matches.opt_present("i") {
        let mut xres: i32 = 1280;
        let mut yres: i32 = 720;
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
                    AssParser::parse(Rule::ass, &str_buf).unwrap_or_else(|e| panic!("{}", e));
                // let tokens: Vec<_> = pairs.flatten().tokens().collect();
                // println!("{} pairs", tokens.len());
                for pair in pairs {
                    let span = pair.clone().into_span();
                    // println!("Rule:    {:?}", pair.as_rule());
                    // println!("Span:    {:?}", span);
                    // println!("Text:    {}", span.as_str());
                    for inner_pair in pair.into_inner() {
                        match inner_pair.as_rule() {
                            Rule::ident => {
                                let node_type = inner_pair.clone().into_span().as_str();
                                print!("{} {{", node_type);
                                let mut iter = span.as_str().split_whitespace();
                                loop {
                                    if let Some(next) = iter.next() {
                                        if next != String::from("}") {
                                            if next == String::from("name") {
                                                if let Some(name) = iter.next() {
                                                    print!(" {} {} ", next, name);
                                                }
                                            } else if node_type == String::from("options") {
                                                if next == String::from("xres") {
                                                    if let Some(xres_str) = iter.next() {
                                                        xres = i32::from_str(xres_str).unwrap();
                                                        print!("\n xres {} ", xres);
                                                    }
                                                } else if next == String::from("yres") {
                                                    if let Some(yres_str) = iter.next() {
                                                        yres = i32::from_str(yres_str).unwrap();
                                                        print!("\n yres {} ", yres);
                                                    }
                                                }
                                            }
                                        } else {
                                            println!("}}");
                                        }
                                    } else {
                                        break;
                                    }
                                }
                            }
                            // WORK
                            _ => println!("TODO: {:?}", inner_pair.as_rule()),
                        }
                    }
                }
            }
            None => panic!("No input file name."),
        }
        let resolution: Point2i = Point2i { x: xres, y: yres };
        println!("resolution = {:?}", resolution);
        return;
    } else if matches.opt_present("v") {
        print_version(&program);
        return;
    } else {
        print_usage(&program, opts);
        return;
    }
}
