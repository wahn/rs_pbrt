// extern crate pbrt;
extern crate getopts;
extern crate num_cpus;
// pest
extern crate pest;
#[macro_use]
extern crate pest_derive;

// parser
use pest::Parser;

// use pbrt::{Bounds2, Point2};

// getopts
use getopts::Options;
// std
use std::env;
use std::fs::File;
use std::io::BufReader;
use std::io::Read;
use std::path::{Path, PathBuf};

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
                    // WORK
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
