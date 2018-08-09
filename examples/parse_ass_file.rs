extern crate getopts;
// pest
extern crate pest;
#[macro_use]
extern crate pest_derive;

// parser
use pest::Parser;

// getopts
use getopts::Options;
// std
use std::env;
use std::fs::File;
use std::io::BufReader;
use std::io::Read;
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[grammar = "../examples/ass.pest"]
struct PbrtParser;

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
