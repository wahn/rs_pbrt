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
use pbrt::core::film::Film;
use pbrt::core::filter::Filter;
use pbrt::core::geometry::{Bounds2f, Point2f, Point2i};
use pbrt::core::paramset::ParamSet;
use pbrt::core::pbrt::Float;
use pbrt::filters::gaussian::GaussianFilter;
// std
use std::env;
use std::fs::File;
use std::io::BufReader;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::sync::Arc;

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
        // default values
        let mut filter_name: String = String::from("box");
        let mut filter_width: Float = 2.0;
        let mut xres: i32 = 1280;
        let mut yres: i32 = 720;
        // input (.ass) file
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
                                            } else if node_type == String::from("gaussian_filter") {
                                                filter_name = String::from("gaussian");
                                                if next == String::from("width") {
                                                    if let Some(filter_width_str) = iter.next() {
                                                        filter_width =
                                                            f32::from_str(filter_width_str)
                                                                .unwrap();
                                                        print!("\n filter_width {} ", filter_width);
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
        println!("filter_name = {:?}", filter_name);
        println!("filter_width = {:?}", filter_width);
        // MakeFilter
        let mut some_filter: Option<Arc<Filter + Sync + Send>> = None;
        if filter_name == String::from("box") {
            println!("TODO: CreateBoxFilter");
        } else if filter_name == String::from("gaussian") {
            let mut filter_params: ParamSet = ParamSet::default();
            filter_params.add_float(String::from("xwidth"), filter_width);
            filter_params.add_float(String::from("ywidth"), filter_width);
            some_filter = Some(GaussianFilter::create(&filter_params));
        } else if filter_name == String::from("mitchell") {
            println!("TODO: CreateMitchellFilter");
        } else if filter_name == String::from("sinc") {
            println!("TODO: CreateSincFilter");
        } else if filter_name == String::from("triangle") {
            println!("TODO: CreateTriangleFilter");
        } else {
            panic!("Filter \"{}\" unknown.", filter_name);
        }
        // MakeFilm
        let resolution: Point2i = Point2i { x: xres, y: yres };
        println!("resolution = {:?}", resolution);
        if let Some(filter) = some_filter {
            let crop: Bounds2f = Bounds2f {
                p_min: Point2f { x: 0.0, y: 0.0 },
                p_max: Point2f { x: 1.0, y: 1.0 },
            };
            let diagonal: Float = 35.0;
            let scale: Float = 1.0;
            let max_sample_luminance: Float = std::f32::INFINITY;
            let film: Arc<Film> = Arc::new(Film::new(
                resolution,
                crop,
                filter,
                diagonal,
                String::from(""),
                scale,
                max_sample_luminance,
            ));
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
