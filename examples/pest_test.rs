#![recursion_limit="2000"]

#[macro_use]
extern crate pest;
extern crate getopts;

// parser
use pest::prelude::*;
// getopts
use getopts::Options;
// std
use std::collections::LinkedList;
use std::env;
use std::fs::File;
use std::io::prelude::*;
use std::io;
use std::io::BufReader;

pub const VERSION: &'static str = env!("CARGO_PKG_VERSION");

#[derive(Debug, PartialEq)]
pub enum Node {
    Sentence(LinkedList<Node>),
    Ident(char),
}

impl_rdp! {
    grammar! {
        pbrt = _{ statement* ~ last_statement }
        statement = { keyword ~ number* }
        keyword = {
            (["Accelerator"] |
             ["ActiveTransform"] |
             ["All"] |
             ["AreaLightSource"] |
             ["AttributeBegin"] |
             ["AttributeEnd"] |
             ["Camera"] |
             ["ConcatTransform"] |
             ["CoordinateSystem"] |
             ["CoordSysTransform"] |
             ["EndTime"] |
             ["Film"] |
             ["Identity"] |
             ["Include"] |
             ["LightSource"] |
             ["LookAt"] |
             ["MakeNamedMedium"] |
             ["MakeNamedMaterial"] |
             ["Material"] |
             ["MediumInterface"] |
             ["NamedMaterial"] |
             ["ObjectBegin"] |
             ["ObjectEnd"] |
             ["ObjectInstance"] |
             ["PixelFilter"] |
             ["ReverseOrientation"] |
             ["Rotate"] |
             ["Sampler"] |
             ["Scale"] |
             ["Shape"] |
             ["StartTime"] |
             ["Integrator"] |
             ["Texture"] |
             ["TransformBegin"] |
             ["TransformEnd"] |
             ["TransformTimes"] |
             ["Transform"] |
             ["Translate"] |
             ["WorldBegin"])
        }
        // NUMBER [-+]?([0-9]+|(([0-9]+\.[0-9]*)|(\.[0-9]+)))([eE][-+]?[0-9]+)?
        number = @{
            (["-"] | ["+"])? ~ // optional sign, followed by
            (
                (
                    (["."] ~ ['0'..'9']+) // dot and digits
                        | // or
                    (['0'..'9']+ ~ ["."] ~ ['0'..'9']*) // digits, dot, and (optional digits)
                )
                    | // or
                ['0'..'9']+ // just digits
            ) ~ ( // followed by (optional)
                (["e"] | ["E"]) ~ // 'e' or 'E', followed by
                (["-"] | ["+"])? ~ // optional sign, followed by
                ['0'..'9']+ // digits
            )?
        }
        last_statement = { ["WorldEnd"] ~ whitespace? }
        whitespace = _{ ([" "] | ["\t"] | ["\r"] | ["\n"]) }
        // // "[" { return LBRACK; }
        // lbrack = { ["["] }
        // // "]" { return RBRACK; }
        // rbrack = { ["]"] }
        // // IDENT [a-zA-Z_][a-zA-Z_0-9]*
        // ident =  { (['a'..'z'] | ['A'..'Z'] | ["_"]) ~
        //            (['a'..'z'] | ['A'..'Z'] | ["_"] | ['0'..'9'])* }
        // string = { ["\""] ~ ident ~ ["\""] }
        // string_two_words = { ["\""] ~ ident ~ whitespaces ~ ident ~ ["\""] }
        // num_list = { number ~ (whitespaces ~ number)* }
        // // num_array: array_init LBRACK num_list RBRACK
        // num_array = { lbrack ~ whitespaces? ~ num_list ~ whitespaces? ~ rbrack }
        // // paramlist_entry: STRING array
        // paramlist_entry = { string_two_words ~ whitespaces ~ num_array }
        // // LOOKAT NUM NUM NUM NUM NUM NUM NUM NUM NUM
        // look_at = { ["LookAt"] ~ whitespaces ~
        //               number ~ whitespaces ~ number ~ whitespaces ~ number ~ whitespaces ~
        //               number ~ whitespaces ~ number ~ whitespaces ~ number ~ whitespaces ~
        //               number ~ whitespaces ~ number ~ whitespaces ~ number
        // }
        // // CAMERA STRING paramlist
        // camera = { ["Camera"] ~ whitespaces ~ string ~ whitespaces ~ paramlist_entry }
    }
}

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
                let f = File::open(x).unwrap();
                let mut reader = BufReader::new(f);
                let mut str_buf: String = String::default();
                reader.read_to_string(&mut str_buf);
                // parser
                let mut parser = Rdp::new(StringInput::new(&str_buf));
                assert!(parser.pbrt());
                assert!(parser.end());
                println!("{:?}", parser.queue());
            }
            None => panic!("no input file name"),
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
