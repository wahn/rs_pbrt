#![recursion_limit="200"]

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
        // IDENT [a-zA-Z_][a-zA-Z_0-9]*
        ident =  { ['a'..'z'] | ['A'..'Z'] | ["_"] ~
                     (['a'..'z'] | ['A'..'Z'] | ["_"] | ['0'..'9'])* }
        // NUMBER [-+]?
        //        ([0-9]+|
        //         (
        //          ([0-9]+\.[0-9]*)|
        //          (\.[0-9]+)
        //         )
        //        )
        //        ([eE][-+]?[0-9]+)?
        number = {
            (["-"] | ["+"])? ~
                (['0'..'9']+ |
                 (
                     (['0'..'9']+ ~ ["."] ~ ['0'..'9']*) |
                     (["."] ~ ['0'..'9']+)
                 )
                ) ~
                (["e"] | ["E"] ~ (["-"] | ["+"])? ~ ['0'..'9']+)?
        }
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
                assert!(parser.number());
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
