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
    Word(LinkedList<Node>),
    Letter(char),
}

impl_rdp! {
    grammar! {
        sentence = _{ word ~ ([" "] ~ word)* }
        word     =  { letter* }
        letter   =  { ['a'..'z'] }
    }

    process! {
        main(&self) -> Node {
            (list: _sentence()) => {
                Node::Sentence(list)
            }
        }

        _sentence(&self) -> LinkedList<Node> {
            (_: word, head: _word(), mut tail: _sentence()) => {
                tail.push_front(Node::Word(head));

                tail
            },
            () => {
                LinkedList::new()
            }
        }

        _word(&self) -> LinkedList<Node> {
            (&head: letter, mut tail: _word()) => {
                tail.push_front(Node::Letter(head.chars().next().unwrap()));

                tail
            },
            () => {
                LinkedList::new()
            }
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
                assert!(parser.sentence());
                println!("{:?}", parser.main());
            }
            None => panic!("no input file name"),
        }
        return;
    } else if matches.opt_present("v") {
        print_version(&program);
        return;
    }
    // parser
    let mut parser = Rdp::new(StringInput::new("abc def"));

    assert!(parser.sentence());
    println!("{:?}", parser.main());
}
