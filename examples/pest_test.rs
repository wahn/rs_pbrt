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
        statement = { look_at | rotate | named_statement | keyword }
        named_statement = { camera |
                            pixel_filter |
                            sampler |
                            film |
                            coord_sys_transform |
                            light_source |
                            texture |
                            material |
                            shape }
        parameter = { float_param |
                      string_param |
                      integer_param |
                      point_param |
                      rgb_param |
                      spectrum_param |
                      texture_param }
        float_param = { ["\"float"] ~ ident ~ ["\""] ~ lbrack ~ number ~ rbrack }
        string_param = { ["\"string"] ~ ident ~ ["\""] ~ lbrack ~ string ~ rbrack }
        integer_param = { ["\"integer"] ~ ident ~ ["\""] ~ lbrack ~ integer ~ rbrack }
        point_param = { ["\"point"] ~ ident ~ ["\""] ~ lbrack ~ number ~ number ~ number ~ rbrack }
        rgb_param = { ["\"rgb"] ~ ident ~ ["\""] ~ lbrack ~ number ~ number ~ number ~ rbrack }
        spectrum_param = { ["\"spectrum\""] ~ string }
        texture_param = { ["\"texture"] ~ ident ~ ["\""] ~ string }
        // Rotate angle x y z
        rotate = { ["Rotate"] ~
                   // followed by 4 numbers:
                   number ~ number ~ number ~ number
        }
        // LookAt eye_x eye_y eye_z look_x look_y look_z up_x up_y up_z
        look_at = { ["LookAt"] ~
                    // followed by 9 numbers:

                    // eye_x eye_y eye_z
                    number ~ number ~ number ~
                    // look_x look_y look_z
                    number ~ number ~ number ~
                    // up_x up_y up_z
                    number ~ number ~ number
        }
        // Camera "perspective" "float fov" [90] ...
        camera = { ["Camera"] ~ string ~ parameter* }
        // PixelFilter "mitchell" "float xwidth" [2] "float ywidth" [2]
        pixel_filter = { ["PixelFilter"] ~ string ~ parameter* }
        // Sampler "halton"
        sampler = { ["Sampler"] ~ string ~ parameter* }
        // Film "image" "string filename" ["..."] ...
        film = { ["Film"] ~ string ~ parameter* }
        // CoordSysTransform "camera"
        coord_sys_transform = { ["CoordSysTransform"] ~ string }
        // LightSource "point" "rgb I" [ .5 .5 .5 ]
        light_source = { ["LightSource"] ~ string ~ parameter* }
        // Texture "mydiffuse" "spectrum" "imagemap" "string filename" "image.tga"
        texture = { ["Texture"] ~ string ~ parameter* }
        // Material "matte" "texture Kd" "mydiffuse"
        material = { ["Material"] ~ string ~ parameter* }
        // Shape "sphere" "float radius" [0.25]
        shape = { ["Shape"] ~ string ~ parameter* }
        keyword = {
            (["Accelerator"] |
             ["ActiveTransform"] |
             ["All"] |
             ["AreaLightSource"] |
             ["AttributeBegin"] |
             ["AttributeEnd"] |
             ["ConcatTransform"] |
             ["CoordinateSystem"] |
             ["CoordSysTransform"] |
             ["EndTime"] |
             ["Identity"] |
             ["Include"] |
             ["MakeNamedMedium"] |
             ["MakeNamedMaterial"] |
             ["MediumInterface"] |
             ["NamedMaterial"] |
             ["ObjectBegin"] |
             ["ObjectEnd"] |
             ["ObjectInstance"] |
             ["ReverseOrientation"] |
             ["Scale"] |
             ["StartTime"] |
             ["Integrator"] |
             ["TransformBegin"] |
             ["TransformEnd"] |
             ["TransformTimes"] |
             ["Transform"] |
             ["Translate"] |
             ["WorldBegin"])
        }
        // IDENT [a-zA-Z_][a-zA-Z_0-9]*
        ident =  { (['a'..'z'] | ['A'..'Z'] | ["_"]) ~
                   (['a'..'z'] | ['A'..'Z'] | ["_"] | ['0'..'9'])* }
        string = { (["\""] ~ ident ~ ["\""]) | (["\""] ~ filename ~ ["\""]) }
        filename = { (['a'..'z'] | ['A'..'Z'] | ["_"]) ~ // TODO: can be a full path
                     (['a'..'z'] | ['A'..'Z'] | ["_"] | ["."] | ['0'..'9'])* }
        // "[" { return LBRACK; }
        lbrack = { ["["] }
        // "]" { return RBRACK; }
        rbrack = { ["]"] }
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
        integer = @{
            (["-"] | ["+"])? ~ // optional sign, followed by
            ['1'..'9'] ~ // at least one non-zero digit, followed by
            ['0'..'9']* // just digits
        }
        last_statement = { ["WorldEnd"] ~ whitespace? }
        whitespace = _{ ([" "] | ["\t"] | ["\r"] | ["\n"]) }
    }
    process! {
        main(&self) -> () {
            (_list: _pbrt()) => {}
        }
        _pbrt(&self) -> () {
            (_head: statement, _tail: _statement()) => { println!("DONE: statement"); },
            (_l: last_statement) => { println!("TODO: last_statement"); },
        }
        _statement(&self) -> () {
            (_head: look_at, _tail: _look_at()) => { println!("DONE: look_at"); },
            (_r: rotate) => { println!("TODO: rotate"); },
            (_n: named_statement) => { println!("TODO: named_statement"); },
            (_k: keyword) => { println!("TODO: keyword"); },
        }
        _look_at(&self) -> () {
            (eye_x: _number(), eye_y: _number(), eye_z: _number(),
             look_x: _number(), look_y: _number(), look_z: _number(),
             up_x: _number(), up_y: _number(), up_z: _number()) => { println!("TODO: look_at"); }
        }
        _number(&self) -> () {
            (&n: number) => {
                    print!("number = ");
                for c in n.chars() {
                    print!("{}", c);
                }
                println!("");
            },
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
                assert!(parser.pbrt());
                assert!(parser.end());
                println!("{:?}", parser.queue());
                println!("do something with created tokens ...");
                parser.main();
                println!("done.");
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
