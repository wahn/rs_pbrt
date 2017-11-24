extern crate ply_rs;

// std
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;
use std::string::String;
use std::sync::Arc;
use std::vec::Vec;
// others
use ply_rs::parser;
use ply_rs::ply;
// pbrt
use core::paramset::ParamSet;
use core::pbrt::Float;
use core::transform::Transform;
use geometry::{Normal3f, Point2f, Point3f};
use shapes::Shape;
use textures::Texture;

struct Vertex {
    p: Point3f,
    n: Normal3f,
    uv: Point2f,
}

impl ply::PropertyAccess for Vertex {
    fn new() -> Vertex {
        Vertex {
            p: Point3f::default(),
            n: Normal3f::default(),
            uv: Point2f::default(),
        }
    }
}

struct Face {
    vertex_indices: Vec<i32>,
}

impl ply::PropertyAccess for Face {
    fn new() -> Face {
        Face { vertex_indices: Vec::new() }
    }
}

pub fn create_ply_mesh(o2w: Transform,
                       w2o: Transform,
                       reverse_orientation: bool,
                       params: &ParamSet,
                       float_textures: HashMap<String, Arc<Texture<Float> + Send + Sync>>,
                       search_directory: Option<&Box<PathBuf>>)
                       -> Vec<Arc<Shape + Send + Sync>> {
    let mut filename: String = params.find_one_string(String::from("filename"), String::new());
    if let Some(ref search_directory) = search_directory {
        let mut path_buf: PathBuf = PathBuf::from("/");
        path_buf.push(search_directory.as_ref());
        path_buf.push(filename);
        filename = String::from(path_buf.to_str().unwrap());
    }
    // header
    let result = File::open(&filename);
    if result.is_err() {
        panic!("Couldn't open PLY file {:?}", filename);
    }
    let f = result.unwrap();
    let mut f = BufReader::new(f);
    let vertex_parser = parser::Parser::<Vertex>::new();
    let face_parser = parser::Parser::<Face>::new();
    let result = vertex_parser.read_header(&mut f);
    if result.is_err() {
        panic!("Unable to read the header of PLY file  {:?}", filename);
    }
    let header = result.unwrap();
    println!("header = {:?}", header);
    // WORK
    Vec::new()
}
