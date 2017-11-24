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

#[derive(Debug)]
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

#[derive(Debug)]
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
    let result = File::open(&filename);
    if result.is_err() {
        panic!("Couldn't open PLY file {:?}", filename);
    }
    let f = result.unwrap();
    let mut buf_reader = BufReader::new(f);
    let p = parser::Parser::<ply::DefaultElement>::new();
    // header
    let result = p.read_header(&mut buf_reader);
    if result.is_err() {
        panic!("Unable to read the header of PLY file  {:?}", filename);
    }
    let header = result.unwrap();
    println!("header = {:?}", header);
    // payload
    let result = p.read_payload(&mut buf_reader, &header);
    if result.is_err() {
        panic!("Unable to read the payload of PLY file  {:?}", filename);
    }
    let payload = result.unwrap();
    println!("payload = {:?}", payload);
    let mut p: Vec<Point3f> = Vec::new();
    for (name, list) in payload.into_iter() {
        println!("name = {:?}", name);
        match name.as_ref() {
            "vertex" => {
                for elem in list.into_iter() {
                    let mut pnt: Point3f = Point3f::default();
                    for (name2, list2) in elem.into_iter() {
                        match name2.as_ref() {
                            "x" => {
                                if let ply::Property::Float(x) = list2 {
                                    pnt.x = x;
                                }
                            }
                            "y" => {
                                if let ply::Property::Float(y) = list2 {
                                    pnt.y = y;
                                }
                            }
                            "z" => {
                                if let ply::Property::Float(z) = list2 {
                                    pnt.z = z;
                                }
                            }
                            _ => unreachable!(),
                        }
                    }
                    p.push(pnt);
                }
            }
            "face" => {
                // println!("list = {:?}", list);
            }
            _ => unreachable!(),
        }
    }
    for i in 0..p.len() {
        println!("{:?}: {:?}", i, p[i]);
    }
    // WORK
    Vec::new()
}
