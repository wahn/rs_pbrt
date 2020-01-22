// std
use std::collections::HashMap;
use std::convert::TryInto;
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
use crate::core::geometry::{Normal3f, Point2f, Point3f, Vector3f};
use crate::core::paramset::ParamSet;
use crate::core::pbrt::Float;
use crate::core::shape::Shape;
use crate::core::texture::Texture;
use crate::core::transform::Transform;
use crate::shapes::triangle::{Triangle, TriangleMesh};
use crate::textures::constant::ConstantTexture;

pub fn create_ply_mesh(
    o2w: &Transform,
    w2o: &Transform,
    reverse_orientation: bool,
    params: &ParamSet,
    float_textures: Arc<HashMap<String, Arc<dyn Texture<Float> + Send + Sync>>>,
    search_directory: Option<&Box<PathBuf>>,
) -> Vec<Arc<Shape>> {
    let mut filename: String = params.find_one_string("filename", String::new());
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
    // println!("header = {:?}", header);
    // payload
    let result = p.read_payload(&mut buf_reader, &header);
    if result.is_err() {
        panic!("Unable to read the payload of PLY file  {:?}", filename);
    }
    let payload = result.unwrap();
    // println!("payload = {:?}", payload);
    let mut p: Vec<Point3f> = Vec::new();
    let mut n: Vec<Normal3f> = Vec::new();
    let mut uvs: Vec<Point2f> = Vec::new();
    let mut has_normals: bool = false;
    let mut has_uvs: bool = false;
    let mut tm_vertex_indices: Vec<u32> = Vec::new();
    for (name, list) in payload.into_iter() {
        match name.as_ref() {
            "vertex" => {
                for elem in list.into_iter() {
                    let mut pnt: Point3f = Point3f::default();
                    let mut nrm: Normal3f = Normal3f::default();
                    let mut pt2: Point2f = Point2f::default();
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
                            "nx" => {
                                has_normals = true;
                                if let ply::Property::Float(x) = list2 {
                                    nrm.x = x;
                                }
                            }
                            "ny" => {
                                has_normals = true;
                                if let ply::Property::Float(y) = list2 {
                                    nrm.y = y;
                                }
                            }
                            "nz" => {
                                has_normals = true;
                                if let ply::Property::Float(z) = list2 {
                                    nrm.z = z;
                                }
                            }
                            "u" | "s" => {
                                has_uvs = true;
                                if let ply::Property::Float(x) = list2 {
                                    pt2.x = x;
                                }
                            }
                            "v" | "t" => {
                                has_uvs = true;
                                if let ply::Property::Float(y) = list2 {
                                    pt2.y = y;
                                }
                            }
                            _ => {
                                println!("name2 = {:?}", name2);
                                unreachable!();
                            }
                        }
                    }
                    p.push(pnt);
                    if has_normals {
                        n.push(nrm);
                    }
                    if has_uvs {
                        uvs.push(pt2);
                    }
                }
            }
            "face" => {
                for elem in list.into_iter() {
                    let mut nrm: Normal3f = Normal3f::default();
                    for (name2, list2) in elem.into_iter() {
                        match name2.as_ref() {
                            "vertex_indices" => {
                                if let ply::Property::ListInt(li) = list2 {
                                    let mut vertex_indices: Vec<usize> = Vec::new();
                                    for i in li.into_iter() {
                                        vertex_indices.push(i as usize);
                                    }
                                    // println!("vertex_indices = {:?}", vertex_indices);
                                    if vertex_indices.len() != 3 {
                                        if vertex_indices.len() == 4 {
                                            // handle quads (split it into 2 triangles)
                                            let v1 = vertex_indices[0];
                                            let v3 = vertex_indices[2];
                                            let v4 = vertex_indices.pop().unwrap();
                                            vertex_indices.push(v4);
                                            vertex_indices.push(v1);
                                            vertex_indices.push(v3);
                                        } else {
                                            panic!("plymesh: Ignoring face with {} vertices (only triangles and quads are supported!)",
                                                   vertex_indices.len());
                                        }
                                    }
                                    // now we can add the indices to the triangle mesh vertex indices
                                    for vi in vertex_indices {
                                        tm_vertex_indices.push(vi.try_into().unwrap());
                                    }
                                } else if let ply::Property::ListUInt(li) = list2 {
                                    let mut vertex_indices: Vec<usize> = Vec::new();
                                    for i in li.into_iter() {
                                        vertex_indices.push(i as usize);
                                    }
                                    // println!("vertex_indices = {:?}", vertex_indices);
                                    if vertex_indices.len() != 3 {
                                        if vertex_indices.len() == 4 {
                                            // handle quads (split it into 2 triangles)
                                            let v1 = vertex_indices[0];
                                            let v3 = vertex_indices[2];
                                            let v4 = vertex_indices.pop().unwrap();
                                            vertex_indices.push(v4);
                                            vertex_indices.push(v1);
                                            vertex_indices.push(v3);
                                        } else {
                                            panic!("plymesh: Ignoring face with {} vertices (only triangles and quads are supported!)",
                                                   vertex_indices.len());
                                        }
                                    }
                                    // now we can add the indices to the triangle mesh vertex indices
                                    for vi in vertex_indices {
                                        tm_vertex_indices.push(vi.try_into().unwrap());
                                    }
                                }
                            }
                            "nx" => {
                                has_normals = true;
                                if let ply::Property::Float(x) = list2 {
                                    nrm.x = x;
                                }
                            }
                            "ny" => {
                                has_normals = true;
                                if let ply::Property::Float(y) = list2 {
                                    nrm.y = y;
                                }
                            }
                            "nz" => {
                                has_normals = true;
                                if let ply::Property::Float(z) = list2 {
                                    nrm.z = z;
                                }
                            }
                            _ => {
                                unreachable!()
                            },
                        }
                    }
                }
            }
            _ => unreachable!(),
        }
    }
    // for i in 0..p.len() {
    //     println!("{:?}: {:?}", i, p[i]);
    // }
    // println!("tm_vertex_indices = {:?}", tm_vertex_indices);
    let mut n_ws: Vec<Normal3f> = Vec::new();
    if !n.is_empty() {
        assert!(n.len() == p.len());
        // transform normals to world space
        let n_normals: usize = n.len();
        for i in 0..n_normals {
            n_ws.push(o2w.transform_normal(&n[i]));
        }
    }
    // transform mesh vertices to world space
    let mut p_ws: Vec<Point3f> = Vec::new();
    let n_vertices: usize = p.len();
    for i in 0..n_vertices {
        p_ws.push(o2w.transform_point(&p[i]));
    }
    let s_ws: Vec<Vector3f> = Vec::new();
    // look up an alpha texture, if applicable
    let mut alpha_tex: Option<Arc<dyn Texture<Float> + Send + Sync>> = None;
    let alpha_tex_name: String = params.find_texture("alpha");
    if alpha_tex_name != String::from("") {
        alpha_tex = match float_textures.get(alpha_tex_name.as_str()) {
            Some(float_texture) => Some(float_texture.clone()),
            None => {
                println!(
                    "Couldn't find float texture {:?} for \"alpha\" parameter",
                    alpha_tex_name.as_str()
                );
                None
            }
        }
    } else {
        if params.find_one_float("alpha", 1.0 as Float) == 0.0 as Float {
            alpha_tex = Some(Arc::new(ConstantTexture::new(0.0 as Float)));
        }
    }
    let mut shadow_alpha_tex: Option<Arc<dyn Texture<Float> + Send + Sync>> = None;
    let shadow_alpha_tex_name: String = params.find_texture("shadowalpha");
    if shadow_alpha_tex_name != String::from("") {
        shadow_alpha_tex = match float_textures.get(shadow_alpha_tex_name.as_str()) {
            Some(float_texture) => Some(float_texture.clone()),
            None => {
                println!(
                    "Couldn't find float texture {:?} for \"shadowalpha\" parameter",
                    shadow_alpha_tex_name.as_str()
                );
                None
            }
        }
    } else {
        if params.find_one_float("shadowalpha", 1.0 as Float) == 0.0 as Float {
            shadow_alpha_tex = Some(Arc::new(ConstantTexture::new(0.0 as Float)));
        }
    }
    let mesh = Arc::new(TriangleMesh::new(
        *o2w,
        *w2o,
        reverse_orientation,
        (tm_vertex_indices.len() / 3).try_into().unwrap(), // n_triangles
        tm_vertex_indices,
        n_vertices.try_into().unwrap(),
        p_ws, // in world space
        s_ws, // in world space
        n_ws, // in world space
        uvs,
        alpha_tex,
        shadow_alpha_tex,
    ));
    let mut shapes: Vec<Arc<Shape>> = Vec::new();
    for id in 0..mesh.n_triangles {
        let triangle = Arc::new(Shape::Trngl(Triangle::new(
            mesh.object_to_world,
            mesh.world_to_object,
            mesh.transform_swaps_handedness,
            mesh.clone(),
            id.try_into().unwrap(),
        )));
        shapes.push(triangle.clone());
    }
    shapes
}
