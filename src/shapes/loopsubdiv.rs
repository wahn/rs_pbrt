// std
use std::sync::Arc;
// pbrt
use core::geometry::{Bounds3f, Normal3f, Point2f, Point3f, Ray, Vector2f, Vector3f};
use core::geometry::{nrm_dot_nrm, nrm_normalize, bnd3_expand, bnd3_union_bnd3, nrm_abs_dot_vec3,
                     nrm_cross_vec3, pnt3_distance, pnt3_distance_squared, pnt3_lerp, vec2_dot,
                     vec3_coordinate_system, vec3_cross_vec3, vec3_normalize};
use core::interaction::{Interaction, InteractionCommon, SurfaceInteraction};
use core::material::Material;
use core::paramset::ParamSet;
use core::pbrt::Float;
use core::pbrt::{clamp_t, float_to_bits, lerp};
use core::shape::Shape;
use core::transform::Transform;

// see loopsubdiv.cpp

#[derive(Debug, Default, Clone)]
struct SDVertex {
    p: Point3f,
    start_face: usize,
}

impl SDVertex {
    pub fn new(p: Point3f) -> Self {
        SDVertex {
            p: p,
            start_face: 0_usize,
        }
    }
}

#[derive(Debug, Default, Clone)]
struct SDFace {
    v: [usize; 3],
    f: [usize; 3],
    children: [usize; 4]
}

pub fn loop_subdivide(
    object_to_world: &Transform,
    world_to_object: &Transform,
    reverse_orientation: bool,
    n_levels: i32,
    vertex_indices: &Vec<i32>,
    p: &Vec<Point3f>,
) -> Vec<Arc<Shape + Send + Sync>> {
    //let mut vertices: Vec<Arc<SDVertex>> = Vec::with_capacity(p.len());
    // allocate _LoopSubdiv_ vertices and faces
    let mut verts: Vec<Arc<SDVertex>> = Vec::with_capacity(p.len());
    for i in 0..p.len() {
        verts.push(Arc::new(SDVertex::new(p[i])));
        //vertices.push(verts[i].clone());
    }
    let n_faces: usize = vertex_indices.len() / 3;
    let mut faces: Vec<Arc<SDFace>> = Vec::with_capacity(n_faces);
    for i in 0..n_faces {
        faces.push(Arc::new(SDFace::default()));
    }
    // set face to vertex pointers
    for i in 0..n_faces {
        let fi = i; // face index
        for j in 0..3_usize {
            let vi: usize = vertex_indices[i * 3 + j] as usize; // vertex index
            if let Some(f) = Arc::get_mut(&mut faces[fi]) {
                f.v[j] = vi;
            } else {
                panic!("Arc::get_mut(&mut faces[{}]) failed", fi);
            }
            if let Some(v) = Arc::get_mut(&mut verts[vi as usize]) {
                v.start_face = fi;
            } else {
                panic!("Arc::get_mut(&mut verts[{}] failed", vi as usize);
            }
        }
    }
    // WORK
    Vec::new()
}
