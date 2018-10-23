// std
use std;
use std::collections::HashSet;
use std::hash::{Hash, Hasher};
use std::sync::Arc;
// pbrt
use core::geometry::{
    bnd3_expand, bnd3_union_bnd3, nrm_abs_dot_vec3, nrm_cross_vec3, nrm_dot_nrm, nrm_normalize,
    pnt3_distance, pnt3_distance_squared, pnt3_lerp, vec2_dot, vec3_coordinate_system,
    vec3_cross_vec3, vec3_normalize,
};
use core::geometry::{Bounds3f, Normal3f, Point2f, Point3f, Ray, Vector2f, Vector3f};
use core::interaction::{Interaction, InteractionCommon, SurfaceInteraction};
use core::material::Material;
use core::paramset::ParamSet;
use core::pbrt::Float;
use core::pbrt::{clamp_t, float_to_bits, lerp};
use core::shape::Shape;
use core::transform::Transform;

// see loopsubdiv.cpp

fn next(i: i32) -> i32 {
    (i + 1) % 3
}

fn prev(i: i32) -> i32 {
    (i + 2) % 3
}

#[derive(Debug, Default, Clone)]
struct SDVertex {
    p: Point3f,
    start_face: i32,
    child: i32,
    regular: bool,
    boundary: bool,
}

impl SDVertex {
    pub fn new(p: Point3f) -> Self {
        SDVertex {
            p: p,
            start_face: -1_i32,
            child: -1_i32,
            regular: false,
            boundary: false,
        }
    }
    pub fn valence(&self, vi: i32, faces: &Vec<Arc<SDFace>>) -> i32 {
        let mut fi: i32 = self.start_face;
        if !self.boundary {
            // compute valence of interior vertex
            let mut nf: i32 = 1;
            loop {
                if fi == -1_i32 {
                    break;
                }
                fi = faces[fi as usize].next_face(vi);
                if fi != self.start_face {
                    nf += 1;
                } else {
                    break;
                }
            }
            nf
        } else {
            // compute valence of boundary vertex
            let mut nf: i32 = 1;
            loop {
                if fi == -1_i32 {
                    break;
                }
                fi = faces[fi as usize].next_face(vi);
                if fi == -1_i32 {
                    break;
                } else {
                    nf += 1;
                }
            }
            fi = self.start_face;
            loop {
                if fi == -1_i32 {
                    break;
                }
                fi = faces[fi as usize].prev_face(vi);
                if fi == -1_i32 {
                    break;
                } else {
                    nf += 1;
                }
            }
            nf + 1
        }
    }
}

#[derive(Debug, Clone)]
struct SDFace {
    v: [i32; 3],
    f: [i32; 3],
    children: [i32; 4],
}

impl SDFace {
    pub fn next_face(&self, vi: i32) -> i32 {
        // see int vnum(SDVertex *vert) const {...}
        let mut fi: i32 = -1;
        for i in 0..3_usize {
            if self.v[i] == vi {
                fi = i as i32;
                break;
            }
        }
        if fi == -1_i32 {
            panic!("next_face({:?}, {})", self, vi);
        }
        self.f[fi as usize]
    }
    pub fn prev_face(&self, vi: i32) -> i32 {
        // see int vnum(SDVertex *vert) const {...}
        let mut fi: i32 = -1;
        for i in 0..3_usize {
            if self.v[i] == vi {
                fi = i as i32;
                break;
            }
        }
        if fi == -1_i32 {
            panic!("next_face({:?}, {})", self, vi);
        }
        self.f[prev(fi) as usize]
    }
}

impl Default for SDFace {
    fn default() -> SDFace {
        SDFace {
            v: [-1_i32, -1_i32, -1_i32],
            f: [-1_i32, -1_i32, -1_i32],
            children: [-1_i32, -1_i32, -1_i32, -1_i32],
        }
    }
}

#[derive(Debug, Default, Clone)]
struct SDEdge {
    v: [i32; 2],
    f: [i32; 2],
    f0_edge_num: i32,
}

impl SDEdge {
    pub fn new(v0: i32, v1: i32) -> Self {
        SDEdge {
            v: [std::cmp::min(v0, v1), std::cmp::max(v0, v1)],
            f: [-1_i32, -1_i32],
            f0_edge_num: -1,
        }
    }
}

impl PartialEq for SDEdge {
    fn eq(&self, other: &SDEdge) -> bool {
        if self.v[0] == other.v[0] {
            self.v[1] == other.v[1]
        } else {
            false
        }
    }
}

impl Eq for SDEdge {}

impl Hash for SDEdge {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.v[0].hash(state);
        self.v[1].hash(state);
    }
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
        let fi = i as i32; // face index
        for j in 0..3_usize {
            let vi: i32 = vertex_indices[i * 3 + j]; // vertex index
            if let Some(f) = Arc::get_mut(&mut faces[fi as usize]) {
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
    // set neighbor pointers in _faces_
    let mut edges: HashSet<SDEdge> = HashSet::new();
    for i in 0..n_faces {
        let fi = i as i32; // face index
        for edge_num in 0..3_usize {
            // update neighbor pointer for _edge_num_
            let v0: usize = edge_num;
            let v1: usize = next(edge_num as i32) as usize;
            let mut e: SDEdge = SDEdge::new(faces[i].v[v0], faces[i].v[v1]);
            if !edges.contains(&e) {
                // handle new edge
                e.f[0] = fi;
                e.f0_edge_num = edge_num as i32;
                edges.insert(e);
            } else {
                // handle previously seen edge
                let e_opt = edges.take(&e);
                if let Some(e) = e_opt {
                    // e.f[0]->f[e.f0_edge_num] = f;
                    if let Some(f) = Arc::get_mut(&mut faces[e.f[0] as usize]) {
                        f.f[e.f0_edge_num as usize] = fi;
                    } else {
                        panic!("Arc::get_mut(&mut faces[{}]) failed", e.f[0]);
                    }
                    // f->f[edge_num] = e.f[0];
                    if let Some(f) = Arc::get_mut(&mut faces[fi as usize]) {
                        f.f[edge_num] = e.f[0];
                    } else {
                        panic!("Arc::get_mut(&mut faces[{}]) failed", fi);
                    }
                }
            }
        }
    }
    // finish vertex initialization
    for vi in 0..p.len() {
        let start_face = verts[vi].start_face;
        let mut fi = start_face;
        loop {
            fi = faces[fi as usize].next_face(vi as i32);
            if fi == -1_i32 || fi == start_face {
                break;
            }
        }
        let valence = verts[vi].valence(vi as i32, &faces);
        if let Some(v) = Arc::get_mut(&mut verts[vi]) {
            v.boundary = (fi == -1_i32);
            if !v.boundary && valence == 6 {
                v.regular = true;
            } else if v.boundary && valence == 4 {
                v.regular = true;
            } else {
                v.regular = false;
            }
        }
    }
    // WORK
    Vec::new()
}
