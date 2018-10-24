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

#[derive(Debug, Clone)]
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
    pub fn one_ring(
        &self,
        p: &mut Vec<Point3f>,
        vi: i32,
        faces: &Vec<Arc<SDFace>>,
        verts: &Vec<Arc<SDVertex>>,
    ) {
        if !self.boundary {
            // get one-ring vertices for interior vertex
            let mut fi: i32 = self.start_face;
            let mut pi = 0_usize;
            loop {
                let nvi = faces[fi as usize].next_vert(vi);
                p[pi] = verts[nvi as usize].p;
                pi += 1_usize;
                fi = faces[fi as usize].next_face(vi);
                if fi == self.start_face {
                    break;
                }
            }
        } else {
            // get one-ring vertices for boundary vertex
            let mut fi: i32 = self.start_face;
            let mut fi2: i32 = -1_i32;
            let mut pi = 0_usize;
            loop {
                fi2 = faces[fi as usize].next_face(vi);
                if fi2 == -1_i32 {
                    break;
                } else {
                    fi = fi2;
                }
            }
            let nvi = faces[fi as usize].next_vert(vi);
            p[pi] = verts[nvi as usize].p;
            pi += 1_usize;
            loop {
                let nvi = faces[fi as usize].prev_vert(vi);
                p[pi] = verts[nvi as usize].p;
                pi += 1_usize;
                fi = faces[fi as usize].prev_face(vi);
                if fi == -1_i32 {
                    break;
                }
            }
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

impl Default for SDVertex {
    fn default() -> SDVertex {
        SDVertex {
            p: Point3f::default(),
            start_face: -1_i32,
            child: -1_i32,
            regular: false,
            boundary: false,
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
    pub fn next_vert(&self, vi: i32) -> i32 {
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
        self.v[next(fi) as usize]
    }
    pub fn prev_vert(&self, vi: i32) -> i32 {
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
        self.v[prev(fi) as usize]
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

fn beta(valence: i32) -> Float {
    if valence == 3_i32 {
        3.0 as Float / 16.0 as Float
    } else {
        3.0 as Float / (8.0 as Float * valence as Float)
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
    // allocate _LoopSubdiv_ vertices and faces
    let mut verts: Vec<Arc<SDVertex>> = Vec::with_capacity(p.len());
    for i in 0..p.len() {
        verts.push(Arc::new(SDVertex::new(p[i])));
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
                    if let Some(f) = Arc::get_mut(&mut faces[e.f[0] as usize]) {
                        f.f[e.f0_edge_num as usize] = fi;
                    } else {
                        panic!("Arc::get_mut(&mut faces[{}]) failed", e.f[0]);
                    }
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
    // refine _LoopSubdiv_ into triangles
    for i in 0..n_levels {
        // update _faces_ and _verts_ for next level of subdivision
        let mut new_faces: Vec<Arc<SDFace>> = Vec::new();
        let mut new_vertices: Vec<Arc<SDVertex>> = Vec::new();
        // allocate next level of children in mesh tree
        for vi in 0..verts.len() {
            if let Some(vertex) = Arc::get_mut(&mut verts[vi]) {
                let ci = new_vertices.len();
                vertex.child = ci as i32;
                new_vertices.push(Arc::new(SDVertex::default()));
                if let Some(child) = Arc::get_mut(&mut new_vertices[ci]) {
                    child.regular = vertex.regular;
                    child.boundary = vertex.boundary;
                }
            }
        }
        for fi in 0..faces.len() {
            if let Some(face) = Arc::get_mut(&mut faces[fi]) {
                for k in 0..4 {
                    let ci = new_faces.len();
                    new_faces.push(Arc::new(SDFace::default()));
                    face.children[k] = ci as i32;
                }
            }
        }
        // update vertex positions for even vertices
        for vi in 0..verts.len() {
            let ci = verts[vi].child as usize;
            if let Some(child) = Arc::get_mut(&mut new_vertices[ci]) {
                if !verts[vi].boundary {
                    // apply one-ring rule for even vertex
                    if (verts[vi].regular) {
                        child.p = weight_one_ring(
                            verts[vi].clone(),
                            1.0 as Float / 16.0 as Float,
                            vi as i32,
                            &faces,
                            &verts,
                        );
                    } else {
                        child.p = weight_one_ring(
                            verts[vi].clone(),
                            beta(verts[vi].valence(vi as i32, &faces)),
                            vi as i32,
                            &faces,
                            &verts,
                        );
                    }
                } else {
                    // apply boundary rule for even vertex
                    child.p = weight_boundary(
                        verts[vi].clone(),
                        1.0 as Float / 8.0 as Float,
                        vi as i32,
                        &faces,
                        &verts,
                    );
                }
            }
        }

        //     // Compute new odd edge vertices
        //     std::map<SDEdge, SDVertex *> edgeVerts;
        //     for (SDFace *face : f) {
        //         for (int k = 0; k < 3; ++k) {
        //             // Compute odd vertex on _k_th edge
        //             SDEdge edge(face->v[k], face->v[NEXT(k)]);
        //             SDVertex *vert = edgeVerts[edge];
        //             if (!vert) {
        //                 // Create and initialize new odd vertex
        //                 vert = arena.Alloc<SDVertex>();
        //                 new_vertices.push_back(vert);
        //                 vert->regular = true;
        //                 vert->boundary = (face->f[k] == nullptr);
        //                 vert->startFace = face->children[3];

        //                 // Apply edge rules to compute new vertex position
        //                 if (vert->boundary) {
        //                     vert->p = 0.5f * edge.v[0]->p;
        //                     vert->p += 0.5f * edge.v[1]->p;
        //                 } else {
        //                     vert->p = 3.f / 8.f * edge.v[0]->p;
        //                     vert->p += 3.f / 8.f * edge.v[1]->p;
        //                     vert->p += 1.f / 8.f *
        //                                face->otherVert(edge.v[0], edge.v[1])->p;
        //                     vert->p +=
        //                         1.f / 8.f *
        //                         face->f[k]->otherVert(edge.v[0], edge.v[1])->p;
        //                 }
        //                 edgeVerts[edge] = vert;
        //             }
        //         }
        //     }

        //     // Update new mesh topology

        //     // Update even vertex face pointers
        //     for (SDVertex *vertex : v) {
        //         int vertNum = vertex->startFace->vnum(vertex);
        //         vertex->child->startFace = vertex->startFace->children[vertNum];
        //     }

        //     // Update face neighbor pointers
        //     for (SDFace *face : f) {
        //         for (int j = 0; j < 3; ++j) {
        //             // Update children _f_ pointers for siblings
        //             face->children[3]->f[j] = face->children[NEXT(j)];
        //             face->children[j]->f[NEXT(j)] = face->children[3];

        //             // Update children _f_ pointers for neighbor children
        //             SDFace *f2 = face->f[j];
        //             face->children[j]->f[j] =
        //                 f2 ? f2->children[f2->vnum(face->v[j])] : nullptr;
        //             f2 = face->f[PREV(j)];
        //             face->children[j]->f[PREV(j)] =
        //                 f2 ? f2->children[f2->vnum(face->v[j])] : nullptr;
        //         }
        //     }

        //     // Update face vertex pointers
        //     for (SDFace *face : f) {
        //         for (int j = 0; j < 3; ++j) {
        //             // Update child vertex pointer to new even vertex
        //             face->children[j]->v[j] = face->v[j]->child;

        //             // Update child vertex pointer to new odd vertex
        //             SDVertex *vert =
        //                 edgeVerts[SDEdge(face->v[j], face->v[NEXT(j)])];
        //             face->children[j]->v[NEXT(j)] = vert;
        //             face->children[NEXT(j)]->v[j] = vert;
        //             face->children[3]->v[j] = vert;
        //         }
        //     }

        // prepare for next level of subdivision
        faces = new_faces.split_off(0);
        verts = new_vertices.split_off(0);
    }

    // // Push vertices to limit surface
    // std::unique_ptr<Point3f[]> pLimit(new Point3f[v.size()]);
    // for (size_t i = 0; i < v.size(); ++i) {
    //     if (v[i]->boundary)
    //         pLimit[i] = weight_boundary(v[i], 1.f / 5.f);
    //     else
    //         pLimit[i] = weight_one_ring(v[i], loopGamma(v[i]->valence()));
    // }
    // for (size_t i = 0; i < v.size(); ++i) v[i]->p = pLimit[i];

    // // Compute vertex tangents on limit surface
    // std::vector<Normal3f> Ns;
    // Ns.reserve(v.size());
    // std::vector<Point3f> pRing(16, Point3f());
    // for (SDVertex *vertex : v) {
    //     Vector3f S(0, 0, 0), T(0, 0, 0);
    //     int valence = vertex->valence();
    //     if (valence > (int)pRing.size()) pRing.resize(valence);
    //     vertex->oneRing(&pRing[0]);
    //     if (!vertex->boundary) {
    //         // Compute tangents of interior face
    //         for (int j = 0; j < valence; ++j) {
    //             S += std::cos(2 * Pi * j / valence) * Vector3f(pRing[j]);
    //             T += std::sin(2 * Pi * j / valence) * Vector3f(pRing[j]);
    //         }
    //     } else {
    //         // Compute tangents of boundary face
    //         S = pRing[valence - 1] - pRing[0];
    //         if (valence == 2)
    //             T = Vector3f(pRing[0] + pRing[1] - 2 * vertex->p);
    //         else if (valence == 3)
    //             T = pRing[1] - vertex->p;
    //         else if (valence == 4)  // regular
    //             T = Vector3f(-1 * pRing[0] + 2 * pRing[1] + 2 * pRing[2] +
    //                          -1 * pRing[3] + -2 * vertex->p);
    //         else {
    //             Float theta = Pi / float(valence - 1);
    //             T = Vector3f(std::sin(theta) * (pRing[0] + pRing[valence - 1]));
    //             for (int k = 1; k < valence - 1; ++k) {
    //                 Float wt = (2 * std::cos(theta) - 2) * std::sin((k)*theta);
    //                 T += Vector3f(wt * pRing[k]);
    //             }
    //             T = -T;
    //         }
    //     }
    //     Ns.push_back(Normal3f(Cross(S, T)));
    // }

    // // Create triangle mesh from subdivision mesh
    // {
    //     size_t ntris = f.size();
    //     std::unique_ptr<int[]> verts(new int[3 * ntris]);
    //     int *vp = verts.get();
    //     size_t totVerts = v.size();
    //     std::map<SDVertex *, int> usedVerts;
    //     for (size_t i = 0; i < totVerts; ++i) usedVerts[v[i]] = i;
    //     for (size_t i = 0; i < ntris; ++i) {
    //         for (int j = 0; j < 3; ++j) {
    //             *vp = usedVerts[f[i]->v[j]];
    //             ++vp;
    //         }
    //     }
    //     return CreateTriangleMesh(ObjectToWorld, WorldToObject,
    //                               reverseOrientation, ntris, verts.get(),
    //                               totVerts, pLimit.get(), nullptr, &Ns[0],
    //                               nullptr, nullptr, nullptr);
    // }
    // WORK
    Vec::new()
}

fn weight_one_ring(
    vert: Arc<SDVertex>,
    beta: Float,
    vi: i32,
    faces: &Vec<Arc<SDFace>>,
    verts: &Vec<Arc<SDVertex>>,
) -> Point3f {
    // put _vert_ one-ring in _p_ring_
    let valence: i32 = vert.valence(vi, faces);
    // Point3f *p_ring = ALLOCA(Point3f, valence);
    let mut p_ring: Vec<Point3f> = Vec::with_capacity(valence as usize);
    for _i in 0..valence as usize {
        p_ring.push(Point3f::default());
    }
    vert.one_ring(&mut p_ring, vi, faces, verts);
    // Point3f p = (1 - valence * beta) * vert->p;
    // for (int i = 0; i < valence; ++i) p += beta * p_ring[i];
    // return p;
    // WORK
    Point3f::default()
}

fn weight_boundary(
    vert: Arc<SDVertex>,
    beta: Float,
    vi: i32,
    faces: &Vec<Arc<SDFace>>,
    verts: &Vec<Arc<SDVertex>>,
) -> Point3f {
    // put _vert_ one-ring in _p_ring_
    let valence: i32 = vert.valence(vi, faces);
    // Point3f *p_ring = ALLOCA(Point3f, valence);
    let mut p_ring: Vec<Point3f> = Vec::with_capacity(valence as usize);
    for _i in 0..valence as usize {
        p_ring.push(Point3f::default());
    }
    vert.one_ring(&mut p_ring, vi, faces, verts);
    // Point3f p = (1 - 2 * beta) * vert->p;
    // p += beta * p_ring[0];
    // p += beta * p_ring[valence - 1];
    // return p;
    // WORK
    Point3f::default()
}
