// std
use std::collections::{HashMap, HashSet};
use std::convert::TryInto;
use std::f32::consts::PI;
use std::hash::{Hash, Hasher};
use std::sync::Arc;
// others
use smallvec::SmallVec;
// pbrt
use crate::core::geometry::vec3_cross_vec3;
use crate::core::geometry::{Normal3f, Point3f, Vector3f};
use crate::core::pbrt::Float;
use crate::core::transform::Transform;
use crate::shapes::triangle::TriangleMesh;

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
            p,
            start_face: -1_i32,
            child: -1_i32,
            regular: false,
            boundary: false,
        }
    }
    pub fn one_ring(
        &self,
        p: &mut SmallVec<[Point3f; 128]>,
        vi: i32,
        faces: &[Arc<SDFace>],
        verts: &[Arc<SDVertex>],
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
            let mut fi2: i32;
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
    pub fn valence(&self, vi: i32, faces: &[Arc<SDFace>]) -> i32 {
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
    pub fn vnum(&self, vi: i32) -> i32 {
        let mut fi: i32 = -1;
        for i in 0..3_usize {
            if self.v[i] == vi {
                fi = i as i32;
                break;
            }
        }
        fi
    }
    pub fn next_face(&self, vi: i32) -> i32 {
        let fi: i32 = self.vnum(vi);
        if fi == -1_i32 {
            panic!("next_face({:?}, {})", self, vi);
        }
        self.f[fi as usize]
    }
    pub fn prev_face(&self, vi: i32) -> i32 {
        let fi: i32 = self.vnum(vi);
        if fi == -1_i32 {
            panic!("prev_face({:?}, {})", self, vi);
        }
        self.f[prev(fi) as usize]
    }
    pub fn next_vert(&self, vi: i32) -> i32 {
        let fi: i32 = self.vnum(vi);
        if fi == -1_i32 {
            panic!("next_vert({:?}, {})", self, vi);
        }
        self.v[next(fi) as usize]
    }
    pub fn prev_vert(&self, vi: i32) -> i32 {
        let fi: i32 = self.vnum(vi);
        if fi == -1_i32 {
            panic!("prev_vert({:?}, {})", self, vi);
        }
        self.v[prev(fi) as usize]
    }
    pub fn other_vert(&self, vi0: i32, vi1: i32) -> i32 {
        let fi: i32 = -1;
        for i in 0..3_usize {
            if self.v[i] != vi0 && self.v[i] != vi1 {
                return self.v[i];
            }
        }
        if fi == -1_i32 {
            panic!("other_vert({:?}, {}, {})", self, vi0, vi1);
        }
        -1_i32
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

fn loop_gamma(valence: i32) -> Float {
    1.0 as Float / (valence as Float + 3.0 as Float / (8.0 as Float * beta(valence)))
}

pub fn loop_subdivide(
    object_to_world: &Transform,
    world_to_object: &Transform,
    reverse_orientation: bool,
    n_levels: i32,
    vertex_indices: &[i32],
    p: &[Point3f],
) -> Arc<TriangleMesh> {
    // allocate _LoopSubdiv_ vertices and faces
    let mut verts: Vec<Arc<SDVertex>> = Vec::with_capacity(p.len());
    for item in p {
        verts.push(Arc::new(SDVertex::new(*item)));
    }
    let n_faces: usize = vertex_indices.len() / 3;
    let mut faces: Vec<Arc<SDFace>> = Vec::with_capacity(n_faces);
    for _i in 0..n_faces {
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
    for (vi, mut vert) in verts.iter_mut().enumerate().take(p.len()) {
        let start_face = vert.start_face;
        let mut fi = start_face;
        loop {
            fi = faces[fi as usize].next_face(vi as i32);
            if fi == -1_i32 || fi == start_face {
                break;
            }
        }
        let valence = vert.valence(vi as i32, &faces);
        if let Some(v) = Arc::get_mut(&mut vert) {
            v.boundary = fi == -1_i32;
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
    for _i in 0..n_levels {
        // update _faces_ and _verts_ for next level of subdivision
        let mut new_faces: Vec<Arc<SDFace>> = Vec::new();
        let mut new_vertices: Vec<Arc<SDVertex>> = Vec::new();
        // allocate next level of children in mesh tree
        for mut vert in &mut verts {
            if let Some(vertex) = Arc::get_mut(&mut vert) {
                let ci = new_vertices.len();
                vertex.child = ci as i32;
                new_vertices.push(Arc::new(SDVertex::default()));
                if let Some(child) = Arc::get_mut(&mut new_vertices[ci]) {
                    child.regular = vertex.regular;
                    child.boundary = vertex.boundary;
                }
            }
        }
        for mut face in &mut faces {
            if let Some(face) = Arc::get_mut(&mut face) {
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
                    if verts[vi].regular {
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
        // compute new odd edge vertices
        let mut edge_verts: HashMap<SDEdge, i32> = HashMap::new();
        for fi in 0..faces.len() {
            for k in 0..3 {
                // compute odd vertex on _k_th edge
                let edge: SDEdge =
                    SDEdge::new(faces[fi].v[k as usize], faces[fi].v[next(k) as usize]);
                let contains_edge: bool = edge_verts.contains_key(&edge);
                if !contains_edge {
                    // create and initialize new odd vertex
                    let nvi = new_vertices.len();
                    new_vertices.push(Arc::new(SDVertex::default()));
                    if let Some(vert) = Arc::get_mut(&mut new_vertices[nvi]) {
                        vert.regular = true;
                        vert.boundary = faces[fi].f[k as usize] == -1_i32;
                        vert.start_face = faces[fi].children[3];
                        // apply edge rules to compute new vertex position
                        if vert.boundary {
                            vert.p = verts[edge.v[0] as usize].p * 0.5 as Float;
                            vert.p += verts[edge.v[1] as usize].p * 0.5 as Float;
                        } else {
                            vert.p = verts[edge.v[0] as usize].p * (3.0 as Float / 8.0 as Float);
                            vert.p += verts[edge.v[1] as usize].p * (3.0 as Float / 8.0 as Float);
                            let vi = faces[fi].other_vert(edge.v[0], edge.v[1]);
                            vert.p += verts[vi as usize].p * (1.0 as Float / 8.0 as Float);
                            let vi = faces[faces[fi].f[k as usize] as usize]
                                .other_vert(edge.v[0], edge.v[1]);
                            vert.p += verts[vi as usize].p * (1.0 as Float / 8.0 as Float);
                        }
                        edge_verts.insert(edge, nvi as i32);
                    }
                }
            }
        }
        // update even vertex face pointers
        for (vi, mut vert) in verts.iter_mut().enumerate() {
            let mut ci = -1_i32;
            let mut face_child = -1_i32;
            if let Some(vertex) = Arc::get_mut(&mut vert) {
                let start_face = vertex.start_face as usize;
                let face = faces[start_face].clone();
                let vert_num: usize = face.vnum(vi as i32) as usize;
                ci = vertex.child;
                face_child = face.children[vert_num];
            }
            if ci != -1_i32 {
                if let Some(child) = Arc::get_mut(&mut new_vertices[ci as usize]) {
                    child.start_face = face_child; // index into new_faces !!!
                }
            }
        }
        // update face neighbor pointers
        for fi in 0..faces.len() {
            let face = faces[fi].clone();
            for j in 0..3 {
                // update children _f_ pointers for siblings
                let ci = face.children[3] as usize;
                if let Some(child) = Arc::get_mut(&mut new_faces[ci]) {
                    child.f[j] = face.children[next(j as i32) as usize];
                }
                let ci = face.children[j] as usize;
                if let Some(child) = Arc::get_mut(&mut new_faces[ci]) {
                    child.f[next(j as i32) as usize] = face.children[3];
                }
                // update children _f_ pointers for neighbor children
                let fi2 = face.f[j];
                if fi2 != -1_i32 {
                    let f2 = faces[fi2 as usize].clone();
                    let ci2 = f2.children[f2.vnum(face.v[j]) as usize];
                    let ci = face.children[j] as usize;
                    if let Some(child) = Arc::get_mut(&mut new_faces[ci]) {
                        child.f[j] = ci2;
                    }
                } else {
                    let ci = face.children[j] as usize;
                    if let Some(child) = Arc::get_mut(&mut new_faces[ci]) {
                        child.f[j] = -1_i32;
                    }
                }
                let fi2 = face.f[prev(j as i32) as usize];
                if fi2 != -1_i32 {
                    let f2 = faces[fi2 as usize].clone();
                    let ci2 = f2.children[f2.vnum(face.v[j]) as usize];
                    let ci = face.children[j] as usize;
                    if let Some(child) = Arc::get_mut(&mut new_faces[ci]) {
                        child.f[prev(j as i32) as usize] = ci2;
                    }
                } else {
                    let ci = face.children[j] as usize;
                    if let Some(child) = Arc::get_mut(&mut new_faces[ci]) {
                        child.f[prev(j as i32) as usize] = -1_i32;
                    }
                }
            }
        }
        // update face vertex pointers
        for face in &faces {
            let mut nvi = -1_i32;
            for j in 0..3_usize {
                // update child vertex pointer to new even vertex
                let ci = face.children[j] as usize;
                if let Some(child) = Arc::get_mut(&mut new_faces[ci]) {
                    let vi = face.v[j] as usize;
                    let vertex = verts[vi].clone();
                    child.v[j] = vertex.child;
                    // update child vertex pointer to new odd vertex
                    let key = SDEdge::new(face.v[j], face.v[next(j as i32) as usize]);
                    let vert_opt = edge_verts.get(&key);
                    if let Some(vi2) = vert_opt {
                        nvi = *vi2;
                        child.v[next(j as i32) as usize] = nvi;
                    }
                }
                let ci = face.children[next(j as i32) as usize] as usize;
                if let Some(child) = Arc::get_mut(&mut new_faces[ci]) {
                    child.v[j] = nvi;
                }
                let ci = face.children[3_usize] as usize;
                if let Some(child) = Arc::get_mut(&mut new_faces[ci]) {
                    child.v[j] = nvi;
                }
            }
        }

        // prepare for next level of subdivision
        faces = new_faces.split_off(0);
        verts = new_vertices.split_off(0);
    }
    // push vertices to limit surface
    let mut p_limit: Vec<Point3f> = Vec::with_capacity(verts.len());
    for i in 0..verts.len() {
        let v = verts[i].clone();
        if v.boundary {
            p_limit.push(weight_boundary(
                v.clone(),
                1.0 as Float / 5.0 as Float,
                i as i32,
                &faces,
                &verts,
            ));
        } else {
            p_limit.push(weight_one_ring(
                v.clone(),
                loop_gamma(v.clone().valence(i as i32, &faces)),
                i as i32,
                &faces,
                &verts,
            ));
        }
    }
    for i in 0..verts.len() {
        if let Some(v) = Arc::get_mut(&mut verts[i]) {
            v.p = p_limit[i];
        }
    }
    // compute vertex tangents on limit surface
    let mut ns: Vec<Normal3f> = Vec::with_capacity(verts.len());
    let mut p_ring: Vec<Point3f> = Vec::new();
    p_ring.resize(16_usize, Point3f::default());
    for vi in 0..verts.len() {
        let vertex = verts[vi].clone();
        let mut s: Vector3f = Vector3f::default();
        let mut t: Vector3f = Vector3f::default();
        let valence: i32 = vertex.valence(vi as i32, &faces);
        if valence as usize > p_ring.len() {
            p_ring.resize(valence as usize, Point3f::default());
        }
        vertex.one_ring(
            &mut SmallVec::from_vec(p_ring.to_vec()),
            vi as i32,
            &faces,
            &verts,
        );
        if !vertex.boundary {
            // compute tangents of interior face
            for (j, item) in p_ring.iter().enumerate().take(valence as usize) {
                s += Vector3f::from(*item)
                    * (2.0 as Float * PI * j as Float / valence as Float).cos();
                t += Vector3f::from(*item)
                    * (2.0 as Float * PI * j as Float / valence as Float).sin();
            }
        } else {
            // compute tangents of boundary face
            s = p_ring[(valence - 1) as usize] - p_ring[0_usize];
            if valence == 2_i32 {
                t = p_ring[0] + p_ring[1] - vertex.p * 2.0 as Float;
            } else if valence == 3_i32 {
                t = p_ring[1] - vertex.p;
            } else if valence == 4_i32 {
                // regular
                t = Vector3f::from(
                    p_ring[0] * -1.0 as Float
                        + p_ring[1] * 2.0 as Float
                        + p_ring[2] * 2.0 as Float
                        + p_ring[3] * -1.0 as Float
                        + vertex.p * -2.0 as Float,
                );
            } else {
                let theta: Float = PI / (valence - 1) as Float;
                t = Vector3f::from((p_ring[0] + p_ring[(valence - 1) as usize]) * theta.sin());
                for (k, item) in p_ring
                    .iter()
                    .enumerate()
                    .take((valence - 1) as usize)
                    .skip(1)
                {
                    let wt: Float =
                        (2.0 as Float * theta.cos() - 2.0 as Float) * (k as Float * theta).sin();
                    t += Vector3f::from(*item * wt);
                }
                t = -t;
            }
        }
        ns.push(Normal3f::from(vec3_cross_vec3(&s, &t)));
    }
    // create triangle mesh from subdivision mesh
    let ntris: usize = faces.len();
    let mut vertex_indices: Vec<u32> = Vec::with_capacity(3 * ntris);
    let tot_verts: usize = verts.len();
    for face in faces.iter().take(ntris) {
        for j in 0..3_usize {
            vertex_indices.push(face.v[j] as u32);
        }
    }
    // transform mesh vertices to world space
    let mut p_ws: Vec<Point3f> = Vec::new();
    let n_vertices: usize = p_limit.len();
    for item in p_limit.iter().take(n_vertices) {
        p_ws.push(object_to_world.transform_point(item));
    }
    // transform normals to world space
    let mut n_ws: Vec<Normal3f> = Vec::new();
    let n_normals: usize = ns.len();
    for item in ns.iter().take(n_normals) {
        n_ws.push(object_to_world.transform_normal(item));
    }
    Arc::new(TriangleMesh::new(
        *object_to_world,
        *world_to_object,
        reverse_orientation,
        ntris.try_into().unwrap(),
        vertex_indices,
        tot_verts.try_into().unwrap(),
        p_ws, // in world space
        Vec::new(),
        n_ws, // in world space
        Vec::new(),
        None,
        None,
    ))
}

fn weight_one_ring(
    vert: Arc<SDVertex>,
    beta: Float,
    vi: i32,
    faces: &[Arc<SDFace>],
    verts: &[Arc<SDVertex>],
) -> Point3f {
    // put _vert_ one-ring in _p_ring_
    let valence: i32 = vert.valence(vi, faces);
    let mut p_ring: SmallVec<[Point3f; 128]> = SmallVec::with_capacity(valence as usize);
    for _i in 0..valence as usize {
        p_ring.push(Point3f::default());
    }
    vert.one_ring(&mut p_ring, vi, faces, verts);
    let mut p: Point3f = vert.p * (1.0 as Float - valence as Float * beta);
    for i in 0..valence as usize {
        p += p_ring[i] * beta;
    }
    p
}

fn weight_boundary(
    vert: Arc<SDVertex>,
    beta: Float,
    vi: i32,
    faces: &[Arc<SDFace>],
    verts: &[Arc<SDVertex>],
) -> Point3f {
    // put _vert_ one-ring in _p_ring_
    let valence: i32 = vert.valence(vi, faces);
    let mut p_ring: SmallVec<[Point3f; 128]> = SmallVec::with_capacity(valence as usize);
    for _i in 0..valence as usize {
        p_ring.push(Point3f::default());
    }
    vert.one_ring(&mut p_ring, vi, faces, verts);
    let mut p: Point3f = vert.p * (1.0 as Float - 2.0 as Float * beta);
    p += p_ring[0] * beta;
    p += p_ring[(valence - 1) as usize] * beta;
    p
}
