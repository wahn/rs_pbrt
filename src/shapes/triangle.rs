// std
use std::mem;
use std::sync::Arc;
// pbrt
use core::geometry::{Bounds3f, Normal3, Normal3f, Point2f, Point3f, Ray, Vector2f, Vector3f};
use core::geometry::{nrm_faceforward_nrm, nrm_normalize, bnd3_union_pnt3, nrm_abs_dot_vec3,
                     pnt3_abs, pnt3_distance_squared, pnt3_permute, vec3_coordinate_system,
                     vec3_cross_nrm, vec3_cross_vec3, vec3_max_component, vec3_max_dimension,
                     vec3_normalize, vec3_permute};
use core::interaction::{Interaction, InteractionCommon, SurfaceInteraction};
use core::material::Material;
use core::pbrt::Float;
use core::pbrt::gamma;
use core::sampling::uniform_sample_triangle;
use core::shape::Shape;
use core::transform::Transform;

// see triangle.h

#[derive(Debug, Clone)]
pub struct TriangleMesh {
    /// the total number of triangles in the mesh
    pub n_triangles: usize,
    /// vector of vertex indices
    pub vertex_indices: Vec<usize>,
    /// the total number of vertices in the mesh
    pub n_vertices: usize,
    /// vector of *n_vertices* vertex positions
    pub p: Vec<Point3f>,
    /// an optional vector of normal vectors (can be empty)
    pub n: Vec<Normal3f>,
    /// an optional vector of tangent vectors (can be empty)
    pub s: Vec<Vector3f>,
    /// an optional vector of paramtric (u, v) values (texture coordinates)
    pub uv: Vec<Point2f>,
    // TODO: std::shared_ptr<Texture<Float>> alphaMask, shadowAlphaMask;
    // inherited from class Shape (see shape.h)
    pub object_to_world: Transform, // TODO: not pub?
    pub world_to_object: Transform, // TODO: not pub?
    pub reverse_orientation: bool,
    pub transform_swaps_handedness: bool, // TODO: not pub?
}

impl TriangleMesh {
    pub fn new(
        object_to_world: Transform,
        world_to_object: Transform,
        reverse_orientation: bool,
        transform_swaps_handedness: bool,
        n_triangles: usize,
        vertex_indices: Vec<usize>,
        n_vertices: usize,
        p: Vec<Point3f>,
        s: Vec<Vector3f>,
        n: Vec<Normal3f>,
        uv: Vec<Point2f>,
    ) -> Self {
        TriangleMesh {
            // Shape
            object_to_world: object_to_world,
            world_to_object: world_to_object,
            reverse_orientation: reverse_orientation,
            transform_swaps_handedness: transform_swaps_handedness,
            // TriangleMesh
            n_triangles: n_triangles,
            vertex_indices: vertex_indices,
            n_vertices: n_vertices,
            p: p,
            n: n,
            s: s,
            uv: uv,
        }
    }
}

#[derive(Clone)]
pub struct Triangle {
    mesh: Arc<TriangleMesh>,
    pub id: usize,
    // inherited from class Shape (see shape.h)
    object_to_world: Transform,
    world_to_object: Transform,
    reverse_orientation: bool,
    transform_swaps_handedness: bool,
    pub material: Option<Arc<Material + Send + Sync>>,
}

impl Triangle {
    pub fn new(
        object_to_world: Transform,
        world_to_object: Transform,
        reverse_orientation: bool,
        mesh: Arc<TriangleMesh>,
        tri_number: usize,
    ) -> Self {
        Triangle {
            mesh: mesh,
            id: tri_number,
            object_to_world: object_to_world,
            world_to_object: world_to_object,
            reverse_orientation: reverse_orientation,
            transform_swaps_handedness: false,
            material: None,
        }
    }
    pub fn get_uvs(&self) -> [Point2f; 3] {
        if self.mesh.uv.is_empty() {
            [
                Point2f { x: 0.0, y: 0.0 },
                Point2f { x: 1.0, y: 0.0 },
                Point2f { x: 1.0, y: 1.0 },
            ]
        } else {
            [
                self.mesh.uv[self.mesh.vertex_indices[self.id * 3 + 0]],
                self.mesh.uv[self.mesh.vertex_indices[self.id * 3 + 1]],
                self.mesh.uv[self.mesh.vertex_indices[self.id * 3 + 2]],
            ]
        }
    }
}

impl Shape for Triangle {
    fn object_bound(&self) -> Bounds3f {
        let p0: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 0]];
        let p1: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 1]];
        let p2: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 2]];
        bnd3_union_pnt3(
            &Bounds3f::new(
                self.world_to_object.transform_point(&p0),
                self.world_to_object.transform_point(&p1),
            ),
            &self.world_to_object.transform_point(&p2),
        )
    }
    fn world_bound(&self) -> Bounds3f {
        let p0: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 0]];
        let p1: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 1]];
        let p2: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 2]];
        bnd3_union_pnt3(&Bounds3f::new(p0, p1), &p2)
    }
    fn intersect(&self, ray: &Ray) -> Option<(SurfaceInteraction, Float)> {
        // get triangle vertices in _p0_, _p1_, and _p2_
        let p0: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 0]];
        let p1: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 1]];
        let p2: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 2]];
        // translate vertices based on ray origin
        let mut p0t: Point3f = p0 - Vector3f {
            x: ray.o.x,
            y: ray.o.y,
            z: ray.o.z,
        };
        let mut p1t: Point3f = p1 - Vector3f {
            x: ray.o.x,
            y: ray.o.y,
            z: ray.o.z,
        };
        let mut p2t: Point3f = p2 - Vector3f {
            x: ray.o.x,
            y: ray.o.y,
            z: ray.o.z,
        };
        // permute components of triangle vertices and ray direction
        let kz: usize = vec3_max_dimension(&ray.d.abs());
        let mut kx: usize = kz + 1;
        if kx == 3 {
            kx = 0;
        }
        let mut ky: usize = kx + 1;
        if ky == 3 {
            ky = 0;
        }
        let d: Vector3f = vec3_permute(&ray.d, kx, ky, kz);
        p0t = pnt3_permute(&p0t, kx, ky, kz);
        p1t = pnt3_permute(&p1t, kx, ky, kz);
        p2t = pnt3_permute(&p2t, kx, ky, kz);
        // apply shear transformation to translated vertex positions
        let sx: Float = -d.x / d.z;
        let sy: Float = -d.y / d.z;
        let sz: Float = 1.0 / d.z;
        p0t.x += sx * p0t.z;
        p0t.y += sy * p0t.z;
        p1t.x += sx * p1t.z;
        p1t.y += sy * p1t.z;
        p2t.x += sx * p2t.z;
        p2t.y += sy * p2t.z;
        // compute edge function coefficients _e0_, _e1_, and _e2_
        let mut e0: Float = p1t.x * p2t.y - p1t.y * p2t.x;
        let mut e1: Float = p2t.x * p0t.y - p2t.y * p0t.x;
        let mut e2: Float = p0t.x * p1t.y - p0t.y * p1t.x;
        // fall back to double precision test at triangle edges
        if mem::size_of::<Float>() == mem::size_of::<f32>() && (e0 == 0.0 || e1 == 0.0 || e2 == 0.0)
        {
            let p2txp1ty: f64 = p2t.x as f64 * p1t.y as f64;
            let p2typ1tx: f64 = p2t.y as f64 * p1t.x as f64;
            e0 = (p2typ1tx - p2txp1ty) as Float;
            let p0txp2ty = p0t.x as f64 * p2t.y as f64;
            let p0typ2tx = p0t.y as f64 * p2t.x as f64;
            e1 = (p0typ2tx - p0txp2ty) as Float;
            let p1txp0ty = p1t.x as f64 * p0t.y as f64;
            let p1typ0tx = p1t.y as f64 * p0t.x as f64;
            e2 = (p1typ0tx - p1txp0ty) as Float;
        }
        // perform triangle edge and determinant tests
        if (e0 < 0.0 || e1 < 0.0 || e2 < 0.0) && (e0 > 0.0 || e1 > 0.0 || e2 > 0.0) {
            return None;
        }
        let det: Float = e0 + e1 + e2;
        if det == 0.0 {
            return None;
        }
        // compute scaled hit distance to triangle and test against ray $t$ range
        p0t.z *= sz;
        p1t.z *= sz;
        p2t.z *= sz;
        let t_scaled: Float = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
        if det < 0.0 && (t_scaled >= 0.0 || t_scaled < ray.t_max * det) {
            return None;
        } else if det > 0.0 && (t_scaled <= 0.0 || t_scaled > ray.t_max * det) {
            return None;
        }
        // compute barycentric coordinates and $t$ value for triangle intersection
        let inv_det: Float = 1.0 / det;
        let b0: Float = e0 * inv_det;
        let b1: Float = e1 * inv_det;
        let b2: Float = e2 * inv_det;
        let t: Float = t_scaled * inv_det;

        // ensure that computed triangle $t$ is conservatively greater than zero

        // compute $\delta_z$ term for triangle $t$ error bounds
        let max_zt: Float = vec3_max_component(&Vector3f {
            x: p0t.z,
            y: p1t.z,
            z: p2t.z,
        }.abs());
        let delta_z: Float = gamma(3_i32) * max_zt;
        // compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
        let max_xt: Float = vec3_max_component(&Vector3f {
            x: p0t.x,
            y: p1t.x,
            z: p2t.x,
        }.abs());
        let max_yt: Float = vec3_max_component(&Vector3f {
            x: p0t.y,
            y: p1t.y,
            z: p2t.y,
        }.abs());
        let delta_x: Float = gamma(5) * (max_xt + max_zt);
        let delta_y: Float = gamma(5) * (max_yt + max_zt);
        // compute $\delta_e$ term for triangle $t$ error bounds
        let delta_e: Float =
            2.0 * (gamma(2) * max_xt * max_yt + delta_y * max_xt + delta_x * max_yt);
        // compute $\delta_t$ term for triangle $t$ error bounds and check _t_
        let max_e: Float = vec3_max_component(&Vector3f {
            x: e0,
            y: e1,
            z: e2,
        }.abs());
        let delta_t: Float =
            3.0 * (gamma(3) * max_e * max_zt + delta_e * max_zt + delta_z * max_e) * inv_det.abs();
        if t <= delta_t {
            return None;
        }
        // compute triangle partial derivatives
        let uv: [Point2f; 3] = self.get_uvs();
        // compute deltas for triangle partial derivatives
        let duv02: Vector2f = uv[0] - uv[2];
        let duv12: Vector2f = uv[1] - uv[2];
        let dp02: Vector3f = p0 - p2;
        let dp12: Vector3f = p1 - p2;
        let determinant: Float = duv02.x * duv12.y - duv02.y * duv12.x;
        let degenerate_uv: bool = determinant.abs() < 1e-8 as Float;
        // Vector3f dpdu, dpdv;
        let mut dpdu: Vector3f = Vector3f::default();
        let mut dpdv: Vector3f = Vector3f::default();
        if !degenerate_uv {
            let invdet: Float = 1.0 / determinant;
            dpdu = (dp02 * duv12.y - dp12 * duv02.y) * invdet;
            dpdv = (dp02 * -duv12.x + dp12 * duv02.x) * invdet;
        }
        if degenerate_uv || vec3_cross_vec3(&dpdu, &dpdv).length_squared() == 0.0 {
            // handle zero determinant for triangle partial derivative matrix
            vec3_coordinate_system(
                &vec3_normalize(&vec3_cross_vec3(&(p2 - p0), &(p1 - p0))),
                &mut dpdu,
                &mut dpdv,
            );
        }
        // compute error bounds for triangle intersection
        let x_abs_sum: Float = (b0 * p0.x).abs() + (b1 * p1.x).abs() + (b2 * p2.x).abs();
        let y_abs_sum: Float = (b0 * p0.y).abs() + (b1 * p1.y).abs() + (b2 * p2.y).abs();
        let z_abs_sum: Float = (b0 * p0.z).abs() + (b1 * p1.z).abs() + (b2 * p2.z).abs();
        let p_error: Vector3f = Vector3f {
            x: x_abs_sum,
            y: y_abs_sum,
            z: z_abs_sum,
        } * gamma(7);
        // interpolate $(u,v)$ parametric coordinates and hit point
        let p_hit: Point3f = p0 * b0 + p1 * b1 + p2 * b2;
        let uv_hit: Point2f = uv[0] * b0 + uv[1] * b1 + uv[2] * b2;
        // TODO: test intersection against alpha texture, if present
        // if (testAlphaTexture && mesh->alphaMask) {
        //     SurfaceInteraction isectLocal(p_hit, Vector3f(0, 0, 0), uv_hit, -ray.d,
        //                                   dpdu, dpdv, Normal3f(0, 0, 0),
        //                                   Normal3f(0, 0, 0), ray.time, this);
        //     if (mesh->alphaMask->Evaluate(isectLocal) == 0) return false;
        // }
        // fill in _SurfaceInteraction_ from triangle hit
        let dndu: Normal3f = Normal3f::default();
        let dndv: Normal3f = Normal3f::default();
        let wo: Vector3f = -ray.d;
        let mut si: SurfaceInteraction = SurfaceInteraction::new(
            &p_hit,
            &p_error,
            &uv_hit,
            &wo,
            &dpdu,
            &dpdv,
            &dndu,
            &dndv,
            ray.time,
            Some(self),
        );
        // override surface normal in _isect_ for triangle
        let surface_normal: Normal3f =
            Normal3f::from(vec3_normalize(&vec3_cross_vec3(&dp02, &dp12)));
        si.n = surface_normal;
        si.shading.n = surface_normal;
        if !self.mesh.n.is_empty() || !self.mesh.s.is_empty() {
            // initialize _Triangle_ shading geometry

            // compute shading normal _ns_ for triangle
            let mut ns: Normal3f;
            if !self.mesh.n.is_empty() {
                let n0 = self.mesh.n[self.mesh.vertex_indices[self.id * 3 + 0]];
                let n1 = self.mesh.n[self.mesh.vertex_indices[self.id * 3 + 1]];
                let n2 = self.mesh.n[self.mesh.vertex_indices[self.id * 3 + 2]];
                ns = Normal3::from(n0) * b0 + Normal3::from(n1) * b1 + Normal3::from(n2) * b2;
                if ns.length_squared() > 0.0 {
                    ns = nrm_normalize(&ns);
                } else {
                    ns = si.n;
                }
            } else {
                ns = si.n;
            }
            // compute shading tangent _ss_ for triangle
            let mut ss: Vector3f;
            if !self.mesh.s.is_empty() {
                let s0 = self.mesh.s[self.mesh.vertex_indices[self.id * 3 + 0]];
                let s1 = self.mesh.s[self.mesh.vertex_indices[self.id * 3 + 1]];
                let s2 = self.mesh.s[self.mesh.vertex_indices[self.id * 3 + 2]];
                ss = s0 * b0 + s1 * b1 + s2 * b2;
                if ss.length_squared() > 0.0 {
                    ss = vec3_normalize(&ss);
                } else {
                    ss = vec3_normalize(&si.dpdu);
                }
            } else {
                ss = vec3_normalize(&si.dpdu);
            }
            // compute shading bitangent _ts_ for triangle and adjust _ss_
            let mut ts: Vector3f = vec3_cross_nrm(&ss, &ns);
            if ts.length_squared() > 0.0 {
                ts = vec3_normalize(&ts);
                ss = vec3_cross_nrm(&ts, &ns);
            } else {
                vec3_coordinate_system(&Vector3f::from(ns), &mut ss, &mut ts);
            }
            // compute $\dndu$ and $\dndv$ for triangle shading geometry
            let dndu: Normal3f;
            let dndv: Normal3f;
            if !self.mesh.n.is_empty() {
                // compute deltas for triangle partial derivatives of normal
                let duv02: Vector2f = uv[0] - uv[2];
                let duv12: Vector2f = uv[1] - uv[2];
                let dn1: Normal3f =
                    Normal3::from(self.mesh.n[self.mesh.vertex_indices[self.id * 3 + 0]])
                        - Normal3::from(self.mesh.n[self.mesh.vertex_indices[self.id * 3 + 2]]);
                let dn2: Normal3f =
                    Normal3::from(self.mesh.n[self.mesh.vertex_indices[self.id * 3 + 1]])
                        - Normal3::from(self.mesh.n[self.mesh.vertex_indices[self.id * 3 + 2]]);
                let determinant: Float = duv02.x * duv12.y - duv02.y * duv12.x;
                let degenerate_uv: bool = determinant.abs() < 1e-8;
                if degenerate_uv {
                    dndu = Normal3f::default();
                    dndv = Normal3f::default();
                } else {
                    let inv_det: Float = 1.0 / determinant;
                    dndu = (dn1 * duv12.y - dn2 * duv02.y) * inv_det;
                    dndv = (dn1 * -duv12.x + dn2 * duv02.x) * inv_det;
                }
            } else {
                dndu = Normal3f::default();
                dndv = Normal3f::default();
            }
            si.set_shading_geometry(&ss, &ts, &dndu, &dndv, true);
        }
        // ensure correct orientation of the geometric normal
        if !self.mesh.n.is_empty() {
            si.n = nrm_faceforward_nrm(&si.n, &si.shading.n);
        } else if self.reverse_orientation ^ self.transform_swaps_handedness {
            si.shading.n = -si.n;
            si.n = -si.n;
        }
        Some((si, t as Float))
    }
    fn intersect_p(&self, ray: &Ray) -> bool {
        // TODO: ProfilePhase p(Prof::TriIntersectP);
        // TODO: ++nTests;
        // get triangle vertices in _p0_, _p1_, and _p2_
        let p0: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 0]];
        let p1: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 1]];
        let p2: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 2]];
        // translate vertices based on ray origin
        let mut p0t: Point3f = p0 - Vector3f {
            x: ray.o.x,
            y: ray.o.y,
            z: ray.o.z,
        };
        let mut p1t: Point3f = p1 - Vector3f {
            x: ray.o.x,
            y: ray.o.y,
            z: ray.o.z,
        };
        let mut p2t: Point3f = p2 - Vector3f {
            x: ray.o.x,
            y: ray.o.y,
            z: ray.o.z,
        };
        // permute components of triangle vertices and ray direction
        let kz: usize = vec3_max_dimension(&ray.d.abs());
        let mut kx: usize = kz + 1;
        if kx == 3 {
            kx = 0;
        }
        let mut ky: usize = kx + 1;
        if ky == 3 {
            ky = 0;
        }
        let d: Vector3f = vec3_permute(&ray.d, kx, ky, kz);
        p0t = pnt3_permute(&p0t, kx, ky, kz);
        p1t = pnt3_permute(&p1t, kx, ky, kz);
        p2t = pnt3_permute(&p2t, kx, ky, kz);
        // apply shear transformation to translated vertex positions
        let sx: Float = -d.x / d.z;
        let sy: Float = -d.y / d.z;
        let sz: Float = 1.0 / d.z;
        p0t.x += sx * p0t.z;
        p0t.y += sy * p0t.z;
        p1t.x += sx * p1t.z;
        p1t.y += sy * p1t.z;
        p2t.x += sx * p2t.z;
        p2t.y += sy * p2t.z;
        // compute edge function coefficients _e0_, _e1_, and _e2_
        let mut e0: Float = p1t.x * p2t.y - p1t.y * p2t.x;
        let mut e1: Float = p2t.x * p0t.y - p2t.y * p0t.x;
        let mut e2: Float = p0t.x * p1t.y - p0t.y * p1t.x;
        // fall back to double precision test at triangle edges
        if mem::size_of::<Float>() == mem::size_of::<f32>() && (e0 == 0.0 || e1 == 0.0 || e2 == 0.0)
        {
            let p2txp1ty: f64 = p2t.x as f64 * p1t.y as f64;
            let p2typ1tx: f64 = p2t.y as f64 * p1t.x as f64;
            e0 = (p2typ1tx - p2txp1ty) as Float;
            let p0txp2ty = p0t.x as f64 * p2t.y as f64;
            let p0typ2tx = p0t.y as f64 * p2t.x as f64;
            e1 = (p0typ2tx - p0txp2ty) as Float;
            let p1txp0ty = p1t.x as f64 * p0t.y as f64;
            let p1typ0tx = p1t.y as f64 * p0t.x as f64;
            e2 = (p1typ0tx - p1txp0ty) as Float;
        }
        // perform triangle edge and determinant tests
        if (e0 < 0.0 || e1 < 0.0 || e2 < 0.0) && (e0 > 0.0 || e1 > 0.0 || e2 > 0.0) {
            return false;
        }
        let det: Float = e0 + e1 + e2;
        if det == 0.0 {
            return false;
        }
        // compute scaled hit distance to triangle and test against ray $t$ range
        p0t.z *= sz;
        p1t.z *= sz;
        p2t.z *= sz;
        let t_scaled: Float = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
        if det < 0.0 && (t_scaled >= 0.0 || t_scaled < ray.t_max * det) {
            return false;
        } else if det > 0.0 && (t_scaled <= 0.0 || t_scaled > ray.t_max * det) {
            return false;
        }
        // compute barycentric coordinates and $t$ value for triangle intersection
        let inv_det: Float = 1.0 / det;
        // let b0: Float = e0 * inv_det;
        // let b1: Float = e1 * inv_det;
        // let b2: Float = e2 * inv_det;
        let t: Float = t_scaled * inv_det;

        // ensure that computed triangle $t$ is conservatively greater than zero

        // compute $\delta_z$ term for triangle $t$ error bounds
        let max_zt: Float = vec3_max_component(&Vector3f {
            x: p0t.z,
            y: p1t.z,
            z: p2t.z,
        }.abs());
        let delta_z: Float = gamma(3_i32) * max_zt;
        // compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
        let max_xt: Float = vec3_max_component(&Vector3f {
            x: p0t.x,
            y: p1t.x,
            z: p2t.x,
        }.abs());
        let max_yt: Float = vec3_max_component(&Vector3f {
            x: p0t.y,
            y: p1t.y,
            z: p2t.y,
        }.abs());
        let delta_x: Float = gamma(5) * (max_xt + max_zt);
        let delta_y: Float = gamma(5) * (max_yt + max_zt);
        // compute $\delta_e$ term for triangle $t$ error bounds
        let delta_e: Float =
            2.0 * (gamma(2) * max_xt * max_yt + delta_y * max_xt + delta_x * max_yt);
        // compute $\delta_t$ term for triangle $t$ error bounds and check _t_
        let max_e: Float = vec3_max_component(&Vector3f {
            x: e0,
            y: e1,
            z: e2,
        }.abs());
        let delta_t: Float =
            3.0 * (gamma(3) * max_e * max_zt + delta_e * max_zt + delta_z * max_e) * inv_det.abs();
        if t <= delta_t {
            return false;
        }
        // TODO: if (testAlphaTexture && (mesh->alphaMask || mesh->shadowAlphaMask)) { ... }
        // TODO: ++nHits;
        true
    }
    fn get_reverse_orientation(&self) -> bool {
        self.reverse_orientation
    }
    fn get_transform_swaps_handedness(&self) -> bool {
        self.transform_swaps_handedness
    }
    fn area(&self) -> Float {
        // get triangle vertices in _p0_, _p1_, and _p2_
        let p0: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 0]];
        let p1: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 1]];
        let p2: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 2]];
        0.5 as Float * vec3_cross_vec3(&(p1 - p0), &(p2 - p0)).length()
    }
    fn sample(&self, u: &Point2f, pdf: &mut Float) -> InteractionCommon {
        let b: Point2f = uniform_sample_triangle(u);
        // get triangle vertices in _p0_, _p1_, and _p2_
        let p0: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 0]];
        let p1: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 1]];
        let p2: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 2]];
        let mut it: InteractionCommon = InteractionCommon::default();
        it.p = p0 * b[0] + p1 * b[1] + p2 * (1.0 as Float - b[0] - b[1]);
        // compute surface normal for sampled point on triangle
        it.n = nrm_normalize(&Normal3f::from(vec3_cross_vec3(&(p1 - p0), &(p2 - p0))));
        // ensure correct orientation of the geometric normal; follow
        // the same approach as was used in Triangle::Intersect().
        if !self.mesh.n.is_empty() {
            let ns: Normal3f = Normal3f::from(
                self.mesh.n[self.mesh.vertex_indices[self.id * 3 + 0]] * b[0]
                    + self.mesh.n[self.mesh.vertex_indices[self.id * 3 + 1]] * b[1]
                    + self.mesh.n[self.mesh.vertex_indices[self.id * 3 + 2]]
                        * (1.0 as Float - b[0] - b[1]),
            );
            it.n = nrm_faceforward_nrm(&it.n, &ns);
        } else if self.reverse_orientation ^ self.transform_swaps_handedness {
            it.n *= -1.0 as Float;
        }
        // compute error bounds for sampled point on triangle
        let p_abs_sum: Point3f = pnt3_abs(&(p0 * b[0])) + pnt3_abs(&(p1 * b[1]))
            + pnt3_abs(&(p2 * (1.0 as Float - b[0] - b[1])));
        it.p_error = Vector3f {
            x: p_abs_sum.x,
            y: p_abs_sum.y,
            z: p_abs_sum.z,
        } * gamma(6);
        *pdf = 1.0 as Float / self.area();
        it
    }
    fn sample_with_ref_point(
        &self,
        iref: &InteractionCommon,
        u: &Point2f,
        pdf: &mut Float,
    ) -> InteractionCommon {
        let intr: InteractionCommon = self.sample(u, pdf);
        let mut wi: Vector3f = intr.p - iref.p;
        if wi.length_squared() == 0.0 as Float {
            *pdf = 0.0 as Float;
        } else {
            wi = vec3_normalize(&wi);
            // convert from area measure, as returned by the Sample()
            // call above, to solid angle measure.
            *pdf *= pnt3_distance_squared(&iref.p, &intr.p) / nrm_abs_dot_vec3(&intr.n, &-wi);
            if (*pdf).is_infinite() {
                *pdf = 0.0 as Float;
            }
        }
        intr
    }
    fn pdf_with_ref_point(&self, iref: &Interaction, wi: &Vector3f) -> Float {
        // intersect sample ray with area light geometry
        let ray: Ray = iref.spawn_ray(wi);
        // ignore any alpha textures used for trimming the shape when
        // performing this intersection. Hack for the "San Miguel"
        // scene, where this is used to make an invisible area light.
        if let Some((isect_light, _t_hit)) = self.intersect(&ray) {
            // convert light sample weight to solid angle measure
            let mut pdf: Float = pnt3_distance_squared(&iref.get_p(), &isect_light.p)
                / (nrm_abs_dot_vec3(&isect_light.n, &-(*wi)) * self.area());
            if pdf.is_infinite() {
                pdf = 0.0 as Float;
            }
            pdf
        } else {
            0.0 as Float
        }
    }
}
