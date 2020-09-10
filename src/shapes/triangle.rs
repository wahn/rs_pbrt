// std
use std::cell::Cell;
use std::mem;
use std::sync::Arc;
// pbrt
use crate::core::geometry::{
    bnd3_union_pnt3f, nrm_abs_dot_vec3f, nrm_faceforward_nrm, pnt3_abs, pnt3_distance_squaredf,
    pnt3_permutef, vec3_coordinate_system, vec3_cross_nrm, vec3_cross_vec3, vec3_max_componentf,
    vec3_max_dimensionf, vec3_permutef,
};
use crate::core::geometry::{
    Bounds3f, Normal3f, Point2f, Point3f, Ray, Vector2f, Vector3f, XYEnum,
};
use crate::core::interaction::{Interaction, InteractionCommon, Shading, SurfaceInteraction};
use crate::core::material::Material;
use crate::core::pbrt::gamma;
use crate::core::pbrt::Float;
// use crate::core::sampling::uniform_sample_triangle;
use crate::core::texture::Texture;
use crate::core::transform::Transform;

// see triangle.h

#[derive(Clone)]
pub struct TriangleMesh {
    /// the total number of triangles in the mesh
    pub n_triangles: u32,
    /// vector of vertex indices
    pub vertex_indices: Vec<u32>,
    /// the total number of vertices in the mesh
    pub n_vertices: u32,
    /// vector of *n_vertices* vertex positions
    pub p: Vec<Point3f>,
    /// an optional vector of normal vectors (can be empty)
    pub n: Vec<Normal3f>,
    /// an optional vector of tangent vectors (can be empty)
    pub s: Vec<Vector3f>,
    /// an optional vector of paramtric (u, v) values (texture coordinates)
    pub uv: Vec<Point2f>,
    pub alpha_mask: Option<Arc<dyn Texture<Float> + Send + Sync>>,
    pub shadow_alpha_mask: Option<Arc<dyn Texture<Float> + Send + Sync>>,
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
        n_triangles: u32,
        vertex_indices: Vec<u32>,
        n_vertices: u32,
        p: Vec<Point3f>,
        s: Vec<Vector3f>,
        n: Vec<Normal3f>,
        uv: Vec<Point2f>,
        alpha_mask: Option<Arc<dyn Texture<Float> + Send + Sync>>,
        shadow_alpha_mask: Option<Arc<dyn Texture<Float> + Send + Sync>>,
    ) -> Self {
        TriangleMesh {
            // Shape
            object_to_world,
            world_to_object,
            reverse_orientation,
            transform_swaps_handedness: object_to_world.swaps_handedness(),
            // TriangleMesh
            n_triangles,
            vertex_indices,
            n_vertices,
            p,
            n,
            s,
            uv,
            alpha_mask,
            shadow_alpha_mask,
        }
    }
}

#[derive(Clone)]
pub struct Triangle {
    mesh: Arc<TriangleMesh>,
    pub id: u32,
    // inherited from class Shape (see shape.h)
    pub object_to_world: Transform,
    pub world_to_object: Transform,
    pub reverse_orientation: bool,
    pub transform_swaps_handedness: bool,
    pub material: Option<Arc<Material>>,
}

impl Triangle {
    pub fn new(
        object_to_world: Transform,
        world_to_object: Transform,
        reverse_orientation: bool,
        mesh: Arc<TriangleMesh>,
        tri_number: u32,
    ) -> Self {
        Triangle {
            mesh,
            id: tri_number,
            object_to_world,
            world_to_object,
            reverse_orientation,
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
                self.mesh.uv[self.mesh.vertex_indices[(self.id * 3) as usize] as usize],
                self.mesh.uv[self.mesh.vertex_indices[(self.id * 3) as usize + 1] as usize],
                self.mesh.uv[self.mesh.vertex_indices[(self.id * 3) as usize + 2] as usize],
            ]
        }
    }
    // Shape
    pub fn object_bound(&self) -> Bounds3f {
        let idx1: usize = (self.id * 3) as usize;
        let idx = &self.mesh.vertex_indices[idx1..(idx1 + 3)];
        let p0: &Point3f = &self.mesh.p[idx[0] as usize];
        let p1: &Point3f = &self.mesh.p[idx[1] as usize];
        let p2: &Point3f = &self.mesh.p[idx[2] as usize];
        bnd3_union_pnt3f(
            &Bounds3f::new(
                self.world_to_object.transform_point(p0),
                self.world_to_object.transform_point(p1),
            ),
            &self.world_to_object.transform_point(p2),
        )
    }
    pub fn world_bound(&self) -> Bounds3f {
        let idx1: usize = (self.id * 3) as usize;
        let idx = &self.mesh.vertex_indices[idx1..(idx1 + 3)];
        let p0: &Point3f = &self.mesh.p[idx[0] as usize];
        let p1: &Point3f = &self.mesh.p[idx[1] as usize];
        let p2: &Point3f = &self.mesh.p[idx[2] as usize];
        bnd3_union_pnt3f(&Bounds3f::new(*p0, *p1), p2)
    }
    pub fn intersect(&self, ray: &Ray, t_hit: &mut Float, isect: &mut SurfaceInteraction) -> bool {
        // get triangle vertices in _p0_, _p1_, and _p2_
        let idx1: usize = (self.id * 3) as usize;
        let idx = &self.mesh.vertex_indices[idx1..(idx1 + 3)];
        let p0: &Point3f = &self.mesh.p[idx[0] as usize];
        let p1: &Point3f = &self.mesh.p[idx[1] as usize];
        let p2: &Point3f = &self.mesh.p[idx[2] as usize];
        // translate vertices based on ray origin
        let mut p0t: Point3f = *p0
            - Vector3f {
                x: ray.o.x,
                y: ray.o.y,
                z: ray.o.z,
            };
        let mut p1t: Point3f = *p1
            - Vector3f {
                x: ray.o.x,
                y: ray.o.y,
                z: ray.o.z,
            };
        let mut p2t: Point3f = *p2
            - Vector3f {
                x: ray.o.x,
                y: ray.o.y,
                z: ray.o.z,
            };
        // permute components of triangle vertices and ray direction
        let kz: usize = vec3_max_dimensionf(&ray.d.abs());
        let mut kx: usize = kz + 1;
        if kx == 3 {
            kx = 0;
        }
        let mut ky: usize = kx + 1;
        if ky == 3 {
            ky = 0;
        }
        let d: Vector3f = vec3_permutef(&ray.d, kx, ky, kz);
        p0t = pnt3_permutef(&p0t, kx, ky, kz);
        p1t = pnt3_permutef(&p1t, kx, ky, kz);
        p2t = pnt3_permutef(&p2t, kx, ky, kz);
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
        let b0: Float = e0 * inv_det;
        let b1: Float = e1 * inv_det;
        let b2: Float = e2 * inv_det;
        let t: Float = t_scaled * inv_det;

        // ensure that computed triangle $t$ is conservatively greater than zero

        // compute $\delta_z$ term for triangle $t$ error bounds
        let max_zt: Float = vec3_max_componentf(
            &Vector3f {
                x: p0t.z,
                y: p1t.z,
                z: p2t.z,
            }
            .abs(),
        );
        let delta_z: Float = gamma(3_i32) * max_zt;
        // compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
        let max_xt: Float = vec3_max_componentf(
            &Vector3f {
                x: p0t.x,
                y: p1t.x,
                z: p2t.x,
            }
            .abs(),
        );
        let max_yt: Float = vec3_max_componentf(
            &Vector3f {
                x: p0t.y,
                y: p1t.y,
                z: p2t.y,
            }
            .abs(),
        );
        let delta_x: Float = gamma(5) * (max_xt + max_zt);
        let delta_y: Float = gamma(5) * (max_yt + max_zt);
        // compute $\delta_e$ term for triangle $t$ error bounds
        let delta_e: Float =
            2.0 * (gamma(2) * max_xt * max_yt + delta_y * max_xt + delta_x * max_yt);
        // compute $\delta_t$ term for triangle $t$ error bounds and check _t_
        let max_e: Float = vec3_max_componentf(
            &Vector3f {
                x: e0,
                y: e1,
                z: e2,
            }
            .abs(),
        );
        let delta_t: Float =
            3.0 * (gamma(3) * max_e * max_zt + delta_e * max_zt + delta_z * max_e) * inv_det.abs();
        if t <= delta_t {
            return false;
        }
        // compute triangle partial derivatives
        let uv: [Point2f; 3] = self.get_uvs();
        // compute deltas for triangle partial derivatives
        let duv02: Vector2f = uv[0] - uv[2];
        let duv12: Vector2f = uv[1] - uv[2];
        let dp02: Vector3f = *p0 - *p2;
        let dp12: Vector3f = *p1 - *p2;
        let determinant: Float = duv02.x * duv12.y - duv02.y * duv12.x;
        let degenerate_uv: bool = determinant.abs() < 1e-8 as Float;
        let mut dpdu: Vector3f = Vector3f::default();
        let mut dpdv: Vector3f = if !degenerate_uv {
            let invdet: Float = 1.0 / determinant;
            dpdu = (dp02 * duv12.y - dp12 * duv02.y) * invdet;
            (dp02 * -duv12.x + dp12 * duv02.x) * invdet
        } else {
            Vector3f::default()
        };
        if degenerate_uv || vec3_cross_vec3(&dpdu, &dpdv).length_squared() == 0.0 {
            // handle zero determinant for triangle partial derivative matrix
            vec3_coordinate_system(
                &vec3_cross_vec3(&(*p2 - *p0), &(*p1 - *p0)).normalize(),
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
        let p_hit: Point3f = *p0 * b0 + *p1 * b1 + *p2 * b2;
        let uv_hit: Point2f = uv[0] * b0 + uv[1] * b1 + uv[2] * b2;
        // test intersection against alpha texture, if present
        // TODO: testAlphaTexture
        if let Some(alpha_mask) = &self.mesh.alpha_mask {
            let wo: Vector3f = -ray.d;
            let isect_local: SurfaceInteraction = SurfaceInteraction::new(
                &p_hit,
                &Vector3f::default(),
                uv_hit,
                &wo,
                &dpdu,
                &dpdv,
                &Normal3f::default(),
                &Normal3f::default(),
                ray.time,
                None,
            );
            if alpha_mask.evaluate(&isect_local) == 0.0 as Float {
                return false;
            }
        }
        // fill in _SurfaceInteraction_ from triangle hit
        let dndu: Normal3f = Normal3f::default();
        let dndv: Normal3f = Normal3f::default();
        let wo: Vector3f = -ray.d;
        // override surface normal in _isect_ for triangle
        let mut surface_normal: Normal3f =
            Normal3f::from(vec3_cross_vec3(&dp02, &dp12).normalize());
        let mut shading: Shading = Shading {
            n: surface_normal,
            dpdu,
            dpdv,
            dndu,
            dndv,
        };
        if !self.mesh.n.is_empty() || !self.mesh.s.is_empty() {
            // initialize _Triangle_ shading geometry

            // compute shading normal _ns_ for triangle
            let mut ns: Normal3f;
            if !self.mesh.n.is_empty() {
                let n0 = self.mesh.n[self.mesh.vertex_indices[(self.id * 3) as usize] as usize];
                let n1 = self.mesh.n[self.mesh.vertex_indices[(self.id * 3) as usize + 1] as usize];
                let n2 = self.mesh.n[self.mesh.vertex_indices[(self.id * 3) as usize + 2] as usize];
                ns = n0 * b0 + n1 * b1 + n2 * b2;
                if ns.length_squared() > 0.0 {
                    ns = ns.normalize();
                } else {
                    ns = surface_normal;
                }
            } else {
                ns = surface_normal;
            }
            // compute shading tangent _ss_ for triangle
            let mut ss: Vector3f;
            if !self.mesh.s.is_empty() {
                let s0 = self.mesh.s[self.mesh.vertex_indices[(self.id * 3) as usize] as usize];
                let s1 = self.mesh.s[self.mesh.vertex_indices[(self.id * 3) as usize + 1] as usize];
                let s2 = self.mesh.s[self.mesh.vertex_indices[(self.id * 3) as usize + 2] as usize];
                ss = s0 * b0 + s1 * b1 + s2 * b2;
                if ss.length_squared() > 0.0 {
                    ss = ss.normalize();
                } else {
                    ss = dpdu.normalize();
                }
            } else {
                ss = dpdu.normalize();
            }
            // compute shading bitangent _ts_ for triangle and adjust _ss_
            let mut ts: Vector3f = vec3_cross_nrm(&ss, &ns);
            if ts.length_squared() > 0.0 {
                ts = ts.normalize();
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
                let dn1: Normal3f = self.mesh.n
                    [self.mesh.vertex_indices[(self.id * 3) as usize] as usize]
                    - self.mesh.n[self.mesh.vertex_indices[(self.id * 3) as usize + 2] as usize];
                let dn2: Normal3f = self.mesh.n
                    [self.mesh.vertex_indices[(self.id * 3) as usize + 1] as usize]
                    - self.mesh.n[self.mesh.vertex_indices[(self.id * 3) as usize + 2] as usize];
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
            shading.n = Normal3f::from(vec3_cross_vec3(&ss, &ts)).normalize();
            surface_normal = nrm_faceforward_nrm(&surface_normal, &shading.n);
            shading.dpdu = ss;
            shading.dpdv = ts;
            shading.dndu = dndu;
            shading.dndv = dndv;
        }
        {
            isect.common.p = p_hit;
            isect.common.time = ray.time;
            isect.common.p_error = p_error;
            isect.common.wo = wo;
            isect.common.n = surface_normal;
            isect.common.medium_interface = None;
        }
        isect.uv = uv_hit;
        isect.dpdu = dpdu;
        isect.dpdv = dpdv;
        isect.dndu = dndu;
        isect.dndv = dndv;
        isect.dpdx = Cell::new(Vector3f::default());
        isect.dpdy = Cell::new(Vector3f::default());
        isect.dudx = Cell::new(0.0 as Float);
        isect.dvdx = Cell::new(0.0 as Float);
        isect.dudy = Cell::new(0.0 as Float);
        isect.dvdy = Cell::new(0.0 as Float);
        isect.primitive = None;
        isect.shading = shading;
        isect.bsdf = None;
        // isect.bssrdf = None;
        isect.shape = None;
        *t_hit = t;
        true
    }
    pub fn intersect_p(&self, ray: &Ray) -> bool {
        // TODO: ProfilePhase p(Prof::TriIntersectP);
        // TODO: ++nTests;
        // get triangle vertices in _p0_, _p1_, and _p2_
        let idx1: usize = (self.id * 3) as usize;
        let idx = &self.mesh.vertex_indices[idx1..(idx1 + 3)];
        let p0: &Point3f = &self.mesh.p[idx[0] as usize];
        let p1: &Point3f = &self.mesh.p[idx[1] as usize];
        let p2: &Point3f = &self.mesh.p[idx[2] as usize];
        // translate vertices based on ray origin
        let mut p0t: Point3f = *p0
            - Vector3f {
                x: ray.o.x,
                y: ray.o.y,
                z: ray.o.z,
            };
        let mut p1t: Point3f = *p1
            - Vector3f {
                x: ray.o.x,
                y: ray.o.y,
                z: ray.o.z,
            };
        let mut p2t: Point3f = *p2
            - Vector3f {
                x: ray.o.x,
                y: ray.o.y,
                z: ray.o.z,
            };
        // permute components of triangle vertices and ray direction
        let kz: usize = vec3_max_dimensionf(&ray.d.abs());
        let mut kx: usize = kz + 1;
        if kx == 3 {
            kx = 0;
        }
        let mut ky: usize = kx + 1;
        if ky == 3 {
            ky = 0;
        }
        let d: Vector3f = vec3_permutef(&ray.d, kx, ky, kz);
        p0t = pnt3_permutef(&p0t, kx, ky, kz);
        p1t = pnt3_permutef(&p1t, kx, ky, kz);
        p2t = pnt3_permutef(&p2t, kx, ky, kz);
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
        let b0: Float = e0 * inv_det;
        let b1: Float = e1 * inv_det;
        let b2: Float = e2 * inv_det;
        let t: Float = t_scaled * inv_det;

        // ensure that computed triangle $t$ is conservatively greater than zero

        // compute $\delta_z$ term for triangle $t$ error bounds
        let max_zt: Float = vec3_max_componentf(
            &Vector3f {
                x: p0t.z,
                y: p1t.z,
                z: p2t.z,
            }
            .abs(),
        );
        let delta_z: Float = gamma(3_i32) * max_zt;
        // compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
        let max_xt: Float = vec3_max_componentf(
            &Vector3f {
                x: p0t.x,
                y: p1t.x,
                z: p2t.x,
            }
            .abs(),
        );
        let max_yt: Float = vec3_max_componentf(
            &Vector3f {
                x: p0t.y,
                y: p1t.y,
                z: p2t.y,
            }
            .abs(),
        );
        let delta_x: Float = gamma(5) * (max_xt + max_zt);
        let delta_y: Float = gamma(5) * (max_yt + max_zt);
        // compute $\delta_e$ term for triangle $t$ error bounds
        let delta_e: Float =
            2.0 * (gamma(2) * max_xt * max_yt + delta_y * max_xt + delta_x * max_yt);
        // compute $\delta_t$ term for triangle $t$ error bounds and check _t_
        let max_e: Float = vec3_max_componentf(
            &Vector3f {
                x: e0,
                y: e1,
                z: e2,
            }
            .abs(),
        );
        let delta_t: Float =
            3.0 * (gamma(3) * max_e * max_zt + delta_e * max_zt + delta_z * max_e) * inv_det.abs();
        if t <= delta_t {
            return false;
        }
        // TODO: if (testAlphaTexture && (mesh->alphaMask || mesh->shadowAlphaMask)) { ... }
        if self.mesh.alpha_mask.is_some() || self.mesh.shadow_alpha_mask.is_some() {
            // compute triangle partial derivatives
            let mut dpdu: Vector3f = Vector3f::default();
            let mut dpdv: Vector3f = Vector3f::default();
            let uv: [Point2f; 3] = self.get_uvs();
            // compute deltas for triangle partial derivatives
            let duv02: Vector2f = uv[0] - uv[2];
            let duv12: Vector2f = uv[1] - uv[2];
            let dp02: Vector3f = *p0 - *p2;
            let dp12: Vector3f = *p1 - *p2;
            let duv02x = duv02[XYEnum::X];
            let duv02y = duv02[XYEnum::Y];
            let duv12x = duv12[XYEnum::X];
            let duv12y = duv12[XYEnum::Y];
            let determinant: Float = duv02x * duv12y - duv02y * duv12x;
            let degenerate_uv: bool = determinant.abs() < 1e-8 as Float;
            if !degenerate_uv {
                let invdet: Float = 1.0 as Float / determinant;
                dpdu = (dp02 * duv12y - dp12 * duv02y) * invdet;
                dpdv = (dp02 * -duv12x + dp12 * duv02x) * invdet;
            }
            if degenerate_uv || vec3_cross_vec3(&dpdu, &dpdv).length_squared() == 0.0 {
                // handle zero determinant for triangle partial derivative matrix
                let ng = vec3_cross_vec3(&(*p2 - *p0), &(*p1 - *p0));
                if ng.length_squared() == 0.0 as Float {
                    // the triangle is actually degenerate; the
                    // intersection is bogus
                    return false;
                }
                vec3_coordinate_system(
                    &vec3_cross_vec3(&(*p2 - *p0), &(*p1 - *p0)).normalize(),
                    &mut dpdu,
                    &mut dpdv,
                );
            }
            // interpolate $(u,v)$ parametric coordinates and hit point
            let p_hit: Point3f = *p0 * b0 + *p1 * b1 + *p2 * b2;
            let uv_hit: Point2f = uv[0] * b0 + uv[1] * b1 + uv[2] * b2;
            let wo: Vector3f = -ray.d;
            let isect_local: SurfaceInteraction = SurfaceInteraction::new(
                &p_hit,
                &Vector3f::default(),
                uv_hit,
                &wo,
                &dpdu,
                &dpdv,
                &Normal3f::default(),
                &Normal3f::default(),
                ray.time,
                None,
            );
            if let Some(alpha_mask) = &self.mesh.alpha_mask {
                if alpha_mask.evaluate(&isect_local) == 0.0 as Float {
                    return false;
                }
            }
            if let Some(shadow_alpha_mask) = &self.mesh.shadow_alpha_mask {
                if shadow_alpha_mask.evaluate(&isect_local) == 0.0 as Float {
                    return false;
                }
            }
        }
        // TODO: ++nHits;
        true
    }
    pub fn get_reverse_orientation(&self) -> bool {
        self.reverse_orientation
    }
    pub fn get_transform_swaps_handedness(&self) -> bool {
        self.transform_swaps_handedness
    }
    pub fn get_object_to_world(&self) -> Transform {
        self.object_to_world
    }
    pub fn area(&self) -> Float {
        let idx1: usize = (self.id * 3) as usize;
        let idx = &self.mesh.vertex_indices[idx1..(idx1 + 3)];
        // get triangle vertices in _p0_, _p1_, and _p2_
        let p0: &Point3f = &self.mesh.p[idx[0] as usize];
        let p1: &Point3f = &self.mesh.p[idx[1] as usize];
        let p2: &Point3f = &self.mesh.p[idx[2] as usize];
        0.5 as Float * vec3_cross_vec3(&(*p1 - *p0), &(*p2 - *p0)).length()
    }
    pub fn sample(&self, u: Point2f, pdf: &mut Float) -> InteractionCommon {
        let idx1: usize = (self.id * 3) as usize;
        let idx = &self.mesh.vertex_indices[idx1..(idx1 + 3)];
        // avoid calling uniform_sample_triangle!!!
        // let b: Point2f = uniform_sample_triangle(u);
        let su0: Float = u[XYEnum::X].sqrt();
        let bx: Float = 1.0 as Float - su0;
        let by: Float = u[XYEnum::Y] * su0;
        // get triangle vertices in _p0_, _p1_, and _p2_
        let p0: &Point3f = &self.mesh.p[idx[0] as usize];
        let p1: &Point3f = &self.mesh.p[idx[1] as usize];
        let p2: &Point3f = &self.mesh.p[idx[2] as usize];
        // let bx = b[XYEnum::X];
        // let by = b[XYEnum::Y];
        let it_p = *p0 * bx + *p1 * by + *p2 * (1.0 as Float - bx - by);
        // compute surface normal for sampled point on triangle
        let mut it_n = Normal3f::from(vec3_cross_vec3(&(*p1 - *p0), &(*p2 - *p0))).normalize();
        // ensure correct orientation of the geometric normal; follow
        // the same approach as was used in Triangle::Intersect().
        if !self.mesh.n.is_empty() {
            let ns: Normal3f = self.mesh.n[idx[0] as usize] * bx
                + self.mesh.n[idx[1] as usize] * by
                + self.mesh.n[idx[2] as usize] * (1.0 as Float - bx - by);
            it_n = nrm_faceforward_nrm(&it_n, &ns);
        } else if self.reverse_orientation ^ self.transform_swaps_handedness {
            it_n *= -1.0 as Float;
        }
        // compute error bounds for sampled point on triangle
        let p_abs_sum: Point3f = pnt3_abs(&(*p0 * bx))
            + pnt3_abs(&(*p1 * by))
            + pnt3_abs(&(*p2 * (1.0 as Float - bx - by)));
        let it_p_error = Vector3f {
            x: p_abs_sum.x,
            y: p_abs_sum.y,
            z: p_abs_sum.z,
        } * gamma(6);
        // avoid calling self.area()!!! *pdf = 1.0 as Float / self.area();
        let area: Float = 0.5 as Float * vec3_cross_vec3(&(*p1 - *p0), &(*p2 - *p0)).length();
        *pdf = 1.0 as Float / area;
        InteractionCommon {
            p: it_p,
            time: 0.0 as Float,
            p_error: it_p_error,
            wo: Vector3f::default(),
            n: it_n,
            medium_interface: None,
        }
    }
    pub fn sample_with_ref_point(
        &self,
        iref: &InteractionCommon,
        u: Point2f,
        pdf: &mut Float,
    ) -> InteractionCommon {
        let intr: InteractionCommon = self.sample(u, pdf);
        let mut wi: Vector3f = intr.p - iref.p;
        if wi.length_squared() == 0.0 as Float {
            *pdf = 0.0 as Float;
        } else {
            wi = wi.normalize();
            // convert from area measure, as returned by the Sample()
            // call above, to solid angle measure.
            *pdf *= pnt3_distance_squaredf(&iref.p, &intr.p) / nrm_abs_dot_vec3f(&intr.n, &-wi);
            if (*pdf).is_infinite() {
                *pdf = 0.0 as Float;
            }
        }
        intr
    }
    pub fn pdf_with_ref_point(&self, iref: &dyn Interaction, wi: &Vector3f) -> Float {
        // intersect sample ray with area light geometry
        let ray: Ray = iref.spawn_ray(wi);
        // ignore any alpha textures used for trimming the shape when
        // performing this intersection. Hack for the "San Miguel"
        // scene, where this is used to make an invisible area light.
        let mut t_hit: Float = 0.0;
        let mut isect_light: SurfaceInteraction = SurfaceInteraction::default();
        if self.intersect(&ray, &mut t_hit, &mut isect_light) {
            // convert light sample weight to solid angle measure
            let mut pdf: Float = pnt3_distance_squaredf(&iref.get_p(), &isect_light.common.p)
                / (nrm_abs_dot_vec3f(&isect_light.common.n, &-(*wi)) * self.area());
            if pdf.is_infinite() {
                pdf = 0.0 as Float;
            }
            pdf
        } else {
            0.0 as Float
        }
    }
}
