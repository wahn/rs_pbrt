//! The geometry of a particular point on a surface is represented by
//! a **SurfaceInteraction**. Having this abstraction lets most of the
//! system work with points on surfaces without needing to consider
//! the particular type of geometric shape the points lie on; the
//! **SurfaceInteraction** abstraction supplies enough information
//! about the surface point to allow the shading and geometric
//! operations in the rest of **pbrt** to be implemented generically.
//!

// std
use std;
use std::sync::Arc;
// pbrt
use core::geometry::{Normal3f, Point2f, Point3f, Ray, Vector3f};
use core::geometry::{nrm_faceforward_nrm, nrm_normalize, pnt3_offset_ray_origin, vec3_cross_vec3,
                     vec3_dot_vec3, vec3_normalize};
use core::pbrt::SHADOW_EPSILON;
use core::pbrt::{Float, Spectrum};
use core::material::TransportMode;
use core::primitive::{GeometricPrimitive, Primitive};
use core::reflection::Bsdf;
use core::shape::Shape;
use core::transform::solve_linear_system_2x2;

// see interaction.h

pub trait Interaction {
    fn is_surface_interaction(&self) -> bool;
    fn is_medium_interaction(&self) -> bool;
    fn spawn_ray(&self, d: &Vector3f) -> Ray;
    fn get_p(&self) -> Point3f;
    fn get_time(&self) -> Float;
    fn get_p_error(&self) -> Vector3f;
    fn get_wo(&self) -> Vector3f;
    fn get_n(&self) -> Normal3f;
}

#[derive(Debug, Default, Copy, Clone)]
pub struct InteractionCommon {
    // Interaction Public Data
    pub p: Point3f,
    pub time: Float,
    pub p_error: Vector3f,
    pub wo: Vector3f,
    pub n: Normal3f,
}

impl InteractionCommon {
    pub fn spawn_ray_to(&self, it: &InteractionCommon) -> Ray {
        let origin: Point3f =
            pnt3_offset_ray_origin(&self.p, &self.p_error, &self.n, &(it.p - self.p));
        let target: Point3f = pnt3_offset_ray_origin(&it.p, &it.p_error, &it.n, &(origin - it.p));
        let d: Vector3f = target - origin;
        Ray {
            o: origin,
            d: d,
            t_max: 1.0 - SHADOW_EPSILON,
            time: self.time,
            differential: None,
        }
    }
}

#[derive(Debug, Default, Copy, Clone)]
pub struct Shading {
    pub n: Normal3f,
    pub dpdu: Vector3f,
    pub dpdv: Vector3f,
    pub dndu: Normal3f,
    pub dndv: Normal3f,
}

#[derive(Default, Clone)]
pub struct SurfaceInteraction<'a, 'b> {
    // Interaction Public Data
    pub p: Point3f,
    pub time: Float,
    pub p_error: Vector3f,
    pub wo: Vector3f,
    pub n: Normal3f,
    // TODO: MediumInterface mediumInterface;
    // SurfaceInteraction Public Data
    pub uv: Point2f,
    pub dpdu: Vector3f,
    pub dpdv: Vector3f,
    pub dndu: Normal3f,
    pub dndv: Normal3f,
    pub dpdx: Vector3f,
    pub dpdy: Vector3f,
    pub dudx: Float,
    pub dvdx: Float,
    pub dudy: Float,
    pub dvdy: Float,
    pub primitive: Option<&'a GeometricPrimitive>,
    pub shading: Shading,
    pub bsdf: Option<Arc<Bsdf>>,
    pub shape: Option<&'b Shape>,
}

impl<'a, 'b> SurfaceInteraction<'a, 'b> {
    pub fn new(
        p: &Point3f,
        p_error: &Vector3f,
        uv: &Point2f,
        wo: &Vector3f,
        dpdu: &Vector3f,
        dpdv: &Vector3f,
        dndu: &Normal3f,
        dndv: &Normal3f,
        time: Float,
        sh: Option<&'b Shape>,
    ) -> Self {
        let nv: Vector3f = vec3_normalize(&vec3_cross_vec3(dpdu, dpdv));
        // TODO: Adjust normal based on orientation and handedness
        let n: Normal3f = Normal3f {
            x: nv.x,
            y: nv.y,
            z: nv.z,
        };
        // initialize shading geometry from true geometry
        let shading: Shading = Shading {
            n: n,
            dpdu: Vector3f::from(*dpdu),
            dpdv: Vector3f::from(*dpdv),
            dndu: *dndu,
            dndv: *dndv,
        };
        SurfaceInteraction {
            p: *p,
            time: time,
            p_error: *p_error,
            wo: vec3_normalize(wo),
            n: n,
            uv: *uv,
            dpdu: *dpdu,
            dpdv: *dpdv,
            dndu: *dndu,
            dndv: *dndv,
            dpdx: Vector3f::default(),
            dpdy: Vector3f::default(),
            dudx: 0.0 as Float,
            dvdx: 0.0 as Float,
            dudy: 0.0 as Float,
            dvdy: 0.0 as Float,
            primitive: None,
            shading: shading,
            bsdf: None,
            shape: sh,
        }
    }
    pub fn set_shading_geometry(
        &mut self,
        dpdus: &Vector3f,
        dpdvs: &Vector3f,
        dndus: &Normal3f,
        dndvs: &Normal3f,
        orientation_is_authoritative: bool,
    ) {
        // compute _shading.n_ for _SurfaceInteraction_
        self.shading.n = nrm_normalize(&Normal3f::from(vec3_cross_vec3(dpdus, dpdvs)));
        if let Some(shape) = self.shape {
            if shape.get_reverse_orientation() ^ shape.get_transform_swaps_handedness() {
                self.shading.n = -self.shading.n;
            }
        }
        if orientation_is_authoritative {
            self.n = nrm_faceforward_nrm(&self.n, &self.shading.n);
        } else {
            self.shading.n = nrm_faceforward_nrm(&self.shading.n, &self.n);
        }
        // initialize _shading_ partial derivative values
        self.shading.dpdu = *dpdus;
        self.shading.dpdv = *dpdvs;
        self.shading.dndu = *dndus;
        self.shading.dndv = *dndvs;
    }
    pub fn compute_scattering_functions(
        &mut self,
        ray: &Ray,
        // arena: &mut Arena,
        allow_multiple_lobes: bool,
        mode: TransportMode,
    ) {
        self.compute_differentials(ray);
        if let Some(primitive) = self.primitive {
            primitive.compute_scattering_functions(
                self, // arena,
                mode,
                allow_multiple_lobes,
            );
        }
    }
    pub fn compute_differentials(&mut self, ray: &Ray) {
        if let Some(ref diff) = ray.differential {
            // estimate screen space change in $\pt{}$ and $(u,v)$

            // compute auxiliary intersection points with plane
            let d: Float = vec3_dot_vec3(
                &Vector3f::from(self.n),
                &Vector3f {
                    x: self.p.x,
                    y: self.p.y,
                    z: self.p.z,
                },
            );
            let tx: Float =
                -(vec3_dot_vec3(&Vector3f::from(self.n), &Vector3f::from(diff.rx_origin)) - d)
                    / vec3_dot_vec3(&Vector3f::from(self.n), &diff.rx_direction);
            if tx.is_nan() {
                self.dudx = 0.0 as Float;
                self.dvdx = 0.0 as Float;
                self.dudy = 0.0 as Float;
                self.dvdy = 0.0 as Float;
                self.dpdx = Vector3f::default();
                self.dpdy = Vector3f::default();
            } else {
                let px: Point3f = diff.rx_origin + diff.rx_direction * tx;
                let ty: Float =
                    -(vec3_dot_vec3(&Vector3f::from(self.n), &Vector3f::from(diff.ry_origin)) - d)
                        / vec3_dot_vec3(&Vector3f::from(self.n), &diff.ry_direction);
                if ty.is_nan() {
                    self.dudx = 0.0 as Float;
                    self.dvdx = 0.0 as Float;
                    self.dudy = 0.0 as Float;
                    self.dvdy = 0.0 as Float;
                    self.dpdx = Vector3f::default();
                    self.dpdy = Vector3f::default();
                } else {
                    let py: Point3f = diff.ry_origin + diff.ry_direction * ty;
                    self.dpdx = px - self.p;
                    self.dpdy = py - self.p;

                    // compute $(u,v)$ offsets at auxiliary points

                    // choose two dimensions to use for ray offset computation
                    let mut dim: [u8; 2] = [0_u8; 2];
                    if self.n.x.abs() > self.n.y.abs() && self.n.x.abs() > self.n.z.abs() {
                        dim[0] = 1;
                        dim[1] = 2;
                    } else if self.n.y.abs() > self.n.z.abs() {
                        dim[0] = 0;
                        dim[1] = 2;
                    } else {
                        dim[0] = 0;
                        dim[1] = 1;
                    }

                    // initialize _a_, _bx_, and _by_ matrices for offset computation
                    let a: [[Float; 2]; 2] = [
                        [self.dpdu[dim[0]], self.dpdv[dim[0]]],
                        [self.dpdu[dim[1]], self.dpdv[dim[1]]],
                    ];
                    let bx: [Float; 2] = [px[dim[0]] - self.p[dim[0]], px[dim[1]] - self.p[dim[1]]];
                    let by: [Float; 2] = [py[dim[0]] - self.p[dim[0]], py[dim[1]] - self.p[dim[1]]];
                    if !solve_linear_system_2x2(a, bx, &mut self.dudx, &mut self.dvdx) {
                        self.dudx = 0.0 as Float;
                        self.dvdx = 0.0 as Float;
                    }
                    if !solve_linear_system_2x2(a, by, &mut self.dudy, &mut self.dvdy) {
                        self.dudy = 0.0 as Float;
                        self.dvdy = 0.0 as Float;
                    }
                }
            }
        } else {
            self.dudx = 0.0 as Float;
            self.dvdx = 0.0 as Float;
            self.dudy = 0.0 as Float;
            self.dvdy = 0.0 as Float;
            self.dpdx = Vector3f::default();
            self.dpdy = Vector3f::default();
        }
    }
    pub fn le(&self, w: &Vector3f) -> Spectrum {
        if let Some(primitive) = self.primitive {
            if let Some(area_light) = primitive.get_area_light() {
                // create InteractionCommon from self
                let interaction: InteractionCommon = InteractionCommon {
                    p: self.p,
                    time: self.time,
                    p_error: self.p_error,
                    wo: self.wo,
                    n: self.n,
                };
                return area_light.l(&interaction, *w);
            }
        }
        Spectrum::default()
    }
}

impl<'a, 'b> Interaction for SurfaceInteraction<'a, 'b> {
    fn is_surface_interaction(&self) -> bool {
        self.n != Normal3f::default()
    }
    fn is_medium_interaction(&self) -> bool {
        !self.is_surface_interaction()
    }
    fn spawn_ray(&self, d: &Vector3f) -> Ray {
        let o: Point3f = pnt3_offset_ray_origin(&self.p, &self.p_error, &self.n, d);
        Ray {
            o: o,
            d: *d,
            t_max: std::f32::INFINITY,
            time: self.time,
            differential: None,
        }
    }
    fn get_p(&self) -> Point3f {
        self.p.clone()
    }
    fn get_time(&self) -> Float {
        self.time
    }
    fn get_p_error(&self) -> Vector3f {
        self.p_error.clone()
    }
    fn get_wo(&self) -> Vector3f {
        self.wo.clone()
    }
    fn get_n(&self) -> Normal3f {
        self.n.clone()
    }
}
