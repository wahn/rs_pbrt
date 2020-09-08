//! The geometry of a particular point on a surface is represented by
//! a **SurfaceInteraction**. Having this abstraction lets most of the
//! system work with points on surfaces without needing to consider
//! the particular type of geometric shape the points lie on; the
//! **SurfaceInteraction** abstraction supplies enough information
//! about the surface point to allow the shading and geometric
//! operations in the rest of **pbrt** to be implemented generically.
//!

// std
use std::cell::Cell;
use std::rc::Rc;
use std::sync::Arc;
// pbrt
use crate::core::bssrdf::TabulatedBssrdf;
use crate::core::geometry::{
    nrm_faceforward_nrm, pnt3_offset_ray_origin, vec3_cross_vec3, vec3_dot_nrm, vec3_dot_vec3,
};
use crate::core::geometry::{Normal3f, Point2f, Point3f, Ray, Vector3f, XYZEnum};
use crate::core::material::TransportMode;
use crate::core::medium::{HenyeyGreenstein, Medium, MediumInterface};
use crate::core::pbrt::SHADOW_EPSILON;
use crate::core::pbrt::{Float, Spectrum};
use crate::core::primitive::Primitive;
use crate::core::reflection::Bsdf;
use crate::core::shape::Shape;
use crate::core::transform::solve_linear_system_2x2;

// see interaction.h

pub trait Interaction {
    fn is_surface_interaction(&self) -> bool;
    fn is_medium_interaction(&self) -> bool;
    fn spawn_ray(&self, d: &Vector3f) -> Ray;
    fn get_common(&self) -> Rc<InteractionCommon>;
    fn get_p(&self) -> &Point3f;
    fn get_time(&self) -> Float;
    fn get_p_error(&self) -> &Vector3f;
    fn get_wo(&self) -> &Vector3f;
    fn get_n(&self) -> &Normal3f;
    fn get_medium_interface(&self) -> Option<Arc<MediumInterface>>;
    fn get_bsdf(&self) -> Option<&Bsdf>;
    fn get_shading_n(&self) -> Option<&Normal3f>;
    fn get_phase(&self) -> Option<Arc<HenyeyGreenstein>>;
}

#[derive(Default, Clone)]
pub struct InteractionCommon {
    // Interaction Public Data
    pub p: Point3f,
    pub time: Float,
    pub p_error: Vector3f,
    pub wo: Vector3f,
    pub n: Normal3f,
    pub medium_interface: Option<Arc<MediumInterface>>,
}

impl InteractionCommon {
    pub fn spawn_ray(&self, d: &Vector3f) -> Ray {
        let o: Point3f = pnt3_offset_ray_origin(&self.p, &self.p_error, &self.n, d);
        Ray {
            o,
            d: *d,
            t_max: std::f32::INFINITY,
            time: self.time,
            differential: None,
            medium: self.get_medium(d),
        }
    }
    pub fn spawn_ray_to_pnt(&self, p2: &Point3f) -> Ray {
        let d: Vector3f = *p2 - self.p;
        let origin: Point3f = pnt3_offset_ray_origin(&self.p, &self.p_error, &self.n, &d);
        Ray {
            o: origin,
            d,
            t_max: 1.0 - SHADOW_EPSILON,
            time: self.time,
            differential: None,
            medium: self.get_medium(&d),
        }
    }
    pub fn spawn_ray_to(&self, it: &InteractionCommon) -> Ray {
        let origin: Point3f =
            pnt3_offset_ray_origin(&self.p, &self.p_error, &self.n, &(it.p - self.p));
        let target: Point3f = pnt3_offset_ray_origin(&it.p, &it.p_error, &it.n, &(origin - it.p));
        let d: Vector3f = target - origin;
        Ray {
            o: origin,
            d,
            t_max: 1.0 - SHADOW_EPSILON,
            time: self.time,
            differential: None,
            medium: self.get_medium(&d),
        }
    }
    pub fn get_medium(&self, w: &Vector3f) -> Option<Arc<Medium>> {
        if vec3_dot_nrm(w, &self.n) > 0.0 as Float {
            if let Some(ref medium_interface_arc) = self.medium_interface {
                medium_interface_arc.get_outside()
            } else {
                None
            }
        } else if let Some(ref medium_interface_arc) = self.medium_interface {
            medium_interface_arc.get_inside()
        } else {
            None
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
pub struct MediumInteraction {
    // Interaction Public Data
    pub common: Rc<InteractionCommon>,
    // MediumInteraction Public Data
    pub phase: Option<Arc<HenyeyGreenstein>>,
}

impl MediumInteraction {
    pub fn new(
        p: &Point3f,
        wo: &Vector3f,
        time: Float,
        medium: Option<Arc<Medium>>,
        phase: Option<Arc<HenyeyGreenstein>>,
    ) -> Self {
        let mut common: InteractionCommon = InteractionCommon::default();
        common.p = *p;
        common.time = time;
        common.p_error = Vector3f::default();
        common.wo = *wo;
        common.n = Normal3f::default();
        if let Some(medium_arc) = medium {
            let inside: Option<Arc<Medium>> = Some(medium_arc.clone());
            let outside: Option<Arc<Medium>> = Some(medium_arc);
            common.medium_interface = Some(Arc::new(MediumInterface::new(inside, outside)));
            MediumInteraction {
                common: Rc::new(common),
                phase,
            }
        } else {
            common.medium_interface = None;
            MediumInteraction {
                common: Rc::new(common),
                phase,
            }
        }
    }
    pub fn get_medium(&self, w: &Vector3f) -> Option<Arc<Medium>> {
        if vec3_dot_nrm(w, &self.get_n()) > 0.0 as Float {
            if let Some(ref medium_interface) = self.get_medium_interface() {
                if let Some(ref outside_arc) = medium_interface.outside {
                    Some(outside_arc.clone())
                } else {
                    None
                }
            } else {
                None
            }
        } else if let Some(ref medium_interface) = self.get_medium_interface() {
            if let Some(ref inside_arc) = medium_interface.inside {
                Some(inside_arc.clone())
            } else {
                None
            }
        } else {
            None
        }
    }
    pub fn is_valid(&self) -> bool {
        if let Some(ref _arc) = self.phase {
            true
        } else {
            false
        }
    }
    pub fn get_phase(&self) -> Option<Arc<HenyeyGreenstein>> {
        if let Some(ref phase) = self.phase {
            Some(phase.clone())
        } else {
            None
        }
    }
}

impl Interaction for MediumInteraction {
    fn is_surface_interaction(&self) -> bool {
        self.common.n != Normal3f::default()
    }
    fn is_medium_interaction(&self) -> bool {
        !self.is_surface_interaction()
    }
    fn spawn_ray(&self, d: &Vector3f) -> Ray {
        let o: Point3f =
            pnt3_offset_ray_origin(&self.common.p, &self.common.p_error, &self.common.n, d);
        Ray {
            o,
            d: *d,
            t_max: std::f32::INFINITY,
            time: self.common.time,
            differential: None,
            medium: self.get_medium(d),
        }
    }
    fn get_common(&self) -> Rc<InteractionCommon> {
        self.common.clone()
    }
    fn get_p(&self) -> &Point3f {
        &self.common.p
    }
    fn get_time(&self) -> Float {
        self.common.time
    }
    fn get_p_error(&self) -> &Vector3f {
        &self.common.p_error
    }
    fn get_wo(&self) -> &Vector3f {
        &self.common.wo
    }
    fn get_n(&self) -> &Normal3f {
        &self.common.n
    }
    fn get_medium_interface(&self) -> Option<Arc<MediumInterface>> {
        if let Some(ref medium_interface) = self.common.medium_interface {
            Some(medium_interface.clone())
        } else {
            None
        }
    }
    fn get_bsdf(&self) -> Option<&Bsdf> {
        None
    }
    fn get_shading_n(&self) -> Option<&Normal3f> {
        None
    }
    fn get_phase(&self) -> Option<Arc<HenyeyGreenstein>> {
        if let Some(ref phase) = self.phase {
            Some(phase.clone())
        } else {
            None
        }
    }
}

#[derive(Default)]
pub struct SurfaceInteraction<'a> {
    // Interaction Public Data
    pub common: Rc<InteractionCommon>,
    // SurfaceInteraction Public Data
    pub uv: Point2f,
    pub dpdu: Vector3f,
    pub dpdv: Vector3f,
    pub dndu: Normal3f,
    pub dndv: Normal3f,
    pub dpdx: Cell<Vector3f>,
    pub dpdy: Cell<Vector3f>,
    pub dudx: Cell<Float>,
    pub dvdx: Cell<Float>,
    pub dudy: Cell<Float>,
    pub dvdy: Cell<Float>,
    pub primitive: Option<*const Primitive>,
    pub shading: Shading,
    pub bsdf: Option<Bsdf>,
    pub bssrdf: Option<TabulatedBssrdf>,
    pub shape: Option<&'a Shape>,
}

impl<'a> SurfaceInteraction<'a> {
    pub fn new(
        p: &Point3f,
        p_error: &Vector3f,
        uv: Point2f,
        wo: &Vector3f,
        dpdu: &Vector3f,
        dpdv: &Vector3f,
        dndu: &Normal3f,
        dndv: &Normal3f,
        time: Float,
        sh: Option<&'a Shape>,
    ) -> Self {
        let nv: Vector3f = vec3_cross_vec3(dpdu, dpdv).normalize();
        let mut n: Normal3f = Normal3f {
            x: nv.x,
            y: nv.y,
            z: nv.z,
        };
        // initialize shading geometry from true geometry
        let mut shading: Shading = Shading {
            n,
            dpdu: *dpdu,
            dpdv: *dpdv,
            dndu: *dndu,
            dndv: *dndv,
        };
        if let Some(ref shape) = sh {
            // adjust normal based on orientation and handedness
            if shape.get_reverse_orientation() ^ shape.get_transform_swaps_handedness() {
                n *= -1.0 as Float;
                shading.n *= -1.0 as Float;
            }
        }
        let mut common: InteractionCommon = InteractionCommon::default();
        common.p = *p;
        common.time = time;
        common.p_error = *p_error;
        common.wo = wo.normalize();
        common.n = n;
        common.medium_interface = None;
        if let Some(ref shape) = sh {
            SurfaceInteraction {
                common: Rc::new(common),
                uv,
                dpdu: *dpdu,
                dpdv: *dpdv,
                dndu: *dndu,
                dndv: *dndv,
                dpdx: Cell::new(Vector3f::default()),
                dpdy: Cell::new(Vector3f::default()),
                dudx: Cell::new(0.0 as Float),
                dvdx: Cell::new(0.0 as Float),
                dudy: Cell::new(0.0 as Float),
                dvdy: Cell::new(0.0 as Float),
                primitive: None,
                shading,
                bsdf: None,
                bssrdf: None,
                shape: Some(shape),
            }
        } else {
            SurfaceInteraction {
                common: Rc::new(common),
                uv,
                dpdu: *dpdu,
                dpdv: *dpdv,
                dndu: *dndu,
                dndv: *dndv,
                dpdx: Cell::new(Vector3f::default()),
                dpdy: Cell::new(Vector3f::default()),
                dudx: Cell::new(0.0 as Float),
                dvdx: Cell::new(0.0 as Float),
                dudy: Cell::new(0.0 as Float),
                dvdy: Cell::new(0.0 as Float),
                primitive: None,
                shading,
                bsdf: None,
                bssrdf: None,
                shape: None,
            }
        }
    }
    pub fn get_medium(&self, w: &Vector3f) -> Option<Arc<Medium>> {
        if vec3_dot_nrm(w, &self.common.n) > 0.0 as Float {
            if let Some(ref medium_interface) = self.common.medium_interface {
                if let Some(ref outside_arc) = medium_interface.outside {
                    Some(outside_arc.clone())
                } else {
                    None
                }
            } else {
                None
            }
        } else if let Some(ref medium_interface) = self.common.medium_interface {
            if let Some(ref inside_arc) = medium_interface.inside {
                Some(inside_arc.clone())
            } else {
                None
            }
        } else {
            None
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
        self.shading.n = Normal3f::from(vec3_cross_vec3(dpdus, dpdvs)).normalize();
        if let Some(ref shape) = self.shape {
            if shape.get_reverse_orientation() ^ shape.get_transform_swaps_handedness() {
                self.shading.n = -self.shading.n;
            }
        }
        if orientation_is_authoritative {
            Rc::get_mut(&mut self.common).unwrap().n =
                nrm_faceforward_nrm(&self.common.n, &self.shading.n);
        } else {
            self.shading.n = nrm_faceforward_nrm(&self.shading.n, &self.common.n);
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
        if let Some(primitive_raw) = self.primitive {
            let primitive = unsafe { &*primitive_raw };
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
            let p = self.common.p;
            let d: Float = vec3_dot_vec3(
                &Vector3f::from(self.common.n),
                &Vector3f {
                    x: p.x,
                    y: p.y,
                    z: p.z,
                },
            );
            let tx: Float = -(vec3_dot_vec3(
                &Vector3f::from(self.common.n),
                &Vector3f::from(diff.rx_origin),
            ) - d)
                / vec3_dot_vec3(&Vector3f::from(self.common.n), &diff.rx_direction);
            if tx.is_infinite() || tx.is_nan() {
                self.dudx.set(0.0 as Float);
                self.dvdx.set(0.0 as Float);
                self.dudy.set(0.0 as Float);
                self.dvdy.set(0.0 as Float);
                self.dpdx.set(Vector3f::default());
                self.dpdy.set(Vector3f::default());
            } else {
                let px: Point3f = diff.rx_origin + diff.rx_direction * tx;
                let ty: Float = -(vec3_dot_vec3(
                    &Vector3f::from(self.common.n),
                    &Vector3f::from(diff.ry_origin),
                ) - d)
                    / vec3_dot_vec3(&Vector3f::from(self.common.n), &diff.ry_direction);
                if ty.is_infinite() || ty.is_nan() {
                    self.dudx.set(0.0 as Float);
                    self.dvdx.set(0.0 as Float);
                    self.dudy.set(0.0 as Float);
                    self.dvdy.set(0.0 as Float);
                    self.dpdx.set(Vector3f::default());
                    self.dpdy.set(Vector3f::default());
                } else {
                    let py: Point3f = diff.ry_origin + diff.ry_direction * ty;
                    self.dpdx.set(px - self.common.p);
                    self.dpdy.set(py - self.common.p);

                    // compute $(u,v)$ offsets at auxiliary points

                    // choose two dimensions to use for ray offset computation
                    let mut dim: [XYZEnum; 2] = [XYZEnum::X; 2];
                    if self.common.n.x.abs() > self.common.n.y.abs()
                        && self.common.n.x.abs() > self.common.n.z.abs()
                    {
                        dim[0] = XYZEnum::Y;
                        dim[1] = XYZEnum::Z;
                    } else if self.common.n.y.abs() > self.common.n.z.abs() {
                        dim[0] = XYZEnum::X;
                        dim[1] = XYZEnum::Z;
                    } else {
                        dim[0] = XYZEnum::X;
                        dim[1] = XYZEnum::Y;
                    }

                    // initialize _a_, _bx_, and _by_ matrices for offset computation
                    let a0: [Float; 2] = [self.dpdu[dim[0]], self.dpdv[dim[0]]];
                    let a1: [Float; 2] = [self.dpdu[dim[1]], self.dpdv[dim[1]]];
                    let a: [[Float; 2]; 2] = [a0, a1];
                    let bx: [Float; 2] = [
                        px[dim[0]] - self.common.p[dim[0]],
                        px[dim[1]] - self.common.p[dim[1]],
                    ];
                    let by: [Float; 2] = [
                        py[dim[0]] - self.common.p[dim[0]],
                        py[dim[1]] - self.common.p[dim[1]],
                    ];
                    if !solve_linear_system_2x2(a, bx, self.dudx.get_mut(), self.dvdx.get_mut()) {
                        self.dudx.set(0.0 as Float);
                        self.dvdx.set(0.0 as Float);
                    }
                    if !solve_linear_system_2x2(a, by, self.dudy.get_mut(), self.dvdy.get_mut()) {
                        self.dudy.set(0.0 as Float);
                        self.dvdy.set(0.0 as Float);
                    }
                }
            }
        } else {
            self.dudx.set(0.0 as Float);
            self.dvdx.set(0.0 as Float);
            self.dudy.set(0.0 as Float);
            self.dvdy.set(0.0 as Float);
            self.dpdx.set(Vector3f::default());
            self.dpdy.set(Vector3f::default());
        }
    }
    pub fn le(&self, w: &Vector3f) -> Spectrum {
        if let Some(primitive_raw) = self.primitive {
            let primitive = unsafe { &*primitive_raw };
            if let Some(area_light) = primitive.get_area_light() {
                // create InteractionCommon from self
                let interaction: InteractionCommon = InteractionCommon {
                    p: self.common.p,
                    time: self.common.time,
                    p_error: self.common.p_error,
                    wo: self.common.wo,
                    n: self.common.n,
                    medium_interface: None,
                };
                return area_light.l(&interaction, w);
            }
        }
        Spectrum::default()
    }
}

impl<'a> Interaction for SurfaceInteraction<'a> {
    fn is_surface_interaction(&self) -> bool {
        self.common.n != Normal3f::default()
    }
    fn is_medium_interaction(&self) -> bool {
        !self.is_surface_interaction()
    }
    fn spawn_ray(&self, d: &Vector3f) -> Ray {
        let o: Point3f =
            pnt3_offset_ray_origin(&self.common.p, &self.common.p_error, &self.common.n, d);
        Ray {
            o,
            d: *d,
            t_max: std::f32::INFINITY,
            time: self.common.time,
            differential: None,
            medium: self.get_medium(d),
        }
    }
    fn get_common(&self) -> Rc<InteractionCommon> {
        self.common.clone()
    }
    fn get_p(&self) -> &Point3f {
        &self.common.p
    }
    fn get_time(&self) -> Float {
        self.common.time
    }
    fn get_p_error(&self) -> &Vector3f {
        &self.common.p_error
    }
    fn get_wo(&self) -> &Vector3f {
        &self.common.wo
    }
    fn get_n(&self) -> &Normal3f {
        &self.common.n
    }
    fn get_medium_interface(&self) -> Option<Arc<MediumInterface>> {
        if let Some(ref medium_interface) = self.common.medium_interface {
            Some(medium_interface.clone())
        } else {
            None
        }
    }
    fn get_bsdf(&self) -> Option<&Bsdf> {
        if let Some(ref bsdf) = self.bsdf {
            Some(bsdf)
        } else {
            None
        }
    }
    fn get_shading_n(&self) -> Option<&Normal3f> {
        Some(&self.shading.n)
    }
    fn get_phase(&self) -> Option<Arc<HenyeyGreenstein>> {
        None
    }
}
