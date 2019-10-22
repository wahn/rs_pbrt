//! Careful abstraction of geometric shapes in a ray tracer is a key
//! component of a clean system design, and shapes are the ideal
//! candidate for an object-oriented approach. All geometric
//! primitives implement a common interface, and the rest of the
//! renderer can use this interface without needing any details about
//! the underlying shape. This makes it possible to separate the
//! geometric and the shading subsystem of pbrt.

// pbrt
use crate::core::geometry::{Bounds3f, Point2f, Ray, Vector3f};
use crate::core::interaction::{Interaction, InteractionCommon, SurfaceInteraction};
use crate::core::pbrt::Float;
use crate::shapes::curve::Curve;
use crate::shapes::cylinder::Cylinder;
use crate::shapes::disk::Disk;
use crate::shapes::sphere::Sphere;
use crate::shapes::triangle::Triangle;

// see shape.h

pub enum Shape {
    Crv(Curve),
    Clndr(Cylinder),
    Dsk(Disk),
    Sphr(Sphere),
    Trngl(Triangle),
}

impl Shape {
    pub fn object_bound(&self) -> Bounds3f {
        match self {
            Shape::Crv(shape) => shape.object_bound(),
            Shape::Clndr(shape) => shape.object_bound(),
            Shape::Dsk(shape) => shape.object_bound(),
            Shape::Sphr(shape) => shape.object_bound(),
            Shape::Trngl(shape) => shape.object_bound(),
        }
    }
    pub fn world_bound(&self) -> Bounds3f {
        match self {
            Shape::Crv(shape) => shape.world_bound(),
            Shape::Clndr(shape) => shape.world_bound(),
            Shape::Dsk(shape) => shape.world_bound(),
            Shape::Sphr(shape) => shape.world_bound(),
            Shape::Trngl(shape) => shape.world_bound(),
        }
    }
    pub fn intersect(&self, r: &Ray) -> Option<(SurfaceInteraction, Float)> {
        match self {
            Shape::Crv(shape) => shape.intersect(r),
            Shape::Clndr(shape) => shape.intersect(r),
            Shape::Dsk(shape) => shape.intersect(r),
            Shape::Sphr(shape) => shape.intersect(r),
            Shape::Trngl(shape) => shape.intersect(r),
        }
    }
    pub fn intersect_p(&self, r: &Ray) -> bool {
        match self {
            Shape::Crv(shape) => shape.intersect_p(r),
            Shape::Clndr(shape) => shape.intersect_p(r),
            Shape::Dsk(shape) => shape.intersect_p(r),
            Shape::Sphr(shape) => shape.intersect_p(r),
            Shape::Trngl(shape) => shape.intersect_p(r),
        }
    }
    pub fn get_reverse_orientation(&self) -> bool {
        match self {
            Shape::Crv(shape) => shape.get_reverse_orientation(),
            Shape::Clndr(shape) => shape.get_reverse_orientation(),
            Shape::Dsk(shape) => shape.get_reverse_orientation(),
            Shape::Sphr(shape) => shape.get_reverse_orientation(),
            Shape::Trngl(shape) => shape.get_reverse_orientation(),
        }
    }
    pub fn get_transform_swaps_handedness(&self) -> bool {
        match self {
            Shape::Crv(shape) => shape.get_transform_swaps_handedness(),
            Shape::Clndr(shape) => shape.get_transform_swaps_handedness(),
            Shape::Dsk(shape) => shape.get_transform_swaps_handedness(),
            Shape::Sphr(shape) => shape.get_transform_swaps_handedness(),
            Shape::Trngl(shape) => shape.get_transform_swaps_handedness(),
        }
    }
    pub fn area(&self) -> Float {
        match self {
            Shape::Crv(shape) => shape.area(),
            Shape::Clndr(shape) => shape.area(),
            Shape::Dsk(shape) => shape.area(),
            Shape::Sphr(shape) => shape.area(),
            Shape::Trngl(shape) => shape.area(),
        }
    }
    pub fn sample(&self, u: &Point2f, pdf: &mut Float) -> InteractionCommon {
        match self {
            Shape::Crv(shape) => shape.sample(u, pdf),
            Shape::Clndr(shape) => shape.sample(u, pdf),
            Shape::Dsk(shape) => shape.sample(u, pdf),
            Shape::Sphr(shape) => shape.sample(u, pdf),
            Shape::Trngl(shape) => shape.sample(u, pdf),
        }
    }
    pub fn pdf(&self, _iref: &InteractionCommon) -> Float {
        1.0 as Float / self.area()
    }
    pub fn sample_with_ref_point(
        &self,
        iref: &InteractionCommon,
        u: &Point2f,
        pdf: &mut Float,
    ) -> InteractionCommon {
        match self {
            Shape::Crv(shape) => shape.sample_with_ref_point(iref, u, pdf),
            Shape::Clndr(shape) => shape.sample_with_ref_point(iref, u, pdf),
            Shape::Dsk(shape) => shape.sample_with_ref_point(iref, u, pdf),
            Shape::Sphr(shape) => shape.sample_with_ref_point(iref, u, pdf),
            Shape::Trngl(shape) => shape.sample_with_ref_point(iref, u, pdf),
        }
    }
    pub fn pdf_with_ref_point(&self, iref: &dyn Interaction, wi: &Vector3f) -> Float {
        match self {
            Shape::Crv(shape) => shape.pdf_with_ref_point(iref, wi),
            Shape::Clndr(shape) => shape.pdf_with_ref_point(iref, wi),
            Shape::Dsk(shape) => shape.pdf_with_ref_point(iref, wi),
            Shape::Sphr(shape) => shape.pdf_with_ref_point(iref, wi),
            Shape::Trngl(shape) => shape.pdf_with_ref_point(iref, wi),
        }
    }
}
