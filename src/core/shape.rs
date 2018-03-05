//! Careful abstraction of geometric shapes in a ray tracer is a key
//! component of a clean system design, and shapes are the ideal
//! candidate for an object-oriented approach. All geometric
//! primitives implement a common interface, and the rest of the
//! renderer can use this interface without needing any details about
//! the underlying shape. This makes it possible to separate the
//! geometric and the shading subsystem of pbrt.

// pbrt
use core::geometry::{Bounds3f, Point2f, Ray, Vector3f};
use core::pbrt::Float;
use core::interaction::{Interaction, InteractionCommon, SurfaceInteraction};

// see shape.h

pub trait Shape {
    fn object_bound(&self) -> Bounds3f;
    fn world_bound(&self) -> Bounds3f;
    fn intersect(&self, r: &Ray) -> Option<(SurfaceInteraction, Float)>;
    fn intersect_p(&self, r: &Ray) -> bool;
    fn get_reverse_orientation(&self) -> bool;
    fn get_transform_swaps_handedness(&self) -> bool;
    fn area(&self) -> Float;
    fn sample(&self, u: &Point2f, pdf: &mut Float) -> InteractionCommon;
    fn pdf(&self, _iref: &InteractionCommon) -> Float {
        1.0 as Float / self.area()
    }
    fn sample_with_ref_point(
        &self,
        iref: &InteractionCommon,
        u: &Point2f,
        pdf: &mut Float,
    ) -> InteractionCommon;
    fn pdf_with_ref_point(&self, iref: &Interaction, wi: &Vector3f) -> Float;
}
