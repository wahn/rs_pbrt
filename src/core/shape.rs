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
    fn sample(&self, u: Point2f, pdf: &mut Float) -> InteractionCommon;
    fn sample_with_ref_point(&self,
                             iref: &InteractionCommon,
                             u: Point2f,
                             pdf: &mut Float)
                             -> InteractionCommon;
    fn pdf(&self, iref: &Interaction, wi: Vector3f) -> Float;
}
