// std
use std::sync::Arc;
// pbrt
use core::geometry::{Bounds3f, Normal3f, Point2f, Point3f, Ray, Vector2f, Vector3f};
use core::geometry::{nrm_dot_nrm, nrm_normalize, bnd3_expand, bnd3_union_bnd3, nrm_abs_dot_vec3,
                     nrm_cross_vec3, pnt3_distance, pnt3_distance_squared, pnt3_lerp, vec2_dot,
                     vec3_coordinate_system, vec3_cross_vec3, vec3_normalize};
use core::interaction::{Interaction, InteractionCommon, SurfaceInteraction};
use core::material::Material;
use core::paramset::ParamSet;
use core::pbrt::Float;
use core::pbrt::{clamp_t, float_to_bits, lerp};
use core::shape::Shape;
use core::transform::Transform;

// see loopsubdiv.cpp

pub fn create_loop_subdiv(
    o2w: &Transform,
    w2o: &Transform,
    reverse_orientation: bool,
    params: &ParamSet,
) -> Vec<Arc<Shape + Send + Sync>> {
    // WORK
    Vec::new()
}
