use crate::core::filter::Filter;
use crate::core::geometry::{Point2f, Vector2f};
use crate::core::paramset::ParamSet;
use crate::core::pbrt::Float;

pub struct LanczosSincFilter {}

impl LanczosSincFilter {}

impl Filter for LanczosSincFilter {
    fn evaluate(&self, p: Point2f) -> Float {
        // WORK
        0.0 as Float
    }

    fn get_radius(&self) -> Vector2f {
        // WORK
        Vector2f::default()
    }
}
