// pbrt
use core::geometry::{Point2f, Vector2f};
use core::pbrt::Float;

// see filter.h

pub trait Filter {
    fn evaluate(&self, p: Point2f) -> Float;
    fn get_radius(&self) -> Vector2f;
}
