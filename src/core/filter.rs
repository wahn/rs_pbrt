//! All filter implementations are derived from an abstract **Filter**
//! class, which provides the interface for the functions used in
//! filtering.

// pbrt
use core::geometry::{Point2f, Vector2f};
use core::pbrt::Float;

// see filter.h

pub trait Filter {
    fn evaluate(&self, p: Point2f) -> Float;
    fn get_radius(&self) -> Vector2f;
}
