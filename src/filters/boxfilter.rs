// pbrt
use core::pbrt::Float;
use filters::Filter;
use geometry::{Point2f, Vector2f};

// see box.h

#[derive(Debug,Default,Copy,Clone)]
pub struct BoxFilter {
    // inherited from Filter (see filter.h)
    pub radius: Vector2f,
    pub inv_radius: Vector2f,
}

impl Filter for BoxFilter {
    fn evaluate(&self, _p: Point2f) -> Float {
        1.0
    }
    fn get_radius(&self) -> Vector2f {
        Vector2f {
            x: self.radius.x,
            y: self.radius.y,
        }
    }
}
