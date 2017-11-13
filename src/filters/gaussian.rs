// pbrt
use core::pbrt::Float;
use filters::Filter;
use geometry::{Point2f, Vector2f};

// see gaussian.h

#[derive(Debug,Default,Copy,Clone)]
pub struct GaussianFilter {
    pub alpha: Float,
    pub exp_x: Float,
    pub exp_y: Float,
    // inherited from Filter (see filter.h)
    pub radius: Vector2f,
    pub inv_radius: Vector2f,
}

impl GaussianFilter {
    pub fn gaussian(&self, d: Float, expv: Float) -> Float {
        (0.0 as Float).max((-self.alpha * d * d).exp() - expv)
    }
}

impl Filter for GaussianFilter {
    fn evaluate(&self, p: Point2f) -> Float {
        self.gaussian(p.x, self.exp_x) * self.gaussian(p.y, self.exp_y)
    }
    fn get_radius(&self) -> Vector2f {
        Vector2f {
            x: self.radius.x,
            y: self.radius.y,
        }
    }
}
