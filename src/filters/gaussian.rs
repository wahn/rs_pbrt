//std
use std::sync::Arc;
// pbrt
use crate::core::filter::Filter;
use crate::core::geometry::{Point2f, Vector2f};
use crate::core::paramset::ParamSet;
use crate::core::pbrt::Float;

// see gaussian.h

#[derive(Debug, Default, Copy, Clone)]
pub struct GaussianFilter {
    pub alpha: Float,
    pub exp_x: Float,
    pub exp_y: Float,
    // inherited from Filter (see filter.h)
    pub radius: Vector2f,
    pub inv_radius: Vector2f,
}

impl GaussianFilter {
    pub fn create(ps: &ParamSet) -> Arc<Filter + Sync + Send> {
        let xw: Float = ps.find_one_float("xwidth", 2.0);
        let yw: Float = ps.find_one_float("ywidth", 2.0);
        let alpha: Float = ps.find_one_float("alpha", 2.0);
        // see gaussian.h (GaussianFilter constructor)
        let exp_x: Float = (-alpha * xw * xw).exp();
        let exp_y: Float = (-alpha * yw * yw).exp();
        let gaussian_filter: Arc<Filter + Sync + Send> = Arc::new(GaussianFilter {
            alpha: alpha,
            exp_x: exp_x,
            exp_y: exp_y,
            radius: Vector2f { x: xw, y: yw },
            inv_radius: Vector2f {
                x: 1.0 / xw,
                y: 1.0 / yw,
            },
        });
        gaussian_filter
    }
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
