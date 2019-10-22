// std
use std::f32::consts::PI;
// pbrt
use crate::core::filter::Filter;
use crate::core::geometry::{Point2f, Vector2f};
use crate::core::paramset::ParamSet;
use crate::core::pbrt::Float;

#[derive(Debug, Default, Copy, Clone)]
pub struct LanczosSincFilter {
    tau: Float,
    // inherited from Filter (see filter.h)
    pub radius: Vector2f,
    pub inv_radius: Vector2f,
}

impl LanczosSincFilter {
    pub fn new(radius: &Vector2f, tau: Float) -> Self {
        LanczosSincFilter {
            tau,
            radius: *radius,
            inv_radius: Vector2f {
                x: 1.0 / radius.x,
                y: 1.0 / radius.y,
            },
        }
    }
    pub fn create(ps: &ParamSet) -> Box<Filter> {
        let xw: Float = ps.find_one_float("xwidth", 4.0);
        let yw: Float = ps.find_one_float("ywidth", 4.0);
        let tau: Float = ps.find_one_float("tau", 3.0);
        let sinc_filter: Box<Filter> = Box::new(Filter::LanczosSinc(LanczosSincFilter::new(
            &Vector2f { x: xw, y: yw },
            tau,
        )));
        sinc_filter
    }
    pub fn sinc(&self, x: Float) -> Float {
        let x = x.abs();
        if x < 1e-5 as Float {
            return 1.0 as Float;
        }
        (PI * x).sin() / (PI * x)
    }
    pub fn windowed_sinc(&self, x: Float, radius: Float) -> Float {
        let x = x.abs();
        if x > radius {
            return 0.0 as Float;
        }
        let lanczos: Float = self.sinc(x / self.tau);
        self.sinc(x) * lanczos
    }
    // Filter
    pub fn evaluate(&self, p: Point2f) -> Float {
        self.windowed_sinc(p.x, self.radius.x) * self.windowed_sinc(p.y, self.radius.y)
    }
    pub fn get_radius(&self) -> Vector2f {
        Vector2f {
            x: self.radius.x,
            y: self.radius.y,
        }
    }
}
