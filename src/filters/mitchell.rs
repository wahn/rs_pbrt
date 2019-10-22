use crate::core::filter::Filter;
use crate::core::geometry::{Point2f, Vector2f};
use crate::core::paramset::ParamSet;
use crate::core::pbrt::Float;

#[derive(Debug, Default, Copy, Clone)]
pub struct MitchellNetravali {
    width: Float,
    height: Float,
    inv_width: Float,
    inv_height: Float,
    b: Float,
    c: Float,
}

impl MitchellNetravali {
    pub fn new(w: Float, h: Float, b: Float, c: Float) -> Self {
        MitchellNetravali {
            width: w,
            height: h,
            inv_width: 1.0 / w,
            inv_height: 1.0 / h,
            b,
            c,
        }
    }
    pub fn mitchell_1d(&self, x: Float) -> Float {
        let fx = x.abs() * 2.0;
        if fx < 1.0 {
            ((12.0 - 9.0 * self.b - 6.0 * self.c) * fx * fx * fx
                + (-18.0 + 12.0 * self.b + 6.0 * self.c) * fx * fx
                + (6.0 - 2.0 * self.b))
                * (1.0 / 6.0)
        } else if fx < 2.0 {
            ((-self.b - 6.0 * self.c) * fx * fx * fx
                + (6.0 * self.b + 30.0 * self.c) * fx * fx
                + (-12.0 * self.b - 48.0 * self.c) * fx
                + (8.0 * self.b + 24.0 * self.c))
                * (1.0 / 6.0)
        } else {
            0.0
        }
    }
    pub fn create(ps: &ParamSet) -> Box<Filter> {
        let xw = ps.find_one_float("xwidth", 2.0);
        let yw = ps.find_one_float("ywidth", 2.0);
        let b = ps.find_one_float("B", 1.0 / 3.0);
        let c = ps.find_one_float("C", 1.0 / 3.0);

        Box::new(Filter::MitchellNetravali(Self::new(xw, yw, b, c)))
    }
    // Filter
    pub fn evaluate(&self, p: Point2f) -> Float {
        self.mitchell_1d(p.x * self.inv_width) * self.mitchell_1d(p.y * self.inv_height)
    }
    pub fn get_radius(&self) -> Vector2f {
        Vector2f {
            x: self.width,
            y: self.height,
        }
    }
}
