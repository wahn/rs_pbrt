//! All filter implementations are derived from an abstract **Filter**
//! class, which provides the interface for the functions used in
//! filtering.

// pbrt
use crate::core::geometry::{Point2f, Vector2f};
use crate::core::pbrt::Float;
use crate::filters::boxfilter::BoxFilter;
use crate::filters::gaussian::GaussianFilter;
use crate::filters::mitchell::MitchellNetravali;
use crate::filters::sinc::LanczosSincFilter;
use crate::filters::triangle::TriangleFilter;

// see filter.h

pub enum Filter {
    Bx(BoxFilter),
    Gaussian(GaussianFilter),
    MitchellNetravali(MitchellNetravali),
    LanczosSinc(LanczosSincFilter),
    Triangle(TriangleFilter),
}

impl Filter {
    pub fn evaluate(&self, p: Point2f) -> Float {
        match self {
            Filter::Bx(filter) => filter.evaluate(p),
            Filter::Gaussian(filter) => filter.evaluate(p),
            Filter::MitchellNetravali(filter) => filter.evaluate(p),
            Filter::LanczosSinc(filter) => filter.evaluate(p),
            Filter::Triangle(filter) => filter.evaluate(p),
        }
    }
    pub fn get_radius(&self) -> Vector2f {
        match self {
            Filter::Bx(filter) => filter.get_radius(),
            Filter::Gaussian(filter) => filter.get_radius(),
            Filter::MitchellNetravali(filter) => filter.get_radius(),
            Filter::LanczosSinc(filter) => filter.get_radius(),
            Filter::Triangle(filter) => filter.get_radius(),
        }
    }
}
