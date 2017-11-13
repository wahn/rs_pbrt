// pbrt
use core::pbrt::Float;
use geometry::Point2f;

// see sampler.h

pub trait Sampler {
    fn get_1d(&mut self) -> Float;
    fn get_2d(&mut self) -> Point2f;
}
