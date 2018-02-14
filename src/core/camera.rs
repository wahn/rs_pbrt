//! The abstract **Camera** base class holds generic camera options
//! and defines the interface that all camera implementations must
//! provide.

// std
use std::sync::Arc;
// pbrt
use core::film::Film;
use core::geometry::{Point2f, Ray};
use core::pbrt::Float;

// see camera.h

pub trait Camera {
    fn generate_ray_differential(&self, sample: &CameraSample, ray: &mut Ray) -> Float;
    fn pdf_we(&self, ray: &Ray) -> (Float, Float);
    fn get_film(&self) -> Arc<Film>;
}

#[derive(Debug,Default,Copy,Clone)]
pub struct CameraSample {
    pub p_film: Point2f,
    pub p_lens: Point2f,
    pub time: Float,
}
