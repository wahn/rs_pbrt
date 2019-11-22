//! The **Sampler** base class not only defines the interface to
//! samplers but also provides some common functionality for use by
//! **Sampler** implementations.

// pbrt
use crate::core::camera::CameraSample;
use crate::core::geometry::{Point2f, Point2i};
use crate::core::pbrt::Float;
use crate::samplers::random::RandomSampler;

// see sampler.h

pub enum Sampler {
    Random(RandomSampler),
}

impl Sampler {
    pub fn start_pixel(&mut self, p: &Point2i) {}
    pub fn get_1d(&mut self) -> Float {
        0.0 as Float
    }
    pub fn get_2d(&mut self) -> Point2f {
        Point2f::default()
    }
    pub fn get_camera_sample(&mut self, p_raster: &Point2i) -> CameraSample {
        let mut cs: CameraSample = CameraSample::default();
        cs.p_film = Point2f {
            x: p_raster.x as Float,
            y: p_raster.y as Float,
        } + self.get_2d();
        cs.time = self.get_1d();
        cs.p_lens = self.get_2d();
        cs
    }
    pub fn request_2d_array(&mut self, n: i32) {}
    pub fn round_count(&self, count: i32) -> i32 {
        0
    }
    pub fn get_2d_array(&mut self, n: i32) -> Option<&[Point2f]> {
        None
    }
    pub fn get_2d_arrays(&mut self, n: i32) -> (Option<&[Point2f]>, Option<&[Point2f]>) {
        (None, None)
    }
    pub fn get_2d_array_vec(&mut self, n: i32) -> Vec<Point2f> {
        Vec::new()
    }
    pub fn start_next_sample(&mut self) -> bool {
        false
    }
    pub fn reseed(&mut self, seed: u64) {}
    pub fn get_current_pixel(&self) -> Point2i {
        Point2i::default()
    }
    pub fn get_current_sample_number(&self) -> i64 {
        0
    }
    pub fn get_samples_per_pixel(&self) -> i64 {
        0
    }
}
