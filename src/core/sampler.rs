//! The **Sampler** base class not only defines the interface to
//! samplers but also provides some common functionality for use by
//! **Sampler** implementations.

// pbrt
use core::camera::CameraSample;
use core::geometry::{Point2f, Point2i};
use core::pbrt::Float;

// see sampler.h

pub trait Sampler: SamplerClone {
    fn start_pixel(&mut self, p: &Point2i);
    fn get_1d(&mut self) -> Float;
    fn get_2d(&mut self) -> Point2f;
    fn get_camera_sample(&mut self, p_raster: &Point2i) -> CameraSample {
        let mut cs: CameraSample = CameraSample::default();
        cs.p_film = Point2f {
            x: p_raster.x as Float,
            y: p_raster.y as Float,
        } + self.get_2d();
        cs.time = self.get_1d();
        cs.p_lens = self.get_2d();
        cs
    }
    fn request_2d_array(&mut self, n: i32);
    fn round_count(&self, count: i32) -> i32;
    fn get_2d_array(&mut self, n: i32) -> Vec<Point2f>;
    fn start_next_sample(&mut self) -> bool;
    fn reseed(&mut self, seed: u64);
    fn get_current_pixel(&self) -> Point2i;
    fn get_current_sample_number(&self) -> i64;
    fn get_samples_per_pixel(&self) -> i64;
}

pub trait PixelSampler: Sampler {}

pub trait GlobalSampler: Sampler {}

pub trait SamplerClone {
    fn box_clone(&self) -> Box<Sampler + Send + Sync>;
}

impl<T> SamplerClone for T
where
    T: 'static + Sampler + Clone + Send + Sync,
{
    fn box_clone(&self) -> Box<Sampler + Send + Sync> {
        Box::new(self.clone())
    }
}

impl Clone for Box<Sampler + Send + Sync> {
    fn clone(&self) -> Box<Sampler + Send + Sync> {
        self.box_clone()
    }
}
