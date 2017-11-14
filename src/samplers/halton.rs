// pbrt
use core::camera::CameraSample;
use core::pbrt::Float;
use core::sampler::Sampler;
use geometry::{Point2f, Point2i};

// see halton.h

#[derive(Debug,Clone)]
pub struct HaltonSampler {
    pub samples_per_pixel: i64,
    pub current_pixel_sample_index: i64,
}

impl Sampler for HaltonSampler {
    fn start_pixel(&mut self, _p: Point2i) {
        // WORK
    }
    fn get_1d(&mut self) -> Float {
        // WORK
        0.0 as Float
    }
    fn get_2d(&mut self) -> Point2f {
        // WORK
        Point2f {
            x: 0.0 as Float,
            y: 0.0 as Float,
        }
    }
    fn request_2d_array(&mut self, _n: i32) {
        // WORK
    }
    fn round_count(&self, _count: i32) -> i32 {
        // WORK
        0_i32
    }
    fn get_2d_array(&mut self, _n: i32) -> Vec<Point2f> {
        let mut samples: Vec<Point2f> = Vec::new();
        // WORK
        samples
    }
    fn start_next_sample(&mut self) -> bool {
        // WORK
        false
    }
    fn get_camera_sample(&mut self, _p_raster: Point2i) -> CameraSample {
        let mut cs: CameraSample = CameraSample::default();
        // WORK
        cs
    }
    fn reseed(&mut self, _seed: u64) {
        // WORK
    }
    fn box_clone(&self) -> Box<Sampler + Send + Sync> {                                                 
        Box::new(self.clone())                                                                          
    }                                                                                                   
    fn get_current_sample_number(&self) -> i64 {
        self.current_pixel_sample_index
    }
    fn get_samples_per_pixel(&self) -> i64 {
        self.samples_per_pixel
    }
}
