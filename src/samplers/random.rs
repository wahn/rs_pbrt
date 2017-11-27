// pbrt
use core::camera::CameraSample;
use core::pbrt::Float;
use core::rng::Rng;
use core::sampler::Sampler;
use geometry::{Point2f, Point2i};

// see random.h

pub struct RandomSampler {
    pub samples_per_pixel: i64,
    pub rng: Rng,
    // inherited from class Sampler (see sampler.h)
    pub current_pixel: Point2i,
    pub current_pixel_sample_index: i64,
    pub samples_1d_array_sizes: Vec<i32>,
    pub samples_2d_array_sizes: Vec<i32>,
    pub sample_array_1d: Vec<Vec<Float>>,
    pub sample_array_2d: Vec<Vec<Point2f>>,
    pub array_1d_offset: usize,
    pub array_2d_offset: usize,
}

impl RandomSampler {
}

impl Sampler for RandomSampler {
    fn start_pixel(&mut self, p: Point2i) {
        // WORK
    }
    fn get_1d(&mut self) -> Float {
        // WORK
        0.0 as Float
    }
    fn get_2d(&mut self) -> Point2f {
        // WORK
        Point2f::default()
    }
    fn request_2d_array(&mut self, n: i32) {
        // WORK
    }
    fn round_count(&self, count: i32) -> i32 {
        // WORK
        0_i32
    }
    fn get_2d_array(&mut self, n: i32) -> Vec<Point2f> {
        // WORK
        Vec::new()
    }
    fn start_next_sample(&mut self) -> bool {
        // WORK
        false
    }
    fn get_camera_sample(&mut self, p_raster: Point2i) -> CameraSample {
        // WORK
        CameraSample::default()
    }                                         
    fn reseed(&mut self, _seed: u64) {
        // do nothing
    }
    fn get_current_sample_number(&self) -> i64 {
        self.current_pixel_sample_index
    }
    fn get_samples_per_pixel(&self) -> i64 {
        self.samples_per_pixel
    }
}

impl Clone for RandomSampler {
    fn clone(&self) -> RandomSampler {
        RandomSampler {
            samples_per_pixel: self.samples_per_pixel,
            rng: Rng::default(),
            current_pixel: self.current_pixel,
            current_pixel_sample_index: self.current_pixel_sample_index,
            samples_1d_array_sizes: self.samples_1d_array_sizes.iter().cloned().collect(),
            samples_2d_array_sizes: self.samples_2d_array_sizes.iter().cloned().collect(),
            sample_array_1d: self.sample_array_1d.iter().cloned().collect(),
            sample_array_2d: self.sample_array_2d.iter().cloned().collect(),
            array_1d_offset: self.array_1d_offset,
            array_2d_offset: self.array_2d_offset,
        }
    }
}
