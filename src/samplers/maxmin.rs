// std
use std::sync::RwLock;
// pbrt
use crate::core::geometry::{Bounds2i, Point2f, Point2i, Vector2i};
use crate::core::pbrt::Float;
use crate::core::rng::Rng;
use crate::core::sampler::Sampler;

pub struct MaxMinDistSampler {
    pub c_pixel: u32,
    // inherited from class PixelSampler (see sampler.h)
    pub samples_1d: Vec<Vec<Float>>,
    pub samples_2d: Vec<Vec<Point2f>>,
    pub current_1d_dimension: i32,
    pub current_2d_dimension: i32,
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

impl MaxMinDistSampler {
    pub fn clone_with_seed(&self, seed: u64) -> Box<Sampler> {
        let mut mmds = MaxMinDistSampler {
            c_pixel: self.c_pixel,
            samples_1d: self.samples_1d.clone(),
            samples_2d: self.samples_2d.clone(),
            current_1d_dimension: self.current_1d_dimension,
            current_2d_dimension: self.current_2d_dimension,
            rng: self.rng.clone(),
            current_pixel: self.current_pixel,
            current_pixel_sample_index: self.current_pixel_sample_index,
            samples_1d_array_sizes: self.samples_1d_array_sizes.iter().cloned().collect(),
            samples_2d_array_sizes: self.samples_2d_array_sizes.iter().cloned().collect(),
            sample_array_1d: self.sample_array_1d.iter().cloned().collect(),
            sample_array_2d: self.sample_array_2d.iter().cloned().collect(),
            array_1d_offset: self.array_1d_offset,
            array_2d_offset: self.array_2d_offset,
        };
        mmds.reseed(seed);
        let sampler = Sampler::MaxMinDist(mmds);
        Box::new(sampler)
    }
    pub fn start_pixel(&mut self, p: &Point2i) {
        // WORK
    }
    pub fn get_1d(&mut self) -> Float {
        // WORK
        0.0 as Float
    }
    pub fn get_2d(&mut self) -> Point2f {
        // WORK
        Point2f::default()
    }
    pub fn request_2d_array(&mut self, n: i32) {
        // WORK
    }
    pub fn round_count(&self, count: i32) -> i32 {
        // WORK
        0_i32
    }
    pub fn get_2d_array(&mut self, n: i32) -> Option<&[Point2f]> {
        // WORK
        None
    }
    pub fn get_2d_arrays(&mut self, n: i32) -> (Option<&[Point2f]>, Option<&[Point2f]>) {
        // WORK
        (None, None)
    }
    pub fn get_2d_array_vec(&mut self, n: i32) -> Vec<Point2f> {
        // WORK
        Vec::new()
    }
    pub fn start_next_sample(&mut self) -> bool {
        // WORK
        false
    }
    pub fn reseed(&mut self, seed: u64) {
        self.rng.set_sequence(seed);
    }
    pub fn get_current_pixel(&self) -> Point2i {
        // WORK
        Point2i::default()
    }
    pub fn get_current_sample_number(&self) -> i64 {
        // WORK
        0_i64
    }
    pub fn get_samples_per_pixel(&self) -> i64 {
        // WORK
        0_i64
    }
}
