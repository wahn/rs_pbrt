//! The **Sampler** base class not only defines the interface to
//! samplers but also provides some common functionality for use by
//! **Sampler** implementations.

// pbrt
use crate::core::camera::CameraSample;
use crate::core::geometry::{Point2f, Point2i};
use crate::core::pbrt::Float;
use crate::integrators::mlt::MLTSampler;
use crate::samplers::halton::HaltonSampler;
use crate::samplers::random::RandomSampler;
use crate::samplers::sobol::SobolSampler;
use crate::samplers::zerotwosequence::ZeroTwoSequenceSampler;

// see sampler.h

pub enum Sampler {
    Halton(HaltonSampler),
    MLT(MLTSampler),
    Random(RandomSampler),
    Sobol(SobolSampler),
    ZeroTwoSequence(ZeroTwoSequenceSampler),
}

impl Sampler {
    pub fn clone_with_seed(&self, seed: u64) -> Box<Sampler> {
        match self {
            Sampler::Halton(sampler) => sampler.clone_with_seed(seed),
            Sampler::MLT(sampler) => sampler.clone_with_seed(seed),
            Sampler::Random(sampler) => sampler.clone_with_seed(seed),
            Sampler::Sobol(sampler) => sampler.clone_with_seed(seed),
            Sampler::ZeroTwoSequence(sampler) => sampler.clone_with_seed(seed),
        }
    }
    pub fn start_pixel(&mut self, p: &Point2i) {
        match self {
            Sampler::Halton(sampler) => sampler.start_pixel(p),
            Sampler::MLT(sampler) => sampler.start_pixel(p),
            Sampler::Random(sampler) => sampler.start_pixel(p),
            Sampler::Sobol(sampler) => sampler.start_pixel(p),
            Sampler::ZeroTwoSequence(sampler) => sampler.start_pixel(p),
        }
    }
    pub fn get_1d(&mut self) -> Float {
        match self {
            Sampler::Halton(sampler) => sampler.get_1d(),
            Sampler::MLT(sampler) => sampler.get_1d(),
            Sampler::Random(sampler) => sampler.get_1d(),
            Sampler::Sobol(sampler) => sampler.get_1d(),
            Sampler::ZeroTwoSequence(sampler) => sampler.get_1d(),
        }
    }
    pub fn get_2d(&mut self) -> Point2f {
        match self {
            Sampler::Halton(sampler) => sampler.get_2d(),
            Sampler::MLT(sampler) => sampler.get_2d(),
            Sampler::Random(sampler) => sampler.get_2d(),
            Sampler::Sobol(sampler) => sampler.get_2d(),
            Sampler::ZeroTwoSequence(sampler) => sampler.get_2d(),
        }
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
    pub fn request_2d_array(&mut self, n: i32) {
        match self {
            Sampler::Halton(sampler) => sampler.request_2d_array(n),
            Sampler::MLT(sampler) => sampler.request_2d_array(n),
            Sampler::Random(sampler) => sampler.request_2d_array(n),
            Sampler::Sobol(sampler) => sampler.request_2d_array(n),
            Sampler::ZeroTwoSequence(sampler) => sampler.request_2d_array(n),
        }
    }
    pub fn round_count(&self, count: i32) -> i32 {
        match self {
            Sampler::Halton(sampler) => sampler.round_count(count),
            Sampler::MLT(sampler) => sampler.round_count(count),
            Sampler::Random(sampler) => sampler.round_count(count),
            Sampler::Sobol(sampler) => sampler.round_count(count),
            Sampler::ZeroTwoSequence(sampler) => sampler.round_count(count),
        }
    }
    pub fn get_2d_array(&mut self, n: i32) -> Option<&[Point2f]> {
        match self {
            Sampler::Halton(sampler) => sampler.get_2d_array(n),
            Sampler::MLT(sampler) => sampler.get_2d_array(n),
            Sampler::Random(sampler) => sampler.get_2d_array(n),
            Sampler::Sobol(sampler) => sampler.get_2d_array(n),
            Sampler::ZeroTwoSequence(sampler) => sampler.get_2d_array(n),
        }
    }
    pub fn get_2d_arrays(&mut self, n: i32) -> (Option<&[Point2f]>, Option<&[Point2f]>) {
        match self {
            Sampler::Halton(sampler) => sampler.get_2d_arrays(n),
            Sampler::MLT(sampler) => sampler.get_2d_arrays(n),
            Sampler::Random(sampler) => sampler.get_2d_arrays(n),
            Sampler::Sobol(sampler) => sampler.get_2d_arrays(n),
            Sampler::ZeroTwoSequence(sampler) => sampler.get_2d_arrays(n),
        }
    }
    pub fn get_2d_array_vec(&mut self, n: i32) -> Vec<Point2f> {
        match self {
            Sampler::Halton(sampler) => sampler.get_2d_array_vec(n),
            Sampler::MLT(sampler) => sampler.get_2d_array_vec(n),
            Sampler::Random(sampler) => sampler.get_2d_array_vec(n),
            Sampler::Sobol(sampler) => sampler.get_2d_array_vec(n),
            Sampler::ZeroTwoSequence(sampler) => sampler.get_2d_array_vec(n),
        }
    }
    pub fn start_next_sample(&mut self) -> bool {
        match self {
            Sampler::Halton(sampler) => sampler.start_next_sample(),
            Sampler::MLT(sampler) => sampler.start_next_sample(),
            Sampler::Random(sampler) => sampler.start_next_sample(),
            Sampler::Sobol(sampler) => sampler.start_next_sample(),
            Sampler::ZeroTwoSequence(sampler) => sampler.start_next_sample(),
        }
    }
    pub fn reseed(&mut self, seed: u64) {
        match self {
            Sampler::Halton(sampler) => sampler.reseed(seed),
            Sampler::MLT(sampler) => sampler.reseed(seed),
            Sampler::Random(sampler) => sampler.reseed(seed),
            Sampler::Sobol(sampler) => sampler.reseed(seed),
            Sampler::ZeroTwoSequence(sampler) => sampler.reseed(seed),
        }
    }
    pub fn get_current_pixel(&self) -> Point2i {
        match self {
            Sampler::Halton(sampler) => sampler.get_current_pixel(),
            Sampler::MLT(sampler) => sampler.get_current_pixel(),
            Sampler::Random(sampler) => sampler.get_current_pixel(),
            Sampler::Sobol(sampler) => sampler.get_current_pixel(),
            Sampler::ZeroTwoSequence(sampler) => sampler.get_current_pixel(),
        }
    }
    pub fn get_current_sample_number(&self) -> i64 {
        match self {
            Sampler::Halton(sampler) => sampler.get_current_sample_number(),
            Sampler::MLT(sampler) => sampler.get_current_sample_number(),
            Sampler::Random(sampler) => sampler.get_current_sample_number(),
            Sampler::Sobol(sampler) => sampler.get_current_sample_number(),
            Sampler::ZeroTwoSequence(sampler) => sampler.get_current_sample_number(),
        }
    }
    pub fn get_samples_per_pixel(&self) -> i64 {
        match self {
            Sampler::Halton(sampler) => sampler.get_samples_per_pixel(),
            Sampler::MLT(sampler) => sampler.get_samples_per_pixel(),
            Sampler::Random(sampler) => sampler.get_samples_per_pixel(),
            Sampler::Sobol(sampler) => sampler.get_samples_per_pixel(),
            Sampler::ZeroTwoSequence(sampler) => sampler.get_samples_per_pixel(),
        }
    }
}
