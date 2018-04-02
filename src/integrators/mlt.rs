// pbrt
use core::geometry::{Point2f, Point2i};
use core::pbrt::Float;
use core::rng::Rng;
use core::sampler::Sampler;

pub const N_SAMPLE_STREAMS: u8 = 3;

pub struct MLTSampler {
    pub samples_per_pixel: i64,
    pub rng: Rng,
    pub sigma: Float,
    pub large_step_probability: Float,
    pub stream_count: i32,
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

impl MLTSampler {
    pub fn new(
        mutations_per_pixel: i64,
        rng_sequence_index: u64,
        sigma: Float,
        large_step_probability: Float,
        stream_count: i32,
    ) -> Self {
        let mut rng: Rng = Rng::default();
        rng.set_sequence(rng_sequence_index);
        MLTSampler {
            samples_per_pixel: mutations_per_pixel,
            rng: rng,
            sigma: sigma,
            large_step_probability: large_step_probability,
            stream_count: stream_count,
            current_pixel: Point2i::default(),
            current_pixel_sample_index: 0_i64,
            samples_1d_array_sizes: Vec::new(),
            samples_2d_array_sizes: Vec::new(),
            sample_array_1d: Vec::new(),
            sample_array_2d: Vec::new(),
            array_1d_offset: 0_usize,
            array_2d_offset: 0_usize,
        }
    }
}

impl Sampler for MLTSampler {
    fn start_pixel(&mut self, p: &Point2i) {}
    fn get_1d(&mut self) -> Float {
        // WORK
        0.0 as Float
    }
    fn get_2d(&mut self) -> Point2f {
        // WORK
        Point2f::default()
    }
    fn reseed(&mut self, seed: u64) {
        // WORK
    }
    fn request_2d_array(&mut self, n: i32) {
        // WORK
    }
    fn round_count(&self, count: i32) -> i32 {
        // WORK
        0_i32
    }
    fn get_2d_array(&mut self, n: i32) -> Vec<Point2f> {
        let mut samples: Vec<Point2f> = Vec::new();
        // WORK
        samples
    }
    fn start_next_sample(&mut self) -> bool {
        // WORK
        false
    }
    fn get_current_pixel(&self) -> Point2i {
        self.current_pixel
    }
    fn get_current_sample_number(&self) -> i64 {
        self.current_pixel_sample_index
    }
    fn get_samples_per_pixel(&self) -> i64 {
        self.samples_per_pixel
    }
}

impl Clone for MLTSampler {
    fn clone(&self) -> MLTSampler {
        MLTSampler {
            samples_per_pixel: self.samples_per_pixel,
            rng: Rng::default(),
            sigma: self.sigma,
            large_step_probability: self.large_step_probability,
            stream_count: self.stream_count,
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

pub struct MLTIntegrator {
    pub max_depth: u32,
    pub n_bootstrap: u32,
    pub n_chains: u32,
    pub mutations_per_pixel: u32,
    pub sigma: Float,
    pub large_step_probability: Float,
}

impl MLTIntegrator {
    pub fn new(
        max_depth: u32,
        n_bootstrap: u32,
        n_chains: u32,
        mutations_per_pixel: u32,
        sigma: Float,
        large_step_probability: Float,
    ) -> Self {
        MLTIntegrator {
            max_depth: max_depth,
            n_bootstrap: n_bootstrap,
            n_chains: n_chains,
            mutations_per_pixel: mutations_per_pixel,
            sigma: sigma,
            large_step_probability: large_step_probability,
        }
    }
}
