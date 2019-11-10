// pbrt
use crate::core::geometry::{Point2f, Point2i};
use crate::core::paramset::ParamSet;
use crate::core::pbrt::Float;
use crate::core::rng::Rng;
use crate::core::sampler::Sampler;

// see random.h

#[derive(Clone)]
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
    pub fn new(samples_per_pixel: i64) -> Self {
        RandomSampler {
            samples_per_pixel,
            rng: Rng::default(),
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
    pub fn create(params: &ParamSet) -> Box<dyn Sampler + Sync + Send> {
        let nsamp: i32 = params.find_one_int("pixelsamples", 4);
        // TODO: if (PbrtOptions.quickRender) nsamp = 1;
        Box::new(RandomSampler::new(nsamp as i64))
    }
}

impl Sampler for RandomSampler {
    fn start_pixel(&mut self, p: &Point2i) {
        // TODO: ProfilePhase _(Prof::StartPixel);
        for i in 0..self.sample_array_1d.len() {
            for j in 0..self.sample_array_1d[i].len() {
                self.sample_array_1d[i][j] = self.rng.uniform_float();
            }
        }
        for i in 0..self.sample_array_2d.len() {
            for j in 0..self.sample_array_2d[i].len() {
                // C++: call y first
                let y = self.rng.uniform_float();
                let x = self.rng.uniform_float();
                self.sample_array_2d[i][j].x = x;
                self.sample_array_2d[i][j].y = y;
            }
        }
        // Sampler::StartPixel(p);
        self.current_pixel = *p;
        self.current_pixel_sample_index = 0_i64;
        self.array_1d_offset = 0_usize;
        self.array_2d_offset = 0_usize;
    }
    fn get_1d(&mut self) -> Float {
        // TODO: ProfilePhase _(Prof::GetSample);
        assert!(self.current_pixel_sample_index < self.samples_per_pixel);
        self.rng.uniform_float()
    }
    fn get_2d(&mut self) -> Point2f {
        // TODO: ProfilePhase _(Prof::GetSample);
        assert!(self.current_pixel_sample_index < self.samples_per_pixel);
        // C++: call x first
        let x = self.rng.uniform_float();
        let y = self.rng.uniform_float();
        Point2f { x, y }
    }
    fn reseed(&mut self, seed: u64) {
        self.rng.set_sequence(seed);
    }
    fn request_2d_array(&mut self, n: i32) {
        assert_eq!(self.round_count(n), n);
        self.samples_2d_array_sizes.push(n);
        let size: usize = (n * self.samples_per_pixel as i32) as usize;
        let additional_points: Vec<Point2f> = vec![Point2f::default(); size];
        self.sample_array_2d.push(additional_points);
    }
    fn round_count(&self, count: i32) -> i32 {
        count
    }
    fn get_2d_array(&mut self, n: i32) -> Option<&[Point2f]> {
        let mut samples: Vec<Point2f> = Vec::new();
        if self.array_2d_offset == self.sample_array_2d.len() {
            return None;
        }
        assert_eq!(self.samples_2d_array_sizes[self.array_2d_offset], n);
        assert!(
            self.current_pixel_sample_index < self.samples_per_pixel,
            "self.current_pixel_sample_index ({}) < self.samples_per_pixel ({})",
            self.current_pixel_sample_index,
            self.samples_per_pixel
        );
        let start: usize = (self.current_pixel_sample_index * n as i64) as usize;
        let end: usize = start + n as usize;
        self.array_2d_offset += 1;
        Some(&self.sample_array_2d[self.array_2d_offset - 1][start..end])
    }
    fn start_next_sample(&mut self) -> bool {
        // reset array offsets for next pixel sample
        self.array_1d_offset = 0_usize;
        self.array_2d_offset = 0_usize;
        self.current_pixel_sample_index += 1_i64;
        self.current_pixel_sample_index < self.samples_per_pixel
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
