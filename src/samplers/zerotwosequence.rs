// pbrt
use core::camera::CameraSample;
use core::lowdiscrepancy::{sobol_2d, van_der_corput};
use core::pbrt::Float;
use core::pbrt::round_up_pow2_32;
use core::rng::Rng;
use core::sampler::Sampler;
use geometry::{Point2i, Point2f};

// see zerotwosequence.h

#[derive(Debug,Clone)]
pub struct ZeroTwoSequenceSampler {
    pub samples_per_pixel: i64,
    pub n_sampled_dimensions: i64,
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
    pub samples_1d_array: Vec<Vec<Float>>,
    pub samples_2d_array: Vec<Vec<Point2f>>,
    pub array_1d_offset: usize,
    pub array_2d_offset: usize,
}

impl Default for ZeroTwoSequenceSampler {
    fn default() -> Self {
        let mut lds: ZeroTwoSequenceSampler = ZeroTwoSequenceSampler {
            samples_per_pixel: 1_i64,
            n_sampled_dimensions: 4_i64,
            samples_1d: Vec::new(),
            samples_2d: Vec::new(),
            current_1d_dimension: 0_i32,
            current_2d_dimension: 0_i32,
            rng: Rng::default(),
            current_pixel: Point2i::default(),
            current_pixel_sample_index: 0_i64,
            samples_1d_array_sizes: Vec::new(),
            samples_2d_array_sizes: Vec::new(),
            samples_1d_array: Vec::new(),
            samples_2d_array: Vec::new(),
            array_1d_offset: 0_usize,
            array_2d_offset: 0_usize,
        };
        for _i in 0..lds.n_sampled_dimensions {
            let additional_1d: Vec<Float> = vec![0.0; lds.samples_per_pixel as usize];
            let additional_2d: Vec<Point2f> =
                vec![Point2f::default(); lds.samples_per_pixel as usize];
            lds.samples_1d.push(additional_1d);
            lds.samples_2d.push(additional_2d);
        }
        lds
    }
}

impl ZeroTwoSequenceSampler {
    pub fn new(samples_per_pixel: i64, n_sampled_dimensions: i64) -> Self {
        let mut lds: ZeroTwoSequenceSampler = ZeroTwoSequenceSampler {
            samples_per_pixel: samples_per_pixel,
            n_sampled_dimensions: n_sampled_dimensions,
            samples_1d: Vec::new(),
            samples_2d: Vec::new(),
            current_1d_dimension: 0_i32,
            current_2d_dimension: 0_i32,
            rng: Rng::default(),
            current_pixel: Point2i::default(),
            current_pixel_sample_index: 0_i64,
            samples_1d_array_sizes: Vec::new(),
            samples_2d_array_sizes: Vec::new(),
            samples_1d_array: Vec::new(),
            samples_2d_array: Vec::new(),
            array_1d_offset: 0_usize,
            array_2d_offset: 0_usize,
        };
        for _i in 0..lds.n_sampled_dimensions {
            let additional_1d: Vec<Float> = vec![0.0; lds.samples_per_pixel as usize];
            let additional_2d: Vec<Point2f> =
                vec![Point2f::default(); lds.samples_per_pixel as usize];
            lds.samples_1d.push(additional_1d);
            lds.samples_2d.push(additional_2d);
        }
        lds
    }
    pub fn get_camera_sample(&mut self, p_raster: Point2i) -> CameraSample {
        let mut cs: CameraSample = CameraSample::default();
        cs.p_film = Point2f {
            x: p_raster.x as Float,
            y: p_raster.y as Float,
        } + self.get_2d();
        cs.time = self.get_1d();
        cs.p_lens = self.get_2d();
        cs
    }
    pub fn start_pixel(&mut self, p: Point2i) {
        // TODO: ProfilePhase _(Prof::StartPixel);
        // generate 1D and 2D pixel sample components using $(0,2)$-sequence
        for samples in &mut self.samples_1d {
            van_der_corput(1, self.samples_per_pixel as i32, samples, &mut self.rng);
        }
        for samples in &mut self.samples_2d {
            sobol_2d(1, self.samples_per_pixel as i32, samples, &mut self.rng);
        }
        // generate 1D and 2D array samples using $(0,2)$-sequence
        for i in 0..self.samples_1d_array_sizes.len() {
            let samples: &mut [Float] = self.samples_1d_array[i].as_mut_slice();
            van_der_corput(self.samples_1d_array_sizes[i],
                           self.samples_per_pixel as i32,
                           samples,
                           &mut self.rng);
        }
        for i in 0..self.samples_2d_array_sizes.len() {
            let samples: &mut [Point2f] = self.samples_2d_array[i].as_mut_slice();
            sobol_2d(self.samples_2d_array_sizes[i],
                     self.samples_per_pixel as i32,
                     samples,
                     &mut self.rng);
        }
        // PixelSampler::StartPixel(p);
        self.current_pixel = p;
        self.current_pixel_sample_index = 0_i64;
        // reset array offsets for next pixel sample
        self.array_1d_offset = 0_usize;
        self.array_2d_offset = 0_usize;
    }
    pub fn round_count(&self, count: i32) -> i32 {
        round_up_pow2_32(count)
    }
    // PixelSampler public methods
    pub fn start_next_sample(&mut self) -> bool {
        self.current_1d_dimension = 0_i32;
        self.current_2d_dimension = 0_i32;
        // Sampler::StartNextSample()
        // reset array offsets for next pixel sample
        self.array_1d_offset = 0_usize;
        self.array_2d_offset = 0_usize;
        self.current_pixel_sample_index += 1_i64;
        self.current_pixel_sample_index < self.samples_per_pixel
    }
    // inherited from class Sampler (see sampler.h)
    pub fn request_2d_array(&mut self, n: i32) {
        assert_eq!(self.round_count(n), n);
        self.samples_2d_array_sizes.push(n);
        let size: usize = (n * self.samples_per_pixel as i32) as usize;
        let additional_points: Vec<Point2f> = vec![Point2f::default(); size];
        self.samples_2d_array.push(additional_points);
    }
    pub fn get_2d_array(&mut self, n: i32) -> Vec<Point2f> {
        let mut samples: Vec<Point2f> = Vec::new();
        if self.array_2d_offset == self.samples_2d_array.len() {
            return samples;
        }
        assert_eq!(self.samples_2d_array_sizes[self.array_2d_offset], n);
        assert!(self.current_pixel_sample_index < self.samples_per_pixel,
                "self.current_pixel_sample_index ({}) < self.samples_per_pixel ({})",
                self.current_pixel_sample_index,
                self.samples_per_pixel);
        let start: usize = (self.current_pixel_sample_index * n as i64) as usize;
        let end: usize = start + n as usize;
        samples = self.samples_2d_array[self.array_2d_offset][start..end].to_vec();
        self.array_2d_offset += 1;
        samples
    }
    pub fn get_current_sample_number(&self) -> i64 {
        self.current_pixel_sample_index
    }
}

impl Sampler for ZeroTwoSequenceSampler {
    fn get_1d(&mut self) -> Float {
        // TODO: ProfilePhase _(Prof::GetSample);
        assert!(self.current_pixel_sample_index < self.samples_per_pixel,
                "current_pixel_sample_index = {}, samples_per_pixel = {}",
                self.current_pixel_sample_index,
                self.samples_per_pixel);
        if self.current_1d_dimension < self.samples_1d.len() as i32 {
            let sample: Float = self.samples_1d[self.current_1d_dimension as usize][self.current_pixel_sample_index as
            usize];
            self.current_1d_dimension += 1;
            sample
        } else {
            self.rng.uniform_float()
        }
    }
    fn get_2d(&mut self) -> Point2f {
        // TODO: ProfilePhase _(Prof::GetSample);
        assert!(self.current_pixel_sample_index < self.samples_per_pixel,
                "current_pixel_sample_index = {}, samples_per_pixel = {}",
                self.current_pixel_sample_index,
                self.samples_per_pixel);
        if self.current_2d_dimension < self.samples_2d.len() as i32 {
            let sample: Point2f = self.samples_2d[self.current_2d_dimension as usize][self.current_pixel_sample_index as
            usize];
            self.current_2d_dimension += 1;
            sample
        } else {
            // C++ call order for Point2f(rng.UniformFloat(), rng.UniformFloat());
            let y = self.rng.uniform_float();
            let x = self.rng.uniform_float();
            Point2f { x: x, y: y }
        }
    }
}
