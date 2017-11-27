// pbrt
use core::lowdiscrepancy::{sobol_interval_to_index, sobol_sample};
use core::pbrt::Float;
use core::pbrt::{clamp_t, is_power_of_2, log_2_int_u32, round_up_pow2_32, round_up_pow2_64};
use core::rng::FLOAT_ONE_MINUS_EPSILON;
use core::sampler::{GlobalSampler, Sampler};
use core::sobolmatrices::NUM_SOBOL_DIMENSIONS;
use geometry::{Bounds2i, Point2f, Point2i, Vector2i};

// see sobol.h

pub struct SobolSampler {
    pub samples_per_pixel: i64,
    pub sample_bounds: Bounds2i,
    pub resolution: i32,
    pub log_2_resolution: i32,
    // inherited from class GlobalSampler (see sampler.h)
    pub dimension: i64,
    pub interval_sample_index: u64,
    pub array_start_dim: i64,
    pub array_end_dim: i64,
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

impl SobolSampler {
    pub fn new(samples_per_pixel: i64, sample_bounds: Bounds2i) -> Self {
        let mut samples_per_pixel: i64 = samples_per_pixel;
        if !is_power_of_2(samples_per_pixel) {
            samples_per_pixel = round_up_pow2_64(samples_per_pixel);
            println!("WARNING: Non power-of-two sample count rounded up to {:?} for SobolSampler.",
                     samples_per_pixel);
        }
        let resolution: i32 =
            round_up_pow2_32(sample_bounds.diagonal().x.max(sample_bounds.diagonal().y));
        let log_2_resolution: i32 = log_2_int_u32(resolution as u32);
        if resolution > 0_i32 {
            assert!(1_i32 << log_2_resolution == resolution);
        }
        SobolSampler {
            samples_per_pixel: samples_per_pixel,
            sample_bounds: sample_bounds,
            resolution: resolution,
            log_2_resolution: log_2_resolution,
            dimension: 0_i64,
            interval_sample_index: 0_u64,
            array_start_dim: 5_i64, // static const int arrayStartDim = 5;
            array_end_dim: 0_i64,
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
    pub fn get_index_for_sample(&self, sample_num: u64) -> u64 {
        let v: Vector2i = self.current_pixel - self.sample_bounds.p_min;
        sobol_interval_to_index(self.log_2_resolution as u32,
                                sample_num,
                                &Point2i { x: v.x, y: v.y })
    }
    pub fn sample_dimension(&self, index: u64, dim: i64) -> Float {
        if dim >= NUM_SOBOL_DIMENSIONS as i64 {
            panic!("SobolSampler can only sample up to {} dimensions! Exiting.",
                   NUM_SOBOL_DIMENSIONS);
        }
        let mut s: Float = sobol_sample(index as i64, dim as i32, 0_u64);
        // remap Sobol$'$ dimensions used for pixel samples
        if dim == 0 || dim == 1 {
            s = s * self.resolution as Float + self.sample_bounds.p_min[dim as u8] as Float;
            s = clamp_t(s - self.current_pixel[dim as u8] as Float,
                        0.0 as Float,
                        FLOAT_ONE_MINUS_EPSILON);
        }
        s
    }
}

impl Sampler for SobolSampler {
    fn start_pixel(&mut self, p: Point2i) {
        // TODO: ProfilePhase _(Prof::StartPixel);
        // Sampler::StartPixel(p);
        self.current_pixel = p;
        self.current_pixel_sample_index = 0_i64;
        self.array_1d_offset = 0_usize;
        self.array_2d_offset = 0_usize;
        // GlobalSampler::StartPixel(p);
        self.dimension = 0_i64;
        self.interval_sample_index = self.get_index_for_sample(0_u64);
        // compute _self.array_end_dim_ for dimensions used for array samples
        self.array_end_dim = self.array_start_dim + self.sample_array_1d.len() as i64 +
                             2_i64 * self.sample_array_2d.len() as i64;
        // compute 1D array samples for _GlobalSampler_
        for i in 0..self.samples_1d_array_sizes.len() {
            let n_samples = self.samples_1d_array_sizes[i] * self.samples_per_pixel as i32;
            for j in 0..n_samples {
                let index: u64 = self.get_index_for_sample(j as u64);
                self.sample_array_1d[i as usize][j as usize] =
                    self.sample_dimension(index, self.array_start_dim + i as i64);
            }
        }
        // compute 2D array samples for _GlobalSampler_
        let mut dim: i64 = self.array_start_dim + self.samples_1d_array_sizes.len() as i64;
        for i in 0..self.samples_2d_array_sizes.len() {
            let n_samples: usize = self.samples_2d_array_sizes[i] as usize *
                                   self.samples_per_pixel as usize;
            for j in 0..n_samples {
                let idx: u64 = self.get_index_for_sample(j as u64);
                self.sample_array_2d[i][j].x = self.sample_dimension(idx, dim);
                self.sample_array_2d[i][j].y = self.sample_dimension(idx, dim + 1_i64);
            }
            dim += 2_i64;
        }
        assert!(self.array_end_dim == dim);
    }
    fn get_1d(&mut self) -> Float {
        // TODO: ProfilePhase _(Prof::GetSample);
        if self.dimension >= self.array_start_dim && self.dimension < self.array_end_dim {
            self.dimension = self.array_end_dim;
        }
        // call first (in C++: return SampleDimension(intervalSampleIndex, dimension++));
        let ret: Float = self.sample_dimension(self.interval_sample_index, self.dimension);
        self.dimension += 1;
        // then return
        ret
    }
    fn get_2d(&mut self) -> Point2f {
        // TODO: ProfilePhase _(Prof::GetSample);
        if self.dimension + 1 >= self.array_start_dim && self.dimension < self.array_end_dim {
            self.dimension = self.array_end_dim;
        }
        // C++: call y first
        let y = self.sample_dimension(self.interval_sample_index, self.dimension + 1);
        let x = self.sample_dimension(self.interval_sample_index, self.dimension);
        let p: Point2f = Point2f { x: x, y: y };
        self.dimension += 2;
        return p;
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
    fn get_2d_array(&mut self, n: i32) -> Vec<Point2f> {
        let mut samples: Vec<Point2f> = Vec::new();
        if self.array_2d_offset == self.sample_array_2d.len() {
            return samples;
        }
        assert_eq!(self.samples_2d_array_sizes[self.array_2d_offset], n);
        assert!(self.current_pixel_sample_index < self.samples_per_pixel,
                "self.current_pixel_sample_index ({}) < self.samples_per_pixel ({})",
                self.current_pixel_sample_index,
                self.samples_per_pixel);
        let start: usize = (self.current_pixel_sample_index * n as i64) as usize;
        let end: usize = start + n as usize;
        samples = self.sample_array_2d[self.array_2d_offset][start..end].to_vec();
        self.array_2d_offset += 1;
        samples
    }
    fn start_next_sample(&mut self) -> bool {
        self.dimension = 0_i64;
        self.interval_sample_index =
            self.get_index_for_sample(self.current_pixel_sample_index as u64 + 1_u64);
        // Sampler::StartNextSample();
        // reset array offsets for next pixel sample
        self.array_1d_offset = 0_usize;
        self.array_2d_offset = 0_usize;
        self.current_pixel_sample_index += 1_i64;
        self.current_pixel_sample_index < self.samples_per_pixel
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

impl GlobalSampler for SobolSampler {
}

impl Clone for SobolSampler {
    fn clone(&self) -> SobolSampler {
        SobolSampler {
            samples_per_pixel: self.samples_per_pixel,
            sample_bounds: self.sample_bounds,
            resolution: self.resolution,
            log_2_resolution: self.log_2_resolution,
            dimension: self.dimension,
            interval_sample_index: self.interval_sample_index,
            array_start_dim: self.array_start_dim,
            array_end_dim: self.array_end_dim,
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
