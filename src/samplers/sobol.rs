// pbrt
use crate::core::geometry::{Bounds2i, Point2f, Point2i, Vector2i, XYEnum};
use crate::core::lowdiscrepancy::{sobol_interval_to_index, sobol_sample};
use crate::core::paramset::ParamSet;
use crate::core::pbrt::Float;
use crate::core::pbrt::{
    clamp_t, is_power_of_2, log_2_int_u32, round_up_pow2_32, round_up_pow2_64,
};
use crate::core::rng::FLOAT_ONE_MINUS_EPSILON;
use crate::core::sampler::Sampler;
use crate::core::sobolmatrices::NUM_SOBOL_DIMENSIONS;

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
    pub fn new(samples_per_pixel: i64, sample_bounds: &Bounds2i) -> Self {
        let mut samples_per_pixel: i64 = samples_per_pixel;
        if !is_power_of_2(samples_per_pixel) {
            samples_per_pixel = round_up_pow2_64(samples_per_pixel);
            println!(
                "WARNING: Non power-of-two sample count rounded up to {:?} for SobolSampler.",
                samples_per_pixel
            );
        }
        let resolution: i32 =
            round_up_pow2_32(sample_bounds.diagonal().x.max(sample_bounds.diagonal().y));
        let log_2_resolution: i32 = log_2_int_u32(resolution as u32);
        if resolution > 0_i32 {
            assert!(1_i32 << log_2_resolution == resolution);
        }
        SobolSampler {
            samples_per_pixel,
            sample_bounds: Bounds2i {
                p_min: Point2i {
                    x: sample_bounds.p_min.x,
                    y: sample_bounds.p_min.y,
                },
                p_max: Point2i {
                    x: sample_bounds.p_max.x,
                    y: sample_bounds.p_max.y,
                },
            },
            resolution,
            log_2_resolution,
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
    pub fn clone_with_seed(&self, _seed: u64) -> Box<Sampler> {
        let sobol_sampler = SobolSampler {
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
            samples_1d_array_sizes: self.samples_1d_array_sizes.to_vec(),
            samples_2d_array_sizes: self.samples_2d_array_sizes.to_vec(),
            sample_array_1d: self.sample_array_1d.to_vec(),
            sample_array_2d: self.sample_array_2d.to_vec(),
            array_1d_offset: self.array_1d_offset,
            array_2d_offset: self.array_2d_offset,
        };
        let sampler = Sampler::Sobol(sobol_sampler);
        Box::new(sampler)
    }
    pub fn create(params: &ParamSet, sample_bounds: &Bounds2i) -> Box<Sampler> {
        let nsamp: i32 = params.find_one_int("pixelsamples", 16);
        // TODO: if (PbrtOptions.quickRender) nsamp = 1;
        Box::new(Sampler::Sobol(SobolSampler::new(
            nsamp as i64,
            sample_bounds,
        )))
    }
    pub fn get_index_for_sample(&self, sample_num: u64) -> u64 {
        let v: Vector2i = self.current_pixel - self.sample_bounds.p_min;
        sobol_interval_to_index(
            self.log_2_resolution as u32,
            sample_num,
            Point2i { x: v.x, y: v.y },
        )
    }
    pub fn sample_dimension(&self, index: u64, dim: i64) -> Float {
        if dim >= NUM_SOBOL_DIMENSIONS as i64 {
            panic!(
                "SobolSampler can only sample up to {} dimensions! Exiting.",
                NUM_SOBOL_DIMENSIONS
            );
        }
        let mut s: Float = sobol_sample(index as i64, dim as i32, 0_u64);
        // remap Sobol$'$ dimensions used for pixel samples
        if dim == 0 || dim == 1 {
            let dim_i: XYEnum = match dim {
                0 => XYEnum::X,
                _ => XYEnum::Y,
            };
            s = s * self.resolution as Float + self.sample_bounds.p_min[dim_i] as Float;
            s = clamp_t(
                s - self.current_pixel[dim_i] as Float,
                0.0 as Float,
                FLOAT_ONE_MINUS_EPSILON,
            );
        }
        s
    }
    // Sampler
    pub fn start_pixel(&mut self, p: Point2i) {
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
        self.array_end_dim = self.array_start_dim
            + self.sample_array_1d.len() as i64
            + 2_i64 * self.sample_array_2d.len() as i64;
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
            let n_samples: usize =
                self.samples_2d_array_sizes[i] as usize * self.samples_per_pixel as usize;
            for j in 0..n_samples {
                let idx: u64 = self.get_index_for_sample(j as u64);
                self.sample_array_2d[i][j].x = self.sample_dimension(idx, dim);
                self.sample_array_2d[i][j].y = self.sample_dimension(idx, dim + 1_i64);
            }
            dim += 2_i64;
        }
        assert!(self.array_end_dim == dim);
    }
    pub fn get_1d(&mut self) -> Float {
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
    pub fn get_2d(&mut self) -> Point2f {
        // TODO: ProfilePhase _(Prof::GetSample);
        if self.dimension + 1 >= self.array_start_dim && self.dimension < self.array_end_dim {
            self.dimension = self.array_end_dim;
        }
        // C++: call y first
        let y = self.sample_dimension(self.interval_sample_index, self.dimension + 1);
        let x = self.sample_dimension(self.interval_sample_index, self.dimension);
        let p: Point2f = Point2f { x, y };
        self.dimension += 2;
        p
    }
    pub fn get_2d_sample(&self, array_idx: usize, idx: usize) -> Point2f {
        self.sample_array_2d[array_idx][idx]
    }
    pub fn request_2d_array(&mut self, n: i32) {
        assert_eq!(self.round_count(n), n);
        self.samples_2d_array_sizes.push(n);
        let size: usize = (n * self.samples_per_pixel as i32) as usize;
        let additional_points: Vec<Point2f> = vec![Point2f::default(); size];
        self.sample_array_2d.push(additional_points);
    }
    pub fn round_count(&self, count: i32) -> i32 {
        count
    }
    pub fn get_2d_array(&mut self, n: i32) -> Option<&[Point2f]> {
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
    pub fn get_2d_arrays(&mut self, n: i32) -> (Option<&[Point2f]>, Option<&[Point2f]>) {
        if self.array_2d_offset == self.sample_array_2d.len() {
            return (None, None);
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
        let ret1 = &self.sample_array_2d[self.array_2d_offset - 1][start..end];
        // repeat code from above
        if self.array_2d_offset == self.sample_array_2d.len() {
            return (None, None);
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
        let ret2 = &self.sample_array_2d[self.array_2d_offset - 1][start..end];
        // return tuple
        (Some(ret1), Some(ret2))
    }
    pub fn get_2d_array_idxs(&mut self, n: i32) -> (bool, usize, usize) {
        if self.array_2d_offset == self.sample_array_2d.len() {
            return (true, 0_usize, 0_usize);
        }
        assert_eq!(self.samples_2d_array_sizes[self.array_2d_offset], n);
        assert!(
            self.current_pixel_sample_index < self.samples_per_pixel,
            "self.current_pixel_sample_index ({}) < self.samples_per_pixel ({})",
            self.current_pixel_sample_index,
            self.samples_per_pixel
        );
        let start: usize = (self.current_pixel_sample_index * n as i64) as usize;
        let idx: usize = self.array_2d_offset;
        self.array_2d_offset += 1;
        (false, idx, start)
    }
    pub fn start_next_sample(&mut self) -> bool {
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
    pub fn reseed(&mut self, _seed: u64) {
        // do nothing
    }
    pub fn get_current_pixel(&self) -> Point2i {
        self.current_pixel
    }
    pub fn get_current_sample_number(&self) -> i64 {
        self.current_pixel_sample_index
    }
    pub fn get_samples_per_pixel(&self) -> i64 {
        self.samples_per_pixel
    }
    // GlobalSampler
    pub fn set_sample_number(&mut self, sample_num: i64) -> bool {
        // GlobalSampler::SetSampleNumber(...)
        self.dimension = 0_i64;
        self.interval_sample_index = self.get_index_for_sample(sample_num as u64);
        // reset array offsets for next pixel sample
        self.array_1d_offset = 0_usize;
        self.array_2d_offset = 0_usize;
        self.current_pixel_sample_index = sample_num;
        self.current_pixel_sample_index < self.samples_per_pixel
    }
}
