// std
use std::sync::RwLock;
// pbrt
use core::camera::CameraSample;
use core::lowdiscrepancy::{PRIME_SUMS, PRIME_TABLE_SIZE};
use core::lowdiscrepancy::{inverse_radical_inverse, scrambled_radical_inverse, radical_inverse};
use core::pbrt::Float;
use core::pbrt::mod_t;
use core::sampler::Sampler;
use geometry::{Point2f, Point2i};

// see halton.h

pub const K_MAX_RESOLUTION: i32 = 128_i32;

pub struct HaltonSampler {
    pub samples_per_pixel: i64,
    pub radical_inverse_permutations: Vec<u16>,
    pub base_scales: Point2i,
    pub base_exponents: Point2i,
    pub sample_stride: u64,
    pub mult_inverse: [i64; 2],
    pub pixel_for_offset: RwLock<Point2i>,
    pub offset_for_current_pixel: RwLock<u64>,
    pub sample_at_pixel_center: bool, // default: false
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

impl HaltonSampler {
    pub fn get_index_for_sample(&self, sample_num: u64) -> u64 {
        let pixel_for_offset: Point2i = *self.pixel_for_offset.read().unwrap();
        if self.current_pixel != pixel_for_offset {
            // compute Halton sample offset for _self.current_pixel_
            let mut offset_for_current_pixel = self.offset_for_current_pixel.write().unwrap();
            *offset_for_current_pixel = 0_u64;
            if self.sample_stride > 1_u64 {
                let pm: Point2i = Point2i {
                    x: mod_t(self.current_pixel[0], K_MAX_RESOLUTION),
                    y: mod_t(self.current_pixel[1], K_MAX_RESOLUTION),
                };
                for i in 0..2 {
                    let dim_offset: u64;
                    if i == 0 {
                        dim_offset =
                            inverse_radical_inverse(2, pm[i] as u64, self.base_exponents[i] as u64);
                    } else {
                        dim_offset =
                            inverse_radical_inverse(3, pm[i] as u64, self.base_exponents[i] as u64);
                    }
                    *offset_for_current_pixel +=
                        dim_offset * (self.sample_stride / self.base_scales[i] as u64) as u64 *
                        self.mult_inverse[i as usize] as u64;
                }
                *offset_for_current_pixel %= self.sample_stride as u64;
            }
            let mut pixel_for_offset = self.pixel_for_offset.write().unwrap();
            *pixel_for_offset = self.current_pixel;
        }
        let offset_for_current_pixel: u64 = *self.offset_for_current_pixel.read().unwrap();
        offset_for_current_pixel + sample_num * self.sample_stride
    }
    pub fn sample_dimension(&self, index: u64, dim: i64) -> Float {
        if self.sample_at_pixel_center && (dim == 0 || dim == 1) {
            return 0.5 as Float;
        }
        if dim == 0 {
            radical_inverse(dim as u16, index >> self.base_exponents[0] as u64)
        } else if dim == 1 {
            radical_inverse(dim as u16, index / self.base_scales[1] as u64)
        } else {
            scrambled_radical_inverse(dim as u16, index, self.permutation_for_dimension(dim))
        }
    }
    fn permutation_for_dimension(&self, dim: i64) -> Vec<u16> {
        if dim >= PRIME_TABLE_SIZE as i64 {
            panic!("FATAL: HaltonSampler can only sample {:?} dimensions.",
                   PRIME_TABLE_SIZE);
        }
        self.radical_inverse_permutations[PRIME_SUMS[dim as usize] as usize..].to_vec()
    }
}

impl Sampler for HaltonSampler {
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
        // int dim = self.array_start_dim + self.samples_1d_array_sizes.len();
        // for (size_t i = 0; i < samples2DArraySizes.len(); ++i) {
        //     int nSamples = samples2DArraySizes[i] * self.samples_per_pixel;
        //     for (int j = 0; j < nSamples; ++j) {
        //         int64_t idx = GetIndexForSample(j);
        //         self.sample_array_2d[i][j].x = SampleDimension(idx, dim);
        //         self.sample_array_2d[i][j].y = SampleDimension(idx, dim + 1);
        //     }
        //     dim += 2;
        // }
        // CHECK_EQ(self.array_end_dim, dim);
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
    fn get_current_sample_number(&self) -> i64 {
        self.current_pixel_sample_index
    }
    fn get_samples_per_pixel(&self) -> i64 {
        self.samples_per_pixel
    }
}

impl Clone for HaltonSampler {
    fn clone(&self) -> HaltonSampler {
        let pixel_for_offset: Point2i = *self.pixel_for_offset.read().unwrap();
        let offset_for_current_pixel: u64 = *self.offset_for_current_pixel.read().unwrap();
        HaltonSampler {
            samples_per_pixel: self.samples_per_pixel,
            radical_inverse_permutations: self.radical_inverse_permutations
                .iter()
                .cloned()
                .collect(),
            base_scales: self.base_scales,
            base_exponents: self.base_exponents,
            sample_stride: self.sample_stride,
            mult_inverse: self.mult_inverse,
            pixel_for_offset: RwLock::new(pixel_for_offset),
            offset_for_current_pixel: RwLock::new(offset_for_current_pixel),
            sample_at_pixel_center: self.sample_at_pixel_center,
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
