// std
use std::sync::RwLock;
// pbrt
use core::camera::CameraSample;
use core::lowdiscrepancy::inverse_radical_inverse;
use core::pbrt::Float;
use core::pbrt::mod_t;
use core::sampler::Sampler;
use geometry::{Point2f, Point2i};

// see halton.h

pub const K_MAX_RESOLUTION: i32 = 128_i32;

pub struct HaltonSampler {
    pub samples_per_pixel: i64,
    pub base_scales: Point2i,
    pub base_exponents: Point2i,
    pub sample_stride: u64,
    pub mult_inverse: [i64; 2],
    pub pixel_for_offset: RwLock<Point2i>,
    pub offset_for_current_pixel: RwLock<u64>,
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
    pub samples_1d_array: Vec<Vec<Float>>,
    pub samples_2d_array: Vec<Vec<Point2f>>,
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
        // // Compute _arrayEndDim_ for dimensions used for array samples
        // arrayEndDim =
        //     arrayStartDim + sampleArray1D.size() + 2 * sampleArray2D.size();

        // // Compute 1D array samples for _GlobalSampler_
        // for (size_t i = 0; i < samples1DArraySizes.size(); ++i) {
        //     int nSamples = samples1DArraySizes[i] * samplesPerPixel;
        //     for (int j = 0; j < nSamples; ++j) {
        //         int64_t index = GetIndexForSample(j);
        //         sampleArray1D[i][j] = SampleDimension(index, arrayStartDim + i);
        //     }
        // }

        // // Compute 2D array samples for _GlobalSampler_
        // int dim = arrayStartDim + samples1DArraySizes.size();
        // for (size_t i = 0; i < samples2DArraySizes.size(); ++i) {
        //     int nSamples = samples2DArraySizes[i] * samplesPerPixel;
        //     for (int j = 0; j < nSamples; ++j) {
        //         int64_t idx = GetIndexForSample(j);
        //         sampleArray2D[i][j].x = SampleDimension(idx, dim);
        //         sampleArray2D[i][j].y = SampleDimension(idx, dim + 1);
        //     }
        //     dim += 2;
        // }
        // CHECK_EQ(arrayEndDim, dim);
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
            base_scales: self.base_scales,
            base_exponents: self.base_exponents,
            sample_stride: self.sample_stride,
            mult_inverse: self.mult_inverse,
            pixel_for_offset: RwLock::new(pixel_for_offset),
            offset_for_current_pixel: RwLock::new(offset_for_current_pixel),
            dimension: self.dimension,
            interval_sample_index: self.interval_sample_index,
            array_start_dim: self.array_start_dim,
            array_end_dim: self.array_end_dim,
            current_pixel: self.current_pixel,
            current_pixel_sample_index: self.current_pixel_sample_index,
            samples_1d_array_sizes: self.samples_1d_array_sizes.iter().cloned().collect(),
            samples_2d_array_sizes: self.samples_2d_array_sizes.iter().cloned().collect(),
            samples_1d_array: self.samples_1d_array.iter().cloned().collect(),
            samples_2d_array: self.samples_2d_array.iter().cloned().collect(),
            array_1d_offset: self.array_1d_offset,
            array_2d_offset: self.array_2d_offset,
        }
    }
}
