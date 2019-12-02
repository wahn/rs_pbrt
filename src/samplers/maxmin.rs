// std
use std::sync::RwLock;
// pbrt
use crate::core::geometry::{Bounds2i, Point2f, Point2i, Vector2i};
use crate::core::lowdiscrepancy::C_MAX_MIN_DIST;
use crate::core::paramset::ParamSet;
use crate::core::pbrt::Float;
use crate::core::pbrt::{is_power_of_2, log_2_int_i64, round_up_pow2_64};
use crate::core::rng::Rng;
use crate::core::sampler::Sampler;

pub struct MaxMinDistSampler {
    pub samples_per_pixel: i64,
    pub c_pixel: [u32; 32],
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
    pub fn new(samples_per_pixel: i64, n_sampled_dimensions: i64) -> Self {
        let mut samples_per_pixel: i64 = samples_per_pixel;
        let c_index: i32 = log_2_int_i64(samples_per_pixel) as i32;
        if c_index as usize >= 17_usize {
            panic!("No more than {} samples per pixel are supported with MaxMinDistSampler.");
            //     Warning(
            //         "No more than %d samples per pixel are supported with "
            //             "MaxMinDistSampler. Rounding down.",
            //         (1 << int(sizeof(CMaxMinDist) / sizeof(CMaxMinDist[0]))) -
            //             1);
            //     spp =
            //         (1 << (sizeof(CMaxMinDist) / sizeof(CMaxMinDist[0]))) - 1;
        }
        if !is_power_of_2(samples_per_pixel) {
            samples_per_pixel = round_up_pow2_64(samples_per_pixel);
            println!(
                "WARNING: Non power-of-two sample count rounded up to {:?} for MaxMinDistSampler.",
                samples_per_pixel
            );
        }
        let c_index: i32 = log_2_int_i64(samples_per_pixel) as i32;
        assert!(c_index >= 0_i32 && c_index < 17);
        let mut mmds: MaxMinDistSampler = MaxMinDistSampler {
            samples_per_pixel: samples_per_pixel,
            c_pixel: C_MAX_MIN_DIST[c_index as usize],
            samples_1d: Vec::new(),
            samples_2d: Vec::new(),
            current_1d_dimension: 0_i32,
            current_2d_dimension: 0_i32,
            rng: Rng::default(),
            current_pixel: Point2i::default(),
            current_pixel_sample_index: 0_i64,
            samples_1d_array_sizes: Vec::new(),
            samples_2d_array_sizes: Vec::new(),
            sample_array_1d: Vec::new(),
            sample_array_2d: Vec::new(),
            array_1d_offset: 0_usize,
            array_2d_offset: 0_usize,
        };
        for _i in 0..n_sampled_dimensions {
            let additional_1d: Vec<Float> = vec![0.0; mmds.samples_per_pixel as usize];
            let additional_2d: Vec<Point2f> =
                vec![Point2f::default(); mmds.samples_per_pixel as usize];
            mmds.samples_1d.push(additional_1d);
            mmds.samples_2d.push(additional_2d);
        }
        mmds
    }
    pub fn clone_with_seed(&self, seed: u64) -> Box<Sampler> {
        let mut mmds = MaxMinDistSampler {
            samples_per_pixel: self.samples_per_pixel,
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
    pub fn create(params: &ParamSet) -> Box<Sampler> {
        let nsamp: i32 = params.find_one_int("pixelsamples", 16);
        let sd: i32 = params.find_one_int("dimensions", 4);
        // TODO: if (PbrtOptions.quickRender) nsamp = 1;
        Box::new(Sampler::MaxMinDist(MaxMinDistSampler::new(
            nsamp as i64,
            sd as i64,
        )))
    }
    // Sampler
    pub fn start_pixel(&mut self, p: &Point2i) {
        // ProfilePhase _(Prof::StartPixel);
        // Float invSPP = (Float)1 / samplesPerPixel;
        // for (int i = 0; i < samplesPerPixel; ++i)
        //     samples2D[0][i] = Point2f(i * invSPP, SampleGeneratorMatrix(CPixel, i));
        // Shuffle(&samples2D[0][0], samplesPerPixel, 1, rng);
        // // Generate remaining samples for _MaxMinDistSampler_
        // for (size_t i = 0; i < samples1D.size(); ++i)
        //     VanDerCorput(1, samplesPerPixel, &samples1D[i][0], rng);

        // for (size_t i = 1; i < samples2D.size(); ++i)
        //     Sobol2D(1, samplesPerPixel, &samples2D[i][0], rng);

        // for (size_t i = 0; i < samples1DArraySizes.size(); ++i) {
        //     int count = samples1DArraySizes[i];
        //     VanDerCorput(count, samplesPerPixel, &sampleArray1D[i][0], rng);
        // }

        // for (size_t i = 0; i < samples2DArraySizes.size(); ++i) {
        //     int count = samples2DArraySizes[i];
        //     Sobol2D(count, samplesPerPixel, &sampleArray2D[i][0], rng);
        // }
        // PixelSampler::StartPixel(p);
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
