// pbrt
use crate::core::geometry::{Point2f, Point2i};
use crate::core::lowdiscrepancy::C_MAX_MIN_DIST;
use crate::core::lowdiscrepancy::{sample_generator_matrix, sobol_2d, van_der_corput};
use crate::core::paramset::ParamSet;
use crate::core::pbrt::Float;
use crate::core::pbrt::{is_power_of_2, log_2_int_i64, round_up_pow2_32, round_up_pow2_64};
use crate::core::rng::Rng;
use crate::core::sampler::Sampler;
use crate::core::sampling::shuffle;

pub struct MaxMinDistSampler {
    pub samples_per_pixel: i64,
    c_pixel: [u32; 32],
    // inherited from class PixelSampler (see sampler.h)
    samples_1d: Vec<Vec<Float>>,
    samples_2d: Vec<Vec<Point2f>>,
    current_1d_dimension: i32,
    current_2d_dimension: i32,
    rng: Rng,
    // inherited from class Sampler (see sampler.h)
    current_pixel: Point2i,
    current_pixel_sample_index: i64,
    samples_1d_array_sizes: Vec<i32>,
    samples_2d_array_sizes: Vec<i32>,
    sample_array_1d: Vec<Vec<Float>>,
    sample_array_2d: Vec<Vec<Point2f>>,
    array_1d_offset: usize,
    array_2d_offset: usize,
}

impl MaxMinDistSampler {
    pub fn new(samples_per_pixel: i64, n_sampled_dimensions: i64) -> Self {
        let mut samples_per_pixel: i64 = samples_per_pixel;
        let c_index: i32 = log_2_int_i64(samples_per_pixel) as i32;
        if c_index as usize >= 17_usize {
            panic!(
                "No more than {} samples per pixel are supported with MaxMinDistSampler.",
                c_index
            );
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
            samples_per_pixel,
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
            samples_1d: self.samples_1d.to_vec(),
            samples_2d: self.samples_2d.to_vec(),
            current_1d_dimension: self.current_1d_dimension,
            current_2d_dimension: self.current_2d_dimension,
            rng: self.rng,
            current_pixel: self.current_pixel,
            current_pixel_sample_index: self.current_pixel_sample_index,
            samples_1d_array_sizes: self.samples_1d_array_sizes.to_vec(),
            samples_2d_array_sizes: self.samples_2d_array_sizes.to_vec(),
            sample_array_1d: self.sample_array_1d.to_vec(),
            sample_array_2d: self.sample_array_2d.to_vec(),
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
    pub fn start_pixel(&mut self, p: Point2i) {
        // TODO: ProfilePhase _(Prof::StartPixel);
        let inv_spp: Float = 1.0 as Float / self.samples_per_pixel as Float;
        for i in 0..self.samples_per_pixel as usize {
            self.samples_2d[0_usize][i] = Point2f {
                x: i as Float * inv_spp,
                y: sample_generator_matrix(&self.c_pixel, i as u32, 0_u32),
            };
        }
        let samples: &mut [Point2f] = self.samples_2d[0].as_mut_slice();
        shuffle(samples, self.samples_per_pixel as i32, 1, &mut self.rng);
        // generate remaining samples for _MaxMinDistSampler_
        for i in 0..self.samples_1d.len() {
            let samples: &mut [Float] = self.samples_1d[i].as_mut_slice();
            van_der_corput(1, self.samples_per_pixel as i32, samples, &mut self.rng);
        }
        for i in 1..self.samples_2d.len() {
            let samples: &mut [Point2f] = self.samples_2d[i].as_mut_slice();
            sobol_2d(1, self.samples_per_pixel as i32, samples, &mut self.rng);
        }
        for i in 0..self.samples_1d_array_sizes.len() {
            let samples: &mut [Float] = self.sample_array_1d[i].as_mut_slice();
            van_der_corput(
                self.samples_1d_array_sizes[i],
                self.samples_per_pixel as i32,
                samples,
                &mut self.rng,
            );
        }
        for i in 0..self.samples_2d_array_sizes.len() {
            let samples: &mut [Point2f] = self.sample_array_2d[i].as_mut_slice();
            sobol_2d(
                self.samples_2d_array_sizes[i],
                self.samples_per_pixel as i32,
                samples,
                &mut self.rng,
            );
        }
        // PixelSampler::StartPixel(p);
        self.current_pixel = p;
        self.current_pixel_sample_index = 0_i64;
        // reset array offsets for next pixel sample
        self.array_1d_offset = 0_usize;
        self.array_2d_offset = 0_usize;
    }
    pub fn get_1d(&mut self) -> Float {
        // TODO: ProfilePhase _(Prof::GetSample);
        assert!(self.current_pixel_sample_index < self.samples_per_pixel);
        if self.current_1d_dimension < self.samples_1d.len() as i32 {
            let sample: Float = self.samples_1d[self.current_1d_dimension as usize]
                [self.current_pixel_sample_index as usize];
            self.current_1d_dimension += 1;
            sample
        } else {
            self.rng.uniform_float()
        }
    }
    pub fn get_2d(&mut self) -> Point2f {
        // TODO: ProfilePhase _(Prof::GetSample);
        assert!(self.current_pixel_sample_index < self.samples_per_pixel);
        if self.current_2d_dimension < self.samples_2d.len() as i32 {
            let sample: Point2f = self.samples_2d[self.current_2d_dimension as usize]
                [self.current_pixel_sample_index as usize];
            self.current_2d_dimension += 1;
            sample
        } else {
            // C++ call order for Point2f(rng.UniformFloat(), rng.UniformFloat());
            let y = self.rng.uniform_float();
            let x = self.rng.uniform_float();
            Point2f { x, y }
        }
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
        round_up_pow2_32(count)
    }
    pub fn get_2d_array(&mut self, n: i32) -> Option<&[Point2f]> {
        if self.array_2d_offset == self.sample_array_2d.len() {
            return None;
        }
        assert_eq!(self.samples_2d_array_sizes[self.array_2d_offset], n);
        assert!(self.current_pixel_sample_index < self.samples_per_pixel);
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
        assert!(self.current_pixel_sample_index < self.samples_per_pixel);
        let start: usize = (self.current_pixel_sample_index * n as i64) as usize;
        let end: usize = start + n as usize;
        self.array_2d_offset += 1;
        let ret1 = &self.sample_array_2d[self.array_2d_offset - 1][start..end];
        // repeat code from above
        if self.array_2d_offset == self.sample_array_2d.len() {
            return (None, None);
        }
        assert_eq!(self.samples_2d_array_sizes[self.array_2d_offset], n);
        assert!(self.current_pixel_sample_index < self.samples_per_pixel);
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
        assert!(self.current_pixel_sample_index < self.samples_per_pixel);
        let start: usize = (self.current_pixel_sample_index * n as i64) as usize;
        let idx: usize = self.array_2d_offset;
        self.array_2d_offset += 1;
        (false, idx, start)
    }
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
    pub fn reseed(&mut self, seed: u64) {
        self.rng.set_sequence(seed);
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
}
