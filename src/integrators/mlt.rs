// std
use std::sync::Arc;
// pbrt
use core::camera::Camera;
use core::geometry::{Bounds2f, Bounds2i, Point2f, Point2i};
use core::pbrt::erf_inv;
use core::pbrt::SQRT_2;
use core::pbrt::{Float, Spectrum};
use core::rng::Rng;
use core::sampler::Sampler;
use core::sampling::Distribution1D;
use core::scene::Scene;
use integrators::bdpt::Vertex;
use integrators::bdpt::{connect_bdpt, generate_camera_subpath, generate_light_subpath};

pub const CAMERA_STREAM_INDEX: u8 = 0;
pub const LIGHT_STREAM_INDEX: u8 = 1;
pub const CONNECTION_STREAM_INDEX: u8 = 2;
pub const N_SAMPLE_STREAMS: u8 = 3;

#[derive(Debug, Default, Copy, Clone)]
pub struct PrimarySample {
    pub value: Float,
    pub last_modification_iteration: i64,
    pub value_backup: Float,
    pub modify_backup: i64,
}

impl PrimarySample {
    pub fn backup(&mut self) {
        self.value_backup = self.value;
        self.modify_backup = self.last_modification_iteration;
    }
    pub fn restore(&mut self) {
        self.value = self.value_backup;
        self.last_modification_iteration = self.modify_backup;
    }
}

pub struct MLTSampler {
    pub samples_per_pixel: i64,
    pub rng: Rng,
    pub sigma: Float,
    pub large_step_probability: Float,
    pub stream_count: i32,
    pub x: Vec<PrimarySample>,
    pub current_iteration: i64,
    pub large_step: bool,
    pub last_large_step_iteration: i64,
    pub stream_index: i32,
    pub sample_index: i32,
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
            x: Vec::new(),
            current_iteration: 0_i64,
            large_step: true,
            last_large_step_iteration: 0_i64,
            stream_index: 0_i32,
            sample_index: 0_i32,
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
    pub fn start_iteration(&mut self) {
        self.current_iteration += 1;
        self.large_step = self.rng.uniform_float() < self.large_step_probability;
    }
    pub fn accept(&mut self) {
        if self.large_step {
            self.last_large_step_iteration = self.current_iteration;
        }
    }
    pub fn reject(&mut self) {
        for mut i in 0..self.x.len() {
            let mut xi: PrimarySample = self.x[i];
            if xi.last_modification_iteration == self.current_iteration {
                xi.restore();
            }
        }
        self.current_iteration -= 1;
    }
    pub fn start_stream(&mut self, index: i32) {
        assert!(index < self.stream_count);
        self.stream_index = index;
        self.sample_index = 0;
    }
    pub fn get_next_index(&mut self) -> i32 {
        let ret = self.stream_index + self.stream_count * self.sample_index;
        self.sample_index += 1;
        ret
    }
    // private
    fn ensure_ready(&mut self, index: i32) {
        // enlarge _MLTSampler::x_ if necessary and get current $\VEC{X}_i$
        if index as usize >= self.x.len() {
            self.x
                .resize((index + 1) as usize, PrimarySample::default());
        }
        if let Some(xi) = self.x.get_mut(index as usize) {
            // reset $\VEC{X}_i$ if a large step took place in the meantime
            if xi.last_modification_iteration < self.last_large_step_iteration {
                xi.value = self.rng.uniform_float();
                xi.last_modification_iteration = self.last_large_step_iteration;
            }
            // apply remaining sequence of mutations to _sample_
            xi.backup();
            if self.large_step {
                xi.value = self.rng.uniform_float();
            } else {
                let n_small: i64 = self.current_iteration - xi.last_modification_iteration;
                // apply _n_small_ small step mutations

                // sample the standard normal distribution $N(0, 1)$
                let normal_sample: Float =
                    SQRT_2 * erf_inv(2.0 as Float * self.rng.uniform_float() - 1.0 as Float);
                // compute the effective standard deviation and apply perturbation to
                // $\VEC{X}_i$
                let eff_sigma: Float = self.sigma * (n_small as Float).sqrt();
                xi.value += normal_sample * eff_sigma;
                xi.value -= xi.value.floor();
            }
            xi.last_modification_iteration = self.current_iteration;
        } else {
            panic!("self.x.get_mut({:?}) failed", index);
        }
    }
}

impl Sampler for MLTSampler {
    fn start_pixel(&mut self, p: &Point2i) {
        // Sampler::StartPixel(p);
        self.current_pixel = *p;
        self.current_pixel_sample_index = 0_i64;
        self.array_1d_offset = 0_usize;
        self.array_2d_offset = 0_usize;
    }
    fn get_1d(&mut self) -> Float {
        // TODO: ProfilePhase _(Prof::GetSample);
        let index: i32 = self.get_next_index();
        self.ensure_ready(index);
        self.x[index as usize].value
    }
    fn get_2d(&mut self) -> Point2f {
        let x: Float = self.get_1d();
        let y: Float = self.get_1d();
        Point2f { x: x, y: y }
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
    fn get_2d_array(&mut self, n: i32) -> Vec<Point2f> {
        let mut samples: Vec<Point2f> = Vec::new();
        if self.array_2d_offset == self.sample_array_2d.len() {
            return samples;
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
        samples = self.sample_array_2d[self.array_2d_offset][start..end].to_vec();
        self.array_2d_offset += 1;
        samples
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

impl Clone for MLTSampler {
    fn clone(&self) -> MLTSampler {
        MLTSampler {
            samples_per_pixel: self.samples_per_pixel,
            rng: self.rng.clone(),
            sigma: self.sigma,
            large_step_probability: self.large_step_probability,
            stream_count: self.stream_count,
            x: self.x.iter().cloned().collect(),
            current_iteration: self.current_iteration,
            large_step: self.large_step,
            last_large_step_iteration: self.last_large_step_iteration,
            stream_index: self.stream_index,
            sample_index: self.sample_index,
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
    pub camera: Arc<Camera + Sync + Send>,
    pub max_depth: u32,
    pub n_bootstrap: u32,
    pub n_chains: u32,
    pub mutations_per_pixel: u32,
    pub sigma: Float,
    pub large_step_probability: Float,
}

impl MLTIntegrator {
    pub fn new(
        camera: Arc<Camera + Sync + Send>,
        max_depth: u32,
        n_bootstrap: u32,
        n_chains: u32,
        mutations_per_pixel: u32,
        sigma: Float,
        large_step_probability: Float,
    ) -> Self {
        MLTIntegrator {
            camera: camera,
            max_depth: max_depth,
            n_bootstrap: n_bootstrap,
            n_chains: n_chains,
            mutations_per_pixel: mutations_per_pixel,
            sigma: sigma,
            large_step_probability: large_step_probability,
        }
    }
    pub fn l(
        &self,
        scene: &Scene,
        light_distr: &Arc<Distribution1D>,
        mlt_sampler: &mut MLTSampler,
        depth: u32,
        p_raster: &mut Point2f,
    ) -> Spectrum {
        mlt_sampler.start_stream(CAMERA_STREAM_INDEX as i32);
        // determine the number of available strategies and pick a specific one
        let mut s: u32 = 0;
        let mut t: u32 = 0;
        let mut n_strategies: u32 = 0;
        if depth == 0_u32 {
            n_strategies = 1;
            s = 0;
            t = 2;
        } else {
            n_strategies = depth + 2;
            s = ((mlt_sampler.get_1d() * n_strategies as Float) as u32).min(n_strategies - 1);
            t = n_strategies - s;
        }
        // generate a camera subpath with exactly _t_ vertices
        let mut camera_vertices: Vec<Vertex> = Vec::with_capacity(t as usize);
        let film = self.camera.get_film();
        let sample_bounds: Bounds2i = film.get_sample_bounds();
        let sample_bounds_f: Bounds2f = Bounds2f {
            p_min: Point2f {
                x: sample_bounds.p_min.x as Float,
                y: sample_bounds.p_min.y as Float,
            },
            p_max: Point2f {
                x: sample_bounds.p_max.x as Float,
                y: sample_bounds.p_max.y as Float,
            },
        };
        *p_raster = sample_bounds_f.lerp(&mlt_sampler.get_2d());
        let mut n_camera;
        let mut p;
        let mut time;
        {
            let (n_camera_new, p_new, time_new) = generate_camera_subpath(
                scene,
                mlt_sampler as &mut Sampler,
                t,
                &self.camera,
                p_raster,
                &mut camera_vertices,
            );
            n_camera = n_camera_new;
            p = p_new;
            time = time_new;
        }
        if n_camera != t as usize {
            return Spectrum::default();
        }
        // generate a light subpath with exactly _s_ vertices
        mlt_sampler.start_stream(LIGHT_STREAM_INDEX as i32);
        let mut light_vertices: Vec<Vertex> = Vec::with_capacity(s as usize);
        let mut n_light;
        {
            n_light = generate_light_subpath(
                scene,
                mlt_sampler as &mut Sampler,
                s,
                time,
                &light_distr,
                // light_to_index,
                &mut light_vertices,
            );
        }
        if n_light != s as usize {
            return Spectrum::default();
        }
        // execute connection strategy and return the radiance estimate
        mlt_sampler.start_stream(CONNECTION_STREAM_INDEX as i32);
        connect_bdpt(
            scene,
            &light_vertices,
            &camera_vertices,
            s as usize,
            t as usize,
            &light_distr,
            // light_to_index,
            &self.camera,
            mlt_sampler as &mut Sampler,
            p_raster,
            None,
        ) * (n_strategies as Float)
    }
}
