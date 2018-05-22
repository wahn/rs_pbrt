// std
use std::sync::Arc;
// pbrt
use core::camera::Camera;
use core::geometry::{Bounds2f, Bounds2i, Point2f, Point2i};
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
        self.stream_index = index;
        self.sample_index = 0;
    }
    pub fn get_next_index(&mut self) -> i32 {
        self.sample_index += 1;
        self.stream_index + self.stream_count * self.sample_index
    }
    // private
    fn ensure_ready(&mut self, index: i32) {
    }
}

impl Sampler for MLTSampler {
    fn start_pixel(&mut self, p: &Point2i) {}
    fn get_1d(&mut self) -> Float {
        // TODO: ProfilePhase _(Prof::GetSample);
        let index: i32 = self.get_next_index();
        self.ensure_ready(index);
        // return X[index].value;
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
        sampler: &mut MLTSampler,
        depth: u32,
        p_raster: &mut Point2f,
    ) -> Spectrum {
        sampler.start_stream(CAMERA_STREAM_INDEX as i32);
        // determine the number of available strategies and pick a specific one
        let mut s: u32 = 0;
        let mut t: u32 = 0;
        let mut n_strategies: u32 = 0;
        if (depth == 0_u32) {
            n_strategies = 1;
            s = 0;
            t = 2;
        } else {
            // n_strategies = depth as u32 + 2;
            // s = std::min((int)(sampler.Get1D() * n_strategies), n_strategies - 1);
            // t = n_strategies - s;
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
        *p_raster = sample_bounds_f.lerp(&sampler.get_2d());
        let mut n_camera;
        {
            let mut sampler: Box<Sampler + Sync + Send> = Box::new(sampler.clone());
            let (n_camera_new, p_new, time_new) = generate_camera_subpath(
                scene,
                &mut sampler,
                t,
                &self.camera,
                p_raster,
                &mut camera_vertices,
            );
            n_camera = n_camera_new;
        }
        if n_camera != t as usize {
            return Spectrum::default();
        }
        // generate a light subpath with exactly _s_ vertices
        sampler.start_stream(LIGHT_STREAM_INDEX as i32);
        let mut light_vertices: Vec<Vertex> = Vec::with_capacity(s as usize);
        let mut n_light;
        {
            let mut sampler: Box<Sampler + Sync + Send> = Box::new(sampler.clone());
            let time: Float = camera_vertices[0].time();
            n_light = generate_light_subpath(
                scene,
                &mut sampler,
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
        sampler.start_stream(CONNECTION_STREAM_INDEX as i32);
        let mut sampler: Box<Sampler + Sync + Send> = Box::new(sampler.clone());
        connect_bdpt(
            scene,
            &light_vertices,
            &camera_vertices,
            s as usize,
            t as usize,
            &light_distr,
            // light_to_index,
            &self.camera,
            &mut sampler,
            p_raster,
            None,
        ) * (n_strategies as Float)
    }
}
