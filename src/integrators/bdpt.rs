// pbrt
use core::camera::{Camera, CameraSample};
use core::geometry::{Bounds2i, Point2f, Ray};
use core::pbrt::Spectrum;
use core::sampler::Sampler;

// see bdpt.h

/// Bidirectional Path Tracing (Global Illumination)
pub struct BDPTIntegrator {
    pub max_depth: u32,
    visualize_strategies: bool,
    visualize_weights: bool,
    pub pixel_bounds: Bounds2i,
    light_sample_strategy: String, // "power"
}

impl BDPTIntegrator {
    pub fn new(
        // TODO: sampler
        // TODO: camera
        max_depth: u32,
        visualize_strategies: bool,
        visualize_weights: bool,
        pixel_bounds: Bounds2i,
        light_sample_strategy: String,
    ) -> Self {
        BDPTIntegrator {
            max_depth: max_depth,
            visualize_strategies: visualize_strategies,
            visualize_weights: visualize_weights,
            pixel_bounds: pixel_bounds,
            light_sample_strategy: light_sample_strategy,
        }
    }
    pub fn get_light_sample_strategy(&self) -> String {
        self.light_sample_strategy.clone()
    }
}

pub fn generate_camera_subpath(
    sampler: &mut Box<Sampler + Send + Sync>,
    max_depth: u32,
    camera: &Box<Camera + Send + Sync>,
    p_film: &Point2f,
) -> u32 {
    if max_depth == 0 {
        return 0_u32;
    }
    // TODO: ProfilePhase _(Prof::BDPTGenerateSubpath);
    // sample initial ray for camera subpath
    let mut camera_sample: CameraSample = CameraSample::default();
    camera_sample.p_film = *p_film;
    camera_sample.time = sampler.get_1d();
    camera_sample.p_lens = sampler.get_2d();
    let mut ray: Ray = Ray::default();
    let beta: Spectrum = Spectrum::new(camera.generate_ray_differential(&camera_sample, &mut ray));
    // WORK
    0_u32
}
