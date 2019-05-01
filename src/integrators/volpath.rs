// std
use std::sync::Arc;
// pbrt
use core::geometry::{Bounds2i, Point2f, Ray, Vector3f};
use core::integrator::SamplerIntegrator;
use core::lightdistrib::create_light_sample_distribution;
use core::lightdistrib::LightDistribution;
use core::pbrt::{Float, Spectrum};
use core::sampler::Sampler;
use core::scene::Scene;

// see volpath.h

/// Accounts for scattering and attenuation from participating media
/// as well as scattering from surfaces
pub struct VolPathIntegrator {
    // inherited from SamplerIntegrator (see integrator.h)
    pixel_bounds: Bounds2i,
    // see volpath.h
    pub max_depth: u32,
    rr_threshold: Float,           // 1.0
    light_sample_strategy: String, // "spatial"
    light_distribution: Option<Arc<LightDistribution + Send + Sync>>,
}

impl VolPathIntegrator {
    pub fn new(
        max_depth: u32,
        pixel_bounds: Bounds2i,
        rr_threshold: Float,
        light_sample_strategy: String,
    ) -> Self {
        VolPathIntegrator {
            pixel_bounds: pixel_bounds,
            max_depth: max_depth,
            rr_threshold: rr_threshold,
            light_sample_strategy: light_sample_strategy,
            light_distribution: None,
        }
    }
}

impl SamplerIntegrator for VolPathIntegrator {
    fn preprocess(&mut self, scene: &Scene, _sampler: &mut Box<Sampler + Send + Sync>) {
        self.light_distribution =
            create_light_sample_distribution(self.light_sample_strategy.clone(), scene);
    }
    fn li(
        &self,
        r: &mut Ray,
        scene: &Scene,
        sampler: &mut Box<Sampler + Send + Sync>,
        // arena: &mut Arena,
        _depth: i32,
    ) -> Spectrum {
        // WORK
        Spectrum::default()
    }
    fn get_pixel_bounds(&self) -> Bounds2i {
        self.pixel_bounds
    }
}
