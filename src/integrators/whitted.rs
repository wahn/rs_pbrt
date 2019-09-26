// pbrt
use crate::core::geometry::{Bounds2i, Ray};
use crate::core::integrator::SamplerIntegrator;
use crate::core::pbrt::Spectrum;
use crate::core::sampler::Sampler;
use crate::core::scene::Scene;

// see whitted.h

pub struct WhittedIntegrator {
    // inherited from SamplerIntegrator (see integrator.h)
    pixel_bounds: Bounds2i,
}

impl WhittedIntegrator {}

impl SamplerIntegrator for WhittedIntegrator {
    fn preprocess(&mut self, _scene: &Scene, _sampler: &mut Box<dyn Sampler + Send + Sync>) {}
    fn li(
        &self,
        r: &mut Ray,
        scene: &Scene,
        sampler: &mut Box<dyn Sampler + Send + Sync>,
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
