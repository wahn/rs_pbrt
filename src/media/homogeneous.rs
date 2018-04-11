// pbrt
use core::geometry::Ray;
use core::interaction::MediumInteraction;
use core::medium::Medium;
use core::pbrt::{Float, Spectrum};
use core::sampler::Sampler;

// see homogeneous.h

pub struct HomogeneousMedium {
    pub sigma_a: Spectrum,
    pub sigma_s: Spectrum,
    pub sigma_t: Spectrum,
    pub g: Float,
}

impl HomogeneousMedium {
    pub fn new(sigma_a: &Spectrum, sigma_s: &Spectrum, g: Float) -> Self {
        HomogeneousMedium {
            sigma_a: *sigma_a,
            sigma_s: *sigma_s,
            sigma_t: *sigma_s + *sigma_a,
            g: g,
        }
    }
}

impl Medium for HomogeneousMedium {
    fn tr(&self, ray: &Ray, sampler: &mut Box<Sampler + Send + Sync>) -> Spectrum {
        // WORK
        Spectrum::default()
    }
    fn sample(
        &self,
        ray: &Ray,
        sampler: &mut Box<Sampler + Send + Sync>,
        mi: &mut MediumInteraction,
    ) -> Spectrum {
        // WORK
        Spectrum::default()
    }
}
