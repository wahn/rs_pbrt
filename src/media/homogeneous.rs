// pbrt
use core::geometry::Ray;
use core::interaction::MediumInteraction;
use core::medium::Medium;
use core::pbrt::{Float, Spectrum};
use core::sampler::Sampler;

// see homogeneous.h

pub struct HomogeneousMedium {}

impl HomogeneousMedium {}

impl Medium for HomogeneousMedium {
    fn tr(ray: &Ray, sampler: &mut Box<Sampler + Send + Sync>) -> Spectrum {
        // WORK
        Spectrum::default()
    }
    fn sample(
        ray: &Ray,
        sampler: &mut Box<Sampler + Send + Sync>,
        mi: &mut MediumInteraction,
    ) -> Spectrum {
        // WORK
        Spectrum::default()
    }
}
