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
        // TODO: ProfilePhase _(Prof::MediumSample);
        // sample a channel and distance along the ray
        let channel: usize = ((sampler.get_1d() * 3.0 as Float) as usize).min(2_usize);
        let dist: Float = (1.0 as Float - sampler.get_1d()).ln() / self.sigma_a[channel];
        let t: Float = (dist / ray.d.length()).min(ray.t_max);
        let sampled_medium: bool = t < ray.t_max;
        if sampled_medium {
            // TODO: *mi = MediumInteraction(...)
        }
        // WORK
        Spectrum::default()
    }
}
