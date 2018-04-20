// std
use std::f32;
use std::sync::Arc;
// pbrt
use core::geometry::Ray;
use core::interaction::MediumInteraction;
use core::medium::{HenyeyGreenstein, Medium};
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
        // TODO: ProfilePhase _(Prof::MediumTr);
        (-self.sigma_t * (ray.t_max * ray.d.length()).min(f32::MAX)).exp()
    }
    fn sample(
        &self,
        ray: &Ray,
        sampler: &mut Box<Sampler + Send + Sync>,
    ) -> (Spectrum, Option<MediumInteraction>) {
        // TODO: ProfilePhase _(Prof::MediumSample);
        // sample a channel and distance along the ray
        let channel: usize = ((sampler.get_1d() * 3.0 as Float) as usize).min(2_usize);
        let dist: Float = (1.0 as Float - sampler.get_1d()).ln() / self.sigma_a[channel];
        let t: Float = (dist / ray.d.length()).min(ray.t_max);
        let sampled_medium: bool = t < ray.t_max;
        let mut mi_opt: Option<MediumInteraction> = None;
        if sampled_medium {
            let mi: MediumInteraction = MediumInteraction::new(
                &ray.position(t),
                &(-ray.d),
                ray.time,
                Some(Arc::new(HomogeneousMedium::new(
                    &self.sigma_a,
                    &self.sigma_s,
                    self.g,
                ))),
                Some(Arc::new(HenyeyGreenstein { g: self.g })),
            );
            mi_opt = Some(mi);
        }
        // compute the transmittance and sampling density
        let tr: Spectrum = (-self.sigma_t * t.min(f32::MAX) * ray.d.length()).exp();
        let density: Spectrum;
        if sampled_medium {
            density = self.sigma_t * tr;
        } else {
            density = tr;
        }
        let mut pdf: Float = 0.0 as Float;
        for i in 0..3 {
            // TODO: Spectrum::nSamples
            pdf += density[i];
        }
        pdf *= 1.0 as Float / 3.0 as Float; // TODO: Spectrum::nSamples
        if pdf == 0.0 as Float {
            assert!(tr.is_black());
            pdf = 1.0 as Float;
        }
        if sampled_medium {
            (tr * self.sigma_s / pdf, mi_opt)
        } else {
            (tr / pdf, mi_opt)
        }
    }
}
