// std
use std::f32;
use std::sync::Arc;
// others
use strum::IntoEnumIterator;
// pbrt
use crate::core::geometry::Ray;
use crate::core::interaction::MediumInteraction;
use crate::core::medium::{HenyeyGreenstein, Medium};
use crate::core::pbrt::{Float, Spectrum};
use crate::core::sampler::Sampler;
use crate::core::spectrum::RGBEnum;

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
            g,
        }
    }
    // Medium
    pub fn tr(&self, ray: &Ray, _sampler: &mut Sampler) -> Spectrum {
        // TODO: ProfilePhase _(Prof::MediumTr);
        (-self.sigma_t * (ray.t_max * ray.d.length()).min(f32::MAX)).exp()
    }
    pub fn sample(
        &self,
        ray: &Ray,
        sampler: &mut Sampler,
    ) -> (Spectrum, Option<MediumInteraction>) {
        // TODO: ProfilePhase _(Prof::MediumSample);
        // sample a channel and distance along the ray
        let channel: usize = ((sampler.get_1d() * 3.0 as Float) as usize).min(2_usize);
        let channel_rgb: RGBEnum = match channel {
            0 => RGBEnum::Red,
            1 => RGBEnum::Green,
            _ => RGBEnum::Blue,
        };
        let dist: Float = -((1.0 as Float - sampler.get_1d()).ln()) / self.sigma_t[channel_rgb];
        let t: Float = (dist / ray.d.length()).min(ray.t_max);
        let sampled_medium: bool = t < ray.t_max;
        let mi_opt = if sampled_medium {
            let mi: MediumInteraction = MediumInteraction::new(
                &ray.position(t),
                &(-ray.d),
                ray.time,
                Some(Arc::new(Medium::Homogeneous(HomogeneousMedium::new(
                    &self.sigma_a,
                    &self.sigma_s,
                    self.g,
                )))),
                Some(Arc::new(HenyeyGreenstein { g: self.g })),
            );
            Some(mi)
        } else {
            None
        };
        // compute the transmittance and sampling density
        let tr: Spectrum = (-self.sigma_t * t.min(f32::MAX) * ray.d.length()).exp();
        let density = if sampled_medium {
            self.sigma_t * tr
        } else {
            tr
        };
        let mut pdf: Float = 0.0 as Float;
        for i in RGBEnum::iter() {
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
