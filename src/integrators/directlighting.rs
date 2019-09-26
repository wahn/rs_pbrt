// pbrt
use crate::core::geometry::{Bounds2i, Ray, Vector3f};
use crate::core::integrator::SamplerIntegrator;
use crate::core::integrator::{uniform_sample_all_lights, uniform_sample_one_light};
use crate::core::material::TransportMode;
use crate::core::pbrt::{Float, Spectrum};
use crate::core::sampler::Sampler;
use crate::core::scene::Scene;

// see directlighting.h

#[derive(Debug, Clone, PartialEq)]
pub enum LightStrategy {
    UniformSampleAll,
    UniformSampleOne,
}

/// Direct Lighting (no Global Illumination)
pub struct DirectLightingIntegrator {
    // inherited from SamplerIntegrator (see integrator.h)
    pixel_bounds: Bounds2i,
    // see directlighting.h
    strategy: LightStrategy,
    max_depth: u32,
    n_light_samples: Vec<i32>,
}

impl DirectLightingIntegrator {
    pub fn new(strategy: LightStrategy, max_depth: u32, pixel_bounds: Bounds2i) -> Self {
        DirectLightingIntegrator {
            pixel_bounds,
            strategy,
            max_depth,
            n_light_samples: Vec::new(),
        }
    }
}

impl SamplerIntegrator for DirectLightingIntegrator {
    fn preprocess(&mut self, scene: &Scene, sampler: &mut Box<dyn Sampler + Send + Sync>) {
        if self.strategy == LightStrategy::UniformSampleAll {
            // compute number of samples to use for each light
            for li in 0..scene.lights.len() {
                let ref light = scene.lights[li];
                self.n_light_samples
                    .push(sampler.round_count(light.get_n_samples()));
            }
            // request samples for sampling all lights
            for _i in 0..self.max_depth {
                for j in 0..scene.lights.len() {
                    sampler.request_2d_array(self.n_light_samples[j]);
                    sampler.request_2d_array(self.n_light_samples[j]);
                }
            }
        }
    }
    fn li(
        &self,
        ray: &mut Ray,
        scene: &Scene,
        sampler: &mut Box<dyn Sampler + Send + Sync>,
        // arena: &mut Arena,
        depth: i32,
    ) -> Spectrum {
        // TODO: ProfilePhase p(Prof::SamplerIntegratorLi);
        let mut l: Spectrum = Spectrum::new(0.0 as Float);
        // find closest ray intersection or return background radiance
        if let Some(mut isect) = scene.intersect(ray) {
            // compute scattering functions for surface interaction
            let mode: TransportMode = TransportMode::Radiance;
            isect.compute_scattering_functions(ray /* arena, */, false, mode);
            // if (!isect.bsdf)
            //     return Li(isect.SpawnRay(ray.d), scene, sampler, arena, depth);
            let wo: Vector3f = isect.wo;
            l += isect.le(&wo);
            if scene.lights.len() > 0 {
                // compute direct lighting for _DirectLightingIntegrator_ integrator
                if self.strategy == LightStrategy::UniformSampleAll {
                    l += uniform_sample_all_lights(
                        &isect,
                        scene,
                        sampler,
                        &self.n_light_samples,
                        false,
                    );
                } else {
                    l += uniform_sample_one_light(&isect, scene, sampler, false, None);
                }
            }
            if ((depth + 1_i32) as u32) < self.max_depth {
                // trace rays for specular reflection and refraction
                l += self.specular_reflect(
                    ray, &isect, scene, sampler, // arena,
                    depth,
                );
                l += self.specular_transmit(
                    ray, &isect, scene, sampler, // arena,
                    depth,
                );
            }
        } else {
            for light in &scene.lights {
                l += light.le(ray);
            }
        }
        l
    }
    fn get_pixel_bounds(&self) -> Bounds2i {
        self.pixel_bounds
    }
}
