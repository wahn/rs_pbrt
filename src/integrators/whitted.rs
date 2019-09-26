// pbrt
use crate::core::geometry::vec3_abs_dot_nrm;
use crate::core::geometry::{Bounds2i, Normal3f, Ray, Vector3f};
use crate::core::integrator::SamplerIntegrator;
use crate::core::interaction::{Interaction, InteractionCommon};
use crate::core::light::VisibilityTester;
use crate::core::material::TransportMode;
use crate::core::pbrt::{Float, Spectrum};
use crate::core::reflection::BxdfType;
use crate::core::sampler::Sampler;
use crate::core::scene::Scene;

// see whitted.h

/// Whitted’s ray-tracing algorithm
pub struct WhittedIntegrator {
    // inherited from SamplerIntegrator (see integrator.h)
    pixel_bounds: Bounds2i,
    // see whitted.h
    max_depth: u32,
}

impl WhittedIntegrator {
    pub fn new(max_depth: u32, pixel_bounds: Bounds2i) -> Self {
        WhittedIntegrator {
            pixel_bounds,
            max_depth,
        }
    }
}

impl SamplerIntegrator for WhittedIntegrator {
    fn preprocess(&mut self, _scene: &Scene, _sampler: &mut Box<dyn Sampler + Send + Sync>) {}
    fn li(
        &self,
        ray: &mut Ray,
        scene: &Scene,
        sampler: &mut Box<dyn Sampler + Send + Sync>,
        // arena: &mut Arena,
        depth: i32,
    ) -> Spectrum {
        let mut l: Spectrum = Spectrum::default();
        // find closest ray intersection or return background radiance
        if let Some(mut isect) = scene.intersect(ray) {
            // compute emitted and reflected light at ray intersection point

            // initialize common variables for Whitted integrator
            let n: Normal3f = isect.shading.n;
            let wo: Vector3f = isect.wo;

            // compute scattering functions for surface interaction
            let mode: TransportMode = TransportMode::Radiance;
            isect.compute_scattering_functions(ray /* arena, */, false, mode);
            // if (!isect.bsdf)
            if let Some(ref _bsdf) = isect.bsdf {
            } else {
                return self.li(&mut isect.spawn_ray(&ray.d), scene, sampler, depth);
            }
            // compute emitted light if ray hit an area light source
            l += isect.le(&wo);

            // add contribution of each light source
            for light in &scene.lights {
                let mut wi: Vector3f = Vector3f::default();
                let mut pdf: Float = 0.0 as Float;
                let mut visibility: VisibilityTester = VisibilityTester::default();
                let it_common: InteractionCommon = InteractionCommon {
                    p: isect.get_p(),
                    time: isect.get_time(),
                    p_error: isect.get_p_error(),
                    wo: isect.get_wo(),
                    n: isect.get_n(),
                    medium_interface: isect.get_medium_interface(),
                };
                let li: Spectrum = light.sample_li(
                    &it_common,
                    &sampler.get_2d(),
                    &mut wi,
                    &mut pdf,
                    &mut visibility,
                );
                if li.is_black() || pdf == 0.0 as Float {
                    continue;
                }
                if let Some(ref bsdf) = isect.bsdf {
                    let bsdf_flags: u8 = BxdfType::BsdfAll as u8;
                    let f: Spectrum = bsdf.f(&wo, &wi, bsdf_flags);
                    if !f.is_black() && visibility.unoccluded(scene) {
                        l += f * li * vec3_abs_dot_nrm(&wi, &n) / pdf;
                    }
                } else {
                    panic!("no isect.bsdf found");
                }
            }
            if depth as u32 + 1 < self.max_depth {
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
            return l;
        } else {
            for light in &scene.lights {
                l += light.le(ray);
            }
            return l;
        }
    }
    fn get_pixel_bounds(&self) -> Bounds2i {
        self.pixel_bounds
    }
}
