// std
use std::borrow::Borrow;
use std::sync::Arc;
// pbrt
use crate::core::camera::Camera;
use crate::core::geometry::{vec3_abs_dot_nrm, vec3_dot_nrm};
use crate::core::geometry::{Bounds2i, Normal3f, Ray, RayDifferential, Vector3f};
use crate::core::integrator::{uniform_sample_all_lights, uniform_sample_one_light};
use crate::core::interaction::{Interaction, SurfaceInteraction};
use crate::core::material::TransportMode;
use crate::core::pbrt::{Float, Spectrum};
use crate::core::reflection::BxdfType;
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
    pub camera: Arc<Camera>,
    pub sampler: Box<Sampler>,
    pixel_bounds: Bounds2i,
    // see directlighting.h
    strategy: LightStrategy,
    max_depth: u32,
    n_light_samples: Vec<i32>,
}

impl DirectLightingIntegrator {
    pub fn new(
        strategy: LightStrategy,
        max_depth: u32,
        camera: Arc<Camera>,
        sampler: Box<Sampler>,
        pixel_bounds: Bounds2i,
    ) -> Self {
        DirectLightingIntegrator {
            camera,
            sampler,
            pixel_bounds,
            strategy,
            max_depth,
            n_light_samples: Vec::new(),
        }
    }
    pub fn preprocess(&mut self, scene: &Scene) {
        if self.strategy == LightStrategy::UniformSampleAll {
            // compute number of samples to use for each light
            for li in 0..scene.lights.len() {
                let light = &scene.lights[li];
                self.n_light_samples
                    .push(self.sampler.round_count(light.get_n_samples()));
            }
            // request samples for sampling all lights
            for _i in 0..self.max_depth {
                for j in 0..scene.lights.len() {
                    self.sampler.request_2d_array(self.n_light_samples[j]);
                    self.sampler.request_2d_array(self.n_light_samples[j]);
                }
            }
        }
    }
    pub fn li(
        &self,
        ray: &mut Ray,
        scene: &Scene,
        sampler: &mut Sampler,
        // arena: &mut Arena,
        depth: i32,
    ) -> Spectrum {
        // TODO: ProfilePhase p(Prof::SamplerIntegratorLi);
        let mut l: Spectrum = Spectrum::new(0.0 as Float);
        // find closest ray intersection or return background radiance
        let mut isect: SurfaceInteraction = SurfaceInteraction::default();
        if scene.intersect(ray, &mut isect) {
            // compute scattering functions for surface interaction
            let mode: TransportMode = TransportMode::Radiance;
            isect.compute_scattering_functions(ray, false, mode);
            // if (!isect.bsdf)
            //     return Li(isect.SpawnRay(ray.d), scene, sampler, arena, depth);
            let wo: Vector3f = isect.common.wo;
            l += isect.le(&wo);
            if !scene.lights.is_empty() {
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
                    let it: &SurfaceInteraction = isect.borrow();
                    l += uniform_sample_one_light(it, scene, sampler, false, None);
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
    pub fn get_camera(&self) -> Arc<Camera> {
        self.camera.clone()
    }
    pub fn get_sampler(&self) -> &Sampler {
        &self.sampler
    }
    pub fn get_pixel_bounds(&self) -> Bounds2i {
        self.pixel_bounds
    }
    pub fn specular_reflect(
        &self,
        ray: &Ray,
        isect: &SurfaceInteraction,
        scene: &Scene,
        sampler: &mut Sampler,
        // arena: &mut Arena,
        depth: i32,
    ) -> Spectrum {
        // compute specular reflection direction _wi_ and BSDF value
        let wo: Vector3f = isect.common.wo;
        let mut wi: Vector3f = Vector3f::default();
        let mut pdf: Float = 0.0 as Float;
        let ns: Normal3f = isect.shading.n;
        let mut sampled_type: u8 = 0_u8;
        let bsdf_flags: u8 = BxdfType::BsdfReflection as u8 | BxdfType::BsdfSpecular as u8;
        let f: Spectrum;
        if let Some(ref bsdf) = isect.bsdf {
            f = bsdf.sample_f(
                &wo,
                &mut wi,
                sampler.get_2d(),
                &mut pdf,
                bsdf_flags,
                &mut sampled_type,
            );
            if pdf > 0.0 as Float && !f.is_black() && vec3_abs_dot_nrm(&wi, &ns) != 0.0 as Float {
                // compute ray differential _rd_ for specular reflection
                let mut rd: Ray = isect.spawn_ray(&wi);
                if let Some(d) = ray.differential.iter().next() {
                    let dndx: Normal3f = isect.shading.dndu * isect.dudx.get()
                        + isect.shading.dndv * isect.dvdx.get();
                    let dndy: Normal3f = isect.shading.dndu * isect.dudy.get()
                        + isect.shading.dndv * isect.dvdy.get();
                    let dwodx: Vector3f = -d.rx_direction - wo;
                    let dwody: Vector3f = -d.ry_direction - wo;
                    let ddndx: Float = vec3_dot_nrm(&dwodx, &ns) + vec3_dot_nrm(&wo, &dndx);
                    let ddndy: Float = vec3_dot_nrm(&dwody, &ns) + vec3_dot_nrm(&wo, &dndy);
                    // compute differential reflected directions
                    let diff: RayDifferential = RayDifferential {
                        rx_origin: isect.common.p + isect.dpdx.get(),
                        ry_origin: isect.common.p + isect.dpdy.get(),
                        rx_direction: wi - dwodx
                            + Vector3f::from(dndx * vec3_dot_nrm(&wo, &ns) + ns * ddndx)
                                * 2.0 as Float,
                        ry_direction: wi - dwody
                            + Vector3f::from(dndy * vec3_dot_nrm(&wo, &ns) + ns * ddndy)
                                * 2.0 as Float,
                    };
                    rd.differential = Some(diff);
                }
                f * self.li(&mut rd, scene, sampler, depth + 1)
                    * Spectrum::new(vec3_abs_dot_nrm(&wi, &ns) / pdf)
            } else {
                Spectrum::new(0.0)
            }
        } else {
            Spectrum::new(0.0)
        }
    }
    pub fn specular_transmit(
        &self,
        ray: &Ray,
        isect: &SurfaceInteraction,
        scene: &Scene,
        sampler: &mut Sampler,
        // arena: &mut Arena,
        depth: i32,
    ) -> Spectrum {
        let wo: Vector3f = isect.common.wo;
        let mut wi: Vector3f = Vector3f::default();
        let mut pdf: Float = 0.0 as Float;
        // let p: Point3f = isect.p;
        let ns: Normal3f = isect.shading.n;
        let mut sampled_type: u8 = 0_u8;
        let bsdf_flags: u8 = BxdfType::BsdfTransmission as u8 | BxdfType::BsdfSpecular as u8;
        let f: Spectrum;
        if let Some(ref bsdf) = isect.bsdf {
            f = bsdf.sample_f(
                &wo,
                &mut wi,
                sampler.get_2d(),
                &mut pdf,
                bsdf_flags,
                &mut sampled_type,
            );
            if pdf > 0.0 as Float && !f.is_black() && vec3_abs_dot_nrm(&wi, &ns) != 0.0 as Float {
                // compute ray differential _rd_ for specular transmission
                let mut rd: Ray = isect.spawn_ray(&wi);
                if let Some(d) = ray.differential.iter().next() {
                    let mut eta: Float = bsdf.eta;
                    let w: Vector3f = -wo;
                    if vec3_dot_nrm(&wo, &ns) < 0.0 as Float {
                        eta = 1.0 / eta;
                    }
                    let dndx: Normal3f = isect.shading.dndu * isect.dudx.get()
                        + isect.shading.dndv * isect.dvdx.get();
                    let dndy: Normal3f = isect.shading.dndu * isect.dudy.get()
                        + isect.shading.dndv * isect.dvdy.get();
                    let dwodx: Vector3f = -d.rx_direction - wo;
                    let dwody: Vector3f = -d.ry_direction - wo;
                    let ddndx: Float = vec3_dot_nrm(&dwodx, &ns) + vec3_dot_nrm(&wo, &dndx);
                    let ddndy: Float = vec3_dot_nrm(&dwody, &ns) + vec3_dot_nrm(&wo, &dndy);
                    let mu: Float = eta * vec3_dot_nrm(&w, &ns) - vec3_dot_nrm(&wi, &ns);
                    let dmudx: Float = (eta
                        - (eta * eta * vec3_dot_nrm(&w, &ns)) / vec3_dot_nrm(&wi, &ns))
                        * ddndx;
                    let dmudy: Float = (eta
                        - (eta * eta * vec3_dot_nrm(&w, &ns)) / vec3_dot_nrm(&wi, &ns))
                        * ddndy;
                    let diff: RayDifferential = RayDifferential {
                        rx_origin: isect.common.p + isect.dpdx.get(),
                        ry_origin: isect.common.p + isect.dpdy.get(),
                        rx_direction: wi + dwodx * eta - Vector3f::from(dndx * mu + ns * dmudx),
                        ry_direction: wi + dwody * eta - Vector3f::from(dndy * mu + ns * dmudy),
                    };
                    rd.differential = Some(diff);
                }
                f * self.li(&mut rd, scene, sampler, depth + 1)
                    * Spectrum::new(vec3_abs_dot_nrm(&wi, &ns) / pdf)
            } else {
                Spectrum::new(0.0)
            }
        } else {
            Spectrum::new(0.0)
        }
    }
}
