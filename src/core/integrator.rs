//! Rendering an image of the scene is handled by an instance of a
//! class that implements the **Integrator** interface.

// std
use std;
use std::sync::Arc;
// pbrt
use crate::core::geometry::{vec3_abs_dot_nrm, vec3_dot_nrm};
use crate::core::geometry::{Bounds2i, Normal3f, Point2f, Ray, RayDifferential, Vector3f};
use crate::core::interaction::{Interaction, InteractionCommon, SurfaceInteraction};
use crate::core::light::is_delta_light;
use crate::core::light::{Light, VisibilityTester};
use crate::core::pbrt::{Float, Spectrum};
use crate::core::reflection::BxdfType;
use crate::core::sampler::Sampler;
use crate::core::sampling::power_heuristic;
use crate::core::sampling::Distribution1D;
use crate::core::scene::Scene;

// see integrator.h

pub trait SamplerIntegrator {
    // TODO: use Sampler trait
    fn preprocess(&mut self, scene: &Scene, sampler: &mut Box<dyn Sampler + Send + Sync>);
    fn render(scene: &Scene);
    /// Returns the incident radiance at the origin of a given
    /// ray. Uses the scene's intersect routine to calculate a
    /// **SurfaceInteraction** and spawns rays if necessary.
    fn li(
        &self,
        ray: &mut Ray,
        scene: &Scene,
        sampler: &mut Box<dyn Sampler + Send + Sync>,
        // arena: &mut Arena,
        depth: i32,
    ) -> Spectrum;
    fn get_pixel_bounds(&self) -> Bounds2i;
    fn specular_reflect(
        &self,
        ray: &Ray,
        isect: &SurfaceInteraction,
        scene: &Scene,
        sampler: &mut Box<dyn Sampler + Send + Sync>,
        // arena: &mut Arena,
        depth: i32,
    ) -> Spectrum {
        // compute specular reflection direction _wi_ and BSDF value
        let wo: Vector3f = isect.wo;
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
                &sampler.get_2d(),
                &mut pdf,
                bsdf_flags,
                &mut sampled_type,
            );
            if pdf > 0.0 as Float && !f.is_black() && vec3_abs_dot_nrm(&wi, &ns) != 0.0 as Float {
                // compute ray differential _rd_ for specular reflection
                let mut rd: Ray = isect.spawn_ray(&wi);
                if let Some(d) = ray.differential.iter().next() {
                    let dudx: Float = *isect.dudx.read().unwrap();
                    let dvdx: Float = *isect.dvdx.read().unwrap();
                    let dndx: Normal3f = isect.shading.dndu * dudx + isect.shading.dndv * dvdx;
                    let dudy: Float = *isect.dudy.read().unwrap();
                    let dvdy: Float = *isect.dvdy.read().unwrap();
                    let dndy: Normal3f = isect.shading.dndu * dudy + isect.shading.dndv * dvdy;
                    let dwodx: Vector3f = -d.rx_direction - wo;
                    let dwody: Vector3f = -d.ry_direction - wo;
                    let ddndx: Float = vec3_dot_nrm(&dwodx, &ns) + vec3_dot_nrm(&wo, &dndx);
                    let ddndy: Float = vec3_dot_nrm(&dwody, &ns) + vec3_dot_nrm(&wo, &dndy);
                    // compute differential reflected directions
                    let dpdx: Vector3f = *isect.dpdx.read().unwrap();
                    let dpdy: Vector3f = *isect.dpdy.read().unwrap();
                    let diff: RayDifferential = RayDifferential {
                        rx_origin: isect.p + dpdx,
                        ry_origin: isect.p + dpdy,
                        rx_direction: wi - dwodx
                            + Vector3f::from(dndx * vec3_dot_nrm(&wo, &ns) + ns * ddndx)
                                * 2.0 as Float,
                        ry_direction: wi - dwody
                            + Vector3f::from(dndy * vec3_dot_nrm(&wo, &ns) + ns * ddndy)
                                * 2.0 as Float,
                    };
                    rd.differential = Some(diff);
                }
                return f
                    * self.li(&mut rd, scene, sampler, depth + 1)
                    * Spectrum::new(vec3_abs_dot_nrm(&wi, &ns) / pdf);
            } else {
                Spectrum::new(0.0)
            }
        } else {
            Spectrum::new(0.0)
        }
    }
    fn specular_transmit(
        &self,
        ray: &Ray,
        isect: &SurfaceInteraction,
        scene: &Scene,
        sampler: &mut Box<dyn Sampler + Send + Sync>,
        // arena: &mut Arena,
        depth: i32,
    ) -> Spectrum {
        let wo: Vector3f = isect.wo;
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
                &sampler.get_2d(),
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
                    let dudx: Float = *isect.dudx.read().unwrap();
                    let dvdx: Float = *isect.dvdx.read().unwrap();
                    let dndx: Normal3f = isect.shading.dndu * dudx + isect.shading.dndv * dvdx;
                    let dudy: Float = *isect.dudy.read().unwrap();
                    let dvdy: Float = *isect.dvdy.read().unwrap();
                    let dndy: Normal3f = isect.shading.dndu * dudy + isect.shading.dndv * dvdy;
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
                    let dpdx: Vector3f = *isect.dpdx.read().unwrap();
                    let dpdy: Vector3f = *isect.dpdy.read().unwrap();
                    let diff: RayDifferential = RayDifferential {
                        rx_origin: isect.p + dpdx,
                        ry_origin: isect.p + dpdy,
                        rx_direction: wi + dwodx * eta - Vector3f::from(dndx * mu + ns * dmudx),
                        ry_direction: wi + dwody * eta - Vector3f::from(dndy * mu + ns * dmudy),
                    };
                    rd.differential = Some(diff);
                }
                return f
                    * self.li(&mut rd, scene, sampler, depth + 1)
                    * Spectrum::new(vec3_abs_dot_nrm(&wi, &ns) / pdf);
            } else {
                Spectrum::new(0.0)
            }
        } else {
            Spectrum::new(0.0)
        }
    }
}

// see integrator.cpp

/// Most basic direct lighting strategy.
pub fn uniform_sample_all_lights(
    it: &SurfaceInteraction,
    scene: &Scene,
    sampler: &mut Box<dyn Sampler + Send + Sync>,
    n_light_samples: &Vec<i32>,
    handle_media: bool,
) -> Spectrum {
    // TODO: ProfilePhase p(Prof::DirectLighting);
    let mut l: Spectrum = Spectrum::new(0.0);
    for j in 0..scene.lights.len() {
        // accumulate contribution of _j_th light to _L_
        let ref light = scene.lights[j];
        let n_samples = n_light_samples[j];
        let u_light_array: Vec<Point2f> = sampler.get_2d_array(n_samples);
        let u_scattering_array: Vec<Point2f> = sampler.get_2d_array(n_samples);
        if u_light_array.is_empty() || u_scattering_array.is_empty() {
            // use a single sample for illumination from _light_
            let u_light: Point2f = sampler.get_2d();
            let u_scattering: Point2f = sampler.get_2d();
            l += estimate_direct(
                it,
                &u_scattering,
                light.clone(),
                &u_light,
                scene,
                sampler,
                handle_media,
                false,
            );
        } else {
            // estimate direct lighting using sample arrays
            let mut ld: Spectrum = Spectrum::new(0.0);
            for k in 0..n_samples {
                ld += estimate_direct(
                    it,
                    &u_scattering_array[k as usize],
                    light.clone(),
                    &u_light_array[k as usize],
                    scene,
                    sampler,
                    handle_media,
                    false,
                );
            }
            l += ld / n_samples as Float;
        }
    }
    l
}

/// Estimate direct lighting for only one randomly chosen light and
/// multiply the result by the number of lights to compensate.
pub fn uniform_sample_one_light(
    it: &dyn Interaction,
    scene: &Scene,
    sampler: &mut Box<dyn Sampler + Send + Sync>,
    handle_media: bool,
    light_distrib: Option<&Distribution1D>,
) -> Spectrum {
    // TODO: ProfilePhase p(Prof::DirectLighting);

    // randomly choose a single light to sample, _light_
    let n_lights: usize = scene.lights.len();
    if n_lights == 0_usize {
        return Spectrum::default();
    }
    let light_num: usize;
    let mut light_pdf: Option<Float> = Some(0.0 as Float);
    let pdf: Float;
    if let Some(light_distribution) = light_distrib {
        // if !light_distrib.is_null() {
        light_num = light_distribution.sample_discrete(sampler.get_1d(), light_pdf.as_mut());
        pdf = light_pdf.unwrap();
        if pdf == 0.0 as Float {
            return Spectrum::default();
        }
    } else {
        light_num = std::cmp::min(
            (sampler.get_1d() * n_lights as Float) as usize,
            n_lights - 1,
        );
        pdf = 1.0 as Float / n_lights as Float;
    }
    let light = &scene.lights[light_num];
    let u_light: Point2f = sampler.get_2d();
    let u_scattering: Point2f = sampler.get_2d();
    estimate_direct(
        it,
        &u_scattering,
        light.clone(),
        &u_light,
        scene,
        sampler,
        handle_media,
        false,
    ) / pdf
}

/// Computes a direct lighting estimate for a single light source sample.
pub fn estimate_direct(
    it: &dyn Interaction,
    u_scattering: &Point2f,
    light: Arc<Light>,
    u_light: &Point2f,
    scene: &Scene,
    sampler: &mut Box<dyn Sampler + Send + Sync>,
    // TODO: arena
    handle_media: bool,
    specular: bool,
) -> Spectrum {
    let mut bsdf_flags: u8 = BxdfType::BsdfAll as u8;
    if !specular {
        // bitwise not in Rust is ! (not the ~ operator like in C)
        bsdf_flags = BxdfType::BsdfAll as u8 & !(BxdfType::BsdfSpecular as u8);
    }
    let mut ld: Spectrum = Spectrum::new(0.0);
    // sample light source with multiple importance sampling
    let mut wi: Vector3f = Vector3f::default();
    let mut light_pdf: Float = 0.0 as Float;
    let mut scattering_pdf: Float = 0.0 as Float;
    let mut visibility: VisibilityTester = VisibilityTester::default();
    let it_common: InteractionCommon = InteractionCommon {
        p: it.get_p(),
        time: it.get_time(),
        p_error: it.get_p_error(),
        wo: it.get_wo(),
        n: it.get_n(),
        medium_interface: it.get_medium_interface(),
    };
    let mut li: Spectrum = light.sample_li(
        &it_common,
        u_light,
        &mut wi,
        &mut light_pdf,
        &mut visibility,
    );
    // TODO: println!("EstimateDirect uLight: {:?} -> Li: {:?}, wi:
    // {:?}, pdf: {:?}", u_light, li, wi, light_pdf);
    if light_pdf > 0.0 as Float && !li.is_black() {
        // compute BSDF or phase function's value for light sample
        let mut f: Spectrum = Spectrum::new(0.0);
        if it.is_surface_interaction() {
            // evaluate BSDF for light sampling strategy
            if let Some(ref bsdf) = it.get_bsdf() {
                if let Some(shading_n) = it.get_shading_n() {
                    f = bsdf.f(&it.get_wo(), &wi, bsdf_flags)
                        * Spectrum::new(vec3_abs_dot_nrm(&wi, &shading_n));
                    scattering_pdf = bsdf.pdf(&it.get_wo(), &wi, bsdf_flags);
                    // TODO: println!("  surf f*dot :{:?}, scatteringPdf: {:?}", f, scattering_pdf);
                }
            }
        } else {
            // evaluate phase function for light sampling strategy
            if let Some(ref phase) = it.get_phase() {
                let p: Float = phase.p(&it.get_wo(), &wi);
                f = Spectrum::new(p);
                scattering_pdf = p;
            }
        }
        if !f.is_black() {
            // compute effect of visibility for light source sample
            if handle_media {
                li *= visibility.tr(scene, sampler);
            // TODO: VLOG(2) << "  after Tr, Li: " << Li;
            } else {
                if !visibility.unoccluded(scene) {
                    // TODO: println!("  shadow ray blocked");
                    li = Spectrum::new(0.0 as Float);
                } else {
                    // TODO: println!("  shadow ray unoccluded");
                }
            }
            // add light's contribution to reflected radiance
            if !li.is_black() {
                if is_delta_light(light.get_flags()) {
                    ld += f * li / light_pdf;
                } else {
                    let weight: Float = power_heuristic(1_u8, light_pdf, 1_u8, scattering_pdf);
                    ld += f * li * Spectrum::new(weight) / light_pdf;
                }
            }
        }
    }
    // sample BSDF with multiple importance sampling
    if !is_delta_light(light.get_flags()) {
        let mut f: Spectrum = Spectrum::new(0.0);
        let mut sampled_specular: bool = false;
        if it.is_surface_interaction() {
            // sample scattered direction for surface interactions
            let mut sampled_type: u8 = 0_u8;
            if let Some(ref bsdf) = it.get_bsdf() {
                if let Some(shading_n) = it.get_shading_n() {
                    f = bsdf.sample_f(
                        &it.get_wo(),
                        &mut wi,
                        u_scattering,
                        &mut scattering_pdf,
                        bsdf_flags,
                        &mut sampled_type,
                    );
                    f *= Spectrum::new(vec3_abs_dot_nrm(&wi, &shading_n));
                    sampled_specular = (sampled_type & BxdfType::BsdfSpecular as u8) != 0_u8;
                }
            } else {
                println!("TODO: if let Some(ref bsdf) = it.get_bsdf() failed");
            }
        } else {
            // sample scattered direction for medium interactions
            if let Some(ref phase) = it.get_phase() {
                let p: Float = phase.sample_p(&it.get_wo(), &mut wi, u_scattering);
                f = Spectrum::new(p);
                scattering_pdf = p;
            }
        }
        // TODO: println!("  BSDF / phase sampling f: {:?}, scatteringPdf: {:?}",
        //          f, scattering_pdf);
        if !f.is_black() && scattering_pdf > 0.0 {
            // account for light contributions along sampled direction _wi_
            let mut weight: Float = 1.0;
            if !sampled_specular {
                light_pdf = light.pdf_li(it, wi);
                if light_pdf == 0.0 {
                    return ld;
                }
                weight = power_heuristic(1, scattering_pdf, 1, light_pdf);
            }
            // find intersection and compute transmittance
            let mut ray: Ray = it.spawn_ray(&wi);
            let mut tr: Spectrum = Spectrum::new(1.0 as Float);
            let mut found_surface_interaction: bool = false;
            // add light contribution from material sampling
            let mut li: Spectrum = Spectrum::default();
            if handle_media {
                let (light_isect_opt, tr_spectrum) = scene.intersect_tr(&mut ray, sampler);
                tr = tr_spectrum; // copy return value
                if let Some(ref light_isect) = light_isect_opt {
                    found_surface_interaction = true;
                    if let Some(primitive) = &light_isect.primitive {
                        if let Some(area_light) = primitive.get_area_light() {
                            let pa = &*area_light as *const _ as *const usize;
                            let pl = &*light as *const _ as *const usize;
                            if pa == pl {
                                li = light_isect.le(&-wi);
                            }
                        }
                    }
                }
            } else {
                if let Some(ref light_isect) = scene.intersect(&mut ray) {
                    found_surface_interaction = true;
                    if let Some(primitive) = &light_isect.primitive {
                        if let Some(area_light) = primitive.get_area_light() {
                            let pa = &*area_light as *const _ as *const usize;
                            let pl = &*light as *const _ as *const usize;
                            if pa == pl {
                                li = light_isect.le(&-wi);
                            }
                        }
                    }
                }
            }
            if !found_surface_interaction {
                li = light.le(&mut ray);
            }
            if !li.is_black() {
                ld += f * li * tr * weight / scattering_pdf;
            }
        }
    }
    ld
}

/// The light to start each photon path from is chosen according to a
/// PDF defined by the lights' respective powers.
pub fn compute_light_power_distribution(scene: &Scene) -> Option<Arc<Distribution1D>> {
    if scene.lights.is_empty() {
        return None;
    }
    let mut light_power: Vec<Float> = Vec::with_capacity(scene.lights.len());
    for li in 0..scene.lights.len() {
        let ref light = scene.lights[li];
        light_power.push(light.power().y());
    }
    Some(Arc::new(Distribution1D::new(light_power)))
}
