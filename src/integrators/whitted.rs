// std
use std::sync::Arc;
// pbrt
use crate::core::camera::Camera;
use crate::core::geometry::{vec3_abs_dot_nrm, vec3_dot_nrm};
use crate::core::geometry::{Bounds2i, Normal3f, Ray, RayDifferential, Vector3f};
use crate::core::interaction::{Interaction, InteractionCommon, SurfaceInteraction};
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
    pub camera: Arc<Camera>,
    pub sampler: Box<Sampler>,
    pixel_bounds: Bounds2i,
    // see whitted.h
    max_depth: u32,
}

impl WhittedIntegrator {
    pub fn new(
        max_depth: u32,
        camera: Arc<Camera>,
        sampler: Box<Sampler>,
        pixel_bounds: Bounds2i,
    ) -> Self {
        WhittedIntegrator {
            camera,
            sampler,
            pixel_bounds,
            max_depth,
        }
    }
    pub fn preprocess(&mut self, _scene: &Scene) {}
    pub fn li(
        &self,
        ray: &mut Ray,
        scene: &Scene,
        sampler: &mut Sampler,
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
            l
        } else {
            for light in &scene.lights {
                l += light.le(ray);
            }
            l
        }
    }
    pub fn get_camera(&self) -> Arc<Camera> {
        self.camera.clone()
    }
    pub fn get_sampler(&self) -> &Box<Sampler> {
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
