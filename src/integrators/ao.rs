// std
use std::sync::Arc;
// pbrt
use crate::core::camera::Camera;
use crate::core::geometry::{nrm_cross_vec3, nrm_faceforward_vec3, vec3_dot_nrm};
use crate::core::geometry::{Bounds2i, Normal3f, Point2f, Ray, Vector3f};
use crate::core::interaction::Interaction;
use crate::core::material::TransportMode;
use crate::core::pbrt::{Float, Spectrum};
use crate::core::sampler::Sampler;
use crate::core::sampling::{
    cosine_hemisphere_pdf, cosine_sample_hemisphere, uniform_hemisphere_pdf,
    uniform_sample_hemisphere,
};
use crate::core::scene::Scene;

// see ao.h

/// Ambient Occlusion
pub struct AOIntegrator {
    // inherited from SamplerIntegrator (see integrator.h)
    pub camera: Arc<Camera>,
    pub sampler: Box<dyn Sampler + Send + Sync>,
    pub pixel_bounds: Bounds2i,
    // see ao.h
    pub cos_sample: bool,
    pub n_samples: i32,
}

impl AOIntegrator {
    pub fn new(
        cos_sample: bool,
        n_samples: i32,
        camera: Arc<Camera>,
        sampler: Box<dyn Sampler + Send + Sync>,
        pixel_bounds: Bounds2i,
    ) -> Self {
        AOIntegrator {
            camera,
            sampler,
            pixel_bounds,
            cos_sample,
            n_samples,
        }
    }
    pub fn preprocess(&mut self, _scene: &Scene) {
        self.sampler.request_2d_array(self.n_samples);
    }
    pub fn li(
        &self,
        r: &mut Ray,
        scene: &Scene,
        sampler: &mut Box<dyn Sampler + Send + Sync>,
        // arena: &mut Arena,
        _depth: i32,
    ) -> Spectrum {
        // TODO: ProfilePhase p(Prof::SamplerIntegratorLi);
        let mut l: Spectrum = Spectrum::default();
        let mut ray: Ray = Ray {
            o: r.o,
            d: r.d,
            t_max: r.t_max,
            time: r.time,
            differential: r.differential,
            medium: r.medium.clone(),
        };
        if let Some(mut isect) = scene.intersect(&mut ray) {
            let mode: TransportMode = TransportMode::Radiance;
            isect.compute_scattering_functions(&mut ray, true, mode);
            // if (!isect.bsdf) {
            //     VLOG(2) << "Skipping intersection due to null bsdf";
            //     ray = isect.SpawnRay(ray.d);
            //     goto retry;
            // }
            // compute coordinate frame based on true geometry, not
            // shading geometry.
            let n: Normal3f = nrm_faceforward_vec3(&isect.n, &-ray.d);
            let s: Vector3f = isect.dpdu.normalize();
            let t: Vector3f = nrm_cross_vec3(&isect.n, &s);
            let u_opt: Option<&[Point2f]> = sampler.get_2d_array(self.n_samples);
            if let Some(u) = u_opt {
                for i in 0..self.n_samples as usize {
                    // Vector3f wi;
                    let mut wi: Vector3f;
                    let pdf: Float;
                    if self.cos_sample {
                        wi = cosine_sample_hemisphere(&u[i]);
                        pdf = cosine_hemisphere_pdf(wi.z.abs());
                    } else {
                        wi = uniform_sample_hemisphere(&u[i]);
                        pdf = uniform_hemisphere_pdf();
                    }
                    // transform wi from local frame to world space.
                    wi = Vector3f {
                        x: s.x * wi.x + t.x * wi.y + n.x * wi.z,
                        y: s.y * wi.x + t.y * wi.y + n.y * wi.z,
                        z: s.z * wi.x + t.z * wi.y + n.z * wi.z,
                    };
                    let mut ray: Ray = isect.spawn_ray(&wi);
                    if !scene.intersect_p(&mut ray) {
                        l += Spectrum::new(vec3_dot_nrm(&wi, &n) / (pdf * self.n_samples as Float));
                    }
                }
            }
        }
        l
    }
    pub fn get_camera(&self) -> Arc<Camera> {
        self.camera.clone()
    }
    pub fn get_sampler(&self) -> &Box<dyn Sampler + Send + Sync> {
        &self.sampler
    }
    pub fn get_pixel_bounds(&self) -> Bounds2i {
        self.pixel_bounds
    }
}
