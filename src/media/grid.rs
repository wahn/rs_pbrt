// std
use std::sync::Arc;
// pbrt
use crate::core::geometry::pnt3i_inside_exclusive;
use crate::core::geometry::{Bounds3f, Bounds3i, Point3f, Point3i, Ray, Vector3f, Vector3i};
use crate::core::interaction::MediumInteraction;
use crate::core::medium::{HenyeyGreenstein, Medium};
use crate::core::pbrt::lerp;
use crate::core::pbrt::{Float, Spectrum};
use crate::core::sampler::Sampler;
use crate::core::spectrum::RGBEnum;
use crate::core::transform::Transform;

// see grid.h

pub struct GridDensityMedium {
    pub sigma_a: Spectrum,
    pub sigma_s: Spectrum,
    pub g: Float,
    pub nx: i32,
    pub ny: i32,
    pub nz: i32,
    pub world_to_medium: Transform,
    pub density: Arc<Vec<Float>>,
    pub sigma_t: Float,
    pub inv_max_density: Float,
}

impl GridDensityMedium {
    pub fn new(
        sigma_a: &Spectrum,
        sigma_s: &Spectrum,
        g: Float,
        nx: i32,
        ny: i32,
        nz: i32,
        medium_to_world: &Transform,
        d: Arc<Vec<Float>>,
    ) -> Self {
        let mut max_density: Float = 0.0;
        for i in 0..(nx * ny * nz) as usize {
            max_density = max_density.max(d[i]);
        }
        GridDensityMedium {
            sigma_a: *sigma_a,
            sigma_s: *sigma_s,
            g,
            nx,
            ny,
            nz,
            world_to_medium: Transform::inverse(medium_to_world),
            density: d,
            sigma_t: (*sigma_s + *sigma_a)[RGBEnum::Red],
            inv_max_density: 1.0 as Float / max_density,
        }
    }
    pub fn d(&self, p: &Point3i) -> Float {
        let sample_bounds: Bounds3i = Bounds3i {
            p_min: Point3i {
                x: 0_i32,
                y: 0_i32,
                z: 0_i32,
            },
            p_max: Point3i {
                x: self.nx,
                y: self.ny,
                z: self.nz,
            },
        };
        if !pnt3i_inside_exclusive(p, &sample_bounds) {
            0.0 as Float
        } else {
            self.density[((p.z * self.ny + p.y) * self.nx + p.x) as usize]
        }
    }
    pub fn density(&self, p: &Point3f) -> Float {
        // compute voxel coordinates and offsets for _p_
        let p_samples: Point3f = Point3f {
            x: p.x * self.nx as Float - 0.5 as Float,
            y: p.y * self.ny as Float - 0.5 as Float,
            z: p.z * self.nz as Float - 0.5 as Float,
        };
        let pi: Point3i = Point3i {
            x: p_samples.x.floor() as i32,
            y: p_samples.y.floor() as i32,
            z: p_samples.z.floor() as i32,
        };
        // Vector3f d = p_samples - (Point3f)pi;
        let d: Vector3f = Vector3f {
            x: p_samples.x - pi.x as Float,
            y: p_samples.y - pi.y as Float,
            z: p_samples.z - pi.z as Float,
        };
        // trilinearly interpolate density values to compute local density
        let d00: Float = lerp(
            d.x,
            self.d(&pi),
            self.d(&(pi
                + Vector3i {
                    x: 1_i32,
                    y: 0_i32,
                    z: 0_i32,
                })),
        );
        let d10: Float = lerp(
            d.x,
            self.d(&(pi
                + Vector3i {
                    x: 0_i32,
                    y: 1_i32,
                    z: 0_i32,
                })),
            self.d(&(pi
                + Vector3i {
                    x: 1_i32,
                    y: 1_i32,
                    z: 0_i32,
                })),
        );
        let d01: Float = lerp(
            d.x,
            self.d(&(pi
                + Vector3i {
                    x: 0_i32,
                    y: 0_i32,
                    z: 1_i32,
                })),
            self.d(&(pi
                + Vector3i {
                    x: 1_i32,
                    y: 0_i32,
                    z: 1_i32,
                })),
        );
        let d11: Float = lerp(
            d.x,
            self.d(&(pi
                + Vector3i {
                    x: 0_i32,
                    y: 1_i32,
                    z: 1_i32,
                })),
            self.d(&(pi
                + Vector3i {
                    x: 1_i32,
                    y: 1_i32,
                    z: 1_i32,
                })),
        );
        let d0: Float = lerp(d.y, d00, d10);
        let d1: Float = lerp(d.y, d01, d11);
        lerp(d.z, d0, d1)
    }
    // Medium
    pub fn tr(&self, r_world: &Ray, sampler: &mut Sampler) -> Spectrum {
        // TODO: ProfilePhase _(Prof::MediumTr);
        // TODO: ++nTrCalls;
        let mut in_ray: Ray = Ray::default();
        in_ray.o = r_world.o;
        in_ray.d = r_world.d.normalize();
        *in_ray.t_max.get_mut() = r_world.t_max.get() * r_world.d.length();
        let ray: Ray = self.world_to_medium.transform_ray(&in_ray);
        // compute $[\tmin, \tmax]$ interval of _ray_'s overlap with medium bounds
        let b: Bounds3f = Bounds3f::new(
            Point3f {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            Point3f {
                x: 1.0,
                y: 1.0,
                z: 1.0,
            },
        );
        let mut t_min: Float = 0.0;
        let mut t_max: Float = 0.0;
        // if (!b.IntersectP(ray, &tMin, &tMax)) return Spectrum(1.f);
        if !b.intersect_b(&ray, &mut t_min, &mut t_max) {
            return Spectrum::new(1.0 as Float);
        }
        // perform ratio tracking to estimate the transmittance value
        let mut tr: Float = 1.0;
        let mut t: Float = t_min;
        loop {
            // TODO: ++nTrSteps;
            t -= (1.0 as Float - sampler.get_1d()).ln() * self.inv_max_density / self.sigma_t;
            if t >= t_max {
                break;
            }
            let density: Float = self.density(&ray.position(t));
            tr *= 1.0 as Float - (0.0 as Float).max(density * self.inv_max_density);
            // added after book publication: when transmittance gets
            // low, start applying Russian roulette to terminate
            // sampling.
            let rr_threshold: Float = 0.1;
            if tr < rr_threshold {
                let q: Float = (0.05 as Float).max(1.0 as Float - tr);
                if sampler.get_1d() < q {
                    return Spectrum::default();
                }
                tr /= 1.0 as Float - q;
            }
        }
        Spectrum::new(tr)
    }
    pub fn sample(
        &self,
        r_world: &Ray,
        sampler: &mut Sampler,
    ) -> (Spectrum, Option<MediumInteraction>) {
        // TODO: ProfilePhase _(Prof::MediumSample);
        let mut in_ray: Ray = Ray::default();
        in_ray.o = r_world.o;
        in_ray.d = r_world.d.normalize();
        *in_ray.t_max.get_mut() = r_world.t_max.get() * r_world.d.length();
        let ray: Ray = self.world_to_medium.transform_ray(&in_ray);
        // compute $[\tmin, \tmax]$ interval of _ray_'s overlap with medium bounds
        let b: Bounds3f = Bounds3f::new(
            Point3f {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            Point3f {
                x: 1.0,
                y: 1.0,
                z: 1.0,
            },
        );
        let mut t_min: Float = 0.0;
        let mut t_max: Float = 0.0;
        if !b.intersect_b(&ray, &mut t_min, &mut t_max) {
            return (Spectrum::new(1.0 as Float), None);
        }
        // run delta-tracking iterations to sample a medium interaction
        let mut t: Float = t_min;
        loop {
            t -= (1.0 as Float - sampler.get_1d()).ln() * self.inv_max_density / self.sigma_t;
            if t >= t_max {
                break;
            }
            if self.density(&ray.position(t)) * self.inv_max_density > sampler.get_1d() {
                let mi_opt: Option<MediumInteraction>;
                // populate _mi_ with medium interaction information and return
                let mi: MediumInteraction = MediumInteraction::new(
                    &r_world.position(t),
                    &(-r_world.d),
                    r_world.time,
                    Some(Arc::new(Medium::GridDensity(GridDensityMedium {
                        sigma_a: self.sigma_a,
                        sigma_s: self.sigma_s,
                        g: self.g,
                        nx: self.nx,
                        ny: self.ny,
                        nz: self.nz,
                        world_to_medium: self.world_to_medium,
                        density: self.density.clone(),
                        sigma_t: self.sigma_t,
                        inv_max_density: self.inv_max_density,
                    }))),
                    Some(Arc::new(HenyeyGreenstein { g: self.g })),
                );
                mi_opt = Some(mi);
                return (self.sigma_s / self.sigma_t, mi_opt);
            }
        }
        (Spectrum::new(1.0 as Float), None)
    }
}
