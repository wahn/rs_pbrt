// pbrt
use core::geometry::{Bounds3f, Point3f, Ray};
use core::interaction::MediumInteraction;
use core::medium::Medium;
use core::pbrt::{Float, Spectrum};
use core::sampler::Sampler;
use core::transform::Transform;

// see grid.h

pub struct GridDensityMedium {
    pub sigma_a: Spectrum,
    pub sigma_s: Spectrum,
    pub g: Float,
    pub nx: i32,
    pub ny: i32,
    pub nz: i32,
    pub world_to_medium: Transform,
    pub density: Vec<Float>,
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
        world_to_medium: &Transform,
        d: Vec<Float>,
    ) -> Self {
        let mut max_density: Float = 0.0;
        for i in 0..(nx * ny * nz) as usize {
            max_density = max_density.max(d[i]);
        }
        GridDensityMedium {
            sigma_a: *sigma_a,
            sigma_s: *sigma_s,
            g: g,
            nx: nx,
            ny: ny,
            nz: nz,
            world_to_medium: *world_to_medium,
            density: d,
            sigma_t: (*sigma_s + *sigma_a)[0],
            inv_max_density: 1.0 as Float / max_density,
        }
    }
}

impl Medium for GridDensityMedium {
    fn tr(&self, ray: &Ray, _sampler: &mut Sampler) -> Spectrum {
        // TODO
        Spectrum::default()
    }
    fn sample(
        &self,
        r_world: &Ray,
        sampler: &mut Sampler,
    ) -> (Spectrum, Option<MediumInteraction>) {
        // TODO: ProfilePhase _(Prof::MediumSample);
        let mut in_ray: Ray = Ray::default();
        in_ray.o = r_world.o;
        in_ray.d = r_world.d.normalize();
        in_ray.t_max = r_world.t_max * r_world.d.length();
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
        // Float tMin, tMax;
        // if (!b.IntersectP(ray, &tMin, &tMax)) return Spectrum(1.f);
        // TODO: we need a new intersection routine (and a name for it)
        // b.intersect_p(&ray, );
        // // Run delta-tracking iterations to sample a medium interaction
        // Float t = tMin;
        // while (true) {
        //     t -= std::log(1 - sampler.Get1D()) * invMaxDensity / sigma_t;
        //     if (t >= tMax) break;
        //     if (Density(ray(t)) * invMaxDensity > sampler.Get1D()) {
        //         // Populate _mi_ with medium interaction information and return
        //         PhaseFunction *phase = ARENA_ALLOC(arena, HenyeyGreenstein)(g);
        //         *mi = MediumInteraction(rWorld(t), -rWorld.d, rWorld.time, this,
        //                                 phase);
        //         return sigma_s / sigma_t;
        //     }
        // }
        (Spectrum::new(1.0 as Float), None)
    }
}
