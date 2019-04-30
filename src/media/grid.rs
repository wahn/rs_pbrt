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
    pub fn density(&self, p: &Point3f) -> Float {
        // TODO
        0.0 as Float
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
        let mut t_min: Float = 0.0;
        let mut t_max: Float = 0.0;
        if !b.intersect_b(&ray, &mut Some(t_min), &mut Some(t_max)) {
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
                // populate _mi_ with medium interaction information and return
                //     PhaseFunction *phase = ARENA_ALLOC(arena, HenyeyGreenstein)(g);
                //     *mi = MediumInteraction(rWorld(t), -rWorld.d, rWorld.time, this,
                //                             phase);
                //     return sigma_s / sigma_t;
            }
        }
        (Spectrum::new(1.0 as Float), None)
    }
}
