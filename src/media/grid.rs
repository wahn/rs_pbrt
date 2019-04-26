// pbrt
use core::geometry::Ray;
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
    fn sample(&self, ray: &Ray, sampler: &mut Sampler) -> (Spectrum, Option<MediumInteraction>) {
        // TODO
        (Spectrum::default(), None)
    }
}
