// pbrt
use core::medium::{HenyeyGreenstein, Medium};
use core::pbrt::{Float, Spectrum};
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
