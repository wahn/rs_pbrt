// pbrt
use core::pbrt::Float;

pub struct MLTSampler {}

pub struct MLTIntegrator {
    pub max_depth: u32,
    pub n_bootstrap: u32,
    pub n_chains: u32,
    pub mutations_per_pixel: u32,
    pub sigma: Float,
    pub large_step_probability: Float,
}

impl MLTIntegrator {
    pub fn new(
        max_depth: u32,
        n_bootstrap: u32,
        n_chains: u32,
        mutations_per_pixel: u32,
        sigma: Float,
        large_step_probability: Float,
    ) -> Self {
        MLTIntegrator {
            max_depth: max_depth,
            n_bootstrap: n_bootstrap,
            n_chains: n_chains,
            mutations_per_pixel: mutations_per_pixel,
            sigma: sigma,
            large_step_probability: large_step_probability,
        }
    }
}
