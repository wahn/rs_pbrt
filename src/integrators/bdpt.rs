// pbrt
use core::geometry::{Bounds2i, Ray, Vector3f};

// see bdpt.h

/// Bidirectional Path Tracing (Global Illumination)
pub struct BDPTIntegrator {
    max_depth: u32,
    visualize_strategies: bool,
    visualize_weights: bool,
    pixel_bounds: Bounds2i,
    light_sample_strategy: String, // "power"
}

impl BDPTIntegrator {
    pub fn new(
        // TODO: sampler
        // TODO: camera
        max_depth: u32,
        visualize_strategies: bool,
        visualize_weights: bool,
        pixel_bounds: Bounds2i,
        light_sample_strategy: String,
    ) -> Self {
        BDPTIntegrator {
            max_depth: max_depth,
            visualize_strategies: visualize_strategies,
            visualize_weights: visualize_weights,
            pixel_bounds: pixel_bounds,
            light_sample_strategy: light_sample_strategy,
        }
    }
}
