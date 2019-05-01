// std
use std::sync::Arc;
// pbrt
use core::geometry::{Bounds2i, Point2f, Ray, Vector3f};
use core::lightdistrib::LightDistribution;
use core::pbrt::{Float, Spectrum};

// see volpath.h

/// Accounts for scattering and attenuation from participating media
/// as well as scattering from surfaces
pub struct VolPathIntegrator {
    // inherited from SamplerIntegrator (see integrator.h)
    pixel_bounds: Bounds2i,
    // see volpath.h
    pub max_depth: u32,
    rr_threshold: Float,           // 1.0
    light_sample_strategy: String, // "spatial"
    light_distribution: Option<Arc<LightDistribution + Send + Sync>>,
}
