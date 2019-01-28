extern crate crossbeam;
extern crate num_cpus;

// std
use std::sync::Arc;
// pbrt
use atomic::{Atomic, Ordering};
use core::camera::Camera;
use core::film::Film;
use core::geometry::{Bounds2f, Bounds2i, Point2f, Point2i};
use core::parallel::AtomicFloat;
use core::pbrt::{Float, Spectrum};
use core::sampler::Sampler;
use core::scene::Scene;

pub struct SPPMIntegrator {
    pub initial_search_radius: Float,
}

impl SPPMIntegrator {
    pub fn new(initial_search_radius: Float) -> Self {
        SPPMIntegrator {
            initial_search_radius: initial_search_radius,
        }
    }
}

#[derive(Debug, Default)]
pub struct VisiblePoint {}

#[derive(Debug, Default)]
pub struct SPPMPixel {
    pub radius: Float,
    pub ld: Spectrum,
    pub vp: VisiblePoint,
    pub phi: [AtomicFloat; 3],
    pub m: Atomic<i32>,
    pub n: Float,
    pub tau: Spectrum,
}

/// **Main function** to **render** a scene multi-threaded (using all
/// available cores) with **Stochastic Progressive Photon Mapping**
/// (SPPM).
pub fn render_mlt(
    scene: &Scene,
    camera: &Arc<Camera + Send + Sync>,
    sampler: &mut Box<Sampler + Send + Sync>,
    integrator: &mut Box<SPPMIntegrator>,
    num_threads: u8,
) {
    // TODO: ProfilePhase p(Prof::IntegratorRender);

    // initialize _pixel_bounds_ and _pixels_ array for SPPM
    let film = camera.get_film();
    let pixel_bounds: Bounds2i = film.cropped_pixel_bounds;
    let n_pixels: i32 = pixel_bounds.area();
    let mut pixels: Vec<SPPMPixel> = Vec::with_capacity(n_pixels as usize);
    for i in 0..n_pixels as usize {
        pixels.push(SPPMPixel::default());
        pixels[i].radius = integrator.initial_search_radius;
    }
    // WORK
}
