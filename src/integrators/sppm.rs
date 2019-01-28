extern crate crossbeam;
extern crate num_cpus;

// std
use std;
use std::sync::Arc;
// pbrt
use atomic::{Atomic, Ordering};
use core::camera::Camera;
use core::film::Film;
use core::geometry::{Bounds2f, Bounds2i, Point2f, Point2i, Vector2i};
use core::integrator::compute_light_power_distribution;
use core::parallel::AtomicFloat;
use core::pbrt::{Float, Spectrum};
use core::sampler::Sampler;
use core::sampling::Distribution1D;
use core::scene::Scene;
use samplers::halton::HaltonSampler;

pub struct SPPMIntegrator {
    pub initial_search_radius: Float,
    pub n_iterations: i32,
}

impl SPPMIntegrator {
    pub fn new(
        camera: Arc<Camera + Send + Sync>,
        n_iterations: i32,
        photons_per_iteration: i32,
        max_depth: i32,
        initial_search_radius: Float,
        write_frequency: i32,
    ) -> Self {
        SPPMIntegrator {
            initial_search_radius: initial_search_radius,
            n_iterations: n_iterations,
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
pub fn render_sppm(
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
    let inv_sqrt_spp: Float = 1.0 as Float / (integrator.n_iterations as Float).sqrt();
    // TODO: let pixel_memory_bytes: usize = n_pixels as usize * std::mem::size_of::<SPPMPixel>();

    // compute _lightDistr_ for sampling lights proportional to power
    let light_distr_opt: Option<Arc<Distribution1D>> = compute_light_power_distribution(scene);
    // perform _nIterations_ of SPPM integration
    let sampler: HaltonSampler =
        HaltonSampler::new(integrator.n_iterations as i64, pixel_bounds, false);
    // compute number of tiles to use for SPPM camera pass
    let pixel_extent: Vector2i = pixel_bounds.diagonal();
    let tile_size: i32 = 16;
    let n_tiles: Point2i = Point2i {
        x: (pixel_extent.x + tile_size - 1) / tile_size,
        y: (pixel_extent.y + tile_size - 1) / tile_size,
    };
    println!("n_tiles = {:?}", n_tiles);
    // int x0 = pixelBounds.pMin.x + tile.x * tileSize;
    // int x1 = std::min(x0 + tileSize, pixelBounds.pMax.x);
    // int y0 = pixelBounds.pMin.y + tile.y * tileSize;
    // int y1 = std::min(y0 + tileSize, pixelBounds.pMax.y);
    // WORK
}
