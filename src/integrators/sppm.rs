extern crate crossbeam;
extern crate num_cpus;

// std
use std;
use std::sync::Arc;
// pbrt
use atomic::{Atomic, Ordering};
use core::camera::{Camera, CameraSample};
use core::film::Film;
use core::geometry::{Bounds2f, Bounds2i, Point2f, Point2i, Point3f, Ray, Vector2i, Vector3f};
use core::integrator::{compute_light_power_distribution, uniform_sample_one_light};
use core::interaction::{Interaction, SurfaceInteraction};
use core::material::{Material, TransportMode};
use core::parallel::AtomicFloat;
use core::pbrt::{Float, Spectrum};
use core::reflection::{Bsdf, BxdfType};
use core::sampler::{GlobalSampler, Sampler, SamplerClone};
use core::sampling::Distribution1D;
use core::scene::Scene;
use samplers::halton::HaltonSampler;

pub struct SPPMIntegrator {
    pub initial_search_radius: Float,
    pub n_iterations: i32,
    pub max_depth: u32,
}

impl SPPMIntegrator {
    pub fn new(
        camera: Arc<Camera + Send + Sync>,
        n_iterations: i32,
        photons_per_iteration: i32,
        max_depth: u32,
        initial_search_radius: Float,
        write_frequency: i32,
    ) -> Self {
        SPPMIntegrator {
            initial_search_radius: initial_search_radius,
            n_iterations: n_iterations,
            max_depth: max_depth,
        }
    }
}

#[derive(Default)]
pub struct VisiblePoint {
    pub p: Point3f,
    pub wo: Vector3f,
    pub bsdf: Option<Arc<Bsdf>>,
    pub beta: Spectrum,
}

#[derive(Default)]
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
    let num_cores: usize;
    if num_threads == 0_u8 {
        num_cores = num_cpus::get();
    } else {
        num_cores = num_threads as usize;
    }
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

    // compute _light_distr_ for sampling lights proportional to power
    let light_distr_opt: Option<Arc<Distribution1D>> = compute_light_power_distribution(scene);
    // perform _n_iterations_ of SPPM integration
    let sampler: Box<HaltonSampler> = Box::new(HaltonSampler::new(
        integrator.n_iterations as i64,
        pixel_bounds,
        false,
    ));
    // compute number of tiles to use for SPPM camera pass
    let pixel_extent: Vector2i = pixel_bounds.diagonal();
    let tile_size: i32 = 16;
    let n_tiles: Point2i = Point2i {
        x: (pixel_extent.x + tile_size - 1) / tile_size,
        y: (pixel_extent.y + tile_size - 1) / tile_size,
    };
    // TODO: ProgressReporter progress(2 * nIterations, "Rendering");
    for iter in 0..integrator.n_iterations {
        // generate SPPM visible points
        {
            // TODO: ProfilePhase _(Prof::SPPMCameraPass);
            // TODO: ParallelFor2D([&](Point2i tile) {
            for y in 0..n_tiles.y {
                for x in 0..n_tiles.x {
                    let tile: Point2i = Point2i { x: x, y: y };
                    // TODO: MemoryArena &arena = perThreadArenas[ThreadIndex];

                    // follow camera paths for _tile_ in image for SPPM
                    let tile_index: i32 = tile.y * n_tiles.x + tile.x;
                    let mut tile_sampler = sampler.clone();
                    // compute _tileBounds_ for SPPM tile
                    let x0: i32 = pixel_bounds.p_min.x + tile.x * tile_size;
                    let x1: i32 = std::cmp::min(x0 + tile_size, pixel_bounds.p_max.x);
                    let y0: i32 = pixel_bounds.p_min.y + tile.y * tile_size;
                    let y1: i32 = std::cmp::min(y0 + tile_size, pixel_bounds.p_max.y);
                    let tile_bounds: Bounds2i =
                        Bounds2i::new(Point2i { x: x0, y: y0 }, Point2i { x: x1, y: y1 });
                    for p_pixel in &tile_bounds {
                        // prepare _tileSampler_ for _p_pixel_
                        tile_sampler.start_pixel(&p_pixel);
                        tile_sampler.set_sample_number(iter as i64);
                        // generate camera ray for pixel for SPPM
                        let camera_sample: CameraSample = tile_sampler.get_camera_sample(&p_pixel);
                        let mut ray: Ray = Ray::default();
                        let beta: Spectrum = Spectrum::new(
                            camera.generate_ray_differential(&camera_sample, &mut ray),
                        );
                        if beta.is_black() {
                            continue;
                        }
                        ray.scale_differentials(inv_sqrt_spp);

                        // follow camera ray path until a visible point is created

                        // get _SPPMPixel_ for _p_pixel_
                        let p_pixel_o: Point2i = Point2i::from(p_pixel - pixel_bounds.p_min);
                        let pixel_offset: i32 = p_pixel_o.x
                            + p_pixel_o.y * (pixel_bounds.p_max.x - pixel_bounds.p_min.x);
                        let pixel: &mut SPPMPixel = &mut pixels[pixel_offset as usize];
                        let specular_bounce: bool = false;
                        for depth in 0..integrator.max_depth {
                            // TODO: ++totalPhotonSurfaceInteractions;
                            if let Some(mut isect) = scene.intersect(&mut ray) {
                                // process SPPM camera ray intersection

                                // compute BSDF at SPPM camera ray intersection
                                let mode: TransportMode = TransportMode::Radiance;
                                isect.compute_scattering_functions(
                                    &mut ray, // arena,
                                    true, mode,
                                );
                                if let Some(ref bsdf) = isect.bsdf {
                                    // accumulate direct illumination
                                    // at SPPM camera ray intersection
                                    let wo: Vector3f = -ray.d;
                                    if depth == 0 || specular_bounce {
                                        pixel.ld += beta * isect.le(&wo);
                                    }
                                    pixel.ld += beta
                                        * uniform_sample_one_light(
                                            &isect,
                                            scene,
                                            &mut tile_sampler.box_clone(),
                                            false,
                                            None,
                                        );
                                    // possibly create visible point and end camera path
                                    let mut bsdf_flags: u8 = BxdfType::BsdfDiffuse as u8
                                        | BxdfType::BsdfReflection as u8
                                        | BxdfType::BsdfTransmission as u8;
                                    let is_diffuse: bool = bsdf.num_components(bsdf_flags) > 0;
                                    bsdf_flags = BxdfType::BsdfGlossy as u8
                                        | BxdfType::BsdfReflection as u8
                                        | BxdfType::BsdfTransmission as u8;
                                    let is_glossy: bool = bsdf.num_components(bsdf_flags) > 0;
                                    if is_diffuse
                                        || (is_glossy && depth == integrator.max_depth - 1)
                                    {
                                        pixel.vp.p = isect.p;
                                        pixel.vp.wo = wo;
                                        pixel.vp.bsdf = Some(bsdf.clone());
                                        pixel.vp.beta = beta;
                                        break;
                                    }
                                // WORK
                                } else {
                                    ray = isect.spawn_ray(&ray.d);
                                    // --depth;
                                    continue;
                                }
                            } else {
                                // accumulate light contributions for
                                // ray with no intersection
                                for light in &scene.lights {
                                    pixel.ld += beta * light.le(&mut ray);
                                }
                                break;
                            }
                        }
                    }
                }
            }
        }
        // create grid of all SPPM visible points
        // WORK
    }
}
