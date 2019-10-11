//! **Integrator** is an abstract base class that defines the
//! **render()** method that must be provided by all integrators.
//!
//! - AOIntegrator
//! - BDPTIntegrator
//! - DirectLightingIntegrator
//! - MLTIntegrator
//! - PathIntegrator
//! - SPPMIntegrator
//! - VolPathIntegrator
//! - WhittedIntegrator
//!
//! ## Ambient Occlusion (AO)
//!
//! Ambient Occlusion is most often calculated by casting rays in
//! every direction from the surface. Rays which reach the background
//! or sky increase the brightness of the surface, whereas a ray which
//! hits any other object contributes no illumination. As a result,
//! points surrounded by a large amount of geometry are rendered dark,
//! whereas points with little geometry on the visible hemisphere
//! appear light.
//!
//! ![Ambient Occlusion](/doc/img/cornell_box_pbrt_rust_ao.png)
//!
//! ## Direct Lighting
//!
//! The **DirectLightingIntegrator** accounts only for direct lighting
//! &mdash; light that has traveled directly from a light source to the
//! point being shaded &mdash; and ignores indirect illumination from
//! objects that are not themselfes emissive, except for basic
//! specular reflection and transmission effects.
//!
//! ![Direct Lighting](/doc/img/cornell_box_pbrt_rust_directlighting.png)
//!
//! ## Path Tracing
//!
//! Path tracing incrementally generates paths of scattering events
//! starting at the camera and ending at light sources in the scene.
//!
//! ![Path Tracing](/doc/img/cornell_box_pbrt_rust_path.png)
//!
//! ## Bidirectional Path Tracing (BDPT)
//!
//! Bidirectional path tracing is a generalization of the standard
//! pathtracing algorithm that can be much more efficient. It
//! constructs paths that start from the camera on one end, from the
//! light on the other end, and connects in the middle with a
//! visibility ray.
//!
//! ![Bidirectional Path
//! Tracing](/doc/img/art_gallery_pbrt_rust_bdpt.png)
//!
//! ## Stochastic Progressive Photon Mapping (SPPM)
//!
//! A photon mapping integrator that uses particles to estimate
//! illumination by interpolating lighting contributions from
//! particles close to but not quite at the point being shaded.
//!
//! ![Stochastic Progressive Photon Mapping](/doc/img/caustic_glass_pbrt_rust_sppm.png)

// std
use std::sync::Arc;
// pbrt
use crate::blockqueue::BlockQueue;
use crate::core::camera::{Camera, CameraSample};
use crate::core::geometry::pnt2_inside_exclusive;
use crate::core::geometry::{Bounds2i, Point2i, Ray, Vector2i};
use crate::core::integrator::SamplerIntegrator;
use crate::core::pbrt::{Float, Spectrum};
use crate::core::sampler::Sampler;
use crate::core::scene::Scene;

pub mod ao;
pub mod bdpt;
pub mod directlighting;
pub mod mlt;
pub mod path;
pub mod sppm;
pub mod volpath;
pub mod whitted;

/// **Main function** to **render** a scene mutli-threaded (using all
/// available cores).
pub fn render(
    scene: &Scene,
    camera: &Arc<dyn Camera + Send + Sync>,
    sampler: &mut Box<dyn Sampler + Send + Sync>,
    integrator: &mut Box<dyn SamplerIntegrator + Send + Sync>,
    num_threads: u8,
) {
    // SamplerIntegrator::Render (integrator.cpp)
    let film = camera.get_film();
    let sample_bounds: Bounds2i = film.get_sample_bounds();
    integrator.preprocess(scene, sampler);
    // use camera below
    let sample_extent: Vector2i = sample_bounds.diagonal();
    let tile_size: i32 = 16;
    let x: i32 = (sample_extent.x + tile_size - 1) / tile_size;
    let y: i32 = (sample_extent.y + tile_size - 1) / tile_size;
    let n_tiles: Point2i = Point2i { x, y };
    // TODO: ProgressReporter reporter(nTiles.x * nTiles.y, "Rendering");
    let num_cores: usize;
    if num_threads == 0_u8 {
        num_cores = num_cpus::get();
    } else {
        num_cores = num_threads as usize;
    }
    println!("Rendering with {:?} thread(s) ...", num_cores);
    {
        let block_queue = BlockQueue::new(
            (
                (n_tiles.x * tile_size) as u32,
                (n_tiles.y * tile_size) as u32,
            ),
            (tile_size as u32, tile_size as u32),
            (0, 0),
        );
        let integrator = &integrator;
        let bq = &block_queue;
        let sampler = sampler;
        let camera = &camera;
        let film = &film;
        let pixel_bounds = integrator.get_pixel_bounds().clone();
        crossbeam::scope(|scope| {
            let (pixel_tx, pixel_rx) = crossbeam_channel::bounded(num_cores);
            // spawn worker threads
            for _ in 0..num_cores {
                let pixel_tx = pixel_tx.clone();
                let mut tile_sampler: Box<dyn Sampler + Send + Sync> = sampler.box_clone();
                scope.spawn(move |_| {
                    while let Some((x, y)) = bq.next() {
                        let tile: Point2i = Point2i {
                            x: x as i32,
                            y: y as i32,
                        };
                        let seed: i32 = tile.y * n_tiles.x + tile.x;
                        tile_sampler.reseed(seed as u64);
                        let x0: i32 = sample_bounds.p_min.x + tile.x * tile_size;
                        let x1: i32 = std::cmp::min(x0 + tile_size, sample_bounds.p_max.x);
                        let y0: i32 = sample_bounds.p_min.y + tile.y * tile_size;
                        let y1: i32 = std::cmp::min(y0 + tile_size, sample_bounds.p_max.y);
                        let tile_bounds: Bounds2i =
                            Bounds2i::new(Point2i { x: x0, y: y0 }, Point2i { x: x1, y: y1 });
                        // println!("Starting image tile {:?}", tile_bounds);
                        let mut film_tile = film.get_film_tile(&tile_bounds);
                        for pixel in &tile_bounds {
                            tile_sampler.start_pixel(&pixel);
                            if !pnt2_inside_exclusive(&pixel, &pixel_bounds) {
                                continue;
                            }
                            let mut done: bool = false;
                            while !done {
                                // let's use the copy_arena crate instead of pbrt's MemoryArena
                                // let mut arena: Arena = Arena::with_capacity(262144); // 256kB

                                // initialize _CameraSample_ for current sample
                                let camera_sample: CameraSample =
                                    tile_sampler.get_camera_sample(&pixel);
                                // generate camera ray for current sample
                                let mut ray: Ray = Ray::default();
                                let ray_weight: Float =
                                    camera.generate_ray_differential(&camera_sample, &mut ray);
                                ray.scale_differentials(
                                    1.0 as Float
                                        / (tile_sampler.get_samples_per_pixel() as Float).sqrt(),
                                );
                                // TODO: ++nCameraRays;
                                // evaluate radiance along camera ray
                                let mut l: Spectrum = Spectrum::new(0.0 as Float);
                                let y: Float = l.y();
                                if ray_weight > 0.0 {
                                    l = integrator.li(
                                        &mut ray,
                                        scene,
                                        &mut tile_sampler, // &mut arena,
                                        0_i32,
                                    );
                                }
                                if l.has_nans() {
                                    println!(
                                        "Not-a-number radiance value returned for pixel \
                                         ({:?}, {:?}), sample {:?}. Setting to black.",
                                        pixel.x,
                                        pixel.y,
                                        tile_sampler.get_current_sample_number()
                                    );
                                    l = Spectrum::new(0.0);
                                } else if y < -10.0e-5 as Float {
                                    println!(
                                        "Negative luminance value, {:?}, returned for pixel \
                                         ({:?}, {:?}), sample {:?}. Setting to black.",
                                        y,
                                        pixel.x,
                                        pixel.y,
                                        tile_sampler.get_current_sample_number()
                                    );
                                    l = Spectrum::new(0.0);
                                } else if y.is_infinite() {
                                    println!(
                                        "Infinite luminance value returned for pixel ({:?}, \
                                         {:?}), sample {:?}. Setting to black.",
                                        pixel.x,
                                        pixel.y,
                                        tile_sampler.get_current_sample_number()
                                    );
                                    l = Spectrum::new(0.0);
                                }
                                // println!("Camera sample: {:?} -> ray: {:?} -> L = {:?}",
                                //          camera_sample, ray, l);
                                // add camera ray's contribution to image
                                film_tile.add_sample(&camera_sample.p_film, &mut l, ray_weight);
                                done = !tile_sampler.start_next_sample();
                            } // arena is dropped here !
                        }
                        // send the tile through the channel to main thread
                        pixel_tx
                            .send(film_tile)
                            .expect(&format!("Failed to send tile"));
                    }
                });
            }
            // spawn thread to collect pixels and render image to file
            scope.spawn(move |_| {
                for _ in pbr::PbIter::new(0..bq.len()) {
                    let film_tile = pixel_rx.recv().unwrap();
                    // merge image tile into _Film_
                    film.merge_film_tile(&film_tile);
                }
            });
        })
        .unwrap();
    }
    film.write_image(1.0 as Float);
}
