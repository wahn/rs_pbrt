//! # pbrt
//!
//! [Rust][rust] crate to implement at least parts of the [PBRT
//! book][book]'s C++ code. You can find a copy of the current code
//! [here][repo].
//!
//! [rust]: https://www.rust-lang.org/en-US
//! [book]: http://www.pbrt.org
//! [repo]: https://github.com/wahn/rs_pbrt
//!

extern crate atomic;
extern crate crossbeam;
#[cfg(feature = "openexr")]
extern crate half;
extern crate image;
extern crate num;
extern crate num_cpus;
#[cfg(feature = "openexr")]
extern crate openexr;
extern crate pbr;
extern crate ply_rs;
extern crate time;
extern crate typed_arena;

// use std::cell::RefCell;
// use std::collections::HashMap;
use std::default::Default;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::sync::mpsc;

pub mod accelerators;
pub mod cameras;
pub mod core;
pub mod filters;
pub mod integrators;
pub mod lights;
pub mod materials;
pub mod samplers;
pub mod shapes;
pub mod textures;

// pbrt
use core::camera::{Camera, CameraSample};
use core::geometry::{Bounds2i, Point2f, Point2i, Ray, Vector2i};
use core::geometry::pnt2_inside_exclusive;
use core::integrator::SamplerIntegrator;
// use core::light::Light;
use core::lightdistrib::create_light_sample_distribution;
use core::pbrt::{Float, Spectrum};
use core::sampler::Sampler;
use core::sampling::Distribution1D;
use core::scene::Scene;
use integrators::bdpt::{BDPTIntegrator, Vertex};
use integrators::bdpt::{connect_bdpt, generate_camera_subpath, generate_light_subpath};

// see github/tray_rust/src/sampler/block_queue.rs

/// The queue of blocks to be worked on shared immutably between worker threads.
pub struct BlockQueue {
    /// The block indices of blocks to work on for the image
    blocks: Vec<(u32, u32)>,
    /// Get the dimensions of an individual block
    dimensions: (u32, u32),
    /// Index of the next block to be worked on
    next: AtomicUsize,
}

impl BlockQueue {
    /// Create a block queue for the image with dimensions `img`.
    /// Panics if the image is not evenly broken into blocks of dimension `dim`
    pub fn new(img: (u32, u32), dim: (u32, u32), select_blocks: (usize, usize)) -> BlockQueue {
        if img.0 % dim.0 != 0 || img.1 % dim.1 != 0 {
            panic!(
                "Image with dimension {:?} not evenly divided by blocks of {:?}",
                img, dim
            );
        }
        let num_blocks = (img.0 / dim.0, img.1 / dim.1);
        // TODO: the .. operator precedence is very low so we need this paren here at the moment
        // once (hopefully) it's raised we can remove the parens
        let mut blocks: Vec<(u32, u32)> = (0..num_blocks.0 * num_blocks.1)
            .map(|i| (i % num_blocks.0, i / num_blocks.0))
            .collect();
        blocks.sort_by(|a, b| morton2(a).cmp(&morton2(b)));
        // If we're only rendering a subset of the blocks then filter our list down
        if select_blocks.1 > 0 {
            blocks = blocks
                .into_iter()
                .skip(select_blocks.0)
                .take(select_blocks.1)
                .collect();
        }
        if blocks.is_empty() {
            println!("Warning: This block queue is empty!");
        }
        BlockQueue {
            blocks: blocks,
            dimensions: dim,
            next: AtomicUsize::new(0),
        }
    }
    /// Get the dimensions of an individual block in the queue
    pub fn block_dim(&self) -> (u32, u32) {
        self.dimensions
    }
    /// Get an iterator to work through the queue
    pub fn iter(&self) -> BlockQueueIterator {
        BlockQueueIterator { queue: self }
    }
    /// Get the next block in the queue or None if the queue is finished
    fn next(&self) -> Option<(u32, u32)> {
        let i = self.next.fetch_add(1, Ordering::AcqRel);
        if i >= self.blocks.len() {
            None
        } else {
            Some(self.blocks[i])
        }
    }
    /// Get the length of the queue
    pub fn len(&self) -> usize {
        self.blocks.len()
    }
    /// Check if the queue is empty
    pub fn is_empty(&self) -> bool {
        self.next.load(Ordering::AcqRel) >= self.blocks.len()
    }
}

/// Iterator to work through the queue safely
pub struct BlockQueueIterator<'a> {
    queue: &'a BlockQueue,
}

impl<'a> Iterator for BlockQueueIterator<'a> {
    type Item = (u32, u32);
    fn next(&mut self) -> Option<(u32, u32)> {
        self.queue.next()
    }
}
// see github/tray_rust/src/sampler/morton.rs

///! Provides utilities for 2D Morton code generation using Fabian
///! Giesen's Morton code decoding functions, see [his post on Morton
///! codes](https://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/)

/// Insert a 0 bit between each of the low 16 bits of x
fn part1_by1(mut x: u32) -> u32 {
    // x = ---- ---- ---- ---- fedc ba98 7654 3210
    x &= 0x0000ffff;
    // x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x ^ (x << 8)) & 0x00ff00ff;
    // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x ^ (x << 4)) & 0x0f0f0f0f;
    // x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x ^ (x << 2)) & 0x33333333;
    // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    (x ^ (x << 1)) & 0x55555555
}
/// Compute the Morton code for the `(x, y)` position.
fn morton2(p: &(u32, u32)) -> u32 {
    (part1_by1(p.1) << 1) + part1_by1(p.0)
}

/// **Main function** to **render** a scene mutli-threaded (using all
/// available cores).
pub fn render(
    scene: &Scene,
    camera: &Box<Camera + Send + Sync>,
    sampler: &mut Box<Sampler + Send + Sync>,
    integrator: &mut Box<SamplerIntegrator + Send + Sync>,
    num_threads: u8,
) {
    // SamplerIntegrator::Render (integrator.cpp)
    let film = camera.get_film();
    let sample_bounds: Bounds2i = film.get_sample_bounds();
    println!("sample_bounds = {:?}", sample_bounds);
    integrator.preprocess(scene, sampler);
    // use camera below
    let sample_extent: Vector2i = sample_bounds.diagonal();
    println!("sample_extent = {:?}", sample_extent);
    let tile_size: i32 = 16;
    let x: i32 = (sample_extent.x + tile_size - 1) / tile_size;
    let y: i32 = (sample_extent.y + tile_size - 1) / tile_size;
    let n_tiles: Point2i = Point2i { x: x, y: y };
    println!("n_tiles = {:?}", n_tiles);
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
        println!("block_queue.len() = {}", block_queue.len());
        let integrator = &integrator;
        let bq = &block_queue;
        let sampler = sampler;
        let camera = &camera;
        let film = &film;
        let pixel_bounds = integrator.get_pixel_bounds().clone();
        crossbeam::scope(|scope| {
            let (pixel_tx, pixel_rx) = mpsc::channel();
            // spawn worker threads
            for _ in 0..num_cores {
                let pixel_tx = pixel_tx.clone();
                let mut tile_sampler: Box<Sampler + Send + Sync> = sampler.box_clone();
                scope.spawn(move || {
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
            scope.spawn(move || {
                for _ in pbr::PbIter::new(0..bq.len()) {
                    let film_tile = pixel_rx.recv().unwrap();
                    // merge image tile into _Film_
                    film.merge_film_tile(&film_tile);
                }
            });
        });
    }
    println!("Rendering finished");
    film.write_image(1.0 as Float);
}

/// **Main function** to **render** a scene mutli-threaded (using all
/// available cores) with **bidirectional** path tracing.
pub fn render_bdpt(
    scene: &Scene,
    camera: &Box<Camera + Send + Sync>,
    sampler: &mut Box<Sampler + Send + Sync>,
    integrator: &mut Box<BDPTIntegrator>,
    num_threads: u8,
) {
    // TODO
    // Compute a reverse mapping from light pointers to offsets into
    // the scene lights vector (and, equivalently, offsets into
    // lightDistr). Added after book text was finalized; this is
    // critical to reasonable performance with 100s+ of light sources.
    // let mut light_to_index = HashMap::new();
    // for li in 0..scene.lights.len() {
    //     let ref light = scene.lights[li];
    //     light_to_index.insert(light, li);
    // }
    // partition the image into tiles
    let film = camera.get_film();
    let sample_bounds: Bounds2i = film.get_sample_bounds();
    println!("sample_bounds = {:?}", sample_bounds);
    let sample_extent: Vector2i = sample_bounds.diagonal();
    println!("sample_extent = {:?}", sample_extent);
    let tile_size: i32 = 16;
    let n_x_tiles: i32 = (sample_extent.x + tile_size - 1) / tile_size;
    let n_y_tiles: i32 = (sample_extent.y + tile_size - 1) / tile_size;
    // TODO: ProgressReporter reporter(nXTiles * nYTiles, "Rendering");
    // TODO: Allocate buffers for debug visualization
    // ...
    // render and write the output image to disk
    if scene.lights.len() > 0 {
        let samples_per_pixel: i64 = sampler.get_samples_per_pixel();
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
                    (n_x_tiles * tile_size) as u32,
                    (n_y_tiles * tile_size) as u32,
                ),
                (tile_size as u32, tile_size as u32),
                (0, 0),
            );
            println!("block_queue.len() = {}", block_queue.len());
            let integrator = &integrator;
            let bq = &block_queue;
            let sampler = sampler;
            let camera = &camera;
            let film = &film;
            // let pixel_bounds = integrator.get_pixel_bounds().clone();
            crossbeam::scope(|scope| {
                let (pixel_tx, pixel_rx) = mpsc::channel();
                // spawn worker threads
                for _ in 0..num_cores {
                    let pixel_tx = pixel_tx.clone();
                    let mut tile_sampler: Box<Sampler + Send + Sync> = sampler.box_clone();
                    scope.spawn(move || {
                        while let Some((x, y)) = bq.next() {
                            let tile: Point2i = Point2i {
                                x: x as i32,
                                y: y as i32,
                            };
                            let seed: i32 = tile.y * n_x_tiles + tile.x;
                            tile_sampler.reseed(seed as u64);
                            let x0: i32 = sample_bounds.p_min.x + tile.x * tile_size;
                            let x1: i32 = std::cmp::min(x0 + tile_size, sample_bounds.p_max.x);
                            let y0: i32 = sample_bounds.p_min.y + tile.y * tile_size;
                            let y1: i32 = std::cmp::min(y0 + tile_size, sample_bounds.p_max.y);
                            let tile_bounds: Bounds2i =
                                Bounds2i::new(Point2i { x: x0, y: y0 }, Point2i { x: x1, y: y1 });
                            // println!("Starting image tile {:?}", tile_bounds);
                            let mut film_tile = film.get_film_tile(&tile_bounds);
                            for p_pixel in &tile_bounds {
                                tile_sampler.start_pixel(&p_pixel);
                                if !pnt2_inside_exclusive(&p_pixel, &integrator.pixel_bounds) {
                                    continue;
                                }
                                let mut done: bool = false;
                                while !done {
                                    // Get a distribution for sampling
                                    // the light at the start of the
                                    // light subpath. Because the
                                    // light path follows multiple
                                    // bounces, basing the sampling
                                    // distribution on any of the
                                    // vertices of the camera path is
                                    // unlikely to be a good
                                    // strategy. We use the
                                    // PowerLightDistribution by
                                    // default here, which doesn't use
                                    // the point passed to it. Now
                                    // trace the light subpath
                                    if let Some(light_distribution) =
                                        create_light_sample_distribution(
                                            integrator.get_light_sample_strategy(),
                                            scene,
                                        ) {
                                        // generate a single sample using BDPT
                                        let p_film: Point2f = Point2f {
                                            x: p_pixel.x as Float,
                                            y: p_pixel.y as Float,
                                        }
                                            + tile_sampler.get_2d();
                                        // trace the camera subpath
                                        let mut camera_vertices: Vec<Vertex> =
                                        Vec::with_capacity((integrator.max_depth + 2) as usize);
                                        let mut n_camera;
                                        let mut p;
                                        let mut time;
                                        {
                                            let (n_camera_new, p_new, time_new) =
                                                generate_camera_subpath(
                                                    scene,
                                                    &mut tile_sampler,
                                                    integrator.max_depth + 2,
                                                    camera,
                                                    &p_film,
                                                    &mut camera_vertices,
                                                );
                                            n_camera = n_camera_new;
                                            p = p_new;
                                            time = time_new;
                                        }
                                        let light_distr: Arc<Distribution1D> =
                                                light_distribution.lookup(&p);
                                        let mut light_vertices: Vec<Vertex> =
                                                Vec::with_capacity((integrator.max_depth + 1)
                                                                   as usize);
                                        let mut n_light;
                                        {
                                            n_light = generate_light_subpath(
                                                scene,
                                                &mut tile_sampler,
                                                integrator.max_depth + 1,
                                                time,
                                                &light_distr,
                                                // light_to_index,
                                                &mut light_vertices,
                                            );
                                        }
                                        // Execute all BDPT connection strategies
                                        let mut l: Spectrum = Spectrum::new(0.0 as Float);
                                        // println!("n_camera = {:?}", n_camera);
                                        // println!("n_light = {:?}", n_light);
                                        for t in 1..n_camera + 1 {
                                            for s in 0..n_light + 1 {
                                                // int depth = t + s - 2;
                                                let depth: isize = (t + s) as isize - 2;
                                                if (s == 1 && t == 1) || depth < 0
                                                    || depth > integrator.max_depth as isize
                                                {
                                                    continue;
                                                }
                                                // execute the $(s, t)$ connection strategy and update _L_
                                                let mut p_film_new: Point2f = Point2f {
                                                    x: p_film.x,
                                                    y: p_film.y,
                                                };
                                                let mut mis_weight: Option<Float> =
                                                    Some(0.0 as Float);
                                                let lpath: Spectrum = connect_bdpt(
                                                    scene,
                                                    &light_vertices,
                                                    &camera_vertices,
                                                    s,
                                                    t,
                                                    &light_distr,
                                                    camera,
                                                    &mut tile_sampler,
                                                    &mut p_film_new,
                                                    mis_weight.as_mut(),
                                                );
                                                // if let Some(mis_weight_flt) = mis_weight {
                                                //     println!("Connect bdpt s: {:?}, t: {:?}, lpath: {:?}, mis_weight: {:?}",
                                                //              s, t, lpath, mis_weight_flt);
                                                // }
                                                // if (visualizeStrategies || visualizeWeights) {
                                                //     Spectrum value;
                                                //     if (visualizeStrategies)
                                                //         value =
                                                //             mis_weight == 0 ? 0 : lpath / mis_weight;
                                                //     if (visualizeWeights) value = lpath;
                                                //     weightFilms[BufferIndex(s, t)]->AddSplat(
                                                //         pFilmNew, value);
                                                // }
                                                if t != 1 {
                                                    l += lpath;
                                                }
                                                else {
                                                    film.add_splat(&p_film_new, &lpath);
                                                }
                                            }
                                        }
                                        // println!(
                                        //     "Add film sample pFilm: {:?}, L: {:?}, (y: {:?})",
                                        //     p_film,
                                        //     l,
                                        //     l.y()
                                        // );
                                        film_tile.add_sample(&p_film, &mut l, 1.0 as Float);
                                        done = !tile_sampler.start_next_sample();
                                    }
                                }
                            }
                            // send the tile through the channel to main thread
                            pixel_tx
                                .send(film_tile)
                                .expect(&format!("Failed to send tile"));
                        }
                    });
                }
                // spawn thread to collect pixels and render image to file
                scope.spawn(move || {
                    for _ in pbr::PbIter::new(0..bq.len()) {
                        let film_tile = pixel_rx.recv().unwrap();
                        // merge image tile into _Film_
                        film.merge_film_tile(&film_tile);
                    }
                });
            });
        }
        println!("Rendering finished");
        film.write_image(1.0 as Float / samples_per_pixel as Float);
        // TODO: Write buffers for debug visualization
    }
}
