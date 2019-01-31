extern crate atom;
extern crate crossbeam;
extern crate num_cpus;

// std
use std;
use std::sync::atomic::Ordering;
use std::sync::{Arc, Barrier};
// others
use atom::*;
use atomic::Atomic;
// pbrt
use core::camera::{Camera, CameraSample};
use core::geometry::{bnd3_expand, bnd3_union_bnd3, vec3_abs_dot_nrm, vec3_max_component};
use core::geometry::{Bounds2i, Bounds3f, Point2i, Point3f, Point3i, Ray, Vector2i, Vector3f};
use core::integrator::{compute_light_power_distribution, uniform_sample_one_light};
use core::interaction::Interaction;
use core::material::TransportMode;
use core::parallel::AtomicFloat;
use core::pbrt::clamp_t;
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

pub struct SPPMPixelListNode {
    pub pixel: Arc<SPPMPixel>,
    pub next: AtomSetOnce<Box<SPPMPixelListNode>>,
}

impl Default for SPPMPixelListNode {
    fn default() -> SPPMPixelListNode {
        SPPMPixelListNode {
            pixel: Arc::default(),
            next: AtomSetOnce::empty(),
        }
    }
}

fn to_grid(p: &Point3f, bounds: &Bounds3f, grid_res: &[i32; 3], pi: &mut Point3i) -> bool {
    let mut in_bounds: bool = true;
    let pg: Vector3f = bounds.offset(p);
    for i in 0..3 as u8 {
        (*pi)[i] = (grid_res[i as usize] as Float * pg[i]) as i32;
        in_bounds &= (*pi)[i] >= 0 && (*pi)[i] < grid_res[i as usize];
        (*pi)[i] = clamp_t((*pi)[i], 0, grid_res[i as usize] - 1);
    }
    in_bounds
}

fn hash(p: &Point3i, hash_size: usize) -> usize {
    (((p.x * 73856093) ^ (p.y * 19349663) ^ (p.z * 83492791)) % hash_size as i32) as usize
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
    let mut pixels: Vec<Arc<SPPMPixel>> = Vec::with_capacity(n_pixels as usize);
    for i in 0..n_pixels as usize {
        let mut pixel = SPPMPixel::default();
        pixel.radius = integrator.initial_search_radius;
        pixels.push(Arc::new(pixel));
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
            // TODO: ParallelFor2D([&](Point2i tile) { ... }, nTiles);
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
                        let mut beta: Spectrum = Spectrum::new(
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
                        if let Some(pixel) = Arc::get_mut(&mut pixels[pixel_offset as usize]) {
                            let mut specular_bounce: bool = false;
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
                                        // spawn ray from SPPM camera path vertex
                                        if depth < integrator.max_depth - 1 {
                                            let mut wi: Vector3f = Vector3f::default();
                                            let mut pdf: Float = 0.0;
                                            let bsdf_flags: u8 = BxdfType::BsdfAll as u8;
                                            let mut sampled_type: u8 = u8::max_value(); // != 0
                                            let f: Spectrum = bsdf.sample_f(
                                                &wo,
                                                &mut wi,
                                                &tile_sampler.get_2d(),
                                                &mut pdf,
                                                bsdf_flags,
                                                &mut sampled_type,
                                            );
                                            if pdf == 0.0 as Float || f.is_black() {
                                                break;
                                            }
                                            specular_bounce = sampled_type
                                                & (BxdfType::BsdfSpecular as u8)
                                                != 0_u8;
                                            beta *=
                                                f * vec3_abs_dot_nrm(&wi, &isect.shading.n) / pdf;
                                            if beta.y() < 0.25 as Float {
                                                let continue_prob: Float =
                                                    (1.0 as Float).min(beta.y());
                                                if tile_sampler.get_1d() > continue_prob {
                                                    break;
                                                }
                                                beta /= continue_prob;
                                            }
                                            ray = isect.spawn_ray(&wi);
                                        }
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
        }
        // create grid of all SPPM visible points
        let mut grid_res: [i32; 3] = [0; 3];
        let mut grid_bounds: Bounds3f = Bounds3f::default();
        // allocate grid for SPPM visible points
        let hash_size: usize = n_pixels as usize;
        let mut grid: Vec<AtomSetOnce<Box<SPPMPixelListNode>>> = Vec::with_capacity(hash_size);
        for _i in 0..hash_size {
            grid.push(AtomSetOnce::empty());
        }
        {
            // TODO: ProfilePhase _(Prof::SPPMGridConstruction);

            // compute grid bounds for SPPM visible points
            let mut max_radius: Float = 0.0 as Float;
            for i in 0..n_pixels as usize {
                if let Some(pixel) = Arc::get_mut(&mut pixels[i]) {
                    if pixel.vp.beta.is_black() {
                        continue;
                    }
                    let vp_bound: Bounds3f = bnd3_expand(
                        &Bounds3f {
                            p_min: pixel.vp.p,
                            p_max: pixel.vp.p,
                        },
                        pixel.radius,
                    );
                    grid_bounds = bnd3_union_bnd3(&grid_bounds, &vp_bound);
                    max_radius = max_radius.max(pixel.radius);
                }
            }
            // compute resolution of SPPM grid in each dimension
            let diag: Vector3f = grid_bounds.diagonal();
            let max_diag: Float = vec3_max_component(&diag);
            let base_grid_res: i32 = (max_diag / max_radius).floor() as i32;
            println!(
                "base_grid_res: {} ({}/{})",
                base_grid_res, max_diag, max_radius
            );
            assert!(base_grid_res > 0_i32);
            for i in 0..3 as usize {
                grid_res[i] =
                    ((base_grid_res as Float * diag[i as u8] / max_diag).floor() as i32).max(1);
            }
            // add visible points to SPPM grid
            // TODO: ParallelFor([&](int pixelIndex) { ... }, nPixels, 4096);
            for pixel_index in 0..n_pixels as usize {
                let pixel: &Arc<SPPMPixel> = &pixels[pixel_index];
                if !pixel.vp.beta.is_black() {
                    // add pixel's visible point to applicable grid cells
                    let radius: Float = pixel.radius;
                    let mut p_min: Point3i = Point3i::default();
                    let mut p_max: Point3i = Point3i::default();
                    to_grid(
                        &(pixel.vp.p
                            - Vector3f {
                                x: radius,
                                y: radius,
                                z: radius,
                            }),
                        &grid_bounds,
                        &grid_res,
                        &mut p_min,
                    );
                    to_grid(
                        &(pixel.vp.p
                            + Vector3f {
                                x: radius,
                                y: radius,
                                z: radius,
                            }),
                        &grid_bounds,
                        &grid_res,
                        &mut p_max,
                    );
                    for z in p_min.z..p_max.z {
                        for y in p_min.y..p_max.y {
                            for x in p_min.x..p_max.x {
                                // add visible point to grid cell $(x, y, z)$
                                let h: usize = hash(&Point3i { x: x, y: y, z: z }, hash_size);
                                let mut node_arc: Arc<SPPMPixelListNode> =
                                    Arc::new(SPPMPixelListNode::default());
                                let pixel_clone: Arc<SPPMPixel> = pixel.clone();
                                if let Some(node) = Arc::get_mut(&mut node_arc) {
                                    node.pixel = pixel_clone;
                                    // atomically add _node_ to the start of _grid[h]_'s linked list
                                    // node.next = grid[h].clone();
                                    // let new = node;
                                    // let mut old = grid[h].load(Ordering::Relaxed);
                                    // loop {
                                    //     match grid[h].compare_exchange_weak(old, new, Ordering::SeqCst, Ordering::Relaxed) {
                                    //         Ok(_) => break,
                                    //         Err(x) => old = x,
                                    //     }
                                    // }
                                    // WORK
                                }
                            }
                        }
                    }
                    // ReportValue(grid_cells_per_visible_point,
                    //             (1 + pMax.x - pMin.x) * (1 + pMax.y - pMin.y) *
                    //                 (1 + pMax.z - pMin.z));
                }
            }
        }
        // trace photons and accumulate contributions
        // update pixel values from this pass's photons
        // periodically store SPPM image in film and write image
        // WORK
    }
}
