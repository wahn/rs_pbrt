extern crate atom;
extern crate crossbeam;
extern crate num_cpus;
extern crate pbr;

// std
use std;
use std::f32::consts::PI;
use std::sync::Arc;
// others
use atom::*;
use atomic::Atomic;
// pbrt
use core::camera::{Camera, CameraSample};
use core::film::Film;
use core::geometry::{
    bnd3_expand, bnd3_union_bnd3, nrm_abs_dot_vec3, pnt3_distance_squared, vec3_abs_dot_nrm,
    vec3_max_component,
};
use core::geometry::{
    Bounds2i, Bounds3f, Normal3f, Point2f, Point2i, Point3f, Point3i, Ray, Vector2i, Vector3f,
};
use core::integrator::{compute_light_power_distribution, uniform_sample_one_light};
use core::interaction::Interaction;
use core::lowdiscrepancy::radical_inverse;
use core::material::TransportMode;
use core::parallel::AtomicFloat;
use core::pbrt::{clamp_t, lerp};
use core::pbrt::{Float, Spectrum};
use core::reflection::{Bsdf, BxdfType};
use core::sampler::{GlobalSampler, Sampler};
use core::scene::Scene;
use samplers::halton::HaltonSampler;

pub struct SPPMIntegrator {
    pub initial_search_radius: Float,
    pub n_iterations: i32,
    pub max_depth: u32,
    pub photons_per_iteration: i32,
    pub write_frequency: i32,
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
        let mut photons_per_iter: i32 = photons_per_iteration;
        if !(photons_per_iter > 0_i32) {
            let film: Arc<Film> = camera.get_film();
            photons_per_iter = film.cropped_pixel_bounds.area();
        }
        SPPMIntegrator {
            initial_search_radius: initial_search_radius,
            n_iterations: n_iterations,
            max_depth: max_depth,
            photons_per_iteration: photons_per_iter,
            write_frequency: write_frequency,
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
    pub next: Atom<Arc<SPPMPixelListNode>>,
}

impl Default for SPPMPixelListNode {
    fn default() -> SPPMPixelListNode {
        SPPMPixelListNode {
            pixel: Arc::default(),
            next: Atom::empty(),
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

fn hash(p: &Point3i, hash_size: i32) -> usize {
    let (x, _overflow) = p.x.overflowing_mul(73856093);
    let (y, _overflow) = p.y.overflowing_mul(19349663);
    let (z, _overflow) = p.z.overflowing_mul(83492791);
    ((x ^ y ^ z) as u32 % hash_size as u32) as usize
}

/// **Main function** to **render** a scene multi-threaded (using all
/// available cores) with **Stochastic Progressive Photon Mapping**
/// (SPPM).
pub fn render_sppm(
    scene: &Scene,
    camera: &Arc<Camera + Send + Sync>,
    _sampler: &mut Box<Sampler + Send + Sync>,
    integrator: &mut Box<SPPMIntegrator>,
    _num_threads: u8,
) {
    // TODO: multi-threading
    // let num_cores: usize;
    // if num_threads == 0_u8 {
    //     num_cores = num_cpus::get();
    // } else {
    //     num_cores = num_threads as usize;
    // }
    // TODO: ProfilePhase p(Prof::IntegratorRender);

    // initialize _pixel_bounds_ and _pixels_ array for SPPM
    let film: Arc<Film> = camera.get_film();
    let pixel_bounds: Bounds2i = film.cropped_pixel_bounds;
    let n_pixels: i32 = pixel_bounds.area();
    let mut pixels: Vec<Arc<SPPMPixel>> = Vec::with_capacity(n_pixels as usize);
    for _i in 0..n_pixels as usize {
        let mut pixel = SPPMPixel::default();
        pixel.radius = integrator.initial_search_radius;
        pixels.push(Arc::new(pixel));
    }
    let inv_sqrt_spp: Float = 1.0 as Float / (integrator.n_iterations as Float).sqrt();
    // TODO: let pixel_memory_bytes: usize = n_pixels as usize * std::mem::size_of::<SPPMPixel>();

    // compute _light_distr_ for sampling lights proportional to power
    if let Some(light_distr) = compute_light_power_distribution(scene) {
        // perform _n_iterations_ of SPPM integration
        let mut sampler: Box<HaltonSampler> = Box::new(HaltonSampler::new(
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
                println!("Generate SPPM visible points ...");
                for y in pbr::PbIter::new(0..n_tiles.y) {
                    for x in 0..n_tiles.x {
                        let tile: Point2i = Point2i { x: x, y: y };
                        // TODO: MemoryArena &arena = perThreadArenas[ThreadIndex];

                        // follow camera paths for _tile_ in image for SPPM
                        // TODO: let tile_index: i32 = tile.y * n_tiles.x + tile.x;
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
                            let camera_sample: CameraSample =
                                tile_sampler.get_camera_sample(&p_pixel);
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
                                                    &mut tile_sampler,
                                                    false,
                                                    None,
                                                );
                                            // possibly create visible point and end camera path
                                            let mut bsdf_flags: u8 = BxdfType::BsdfDiffuse as u8
                                                | BxdfType::BsdfReflection as u8
                                                | BxdfType::BsdfTransmission as u8;
                                            let is_diffuse: bool =
                                                bsdf.num_components(bsdf_flags) > 0;
                                            bsdf_flags = BxdfType::BsdfGlossy as u8
                                                | BxdfType::BsdfReflection as u8
                                                | BxdfType::BsdfTransmission as u8;
                                            let is_glossy: bool =
                                                bsdf.num_components(bsdf_flags) > 0;
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
                                                beta *= f * vec3_abs_dot_nrm(&wi, &isect.shading.n)
                                                    / pdf;
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
            let mut grid: Vec<Atom<Arc<SPPMPixelListNode>>> = Vec::with_capacity(hash_size);
            for _i in 0..hash_size {
                grid.push(Atom::empty());
            }
            {
                // TODO: ProfilePhase _(Prof::SPPMGridConstruction);

                // compute grid bounds for SPPM visible points
                let mut max_radius: Float = 0.0 as Float;
                println!("Compute grid bounds for SPPM visible points ...");
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
                assert!(base_grid_res > 0_i32);
                for i in 0..3 as usize {
                    grid_res[i] =
                        ((base_grid_res as Float * diag[i as u8] / max_diag).floor() as i32).max(1);
                }
                // add visible points to SPPM grid
                // TODO: ParallelFor([&](int pixelIndex) { ... }, nPixels, 4096);
                println!("Add visible points to SPPM grid ...");
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
                        for z in p_min.z..(p_max.z + 1) {
                            for y in p_min.y..(p_max.y + 1) {
                                for x in p_min.x..(p_max.x + 1) {
                                    // add visible point to grid cell $(x, y, z)$
                                    let h: usize =
                                        hash(&Point3i { x: x, y: y, z: z }, hash_size as i32);
                                    let mut node_arc: Arc<SPPMPixelListNode> =
                                        Arc::new(SPPMPixelListNode::default());
                                    let pixel_clone: Arc<SPPMPixel> = pixel.clone();
                                    if let Some(node) = Arc::get_mut(&mut node_arc) {
                                        node.pixel = pixel_clone;
                                        if !grid[h].is_none() {
                                            node.next.set_if_none(grid[h].take().unwrap());
                                        }
                                    }
                                    // atomically add _node_ to the start of _grid[h]_'s linked list
                                    grid[h].set_if_none(node_arc);
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
            {
                // TODO: ProfilePhase _(Prof::SPPMPhotonPass);
                // TODO: ParallelFor([&](int photon_index) { ... }, photonsPerIteration, 8192);
                println!("Trace photons and accumulate contributions ...");
                for photon_index in 0..integrator.photons_per_iteration as usize {
                    // MemoryArena &arena = photonShootArenas[ThreadIndex];
                    // follow photon path for _photon_index_
                    let halton_index: u64 =
                        iter as u64 * integrator.photons_per_iteration as u64 + photon_index as u64;
                    let mut halton_dim: i32 = 0;
                    // choose light to shoot photon from
                    let mut light_pdf_opt: Option<Float> = Some(0.0 as Float);
                    let light_sample: Float = radical_inverse(halton_dim as u16, halton_index);
                    halton_dim += 1;
                    let light_num: usize =
                        light_distr.sample_discrete(light_sample, light_pdf_opt.as_mut());
                    if let Some(light_pdf) = light_pdf_opt {
                        let ref light = scene.lights[light_num];
                        // compute sample values for photon ray leaving light source
                        let u_light_0: Point2f = Point2f {
                            x: radical_inverse(halton_dim as u16, halton_index),
                            y: radical_inverse((halton_dim + 1) as u16, halton_index),
                        };
                        let u_light_1: Point2f = Point2f {
                            x: radical_inverse((halton_dim + 2) as u16, halton_index),
                            y: radical_inverse((halton_dim + 3) as u16, halton_index),
                        };
                        let u_light_time: Float = lerp(
                            radical_inverse((halton_dim + 4) as u16, halton_index),
                            camera.get_shutter_open(),
                            camera.get_shutter_close(),
                        );
                        halton_dim += 5;
                        // generate _photon_ray_ from light source and initialize _beta_
                        // RayDifferential photon_ray;
                        let mut photon_ray: Ray = Ray::default();
                        let mut n_light: Normal3f = Normal3f::default();
                        // Float pdf_pos, pdf_dir;
                        let mut pdf_pos: Float = 0.0;
                        let mut pdf_dir: Float = 0.0;
                        let le: Spectrum = light.sample_le(
                            &u_light_0,
                            &u_light_1,
                            u_light_time,
                            &mut photon_ray,
                            &mut n_light,
                            &mut pdf_pos,
                            &mut pdf_dir,
                        );
                        if pdf_pos == 0.0 as Float || pdf_dir == 0.0 as Float || le.is_black() {
                            println!(
                                "light[{}]: pdf_pos = {}, pdf_dir = {}, le = {:?}",
                                light_num, pdf_pos, pdf_dir, le
                            );
                            return;
                        }
                        let mut beta: Spectrum = (le * nrm_abs_dot_vec3(&n_light, &photon_ray.d))
                            / (light_pdf * pdf_pos * pdf_dir);
                        if beta.is_black() {
                            println!("light[{}]: beta = {:?}", light_num, beta);
                            return;
                        }
                        // follow photon path through scene and record intersections
                        for depth in 0..integrator.max_depth {
                            if let Some(mut isect) = scene.intersect(&mut photon_ray) {
                                // TODO: ++totalPhotonSurfaceInteractions;
                                if depth > 0 {
                                    // add photon contribution to nearby visible points
                                    let mut photon_grid_index: Point3i = Point3i::default();
                                    if to_grid(
                                        &isect.p,
                                        &grid_bounds,
                                        &grid_res,
                                        &mut photon_grid_index,
                                    ) {
                                        let h: usize = hash(&photon_grid_index, hash_size as i32);
                                        // add photon contribution to visible points in _grid[h]_
                                        assert!(
                                            h < hash_size,
                                            "hash({:?}, {:?})",
                                            photon_grid_index,
                                            hash_size
                                        );
                                        if !grid[h].is_none() {
                                            let mut node = grid[h].take().unwrap();
                                            loop {
                                                // TODO: ++visiblePointsChecked;
                                                let pixel = node.pixel.clone();
                                                let radius: Float = pixel.radius;
                                                if pnt3_distance_squared(&pixel.vp.p, &isect.p)
                                                    > radius * radius
                                                {
                                                    if !node.next.is_none() {
                                                        node = node.next.take().unwrap();
                                                    } else {
                                                        break;
                                                    }
                                                    continue;
                                                }
                                                // update _pixel_ $\phi$ and $m$ for nearby photon
                                                let wi: Vector3f = -photon_ray.d;
                                                if let Some(ref bsdf) = pixel.vp.bsdf {
                                                    let bsdf_flags: u8 = BxdfType::BsdfAll as u8;
                                                    let phi: Spectrum = beta
                                                        * bsdf.f(&pixel.vp.wo, &wi, bsdf_flags);
                                                    for i in 0..3 {
                                                        pixel.phi[i].add(phi[i]);
                                                    }
                                                    pixel.m.fetch_add(
                                                        1_i32,
                                                        atomic::Ordering::Relaxed,
                                                    );
                                                }
                                                if !node.next.is_none() {
                                                    node = node.next.take().unwrap();
                                                } else {
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                }
                                // sample new photon ray direction

                                // compute BSDF at photon intersection point
                                let mode: TransportMode = TransportMode::Importance;
                                isect.compute_scattering_functions(
                                    &mut photon_ray, // arena,
                                    true,
                                    mode,
                                );
                                if let Some(ref photon_bsdf) = isect.bsdf {
                                    // sample BSDF _fr_ and direction _wi_ for reflected photon
                                    let mut wi: Vector3f = Vector3f::default();
                                    let wo: Vector3f = -photon_ray.d;
                                    let mut pdf: Float = 0.0;
                                    let bsdf_flags: u8 = BxdfType::BsdfAll as u8;
                                    let mut sampled_type: u8 = u8::max_value();
                                    // generate _bsdf_sample_ for outgoing photon sample
                                    let bsdf_sample: Point2f = Point2f {
                                        x: radical_inverse(halton_dim as u16, halton_index),
                                        y: radical_inverse((halton_dim + 1) as u16, halton_index),
                                    };
                                    halton_dim += 2;
                                    let fr: Spectrum = photon_bsdf.sample_f(
                                        &wo,
                                        &mut wi,
                                        &bsdf_sample,
                                        &mut pdf,
                                        bsdf_flags,
                                        &mut sampled_type,
                                    );
                                    if fr.is_black() || pdf == 0.0 as Float {
                                        break;
                                    }
                                    let bnew: Spectrum =
                                        beta * fr * vec3_abs_dot_nrm(&wi, &isect.shading.n) / pdf;
                                    // possibly terminate photon path with Russian roulette
                                    let q: Float =
                                        (0.0 as Float).max(1.0 as Float - bnew.y() / beta.y());
                                    if radical_inverse(halton_dim as u16, halton_index) < q {
                                        break;
                                    } else {
                                        halton_dim += 1;
                                    }
                                    beta = bnew / (1.0 as Float - q);
                                    photon_ray = isect.spawn_ray(&wi);
                                } else {
                                    photon_ray = isect.spawn_ray(&photon_ray.d);
                                    // --depth;
                                    continue;
                                }
                            } else {
                                break;
                            }
                        }
                        // arena.Reset();
                    }
                }
            }
            // update pixel values from this pass's photons
            {
                // TODO: ProfilePhase _(Prof::SPPMStatsUpdate);
                println!("Update pixel values from this pass's photons ...");
                // TODO: ParallelFor([&](int i) { ... }, nPixels, 4096);
                for i in 0..n_pixels as usize {
                    // copy immutable data ...
                    let p_n: Float;
                    let p_tau: Spectrum;
                    let p_radius: Float;
                    let mut p_phi: [Float; 3] = [0.0 as Float; 3];
                    let p_vp_beta: Spectrum;
                    {
                        let p: &Arc<SPPMPixel> = &pixels[i];
                        p_n = p.n;
                        p_tau = p.tau;
                        p_radius = p.radius;
                        for j in 0..3 {
                            p_phi[j] = Float::from(&p.phi[j]);
                        }
                        p_vp_beta = p.vp.beta;
                    }
                    // ... before borrowing mutably
                    let p: &mut Arc<SPPMPixel> = &mut pixels[i];
                    let p_m = p.m.load(atomic::Ordering::Relaxed);
                    if let Some(p_mut) = Arc::<SPPMPixel>::get_mut(p) {
                        if p_m > 0_i32 {
                            // update pixel photon count, search radius, and $\tau$ from photons
                            let gamma: Float = 2.0 as Float / 3.0 as Float;
                            let n_new: Float = p_n + gamma * p_m as Float;
                            let r_new: Float = p_radius * (n_new / (p_n + p_m as Float)).sqrt();
                            let mut phi: Spectrum = Spectrum::default();
                            for j in 0..3 {
                                phi[j] = p_phi[j];
                            }
                            p_mut.tau =
                                (p_tau + p_vp_beta * phi) * (r_new * r_new) / (p_radius * p_radius);
                            p_mut.n = n_new;
                            p_mut.radius = r_new;
                            p_mut.m.store(0, atomic::Ordering::Relaxed);
                            for j in 0..3 {
                                p_mut.phi[j] = AtomicFloat::new(0.0 as Float);
                            }
                        }
                        // reset _VisiblePoint_ in pixel
                        p_mut.vp.beta = Spectrum::default();
                        p_mut.vp.bsdf = None;
                    }
                }
            }
            // periodically store SPPM image in film and write image
            if iter + 1 == integrator.n_iterations || ((iter + 1) % integrator.write_frequency) == 0
            {
                let x0: i32 = pixel_bounds.p_min.x;
                let x1: i32 = pixel_bounds.p_max.x;
                let np: u64 = (iter + 1) as u64 * integrator.photons_per_iteration as u64;
                let mut image: Vec<Spectrum> = Vec::with_capacity(pixel_bounds.area() as usize);
                for y in (pixel_bounds.p_min.y as usize)..(pixel_bounds.p_max.y as usize) {
                    for x in (x0 as usize)..(x1 as usize) {
                        // compute radiance _L_ for SPPM pixel _pixel_
                        let pixel: &Arc<SPPMPixel> = &pixels[(y - pixel_bounds.p_min.y as usize)
                            * (x1 as usize - x0 as usize)
                            + (x - x0 as usize)];
                        let mut l: Spectrum = pixel.ld / (iter + 1) as Float;
                        l += pixel.tau / (np as Float * PI * pixel.radius * pixel.radius);
                        image.push(l);
                    }
                }
                film.set_image(&image[..]);
                film.write_image(1.0 as Float);
                // TODO: write SPPM radius image, if requested
                // if (getenv("SPPM_RADIUS")) {
                //     std::unique_ptr<Float[]> rimg(
                //         new Float[3 * pixel_bounds.area()]);
                //     Float minrad = 1e30f, maxrad = 0;
                //     for (int y = pixel_bounds.p_min.y; y < pixel_bounds.p_max.y; ++y) {
                //         for (int x = x0; x < x1; ++x) {
                //             const SPPMPixel &p =
                //                 pixels[(y - pixel_bounds.p_min.y) * (x1 - x0) +
                //                        (x - x0)];
                //             minrad = std::min(minrad, p.radius);
                //             maxrad = std::max(maxrad, p.radius);
                //         }
                //     }
                //     fprintf(stderr,
                //             "iterations: %d (%.2f s) radius range: %f - %f\n",
                //             iter + 1, progress.ElapsedMS() / 1000., minrad, maxrad);
                //     int offset = 0;
                //     for (int y = pixel_bounds.p_min.y; y < pixel_bounds.p_max.y; ++y) {
                //         for (int x = x0; x < x1; ++x) {
                //             const SPPMPixel &p =
                //                 pixels[(y - pixel_bounds.p_min.y) * (x1 - x0) +
                //                        (x - x0)];
                //             Float v = 1.f - (p.radius - minrad) / (maxrad - minrad);
                //             rimg[offset++] = v;
                //             rimg[offset++] = v;
                //             rimg[offset++] = v;
                //         }
                //     }
                //     Point2i res(pixel_bounds.p_max.x - pixel_bounds.p_min.x,
                //                 pixel_bounds.p_max.y - pixel_bounds.p_min.y);
                //     WriteImage("sppm_radius.png", rimg.get(), pixel_bounds, res);
                // }
            }
        }
        // TODO: progress.Done();
    }
}
