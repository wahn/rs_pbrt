// std
use std::borrow::Borrow;
use std::f32::consts::PI;
use std::sync::Arc;
// others
use atom::*;
use atomic::Atomic;
use strum::IntoEnumIterator;
// pbrt
use crate::blockqueue::BlockQueue;
use crate::core::camera::{Camera, CameraSample};
use crate::core::film::Film;
use crate::core::geometry::{
    bnd3_expand, bnd3_union_bnd3f, nrm_abs_dot_vec3f, pnt3_distance_squaredf, vec3_abs_dot_nrmf,
    vec3_max_componentf,
};
use crate::core::geometry::{
    Bounds2i, Bounds3f, Normal3f, Point2f, Point2i, Point3f, Point3i, Ray, Vector2i, Vector3f,
    XYZEnum,
};
use crate::core::integrator::{compute_light_power_distribution, uniform_sample_one_light};
use crate::core::interaction::{Interaction, SurfaceInteraction};
use crate::core::lowdiscrepancy::radical_inverse;
use crate::core::material::TransportMode;
use crate::core::parallel::AtomicFloat;
use crate::core::pbrt::{clamp_t, lerp};
use crate::core::pbrt::{Float, Spectrum};
use crate::core::reflection::{Bsdf, BxdfType};
use crate::core::scene::Scene;
use crate::core::spectrum::RGBEnum;
use crate::samplers::halton::HaltonSampler;

/// Stochastic Progressive Photon Mapping
pub struct SPPMIntegrator {
    pub camera: Arc<Camera>,
    pub initial_search_radius: Float,
    pub n_iterations: i32,
    pub max_depth: u32,
    pub photons_per_iteration: i32,
    pub write_frequency: i32,
}

impl SPPMIntegrator {
    pub fn new(
        camera: Arc<Camera>,
        n_iterations: i32,
        photons_per_iteration: i32,
        max_depth: u32,
        initial_search_radius: Float,
        write_frequency: i32,
    ) -> Self {
        let photons_per_iteration = if photons_per_iteration <= 0_i32 {
            let film: Arc<Film> = camera.get_film();
            film.cropped_pixel_bounds.area()
        } else {
            photons_per_iteration
        };
        SPPMIntegrator {
            camera,
            initial_search_radius,
            n_iterations,
            max_depth,
            photons_per_iteration,
            write_frequency,
        }
    }
    pub fn render(&self, scene: &Scene, num_threads: u8) {
        let num_cores = if num_threads == 0_u8 {
            num_cpus::get()
        } else {
            num_threads as usize
        };
        println!("Rendering with {:?} thread(s) ...", num_cores);
        // TODO: ProfilePhase p(Prof::IntegratorRender);

        // initialize _pixel_bounds_ and _pixels_ array for SPPM
        let film: Arc<Film> = self.get_camera().get_film();
        let pixel_bounds: Bounds2i = film.cropped_pixel_bounds;
        let n_pixels: i32 = pixel_bounds.area();
        let mut pixels: Vec<SPPMPixel> = Vec::with_capacity(n_pixels as usize);
        for _i in 0..n_pixels as usize {
            let mut pixel = SPPMPixel::default();
            pixel.radius = self.initial_search_radius;
            pixels.push(pixel);
        }
        let inv_sqrt_spp: Float = 1.0 as Float / (self.n_iterations as Float).sqrt();
        // TODO: let pixel_memory_bytes: usize = n_pixels as usize * std::mem::size_of::<SPPMPixel>();

        // compute _light_distr_ for sampling lights proportional to power
        if let Some(light_distr) = compute_light_power_distribution(scene) {
            // perform _n_iterations_ of SPPM integration
            let sampler: Box<HaltonSampler> = Box::new(HaltonSampler::new(
                self.n_iterations as i64,
                &pixel_bounds,
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
            for iteration in pbr::PbIter::new(0..self.n_iterations) {
                // generate SPPM visible points
                {
                    // TODO: ProfilePhase _(Prof::SPPMCameraPass);
                    // println!("Generate SPPM visible points ...");
                    {
                        let block_queue = BlockQueue::new(
                            (
                                (n_tiles.x * tile_size) as u32,
                                (n_tiles.y * tile_size) as u32,
                            ),
                            (tile_size as u32, tile_size as u32),
                            (0, 0),
                        );
                        let integrator = &self;
                        let bq = &block_queue;
                        let sampler = &sampler;
                        let pixels = &mut pixels;
                        crossbeam::scope(|scope| {
                            let (pixel_tx, pixel_rx) = crossbeam_channel::bounded(num_cores);
                            // spawn worker threads
                            for _ in 0..num_cores {
                                let pixel_tx = pixel_tx.clone();
                                scope.spawn(move |_| {
                                    while let Some((x, y)) = bq.next() {
                                        let tile: Point2i = Point2i {
                                            x: x as i32,
                                            y: y as i32,
                                        };
                                        let mut tile_bq: Vec<(i32, Spectrum, VisiblePoint)> =
                                            Vec::new();
                                        // TODO: MemoryArena &arena = perThreadArenas[ThreadIndex];

                                        // follow camera paths for _tile_ in image for SPPM
                                        // TODO: let tile_index: i32 = tile.y * n_tiles.x + tile.x;
                                        let mut tile_sampler = sampler.clone_with_seed(0_u64);
                                        // compute _tileBounds_ for SPPM tile
                                        let x0: i32 = pixel_bounds.p_min.x + tile.x * tile_size;
                                        let x1: i32 =
                                            std::cmp::min(x0 + tile_size, pixel_bounds.p_max.x);
                                        let y0: i32 = pixel_bounds.p_min.y + tile.y * tile_size;
                                        let y1: i32 =
                                            std::cmp::min(y0 + tile_size, pixel_bounds.p_max.y);
                                        let tile_bounds: Bounds2i = Bounds2i::new(
                                            Point2i { x: x0, y: y0 },
                                            Point2i { x: x1, y: y1 },
                                        );
                                        for p_pixel in &tile_bounds {
                                            // prepare _tileSampler_ for _p_pixel_
                                            tile_sampler.start_pixel(p_pixel);
                                            tile_sampler.set_sample_number(iteration as i64);
                                            // generate camera ray for pixel for SPPM
                                            let camera_sample: CameraSample =
                                                tile_sampler.get_camera_sample(p_pixel);
                                            let mut ray: Ray = Ray::default();
                                            let mut beta: Spectrum = Spectrum::new(
                                                self.get_camera().generate_ray_differential(
                                                    &camera_sample,
                                                    &mut ray,
                                                ),
                                            );
                                            if beta.is_black() {
                                                continue;
                                            }
                                            ray.scale_differentials(inv_sqrt_spp);

                                            // follow camera ray path until a visible point is created

                                            // get _SPPMPixel_ for _p_pixel_
                                            let p_pixel_o: Point2i =
                                                Point2i::from(p_pixel - pixel_bounds.p_min);
                                            let pixel_offset: i32 = p_pixel_o.x
                                                + p_pixel_o.y
                                                    * (pixel_bounds.p_max.x - pixel_bounds.p_min.x);
                                            // let mut pixel = &mut pixels[pixel_offset as usize];
                                            let mut pixel = (
                                                pixel_offset,
                                                Spectrum::default(),
                                                VisiblePoint::default(),
                                            );
                                            let mut specular_bounce: bool = false;
                                            for depth in 0..integrator.max_depth {
                                                // TODO: ++totalPhotonSurfaceInteractions;
                                                let mut isect: SurfaceInteraction =
                                                    SurfaceInteraction::default();
                                                if scene.intersect(&mut ray, &mut isect) {
                                                    // process SPPM camera ray intersection

                                                    // compute BSDF at SPPM camera ray intersection
                                                    let mode: TransportMode =
                                                        TransportMode::Radiance;
                                                    isect.compute_scattering_functions(
                                                        &ray, true, mode,
                                                    );
                                                    if let Some(bsdf) = &isect.bsdf {
                                                        // accumulate direct illumination
                                                        // at SPPM camera ray intersection
                                                        let wo: Vector3f = -ray.d;
                                                        if depth == 0 || specular_bounce {
                                                            pixel.1 += beta * isect.le(&wo);
                                                        }
                                                        let it: &SurfaceInteraction =
                                                            isect.borrow();
                                                        pixel.1 += beta
                                                            * uniform_sample_one_light(
                                                                it,
                                                                scene,
                                                                &mut tile_sampler,
                                                                false,
                                                                None,
                                                            );
                                                        // possibly create visible point and end camera path
                                                        let mut bsdf_flags: u8 =
                                                            BxdfType::BsdfDiffuse as u8
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
                                                            || (is_glossy
                                                                && depth
                                                                    == integrator.max_depth - 1)
                                                        {
                                                            pixel.2.p = isect.common.p;
                                                            pixel.2.wo = wo;
                                                            pixel.2.bsdf = Some(bsdf.clone());
                                                            pixel.2.beta = beta;
                                                            break;
                                                        }
                                                        // spawn ray from SPPM camera path vertex
                                                        if depth < integrator.max_depth - 1 {
                                                            let mut wi: Vector3f =
                                                                Vector3f::default();
                                                            let mut pdf: Float = 0.0;
                                                            let bsdf_flags: u8 =
                                                                BxdfType::BsdfAll as u8;
                                                            let mut sampled_type: u8 =
                                                                u8::max_value(); // != 0
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
                                                            beta *= f * vec3_abs_dot_nrmf(
                                                                &wi,
                                                                &isect.shading.n,
                                                            ) / pdf;
                                                            if beta.y() < 0.25 as Float {
                                                                let continue_prob: Float =
                                                                    (1.0 as Float).min(beta.y());
                                                                if tile_sampler.get_1d()
                                                                    > continue_prob
                                                                {
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
                                                        pixel.1 += beta * light.le(&mut ray);
                                                    }
                                                    break;
                                                }
                                            }
                                            tile_bq.push(pixel);
                                        }
                                        // send progress through the channel to main thread
                                        pixel_tx
                                            .send(tile_bq)
                                            .unwrap_or_else(|_| panic!("Failed to send progress"));
                                    }
                                });
                            }
                            // spawn thread to collect
                            scope.spawn(move |_| {
                                for _ in 0..bq.len() {
                                    let tile = pixel_rx.recv().unwrap();
                                    for (pixel_offset, ld, vp) in tile {
                                        let mut pixel = &mut pixels[pixel_offset as usize];
                                        pixel.ld += ld;
                                        pixel.vp.p = vp.p;
                                        pixel.vp.wo = vp.wo;
                                        pixel.vp.bsdf = vp.bsdf;
                                        pixel.vp.beta = vp.beta;
                                    }
                                }
                            });
                        })
                        .unwrap();
                    }
                }
                // create grid of all SPPM visible points
                let mut grid_res: [i32; 3] = [0; 3];
                let mut grid_bounds: Bounds3f = Bounds3f::default();
                // allocate grid for SPPM visible points
                let hash_size: usize = n_pixels as usize;
                let mut grid: Vec<Atom<Arc<SPPMPixelListNode>>> = Vec::with_capacity(hash_size);
                let mut grid_once: Vec<AtomSetOnce<Arc<SPPMPixelListNode>>> =
                    Vec::with_capacity(hash_size);
                for _i in 0..hash_size {
                    grid.push(Atom::empty());
                    grid_once.push(AtomSetOnce::empty());
                }
                {
                    // TODO: ProfilePhase _(Prof::SPPMGridConstruction);

                    // compute grid bounds for SPPM visible points
                    let mut max_radius: Float = 0.0 as Float;
                    // println!("Compute grid bounds for SPPM visible points ...");
                    for pixel in pixels.iter().take(n_pixels as usize) {
                        if !pixel.vp.beta.is_black() {
                            let vp_bound: Bounds3f = bnd3_expand(
                                &Bounds3f {
                                    p_min: pixel.vp.p,
                                    p_max: pixel.vp.p,
                                },
                                pixel.radius,
                            );
                            grid_bounds = bnd3_union_bnd3f(&grid_bounds, &vp_bound);
                            max_radius = max_radius.max(pixel.radius);
                        }
                    }
                    // compute resolution of SPPM grid in each dimension
                    let diag: Vector3f = grid_bounds.diagonal();
                    let max_diag: Float = vec3_max_componentf(&diag);
                    let base_grid_res: i32 = (max_diag / max_radius).floor() as i32;
                    assert!(base_grid_res > 0_i32);
                    for i in XYZEnum::iter() {
                        grid_res[i as usize] =
                            ((base_grid_res as Float * diag[i] / max_diag).floor() as i32).max(1);
                    }
                    // add visible points to SPPM grid
                    // println!("Add visible points to SPPM grid ...");
                    let chunk_size: usize = (n_pixels / num_cores as i32) as usize;
                    {
                        let bands: Vec<&mut [SPPMPixel]> = pixels.chunks_mut(chunk_size).collect();
                        let grid = &grid;
                        crossbeam::scope(|scope| {
                            let (band_tx, band_rx) = crossbeam_channel::bounded(num_cores);
                            // spawn worker threads
                            for (b, band) in bands.into_iter().enumerate() {
                                let band_tx = band_tx.clone();
                                scope.spawn(move |_| {
                                    for pixel in band.iter_mut() {
                                        // for pixel_index in 0..n_pixels as usize {
                                        // let pixel = &pixels[pixel_index];
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
                                            for z in p_min.z..=p_max.z {
                                                for y in p_min.y..=p_max.y {
                                                    for x in p_min.x..=p_max.x {
                                                        // add visible point to grid cell $(x, y, z)$
                                                        let h: usize = hash(
                                                            &Point3i { x, y, z },
                                                            hash_size as i32,
                                                        );
                                                        let node_arc =
                                                            Arc::new(SPPMPixelListNode::new(pixel));
                                                        let old_opt =
                                                            grid[h].swap(node_arc.clone());
                                                        if let Some(old) = old_opt {
                                                            node_arc.next.set_if_none(old);
                                                        }
                                                    }
                                                }
                                            }
                                            // ReportValue(grid_cells_per_visible_point,
                                            //             (1 + pMax.x - pMin.x) * (1 + pMax.y - pMin.y) *
                                            //                 (1 + pMax.z - pMin.z));
                                        }
                                    }
                                });
                                // send progress through the channel to main thread
                                band_tx
                                    .send(b)
                                    .unwrap_or_else(|_| panic!("Failed to send progress"));
                            }
                            // spawn thread to report progress
                            scope.spawn(move |_| {
                                for _ in 0..num_cores {
                                    band_rx.recv().unwrap();
                                }
                            });
                        })
                        .unwrap();
                    }
                }
                // trace photons and accumulate contributions
                for h in 0..hash_size {
                    // take
                    let opt = grid[h].take();
                    if let Some(p) = opt {
                        grid_once[h].set_if_none(p);
                    }
                }
                std::mem::drop(grid);
                {
                    // TODO: ProfilePhase _(Prof::SPPMPhotonPass);
                    // println!("Trace photons and accumulate contributions ...");
                    let chunk_size: usize =
                        (self.photons_per_iteration / num_cores as i32) as usize;
                    {
                        let photons_vec: Vec<i32> = (0..self.photons_per_iteration).collect();
                        let bands: Vec<&[i32]> = photons_vec.chunks(chunk_size).collect();
                        let grid_once = &grid_once;
                        let integrator = &self;
                        let light_distr = &light_distr;
                        crossbeam::scope(|scope| {
                        let (band_tx, band_rx) = crossbeam_channel::bounded(num_cores);
                        // spawn worker threads
                        for (b, band) in bands.into_iter().enumerate() {
                            let band_tx = band_tx.clone();
                            scope.spawn(move |_| {
                                for photon_index in band.iter() {
                                    // for photon_index in 0..integrator.photons_per_iteration as usize {
                                    // MemoryArena &arena = photonShootArenas[ThreadIndex];
                                    // follow photon path for _photon_index_
                                    let halton_index: u64 = iteration as u64
                                        * integrator.photons_per_iteration as u64
                                        + *photon_index as u64;
                                    let mut halton_dim: i32 = 0;
                                    // choose light to shoot photon from
                                    let mut light_pdf_opt: Option<Float> = Some(0.0 as Float);
                                    let light_sample: Float =
                                        radical_inverse(halton_dim as u16, halton_index);
                                    halton_dim += 1;
                                    let light_num: usize = light_distr
                                        .sample_discrete(light_sample, light_pdf_opt.as_mut());
                                    if let Some(light_pdf) = light_pdf_opt {
                                        let light = &scene.lights[light_num];
                                        // compute sample values for photon ray leaving light source
                                        let u_light_0: Point2f = Point2f {
                                            x: radical_inverse(halton_dim as u16, halton_index),
                                            y: radical_inverse(
                                                (halton_dim + 1) as u16,
                                                halton_index,
                                            ),
                                        };
                                        let u_light_1: Point2f = Point2f {
                                            x: radical_inverse(
                                                (halton_dim + 2) as u16,
                                                halton_index,
                                            ),
                                            y: radical_inverse(
                                                (halton_dim + 3) as u16,
                                                halton_index,
                                            ),
                                        };
                                        let u_light_time: Float = lerp(
                                            radical_inverse((halton_dim + 4) as u16, halton_index),
                                            self.get_camera().get_shutter_open(),
                                            self.get_camera().get_shutter_close(),
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
                                            u_light_0,
                                            u_light_1,
                                            u_light_time,
                                            &mut photon_ray,
                                            &mut n_light,
                                            &mut pdf_pos,
                                            &mut pdf_dir,
                                        );
                                        if pdf_pos == 0.0 as Float
                                            || pdf_dir == 0.0 as Float
                                            || le.is_black()
                                        {
                                            // println!(
                                            //     "light[{}]: pdf_pos = {}, pdf_dir = {}, le = {:?}",
                                            //     light_num, pdf_pos, pdf_dir, le
                                            // );
                                            // C++: return; (from ParallelFor(...{}, photonsPerIteration, 8192);)
                                            break;
                                        }
                                        let mut beta: Spectrum = (le
                                            * nrm_abs_dot_vec3f(&n_light, &photon_ray.d))
                                            / (light_pdf * pdf_pos * pdf_dir);
                                        if beta.is_black() {
                                            // println!("light[{}]: beta = {:?}", light_num, beta);
                                            // C++:  return; (from ParallelFor(...{}, photonsPerIteration, 8192);)
                                            break;
                                        }
                                        // follow photon path through scene and record intersections
                                        for depth in 0..integrator.max_depth {
					    let mut isect: SurfaceInteraction = SurfaceInteraction::default();
					    if scene.intersect(&mut photon_ray, &mut isect) {
                                                // TODO: ++totalPhotonSurfaceInteractions;
                                                if depth > 0 {
                                                    // add photon contribution to nearby visible points
                                                    let mut photon_grid_index: Point3i =
                                                        Point3i::default();
                                                    if to_grid(
                                                        &isect.common.p,
                                                        &grid_bounds,
                                                        &grid_res,
                                                        &mut photon_grid_index,
                                                    ) {
                                                        let h: usize = hash(
                                                            &photon_grid_index,
                                                            hash_size as i32,
                                                        );
                                                        // add photon contribution to visible points in _grid[h]_
                                                        assert!(
                                                            h < hash_size,
                                                            "hash({:?}, {:?})",
                                                            photon_grid_index,
                                                            hash_size
                                                        );
                                                        if !grid_once[h].is_none() {
                                                            let mut opt = grid_once[h].get();
                                                            while let Some(node) = opt {
                                                                // deal with linked list
                                                                let pixel = node.pixel;
                                                                let radius: Float = pixel.radius;
                                                                    if pnt3_distance_squaredf(
                                                                        &pixel.vp.p,
                                                                        &isect.common.p,
                                                                    ) > radius * radius
                                                                    {
                                                                        // update opt
                                                                        opt = node.next.get();
                                                                    } else {
                                                                        // update
                                                                        // _pixel_
                                                                        // $\phi$
                                                                        // and
                                                                        // $m$
                                                                        // for
                                                                        // nearby
                                                                        // photon
                                                                        let wi: Vector3f =
                                                                            -photon_ray.d;
                                                                        if let Some(ref bsdf) =
                                                                            pixel.vp.bsdf
                                                                        {
                                                                            let bsdf_flags: u8 =
                                                                                BxdfType::BsdfAll
                                                                                    as u8;
                                                                            let phi: Spectrum = beta
                                                                                * bsdf.f(
                                                                                    &pixel.vp.wo,
                                                                                    &wi,
                                                                                    bsdf_flags,
                                                                                );
                                                                            for i in 0..3 {
                                                                                let rgb_i: RGBEnum =
                                                                                    match i {
                                                                                        0 =>
                                                                                            RGBEnum::Red,
                                                                                        1 =>
                                                                                            RGBEnum::Green,
                                                                                        _ =>
                                                                                            RGBEnum::Blue,
                                                                                    };
                                                                                let phi_i: Float = phi[rgb_i];
                                                                                pixel.phi[i]
                                                                                    .add(phi_i);
                                                                            }
                                                                            pixel.m.fetch_add(
                                                                                1_i32,
                                                                                atomic::Ordering::Relaxed,
                                                                            );
                                                                        }
                                                                        // update opt
                                                                        opt = node.next.get();
                                                                    }
                                                        }
                                                        }
                                                    }
                                                }
                                                // sample new photon ray direction

                                                // compute BSDF at photon intersection point
                                                let mode: TransportMode = TransportMode::Importance;
						isect.compute_scattering_functions(&photon_ray, true, mode);
                                                if let Some(ref photon_bsdf) = isect.bsdf {
                                                    // sample BSDF _fr_ and direction _wi_ for reflected photon
                                                    let mut wi: Vector3f = Vector3f::default();
                                                    let wo: Vector3f = -photon_ray.d;
                                                    let mut pdf: Float = 0.0;
                                                    let bsdf_flags: u8 = BxdfType::BsdfAll as u8;
                                                    let mut sampled_type: u8 = u8::max_value();
                                                    // generate _bsdf_sample_ for outgoing photon sample
                                                    let bsdf_sample: Point2f = Point2f {
                                                        x: radical_inverse(
                                                            halton_dim as u16,
                                                            halton_index,
                                                        ),
                                                        y: radical_inverse(
                                                            (halton_dim + 1) as u16,
                                                            halton_index,
                                                        ),
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
                                                    let bnew: Spectrum = beta
                                                        * fr
                                                        * vec3_abs_dot_nrmf(&wi, &isect.shading.n)
                                                        / pdf;
                                                    // possibly terminate photon path with Russian roulette
                                                    let q: Float = (0.0 as Float)
                                                        .max(1.0 as Float - bnew.y() / beta.y());
                                                    if radical_inverse(
                                                        halton_dim as u16,
                                                        halton_index,
                                                    ) < q
                                                    {
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
                                    }
                                }
                            });
                            // send progress through the channel to main thread
                            band_tx.send(b).unwrap_or_else(|_| panic!("Failed to send progress"));
                        }
                        // spawn thread to report progress
                        scope.spawn(move |_| {
                            for _ in 0..num_cores {
                                band_rx.recv().unwrap();
                            }
                        });
                    })
                    .unwrap();
                    }
                }
                // update pixel values from this pass's photons
                {
                    // TODO: ProfilePhase _(Prof::SPPMStatsUpdate);
                    // println!("Update pixel values from this pass's photons ...");
                    let chunk_size: usize = (n_pixels / num_cores as i32) as usize;
                    {
                        let bands: Vec<&mut [SPPMPixel]> = pixels.chunks_mut(chunk_size).collect();
                        crossbeam::scope(|scope| {
                            let (band_tx, band_rx) = crossbeam_channel::bounded(num_cores);
                            // spawn worker threads
                            for (b, band) in bands.into_iter().enumerate() {
                                let band_tx = band_tx.clone();
                                scope.spawn(move |_| {
                                    for p in band.iter_mut() {
                                        // let mut p = &mut pixels[i];
                                        let p_m = p.m.load(atomic::Ordering::Relaxed);
                                        if p_m > 0_i32 {
                                            // update pixel photon count, search radius, and $\tau$ from photons
                                            let gamma: Float = 2.0 as Float / 3.0 as Float;
                                            let n_new: Float = p.n + gamma * p_m as Float;
                                            let r_new: Float =
                                                p.radius * (n_new / (p.n + p_m as Float)).sqrt();
                                            let mut phi: Spectrum = Spectrum::default();
                                            for j in 0..3 {
                                                match j {
                                                    0 => {
                                                        phi[RGBEnum::Red] = Float::from(&p.phi[j]);
                                                    }
                                                    1 => {
                                                        phi[RGBEnum::Green] =
                                                            Float::from(&p.phi[j]);
                                                    }
                                                    _ => {
                                                        phi[RGBEnum::Blue] = Float::from(&p.phi[j]);
                                                    }
                                                }
                                            }
                                            p.tau = (p.tau + p.vp.beta * phi) * (r_new * r_new)
                                                / (p.radius * p.radius);
                                            p.n = n_new;
                                            p.radius = r_new;
                                            p.m.store(0, atomic::Ordering::Relaxed);
                                            for j in 0..3 {
                                                p.phi[j] = AtomicFloat::new(0.0 as Float);
                                            }
                                        }
                                        // reset _VisiblePoint_ in pixel
                                        p.vp.beta = Spectrum::default();
                                        p.vp.bsdf = None;
                                    }
                                });
                                // send progress through the channel to main thread
                                band_tx
                                    .send(b)
                                    .unwrap_or_else(|_| panic!("Failed to send progress"));
                            }
                            // spawn thread to report progress
                            scope.spawn(move |_| {
                                for _ in 0..num_cores {
                                    band_rx.recv().unwrap();
                                }
                            });
                        })
                        .unwrap();
                    }
                }
                // periodically store SPPM image in film and write image
                if iteration + 1 == self.n_iterations
                    || ((iteration + 1) % self.write_frequency) == 0
                {
                    let x0: i32 = pixel_bounds.p_min.x;
                    let x1: i32 = pixel_bounds.p_max.x;
                    let np: u64 = (iteration + 1) as u64 * self.photons_per_iteration as u64;
                    let mut image: Vec<Spectrum> = Vec::with_capacity(pixel_bounds.area() as usize);
                    for y in (pixel_bounds.p_min.y as usize)..(pixel_bounds.p_max.y as usize) {
                        for x in (x0 as usize)..(x1 as usize) {
                            // compute radiance _L_ for SPPM pixel _pixel_
                            let pixel = &pixels[(y - pixel_bounds.p_min.y as usize)
                                * (x1 as usize - x0 as usize)
                                + (x - x0 as usize)];
                            let mut l: Spectrum = pixel.ld / (iteration + 1) as Float;
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
    pub fn get_camera(&self) -> Arc<Camera> {
        self.camera.clone()
    }
}

#[derive(Default)]
pub struct VisiblePoint {
    pub p: Point3f,
    pub wo: Vector3f,
    pub bsdf: Option<Bsdf>,
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

pub struct SPPMPixelListNode<'p> {
    pub pixel: &'p SPPMPixel,
    pub next: AtomSetOnce<Arc<SPPMPixelListNode<'p>>>,
}

impl<'p> SPPMPixelListNode<'p> {
    pub fn new(pixel: &'p SPPMPixel) -> Self {
        SPPMPixelListNode {
            pixel,
            next: AtomSetOnce::empty(),
        }
    }
}

fn to_grid(p: &Point3f, bounds: &Bounds3f, grid_res: &[i32; 3], pi: &mut Point3i) -> bool {
    let mut in_bounds: bool = true;
    let pg: Vector3f = bounds.offset(p);
    for i in XYZEnum::iter() {
        (*pi)[i] = (grid_res[i as usize] as Float * pg[i]) as i32;
        in_bounds &= (*pi)[i] >= 0 && (*pi)[i] < grid_res[i as usize];
        (*pi)[i] = clamp_t((*pi)[i], 0, grid_res[i as usize] - 1);
    }
    in_bounds
}

fn hash(p: &Point3i, hash_size: i32) -> usize {
    let (x, _overflow) = p.x.overflowing_mul(73_856_093);
    let (y, _overflow) = p.y.overflowing_mul(19_349_663);
    let (z, _overflow) = p.z.overflowing_mul(83_492_791);
    ((x ^ y ^ z) as u32 % hash_size as u32) as usize
}
