//! Rendering an image of the scene is handled by an instance of a
//! class that implements the **Integrator** interface.

// std
use std::sync::Arc;
// pbrt
use crate::blockqueue::BlockQueue;
use crate::core::camera::{Camera, CameraSample};
use crate::core::geometry::{pnt2_inside_exclusive, vec3_abs_dot_nrm};
use crate::core::geometry::{Bounds2i, Point2f, Point2i, Ray, Vector2i, Vector3f};
use crate::core::interaction::{Interaction, InteractionCommon, SurfaceInteraction};
use crate::core::light::is_delta_light;
use crate::core::light::{Light, VisibilityTester};
use crate::core::pbrt::{Float, Spectrum};
use crate::core::reflection::BxdfType;
use crate::core::sampler::Sampler;
use crate::core::sampling::power_heuristic;
use crate::core::sampling::Distribution1D;
use crate::core::scene::Scene;
use crate::integrators::ao::AOIntegrator;
use crate::integrators::bdpt::BDPTIntegrator;
use crate::integrators::directlighting::DirectLightingIntegrator;
use crate::integrators::mlt::MLTIntegrator;
use crate::integrators::path::PathIntegrator;
use crate::integrators::sppm::SPPMIntegrator;
use crate::integrators::volpath::VolPathIntegrator;
use crate::integrators::whitted::WhittedIntegrator;

// see integrator.h

pub enum Integrator {
    BDPT(BDPTIntegrator),
    MLT(MLTIntegrator),
    SPPM(SPPMIntegrator),
    Sampler(SamplerIntegrator),
}

impl Integrator {
    pub fn render(&mut self, scene: &Scene, num_threads: u8) {
        match self {
            Integrator::BDPT(integrator) => integrator.render(scene, num_threads),
            Integrator::MLT(integrator) => integrator.render(scene, num_threads),
            Integrator::SPPM(integrator) => integrator.render(scene, num_threads),
            Integrator::Sampler(integrator) => integrator.render(scene, num_threads),
        }
    }
}

pub enum SamplerIntegrator {
    AO(AOIntegrator),
    DirectLighting(DirectLightingIntegrator),
    Path(PathIntegrator),
    VolPath(VolPathIntegrator),
    Whitted(WhittedIntegrator),
}

impl SamplerIntegrator {
    pub fn preprocess(&mut self, scene: &Scene) {
        match self {
            SamplerIntegrator::AO(integrator) => integrator.preprocess(scene),
            SamplerIntegrator::DirectLighting(integrator) => integrator.preprocess(scene),
            SamplerIntegrator::Path(integrator) => integrator.preprocess(scene),
            SamplerIntegrator::VolPath(integrator) => integrator.preprocess(scene),
            SamplerIntegrator::Whitted(integrator) => integrator.preprocess(scene),
        }
    }
    pub fn render(&mut self, scene: &Scene, num_threads: u8) {
        match self {
            _ => {
                let film = self.get_camera().get_film();
                let sample_bounds: Bounds2i = film.get_sample_bounds();
                self.preprocess(scene);
                let sample_extent: Vector2i = sample_bounds.diagonal();
                let tile_size: i32 = 16;
                let x: i32 = (sample_extent.x + tile_size - 1) / tile_size;
                let y: i32 = (sample_extent.y + tile_size - 1) / tile_size;
                let n_tiles: Point2i = Point2i { x, y };
                // TODO: ProgressReporter reporter(nTiles.x * nTiles.y, "Rendering");
                let num_cores = if num_threads == 0_u8 {
                    num_cpus::get()
                } else {
                    num_threads as usize
                };
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
                    let integrator = &self;
                    let bq = &block_queue;
                    let sampler = &self.get_sampler();
                    let camera = &self.get_camera();
                    let film = &film;
                    let pixel_bounds = &self.get_pixel_bounds();
                    crossbeam::scope(|scope| {
                        let (pixel_tx, pixel_rx) = crossbeam_channel::bounded(num_cores);
                        // spawn worker threads
                        for _ in 0..num_cores {
                            let pixel_tx = pixel_tx.clone();
                            let mut tile_sampler: Arc<Sampler> =
                                sampler.clone_with_seed(0_u64);
                            scope.spawn(move |_| {
                                while let Some((x, y)) = bq.next() {
                                    let tile: Point2i = Point2i {
                                        x: x as i32,
                                        y: y as i32,
                                    };
                                    let seed: i32 = tile.y * n_tiles.x + tile.x;
                                    Arc::get_mut(&mut tile_sampler).unwrap().reseed(seed as u64);
                                    let x0: i32 = sample_bounds.p_min.x + tile.x * tile_size;
                                    let x1: i32 =
                                        std::cmp::min(x0 + tile_size, sample_bounds.p_max.x);
                                    let y0: i32 = sample_bounds.p_min.y + tile.y * tile_size;
                                    let y1: i32 =
                                        std::cmp::min(y0 + tile_size, sample_bounds.p_max.y);
                                    let tile_bounds: Bounds2i = Bounds2i::new(
                                        Point2i { x: x0, y: y0 },
                                        Point2i { x: x1, y: y1 },
                                    );
                                    // println!("Starting image tile {:?}", tile_bounds);
                                    let mut film_tile = film.get_film_tile(&tile_bounds);
                                    for pixel in &tile_bounds {
                                        Arc::get_mut(&mut tile_sampler).unwrap().start_pixel(pixel);
                                        if !pnt2_inside_exclusive(pixel, &pixel_bounds) {
                                            continue;
                                        }
                                        let mut done: bool = false;
                                        while !done {
                                            // let's use the copy_arena crate instead of pbrt's MemoryArena
                                            // let mut arena: Arena = Arena::with_capacity(262144); // 256kB

                                            // initialize _CameraSample_ for current sample
                                            let camera_sample: CameraSample =
                                                Arc::get_mut(&mut tile_sampler).unwrap().get_camera_sample(pixel);
                                            // generate camera ray for current sample
                                            let mut ray: Ray = Ray::default();
                                            let ray_weight: Float = camera
                                                .generate_ray_differential(
                                                    &camera_sample,
                                                    &mut ray,
                                                );
                                            ray.scale_differentials(
                                                1.0 as Float
                                                    / (tile_sampler.get_samples_per_pixel()
                                                        as Float)
                                                        .sqrt(),
                                            );
                                            // TODO: ++nCameraRays;
                                            // evaluate radiance along camera ray
                                            let mut l: Spectrum = Spectrum::new(0.0 as Float);
                                            let y: Float = l.y();
                                            if ray_weight > 0.0 {
                                                l = integrator.li(
                                                    &mut ray,
                                                    scene,
                                                    Arc::get_mut(&mut tile_sampler).unwrap(), // &mut arena,
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
                                            film_tile.add_sample(
                                                camera_sample.p_film,
                                                &mut l,
                                                ray_weight,
                                            );
                                            done = !Arc::get_mut(&mut tile_sampler).unwrap().start_next_sample();
                                        } // arena is dropped here !
                                    }
                                    // send the tile through the channel to main thread
                                    pixel_tx
                                        .send(film_tile)
                                        .unwrap_or_else(|_| panic!("Failed to send tile"));
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
        }
    }
    pub fn li(&self, ray: &mut Ray, scene: &Scene, sampler: &mut Sampler, depth: i32) -> Spectrum {
        match self {
            SamplerIntegrator::AO(integrator) => integrator.li(ray, scene, sampler, depth),
            SamplerIntegrator::DirectLighting(integrator) => {
                integrator.li(ray, scene, sampler, depth)
            }
            SamplerIntegrator::Path(integrator) => integrator.li(ray, scene, sampler, depth),
            SamplerIntegrator::VolPath(integrator) => integrator.li(ray, scene, sampler, depth),
            SamplerIntegrator::Whitted(integrator) => integrator.li(ray, scene, sampler, depth),
        }
    }
    pub fn get_camera(&self) -> Arc<Camera> {
        match self {
            SamplerIntegrator::AO(integrator) => integrator.get_camera(),
            SamplerIntegrator::DirectLighting(integrator) => integrator.get_camera(),
            SamplerIntegrator::Path(integrator) => integrator.get_camera(),
            SamplerIntegrator::VolPath(integrator) => integrator.get_camera(),
            SamplerIntegrator::Whitted(integrator) => integrator.get_camera(),
        }
    }
    pub fn get_sampler(&self) -> Arc<Sampler> {
        match self {
            SamplerIntegrator::AO(integrator) => integrator.get_sampler(),
            SamplerIntegrator::DirectLighting(integrator) => integrator.get_sampler(),
            SamplerIntegrator::Path(integrator) => integrator.get_sampler(),
            SamplerIntegrator::VolPath(integrator) => integrator.get_sampler(),
            SamplerIntegrator::Whitted(integrator) => integrator.get_sampler(),
        }
    }
    pub fn get_pixel_bounds(&self) -> Bounds2i {
        match self {
            SamplerIntegrator::AO(integrator) => integrator.get_pixel_bounds(),
            SamplerIntegrator::DirectLighting(integrator) => integrator.get_pixel_bounds(),
            SamplerIntegrator::Path(integrator) => integrator.get_pixel_bounds(),
            SamplerIntegrator::VolPath(integrator) => integrator.get_pixel_bounds(),
            SamplerIntegrator::Whitted(integrator) => integrator.get_pixel_bounds(),
        }
    }
    pub fn specular_reflect(
        &self,
        ray: &Ray,
        isect: &SurfaceInteraction,
        scene: &Scene,
        sampler: &mut Sampler,
        depth: i32,
    ) -> Spectrum {
        match self {
            SamplerIntegrator::DirectLighting(integrator) => {
                integrator.specular_reflect(ray, isect, scene, sampler, depth)
            }
            SamplerIntegrator::Whitted(integrator) => {
                integrator.specular_reflect(ray, isect, scene, sampler, depth)
            }
            _ => Spectrum::default(),
        }
    }
    pub fn specular_transmit(
        &self,
        ray: &Ray,
        isect: &SurfaceInteraction,
        scene: &Scene,
        sampler: &mut Sampler,
        depth: i32,
    ) -> Spectrum {
        match self {
            SamplerIntegrator::DirectLighting(integrator) => {
                integrator.specular_transmit(ray, isect, scene, sampler, depth)
            }
            SamplerIntegrator::Whitted(integrator) => {
                integrator.specular_transmit(ray, isect, scene, sampler, depth)
            }
            _ => Spectrum::default(),
        }
    }
}

// see integrator.cpp

/// Most basic direct lighting strategy.
pub fn uniform_sample_all_lights(
    it: &SurfaceInteraction,
    scene: &Scene,
    sampler: &mut Sampler,
    n_light_samples: &[i32],
    handle_media: bool,
) -> Spectrum {
    // TODO: ProfilePhase p(Prof::DirectLighting);
    let mut l: Spectrum = Spectrum::new(0.0);
    for (j, n_samples) in n_light_samples.iter().enumerate().take(scene.lights.len()) {
        // accumulate contribution of _j_th light to _L_
        let light = &scene.lights[j];
        let (u_light_array_is_empty, u_light_array_idx, u_light_array_start) =
            sampler.get_2d_array_idxs(*n_samples);
        let (u_scattering_array_is_empty, u_scattering_array_idx, u_scattering_array_start) =
            sampler.get_2d_array_idxs(*n_samples);
        if u_light_array_is_empty || u_scattering_array_is_empty {
            // use a single sample for illumination from _light_
            let u_light: Point2f = sampler.get_2d();
            let u_scattering: Point2f = sampler.get_2d();
            l += estimate_direct(
                it,
                u_scattering,
                light.clone(),
                u_light,
                scene,
                sampler,
                handle_media,
                false,
            );
        } else {
            // estimate direct lighting using sample arrays
            let mut ld: Spectrum = Spectrum::new(0.0);
            for k in 0..*n_samples {
                let u_scattering_array_sample: Point2f = sampler.get_2d_sample(
                    u_scattering_array_idx,
                    u_scattering_array_start + k as usize,
                );
                let u_light_array_sample: Point2f =
                    sampler.get_2d_sample(u_light_array_idx, u_light_array_start + k as usize);
                ld += estimate_direct(
                    it,
                    u_scattering_array_sample,
                    light.clone(),
                    u_light_array_sample,
                    scene,
                    sampler,
                    handle_media,
                    false,
                );
            }
            l += ld / *n_samples as Float;
        }
    }
    l
}

/// Estimate direct lighting for only one randomly chosen light and
/// multiply the result by the number of lights to compensate.
pub fn uniform_sample_one_light(
    it: &dyn Interaction,
    scene: &Scene,
    sampler: &mut Sampler,
    handle_media: bool,
    light_distrib: Option<&Distribution1D>,
) -> Spectrum {
    // TODO: ProfilePhase p(Prof::DirectLighting);

    // randomly choose a single light to sample, _light_
    let n_lights: usize = scene.lights.len();
    if n_lights == 0_usize {
        return Spectrum::default();
    }
    let light_num: usize;
    let mut light_pdf: Option<Float> = Some(0.0 as Float);
    let pdf: Float;
    if let Some(light_distribution) = light_distrib {
        // if !light_distrib.is_null() {
        light_num = light_distribution.sample_discrete(sampler.get_1d(), light_pdf.as_mut());
        pdf = light_pdf.unwrap();
        if pdf == 0.0 as Float {
            return Spectrum::default();
        }
    } else {
        light_num = std::cmp::min(
            (sampler.get_1d() * n_lights as Float) as usize,
            n_lights - 1,
        );
        pdf = 1.0 as Float / n_lights as Float;
    }
    let light = &scene.lights[light_num];
    let u_light: Point2f = sampler.get_2d();
    let u_scattering: Point2f = sampler.get_2d();
    estimate_direct(
        it,
        u_scattering,
        light.clone(),
        u_light,
        scene,
        sampler,
        handle_media,
        false,
    ) / pdf
}

/// Computes a direct lighting estimate for a single light source sample.
pub fn estimate_direct(
    it: &dyn Interaction,
    u_scattering: Point2f,
    light: Arc<Light>,
    u_light: Point2f,
    scene: &Scene,
    sampler: &mut Sampler,
    // TODO: arena
    handle_media: bool,
    specular: bool,
) -> Spectrum {
    let bsdf_flags = if !specular {
        // bitwise not in Rust is ! (not the ~ operator like in C)
        BxdfType::BsdfAll as u8 & !(BxdfType::BsdfSpecular as u8)
    } else {
        BxdfType::BsdfAll as u8
    };
    let mut ld: Spectrum = Spectrum::new(0.0);
    // sample light source with multiple importance sampling
    let mut wi: Vector3f = Vector3f::default();
    let mut light_pdf: Float = 0.0 as Float;
    let mut scattering_pdf: Float = 0.0 as Float;
    let mut visibility: VisibilityTester = VisibilityTester::default();
    let it_common: InteractionCommon = InteractionCommon {
        p: it.get_p(),
        time: it.get_time(),
        p_error: it.get_p_error(),
        wo: it.get_wo(),
        n: it.get_n(),
        medium_interface: it.get_medium_interface(),
    };
    let mut li: Spectrum = light.sample_li(
        &it_common,
        u_light,
        &mut wi,
        &mut light_pdf,
        &mut visibility,
    );
    // TODO: println!("EstimateDirect uLight: {:?} -> Li: {:?}, wi:
    // {:?}, pdf: {:?}", u_light, li, wi, light_pdf);
    if light_pdf > 0.0 as Float && !li.is_black() {
        // compute BSDF or phase function's value for light sample
        let mut f: Spectrum = Spectrum::new(0.0);
        if it.is_surface_interaction() {
            // evaluate BSDF for light sampling strategy
            if let Some(ref bsdf) = it.get_bsdf() {
                if let Some(shading_n) = it.get_shading_n() {
                    f = bsdf.f(&it.get_wo(), &wi, bsdf_flags)
                        * Spectrum::new(vec3_abs_dot_nrm(&wi, &shading_n));
                    scattering_pdf = bsdf.pdf(&it.get_wo(), &wi, bsdf_flags);
                    // TODO: println!("  surf f*dot :{:?}, scatteringPdf: {:?}", f, scattering_pdf);
                }
            }
        } else {
            // evaluate phase function for light sampling strategy
            if let Some(ref phase) = it.get_phase() {
                let p: Float = phase.p(&it.get_wo(), &wi);
                f = Spectrum::new(p);
                scattering_pdf = p;
            }
        }
        if !f.is_black() {
            // compute effect of visibility for light source sample
            if handle_media {
                li *= visibility.tr(scene, sampler);
            } else if !visibility.unoccluded(scene) {
                li = Spectrum::new(0.0 as Float);
            }
            // add light's contribution to reflected radiance
            if !li.is_black() {
                if is_delta_light(light.get_flags()) {
                    ld += f * li / light_pdf;
                } else {
                    let weight: Float = power_heuristic(1_u8, light_pdf, 1_u8, scattering_pdf);
                    ld += f * li * Spectrum::new(weight) / light_pdf;
                }
            }
        }
    }
    // sample BSDF with multiple importance sampling
    if !is_delta_light(light.get_flags()) {
        let mut f: Spectrum = Spectrum::new(0.0);
        let mut sampled_specular: bool = false;
        if it.is_surface_interaction() {
            // sample scattered direction for surface interactions
            let mut sampled_type: u8 = 0_u8;
            if let Some(ref bsdf) = it.get_bsdf() {
                if let Some(shading_n) = it.get_shading_n() {
                    f = bsdf.sample_f(
                        &it.get_wo(),
                        &mut wi,
                        u_scattering,
                        &mut scattering_pdf,
                        bsdf_flags,
                        &mut sampled_type,
                    );
                    f *= Spectrum::new(vec3_abs_dot_nrm(&wi, &shading_n));
                    sampled_specular = (sampled_type & BxdfType::BsdfSpecular as u8) != 0_u8;
                }
            } else {
                println!("TODO: if let Some(ref bsdf) = it.get_bsdf() failed");
            }
        } else {
            // sample scattered direction for medium interactions
            if let Some(ref phase) = it.get_phase() {
                let p: Float = phase.sample_p(&it.get_wo(), &mut wi, u_scattering);
                f = Spectrum::new(p);
                scattering_pdf = p;
            }
        }
        // TODO: println!("  BSDF / phase sampling f: {:?}, scatteringPdf: {:?}",
        //          f, scattering_pdf);
        if !f.is_black() && scattering_pdf > 0.0 {
            // account for light contributions along sampled direction _wi_
            let weight = if !sampled_specular {
                light_pdf = light.pdf_li(it, wi);
                if light_pdf == 0.0 {
                    return ld;
                }
                power_heuristic(1, scattering_pdf, 1, light_pdf)
            } else {
                1.0
            };
            // find intersection and compute transmittance
            let mut ray: Ray = it.spawn_ray(&wi);
            let mut tr: Spectrum = Spectrum::new(1.0 as Float);
            let mut found_surface_interaction: bool = false;
            // add light contribution from material sampling
            let mut li: Spectrum = Spectrum::default();
            let mut light_isect: SurfaceInteraction = SurfaceInteraction::default();
            let mut tr_spectrum: Spectrum = Spectrum::default();
            if handle_media {
                let hit_surface: bool =
                    scene.intersect_tr(&mut ray, sampler, &mut light_isect, &mut tr_spectrum);
                tr = tr_spectrum; // copy return value
                if hit_surface {
                    found_surface_interaction = true;
                    if let Some(primitive_raw) = light_isect.primitive {
                        let primitive = unsafe { &*primitive_raw };
                        if let Some(area_light) = primitive.get_area_light() {
                            let pa = &*area_light as *const _ as *const usize;
                            let pl = &*light as *const _ as *const usize;
                            if pa == pl {
                                li = light_isect.le(&-wi);
                            }
                        }
                    }
                }
            } else if scene.intersect(&mut ray, &mut light_isect) {
                found_surface_interaction = true;
                if let Some(primitive_raw) = light_isect.primitive {
                    let primitive = unsafe { &*primitive_raw };
                    if let Some(area_light) = primitive.get_area_light() {
                        let pa = &*area_light as *const _ as *const usize;
                        let pl = &*light as *const _ as *const usize;
                        if pa == pl {
                            li = light_isect.le(&-wi);
                        }
                    }
                }
            }
            if !found_surface_interaction {
                li = light.le(&mut ray);
            }
            if !li.is_black() {
                ld += f * li * tr * weight / scattering_pdf;
            }
        }
    }
    ld
}

/// The light to start each photon path from is chosen according to a
/// PDF defined by the lights' respective powers.
pub fn compute_light_power_distribution(scene: &Scene) -> Option<Arc<Distribution1D>> {
    if scene.lights.is_empty() {
        return None;
    }
    let mut light_power: Vec<Float> = Vec::with_capacity(scene.lights.len());
    for li in 0..scene.lights.len() {
        let light = &scene.lights[li];
        light_power.push(light.power().y());
    }
    Some(Arc::new(Distribution1D::new(light_power)))
}
