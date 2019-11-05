// std
use std::sync::Arc;
// pbrt
use crate::blockqueue::BlockQueue;
use crate::core::camera::{Camera, CameraSample};
use crate::core::geometry::{pnt2_inside_exclusive, vec3_abs_dot_nrm, vec3_dot_nrm};
use crate::core::geometry::{
    Bounds2i, Normal3f, Point2i, Ray, RayDifferential, Vector2i, Vector3f,
};
use crate::core::integrator::{uniform_sample_all_lights, uniform_sample_one_light};
use crate::core::interaction::{Interaction, SurfaceInteraction};
use crate::core::material::TransportMode;
use crate::core::pbrt::{Float, Spectrum};
use crate::core::reflection::BxdfType;
use crate::core::sampler::Sampler;
use crate::core::scene::Scene;

// see directlighting.h

#[derive(Debug, Clone, PartialEq)]
pub enum LightStrategy {
    UniformSampleAll,
    UniformSampleOne,
}

/// Direct Lighting (no Global Illumination)
pub struct DirectLightingIntegrator {
    // inherited from SamplerIntegrator (see integrator.h)
    pub camera: Arc<Camera>,
    pub sampler: Box<dyn Sampler + Send + Sync>,
    pixel_bounds: Bounds2i,
    // see directlighting.h
    strategy: LightStrategy,
    max_depth: u32,
    n_light_samples: Vec<i32>,
}

impl DirectLightingIntegrator {
    pub fn new(
        strategy: LightStrategy,
        max_depth: u32,
        camera: Arc<Camera>,
        sampler: Box<dyn Sampler + Send + Sync>,
        pixel_bounds: Bounds2i,
    ) -> Self {
        DirectLightingIntegrator {
            camera,
            sampler,
            pixel_bounds,
            strategy,
            max_depth,
            n_light_samples: Vec::new(),
        }
    }
    pub fn render(&mut self, scene: &Scene, num_threads: u8) {
        let film = self.camera.get_film();
        let sample_bounds: Bounds2i = film.get_sample_bounds();
        self.preprocess(scene);
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
            let integrator = &self;
            let bq = &block_queue;
            let sampler = &self.sampler;
            let camera = &self.camera;
            let film = &film;
            let pixel_bounds = self.get_pixel_bounds().clone();
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
                                            / (tile_sampler.get_samples_per_pixel() as Float)
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
    pub fn preprocess(&mut self, scene: &Scene) {
        if self.strategy == LightStrategy::UniformSampleAll {
            // compute number of samples to use for each light
            for li in 0..scene.lights.len() {
                let ref light = scene.lights[li];
                self.n_light_samples
                    .push(self.sampler.round_count(light.get_n_samples()));
            }
            // request samples for sampling all lights
            for _i in 0..self.max_depth {
                for j in 0..scene.lights.len() {
                    self.sampler.request_2d_array(self.n_light_samples[j]);
                    self.sampler.request_2d_array(self.n_light_samples[j]);
                }
            }
        }
    }
    pub fn li(
        &self,
        ray: &mut Ray,
        scene: &Scene,
        sampler: &mut Box<dyn Sampler + Send + Sync>,
        // arena: &mut Arena,
        depth: i32,
    ) -> Spectrum {
        // TODO: ProfilePhase p(Prof::SamplerIntegratorLi);
        let mut l: Spectrum = Spectrum::new(0.0 as Float);
        // find closest ray intersection or return background radiance
        if let Some(mut isect) = scene.intersect(ray) {
            // compute scattering functions for surface interaction
            let mode: TransportMode = TransportMode::Radiance;
            isect.compute_scattering_functions(ray /* arena, */, false, mode);
            // if (!isect.bsdf)
            //     return Li(isect.SpawnRay(ray.d), scene, sampler, arena, depth);
            let wo: Vector3f = isect.wo;
            l += isect.le(&wo);
            if scene.lights.len() > 0 {
                // compute direct lighting for _DirectLightingIntegrator_ integrator
                if self.strategy == LightStrategy::UniformSampleAll {
                    l += uniform_sample_all_lights(
                        &isect,
                        scene,
                        sampler,
                        &self.n_light_samples,
                        false,
                    );
                } else {
                    l += uniform_sample_one_light(&isect, scene, sampler, false, None);
                }
            }
            if ((depth + 1_i32) as u32) < self.max_depth {
                // trace rays for specular reflection and refraction
                l += self.specular_reflect(
                    ray, &isect, scene, sampler, // arena,
                    depth,
                );
                l += self.specular_transmit(
                    ray, &isect, scene, sampler, // arena,
                    depth,
                );
            }
        } else {
            for light in &scene.lights {
                l += light.le(ray);
            }
        }
        l
    }
    pub fn get_pixel_bounds(&self) -> Bounds2i {
        self.pixel_bounds
    }
    pub fn specular_reflect(
        &self,
        ray: &Ray,
        isect: &SurfaceInteraction,
        scene: &Scene,
        sampler: &mut Box<dyn Sampler + Send + Sync>,
        // arena: &mut Arena,
        depth: i32,
    ) -> Spectrum {
        // compute specular reflection direction _wi_ and BSDF value
        let wo: Vector3f = isect.wo;
        let mut wi: Vector3f = Vector3f::default();
        let mut pdf: Float = 0.0 as Float;
        let ns: Normal3f = isect.shading.n;
        let mut sampled_type: u8 = 0_u8;
        let bsdf_flags: u8 = BxdfType::BsdfReflection as u8 | BxdfType::BsdfSpecular as u8;
        let f: Spectrum;
        if let Some(ref bsdf) = isect.bsdf {
            f = bsdf.sample_f(
                &wo,
                &mut wi,
                &sampler.get_2d(),
                &mut pdf,
                bsdf_flags,
                &mut sampled_type,
            );
            if pdf > 0.0 as Float && !f.is_black() && vec3_abs_dot_nrm(&wi, &ns) != 0.0 as Float {
                // compute ray differential _rd_ for specular reflection
                let mut rd: Ray = isect.spawn_ray(&wi);
                if let Some(d) = ray.differential.iter().next() {
                    let dudx: Float = *isect.dudx.read().unwrap();
                    let dvdx: Float = *isect.dvdx.read().unwrap();
                    let dndx: Normal3f = isect.shading.dndu * dudx + isect.shading.dndv * dvdx;
                    let dudy: Float = *isect.dudy.read().unwrap();
                    let dvdy: Float = *isect.dvdy.read().unwrap();
                    let dndy: Normal3f = isect.shading.dndu * dudy + isect.shading.dndv * dvdy;
                    let dwodx: Vector3f = -d.rx_direction - wo;
                    let dwody: Vector3f = -d.ry_direction - wo;
                    let ddndx: Float = vec3_dot_nrm(&dwodx, &ns) + vec3_dot_nrm(&wo, &dndx);
                    let ddndy: Float = vec3_dot_nrm(&dwody, &ns) + vec3_dot_nrm(&wo, &dndy);
                    // compute differential reflected directions
                    let dpdx: Vector3f = *isect.dpdx.read().unwrap();
                    let dpdy: Vector3f = *isect.dpdy.read().unwrap();
                    let diff: RayDifferential = RayDifferential {
                        rx_origin: isect.p + dpdx,
                        ry_origin: isect.p + dpdy,
                        rx_direction: wi - dwodx
                            + Vector3f::from(dndx * vec3_dot_nrm(&wo, &ns) + ns * ddndx)
                                * 2.0 as Float,
                        ry_direction: wi - dwody
                            + Vector3f::from(dndy * vec3_dot_nrm(&wo, &ns) + ns * ddndy)
                                * 2.0 as Float,
                    };
                    rd.differential = Some(diff);
                }
                return f
                    * self.li(&mut rd, scene, sampler, depth + 1)
                    * Spectrum::new(vec3_abs_dot_nrm(&wi, &ns) / pdf);
            } else {
                Spectrum::new(0.0)
            }
        } else {
            Spectrum::new(0.0)
        }
    }
    pub fn specular_transmit(
        &self,
        ray: &Ray,
        isect: &SurfaceInteraction,
        scene: &Scene,
        sampler: &mut Box<dyn Sampler + Send + Sync>,
        // arena: &mut Arena,
        depth: i32,
    ) -> Spectrum {
        let wo: Vector3f = isect.wo;
        let mut wi: Vector3f = Vector3f::default();
        let mut pdf: Float = 0.0 as Float;
        // let p: Point3f = isect.p;
        let ns: Normal3f = isect.shading.n;
        let mut sampled_type: u8 = 0_u8;
        let bsdf_flags: u8 = BxdfType::BsdfTransmission as u8 | BxdfType::BsdfSpecular as u8;
        let f: Spectrum;
        if let Some(ref bsdf) = isect.bsdf {
            f = bsdf.sample_f(
                &wo,
                &mut wi,
                &sampler.get_2d(),
                &mut pdf,
                bsdf_flags,
                &mut sampled_type,
            );
            if pdf > 0.0 as Float && !f.is_black() && vec3_abs_dot_nrm(&wi, &ns) != 0.0 as Float {
                // compute ray differential _rd_ for specular transmission
                let mut rd: Ray = isect.spawn_ray(&wi);
                if let Some(d) = ray.differential.iter().next() {
                    let mut eta: Float = bsdf.eta;
                    let w: Vector3f = -wo;
                    if vec3_dot_nrm(&wo, &ns) < 0.0 as Float {
                        eta = 1.0 / eta;
                    }
                    let dudx: Float = *isect.dudx.read().unwrap();
                    let dvdx: Float = *isect.dvdx.read().unwrap();
                    let dndx: Normal3f = isect.shading.dndu * dudx + isect.shading.dndv * dvdx;
                    let dudy: Float = *isect.dudy.read().unwrap();
                    let dvdy: Float = *isect.dvdy.read().unwrap();
                    let dndy: Normal3f = isect.shading.dndu * dudy + isect.shading.dndv * dvdy;
                    let dwodx: Vector3f = -d.rx_direction - wo;
                    let dwody: Vector3f = -d.ry_direction - wo;
                    let ddndx: Float = vec3_dot_nrm(&dwodx, &ns) + vec3_dot_nrm(&wo, &dndx);
                    let ddndy: Float = vec3_dot_nrm(&dwody, &ns) + vec3_dot_nrm(&wo, &dndy);
                    let mu: Float = eta * vec3_dot_nrm(&w, &ns) - vec3_dot_nrm(&wi, &ns);
                    let dmudx: Float = (eta
                        - (eta * eta * vec3_dot_nrm(&w, &ns)) / vec3_dot_nrm(&wi, &ns))
                        * ddndx;
                    let dmudy: Float = (eta
                        - (eta * eta * vec3_dot_nrm(&w, &ns)) / vec3_dot_nrm(&wi, &ns))
                        * ddndy;
                    let dpdx: Vector3f = *isect.dpdx.read().unwrap();
                    let dpdy: Vector3f = *isect.dpdy.read().unwrap();
                    let diff: RayDifferential = RayDifferential {
                        rx_origin: isect.p + dpdx,
                        ry_origin: isect.p + dpdy,
                        rx_direction: wi + dwodx * eta - Vector3f::from(dndx * mu + ns * dmudx),
                        ry_direction: wi + dwody * eta - Vector3f::from(dndy * mu + ns * dmudy),
                    };
                    rd.differential = Some(diff);
                }
                return f
                    * self.li(&mut rd, scene, sampler, depth + 1)
                    * Spectrum::new(vec3_abs_dot_nrm(&wi, &ns) / pdf);
            } else {
                Spectrum::new(0.0)
            }
        } else {
            Spectrum::new(0.0)
        }
    }
}
