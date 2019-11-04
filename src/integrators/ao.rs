// std
use std::sync::Arc;
// pbrt
use crate::core::camera::{Camera, CameraSample};
use crate::core::geometry::{nrm_cross_vec3, nrm_faceforward_vec3, vec3_dot_nrm};
use crate::core::geometry::{Bounds2i, Normal3f, Point2f, Point2i, Ray, Vector2i, Vector3f};
// use crate::core::integrator::SamplerIntegrator;
use crate::blockqueue::BlockQueue;
use crate::core::geometry::pnt2_inside_exclusive;
use crate::core::interaction::Interaction;
use crate::core::material::TransportMode;
use crate::core::pbrt::{Float, Spectrum};
use crate::core::sampler::Sampler;
use crate::core::sampling::{
    cosine_hemisphere_pdf, cosine_sample_hemisphere, uniform_hemisphere_pdf,
    uniform_sample_hemisphere,
};
use crate::core::scene::Scene;

// see ao.h

/// Ambient Occlusion
pub struct AOIntegrator {
    // inherited from SamplerIntegrator (see integrator.h)
    pub camera: Arc<Camera>,
    pub sampler: Box<dyn Sampler + Send + Sync>,
    pub pixel_bounds: Bounds2i,
    // see ao.h
    pub cos_sample: bool,
    pub n_samples: i32,
}

impl AOIntegrator {
    pub fn new(
        cos_sample: bool,
        n_samples: i32,
        camera: Arc<Camera>,
        sampler: Box<dyn Sampler + Send + Sync>,
        pixel_bounds: Bounds2i,
    ) -> Self {
        AOIntegrator {
            camera,
            sampler,
            pixel_bounds,
            cos_sample,
            n_samples,
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
    pub fn preprocess(&mut self, _scene: &Scene) {
        self.sampler.request_2d_array(self.n_samples);
    }
    pub fn li(
        &self,
        r: &mut Ray,
        scene: &Scene,
        sampler: &mut Box<dyn Sampler + Send + Sync>,
        // arena: &mut Arena,
        _depth: i32,
    ) -> Spectrum {
        // TODO: ProfilePhase p(Prof::SamplerIntegratorLi);
        let mut l: Spectrum = Spectrum::default();
        let mut ray: Ray = Ray {
            o: r.o,
            d: r.d,
            t_max: r.t_max,
            time: r.time,
            differential: r.differential,
            medium: r.medium.clone(),
        };
        if let Some(mut isect) = scene.intersect(&mut ray) {
            let mode: TransportMode = TransportMode::Radiance;
            isect.compute_scattering_functions(&mut ray, true, mode);
            // if (!isect.bsdf) {
            //     VLOG(2) << "Skipping intersection due to null bsdf";
            //     ray = isect.SpawnRay(ray.d);
            //     goto retry;
            // }
            // compute coordinate frame based on true geometry, not
            // shading geometry.
            let n: Normal3f = nrm_faceforward_vec3(&isect.n, &-ray.d);
            let s: Vector3f = isect.dpdu.normalize();
            let t: Vector3f = nrm_cross_vec3(&isect.n, &s);
            let u: Vec<Point2f> = sampler.get_2d_array(self.n_samples);
            for i in 0..self.n_samples as usize {
                // Vector3f wi;
                let mut wi: Vector3f;
                let pdf: Float;
                if self.cos_sample {
                    wi = cosine_sample_hemisphere(&u[i]);
                    pdf = cosine_hemisphere_pdf(wi.z.abs());
                } else {
                    wi = uniform_sample_hemisphere(&u[i]);
                    pdf = uniform_hemisphere_pdf();
                }
                // transform wi from local frame to world space.
                wi = Vector3f {
                    x: s.x * wi.x + t.x * wi.y + n.x * wi.z,
                    y: s.y * wi.x + t.y * wi.y + n.y * wi.z,
                    z: s.z * wi.x + t.z * wi.y + n.z * wi.z,
                };
                let mut ray: Ray = isect.spawn_ray(&wi);
                if !scene.intersect_p(&mut ray) {
                    l += Spectrum::new(vec3_dot_nrm(&wi, &n) / (pdf * self.n_samples as Float));
                }
            }
        }
        l
    }
    pub fn get_pixel_bounds(&self) -> Bounds2i {
        self.pixel_bounds
    }
}
