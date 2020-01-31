// std
use std;
use std::path::PathBuf;
use std::sync::Arc;
// pbrt
use crate::core::camera::{Camera, CameraSample};
use crate::core::film::Film;
use crate::core::floatfile::read_float_file;
use crate::core::geometry::{bnd2_expand, bnd2_union_pnt2, nrm_faceforward_vec3, pnt2_inside_bnd2};
use crate::core::geometry::{Bounds2f, Normal3f, Point2f, Point3f, Ray, RayDifferential, Vector3f};
use crate::core::interaction::InteractionCommon;
use crate::core::light::VisibilityTester;
use crate::core::lowdiscrepancy::radical_inverse;
use crate::core::medium::Medium;
use crate::core::paramset::ParamSet;
use crate::core::pbrt::{lerp, quadratic};
use crate::core::pbrt::{Float, Spectrum};
use crate::core::reflection::refract;
use crate::core::transform::{AnimatedTransform, Transform};

// see realistic.h

#[derive(Debug, Default, Copy, Clone)]
pub struct LensElementInterface {
    pub curvature_radius: Float,
    pub thickness: Float,
    pub eta: Float,
    pub aperture_radius: Float,
}

#[derive(Clone)]
pub struct RealisticCamera {
    // inherited from Camera (see camera.h)
    pub camera_to_world: AnimatedTransform,
    pub shutter_open: Float,
    pub shutter_close: Float,
    pub film: Arc<Film>,
    pub medium: Option<Arc<Medium>>,
    // private data (see realistic.h)
    pub simple_weighting: bool,
    pub element_interfaces: Vec<LensElementInterface>,
    pub exit_pupil_bounds: Vec<Bounds2f>,
}

impl RealisticCamera {
    pub fn new(
        camera_to_world: AnimatedTransform,
        shutter_open: Float,
        shutter_close: Float,
        aperture_diameter: Float,
        focus_distance: Float,
        simple_weighting: bool,
        lens_data: &[Float],
        film: Arc<Film>,
        medium: Option<Arc<Medium>>,
    ) -> Self {
        let mut element_interfaces: Vec<LensElementInterface> = Vec::new();
        for i in (0..lens_data.len()).step_by(4) {
            let mut diameter: Float = lens_data[i + 3];
            if lens_data[i] == 0.0 as Float {
                if aperture_diameter > lens_data[i + 3] {
                    println!("Specified aperture diameter {} is greater than maximum possible {}.  Clamping it.",
                             aperture_diameter,
                             lens_data[i + 3]);
                } else {
                    diameter = aperture_diameter;
                }
            }
            element_interfaces.push(LensElementInterface {
                curvature_radius: lens_data[i] * 0.001 as Float,
                thickness: lens_data[i + 1] * 0.001 as Float,
                eta: lens_data[i + 2],
                aperture_radius: diameter * 0.001 as Float / 2.0 as Float,
            });
            // println!("{:?}", element_interfaces[i / 4]);
        }
        let mut camera = RealisticCamera {
            camera_to_world,
            shutter_open,
            shutter_close,
            film: film.clone(),
            medium,
            simple_weighting,
            element_interfaces,
            exit_pupil_bounds: Vec::new(),
        };
        // compute lens--film distance for given focus distance
        let _fb: Float = camera.focus_binary_search(focus_distance);
        // LOG(INFO) << StringPrintf("Binary search focus: %f -> %f\n", fb,
        //                           camera.focus_distance(fb));
        camera.element_interfaces.last_mut().unwrap().thickness =
            camera.focus_thick_lens(focus_distance);
        // LOG(INFO) << StringPrintf("Thick lens focus: %f -> %f\n",
        //                           camera.element_interfaces.last().unwrap().thickness,
        //                           camera.focus_distance(elementInterfaces.back().thickness));
        // compute exit pupil bounds at sampled points on the film
        let n_samples: usize = 64;
        let mut exit_pupil_bounds: Vec<Bounds2f> = Vec::new();
        exit_pupil_bounds.resize(n_samples, Bounds2f::default());
        let num_cores: usize = num_cpus::get();
        let chunk_size: usize = n_samples / num_cores;
        {
            let bands: Vec<&mut [Bounds2f]> = exit_pupil_bounds.chunks_mut(chunk_size).collect();
            let camera = &camera;
            let film = &film;
            crossbeam::scope(|scope| {
                let (band_tx, band_rx) = crossbeam_channel::bounded(num_cores);
                // spawn worker threads
                for (b, band) in bands.into_iter().enumerate() {
                    let band_tx = band_tx.clone();
                    scope.spawn(move |_| {
                        for (index, bound) in band.iter_mut().enumerate() {
                            let i: usize = (b * chunk_size) + index;
                            let r0: Float =
                                i as Float / n_samples as Float * film.diagonal / 2.0 as Float;
                            let r1: Float = (i + 1) as Float / n_samples as Float * film.diagonal
                                / 2.0 as Float;
                            *bound = camera.bound_exit_pupil(r0, r1);
                        }
                    });
                    // send progress through the channel to main thread
                    band_tx
                        .send(b)
                        .unwrap_or_else(|_| panic!("Failed to send progress"));
                }
                // spawn thread to report progress
                scope.spawn(move |_| {
                    for _ in pbr::PbIter::new(0..num_cores) {
                        band_rx.recv().unwrap();
                    }
                });
            })
            .unwrap();
        }
        camera.exit_pupil_bounds = exit_pupil_bounds;
        if camera.simple_weighting {
            println!("WARNING: \"simpleweighting\" option with RealisticCamera no longer necessarily matches regular camera images. Further, pixel values will vary a bit depending on the aperture size. See this discussion for details: https://github.com/mmp/pbrt-v3/issues/162#issuecomment-348625837");
        }
        camera
    }
    pub fn create(
        params: &ParamSet,
        cam2world: AnimatedTransform,
        film: Arc<Film>,
        medium: Option<Arc<Medium>>,
        search_directory: Option<&PathBuf>,
    ) -> Arc<Camera> {
        let shutteropen: Float = params.find_one_float("shutteropen", 0.0);
        let shutterclose: Float = params.find_one_float("shutterclose", 1.0);
        // TODO: std::swap(shutterclose, shutteropen);
        assert!(shutterclose >= shutteropen);
        // realistic camera-specific parameters
        let mut lens_file: String = params.find_one_filename("lensfile", String::from(""));
        if lens_file != "" {
            if let Some(ref search_directory) = search_directory {
                let mut path_buf: PathBuf = PathBuf::from("/");
                path_buf.push(search_directory);
                path_buf.push(lens_file);
                lens_file = String::from(path_buf.to_str().unwrap());
            }
        }
        if lens_file == "" {
            println!("ERROR: No lens description file supplied!");
        } else {
            println!("lens_file = {:?}", lens_file);
        }
        let aperture_diameter: Float = params.find_one_float("aperturediameter", 1.0);
        let focus_distance: Float = params.find_one_float("focusdistance", 10.0);
        let simple_weighting: bool = params.find_one_bool("simpleweighting", true);
        let mut lens_data: Vec<Float> = Vec::new();
        if !read_float_file(&lens_file, &mut lens_data) {
            println!(
                "ERROR: Error reading lens specification file {:?}.",
                lens_file
            );
        }
        if lens_data.len() % 4_usize != 0_usize {
            println!("ERROR: Excess values in lens specification file {:?}; must be multiple-of-four values, read {}.",
                     lens_file, lens_data.len());
        }
        // println!("lens_data = {:?}", lens_data);
        Arc::new(Camera::Realistic(Box::new(RealisticCamera::new(
            cam2world,
            shutteropen,
            shutterclose,
            aperture_diameter,
            focus_distance,
            simple_weighting,
            &lens_data,
            film,
            medium,
        ))))
    }
    pub fn generate_ray(&self, sample: &CameraSample, ray: &mut Ray) -> Float {
        // TODO: ProfilePhase prof(Prof::GenerateCameraRay);
        // ++totalRays;
        // find point on film, _p_film_, corresponding to _sample.p_film_
        let s: Point2f = Point2f {
            x: sample.p_film.x / self.film.full_resolution.x as Float,
            y: sample.p_film.y / self.film.full_resolution.y as Float,
        };
        let p_film2: Point2f = self.film.get_physical_extent().lerp(s);
        let p_film: Point3f = Point3f {
            x: -p_film2.x,
            y: p_film2.y,
            z: 0.0 as Float,
        };
        // trace ray from _p_film_ through lens system
        let mut exit_pupil_bounds_area: Float = 0.0 as Float;
        let p_rear: Point3f = self.sample_exit_pupil(
            Point2f {
                x: p_film.x,
                y: p_film.y,
            },
            sample.p_lens,
            &mut exit_pupil_bounds_area,
        );
        let mut r_film: Ray = Ray::default();
        r_film.o = p_film;
        r_film.d = p_rear - p_film;
        r_film.t_max = std::f32::INFINITY;
        r_film.time = lerp(sample.time, self.shutter_open, self.shutter_close);
        if !self.trace_lenses_from_film(&r_film, Some(ray)) {
            // ++vignettedRays;
            return 0.0 as Float;
        }
        // finish initialization of _RealisticCamera_ ray
        *ray = self.camera_to_world.transform_ray(&ray);
        ray.d = ray.d.normalize();
        if let Some(ref medium_arc) = self.medium {
            ray.medium = Some(medium_arc.clone());
        } else {
            ray.medium = None;
        }
        // return weighting for _RealisticCamera_ ray
        let cos_theta: Float = r_film.d.normalize().z;
        let cos_2_theta: Float = cos_theta * cos_theta;
        let cos_4_theta: Float = cos_2_theta * cos_2_theta;
        if self.simple_weighting {
            cos_4_theta * exit_pupil_bounds_area / self.exit_pupil_bounds[0].area()
        } else {
            (self.shutter_close - self.shutter_open) * (cos_4_theta * exit_pupil_bounds_area)
                / (self.lens_rear_z() * self.lens_rear_z())
        }
    }
    pub fn lens_rear_z(&self) -> Float {
        self.element_interfaces.last().unwrap().thickness
    }
    pub fn lens_front_z(&self) -> Float {
        let mut z_sum = 0.0;
        for i in 0..self.element_interfaces.len() {
            let element = self.element_interfaces[i];
            z_sum += element.thickness
        }
        z_sum
    }
    pub fn rear_element_radius(&self) -> Float {
        self.element_interfaces.last().unwrap().aperture_radius
    }
    pub fn trace_lenses_from_film(&self, r_camera: &Ray, r_out: Option<&mut Ray>) -> bool {
        let mut element_z: Float = 0.0 as Float;
        // transform _rCamera_ from camera to lens system space
        let camera_to_lens: Transform = Transform::scale(1.0 as Float, 1.0 as Float, -1.0 as Float);
        let mut r_lens: Ray = camera_to_lens.transform_ray(r_camera);
        let ei_len = self.element_interfaces.len();
        for idx in 0..ei_len {
            let i = ei_len - 1 - idx;
            let element = self.element_interfaces[i];
            // update ray from film accounting for interaction with _element_
            element_z -= element.thickness;
            // compute intersection of ray with lens element
            let mut t: Float = 0.0 as Float;
            let mut n: Normal3f = Normal3f::default();
            let is_stop: bool = element.curvature_radius == 0.0 as Float;
            if is_stop {
                // The refracted ray computed in the previous lens
                // element interface may be pointed towards film
                // plane(+z) in some extreme situations; in such
                // cases, 't' becomes negative.
                if r_lens.d.z >= 0.0 as Float {
                    return false;
                }
                t = (element_z - r_lens.o.z) / r_lens.d.z;
            } else {
                let radius: Float = element.curvature_radius;
                let z_center: Float = element_z + element.curvature_radius;
                if !self.intersect_spherical_element(radius, z_center, &r_lens, &mut t, &mut n) {
                    return false;
                }
            }
            assert!(t >= 0.0 as Float);
            // test intersection point against element aperture
            let p_hit: Point3f = r_lens.position(t);
            let r2: Float = p_hit.x * p_hit.x + p_hit.y * p_hit.y;
            if r2 > element.aperture_radius * element.aperture_radius {
                return false;
            }
            r_lens.o = p_hit;
            // update ray path for element interface interaction
            if !is_stop {
                let mut w: Vector3f = Vector3f::default();
                let eta_i: Float = element.eta;
                let eta_t = if i > 0_usize && self.element_interfaces[i - 1].eta != 0.0 as Float {
                    self.element_interfaces[i - 1].eta
                } else {
                    1.0 as Float
                };
                if !refract(&(-r_lens.d).normalize(), &n, eta_i / eta_t, &mut w) {
                    return false;
                }
                r_lens.d = w;
            }
        }
        // transform _r_lens_ from lens system space back to camera space
        if let Some(r_out) = r_out {
            let lens_to_camera: Transform =
                Transform::scale(1.0 as Float, 1.0 as Float, -1.0 as Float);
            *r_out = lens_to_camera.transform_ray(&r_lens);
        }
        true
    }
    pub fn intersect_spherical_element(
        &self,
        radius: Float,
        z_center: Float,
        ray: &Ray,
        t: &mut Float,
        n: &mut Normal3f,
    ) -> bool {
        // compute _t0_ and _t1_ for ray--element intersection
        let o: Point3f = ray.o
            - Vector3f {
                x: 0.0 as Float,
                y: 0.0 as Float,
                z: z_center,
            };
        let a: Float = ray.d.x * ray.d.x + ray.d.y * ray.d.y + ray.d.z * ray.d.z;
        let b: Float = 2.0 as Float * (ray.d.x * o.x + ray.d.y * o.y + ray.d.z * o.z);
        let c: Float = o.x * o.x + o.y * o.y + o.z * o.z - radius * radius;
        let mut t0: Float = 0.0 as Float;
        let mut t1: Float = 0.0 as Float;
        if !quadratic(a, b, c, &mut t0, &mut t1) {
            return false;
        }
        // select intersection $t$ based on ray direction and element curvature
        let use_closer_t: bool = (ray.d.z > 0.0 as Float) ^ (radius < 0.0 as Float);
        if use_closer_t {
            *t = t0.min(t1);
        } else {
            *t = t0.max(t1);
        }
        if *t < 0.0 as Float {
            return false;
        }
        // compute surface normal of element at ray intersection point
        *n = Normal3f::from(Vector3f::from(o + ray.d * *t));
        *n = nrm_faceforward_vec3(&n.normalize(), &-ray.d);
        true
    }
    pub fn trace_lenses_from_scene(&self, r_camera: &Ray, r_out: Option<&mut Ray>) -> bool {
        let mut element_z: Float = -self.lens_front_z();
        // transform _r_camera_ from camera to lens system space
        let camera_to_lens: Transform = Transform::scale(1.0 as Float, 1.0 as Float, -1.0 as Float);
        let mut r_lens: Ray = camera_to_lens.transform_ray(r_camera);
        for i in 0..self.element_interfaces.len() {
            let element = self.element_interfaces[i];
            // compute intersection of ray with lens element
            let mut t: Float = 0.0 as Float;
            let mut n: Normal3f = Normal3f::default();
            let is_stop: bool = element.curvature_radius == 0.0 as Float;
            if is_stop {
                t = (element_z - r_lens.o.z) / r_lens.d.z;
            } else {
                let radius: Float = element.curvature_radius;
                let z_center: Float = element_z + element.curvature_radius;
                if !self.intersect_spherical_element(radius, z_center, &r_lens, &mut t, &mut n) {
                    return false;
                }
            }
            assert!(t >= 0.0 as Float);
            // test intersection point against element aperture
            let p_hit: Point3f = r_lens.position(t);
            let r2: Float = p_hit.x * p_hit.x + p_hit.y * p_hit.y;
            if r2 > element.aperture_radius * element.aperture_radius {
                return false;
            }
            r_lens.o = p_hit;
            // update ray path for from-scene element interface interaction
            if !is_stop {
                let mut wt: Vector3f = Vector3f::default();
                let eta_i = if i == 0 || self.element_interfaces[i - 1].eta == 0.0 as Float {
                    1.0 as Float
                } else {
                    self.element_interfaces[i - 1].eta
                };
                let eta_t = if self.element_interfaces[i].eta != 0.0 as Float {
                    self.element_interfaces[i].eta
                } else {
                    1.0 as Float
                };
                if !refract(&(-r_lens.d).normalize(), &n, eta_i / eta_t, &mut wt) {
                    return false;
                }
                r_lens.d = wt;
            }
            element_z += element.thickness;
        }
        // transform _r_lens_ from lens system space back to camera space
        if let Some(r_out) = r_out {
            let lens_to_camera: Transform =
                Transform::scale(1.0 as Float, 1.0 as Float, -1.0 as Float);
            *r_out = lens_to_camera.transform_ray(&r_lens);
        }
        true
    }
    pub fn draw_lens_system(&self) {
        println!("TODO: RealisticCamera::draw_lens_system()");
    }
    pub fn draw_ray_path_from_film(&self, _r: &Ray, _arrow: bool, _to_optical_intercept: bool) {
        println!("TODO: RealisticCamera::draw_ray_path_from_film()");
    }
    pub fn draw_ray_path_from_scene(&self, _r: &Ray, _arrow: bool, _to_optical_intercept: bool) {
        println!("TODO: RealisticCamera::draw_ray_path_from_scene()");
    }
    pub fn compute_cardinal_points(
        &self,
        r_in: &Ray,
        r_out: &Ray,
        idx: usize,
        pz: &mut [Float; 2],
        fz: &mut [Float; 2],
    ) {
        let tf: Float = -r_out.o.x / r_out.d.x;
        fz[idx] = -r_out.position(tf).z;
        let tp: Float = (r_in.o.x - r_out.o.x) / r_out.d.x;
        pz[idx] = -r_out.position(tp).z;
    }
    pub fn compute_thick_lens_approximation(&self, pz: &mut [Float; 2], fz: &mut [Float; 2]) {
        // find height $x$ from optical axis for parallel rays
        let x: Float = 0.001 as Float * self.film.diagonal;
        // compute cardinal points for film side of lens system
        let mut r_scene: Ray = Ray {
            o: Point3f {
                x,
                y: 0.0 as Float,
                z: self.lens_front_z() + 1.0 as Float,
            },
            d: Vector3f {
                x: 0.0 as Float,
                y: 0.0 as Float,
                z: -1.0 as Float,
            },
            t_max: std::f32::INFINITY,
            time: 0.0 as Float,
            medium: None,
            differential: None,
        };
        let mut r_film: Ray = Ray::default();
        assert!(self.trace_lenses_from_scene(&r_scene, Some(&mut r_film)),
                "Unable to trace ray from scene to film for thick lens approximation. Is aperture stop extremely small?");
        self.compute_cardinal_points(&r_scene, &r_film, 0, pz, fz);
        // compute cardinal points for scene side of lens system
        r_film.o = Point3f {
            x,
            y: 0.0 as Float,
            z: self.lens_rear_z() - 1.0 as Float,
        };
        r_film.d = Vector3f {
            x: 0.0 as Float,
            y: 0.0 as Float,
            z: 1.0 as Float,
        };
        assert!(self.trace_lenses_from_film(&r_film, Some(&mut r_scene)),
                "Unable to trace ray from film to scene for thick lens approximation. Is aperture stop extremely small?");
        self.compute_cardinal_points(&r_film, &r_scene, 1, pz, fz);
    }
    pub fn focus_thick_lens(&self, focus_distance: Float) -> Float {
        let mut pz: [Float; 2] = [0.0 as Float; 2];
        let mut fz: [Float; 2] = [0.0 as Float; 2];
        self.compute_thick_lens_approximation(&mut pz, &mut fz);
        // LOG(INFO) << StringPrintf("Cardinal points: p' = %f f' = %f, p = %f f = %f.\n",
        //                           pz[0], fz[0], pz[1], fz[1]);
        // LOG(INFO) << StringPrintf("Effective focal length %f\n", fz[0] - pz[0]);
        // compute translation of lens, _delta_, to focus at _focus_distance_
        let f: Float = fz[0] - pz[0];
        let z: Float = -focus_distance;
        let c: Float = (pz[1] - z - pz[0]) * (pz[1] - z - 4.0 as Float * f - pz[0]);
        assert!(c > 0.0 as Float,
                "Coefficient must be positive. It looks focus_distance: {} is too short for a given lenses configuration",
                focus_distance);
        let delta: Float = 0.5 as Float * (pz[1] - z + pz[0] - c.sqrt());
        self.element_interfaces.last().unwrap().thickness + delta
    }
    pub fn focus_binary_search(&self, focus_distance: Float) -> Float {
        // find _film_distance_lower_, _film_distance_upper_ that bound focus distance
        let mut film_distance_upper: Float = self.focus_thick_lens(focus_distance);
        let mut film_distance_lower: Float = film_distance_upper;
        while self.focus_distance(film_distance_lower) > focus_distance {
            film_distance_lower *= 1.005 as Float;
        }
        while self.focus_distance(film_distance_upper) < focus_distance {
            film_distance_upper /= 1.005 as Float;
        }
        // do binary search on film distances to focus
        for _i in 0..20 {
            let fmid: Float = 0.5 as Float * (film_distance_lower + film_distance_upper);
            let mid_focus: Float = self.focus_distance(fmid);
            if mid_focus < focus_distance {
                film_distance_lower = fmid;
            } else {
                film_distance_upper = fmid;
            }
        }
        0.5 as Float * (film_distance_lower + film_distance_upper)
    }
    pub fn focus_distance(&self, film_distance: Float) -> Float {
        // find offset ray from film center through lens
        let bounds: Bounds2f =
            self.bound_exit_pupil(0.0 as Float, 0.001 as Float * self.film.diagonal);
        let scale_factors: [Float; 3] = [0.1 as Float, 0.01 as Float, 0.001 as Float];
        let mut lu: Float = 0.0;
        let mut ray: Ray = Ray::default();
        // Try some different and decreasing scaling factor to find
        // focus ray more quickly when `aperturediameter` is too
        // small.  (e.g. 2 [mm] for `aperturediameter` with
        // wide.22mm.dat),
        let mut found_focus_ray: bool = false;
        for scale in scale_factors.iter() {
            lu = scale * bounds.p_max[0];
            if self.trace_lenses_from_film(
                &Ray {
                    o: Point3f {
                        x: 0.0 as Float,
                        y: 0.0 as Float,
                        z: self.lens_rear_z() - film_distance,
                    },
                    d: Vector3f {
                        x: lu,
                        y: 0.0 as Float,
                        z: film_distance,
                    },
                    t_max: std::f32::INFINITY,
                    time: 0.0 as Float,
                    medium: None,
                    differential: None,
                },
                Some(&mut ray),
            ) {
                found_focus_ray = true;
                break;
            }
        }
        if !found_focus_ray {
            println!(
                "ERROR: Focus ray at lens pos({},0) didn't make it through the lenses with film distance {}?!??",
                lu, film_distance);
            return std::f32::INFINITY;
        }
        // compute distance _zFocus_ where ray intersects the principal axis
        let t_focus: Float = -ray.o.x / ray.d.x;
        let mut z_focus: Float = ray.position(t_focus).z;
        if z_focus < 0.0 as Float {
            z_focus = std::f32::INFINITY;
        }
        z_focus
    }
    pub fn bound_exit_pupil(&self, p_film_x0: Float, p_film_x1: Float) -> Bounds2f {
        let mut pupil_bounds: Bounds2f = Bounds2f::default();
        // sample a collection of points on the rear lens to find exit pupil
        let n_samples: i32 = 1024 * 1024;
        let mut n_exiting_rays: i32 = 0;
        // compute bounding box of projection of rear element on sampling plane
        let rear_radius: Float = self.rear_element_radius();
        let proj_rear_bounds: Bounds2f = Bounds2f {
            p_min: Point2f {
                x: -1.5 as Float * rear_radius,
                y: -1.5 as Float * rear_radius,
            },
            p_max: Point2f {
                x: 1.5 as Float * rear_radius,
                y: 1.5 as Float * rear_radius,
            },
        };
        for i in 0..n_samples {
            // find location of sample points on $x$ segment and rear lens element
            let p_film: Point3f = Point3f {
                x: lerp(
                    (i as Float + 0.5 as Float) / n_samples as Float,
                    p_film_x0,
                    p_film_x1,
                ),
                y: 0.0 as Float,
                z: 0.0 as Float,
            };
            let u: [Float; 2] = [
                radical_inverse(0 as u16, i as u64),
                radical_inverse(1 as u16, i as u64),
            ];
            let p_rear: Point3f = Point3f {
                x: lerp(u[0], proj_rear_bounds.p_min.x, proj_rear_bounds.p_max.x),
                y: lerp(u[1], proj_rear_bounds.p_min.y, proj_rear_bounds.p_max.y),
                z: self.lens_rear_z(),
            };
            // expand pupil bounds if ray makes it through the lens system
            if pnt2_inside_bnd2(
                Point2f {
                    x: p_rear.x,
                    y: p_rear.y,
                },
                &pupil_bounds,
            ) || self.trace_lenses_from_film(
                &Ray {
                    o: p_film,
                    d: p_rear - p_film,
                    t_max: std::f32::INFINITY,
                    time: 0.0 as Float,
                    medium: None,
                    differential: None,
                },
                None,
            ) {
                pupil_bounds = bnd2_union_pnt2(
                    &pupil_bounds,
                    Point2f {
                        x: p_rear.x,
                        y: p_rear.y,
                    },
                );
                n_exiting_rays += 1;
            }
        }
        // return entire element bounds if no rays made it through the lens system
        if n_exiting_rays == 0_i32 {
            // println!(
            //     "Unable to find exit pupil in x = [{},{}] on film.",
            //     p_film_x0, p_film_x1
            // );
            return proj_rear_bounds;
        }
        // expand bounds to account for sample spacing
        pupil_bounds = bnd2_expand(
            &pupil_bounds,
            2.0 as Float * proj_rear_bounds.diagonal().length() / (n_samples as Float).sqrt(),
        );
        pupil_bounds
    }
    pub fn render_exit_pupil(&self, _sx: Float, _sy: Float, _filename: String) {
        println!("TODO: RealisticCamera::render_exit_pupil()");
    }
    pub fn sample_exit_pupil(
        &self,
        p_film: Point2f,
        lens_sample: Point2f,
        sample_bounds_area: &mut Float,
    ) -> Point3f {
        // find exit pupil bound for sample distance from film center
        let r_film: Float = (p_film.x * p_film.x + p_film.y * p_film.y).sqrt();
        let mut r_index: usize = (r_film / (self.film.diagonal / 2.0 as Float)
            * self.exit_pupil_bounds.len() as Float)
            .floor() as usize;
        r_index = (self.exit_pupil_bounds.len() - 1).min(r_index);
        let pupil_bounds: Bounds2f = self.exit_pupil_bounds[r_index];
        *sample_bounds_area = pupil_bounds.area();
        // generate sample point inside exit pupil bound
        let p_lens: Point2f = pupil_bounds.lerp(lens_sample);
        // return sample point rotated by angle of _p_film_ with $+x$ axis
        let sin_theta = if r_film != 0.0 as Float {
            p_film.y / r_film
        } else {
            0.0 as Float
        };
        let cos_theta = if r_film != 0.0 as Float {
            p_film.x / r_film
        } else {
            1.0 as Float
        };
        Point3f {
            x: cos_theta * p_lens.x - sin_theta * p_lens.y,
            y: sin_theta * p_lens.x + cos_theta * p_lens.y,
            z: self.lens_rear_z(),
        }
    }
    pub fn test_exit_pupil_bounds(&self) {
        println!("TODO: RealisticCamera::test_exit_pupil_bounds()");
    }
    // Camera
    pub fn generate_ray_differential(&self, sample: &CameraSample, ray: &mut Ray) -> Float {
        let wt: Float = self.generate_ray(sample, ray);
        if wt == 0.0 as Float {
            return 0.0 as Float;
        }
        let mut rd = RayDifferential::default();
        // find camera ray after shifting a fraction of a pixel in the $x$ direction
        let mut wtx: Float = 0.0 as Float;
        let eps_values: [Float; 2] = [0.05 as Float, -0.05 as Float];
        for eps in eps_values.iter() {
            let mut sshift: CameraSample = *sample;
            sshift.p_film.x += eps;
            let mut rx: Ray = Ray::default();
            wtx = self.generate_ray(&sshift, &mut rx);
            rd.rx_origin = ray.o + (rx.o - ray.o) / *eps;
            rd.rx_direction = ray.d + (rx.d - ray.d) / *eps;
            if wtx != 0.0 as Float {
                break;
            }
        }
        if wtx == 0.0 as Float {
            return 0.0 as Float;
        }
        // find camera ray after shifting a fraction of a pixel in the $y$ direction
        let mut wty: Float = 0.0 as Float;
        for eps in eps_values.iter() {
            let mut sshift: CameraSample = *sample;
            sshift.p_film.y += eps;
            let mut ry: Ray = Ray::default();
            wty = self.generate_ray(&sshift, &mut ry);
            rd.ry_origin = ray.o + (ry.o - ray.o) / *eps;
            rd.ry_direction = ray.d + (ry.d - ray.d) / *eps;
            if wty != 0.0 as Float {
                break;
            }
        }
        if wty == 0.0 as Float {
            return 0.0 as Float;
        }
        // rd.has_differentials = true;
        ray.differential = Some(rd);
        wt
    }
    pub fn we(&self, _ray: &Ray, _p_raster2: Option<&mut Point2f>) -> Spectrum {
        panic!("camera::we() is not implemented!");
        // Spectrum::default()
    }
    pub fn pdf_we(&self, _ray: &Ray) -> (Float, Float) {
        // let mut pdf_pos: Float = 0.0;
        // let mut pdf_dir: Float = 0.0;
        panic!("camera::pdf_we() is not implemented!");
        // (pdf_pos, pdf_dir)
    }
    pub fn sample_wi(
        &self,
        _iref: &InteractionCommon,
        _u: Point2f,
        _wi: &mut Vector3f,
        _pdf: &mut Float,
        _p_raster: &mut Point2f,
        _vis: &mut VisibilityTester,
    ) -> Spectrum {
        panic!("camera::sample_wi() is not implemented!");
        // Spectrum::default()
    }
    pub fn get_shutter_open(&self) -> Float {
        self.shutter_open
    }
    pub fn get_shutter_close(&self) -> Float {
        self.shutter_close
    }
    pub fn get_film(&self) -> Arc<Film> {
        self.film.clone()
    }
}
