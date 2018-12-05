// std
use std;
use std::path::PathBuf;
use std::sync::Arc;
// pbrt
use core::camera::{Camera, CameraSample};
use core::film::Film;
use core::floatfile::read_float_file;
use core::geometry::{nrm_faceforward_vec3, nrm_normalize, vec3_normalize};
use core::geometry::{Bounds2f, Normal3f, Point2f, Point3f, Ray, Vector3f};
use core::interaction::InteractionCommon;
use core::light::VisibilityTester;
use core::medium::Medium;
use core::paramset::ParamSet;
use core::pbrt::quadratic;
use core::pbrt::{Float, Spectrum};
use core::reflection::refract;
use core::transform::{AnimatedTransform, Transform};

// see realistic.h

#[derive(Debug, Default, Copy, Clone)]
pub struct LensElementInterface {
    pub curvature_radius: Float,
    pub thickness: Float,
    pub eta: Float,
    pub aperture_radius: Float,
}

pub struct RealisticCamera {
    // inherited from Camera (see camera.h)
    pub camera_to_world: AnimatedTransform,
    pub shutter_open: Float,
    pub shutter_close: Float,
    pub film: Arc<Film>,
    pub medium: Option<Arc<Medium + Send + Sync>>,
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
        lens_data: &Vec<Float>,
        film: Arc<Film>,
        medium: Option<Arc<Medium + Send + Sync>>,
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
            println!("{:?}", element_interfaces[i / 4]);
        }
        let camera = RealisticCamera {
            camera_to_world: camera_to_world,
            shutter_open: shutter_open,
            shutter_close: shutter_close,
            film: film,
            medium: medium,
            simple_weighting: simple_weighting,
            element_interfaces: element_interfaces,
            exit_pupil_bounds: Vec::new(),
        };
        // compute lens--film distance for given focus distance
        camera.focus_binary_search(focus_distance);
        // WORK
        camera
    }
    pub fn create(
        params: &ParamSet,
        cam2world: AnimatedTransform,
        film: Arc<Film>,
        medium: Option<Arc<Medium + Send + Sync>>,
        search_directory: Option<&Box<PathBuf>>,
    ) -> Arc<Camera + Send + Sync> {
        let shutteropen: Float = params.find_one_float("shutteropen", 0.0);
        let shutterclose: Float = params.find_one_float("shutterclose", 1.0);
        // TODO: std::swap(shutterclose, shutteropen);
        assert!(shutterclose >= shutteropen);
        // realistic camera-specific parameters
        let mut lens_file: String = params.find_one_filename("lensfile", String::from(""));
        if lens_file != String::from("") {
            if let Some(ref search_directory) = search_directory {
                let mut path_buf: PathBuf = PathBuf::from("/");
                path_buf.push(search_directory.as_ref());
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
        println!("lens_data = {:?}", lens_data);
        let camera = Arc::new(RealisticCamera::new(
            cam2world,
            shutteropen,
            shutterclose,
            aperture_diameter,
            focus_distance,
            simple_weighting,
            &lens_data,
            film,
            medium,
        ));
        camera
    }
    pub fn lens_rear_z(&self) -> Float {
        // WORK
        0.0
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
        // WORK
        0.0
    }
    pub fn trace_lenses_from_film(&self, r_camera: &Ray, r_out: &mut Ray) -> bool {
        // WORK
        false
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
        let o: Point3f = ray.o - Vector3f {
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
        *n = nrm_faceforward_vec3(&nrm_normalize(&*n), &-ray.d);
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
                if self.intersect_spherical_element(radius, z_center, &r_lens, &mut t, &mut n) {
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
                let eta_i: Float;
                if i == 0 || self.element_interfaces[i - 1].eta == 0.0 as Float {
                    eta_i = 0.0 as Float;
                } else {
                    eta_i = self.element_interfaces[i - 1].eta;
                }
                let eta_t: Float;
                if self.element_interfaces[i].eta != 0.0 as Float {
                    eta_t = self.element_interfaces[i].eta;
                } else {
                    eta_t = 1.0 as Float;
                }
                if !refract(&vec3_normalize(&-r_lens.d), &n, eta_i / eta_t, &mut wt) {
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
        // WORK
    }
    pub fn draw_ray_path_from_film(&self, r: &Ray, arrow: bool, to_optical_intercept: bool) {
        // WORK
    }
    pub fn draw_ray_path_from_scene(&self, r: &Ray, arrow: bool, to_optical_intercept: bool) {
        // WORK
    }
    pub fn compute_cardinal_points(&self, r_in: &Ray, r_out: &Ray, p: &mut Float, f: &mut Float) {
        // WORK
    }
    pub fn compute_thick_lens_approximation(&self, pz: &mut [Float; 2], f: &mut [Float; 2]) {
        // find height $x$ from optical axis for parallel rays
        let x: Float = 0.001 as Float * self.film.diagonal;
        // compute cardinal points for film side of lens system
        let r_scene: Ray = Ray {
            o: Point3f {
                x: x,
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
        // ComputeCardinalPoints(r_scene, r_film, &pz[0], &fz[0]);

        // // Compute cardinal points for scene side of lens system
        // r_film = Ray(Point3f(x, 0, LensRearZ() - 1), Vector3f(0, 0, 1));
        // CHECK(TraceLensesFromFilm(r_film, &r_scene))
        //     << "Unable to trace ray from film to scene for thick lens "
        //        "approximation. Is aperture stop extremely small?";
        // ComputeCardinalPoints(r_film, r_scene, &pz[1], &fz[1]);
        // WORK
    }
    pub fn focus_thick_lens(&self, focus_distance: Float) -> Float {
        // WORK
        0.0
    }
    pub fn focus_binary_search(&self, focus_distance: Float) -> Float {
        // Float pz[2], fz[2];
        let mut pz: [Float; 2] = [0.0 as Float; 2];
        let mut fz: [Float; 2] = [0.0 as Float; 2];
        // ComputeThickLensApproximation(pz, fz);
        self.compute_thick_lens_approximation(&mut pz, &mut fz);
        // LOG(INFO) << StringPrintf("Cardinal points: p' = %f f' = %f, p = %f f = %f.\n",
        //                           pz[0], fz[0], pz[1], fz[1]);
        // LOG(INFO) << StringPrintf("Effective focal length %f\n", fz[0] - pz[0]);
        // // Compute translation of lens, _delta_, to focus at _focusDistance_
        // Float f = fz[0] - pz[0];
        // Float z = -focusDistance;
        // Float c = (pz[1] - z - pz[0]) * (pz[1] - z - 4 * f - pz[0]);
        // CHECK_GT(c, 0) << "Coefficient must be positive. It looks focusDistance: " << focusDistance << " is too short for a given lenses configuration";
        // Float delta =
        //     0.5f * (pz[1] - z + pz[0] - std::sqrt(c));
        // return elementInterfaces.back().thickness + delta;
        // WORK
        0.0
    }
    pub fn focus_distance(&self, film_dist: Float) -> Float {
        // WORK
        0.0
    }
    pub fn bound_exit_pupil(&self, p_film_x0: Float, p_film_x1: Float) -> Bounds2f {
        // WORK
        Bounds2f::default()
    }
    pub fn render_exit_pupil(&self, sx: Float, sy: Float, filename: String) {
        // WORK
    }
    pub fn sample_exit_pupil(
        &self,
        p_film: &Point2f,
        lens_sample: &Point2f,
        sample_bounds_area: &mut Float,
    ) -> Point3f {
        // WORK
        Point3f::default()
    }
    pub fn test_exit_pupil_bounds(&self) {
        // WORK
    }
}

impl Camera for RealisticCamera {
    fn generate_ray_differential(&self, sample: &CameraSample, ray: &mut Ray) -> Float {
        // WORK
        0.0
    }
    fn we(&self, _ray: &Ray, _p_raster2: Option<&mut Point2f>) -> Spectrum {
        panic!("camera::we() is not implemented!");
        // Spectrum::default()
    }
    fn pdf_we(&self, _ray: &Ray) -> (Float, Float) {
        // let mut pdf_pos: Float = 0.0;
        // let mut pdf_dir: Float = 0.0;
        panic!("camera::pdf_we() is not implemented!");
        // (pdf_pos, pdf_dir)
    }
    fn sample_wi(
        &self,
        _iref: &InteractionCommon,
        _u: &Point2f,
        _wi: &mut Vector3f,
        _pdf: &mut Float,
        _p_raster: &mut Point2f,
        _vis: &mut VisibilityTester,
    ) -> Spectrum {
        panic!("camera::sample_wi() is not implemented!");
        // Spectrum::default()
    }
    fn get_film(&self) -> Arc<Film> {
        self.film.clone()
    }
}
