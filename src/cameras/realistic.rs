// std
use std;
use std::f32::consts::PI;
use std::path::PathBuf;
use std::sync::Arc;
// pbrt
use core::camera::{Camera, CameraSample};
use core::film::Film;
use core::geometry::{nrm_abs_dot_vec3, vec3_dot_vec3, vec3_normalize};
use core::geometry::{
    Bounds2f, Bounds2i, Normal3f, Point2f, Point2i, Point3f, Ray, RayDifferential, Vector3f,
};
use core::interaction::InteractionCommon;
use core::light::VisibilityTester;
use core::medium::{Medium, MediumInterface};
use core::paramset::ParamSet;
use core::pbrt::lerp;
use core::pbrt::{Float, Spectrum};
use core::sampling::concentric_sample_disk;
use core::transform::{AnimatedTransform, Transform};

// see realistic.h

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
        RealisticCamera {
            camera_to_world: camera_to_world,
            shutter_open: shutter_open,
            shutter_close: shutter_close,
            film: film,
            medium: medium,
            // TODO
            simple_weighting: false,
            element_interfaces: Vec::new(),
            exit_pupil_bounds: Vec::new(),
        }
    }
    pub fn create(
        params: &ParamSet,
        cam2world: AnimatedTransform,
        film: Arc<Film>,
        medium: Option<Arc<Medium + Send + Sync>>,
        search_directory: Option<&Box<PathBuf>>,
    ) -> Arc<Camera + Send + Sync> {
        let shutteropen: Float = params.find_one_float(String::from("shutteropen"), 0.0);
        let shutterclose: Float = params.find_one_float(String::from("shutterclose"), 1.0);
        // TODO: std::swap(shutterclose, shutteropen);
        assert!(shutterclose >= shutteropen);
        // realistic camera-specific parameters
        let mut lens_file: String =
            params.find_one_filename(String::from("lensfile"), String::from(""));
        if lens_file != String::from("") {
            if let Some(ref search_directory) = search_directory {
                let mut path_buf: PathBuf = PathBuf::from("/");
                path_buf.push(search_directory.as_ref());
                path_buf.push(lens_file);
                lens_file = String::from(path_buf.to_str().unwrap());
            }
        }
        println!("lens_file = {:?}", lens_file);
        // WORK
        let camera = Arc::new(RealisticCamera::new(
            cam2world,
            shutteropen,
            shutterclose,
            // TODO
            0.0,         // aperture_diameter
            0.0,         // focus_distance
            false,       // simple_weighting
            &Vec::new(), // lens_data
            // TODO
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
        // WORK
        z_sum
    }
    pub fn rear_element_radius(&self) -> Float {
        // WORK
        0.0
    }
    pub fn trace_lenses_from_film(&self, ray: &Ray, r_out: &mut Ray) -> bool {
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
        // WORK
        false
    }
    pub fn trace_lenses_from_scene(&self, r_camera: &Ray, r_out: &mut Ray) -> bool {
        // WORK
        false
    }
    pub fn draw_lens_system(&self) {
        // WORK
    }
    pub fn draw_ray_path_from_film(r: &Ray, arrow: bool, to_optical_intercept: bool) {
        // WORK
    }
    pub fn draw_ray_path_from_scene(r: &Ray, arrow: bool, to_optical_intercept: bool) {
        // WORK
    }
    pub fn compute_cardinal_points(r_in: &Ray, r_out: &Ray, p: &mut Float, f: &mut Float) {
        // WORK
    }
    pub fn compute_thick_lens_approximation(pz: [Float; 2], f: [Float; 2]) {
        // WORK
    }
    pub fn focus_thick_lens(focus_distance: Float) -> Float {
        // WORK
        0.0
    }
    pub fn focus_binary_search(focus_distance: Float) -> Float {
        // WORK
        0.0
    }
    pub fn focus_distance(film_dist: Float) -> Float {
        // WORK
        0.0
    }
    pub fn bound_exit_pupil(p_film_x0: Float, p_film_x1: Float) -> Bounds2f {
        // WORK
        Bounds2f::default()
    }
    pub fn render_exit_pupil(sx: Float, sy: Float, filename: String) {
        // WORK
    }
    pub fn sample_exit_pupil(
        p_film: &Point2f,
        lens_sample: &Point2f,
        sample_bounds_area: &mut Float,
    ) -> Point3f {
        // WORK
        Point3f::default()
    }
    pub fn test_exit_pupil_bounds() {
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
