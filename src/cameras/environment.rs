// std
use std;
use std::f32::consts::PI;
use std::sync::Arc;
// pbrt
use core::camera::{Camera, CameraSample};
use core::film::Film;
use core::geometry::{Bounds2f, Bounds2i, Normal3f, Point2f, Point2i, Point3f, Ray,
                     RayDifferential, Vector3f};
use core::geometry::{nrm_abs_dot_vec3, vec3_dot_vec3, vec3_normalize};
use core::interaction::InteractionCommon;
use core::light::VisibilityTester;
use core::medium::{Medium, MediumInterface};
use core::paramset::ParamSet;
use core::pbrt::lerp;
use core::pbrt::{Float, Spectrum};
use core::sampling::concentric_sample_disk;
use core::transform::{AnimatedTransform, Transform};

// see environment.h

pub struct EnvironmentCamera {
    // inherited from Camera (see camera.h)
    pub camera_to_world: AnimatedTransform,
    pub shutter_open: Float,
    pub shutter_close: Float,
    pub film: Arc<Film>,
    pub medium: Option<Arc<Medium + Send + Sync>>,
}

impl EnvironmentCamera {
    pub fn new(
        camera_to_world: AnimatedTransform,
        shutter_open: Float,
        shutter_close: Float,
        film: Arc<Film>,
        medium: Option<Arc<Medium + Send + Sync>>,
    ) -> Self {
        EnvironmentCamera {
            camera_to_world: camera_to_world,
            shutter_open: shutter_open,
            shutter_close: shutter_close,
            film: film,
            medium: medium,
        }
    }
    pub fn create(
        params: &ParamSet,
        cam2world: AnimatedTransform,
        film: Arc<Film>,
        medium: Option<Arc<Medium + Send + Sync>>,
    ) -> Arc<Camera + Send + Sync> {
        let shutteropen: Float = params.find_one_float(String::from("shutteropen"), 0.0);
        let shutterclose: Float = params.find_one_float(String::from("shutterclose"), 1.0);
        // TODO: std::swap(shutterclose, shutteropen);
        assert!(shutterclose >= shutteropen);
        let lensradius: Float = params.find_one_float(String::from("lensradius"), 0.0);
        let focaldistance: Float = params.find_one_float(String::from("focaldistance"), 1e30);
        let frame: Float = params.find_one_float(
            String::from("frameaspectratio"),
            (film.full_resolution.x as Float) / (film.full_resolution.y as Float),
        );
        let mut screen: Bounds2f = Bounds2f::default();
        if frame > 1.0 {
            screen.p_min.x = -frame;
            screen.p_max.x = frame;
            screen.p_min.y = -1.0;
            screen.p_max.y = 1.0;
        } else {
            screen.p_min.x = -1.0;
            screen.p_max.x = 1.0;
            screen.p_min.y = -1.0 / frame;
            screen.p_max.y = 1.0 / frame;
        }
        let sw: Vec<Float> = params.find_float(String::from("screenwindow"));
        if sw.len() > 0_usize {
            if sw.len() == 4 {
                screen.p_min.x = sw[0];
                screen.p_max.x = sw[1];
                screen.p_min.y = sw[2];
                screen.p_max.y = sw[3];
            } else {
                panic!("\"screenwindow\" should have four values");
            }
        }
        let camera = Arc::new(EnvironmentCamera::new(
            cam2world,
            shutteropen,
            shutterclose,
            film,
            medium,
        ));
        camera
    }
}


impl Camera for EnvironmentCamera {
    fn generate_ray_differential(&self, sample: &CameraSample, ray: &mut Ray) -> Float {
        // WORK
        0.0 as Float
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
