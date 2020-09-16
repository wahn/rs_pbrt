// std
use std::cell::Cell;
use std::f32::consts::PI;
use std::sync::Arc;
// pbrt
use crate::core::camera::{Camera, CameraSample};
use crate::core::film::Film;
use crate::core::geometry::{Bounds2f, Point2f, Point3f, Ray, Vector3f};
use crate::core::interaction::InteractionCommon;
use crate::core::light::VisibilityTester;
use crate::core::medium::Medium;
use crate::core::paramset::ParamSet;
use crate::core::pbrt::lerp;
use crate::core::pbrt::{Float, Spectrum};
use crate::core::transform::AnimatedTransform;

// see environment.h

pub struct EnvironmentCamera {
    // inherited from Camera (see camera.h)
    pub camera_to_world: AnimatedTransform,
    pub shutter_open: Float,
    pub shutter_close: Float,
    pub film: Arc<Film>,
    pub medium: Option<Arc<Medium>>,
}

impl EnvironmentCamera {
    pub fn new(
        camera_to_world: AnimatedTransform,
        shutter_open: Float,
        shutter_close: Float,
        film: Arc<Film>,
        medium: Option<Arc<Medium>>,
    ) -> Self {
        EnvironmentCamera {
            camera_to_world,
            shutter_open,
            shutter_close,
            film,
            medium,
        }
    }
    pub fn create(
        params: &ParamSet,
        cam2world: AnimatedTransform,
        film: Arc<Film>,
        medium: Option<Arc<Medium>>,
    ) -> Arc<Camera> {
        let shutteropen: Float = params.find_one_float("shutteropen", 0.0);
        let shutterclose: Float = params.find_one_float("shutterclose", 1.0);
        // TODO: std::swap(shutterclose, shutteropen);
        assert!(shutterclose >= shutteropen);
        // let lensradius: Float = params.find_one_float(String::from("lensradius"), 0.0);
        // let focaldistance: Float = params.find_one_float(String::from("focaldistance"), 1e30);
        let frame: Float = params.find_one_float(
            "frameaspectratio",
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
        let sw: Vec<Float> = params.find_float("screenwindow");
        if !sw.is_empty() {
            if sw.len() == 4 {
                screen.p_min.x = sw[0];
                screen.p_max.x = sw[1];
                screen.p_min.y = sw[2];
                screen.p_max.y = sw[3];
            } else {
                panic!("\"screenwindow\" should have four values");
            }
        }
        Arc::new(Camera::Environment(Box::new(EnvironmentCamera::new(
            cam2world,
            shutteropen,
            shutterclose,
            film,
            medium,
        ))))
    }
    // Camera
    pub fn generate_ray_differential(&self, sample: &CameraSample, ray: &mut Ray) -> Float {
        let theta: Float = PI * sample.p_film.y / self.film.full_resolution.y as Float;
        let phi: Float = 2.0 as Float * PI * sample.p_film.x / self.film.full_resolution.x as Float;
        let dir: Vector3f = Vector3f {
            x: theta.sin() * phi.cos(),
            y: theta.cos(),
            z: theta.sin() * phi.sin(),
        };
        let mut in_ray: Ray = Ray {
            o: Point3f::default(),
            d: dir,
            t_max: Cell::new(std::f32::INFINITY),
            time: lerp(sample.time, self.shutter_open, self.shutter_close),
            medium: None,
            differential: None,
        };
        // ray->medium = medium;
        if let Some(ref medium_arc) = self.medium {
            in_ray.medium = Some(medium_arc.clone());
        } else {
            in_ray.medium = None;
        }
        *ray = self.camera_to_world.transform_ray(&in_ray);
        1.0
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
        _lens_intr: &mut InteractionCommon,
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
