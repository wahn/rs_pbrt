// std
use std;
use std::sync::Arc;
// pbrt
use core::camera::{Camera, CameraSample};
use core::film::Film;
use core::geometry::{Bounds2f, Point2f, Point2i, Point3f, Ray, RayDifferential, Vector3f};
use core::geometry::vec3_normalize;
use core::paramset::ParamSet;
use core::pbrt::Float;
use core::pbrt::lerp;
use core::sampling::concentric_sample_disk;
use core::transform::{AnimatedTransform, Transform};

// see perspective.h

pub struct PerspectiveCamera {
    // inherited from Camera (see camera.h)
    pub camera_to_world: AnimatedTransform,
    pub shutter_open: Float,
    pub shutter_close: Float,
    pub film: Arc<Film>,
    // TODO: const Medium *medium;
    // inherited from ProjectiveCamera (see camera.h)
    // camera_to_screen: Transform,
    pub raster_to_camera: Transform,
    // screen_to_raster: Transform,
    // raster_to_screen: Transform,
    pub lens_radius: Float,
    pub focal_distance: Float,
    // private data (see perspective.h)
    pub dx_camera: Vector3f,
    pub dy_camera: Vector3f,
    // a: Float,
}

impl PerspectiveCamera {
    pub fn new(
        camera_to_world: AnimatedTransform,
        screen_window: Bounds2f,
        shutter_open: Float,
        shutter_close: Float,
        lens_radius: Float,
        focal_distance: Float,
        fov: Float,
        film: Arc<Film>, /* const Medium *medium */
    ) -> Self {
        // see perspective.cpp
        let camera_to_screen: Transform = Transform::perspective(fov, 1e-2, 1000.0);
        // see camera.h
        // compute projective camera screen transformations
        let scale1 = Transform::scale(
            film.full_resolution.x as Float,
            film.full_resolution.y as Float,
            1.0,
        );
        let scale2 = Transform::scale(
            1.0 / (screen_window.p_max.x - screen_window.p_min.x),
            1.0 / (screen_window.p_min.y - screen_window.p_max.y),
            1.0,
        );
        let translate = Transform::translate(&Vector3f {
            x: -screen_window.p_min.x,
            y: -screen_window.p_max.y,
            z: 0.0,
        });
        let screen_to_raster = scale1 * scale2 * translate;
        let raster_to_screen = Transform::inverse(&screen_to_raster);
        let raster_to_camera = Transform::inverse(&camera_to_screen) * raster_to_screen;
        // see perspective.cpp
        // compute differential changes in origin for perspective camera rays
        let dx_camera: Vector3f = raster_to_camera.transform_point(&Point3f {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        }) - raster_to_camera.transform_point(&Point3f {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        });
        let dy_camera: Vector3f = raster_to_camera.transform_point(&Point3f {
            x: 0.0,
            y: 1.0,
            z: 0.0,
        }) - raster_to_camera.transform_point(&Point3f {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        });
        // compute image plane bounds at $z=1$ for _PerspectiveCamera_
        let res: Point2i = film.full_resolution;
        let mut p_min: Point3f = raster_to_camera.transform_point(&Point3f {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        });
        // Point3f p_max = RasterToCamera(Point3f(res.x, res.y, 0));
        let mut p_max: Point3f = raster_to_camera.transform_point(&Point3f {
            x: res.x as Float,
            y: res.y as Float,
            z: 0.0,
        });
        p_min /= p_min.z;
        p_max /= p_max.z;
        // let a: Float = ((p_max.x - p_min.x) * (p_max.y - p_min.y)).abs();

        PerspectiveCamera {
            camera_to_world: camera_to_world,
            shutter_open: shutter_open,
            shutter_close: shutter_close,
            film: film,
            // camera_to_screen: camera_to_screen,
            raster_to_camera: raster_to_camera,
            // screen_to_raster: screen_to_raster,
            // raster_to_screen: raster_to_screen,
            lens_radius: lens_radius,
            focal_distance: focal_distance,
            dx_camera: dx_camera,
            dy_camera: dy_camera,
            // a: a,
        }
    }
    pub fn create(
        params: &ParamSet,
        cam2world: AnimatedTransform,
        film: Arc<Film>, /* const Medium *medium */
    ) -> Box<Camera + Send + Sync> {
        let shutteropen: Float = params.find_one_float(String::from("shutteropen"), 0.0);
        let shutterclose: Float = params.find_one_float(String::from("shutterclose"), 1.0);
        // TODO: std::swap(shutterclose, shutteropen);
        assert!(shutterclose >= shutteropen);
        let lensradius: Float = params.find_one_float(String::from("lensradius"), 0.0);
        let focaldistance: Float = params.find_one_float(String::from("focaldistance"), 1e6);
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
        // TODO: const Float *sw = params.FindFloat("screenwindow", &swi);
        let fov: Float = params.find_one_float(String::from("fov"), 90.0);
        // let halffov: Float =
        //     params.find_one_float(String::from("halffov"), -1.0);
        // TODO: if (halffov > 0.f)
        // TODO: let perspective_camera: Arc<Camera + Sync + Send> =
        let camera = Box::new(PerspectiveCamera::new(
            cam2world,
            screen,
            shutteropen,
            shutterclose,
            lensradius,
            focaldistance,
            fov,
            film,
        ));
        camera
    }
}

impl Camera for PerspectiveCamera {
    fn generate_ray_differential(&self, sample: &CameraSample, ray: &mut Ray) -> Float {
        // TODO: ProfilePhase prof(Prof::GenerateCameraRay);
        // compute raster and camera sample positions
        let p_film: Point3f = Point3f {
            x: sample.p_film.x,
            y: sample.p_film.y,
            z: 0.0,
        };
        let p_camera: Point3f = self.raster_to_camera.transform_point(&p_film);
        let dir: Vector3f = vec3_normalize(&Vector3f {
            x: p_camera.x,
            y: p_camera.y,
            z: p_camera.z,
        });
        let mut diff: RayDifferential = RayDifferential {
            rx_origin: ray.o,
            ry_origin: ray.o,
            rx_direction: vec3_normalize(
                &(Vector3f {
                    x: p_camera.x,
                    y: p_camera.y,
                    z: p_camera.z,
                } + self.dx_camera),
            ),
            ry_direction: vec3_normalize(
                &(Vector3f {
                    x: p_camera.x,
                    y: p_camera.y,
                    z: p_camera.z,
                } + self.dy_camera),
            ),
        };
        // *ray = RayDifferential(Point3f(0, 0, 0), dir);
        let mut in_ray: Ray = Ray {
            o: Point3f::default(),
            d: dir,
            t_max: std::f32::INFINITY,
            time: lerp(sample.time, self.shutter_open, self.shutter_close),
            differential: Some(diff),
        };
        // modify ray for depth of field
        if self.lens_radius > 0.0 as Float {
            // sample point on lens
            let p_lens: Point2f = concentric_sample_disk(sample.p_lens) * self.lens_radius;
            // compute point on plane of focus
            let ft: Float = self.focal_distance / in_ray.d.z;
            let p_focus: Point3f = in_ray.position(ft);
            // update ray for effect of lens
            in_ray.o = Point3f {
                x: p_lens.x,
                y: p_lens.y,
                z: 0.0 as Float,
            };
            in_ray.d = vec3_normalize(&(p_focus - in_ray.o));
        }
        // compute offset rays for _PerspectiveCamera_ ray differentials
        if self.lens_radius > 0.0 as Float {
            // compute _PerspectiveCamera_ ray differentials accounting for lens

            // sample point on lens
            let p_lens: Point2f = concentric_sample_disk(sample.p_lens) * self.lens_radius;
            let dx: Vector3f = vec3_normalize(&Vector3f::from(p_camera + self.dx_camera));
            let ft: Float = self.focal_distance / dx.z;
            let p_focus: Point3f = Point3f::default() + (dx * ft);
            diff.rx_origin = Point3f {
                x: p_lens.x,
                y: p_lens.y,
                z: 0.0 as Float,
            };
            diff.rx_direction = vec3_normalize(&(p_focus - diff.rx_origin));
            let dy: Vector3f = vec3_normalize(&Vector3f::from(p_camera + self.dy_camera));
            let ft: Float = self.focal_distance / dy.z;
            let p_focus: Point3f = Point3f::default() + (dy * ft);
            diff.ry_origin = Point3f {
                x: p_lens.x,
                y: p_lens.y,
                z: 0.0 as Float,
            };
            diff.ry_direction = vec3_normalize(&(p_focus - diff.ry_origin));
            // replace differential
            in_ray.differential = Some(diff);
        }
        // TODO: ray->medium = medium;
        *ray = self.camera_to_world.transform_ray(&in_ray);
        1.0
    }
    fn get_film(&self) -> Arc<Film> {
        self.film.clone()
    }
}
