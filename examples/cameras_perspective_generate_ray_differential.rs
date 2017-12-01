extern crate pbrt;

use pbrt::cameras::perspective::PerspectiveCamera;
use pbrt::core::camera::{Camera, CameraSample};
use pbrt::core::film::Film;
use pbrt::core::filter::Filter;
use pbrt::core::geometry::{Bounds2f, Point2f, Point2i, Point3f, Ray, Vector2f, Vector3f};
use pbrt::core::pbrt::Float;
use pbrt::core::transform::{AnimatedTransform, Transform};
use pbrt::filters::boxfilter::BoxFilter;
use std::sync::Arc;

fn main() {
    // CameraSample
    let p_film: Point2f = Point2f {
        x: 0.0453399941,
        y: 0.304356605,
    };
    let p_lens: Point2f = Point2f {
        x: 0.295576781,
        y: 0.392422706,
    };
    let time: Float = 0.41321516;
    let camera_sample: CameraSample = CameraSample {
        p_film: p_film,
        p_lens: p_lens,
        time: time,
    };
    println!("camera_sample = {:?}", camera_sample);
    // PerspectiveCamera
    let xw: Float = 0.5;
    let yw: Float = 0.5;
    let filter: Arc<Filter + Sync + Send> = Arc::new(BoxFilter {
                                                         radius: Vector2f { x: xw, y: yw },
                                                         inv_radius: Vector2f {
                                                             x: 1.0 / xw,
                                                             y: 1.0 / yw,
                                                         },
                                                     });
    let filename: String = String::from("spheres-differentials-texfilt.exr");
    let xres = 1000;
    let yres = 500;
    let crop: Bounds2f = Bounds2f {
        p_min: Point2f { x: 0.0, y: 0.0 },
        p_max: Point2f { x: 1.0, y: 1.0 },
    };
    let film: Arc<Film> = Arc::new(Film::new(Point2i { x: xres, y: yres },
                                             crop,
                                             filter,
                                             35.0,
                                             filename,
                                             1.0,
                                             std::f32::INFINITY));
    let pos = Point3f {
        x: 2.0,
        y: 2.0,
        z: 5.0,
    };
    let look = Point3f {
        x: 0.0,
        y: -0.4,
        z: 0.0,
    };
    let up = Vector3f {
        x: 0.0,
        y: 1.0,
        z: 0.0,
    };
    let t: Transform = Transform::look_at(pos, look, up);
    let it: Transform = Transform {
        m: t.m_inv.clone(),
        m_inv: t.m.clone(),
    };
    let animated_cam_to_world: AnimatedTransform = AnimatedTransform::new(&it, 0.0, &it, 1.0);
    let shutteropen: Float = 0.0;
    let shutterclose: Float = 1.0;
    let lensradius: Float = 0.0;
    let focaldistance: Float = 1e6;
    let frame: Float = xres as Float / yres as Float;
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
    let fov: Float = 30.0;
    let perspective_camera: PerspectiveCamera = PerspectiveCamera::new(animated_cam_to_world,
                                                                       screen,
                                                                       shutteropen,
                                                                       shutterclose,
                                                                       lensradius,
                                                                       focaldistance,
                                                                       fov,
                                                                       film);
    // println!("perspective_camera = {:?}", perspective_camera);
    let mut ray: Ray = Ray::default();
    let _ray_weight: Float = perspective_camera.generate_ray_differential(&camera_sample, &mut ray);
    println!("ray = {:?}", ray);
}
