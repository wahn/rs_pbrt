extern crate pbrt;

use pbrt::{AnimatedTransform, Bounds2f, Bounds2i, BoxFilter, Film, Float, PerspectiveCamera,
           Point2i, Point2f, Point3f, Transform, Vector2f, Vector3f};

fn main() {
    // see  perspective.cpp CreatePerspectiveCamera()

    let shutteropen: Float = 0.0;
    let shutterclose: Float = 1.0;
    // TODO: if (shutterclose < shutteropen) {}
    let lensradius: Float = 0.0;
    let focaldistance: Float = 1e6;
    let frame: Float = 1280.0 / 720.0;
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
    let fov: Float = 30.0;
    // TODO: Float halffov ...
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
    let cam_to_world: AnimatedTransform = AnimatedTransform::new(&it, 0.0, &it, 1.0);
    let film: Film = Film::new(Point2i { x: 1280, y: 720 },
                               Bounds2f {
                                   p_min: Point2f { x: 0.0, y: 0.0 },
                                   p_max: Point2f { x: 1.0, y: 1.0 },
                               },
                               BoxFilter {
                                   radius: Vector2f { x: 0.5, y: 0.5 },
                                   inv_radius: Vector2f {
                                       x: 1.0 / 0.5,
                                       y: 1.0 / 0.5,
                                   },
                               },
                               35.0,
                               String::from("pbrt.exr"),
                               1.0,
                               std::f64::INFINITY);
    let camera_to_screen: Transform = Transform::perspective(fov, 1e-2, 1000.0);
    let perspective_camera: PerspectiveCamera = PerspectiveCamera::new(cam_to_world,
                                                                       camera_to_screen,
                                                                       screen,
                                                                       shutteropen,
                                                                       shutterclose,
                                                                       lensradius,
                                                                       focaldistance,
                                                                       fov,
                                                                       film /* ,
                                                                             * medium */);
    println!("cam_to_world = {:?}", cam_to_world);
    println!("camera_to_screen = {:?}", camera_to_screen);
    println!("perspective_camera = {:?}", perspective_camera);
}
