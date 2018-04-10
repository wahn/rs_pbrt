extern crate pbrt;

use pbrt::core::geometry::{Point3f, Ray, Vector3f};
use pbrt::core::transform::Transform;

fn main() {
    let t: Transform = Transform::translate(&Vector3f {
        x: -1.3,
        y: 0.0,
        z: 0.0,
    });
    let o: Point3f = Point3f {
        x: 2.0,
        y: 1.99999988,
        z: 4.99999905,
    };
    let d: Vector3f = Vector3f {
        x: -0.0607556403,
        y: -0.164096087,
        z: -0.984571517,
    };
    let r: Ray = Ray {
        o: o,
        d: d,
        t_max: std::f32::INFINITY,
        time: 0.0,
        medium: None,
        differential: None,
    };
    let mut o_error: Vector3f = Vector3f::default();
    let mut d_error: Vector3f = Vector3f::default();
    let tr: Ray = t.transform_ray_with_error(&r, &mut o_error, &mut d_error);
}
