extern crate pbrt;

use pbrt::core::pbrt::Float;
use pbrt::core::transform::Transform;
use pbrt::geometry::{Point3f, Ray, Vector3f};
use pbrt::shapes::Shape;
use pbrt::shapes::sphere::Sphere;

fn main() {
    // see CreateSphereShape() in sphere.cpp
    let radius: Float = 1.0;
    let z_min: Float = -radius;
    let z_max: Float = radius;
    let phi_max: Float = 360.0;
    let translate: Transform = Transform::translate(Vector3f {
        x: -1.3,
        y: 0.0,
        z: 0.0,
    });
    let inverse: Transform = Transform::inverse(translate);
    let sphere: Sphere = Sphere::new(translate,
                                     inverse,
                                     false,
                                     false,
                                     radius,
                                     z_min,
                                     z_max,
                                     phi_max);
    // see Sphere::Intersect() in sphere.cpp
    let o: Point3f = Point3f {
        x: 1.99999952,
        y: 1.99999988,
        z: 4.99999905,
    };
    let d: Vector3f = Vector3f {
        x: -0.505975366,
        y: -0.168564394,
        z: -0.84591651,
    };
    let r: Ray = Ray {
        o: o,
        d: d,
        t_max: 17.7973537,
        time: 0.0,
        differential: None,
    };
    println!("translate = {:?}", translate);
    println!("r = {:?}", r);
    if let Some((_isect, t_hit)) = Shape::intersect(&sphere, &r) {
        println!("sphere.intersect(r) = (isect, {:?})", t_hit);
    }
}
