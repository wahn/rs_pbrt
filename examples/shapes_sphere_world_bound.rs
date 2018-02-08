extern crate pbrt;

use pbrt::core::geometry::{Bounds3f, Vector3f};
use pbrt::core::pbrt::Float;
use pbrt::core::shape::Shape;
use pbrt::core::transform::Transform;
use pbrt::shapes::sphere::Sphere;

fn main() {
    let radius: Float = 1.0;
    let z_min: Float = -radius;
    let z_max: Float = radius;
    let phi_max: Float = 360.0;
    let translate: Transform = Transform::translate(&Vector3f {
        x: -1.3,
        y: 0.0,
        z: 0.0,
    });
    let inverse: Transform = Transform::inverse(&translate);
    let sphere: Sphere = Sphere::new(
        translate,
        inverse,
        false,
        false,
        radius,
        z_min,
        z_max,
        phi_max,
    );
    let world_bound: Bounds3f = sphere.world_bound(); // Primitive
                                                      // println!("sphere() = {:?}", sphere);
    println!("world_bound() = {:?}", world_bound);
    let translate: Transform = Transform::translate(&Vector3f {
        x: 1.3,
        y: 0.0,
        z: 0.0,
    });
    let inverse: Transform = Transform::inverse(&translate);
    let sphere: Sphere = Sphere::new(
        translate,
        inverse,
        false,
        false,
        radius,
        z_min,
        z_max,
        phi_max,
    );
    let world_bound: Bounds3f = sphere.world_bound(); // Primitive
                                                      // println!("sphere() = {:?}", sphere);
    println!("world_bound() = {:?}", world_bound);
}
