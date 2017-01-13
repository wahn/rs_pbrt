extern crate pbrt;

use pbrt::{Sphere, Transform, Vector3f};

fn main() {
    let translate: Transform = Transform::translate(Vector3f {
        x: -1.3,
        y: 0.0,
        z: 0.0,
    });
    let inverse: Transform = Transform::inverse(translate);
    let sphere: Sphere = Sphere::new(translate, inverse, false, false, 2.0, -0.5, 0.75, 270.0);

    println!("sphere = {:?}", sphere);
}
