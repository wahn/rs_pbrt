extern crate pbrt;

use pbrt::{Float, Transform, Vector3f};

fn main() {
    let theta: Float = 30.0;
    let axis = Vector3f {
        x: 1.0,
        y: 2.0,
        z: 3.0,
    };
    let t: Transform = Transform::rotate(theta, axis);

    println!("Transform::rotate({}, {:?}) = {:?}", theta, axis, t);
}
