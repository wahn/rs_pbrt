extern crate pbrt;

use pbrt::core::geometry::Vector3f;
use pbrt::core::pbrt::Float;
use pbrt::core::transform::Transform;

fn main() {
    let theta: Float = 30.0;
    let axis = Vector3f {
        x: 1.0,
        y: 2.0,
        z: 3.0,
    };
    let t: Transform = Transform::rotate(theta, &axis);

    println!("Transform::rotate({}, {:?}) = {:?}", theta, axis, t);
}
