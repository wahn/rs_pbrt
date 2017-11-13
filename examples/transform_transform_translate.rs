extern crate pbrt;

use pbrt::core::transform::Transform;
use pbrt::geometry::Vector3f;

fn main() {
    let t: Transform = Transform::translate(Vector3f {
        x: -1.25,
        y: 3.5,
        z: 7.875,
    });

    println!("t = {:?}", t);
}
