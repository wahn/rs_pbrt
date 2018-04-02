extern crate pbrt;

use pbrt::core::transform::Transform;

fn main() {
    let t: Transform = Transform::new(
        2.0,
        0.0,
        0.0,
        -1.25,
        0.0,
        -8.0,
        0.0,
        3.5,
        0.0,
        0.0,
        1.0,
        7.875,
        0.0,
        0.0,
        0.0,
        1.0,
    );

    println!("t = {:?}", t);
}
