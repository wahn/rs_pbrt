extern crate pbrt;

use pbrt::{Float, Transform};

fn main() {
    let theta_x: Float = 30.0;
    let t_x: Transform = Transform::rotate_x(theta_x);

    println!("Transform::rotate_x({}) = {:?}", theta_x, t_x);
}
