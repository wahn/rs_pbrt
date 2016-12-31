extern crate pbrt;

use pbrt::{Float, Transform};

fn main() {
    let theta_z: Float = 60.0;
    let t_z: Transform = Transform::rotate_z(theta_z);

    println!("Transform::rotate_z({}) = {:?}", theta_z, t_z);
}
