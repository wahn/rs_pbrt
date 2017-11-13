extern crate pbrt;

use pbrt::core::pbrt::Float;
use pbrt::core::transform::Transform;

fn main() {
    let theta_y: Float = 45.0;
    let t_y: Transform = Transform::rotate_y(theta_y);

    println!("Transform::rotate_y({}) = {:?}", theta_y, t_y);
}
