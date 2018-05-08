extern crate pbrt;

use pbrt::core::geometry::spherical_direction_vec3;
use pbrt::core::geometry::Vector3f;
use pbrt::core::pbrt::Float;

fn main() {
    let sin_theta: Float = 0.992971241;
    let cos_theta: Float = 0.11835587;
    let phi: Float = 0.790210605;
    let x: Vector3f = Vector3f {
        x: 0.0 as Float,
        y: -0.808129907 as Float,
        z: 0.589004219 as Float,
    };
    let y: Vector3f = Vector3f {
        x: -0.874144614 as Float,
        y: -0.286059141 as Float,
        z: -0.392480969 as Float,
    };
    let z: Vector3f = Vector3f {
        x: -0.485665679 as Float,
        y: 0.514874935 as Float,
        z: 0.706422448 as Float,
    };
    let sp: Vector3f = spherical_direction_vec3(sin_theta, cos_theta, phi, &x, &y, &z);
    println!("sp = {:?}", sp);
}
