extern crate pbrt;

use pbrt::core::pbrt::Float;
use pbrt::core::transform::Transform;
use pbrt::core::geometry::Vector3f;
use pbrt::shapes::sphere::Sphere;

fn main() {
    let translate: Transform = Transform::translate(Vector3f {
                                                        x: -1.3,
                                                        y: 0.0,
                                                        z: 0.0,
                                                    });
    let inverse: Transform = Transform::inverse(translate);
    let radius: Float = 1.0;
    let z_min: Float = -1.0;
    let z_max: Float = 1.0;
    let phi_max: Float = 360.0;
    let sphere = Sphere::new(translate,
                             inverse,
                             false,
                             false,
                             radius,
                             z_min,
                             z_max,
                             phi_max);
    println!("translate = {:?}", translate);
    println!("inverse = {:?}", inverse);
    println!("sphere.radius = {:?}", sphere.radius);
    println!("sphere.z_min = {:?}", sphere.z_min);
    println!("sphere.z_max = {:?}", sphere.z_max);
    println!("sphere.theta_min = {:?}", sphere.theta_min);
    println!("sphere.theta_max = {:?}", sphere.theta_max);
    println!("sphere.phi_max = {:?}", sphere.phi_max);
}
