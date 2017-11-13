extern crate pbrt;

use pbrt::core::pbrt::Float;
use pbrt::core::transform::Transform;
use pbrt::geometry::Vector3f;
use pbrt::shapes::cylinder::Cylinder;

fn main() {
    let translate: Transform = Transform::translate(Vector3f {
                                                        x: 2.0,
                                                        y: 1.0,
                                                        z: 1.6,
                                                    });
    let inverse: Transform = Transform::inverse(translate);
    let radius: Float = 1.0;
    let z_min: Float = 0.0;
    let z_max: Float = 1.0;
    let phi_max: Float = 360.0;
    let cylinder = Cylinder::new(translate,
                                 inverse,
                                 false,
                                 radius,
                                 z_min,
                                 z_max,
                                 phi_max);
    println!("translate = {:?}", translate);
    println!("inverse = {:?}", inverse);
    println!("cylinder.radius = {:?}", cylinder.radius);
    println!("cylinder.z_min = {:?}", cylinder.z_min);
    println!("cylinder.z_max = {:?}", cylinder.z_max);
    println!("cylinder.phi_max = {:?}", cylinder.phi_max);
}
