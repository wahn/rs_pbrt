extern crate pbrt;

use pbrt::core::geometry::Vector3f;
use pbrt::core::pbrt::Float;
use pbrt::core::transform::Transform;
use pbrt::shapes::disk::Disk;

fn main() {
    let translate: Transform = Transform::translate(&Vector3f {
        x: 0.0,
        y: 0.0,
        z: -0.01,
    });
    let inverse: Transform = Transform::inverse(&translate);
    let height: Float = 0.0;
    let radius: Float = 30.0;
    let inner_radius: Float = 0.0;
    let phi_max: Float = 360.0;
    let disk = Disk::new(
        translate,
        inverse,
        false,
        false,
        height,
        radius,
        inner_radius,
        phi_max,
    );
    println!("translate = {:?}", translate);
    println!("inverse = {:?}", inverse);
    println!("disk.height = {:?}", disk.height);
    println!("disk.radius = {:?}", disk.radius);
    println!("disk.inner_radius = {:?}", disk.inner_radius);
    println!("disk.phi_max = {:?}", disk.phi_max);
}
