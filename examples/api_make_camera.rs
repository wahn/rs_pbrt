extern crate pbrt;

use pbrt::core::geometry::{Point3f, Vector3f};
use pbrt::core::transform::{AnimatedTransform, Transform};

fn main() {
    // LookAt 2 2 5  0 -.4 0  0 1 0 (see spheres-differentials-texfilt.pbrt)
    let pos = Point3f {
        x: 2.0,
        y: 2.0,
        z: 5.0,
    };
    let look = Point3f {
        x: 0.0,
        y: -0.4,
        z: 0.0,
    };
    let up = Vector3f {
        x: 0.0,
        y: 1.0,
        z: 0.0,
    };
    let t: Transform = Transform::look_at(pos, look, up);
    let it: Transform = Transform {
        m: t.m_inv.clone(),
        m_inv: t.m.clone(),
    };
    let animated_cam_2_world: AnimatedTransform = AnimatedTransform::new(&it, 0.0, &it, 1.0);

    println!("Transform::look_at({:?}, {:?}, {:?}) = {:?}",
             pos,
             look,
             up,
             t);
    println!("it = {:?}", it);
    println!("animated_cam_2_world = {:?}", animated_cam_2_world);
    // TODO: check results vs. C++ version
}
