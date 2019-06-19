use pbrt::core::geometry::{Point3f, Vector3f};
use pbrt::core::transform::Transform;

fn main() {
    let t: Transform = Transform::translate(&Vector3f {
        x: -1.3,
        y: 0.0,
        z: 0.0,
    });
    let p: Point3f = Point3f {
        x: 2.0,
        y: 1.99999988,
        z: 4.99999905,
    };
    let mut p_error: Vector3f = Vector3f::default();
    let tp: Point3f = t.transform_point_with_error(&p, &mut p_error);

    println!("t = {:?}", t);
    println!("p = {:?}", p);
    println!(
        "tp = transform_point_with_error(p, {:?}) = {:?}",
        p_error, tp
    );
}
