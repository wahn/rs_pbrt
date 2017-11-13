extern crate pbrt;

use pbrt::core::transform::Transform;
use pbrt::geometry::Vector3f;

fn main() {
    let t: Transform = Transform::translate(Vector3f {
        x: -1.3,
        y: 0.0,
        z: 0.0,
    });
    let v: Vector3f = Vector3f {
        x: -0.0607556403,
        y: -0.164096087,
        z: -0.984571517,
    };
    let mut abs_error: Vector3f = Vector3f::default();
    let tv: Vector3f = t.transform_vector_with_error(v, &mut abs_error);

    println!("t = {:?}", t);
    println!("v = {:?}", v);
    println!("tv = transform_vector_with_error(v, {:?}) = {:?}", abs_error, tv);
}
