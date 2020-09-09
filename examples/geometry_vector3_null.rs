use pbrt::core::geometry::{Vector3f, Vector3i};

fn main() {
    let int_null = Vector3i { x: 0, y: 0, z: 0 };
    let float_null = Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };

    println!("int   {:?}", int_null);
    println!("float {:?}", float_null);
}
