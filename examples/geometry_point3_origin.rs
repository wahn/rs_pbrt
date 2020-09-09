use pbrt::core::geometry::{Point3f, Point3i};

fn main() {
    let int_origin = Point3i { x: 0, y: 0, z: 0 };
    let float_origin = Point3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };

    println!("int   {:?}", int_origin);
    println!("float {:?}", float_origin);
}
