use pbrt::core::geometry::{Normal3f, Normal3i};

fn main() {
    let int_null = Normal3i { x: 0, y: 0, z: 0 };
    let float_null = Normal3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };

    println!("int   {:?}", int_null);
    println!("float {:?}", float_null);
}
