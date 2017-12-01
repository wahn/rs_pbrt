extern crate pbrt;

use pbrt::core::geometry::Normal3;

fn main() {
    let int_null = Normal3 { x: 0, y: 0, z: 0 };
    let float_null = Normal3 {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };

    println!("int   {:?}", int_null);
    println!("float {:?}", float_null);
}
