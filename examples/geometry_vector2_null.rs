extern crate pbrt;

use pbrt::core::geometry::Vector2;

fn main() {
    let int_null = Vector2 { x: 0, y: 0 };
    let float_null = Vector2 { x: 0.0, y: 0.0 };

    println!("int   {:?}", int_null);
    println!("float {:?}", float_null);
}
