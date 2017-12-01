extern crate pbrt;

use pbrt::core::geometry::Point2;

fn main() {
    let int_origin = Point2 { x: 0, y: 0 };
    let float_origin = Point2 { x: 0.0, y: 0.0 };

    println!("int   {:?}", int_origin);
    println!("float {:?}", float_origin);
}
