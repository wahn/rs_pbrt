use pbrt::core::geometry::{Point2f, Point2i};

fn main() {
    let int_origin = Point2i { x: 0, y: 0 };
    let float_origin = Point2f { x: 0.0, y: 0.0 };

    println!("int   {:?}", int_origin);
    println!("float {:?}", float_origin);
}
