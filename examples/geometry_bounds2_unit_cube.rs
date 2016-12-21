extern crate pbrt;

use pbrt::{Bounds2, Point2};

fn main() {
    let int_origin = Point2 { x: 0, y: 0 };
    let int_xy11 = Point2 { x: 1, y: 1 };
    let float_origin = Point2 { x: 0.0, y: 0.0 };
    let float_xy11 = Point2 { x: 1.0, y: 1.0 };
    let int_unit_cube = Bounds2 {
        p_min: int_origin,
        p_max: int_xy11,
    };
    let float_unit_cube = Bounds2 {
        p_min: float_origin,
        p_max: float_xy11,
    };

    println!("int   {:?}", int_unit_cube);
    println!("float {:?}", float_unit_cube);
}
