use pbrt::core::geometry::{Bounds2f, Bounds2i, Point2f, Point2i};

fn main() {
    let int_origin = Point2i { x: 0, y: 0 };
    let int_xy11 = Point2i { x: 1, y: 1 };
    let float_origin = Point2f { x: 0.0, y: 0.0 };
    let float_xy11 = Point2f { x: 1.0, y: 1.0 };
    let int_unit_cube = Bounds2i {
        p_min: int_origin,
        p_max: int_xy11,
    };
    let float_unit_cube = Bounds2f {
        p_min: float_origin,
        p_max: float_xy11,
    };

    println!("int   {:?}", int_unit_cube);
    println!("float {:?}", float_unit_cube);
}
