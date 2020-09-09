use pbrt::core::geometry::{Bounds3f, Bounds3i, Point3f, Point3i};

fn main() {
    let int_origin = Point3i { x: 0, y: 0, z: 0 };
    let int_xyz111 = Point3i { x: 1, y: 1, z: 1 };
    let float_origin = Point3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };
    let float_xyz111 = Point3f {
        x: 1.0,
        y: 1.0,
        z: 1.0,
    };
    let int_unit_cube = Bounds3i {
        p_min: int_origin,
        p_max: int_xyz111,
    };
    let float_unit_cube = Bounds3f {
        p_min: float_origin,
        p_max: float_xyz111,
    };

    println!("int   {:?}", int_unit_cube);
    println!("float {:?}", float_unit_cube);
}
