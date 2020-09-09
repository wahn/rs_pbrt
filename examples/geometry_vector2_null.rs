use pbrt::core::geometry::{Vector2f, Vector2i};

fn main() {
    let int_null = Vector2i { x: 0, y: 0 };
    let float_null = Vector2f { x: 0.0, y: 0.0 };

    println!("int   {:?}", int_null);
    println!("float {:?}", float_null);
}
