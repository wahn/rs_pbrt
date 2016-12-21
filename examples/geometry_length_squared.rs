extern crate pbrt;

use pbrt::{Vector2, Vector3, Normal3};

fn main() {
    let int_vec2 = Vector2 { x: 1, y: 2 };

    println!("int_vec2 = {:?}", int_vec2);
    println!("int_vec2.length_squared() = {:?}",
             int_vec2.length_squared());

    let float_vec2 = Vector2 { x: 1.2, y: 2.3 };

    println!("float_vec2 = {:?}", float_vec2);
    println!("float_vec2.length_squared() = {:?}",
             float_vec2.length_squared());

    let int_vec3 = Vector3 { x: 1, y: 2, z: 3 };

    println!("int_vec3 = {:?}", int_vec3);
    println!("int_vec3.length_squared() = {:?}",
             int_vec3.length_squared());

    let float_vec3 = Vector3 {
        x: 1.2,
        y: 2.3,
        z: 3.4,
    };

    println!("float_vec3 = {:?}", float_vec3);
    println!("float_vec3.length_squared() = {:?}",
             float_vec3.length_squared());

    let int_normal3 = Normal3 { x: 2, y: 3, z: 4 };

    println!("int_normal3 = {:?}", int_normal3);
    println!("int_normal3.length_squared() = {:?}",
             int_normal3.length_squared());

    let float_normal3 = Normal3 {
        x: 2.2,
        y: 3.3,
        z: 4.4,
    };

    println!("float_normal3 = {:?}", float_normal3);
    println!("float_normal3.length_squared() = {:?}",
             float_normal3.length_squared());
}
