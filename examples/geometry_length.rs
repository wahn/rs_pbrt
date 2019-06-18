use pbrt::core::geometry::{Normal3, Vector2, Vector3};

fn main() {
    let float_vec2_f32 = Vector2 {
        x: 1.2_f32,
        y: 2.3_f32,
    };
    let float_vec2_f64 = Vector2 {
        x: 1.2_f64,
        y: 2.3_f64,
    };

    println!("float_vec2_f32 = {:?}", float_vec2_f32);
    println!("float_vec2_f32.length() = {:?}", float_vec2_f32.length());
    println!("float_vec2_f64 = {:?}", float_vec2_f64);
    println!("float_vec2_f64.length() = {:?}", float_vec2_f64.length());

    let float_vec3_f32 = Vector3 {
        x: 1.2_f32,
        y: 2.3_f32,
        z: 3.4_f32,
    };
    let float_vec3_f64 = Vector3 {
        x: 1.2_f64,
        y: 2.3_f64,
        z: 3.4_f64,
    };

    println!("float_vec3_f32 = {:?}", float_vec3_f32);
    println!("float_vec3_f32.length() = {:?}", float_vec3_f32.length());
    println!("float_vec3_f64 = {:?}", float_vec3_f64);
    println!("float_vec3_f64.length() = {:?}", float_vec3_f64.length());

    let float_normal3_f32 = Normal3 {
        x: 2.2_f32,
        y: 3.3_f32,
        z: 4.4_f32,
    };
    let float_normal3_f64 = Normal3 {
        x: 2.2_f64,
        y: 3.3_f64,
        z: 4.4_f64,
    };

    println!("float_normal3_f32 = {:?}", float_normal3_f32);
    println!(
        "float_normal3_f32.length() = {:?}",
        float_normal3_f32.length()
    );
    println!("float_normal3_f64 = {:?}", float_normal3_f64);
    println!(
        "float_normal3_f64.length() = {:?}",
        float_normal3_f64.length()
    );
}
