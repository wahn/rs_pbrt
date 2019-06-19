use pbrt::core::transform::Matrix4x4;

fn main() {
    let m: Matrix4x4 = Matrix4x4::new(
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0,
    );

    println!("          m  = {:?}", m);
    println!("transpose(m) = {:?}", Matrix4x4::transpose(&m));
}
