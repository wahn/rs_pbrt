use pbrt::core::transform::Matrix4x4;

fn main() {
    let identity: Matrix4x4 = Matrix4x4::default();

    println!("identity matrix = {:?}", identity);
}
