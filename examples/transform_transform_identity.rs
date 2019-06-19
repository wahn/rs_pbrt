use pbrt::core::transform::Transform;

fn main() {
    let identity: Transform = Transform::default();

    println!("identity transform = {:?}", identity);
}
