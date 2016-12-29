extern crate pbrt;

use pbrt::Transform;

fn main() {
    let identity: Transform = Transform::new();

    println!("identity matrix = {:?}", identity);
}
