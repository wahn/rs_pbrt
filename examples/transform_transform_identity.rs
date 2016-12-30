extern crate pbrt;

use pbrt::Transform;

fn main() {
    let identity: Transform = Transform::default();

    println!("identity matrix = {:?}", identity);
}
