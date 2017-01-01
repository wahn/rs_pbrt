extern crate pbrt;

use pbrt::Sphere;

fn main() {
    let sphere: Sphere = Sphere::new(2.0, -0.5, 0.75, 270.0);

    println!("sphere = {:?}", sphere);
}
