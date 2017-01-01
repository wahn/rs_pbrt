extern crate pbrt;

use pbrt::Sphere;

fn main() {
    let default_sphere: Sphere = Sphere::default();

    println!("default sphere = {:?}", default_sphere);
}
