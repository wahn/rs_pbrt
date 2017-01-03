extern crate pbrt;

use pbrt::Quaternion;

fn main() {
    let default_quaternion: Quaternion = Quaternion::default();

    println!("default quaternion = {:?}", default_quaternion);
}
