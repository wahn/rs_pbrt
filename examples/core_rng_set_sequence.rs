extern crate pbrt;

use pbrt::core::rng::Rng;

fn main() {
    let mut rng: Rng = Rng::new();
    rng.set_sequence(0_u64);
    println!("rng = {:?}", rng);
}
