extern crate pbrt;

use pbrt::radical_inverse;

fn main() {
    let a: u64 = 0;
    let mut base_index: u16 = 0;
    println!("radical_inverse({:?}, {:?}) = {:?}", base_index, a, radical_inverse(base_index, a));
    base_index = 1;
    println!("radical_inverse({:?}, {:?}) = {:?}", base_index, a, radical_inverse(base_index, a));
    base_index = 2;
    println!("radical_inverse({:?}, {:?}) = {:?}", base_index, a, radical_inverse(base_index, a));
}
