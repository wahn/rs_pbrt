extern crate pbrt;

use pbrt::core::lowdiscrepancy::radical_inverse;

fn main() {
    for a in 0..128 {
        let mut base_index: u16 = 2;
        println!("radical_inverse({:?}, {:?}) = {:?}", base_index, a, radical_inverse(base_index, a));
        base_index = 1;
        println!("radical_inverse({:?}, {:?}) = {:?}", base_index, a, radical_inverse(base_index, a));
        base_index = 0;
        println!("radical_inverse({:?}, {:?}) = {:?}", base_index, a, radical_inverse(base_index, a));
    }
}
