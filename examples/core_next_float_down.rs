extern crate pbrt;

use pbrt::next_float_down;

fn main() {
    let v: f32 = 0.699999392;
    let vu: f32 = next_float_down(v);
    println!("next_float_down({:?}) = {:?})", v, vu);
}
