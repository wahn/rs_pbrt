use pbrt::core::pbrt::next_float_up;

fn main() {
    let v: f32 = -0.999999583;
    let vu: f32 = next_float_up(v);
    println!("next_float_up({:?}) = {:?})", v, vu);
}
