use pbrt::core::quaternion::Quaternion;

fn main() {
    let default_quaternion: Quaternion = Quaternion::default();

    println!("default quaternion = {:?}", default_quaternion);
}
