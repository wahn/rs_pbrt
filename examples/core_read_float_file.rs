extern crate pbrt;

use pbrt::core::floatfile::read_float_file;
use pbrt::core::pbrt::Float;

fn main() {
    {
        let filename: String =
            String::from("/home/jan/git/self_hosted/Rust/pbrt/assets/spds/Al.k.spd");
        let mut values: Vec<Float> = Vec::new();
        println!(
            "read_float_file({:?}, ...) returns {}.",
            filename,
            read_float_file(&filename, &mut values)
        );
        println!("{:?}", values);
    }
    {
        let filename: String =
            String::from("/home/jan/git/self_hosted/Rust/pbrt/assets/spds/Al.eta.spd");
        let mut values: Vec<Float> = Vec::new();
        println!(
            "read_float_file({:?}, ...) returns {}.",
            filename,
            read_float_file(&filename, &mut values)
        );
        println!("{:?}", values);
    }
}
