extern crate pbrt;

use pbrt::{PointLight, Spectrum, Transform};

fn main() {
    let i: Spectrum = Spectrum::new(50.0);
    let light_to_world: Transform = Transform::default();
    let point_light: PointLight = PointLight::new(&light_to_world, &i);
    println!("point_light = {:?}", point_light);
}
