extern crate pbrt;

use pbrt::{DistantLight, Point3f, Spectrum, Transform, Vector3f};

fn main() {
    let l: Spectrum = Spectrum::new(3.141593);
    let sc: Spectrum = Spectrum::new(1.0);
    let from: Point3f = Point3f {
        x: 0.0,
        y: 10.0,
        z: 0.0,
    };
    let to: Point3f = Point3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };
    let dir: Vector3f = from - to;
    let light_to_world: Transform = Transform::default();
    let lsc: Spectrum = l * sc;
    let distant_light: DistantLight = DistantLight::new(&light_to_world, &lsc, &dir);
    println!("distant_light = {:?}", distant_light);
}
