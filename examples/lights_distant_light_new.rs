use pbrt::core::geometry::{Point3f, Vector3f};
use pbrt::core::pbrt::Spectrum;
use pbrt::core::transform::Transform;
use pbrt::lights::distant::DistantLight;

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
    let _distant_light: DistantLight = DistantLight::new(&light_to_world, &lsc, &dir);
}
