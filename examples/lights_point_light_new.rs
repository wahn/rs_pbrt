extern crate pbrt;

use pbrt::core::medium::MediumInterface;
use pbrt::core::pbrt::Spectrum;
use pbrt::core::transform::Transform;
use pbrt::lights::point::PointLight;

fn main() {
    let i: Spectrum = Spectrum::new(50.0);
    let light_to_world: Transform = Transform::default();
    let _point_light: PointLight =
        PointLight::new(&light_to_world, &MediumInterface::default(), &i);
}
