extern crate pbrt;

use pbrt::core::pbrt::Spectrum;
use pbrt::core::transform::Transform;
use pbrt::lights::infinite::InfiniteAreaLight;

fn main() {
    let l: Spectrum = Spectrum::new(50.0);
    let light_to_world: Transform = Transform::default();
    let n_samples: i32 = 1;
    let texmap: String = String::from("");
    let infinite_light: InfiniteAreaLight =
        InfiniteAreaLight::new(&light_to_world, &l, n_samples, texmap);
}
