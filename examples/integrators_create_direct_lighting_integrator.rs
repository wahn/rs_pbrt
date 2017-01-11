extern crate pbrt;

use pbrt::DirectLightingIntegrator;

fn main() {
    // see directlighting.cpp CreateDirectLightingIntegrator()
    let integrator = DirectLightingIntegrator::default();

    println!("integrator = {:?}", integrator);
}
