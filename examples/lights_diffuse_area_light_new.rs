use pbrt::core::geometry::Vector3f;
use pbrt::core::medium::MediumInterface;
use pbrt::core::pbrt::{Float, Spectrum};
use pbrt::core::shape::Shape;
use pbrt::core::transform::Transform;
use pbrt::lights::diffuse::DiffuseAreaLight;
use pbrt::shapes::disk::Disk;
use std::sync::Arc;

fn main() {
    let t: Transform = Transform::translate(&Vector3f {
        x: 2.0,
        y: -4.0,
        z: 4.0,
    });
    let theta: Float = -120.0;
    let axis = Vector3f {
        x: 1.0,
        y: 0.0,
        z: 0.0,
    };
    let r: Transform = Transform::rotate(theta, &axis);
    let light_to_world: Transform = Transform::default() * r * t;
    let inverse: Transform = Transform::inverse(&light_to_world);
    let l_emit: Spectrum = Spectrum::new(8.0);
    let n_samples: i32 = 16;
    let height: Float = 0.0;
    let radius: Float = 2.0;
    let inner_radius: Float = 0.0;
    let phi_max: Float = 360.0;
    let shape: Arc<dyn Shape + Send + Sync> = Arc::new(Disk::new(
        light_to_world,
        inverse,
        false,
        height,
        radius,
        inner_radius,
        phi_max,
    ));
    let two_sided: bool = false;
    let diffuse_area_light: DiffuseAreaLight = DiffuseAreaLight::new(
        &light_to_world,
        &MediumInterface::default(),
        &l_emit,
        n_samples,
        shape,
        two_sided,
    );
    println!(
        "diffuse_area_light.l_emit = {:?}",
        diffuse_area_light.l_emit
    );
    println!(
        "diffuse_area_light.two_sided = {:?}",
        diffuse_area_light.two_sided
    );
    println!("diffuse_area_light.area = {:?}", diffuse_area_light.area);
}
