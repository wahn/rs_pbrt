extern crate pbrt;

use pbrt::{Float, Point3f, Ray, Sphere, SurfaceInteraction, Transform, Vector3f};

fn main() {
    // see CreateSphereShape() in sphere.cpp
    let radius: Float = 1.0;
    let z_min: Float = -radius;
    let z_max: Float = radius;
    let phi_max: Float = 360.0;
    let translate: Transform = Transform::translate(Vector3f {
        x: -1.3,
        y: 0.0,
        z: 0.0,
    });
    let inverse: Transform = Transform::inverse(translate);
    let sphere: Sphere = Sphere::new(translate,
                                     inverse,
                                     false,
                                     false,
                                     radius,
                                     z_min,
                                     z_max,
                                     phi_max);
    // see Sphere::Intersect() in sphere.cpp
    let o: Point3f = Point3f {
        x: 1.99999952,
        y: 1.99999988,
        z: 4.99999905,
    };
    let d: Vector3f = Vector3f {
        x: -0.505975366,
        y: -0.168564394,
        z: -0.84591651,
    };
    let r: Ray = Ray {
        o: o,
        d: d,
        t_max: 17.7973537,
        time: 0.0,
        differential: None,
    };
    let mut t_hit: Float = 0.0;
    let mut isect: SurfaceInteraction = SurfaceInteraction::default();
    let did_ray_interesect: bool = sphere.intersect_hit(&r, &mut t_hit, &mut isect); // Primitive

    println!("translate = {:?}", translate);
    println!("r = {:?}", r);
    println!("sphere.intersect(r, {:?}) = {:?}",
             t_hit,
             did_ray_interesect);
}
