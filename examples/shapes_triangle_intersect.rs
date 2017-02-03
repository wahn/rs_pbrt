extern crate pbrt;

use pbrt::{Float, Point3f, Ray, Triangle, Vector3f};

fn main() {
    // first ray
    let triangle: Triangle = Triangle {};
    let o: Point3f = Point3f {
        x: 2.0,
        y: 1.99999988,
        z: 4.99999905,
    };
    let d: Vector3f = Vector3f {
        x: 0.101090848,
        y: -0.139514774,
        z: -0.985046268,
    };
    let r: Ray = Ray {
        o: o,
        d: d,
        t_max: std::f64::INFINITY,
    };
    let mut t_hit: Float = 0.0;
    let did_ray_interesect: bool = triangle.intersect(&r, &mut t_hit);

    println!("r = {:?}", r);
    println!("sphere.intersect(r, {:?}) = {:?}",
             t_hit,
             did_ray_interesect);
    // second ray (same origin)
    let d: Vector3f = Vector3f {
        x: 0.0997664928,
        y: -0.139622703,
        z: -0.985166073
    };
    let r: Ray = Ray {
        o: o,
        d: d,
        t_max: std::f64::INFINITY,
    };
    let mut t_hit: Float = 0.0;
    let did_ray_interesect: bool = triangle.intersect(&r, &mut t_hit);

    println!("r = {:?}", r);
    println!("sphere.intersect(r, {:?}) = {:?}",
             t_hit,
             did_ray_interesect);
    // WORK
}
