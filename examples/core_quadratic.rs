extern crate pbrt;

use pbrt::{quadratic, Float, Point3f, Ray, Vector3f};

fn main() {
    // see bool Sphere::Intersect(...) in sphere.cpp
    let radius: Float = 1.0;
    let ray: Ray = Ray {
        o: Point3f {
            x: 0.699999988,
            y: 1.99999976,
            z: 4.99999809,
        },
        d: Vector3f {
            x: -0.0607556403,
            y: -0.164096087,
            z: -0.984571517,
        },
        t_max: std::f64::INFINITY,
    };
    let ox = ray.o.x;
    let oy = ray.o.y;
    let oz = ray.o.z;
    let dx = ray.d.x;
    let dy = ray.d.y;
    let dz = ray.d.z;
    let a = dx * dx + dy * dy + dz * dz;
    let b = 2.0 * (dx * ox + dy * oy + dz * oz);
    let c = ox * ox + oy * oy + oz * oz - radius * radius;
    // return values of quadratic
    let mut t0: Float = 0.0;
    let mut t1: Float = 0.0;
    // expected to return false
    let is_quadratic = quadratic(a, b, c, &mut t0, &mut t1);
    println!("ray = {:?}", ray);
    if is_quadratic {
        println!("quadratic({:?}, {:?}, {:?}, {:?}, {:?}) = {:?}",
                 a,
                 b,
                 c,
                 t0,
                 t1,
                 is_quadratic);
    } else {
        // values of _t0_ and _t1_ therefore don't matter
        println!("quadratic({:?}, {:?}, {:?}, ...) = {:?}",
                 a,
                 b,
                 c,
                 is_quadratic);
    }
    let radius: Float = 1.0;
    let ray: Ray = Ray {
        o: Point3f {
            x: 3.299999,
            y: 1.99999964,
            z: 4.99999809,
        },
        d: Vector3f {
            x: -0.505975366,
            y: -0.168564394,
            z: -0.84591651,
        },
        t_max: std::f64::INFINITY,
    };
    let ox = ray.o.x;
    let oy = ray.o.y;
    let oz = ray.o.z;
    let dx = ray.d.x;
    let dy = ray.d.y;
    let dz = ray.d.z;
    let a = dx * dx + dy * dy + dz * dz;
    let b = 2.0 * (dx * ox + dy * oy + dz * oz);
    let c = ox * ox + oy * oy + oz * oz - radius * radius;
    // return values of quadratic
    let mut t0: Float = 0.0;
    let mut t1: Float = 0.0;
    // expected to return true
    let is_quadratic = quadratic(a, b, c, &mut t0, &mut t1);
    // (gdb) p t0
    // $8 = {v = 6.18103886, low = 6.18102741, high = 6.18105078, vPrecise = 6.1810388691983602651135176753172118}
    // (gdb) p t1
    // $9 = {v = 6.29181957, low = 6.29180765, high = 6.29183102, vPrecise = 6.2918196725310140370035494328249115}
    println!("ray = {:?}", ray);
    if is_quadratic {
        // compare vs. the gdb output above
        println!("quadratic({:?}, {:?}, {:?}, {:?}, {:?}) = {:?}",
                 a,
                 b,
                 c,
                 t0,
                 t1,
                 is_quadratic);
    } else {
        println!("quadratic({:?}, {:?}, {:?}, ...) = {:?}",
                 a,
                 b,
                 c,
                 is_quadratic);
    }
}
