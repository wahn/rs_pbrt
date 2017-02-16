extern crate pbrt;

use pbrt::{Float, Point2f, Point3f, Primitive, Transform, Ray, SurfaceInteraction, Triangle,
           TriangleMesh, Vector3f};

fn main() {
    let vertex_indices: Vec<usize> = vec![0_usize, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: -100.0,
                                   y: -1.0,
                                   z: -100.0,
                               },
                               Point3f {
                                   x: 400.0,
                                   y: -1.0,
                                   z: -100.0,
                               },
                               Point3f {
                                   x: 400.0,
                                   y: -1.0,
                                   z: 400.0,
                               },
                               Point3f {
                                   x: -100.0,
                                   y: -1.0,
                                   z: 400.0,
                               }];
    let uv: Vec<Point2f> = vec![Point2f { x: 0.0, y: 0.0 },
                                Point2f { x: 1.0, y: 0.0 },
                                Point2f { x: 0.0, y: 1.0 },
                                Point2f { x: 1.0, y: 1.0 }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.25,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let mut p_ws: Vec<Point3f> = Vec::new();
    let n_vertices: usize = p.len();
    for i in 0..n_vertices {
        p_ws.push(object_to_world.transform_point(p[i]));
    }
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let triangle_mesh: TriangleMesh = TriangleMesh::new(object_to_world,
                                                        world_to_object,
                                                        false,
                                                        false,
                                                        n_triangles,
                                                        vertex_indices,
                                                        n_vertices,
                                                        p_ws, // in world space
                                                        s, // empty
                                                        n, // empty
                                                        uv);
    let mut tris: Vec<Triangle> = Vec::new();
    for i in 0..n_triangles {
        tris.push(Triangle::new(object_to_world, world_to_object, false, &triangle_mesh, i));
    }
    // first ray
    let ref triangle: Triangle = tris[0];
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
        time: 0.0,
    };
    let mut t_hit: Float = 0.0;
    let mut isect: SurfaceInteraction = SurfaceInteraction::default();
    let did_ray_interesect: bool = triangle.intersect(&r, &mut t_hit, &mut isect); // Primitive

    println!("r = {:?}", r);
    println!("sphere.intersect(r, {:?}) = {:?}",
             t_hit,
             did_ray_interesect);
    // second ray (same origin)
    let d: Vector3f = Vector3f {
        x: 0.0997664928,
        y: -0.139622703,
        z: -0.985166073,
    };
    let r: Ray = Ray {
        o: o,
        d: d,
        t_max: std::f64::INFINITY,
        time: 0.0,
    };
    let mut t_hit: Float = 0.0;
    let did_ray_interesect: bool = triangle.intersect(&r, &mut t_hit, &mut isect); // Primitive

    println!("r = {:?}", r);
    println!("sphere.intersect(r, {:?}) = {:?}",
             t_hit,
             did_ray_interesect);
}
