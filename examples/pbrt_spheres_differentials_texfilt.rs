extern crate pbrt;

use pbrt::{Bounds2f, BoxFilter, Film, Float, Point2f, Point2i, Point3f, Sphere, Transform,
           Triangle, TriangleMesh, Vector2f, Vector3f};
use std::string::String;

fn main() {
    // pbrt::MakeShapes

    // trianglemesh

    // Shape "trianglemesh"  "integer indices" [0 2 1 0 3 2 ]
    let vertex_indices: Vec<usize> = vec![0_usize, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    // "point P" [-100 -1 -100 400 -1 -100 400 -1 400 -100 -1 400 ]
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
    // "float st" [ 0 0 1 0 0 1 1 1]
    let uv: Vec<Point2f> = vec![Point2f { x: 0.0, y: 0.0 },
                                Point2f { x: 1.0, y: 0.0 },
                                Point2f { x: 0.0, y: 1.0 },
                                Point2f { x: 1.0, y: 1.0 }];

    // CreateTriangleMeshShape()
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.25,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    // reverseOrientation = false
    // p graphicsState.floatTextures
    // $7 = std::map with 0 elements

    // transform mesh vertices to world space
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
    println!("triangle_mesh = {:?}", triangle_mesh);
    let mut tris: Vec<Triangle> = Vec::new();
    for i in 0..n_triangles {
        tris.push(Triangle::new(object_to_world, world_to_object, false, &triangle_mesh, i));
    }
    println!("vertex_indices = {:?}", triangle_mesh.vertex_indices);
    for i in 0..n_triangles {
        let uv = tris[i].get_uvs();
        println!("uvs[{:?}] = {:?}", i, uv);
    }

    // sphere

    // Translate -1.3 0 0
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: -1.3,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);

    // Shape "sphere"
    let radius: Float = 1.0;
    let z_min: Float = -1.0;
    let z_max: Float = 1.0;
    let phi_max: Float = 360.0;
    let sphere1: Sphere = Sphere::new(object_to_world,
                                      world_to_object,
                                      false,
                                      false,
                                      radius,
                                      z_min,
                                      z_max,
                                      phi_max);
    println!("sphere1 = {:?}", sphere1);

    // sphere

    // Translate 2.6 0 0 (not protected by Attribute[Begin|End])
    let object_to_world: Transform = object_to_world *
                                     Transform::translate(Vector3f {
        x: 2.6,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);

    // Shape "sphere"
    let sphere2: Sphere = Sphere::new(object_to_world,
                                      world_to_object,
                                      false,
                                      false,
                                      radius,
                                      z_min,
                                      z_max,
                                      phi_max);
    println!("sphere2 = {:?}", sphere2);

    // pbrtWorldEnd

    // RenderOptions::MakeIntegrator
    let xw: Float = 0.5;
    let yw: Float = 0.5;
    let box_filter = BoxFilter {
        radius: Vector2f { x: xw, y: yw },
        inv_radius: Vector2f {
            x: 1.0 / xw,
            y: 1.0 / yw,
        },
    };
    let filename: String = String::from("spheres-differentials-texfilt.exr");
    let xres = 1000;
    let yres = 500;
    let crop: Bounds2f = Bounds2f {
        p_min: Point2f { x: 0.0, y: 0.0 },
        p_max: Point2f { x: 1.0, y: 1.0 },
    };
    println!("crop = {:?}", crop);
    let film: Film = Film::new(Point2i { x: xres, y: yres },
                               crop,
                               box_filter,
                               35.0,
                               filename,
                               1.0,
                               std::f64::INFINITY);

    println!("film = {:?}", film);
}
