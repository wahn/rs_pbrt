extern crate pbrt;

use pbrt::{Point2f, Point3f, Transform, Triangle, TriangleMesh, Vector3f};
use std::rc::Rc;

fn main() {
    // pbrt::MakeShapes

    // trianglemesh

    // Shape "trianglemesh"  "integer indices" [0 2 1 0 3 2 ]
    let vertex_indices: Vec<usize> = vec![0_usize, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    println!("n_triangles = {:?}", n_triangles);
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
    println!("object_to_world = {:?}", object_to_world);
    println!("world_to_object = {:?}", world_to_object);
    // reverseOrientation = false
    // p graphicsState.floatTextures
    // $7 = std::map with 0 elements

    // transform mesh vertices to world space
    let mut p_ws: Vec<Point3f> = Vec::new();
    let n_vertices: usize = p.len();
    println!("n_vertices = {:?}", n_vertices);
    for i in 0..n_vertices {
        p_ws.push(object_to_world.transform_point(p[i]));
    }
    println!("p_ws = {:?}", p_ws);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let triangle_mesh: Rc<TriangleMesh> = Rc::new(TriangleMesh::new(object_to_world,
                                                                    world_to_object,
                                                                    false,
                                                                    false,
                                                                    n_triangles,
                                                                    vertex_indices,
                                                                    n_vertices,
                                                                    p_ws, // in world space
                                                                    s, // empty
                                                                    n, // empty
                                                                    uv));
    let mut tris: Vec<Triangle> = Vec::new();
    for i in 0..n_triangles {
        tris.push(Triangle::new(object_to_world,
                                world_to_object,
                                false,
                                triangle_mesh.clone(),
                                i));
    }
    println!("tris = {:?}", tris);

    // sphere

    // Shape "sphere"

    // sphere

    // Shape "sphere"
}
