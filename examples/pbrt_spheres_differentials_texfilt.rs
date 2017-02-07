extern crate pbrt;

use pbrt::{Point3f, Transform, Vector3f};

fn main() {
    // pbrt::MakeShapes

    // trianglemesh

    // Shape "trianglemesh"  "integer indices" [0 2 1 0 3 2 ]
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

    // CreateTriangleMeshShape()
    // p *object2world
    // $4 = {m = {m = {{1, 0, 0, 0.25}, {0, 1, 0, 0}, {0, 0, 1, 0},
    // {0, 0, 0, 1}}}, mInv = {m = {{1, 0, 0, -0.25}, {0, 1, 0, 0},
    // {0, 0, 1, 0}, {0, 0, 0, 1}}}}
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.25,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    println!("object_to_world = {:?}", object_to_world);
    println!("world_to_object = {:?}", world_to_object);
    // p *world2object
    // $5 = {m = {m = {{1, 0, 0, -0.25}, {0, 1, 0, 0}, {0, 0, 1, 0},
    // {0, 0, 0, 1}}}, mInv = {m = {{1, 0, 0, 0.25}, {0, 1, 0, 0}, {0,
    // 0, 1, 0}, {0, 0, 0, 1}}}}
    // reverseOrientation = false
    // p graphicsState.floatTextures
    // $7 = std::map with 0 elements

    // transform mesh vertices to world space
    let mut p_ws: Vec<Point3f> = Vec::new();
    let p_len: usize = p.len();
    println!("p_len = {:?}", p_len);
    for i in 0..p_len {
        p_ws.push(object_to_world.transform_point(p[i]));
    }
    println!("p_ws = {:?}", p_ws);
    // sphere

    // Shape "sphere"

    // sphere

    // Shape "sphere"
}
