extern crate pbrt;

use pbrt::{Bounds3f, Point2f, Point3f, Primitive, Transform, Triangle, TriangleMesh, Vector3f};

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
    println!("tris[0] = {:?}", tris[0]);
    let p0: Point3f = triangle_mesh.p[triangle_mesh.vertex_indices[0]];
    let p1: Point3f = triangle_mesh.p[triangle_mesh.vertex_indices[1]];
    let p2: Point3f = triangle_mesh.p[triangle_mesh.vertex_indices[2]];
    println!("p0 = {:?}", p0);
    println!("p1 = {:?}", p1);
    println!("p2 = {:?}", p2);
    let world_bound: Bounds3f = tris[1].world_bound(); // Primitive
    println!("tris[1].world_bound() = {:?}", world_bound);
    println!("tris[1] = {:?}", tris[1]);
    let p0: Point3f = triangle_mesh.p[triangle_mesh.vertex_indices[3]];
    let p1: Point3f = triangle_mesh.p[triangle_mesh.vertex_indices[4]];
    let p2: Point3f = triangle_mesh.p[triangle_mesh.vertex_indices[5]];
    println!("p0 = {:?}", p0);
    println!("p1 = {:?}", p1);
    println!("p2 = {:?}", p2);
    let world_bound: Bounds3f = tris[1].world_bound(); // Primitive
    println!("tris[1].world_bound() = {:?}", world_bound);
}
