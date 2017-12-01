extern crate pbrt;

use pbrt::core::geometry::{Normal3f, Point2f, Point3f, Vector3f};
use pbrt::core::transform::Transform;
use pbrt::shapes::triangle::{Triangle, TriangleMesh};
use std::sync::Arc;

fn main() {
    let translate: Transform = Transform::translate(Vector3f {
                                                        x: 0.25,
                                                        y: 0.0,
                                                        z: 0.0,
                                                    });
    let inverse: Transform = Transform::inverse(translate);
    let n_triangles: usize = 2;
    let vertex_indices: Vec<usize> = vec![0_usize, 2, 1, 0, 3, 2];
    let n_vertices: usize = 4;
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
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Normal3f> = Vec::new();
    let uv: Vec<Point2f> = vec![Point2f { x: 0.0, y: 0.0 },
                                Point2f { x: 1.0, y: 0.0 },
                                Point2f { x: 0.0, y: 1.0 },
                                Point2f { x: 1.0, y: 1.0 }];
    let triangle_mesh: TriangleMesh = TriangleMesh::new(translate,
                                                        inverse,
                                                        false,
                                                        false,
                                                        n_triangles,
                                                        vertex_indices,
                                                        n_vertices,
                                                        p,
                                                        s,
                                                        n,
                                                        uv);
    println!("translate = {:?}", translate);
    println!("inverse = {:?}", inverse);
    println!("triangle_mesh = {:?}", triangle_mesh);
    for id in 0..triangle_mesh.n_triangles {
        let triangle = Triangle::new(triangle_mesh.object_to_world,
                                     triangle_mesh.world_to_object,
                                     triangle_mesh.transform_swaps_handedness,
                                     Arc::new(triangle_mesh.clone()),
                                     id);
        println!("triangle.id = {:?}", triangle.id);
    }
}
