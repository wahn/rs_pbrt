extern crate pbrt;

use pbrt::{Point3f, Transform, TriangleMesh, Vector3f};

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
    let triangle_mesh: TriangleMesh = TriangleMesh::new(translate,
                                                        inverse,
                                                        false,
                                                        false,
                                                        n_triangles,
                                                        vertex_indices,
                                                        n_vertices,
                                                        p);

    println!("translate = {:?}", translate);
    println!("inverse = {:?}", inverse);
    println!("triangle_mesh = {:?}", triangle_mesh);
}
