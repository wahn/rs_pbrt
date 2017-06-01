extern crate pbrt;

use pbrt::{AnimatedTransform, Bounds2f, BoxFilter, BVHAccel, Checkerboard2DTexture,
           ConstantTexture, DistantLight, Film, Float, GeometricPrimitive, GlassMaterial,
           ImageTexture, ImageWrap, MatteMaterial, MirrorMaterial, PerspectiveCamera,
           PlanarMapping2D, Point2f, Point2i, Point3f, Primitive, Scene, Spectrum, Disk,
           SplitMethod, Transform, Triangle, TriangleMesh, UVMapping2D, Vector2f, Vector3f};
use std::sync::Arc;

struct SceneDescription {
    meshes: Vec<Arc<TriangleMesh>>,
    disks: Vec<Arc<Disk>>,
    lights: Vec<DistantLight>,
}

struct SceneDescriptionBuilder {
    meshes: Vec<Arc<TriangleMesh>>,
    disks: Vec<Arc<Disk>>,
    lights: Vec<DistantLight>,
}

impl SceneDescriptionBuilder {
    fn new() -> SceneDescriptionBuilder {
        SceneDescriptionBuilder {
            meshes: Vec::new(),
            disks: Vec::new(),
            lights: Vec::new(),
        }
    }
    fn add_light(&mut self,
                 light_to_world: &Transform,
                 l: &Spectrum,
                 w_light: &Vector3f)
                 -> &mut SceneDescriptionBuilder {
        let distant_light: DistantLight = DistantLight::new(light_to_world, &l, &w_light);
        self.lights.push(distant_light);
        self
    }
    fn add_mesh(&mut self,
                object_to_world: Transform,
                world_to_object: Transform,
                n_triangles: usize,
                vertex_indices: Vec<usize>,
                n_vertices: usize,
                p_ws: Vec<Point3f>,
                s: Vec<Vector3f>,
                n: Vec<Vector3f>,
                uv: Vec<Point2f>)
                -> &mut SceneDescriptionBuilder {
        let triangle_mesh = Arc::new(TriangleMesh::new(object_to_world,
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
        self.meshes.push(triangle_mesh);
        self
    }
    fn add_disk(&mut self,
                object_to_world: Transform,
                world_to_object: Transform,
                height: Float,
                radius: Float,
                inner_radius: Float,
                phi_max: Float)
                -> &mut SceneDescriptionBuilder {
        let disk = Arc::new(Disk::new(object_to_world,
                                      world_to_object,
                                      false,
                                      false,
                                      height,
                                      radius,
                                      inner_radius,
                                      phi_max));
        // println!("disk = {:?}", disk);
        self.disks.push(disk);
        self
    }
    fn finalize(self) -> SceneDescription {
        SceneDescription {
            meshes: self.meshes,
            disks: self.disks,
            lights: self.lights,
        }
    }
}

struct RenderOptions {
    primitives: Vec<Arc<Primitive + Sync + Send>>,
    triangles: Vec<Arc<Triangle>>,
    disks: Vec<Arc<Disk>>,
    lights: Vec<DistantLight>,
}

impl RenderOptions {
    fn new(scene: SceneDescription) -> RenderOptions {
        let primitives: Vec<Arc<Primitive + Sync + Send>> = Vec::new();
        let mut triangles: Vec<Arc<Triangle>> = Vec::new();
        let mut disks: Vec<Arc<Disk>> = Vec::new();
        let mut lights: Vec<DistantLight> = Vec::new();
        // lights
        for light in &scene.lights {
            lights.push(*light);
        }
        // disks
        for disk in scene.disks {
            disks.push(disk);
        }
        // meshes
        for mesh in scene.meshes {
            // create individual triangles
            for id in 0..mesh.n_triangles {
                let triangle = Arc::new(Triangle::new(mesh.object_to_world,
                                                      mesh.world_to_object,
                                                      mesh.transform_swaps_handedness,
                                                      mesh.clone(),
                                                      id));
                triangles.push(triangle);
            }
        }
        RenderOptions {
            primitives: primitives,
            triangles: triangles,
            disks: disks,
            lights: lights,
        }
    }
}


fn main() {
    let mut builder: SceneDescriptionBuilder = SceneDescriptionBuilder::new();
    // pbrt::MakeLight
    let l: Spectrum = Spectrum::new(3.141593);
    let sc: Spectrum = Spectrum::new(1.0);
    let from: Point3f = Point3f {
        x: 50.0,
        y: 50.0,
        z: 50.0,
    };
    let to: Point3f = Point3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };
    let dir: Vector3f = from - to;
    let light_to_world: Transform = Transform::default();
    let lsc: Spectrum = l * sc;
    builder.add_light(&light_to_world, &lsc, &dir);

    // pbrt::MakeShapes

    // disk

    // Translate 2.000000 -4.000000 4.000000
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 2.0,
        y: -4.0,
        z: 4.0,
    });
    // Rotate -120.000000 1.000000 0.000000 0.000000
    let theta: Float = -120.0;
    let object_to_world: Transform = object_to_world *
                                     Transform::rotate(theta,
                                                       Vector3f {
                                                           x: 1.0,
                                                           y: 0.0,
                                                           z: 0.0,
                                                       });
    let world_to_object: Transform = Transform::inverse(object_to_world);

    // Shape "disk"
    let height: Float = 0.0;
    let radius: Float = 2.0;
    let inner_radius: Float = 0.0;
    let phi_max: Float = 360.0;
    // builder.add_disk(object_to_world,
    //                  world_to_object,
    //                  height,
    //                  radius,
    //                  inner_radius,
    //                  phi_max);

    // we need camera transformation below to transform triangles
    // INFO: The order in PBRT file is different !!!
    let t: Transform = 
                       // Transform::new(0.828849,
                       //                -0.295370,
                       //                -0.475149,
                       //                0.000000,
                       //                -0.559473,
                       //                -0.437585,
                       //                -0.703924,
                       //                0.000000,
                       //                -0.000000,
                       //                0.849280,
                       //                -0.527943,
                       //                0.000000,
                       //                0.000000,
                       //                0.000000,
                       //                0.000000,
                       //                1.000000);
                       Transform::new(0.828849,
                                      -0.559473,
                                      -0.000000,
                                      0.000000,
                                      -0.295370,
                                      -0.437585,
                                      0.849280,
                                      0.000000,
                                      -0.475149,
                                      -0.703924,
                                      -0.527943,
                                      0.000000,
                                      0.000000,
                                      0.000000,
                                      0.000000,
                                      1.000000);
    let t: Transform = t *
        Transform::translate(Vector3f {
            x: -4.86,
            y: -7.2,
            z: -5.4,
        });

    // trianglemeshes
    // TODO: insert triangles below
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: -4.0,
                                   y: -1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -2.66667,
                                   y: -1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -2.66667,
                                   y: -2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -4.0,
                                   y: -2.66667,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: -4.0,
                                   y: 1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -2.66667,
                                   y: 1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -2.66667,
                                   y: 0.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -4.0,
                                   y: 0.0,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: -4.0,
                                   y: 4.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -2.66667,
                                   y: 4.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -2.66667,
                                   y: 2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -4.0,
                                   y: 2.66667,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: -2.66667,
                                   y: -2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -1.33333,
                                   y: -2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -1.33333,
                                   y: -4.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -2.66667,
                                   y: -4.0,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: -2.66667,
                                   y: 0.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -1.33333,
                                   y: 0.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -1.33333,
                                   y: -1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -2.66667,
                                   y: -1.33333,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: -2.66667,
                                   y: 2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -1.33333,
                                   y: 2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -1.33333,
                                   y: 1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -2.66667,
                                   y: 1.33333,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: -1.33333,
                                   y: -1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 0.0,
                                   y: -1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 0.0,
                                   y: -2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -1.33333,
                                   y: -2.66667,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: -1.33333,
                                   y: 1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 0.0,
                                   y: 1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 0.0,
                                   y: 0.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -1.33333,
                                   y: 0.0,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: -1.33333,
                                   y: 4.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 0.0,
                                   y: 4.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 0.0,
                                   y: 2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -1.33333,
                                   y: 2.66667,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: 0.0,
                                   y: -2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 1.33333,
                                   y: -2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 1.33333,
                                   y: -4.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 0.0,
                                   y: -4.0,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: 0.0,
                                   y: 0.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 1.33333,
                                   y: 0.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 1.33333,
                                   y: -1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 0.0,
                                   y: -1.33333,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: 0.0,
                                   y: 2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 1.33333,
                                   y: 2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 1.33333,
                                   y: 1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 0.0,
                                   y: 1.33333,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: 1.33333,
                                   y: -1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 2.66667,
                                   y: -1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 2.66667,
                                   y: -2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 1.33333,
                                   y: -2.66667,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: 1.33333,
                                   y: 1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 2.66667,
                                   y: 1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 2.66667,
                                   y: 0.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 1.33333,
                                   y: 0.0,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: 1.33333,
                                   y: 4.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 2.66667,
                                   y: 4.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 2.66667,
                                   y: 2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 1.33333,
                                   y: 2.66667,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: 2.66667,
                                   y: -2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 4.0,
                                   y: -2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 4.0,
                                   y: -4.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 2.66667,
                                   y: -4.0,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: 2.66667,
                                   y: 0.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 4.0,
                                   y: 0.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 4.0,
                                   y: -1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 2.66667,
                                   y: -1.33333,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: 2.66667,
                                   y: 2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 4.0,
                                   y: 2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 4.0,
                                   y: 1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 2.66667,
                                   y: 1.33333,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: -4.0,
                                   y: -2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -2.66667,
                                   y: -2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -2.66667,
                                   y: -4.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -4.0,
                                   y: -4.0,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: -4.0,
                                   y: 0.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -2.66667,
                                   y: 0.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -2.66667,
                                   y: -1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -4.0,
                                   y: -1.33333,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: -4.0,
                                   y: 2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -2.66667,
                                   y: 2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -2.66667,
                                   y: 1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -4.0,
                                   y: 1.33333,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: -2.66667,
                                   y: -1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -1.33333,
                                   y: -1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -1.33333,
                                   y: -2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -2.66667,
                                   y: -2.66667,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: -2.66667,
                                   y: 1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -1.33333,
                                   y: 1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -1.33333,
                                   y: 0.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -2.66667,
                                   y: 0.0,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: -2.66667,
                                   y: 4.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -1.33333,
                                   y: 4.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -1.33333,
                                   y: 2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -2.66667,
                                   y: 2.66667,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: -1.33333,
                                   y: -2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 0.0,
                                   y: -2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 0.0,
                                   y: -4.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -1.33333,
                                   y: -4.0,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: -1.33333,
                                   y: 0.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 0.0,
                                   y: 0.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 0.0,
                                   y: -1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -1.33333,
                                   y: -1.33333,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: -1.33333,
                                   y: 2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 0.0,
                                   y: 2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 0.0,
                                   y: 1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: -1.33333,
                                   y: 1.33333,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: 0.0,
                                   y: -1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 1.33333,
                                   y: -1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 1.33333,
                                   y: -2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 0.0,
                                   y: -2.66667,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: 0.0,
                                   y: 1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 1.33333,
                                   y: 1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 1.33333,
                                   y: 0.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 0.0,
                                   y: 0.0,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: 0.0,
                                   y: 4.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 1.33333,
                                   y: 4.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 1.33333,
                                   y: 2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 0.0,
                                   y: 2.66667,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: 1.33333,
                                   y: -2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 2.66667,
                                   y: -2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 2.66667,
                                   y: -4.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 1.33333,
                                   y: -4.0,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: 1.33333,
                                   y: 0.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 2.66667,
                                   y: 0.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 2.66667,
                                   y: -1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 1.33333,
                                   y: -1.33333,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: 1.33333,
                                   y: 2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 2.66667,
                                   y: 2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 2.66667,
                                   y: 1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 1.33333,
                                   y: 1.33333,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: 2.66667,
                                   y: -1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 4.0,
                                   y: -1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 4.0,
                                   y: -2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 2.66667,
                                   y: -2.66667,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: 2.66667,
                                   y: 1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 4.0,
                                   y: 1.33333,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 4.0,
                                   y: 0.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 2.66667,
                                   y: 0.0,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    {
    let vertex_indices: Vec<usize> = vec![0, 2, 1, 0, 3, 2];
    let n_triangles: usize = vertex_indices.len() / 3;
    let p: Vec<Point3f> = vec![Point3f {
                                   x: 2.66667,
                                   y: 4.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 4.0,
                                   y: 4.0,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 4.0,
                                   y: 2.66667,
                                   z: 0.0,
                               },
                               Point3f {
                                   x: 2.66667,
                                   y: 2.66667,
                                   z: 0.0,
                               }];
    let object_to_world: Transform = Transform::translate(Vector3f {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    });
    let world_to_object: Transform = Transform::inverse(object_to_world);
    let s: Vec<Vector3f> = Vec::new();
    let n: Vec<Vector3f> = Vec::new();
    let uv: Vec<Point2f> = Vec::new();
    let n_vertices: usize = p.len();
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p,
                     s,
                     n,
                     uv);
    }
    // TODO: insert triangles above

    let scene_description: SceneDescription = builder.finalize();

    // pbrtWorldEnd

    // TMP: process SceneDescription before handing primitives to BVHAccel
    let mut render_options: RenderOptions = RenderOptions::new(scene_description);
    // add triangles created above (not meshes)
    let kd = Arc::new(ConstantTexture::new(Spectrum::new(0.5)));
    let matte = Arc::new(MatteMaterial::new(kd, 0.0 as Float));
    for triangle in render_options.triangles {
        let geo_prim = Arc::new(GeometricPrimitive::new(triangle, matte.clone()));
        render_options.primitives.push(geo_prim.clone());
    }
    for disk in render_options.disks {
        let geo_prim = Arc::new(GeometricPrimitive::new(disk, matte.clone()));
        render_options.primitives.push(geo_prim.clone());
    }
    // TMP: process SceneDescription before handing primitives to BVHAccel
    // pbrt::RenderOptions::MakeScene
    let accelerator = Arc::new(BVHAccel::new(render_options.primitives, 4, SplitMethod::SAH));
    // SamplerIntegrator::Render (integrator.cpp)
    let scene: Scene = Scene::new(accelerator.clone(), render_options.lights);
    // create camera
    let it: Transform = Transform {
        m: t.m_inv.clone(),
        m_inv: t.m.clone(),
    };
    let animated_cam_to_world: AnimatedTransform = AnimatedTransform::new(&it, 0.0, &it, 1.0);
    let fov: Float = 45.0;
    let camera_to_screen: Transform = Transform::perspective(fov, 1e-2, 1000.0);
    let xres = 256;
    let yres = 256;
    let frame: Float = xres as Float / yres as Float;
    let mut screen: Bounds2f = Bounds2f::default();
    if frame > 1.0 {
        screen.p_min.x = -frame;
        screen.p_max.x = frame;
        screen.p_min.y = -1.0;
        screen.p_max.y = 1.0;
    } else {
        screen.p_min.x = -1.0;
        screen.p_max.x = 1.0;
        screen.p_min.y = -1.0 / frame;
        screen.p_max.y = 1.0 / frame;
    }
    let shutteropen: Float = 0.0;
    let shutterclose: Float = 1.0;
    let lensradius: Float = 0.0;
    let focaldistance: Float = 1e6;
    let crop: Bounds2f = Bounds2f {
        p_min: Point2f { x: 0.0, y: 0.0 },
        p_max: Point2f { x: 1.0, y: 1.0 },
    };
    let xw: Float = 0.5;
    let yw: Float = 0.5;
    let box_filter = BoxFilter {
        radius: Vector2f { x: xw, y: yw },
        inv_radius: Vector2f {
            x: 1.0 / xw,
            y: 1.0 / yw,
        },
    };
    let filename: String = String::from("teapot_area_light.exr");
    let film: Film = Film::new(Point2i { x: xres, y: yres },
                               crop,
                               box_filter,
                               35.0,
                               filename,
                               1.0,
                               std::f32::INFINITY);
    let perspective_camera: PerspectiveCamera = PerspectiveCamera::new(animated_cam_to_world,
                                                                       camera_to_screen,
                                                                       screen,
                                                                       shutteropen,
                                                                       shutterclose,
                                                                       lensradius,
                                                                       focaldistance,
                                                                       fov,
                                                                       film);
    pbrt::render(&scene, &perspective_camera);
}
