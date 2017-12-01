extern crate pbrt;

use pbrt::accelerators::bvh::{BVHAccel, SplitMethod};
use pbrt::cameras::perspective::PerspectiveCamera;
use pbrt::core::filter::Filter;
use pbrt::core::geometry::{Bounds2f, Bounds2i, Normal3f, Point2f, Point2i, Point3f, Vector2f, Vector3f};
use pbrt::core::integrator::SamplerIntegrator;
use pbrt::core::light::Light;
use pbrt::core::pbrt::{Float, Spectrum};
use pbrt::core::primitive::{GeometricPrimitive, Primitive};
use pbrt::core::transform::{AnimatedTransform, Transform};
use pbrt::core::film::Film;
use pbrt::core::sampler::Sampler;
use pbrt::core::scene::Scene;
use pbrt::filters::boxfilter::BoxFilter;
use pbrt::integrators::directlighting::{DirectLightingIntegrator, LightStrategy};
use pbrt::materials::matte::MatteMaterial;
use pbrt::materials::plastic::PlasticMaterial;
use pbrt::lights::point::PointLight;
use pbrt::samplers::zerotwosequence::ZeroTwoSequenceSampler;
use pbrt::shapes::disk::Disk;
use pbrt::shapes::triangle::{Triangle, TriangleMesh};
use pbrt::textures::constant::ConstantTexture;
use std::sync::Arc;

struct SceneDescription {
    meshes: Vec<Arc<TriangleMesh>>,
    disks: Vec<Arc<Disk>>,
    lights: Vec<Arc<Light + Sync + Send>>,
}

struct SceneDescriptionBuilder {
    meshes: Vec<Arc<TriangleMesh>>,
    disks: Vec<Arc<Disk>>,
    lights: Vec<Arc<Light + Sync + Send>>,
}

impl SceneDescriptionBuilder {
    fn new() -> SceneDescriptionBuilder {
        SceneDescriptionBuilder {
            meshes: Vec::new(),
            disks: Vec::new(),
            lights: Vec::new(),
        }
    }
    fn add_point_light(&mut self,
                       light_to_world: &Transform,
                       i: &Spectrum)
                       -> &mut SceneDescriptionBuilder {
        let point_light = Arc::new(PointLight::new(light_to_world, &i));
        self.lights.push(point_light);
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
                n: Vec<Normal3f>,
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
    lights: Vec<Arc<Light + Sync + Send>>,
}

impl RenderOptions {
    fn new(scene: SceneDescription) -> RenderOptions {
        let primitives: Vec<Arc<Primitive + Sync + Send>> = Vec::new();
        let mut triangles: Vec<Arc<Triangle>> = Vec::new();
        let mut disks: Vec<Arc<Disk>> = Vec::new();
        let mut lights: Vec<Arc<Light + Sync + Send>> = Vec::new();
        // lights
        for light in &scene.lights {
            lights.push(light.clone());
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
    // we need camera transformation below to transform triangles
    // INFO: The order in PBRT file is different !!!
    let t: Transform = Transform::new(0.828849,
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
    let t_inv: Transform = Transform {
        m: t.m_inv,
        m_inv: t.m,
    };
    // LightSource "point" "color I" [ 50 50 50 ]
    let i: Spectrum = Spectrum::new(50.0);
    builder.add_point_light(&t_inv, &i);

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
    builder.add_disk(object_to_world,
                     world_to_object,
                     height,
                     radius,
                     inner_radius,
                     phi_max);

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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
        let n: Vec<Normal3f> = Vec::new();
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
    let kd = Arc::new(ConstantTexture::new(Spectrum::rgb(0.5, 0.3, 0.8)));
    let ks = Arc::new(ConstantTexture::new(Spectrum::rgb(0.2, 0.2, 0.2)));
    let roughness = Arc::new(ConstantTexture::new(1.0 as Float));
    let plastic1 = Arc::new(PlasticMaterial::new(kd, ks.clone(), roughness.clone(), true));
    let kd = Arc::new(ConstantTexture::new(Spectrum::rgb(0.8, 0.5, 0.1)));
    let plastic2 = Arc::new(PlasticMaterial::new(kd, ks.clone(), roughness.clone(), true));
    let mut triangle_count: usize = 0;
    for triangle in render_options.triangles {
        if triangle_count < 72 {
            let geo_prim = Arc::new(GeometricPrimitive::new(triangle, plastic1.clone(), None));
            render_options.primitives.push(geo_prim.clone());
        } else {
            let geo_prim = Arc::new(GeometricPrimitive::new(triangle, plastic2.clone(), None));
            render_options.primitives.push(geo_prim.clone());
        }
        triangle_count += 1;
    }
    println!("triangle_count = {}", triangle_count);
    let kd = Arc::new(ConstantTexture::new(Spectrum::new(0.0)));
    let matte = Arc::new(MatteMaterial::new(kd, 0.0 as Float));
    for disk in render_options.disks {
        let geo_prim = Arc::new(GeometricPrimitive::new(disk, matte.clone(), None));
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
    let filter: Arc<Filter + Sync + Send> = Arc::new(BoxFilter {
                                                         radius: Vector2f { x: xw, y: yw },
                                                         inv_radius: Vector2f {
                                                             x: 1.0 / xw,
                                                             y: 1.0 / yw,
                                                         },
                                                     });
    let filename: String = String::from("teapot_area_light.exr");
    let film: Arc<Film> = Arc::new(Film::new(Point2i { x: xres, y: yres },
                                             crop,
                                             filter,
                                             35.0,
                                             filename,
                                             1.0,
                                             std::f32::INFINITY));
    let camera = Box::new(PerspectiveCamera::new(animated_cam_to_world,
                                                 screen,
                                                 shutteropen,
                                                 shutterclose,
                                                 lensradius,
                                                 focaldistance,
                                                 fov,
                                                 film.clone()));
    let mut sampler: Box<Sampler + Sync + Send> = Box::new(ZeroTwoSequenceSampler::default());
    let sample_bounds: Bounds2i = film.get_sample_bounds();
    let mut integrator: Box<SamplerIntegrator + Send + Sync> =
        Box::new(DirectLightingIntegrator::new(LightStrategy::UniformSampleAll, 10, sample_bounds));
    pbrt::render(&scene, camera, &mut sampler, &mut integrator, 0_u8);
}
