extern crate pbrt;

use pbrt::{BVHAccel, Checkerboard2DTexture, ConstantTexture, DistantLight, Float,
           GeometricPrimitive, GlassMaterial, ImageTexture, ImageWrap, MatteMaterial,
           MirrorMaterial, PlanarMapping2D, Point2f, Point3f, Primitive, Scene, Spectrum, Sphere,
           SplitMethod, Transform, Triangle, TriangleMesh, UVMapping2D, Vector3f};
use std::string::String;
use std::sync::Arc;

struct SceneDescription {
    meshes: Vec<Arc<TriangleMesh>>,
    spheres: Vec<Arc<Sphere>>,
    lights: Vec<DistantLight>,
}

struct SceneDescriptionBuilder {
    meshes: Vec<Arc<TriangleMesh>>,
    spheres: Vec<Arc<Sphere>>,
    lights: Vec<DistantLight>,
}

impl SceneDescriptionBuilder {
    fn new() -> SceneDescriptionBuilder {
        SceneDescriptionBuilder {
            meshes: Vec::new(),
            spheres: Vec::new(),
            lights: Vec::new(),
        }
    }
    fn add_light(&mut self,
                 light_to_world: &Transform,
                 l: &Spectrum,
                 w_light: &Vector3f)
                 -> &mut SceneDescriptionBuilder {
        let distant_light: DistantLight = DistantLight::new(light_to_world, &l, &w_light);
        println!("distant_light = {:?}", distant_light);
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
        println!("triangle_mesh = {:?}", triangle_mesh);
        println!("vertex_indices = {:?}", triangle_mesh.vertex_indices);
        self.meshes.push(triangle_mesh);
        self
    }
    fn add_sphere(&mut self,
                  object_to_world: Transform,
                  world_to_object: Transform,
                  radius: Float,
                  z_min: Float,
                  z_max: Float,
                  phi_max: Float)
                  -> &mut SceneDescriptionBuilder {
        let sphere = Arc::new(Sphere::new(object_to_world,
                                          world_to_object,
                                          false,
                                          false,
                                          radius,
                                          z_min,
                                          z_max,
                                          phi_max));
        // println!("sphere = {:?}", sphere);
        self.spheres.push(sphere);
        self
    }
    fn finalize(self) -> SceneDescription {
        SceneDescription {
            meshes: self.meshes,
            spheres: self.spheres,
            lights: self.lights,
        }
    }
}

struct RenderOptions {
    primitives: Vec<Arc<Primitive + Sync + Send>>,
    triangles: Vec<Arc<Triangle>>,
    spheres: Vec<Arc<Sphere>>,
    lights: Vec<DistantLight>,
}

impl RenderOptions {
    fn new(scene: SceneDescription) -> RenderOptions {
        let primitives: Vec<Arc<Primitive + Sync + Send>> = Vec::new();
        let mut triangles: Vec<Arc<Triangle>> = Vec::new();
        let mut spheres: Vec<Arc<Sphere>> = Vec::new();
        let mut lights: Vec<DistantLight> = Vec::new();
        // lights
        for light in &scene.lights {
            lights.push(*light);
        }
        // spheres
        for sphere in scene.spheres {
            spheres.push(sphere);
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
            spheres: spheres,
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
        x: 0.0,
        y: 10.0,
        z: 0.0,
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
    println!("########");
    println!("# mesh #");
    println!("########");
    builder.add_mesh(object_to_world,
                     world_to_object,
                     n_triangles,
                     vertex_indices,
                     n_vertices,
                     p_ws, // in world space
                     s, // empty
                     n, // empty
                     uv);

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
    println!("############");
    println!("# sphere 1 #");
    println!("############");
    builder.add_sphere(object_to_world,
                       world_to_object,
                       radius,
                       z_min,
                       z_max,
                       phi_max);

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
    println!("############");
    println!("# sphere 2 #");
    println!("############");
    builder.add_sphere(object_to_world,
                       world_to_object,
                       radius,
                       z_min,
                       z_max,
                       phi_max);

    let scene_description: SceneDescription = builder.finalize();

    // pbrtWorldEnd

    // TMP: process SceneDescription before handing primitives to BVHAccel
    let mut render_options: RenderOptions = RenderOptions::new(scene_description);
    // add triangles created above (not meshes)
    // let tex1 = Arc::new(ConstantTexture { value: Spectrum::new(0.0) });
    // let tex2 = Arc::new(ConstantTexture { value: Spectrum::new(1.0) });
    // let mapping = Box::new(PlanarMapping2D {
    //     vs: Vector3f {
    //         x: 1.0,
    //         y: 0.0,
    //         z: 0.0,
    //     },
    //     vt: Vector3f {
    //         x: 0.0,
    //         y: 0.0,
    //         z: 1.0,
    //     },
    //     ds: 0.0 as Float,
    //     dt: 1.0 as Float,
    // });
    let uscale: Float = 100.0;
    let vscale: Float = 100.0;
    let udelta: Float = 0.0;
    let vdelta: Float = 0.0;
    let mapping = Box::new(UVMapping2D {
        su: uscale,
        sv: vscale,
        du: udelta,
        dv: vdelta,
    });
    let filename: String = String::from("assets/textures/lines.png");
    let do_trilinear: bool = false;
    let max_aniso: Float = 8.0;
    let wrap_mode: ImageWrap = ImageWrap::Repeat;
    let scale: Float = 1.0;
    let gamma: bool = true;
    let lines_tex = Arc::new(ImageTexture::new(mapping,
                                               filename,
                                               do_trilinear,
                                               max_aniso,
                                               wrap_mode,
                                               scale,
                                               gamma));
    // let checker = Arc::new(Checkerboard2DTexture::new(mapping, tex1, tex2));
    // let kd = Arc::new(ConstantTexture::new(Spectrum::new(0.5)));
    // let matte = Arc::new(MatteMaterial::new(kd, 0.0 as Float));
    // let matte = Arc::new(MatteMaterial::new(checker, 0.0 as Float));
    let matte = Arc::new(MatteMaterial::new(lines_tex, 0.0 as Float));
    let mirror = Arc::new(MirrorMaterial { kr: Spectrum::new(0.9) });
    let glass = Arc::new(GlassMaterial {
        kr: Spectrum::new(1.0),
        kt: Spectrum::new(1.0),
        u_roughness: 0.0 as Float,
        v_roughness: 0.0 as Float,
        index: 0.0 as Float,
        remap_roughness: true,
    });
    for triangle in render_options.triangles {
        let geo_prim = Arc::new(GeometricPrimitive::new(triangle, matte.clone()));
        render_options.primitives.push(geo_prim.clone());
    }
    let mut sphere_counter: u8 = 0;
    for sphere in render_options.spheres {
        if sphere_counter == 0 {
            let geo_prim = Arc::new(GeometricPrimitive::new(sphere, mirror.clone()));
            render_options.primitives.push(geo_prim.clone());
        } else {
            let geo_prim = Arc::new(GeometricPrimitive::new(sphere, glass.clone()));
            render_options.primitives.push(geo_prim.clone());
        }
        sphere_counter += 1;
    }
    // TMP: process SceneDescription before handing primitives to BVHAccel
    // pbrt::RenderOptions::MakeScene
    let accelerator = Arc::new(BVHAccel::new(render_options.primitives, 4, SplitMethod::SAH));
    println!("###############");
    println!("# accelerator #");
    println!("###############");
    // SamplerIntegrator::Render (integrator.cpp)
    let scene: Scene = Scene::new(accelerator.clone(), render_options.lights);
    pbrt::render(&scene);
}
