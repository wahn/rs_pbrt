//! # pbrt
//!
//! [Rust][rust] crate to implement at least parts of the [PBRT
//! book][book]'s C++ code. You can find a copy of the current code
//! [here][repo].
//!
//! [rust]: https://www.rust-lang.org/en-US
//! [book]: http://www.pbrt.org
//! [repo]: https://github.com/wahn/rs_pbrt
//!
//! ## Scene
//!
//! As the scene file is parsed, objects are created that represent
//! the lights and geometric primitives in the scene. These are all
//! stored in the **Scene** object.
//!
//! ## Vectors
//!
//! **pbrt** provides both 2D and 3D **vector** classes. Both are
//! parameterized by the type of the underlying vector element, thus
//! making it easy to instantiate vectors of both integer and
//! floating-point types.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::Vector3;
//!
//! fn main() {
//!     let int_null = Vector3 { x: 0, y: 0, z: 0 };
//!     let float_null = Vector3 {
//!         x: 0.0,
//!         y: 0.0,
//!         z: 0.0,
//!     };
//!
//!     println!("int   {:?}", int_null);
//!     println!("float {:?}", float_null);
//! }
//! ```
//!
//! ## Points
//!
//! A **point** is a zero-dimensional location in 2D or 3D space. The
//! **Point2** and **Point3** classes in **pbrt** represent points in
//! the obvious way: using x, y, z (in 3D) coordinates with respect to
//! a coordinate system. Although the same representation is used for
//! vectors, the fact that a point represents a position whereas a
//! vector represents a direction leads to a number of important
//! differences in how they are treated.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::Point3;
//!
//! fn main() {
//!     let int_origin = Point3 { x: 0, y: 0, z: 0 };
//!     let float_origin = Point3 {
//!         x: 0.0,
//!         y: 0.0,
//!         z: 0.0,
//!     };
//!
//!     println!("int   {:?}", int_origin);
//!     println!("float {:?}", float_origin);
//! }
//! ```
//!
//! ## Normals
//!
//! A surface **normal** (or just normal) is a vector that is
//! perpendicular to a surface at a particular position. It can be
//! defined as the cross product of any two nonparallel vectors that
//! are tangent to the surface at a point. Although normals are
//! superficially similar to vectors, it is important to distinguish
//! between the two of them: because normals are defined in terms of
//! their relationship to a particular surface, they behave
//! differently than vectors in some situations, particularly when
//! applying transformations.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::Normal3;
//!
//! fn main() {
//!     let int_null = Normal3 { x: 0, y: 0, z: 0 };
//!     let float_null = Normal3 {
//!         x: 0.0,
//!         y: 0.0,
//!         z: 0.0,
//!     };
//!
//!     println!("int   {:?}", int_null);
//!     println!("float {:?}", float_null);
//! }
//! ```
//!
//! ## Rays
//!
//! A **ray** is a semi-infinite line specified by its origin and
//! direction. **pbrt** represents a **Ray** with a **Point3f** for
//! the origin and a **Vector3f** for the direction. We only need rays
//! with floating-point origins and directions, so **Ray** isn't a
//! template class parameterized by an arbitrary type, as points,
//! vectors, and normals were.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::{Ray, Point3f, Vector3f};
//!
//! fn main() {
//!     let origin = Point3f {
//!         x: -5.5,
//!         y: 2.75,
//!         z: 0.0,
//!     };
//!     let direction = Vector3f {
//!         x: 1.0,
//!         y: -8.75,
//!         z: 2.25,
//!     };
//!     let ray = Ray {
//!         o: origin,
//!         d: direction,
//!         t_max: std::f32::INFINITY,
//!         time: 0.0,
//!         differential: None,
//!     };
//!
//!     println!("{:?}", ray);
//! }
//! ```
//!
//! ### RayDifferentials
//!
//! **RayDifferential** is a subclass of **Ray** that contains
//! additional information about two auxiliary rays. These extra rays
//! represent camera rays offset by one sample in the *x* and *y*
//! direction from the main ray on the film plane. By determining the
//! area that these three rays project on an object being shaded, the
//! **Texture** can estimate an area to average over for proper
//! antialiasing.
//!
//! In Rust we don't have inheritance, therefore we use an
//! **Option<RayDifferential>** in the **Ray** struct, which means the
//! additional information can be present (or not).
//!
//! ## Bounding Boxes
//!
//! Many parts of the system operate on axis-aligned regions of
//! space. For example, multi-threading in **pbrt** is implemented by
//! subdividing the image into rectangular tiles that can be processed
//! independently, and the bounding volume hierarchy uses 3D boxes to
//! bound geometric primitives in the scene. The **Bounds2** and
//! **Bounds3** template classes are used to represent the extent of
//! these sort of regions. Both are parameterized by a type T that is
//! used to represent the coordinates of its extents.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::{Bounds3, Point3};
//!
//! fn main() {
//!     let int_origin = Point3 { x: 0, y: 0, z: 0 };
//!     let int_xyz111 = Point3 { x: 1, y: 1, z: 1 };
//!     let float_origin = Point3 {
//!         x: 0.0,
//!         y: 0.0,
//!         z: 0.0,
//!     };
//!     let float_xyz111 = Point3 {
//!         x: 1.0,
//!         y: 1.0,
//!         z: 1.0,
//!     };
//!     let int_unit_cube = Bounds3 {
//!         p_min: int_origin,
//!         p_max: int_xyz111,
//!     };
//!     let float_unit_cube = Bounds3 {
//!         p_min: float_origin,
//!         p_max: float_xyz111,
//!     };
//!
//!     println!("int   {:?}", int_unit_cube);
//!     println!("float {:?}", float_unit_cube);
//! }
//! ```
//!
//! ## 4 x 4 Matrices
//!
//! The **Matrix4x4** structure provides a low-level representation of
//! 4 x 4 matrices. It is an integral part of the **Transform** class.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::Matrix4x4;
//!
//! fn main() {
//!     let identity: Matrix4x4 = Matrix4x4::default();
//!     let m: Matrix4x4 = Matrix4x4::new(0.0,
//!                                       0.1,
//!                                       0.2,
//!                                       0.3,
//!                                       1.0,
//!                                       1.1,
//!                                       1.2,
//!                                       1.3,
//!                                       2.0,
//!                                       2.1,
//!                                       2.2,
//!                                       2.3,
//!                                       3.0,
//!                                       3.1,
//!                                       3.2,
//!                                       3.3);
//!
//!     println!("identity matrix = {:?}", identity);
//!     println!("m = {:?}", m);
//! }
//! ```
//!
//! ## Transformations
//!
//! In general a transformation is a mapping from points to points and
//! from vectors to vectors. When a new **Transform** is created, it
//! defaults to the *identity transformation* - the transformation
//! that maps each point and each vector to itself.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::Transform;
//!
//! fn main() {
//!     let identity: Transform = Transform::default();
//!     let t: Transform = Transform::new(2.0,
//!                                       0.0,
//!                                       0.0,
//!                                       -1.25,
//!                                       0.0,
//!                                       -8.0,
//!                                       0.0,
//!                                       3.5,
//!                                       0.0,
//!                                       0.0,
//!                                       1.0,
//!                                       7.875,
//!                                       0.0,
//!                                       0.0,
//!                                       0.0,
//!                                       1.0);
//!
//!     println!("identity matrix = {:?}", identity);
//!     println!("t = {:?}", t);
//! }
//! ```
//!
//! ### Translations
//!
//! One of the simplest transformations is the translation
//! transformation. Translations only affect points, leaving vectors
//! unchanged.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::{Transform, Vector3f};
//!
//! fn main() {
//!     let t: Transform = Transform::translate(Vector3f {
//!         x: -1.25,
//!         y: 3.5,
//!         z: 7.875,
//!     });
//!
//!     println!("t = {:?}", t);
//! }
//! ```
//!
//! ### Scaling
//!
//! Another basic transformations is the scale transformation. We can
//! differentiate between **uniform** scaling, where all three scale
//! factors have the same value, and **nonuniform** scaling, where
//! they may have different values.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::Transform;
//!
//! fn main() {
//!     let t: Transform = Transform::scale(2.0, -8.0, 1.0);
//!
//!     println!("t = {:?}", t);
//! }
//! ```
//!
//! ### X, Y, And Z Axis Rotations
//!
//! Another useful type of transformation is the rotation
//! transformation.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::{Float, Transform};
//!
//! fn main() {
//!     let theta_x: Float = 30.0;
//!     let theta_y: Float = 45.0;
//!     let theta_z: Float = 60.0;
//!     let t_x: Transform = Transform::rotate_x(theta_x);
//!     let t_y: Transform = Transform::rotate_y(theta_y);
//!     let t_z: Transform = Transform::rotate_z(theta_z);
//!
//!     println!("Transform::rotate_x({}) = {:?}", theta_x, t_x);
//!     println!("Transform::rotate_y({}) = {:?}", theta_y, t_y);
//!     println!("Transform::rotate_z({}) = {:?}", theta_z, t_z);
//! }
//! ```
//!
//! ### Rotation Around an Arbitrary Axis
//!
//! We also provide a routine to compute the transformation that
//! represents rotation around an arbitrary axis.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::{Float, Transform, Vector3f};
//!
//! fn main() {
//!     let theta: Float = 30.0;
//!     let axis = Vector3f {
//!         x: 1.0,
//!         y: 2.0,
//!         z: 3.0,
//!     };
//!     let t: Transform = Transform::rotate(theta, axis);
//!
//!     println!("Transform::rotate({}, {:?}) = {:?}", theta, axis, t);
//! }
//! ```
//!
//! ### The Look-At Transformation
//!
//! The *look-at* transformation is particularly useful for placing a
//! camera in the scene. The caller specifies the desired position of
//! the camera, a point the camera is looking at, and an "up" vector
//! that orients the camera along the viewing direction implied by the
//! first two parameters. All of these values are given in world space
//! coordinates. The look-at construction then gives a transformation
//! between camera space and world space.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::{Point3f, Transform, Vector3f};
//!
//! fn main() {
//!     // LookAt 2 2 5  0 -.4 0  0 1 0 (see spheres-differentials-texfilt.pbrt)
//!     let pos = Point3f {
//!         x: 2.0,
//!         y: 2.0,
//!         z: 5.0,
//!     };
//!     let look = Point3f {
//!         x: 0.0,
//!         y: -0.4,
//!         z: 0.0,
//!     };
//!     let up = Vector3f {
//!         x: 0.0,
//!         y: 1.0,
//!         z: 0.0,
//!     };
//!     let t: Transform = Transform::look_at(pos, look, up);
//!
//!     println!("Transform::look_at({:?}, {:?}, {:?}) = {:?}",
//!              pos,
//!              look,
//!              up,
//!              t);
//! }
//! ```
//!
//! ## Quaternions
//!
//! Quaternions were originally invented by Sir William Hamilton in
//! 1843 as a generalization of complex numbers (2 dimensions) to four
//! dimensions.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::Quaternion;
//!
//! fn main() {
//!     let default_quaternion: Quaternion = Quaternion::default();
//!
//!     println!("default quaternion = {:?}", default_quaternion);
//! }
//! ```
//!
//! ## Surface Interaction
//!
//! The geometry of a particular point on a surface is represented by
//! a **SurfaceInteraction**. Having this abstraction lets most of the
//! system work with points on surfaces without needing to consider
//! the particular type of geometric shape the points lie on; the
//! **SurfaceInteraction** abstraction supplies enough information
//! about the surface point to allow the shading and geometric
//! operations in the rest of **pbrt** to be implemented generically.
//!
//! ## Shapes
//!
//! Careful abstraction of geometric shapes in a ray tracer is a key
//! component of a clean system design, and shapes are the ideal
//! candidate for an object-oriented approach. All geometric
//! primitives implement a common interface, and the rest of the
//! renderer can use this interface without needing any details about
//! the underlying shape. This makes it possible to separate the
//! geometric and the shading subsystem of pbrt.
//!
//! ### Spheres
//!
//! Spheres are a special case of a general type of surfaces called
//! quadrics. They are the simplest type of curved surfaces that is
//! useful to a ray tracer and are a good starting point for general
//! ray intersection routines.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::{Float, Sphere, Transform, Vector3f};
//!
//! fn main() {
//!     let translate: Transform = Transform::translate(Vector3f {
//!                                                         x: -1.3,
//!                                                         y: 0.0,
//!                                                         z: 0.0,
//!                                                     });
//!     let inverse: Transform = Transform::inverse(translate);
//!     let radius: Float = 1.0;
//!     let z_min: Float = -1.0;
//!     let z_max: Float = 1.0;
//!     let phi_max: Float = 360.0;
//!     let sphere = Sphere::new(translate,
//!                              inverse,
//!                              false,
//!                              false,
//!                              radius,
//!                              z_min,
//!                              z_max,
//!                              phi_max);
//!     println!("translate = {:?}", translate);
//!     println!("inverse = {:?}", inverse);
//!     println!("sphere.radius = {:?}", sphere.radius);
//!     println!("sphere.z_min = {:?}", sphere.z_min);
//!     println!("sphere.z_max = {:?}", sphere.z_max);
//!     println!("sphere.theta_min = {:?}", sphere.theta_min);
//!     println!("sphere.theta_max = {:?}", sphere.theta_max);
//!     println!("sphere.phi_max = {:?}", sphere.phi_max);
//! }
//! ```
//!
//! ### Triangle Meshes
//!
//! While a natural representation would be to have a **Triangle**
//! shape implementation where each triangle stored the positions of
//! its three vertices, a more memory-efficient representation is to
//! separately store entire triangle meshes with an array of vertex
//! positions where each individual triangle just stores three offsets
//! into this array for its three vertices.
//!
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::{Normal3f, Point2f, Point3f, Transform, Triangle, TriangleMesh, Vector3f};
//! use std::sync::Arc;
//!
//! fn main() {
//!     let translate: Transform = Transform::translate(Vector3f {
//!         x: 0.25,
//!         y: 0.0,
//!         z: 0.0,
//!     });
//!     let inverse: Transform = Transform::inverse(translate);
//!     let n_triangles: usize = 2;
//!     let vertex_indices: Vec<usize> = vec![0_usize, 2, 1, 0, 3, 2];
//!     let n_vertices: usize = 4;
//!     let p: Vec<Point3f> = vec![Point3f {
//!                                    x: -100.0,
//!                                    y: -1.0,
//!                                    z: -100.0,
//!                                },
//!                                Point3f {
//!                                    x: 400.0,
//!                                    y: -1.0,
//!                                    z: -100.0,
//!                                },
//!                                Point3f {
//!                                    x: 400.0,
//!                                    y: -1.0,
//!                                    z: 400.0,
//!                                },
//!                                Point3f {
//!                                    x: -100.0,
//!                                    y: -1.0,
//!                                    z: 400.0,
//!                                }];
//!     let s: Vec<Vector3f> = Vec::new();
//!     let n: Vec<Normal3f> = Vec::new();
//!     let uv: Vec<Point2f> = vec![Point2f { x: 0.0, y: 0.0 },
//!                                 Point2f { x: 1.0, y: 0.0 },
//!                                 Point2f { x: 0.0, y: 1.0 },
//!                                 Point2f { x: 1.0, y: 1.0 }];
//!     let triangle_mesh: TriangleMesh = TriangleMesh::new(translate,
//!                                                         inverse,
//!                                                         false,
//!                                                         false,
//!                                                         n_triangles,
//!                                                         vertex_indices,
//!                                                         n_vertices,
//!                                                         p,
//!                                                         s,
//!                                                         n,
//!                                                         uv);
//!     println!("translate = {:?}", translate);
//!     println!("inverse = {:?}", inverse);
//!     println!("triangle_mesh = {:?}", triangle_mesh);
//!     for id in 0..triangle_mesh.n_triangles {
//!         let triangle = Triangle::new(triangle_mesh.object_to_world,
//!                                      triangle_mesh.world_to_object,
//!                                      triangle_mesh.transform_swaps_handedness,
//!                                      Arc::new(triangle_mesh.clone()),
//!                                      id);
//!         println!("triangle.id = {:?}", triangle.id);
//!     }
//! }
//! ```
//!
//! ### Cones
//!
//! TODO
//!
//! ### Disks
//!
//! The disk is an interesting quadric since it has a particularly
//! straightforward intersection routine that avoids solving the
//! quadric equation.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::{Disk, Float, Transform, Vector3f};
//!
//! fn main() {
//!     let translate: Transform = Transform::translate(Vector3f {
//!                                                         x: 0.0,
//!                                                         y: 0.0,
//!                                                         z: -0.01,
//!                                                     });
//!     let inverse: Transform = Transform::inverse(translate);
//!     let height: Float = 0.0;
//!     let radius: Float = 30.0;
//!     let inner_radius: Float = 0.0;
//!     let phi_max: Float = 360.0;
//!     let disk = Disk::new(translate,
//!                          inverse,
//!                          false,
//!                          false,
//!                          height,
//!                          radius,
//!                          inner_radius,
//!                          phi_max);
//!     println!("translate = {:?}", translate);
//!     println!("inverse = {:?}", inverse);
//!     println!("disk.height = {:?}", disk.height);
//!     println!("disk.radius = {:?}", disk.radius);
//!     println!("disk.inner_radius = {:?}", disk.inner_radius);
//!     println!("disk.phi_max = {:?}", disk.phi_max);
//! }
//! ```
//!
//! ### Cylinders
//!
//! Another useful quadric is the cylinder. Cylinder shapes are
//! centered around the z axis.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::{Cylinder, Float, Transform, Vector3f};
//!
//! fn main() {
//!     let translate: Transform = Transform::translate(Vector3f {
//!                                                         x: 2.0,
//!                                                         y: 1.0,
//!                                                         z: 1.6,
//!                                                     });
//!     let inverse: Transform = Transform::inverse(translate);
//!     let radius: Float = 1.0;
//!     let z_min: Float = 0.0;
//!     let z_max: Float = 1.0;
//!     let phi_max: Float = 360.0;
//!     let cylinder = Cylinder::new(translate,
//!                                  inverse,
//!                                  false,
//!                                  radius,
//!                                  z_min,
//!                                  z_max,
//!                                  phi_max);
//!     println!("translate = {:?}", translate);
//!     println!("inverse = {:?}", inverse);
//!     println!("cylinder.radius = {:?}", cylinder.radius);
//!     println!("cylinder.z_min = {:?}", cylinder.z_min);
//!     println!("cylinder.z_max = {:?}", cylinder.z_max);
//!     println!("cylinder.phi_max = {:?}", cylinder.phi_max);
//! }
//! ```
//!
//! ### Hyperboloids
//!
//! TODO
//!
//! ### Paraboloids
//!
//! TODO
//!
//! ## Filters
//!
//! ### Box Filter
//!
//! One of the most commonly used filters in graphics is the box
//! filter. The box filter equally weights all samples within a square
//! region of the image. Although computational efficient, it's just
//! about the worst filter possible.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::{BoxFilter, Float, Vector2f};
//!
//! fn main() {
//!     let xw: Float = 0.5;
//!     let yw: Float = 0.5;
//!     let box_filter = BoxFilter {
//!         radius: Vector2f { x: xw, y: yw },
//!         inv_radius: Vector2f {
//!             x: 1.0 / xw,
//!             y: 1.0 / yw,
//!         },
//!     };
//!
//!     println!("box_filter = {:?}", box_filter);
//! }
//! ```
//! ### TriangleFilter
//!
//! TODO
//!
//! ### Gaussian Filter
//!
//! Unlike the box and triangle filters, the Gaussian filter gives a
//! reasonably good result in practice. The Gaussian filter does tend
//! to cause slight blurring of the final image compared to some of
//! the other filters, but this blurring can actually help mask any
//! remaining aliasing in the image.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::{GaussianFilter, Float, Vector2f};
//!
//! fn main() {
//!     let xw: Float = 2.0;
//!     let yw: Float = 2.0;
//!     let alpha: Float = 2.0;
//!     let exp_x: Float = (-alpha * xw * xw).exp();
//!     let exp_y: Float = (-alpha * yw * yw).exp();
//!     let gaussian_filter = GaussianFilter {
//!         alpha: alpha,
//!         exp_x: exp_x,
//!         exp_y: exp_y,
//!         radius: Vector2f { x: xw, y: yw },
//!         inv_radius: Vector2f {
//!             x: 1.0 / xw,
//!             y: 1.0 / yw,
//!         },
//!     };
//!
//!     println!("gaussian_filter = {:?}", gaussian_filter);
//! }
//! ```
//!
//! ### MitchellFilter
//!
//! TODO
//!
//! ### LanczosSincFilter
//!
//! TODO
//!
//! ## Film
//!
//! The type of film or sensor in a camera has a dramatic effect on
//! the way that incident light is transformed into colors in an
//! image. In **pbrt**, the **Film** class models the sensing device
//! in the simulated camera. After the radiance is found for each
//! camera ray, the **Film** implementation determines the sample's
//! contribution to the pixel around the point on the film plane where
//! the camera ray began and updates its representation of the
//! image. When the main rendering loop exits, the **Film** writes the
//! final image to file.
//!
//! ## Light Sources
//!
//! ### Point Lights
//!
//! Isotropic point light source that emits the same amount of light
//! in all directions.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::{PointLight, Spectrum, Transform};
//!
//! fn main() {
//!     let i: Spectrum = Spectrum::new(50.0);
//!     let light_to_world: Transform = Transform::default();
//!     let point_light: PointLight = PointLight::new(&light_to_world, &i);
//!     println!("point_light = {:?}", point_light);
//! }
//! ```
//!
//! #### Spotlights
//!
//! TODO
//!
//! #### Texture Projection Lights
//!
//! TODO
//!
//! #### Goniophotometric Diagram Lights
//!
//! TODO
//!
//! ### Distant Lights
//!
//! A distant light, also known as directional light, describes an
//! emitter that deposits illumination from the same direction at
//! every point in space.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::{DistantLight, Point3f, Spectrum, Transform, Vector3f};
//!
//! fn main() {
//!     let l: Spectrum = Spectrum::new(3.141593);
//!     let sc: Spectrum = Spectrum::new(1.0);
//!     let from: Point3f = Point3f {
//!         x: 0.0,
//!         y: 10.0,
//!         z: 0.0,
//!     };
//!     let to: Point3f = Point3f {
//!         x: 0.0,
//!         y: 0.0,
//!         z: 0.0,
//!     };
//!     let dir: Vector3f = from - to;
//!     let light_to_world: Transform = Transform::default();
//!     let lsc: Spectrum = l * sc;
//!     let distant_light: DistantLight = DistantLight::new(&light_to_world, &lsc, &dir);
//!     println!("distant_light = {:?}", distant_light);
//! }
//! ```
//!
//! ### Area Lights
//!
//! Area lights are light sources defined by one or more **Shapes**
//! that emit light from their surface, with some directional
//! distribution of radiance at each point on the surface.
//!
//! ### Infinite Area Lights
//!
//! TODO
//!
//! #### Diffuse Area Lights
//!
//! **DiffuseAreaLight** implements a basic area light source with a
//! uniform spatial and directional radiance distribution. The surface
//! it emits from is defined by a **Shape**. It only emits light on
//! the side of the surface with outward-facing surface normal; there
//! is no emission from the other side.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::{DiffuseAreaLight, Disk, Float, Shape, Spectrum, Transform, Vector3f};
//! use std::sync::Arc;
//!
//! fn main() {
//!     let t: Transform = Transform::translate(Vector3f {
//!                                                 x: 2.0,
//!                                                 y: -4.0,
//!                                                 z: 4.0,
//!                                             });
//!     let theta: Float = -120.0;
//!     let axis = Vector3f {
//!         x: 1.0,
//!         y: 0.0,
//!         z: 0.0,
//!     };
//!     let r: Transform = Transform::rotate(theta, axis);
//!     let light_to_world: Transform = Transform::default() * r * t;
//!     let inverse: Transform = Transform::inverse(light_to_world);
//!     let l_emit: Spectrum = Spectrum::new(8.0);
//!     let n_samples: i32 = 16;
//!     let height: Float = 0.0;
//!     let radius: Float = 2.0;
//!     let inner_radius: Float = 0.0;
//!     let phi_max: Float = 360.0;
//!     let shape: Arc<Shape + Send + Sync> = Arc::new(Disk::new(light_to_world,
//!                                                              inverse,
//!                                                              false,
//!                                                              false,
//!                                                              height,
//!                                                              radius,
//!                                                              inner_radius,
//!                                                              phi_max));
//!     let two_sided: bool = false;
//!     let diffuse_area_light: DiffuseAreaLight =
//!         DiffuseAreaLight::new(&light_to_world, &l_emit, n_samples, shape, two_sided);
//!     println!("diffuse_area_light.l_emit = {:?}", diffuse_area_light.l_emit);
//!     println!("diffuse_area_light.two_sided = {:?}", diffuse_area_light.two_sided);
//!     println!("diffuse_area_light.area = {:?}", diffuse_area_light.area);
//! }
//! ```
//!
//! ## Direct Lighting
//!
//! The **DirectLightingIntegrator** accounts only for direct lighting
//! &mdash; light that has traveled directly from a light source to the
//! point being shaded &mdash; and ignores indirect illumination from
//! objects that are not themselfes emissive, except for basic
//! specular reflection and transmission effects.
#![feature(integer_atomics)]

extern crate crossbeam;
extern crate half;
extern crate image;
extern crate num;
extern crate num_cpus;
extern crate openexr;
extern crate pbr;
extern crate time;
extern crate typed_arena;

// use std::cell::RefCell;
use std::borrow::Borrow;
use std::cmp::PartialEq;
use std::collections::HashMap;
use std::default::Default;
use std::f32::consts::PI;
use std::mem;
use std::ops::{BitAnd, Add, AddAssign, Sub, Mul, MulAssign, Div, DivAssign, Neg, Index, IndexMut};
use std::path::Path;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, RwLock};
use std::sync::mpsc;
// use copy_arena::{Arena, Allocator};
use typed_arena::Arena;
use num::Zero;
use image::{ImageResult, DynamicImage};
use openexr::{FrameBuffer, FrameBufferMut, Header, InputFile, PixelType, ScanlineOutputFile};
use half::f16;
use time::PreciseTime;

pub type Float = f32;

// see https://stackoverflow.com/questions/36008434/how-can-i-decode-f16-to-f32-using-only-the-stable-standard-library
#[inline]
fn decode_f16(half: u16) -> f32 {
    let exp: u16 = half >> 10 & 0x1f;
    let mant: u16 = half & 0x3ff;
    let val: f32 = if exp == 0 {
        (mant as f32) * (2.0f32).powi(-24)
    } else if exp != 31 {
        (mant as f32 + 1024f32) * (2.0f32).powi(exp as i32 - 25)
    } else if mant == 0 {
        ::std::f32::INFINITY
    } else {
        ::std::f32::NAN
    };
    if half & 0x8000 != 0 {
        -val
    } else {
        val
    }
}

// see scene.h

#[derive(Clone)]
pub struct Scene {
    pub lights: Vec<Arc<Light + Sync + Send>>,
    pub infinite_lights: Vec<Arc<Light + Sync + Send>>,
    aggregate: Arc<BVHAccel>, // TODO: Primitive,
    world_bound: Bounds3f,
}

impl Scene {
    pub fn new(aggregate: Arc<BVHAccel>, lights: Vec<Arc<Light + Sync + Send>>) -> Self {
        let world_bound: Bounds3f = aggregate.world_bound();
        let scene: Scene = Scene {
            lights: Vec::new(),
            infinite_lights: Vec::new(),
            aggregate: aggregate.clone(),
            world_bound: world_bound,
        };
        let mut changed_lights = Vec::new();
        let mut infinite_lights = Vec::new();
        for light in lights {
            light.preprocess(&scene);
            changed_lights.push(light.clone());
            let check: u8 = light.get_flags() & LightFlags::Infinite as u8;
            if check == LightFlags::Infinite as u8 {
                infinite_lights.push(light);
            }
        }
        Scene {
            lights: changed_lights,
            infinite_lights: infinite_lights,
            aggregate: aggregate,
            world_bound: world_bound,
        }
    }
    pub fn world_bound(&self) -> Bounds3f {
        self.world_bound
    }
    pub fn intersect(&self, ray: &mut Ray) -> Option<SurfaceInteraction> {
        // TODO: ++nIntersectionTests;
        assert_ne!(ray.d,
                   Vector3f {
                       x: 0.0,
                       y: 0.0,
                       z: 0.0,
                   });
        self.aggregate.intersect(ray)
    }
    pub fn intersect_p(&self, ray: &mut Ray) -> bool {
        // TODO: ++nShadowTests;
        assert_ne!(ray.d,
                   Vector3f {
                       x: 0.0,
                       y: 0.0,
                       z: 0.0,
                   });
        self.aggregate.intersect_p(ray)
    }
}

// see pbrt.h

pub type Spectrum = RGBSpectrum;

const MACHINE_EPSILON: Float = std::f32::EPSILON * 0.5;
const SHADOW_EPSILON: Float = 0.0001;
const INV_PI: Float = 0.31830988618379067154;
const INV_2_PI: Float = 0.15915494309189533577;
const PI_OVER_2: Float = 1.57079632679489661923;
const PI_OVER_4: Float = 0.78539816339744830961;

/// Use **unsafe**
/// [std::mem::transmute_copy][transmute_copy]
/// to convert *f32* to *u32*.
///
/// [transmute_copy]: https://doc.rust-lang.org/std/mem/fn.transmute_copy.html
pub fn float_to_bits(f: f32) -> u32 {
    // uint64_t ui;
    // memcpy(&ui, &f, sizeof(double));
    // return ui;
    let rui: u32;
    unsafe {
        let ui: u32 = std::mem::transmute_copy(&f);
        rui = ui;
    }
    rui
}

/// Use **unsafe**
/// [std::mem::transmute_copy][transmute_copy]
/// to convert *u32* to *f32*.
///
/// [transmute_copy]: https://doc.rust-lang.org/std/mem/fn.transmute_copy.html
pub fn bits_to_float(ui: u32) -> f32 {
    // float f;
    // memcpy(&f, &ui, sizeof(uint32_t));
    // return f;
    let rf: f32;
    unsafe {
        let f: f32 = std::mem::transmute_copy(&ui);
        rf = f;
    }
    rf
}

/// Bump a floating-point value up to the next greater representable
/// floating-point value.
pub fn next_float_up(v: f32) -> f32 {
    if v.is_infinite() && v > 0.0 {
        v
    } else {
        let new_v: f32;
        if v == -0.0 {
            new_v = 0.0;
        } else {
            new_v = v;
        }
        let mut ui: u32 = float_to_bits(new_v);
        if new_v >= 0.0 {
            ui += 1;
        } else {
            ui -= 1;
        }
        bits_to_float(ui)
    }
}

/// Bump a floating-point value down to the next smaller representable
/// floating-point value.
pub fn next_float_down(v: f32) -> f32 {
    if v.is_infinite() && v < 0.0 {
        v
    } else {
        let new_v: f32;
        if v == 0.0 {
            new_v = -0.0;
        } else {
            new_v = v;
        }
        let mut ui: u32 = float_to_bits(new_v);
        if new_v > 0.0 {
            ui -= 1;
        } else {
            ui += 1;
        }
        bits_to_float(ui)
    }
}

/// Error propagation.
pub fn gamma(n: i32) -> Float {
    (n as Float * MACHINE_EPSILON) / (1.0 - n as Float * MACHINE_EPSILON)
}

/// Is used to write sRGB-compatible 8-bit image files.
pub fn gamma_correct(value: Float) -> Float {
    if value <= 0.0031308 {
        12.92 * value
    } else {
        1.055 as Float * value.powf((1.0 / 2.4) as Float) - 0.055
    }
}

/// Clamp the given value *val* to lie between the values *low* and *high*.
pub fn clamp_t<T>(val: T, low: T, high: T) -> T
    where T: PartialOrd
{
    let r: T;
    if val < low {
        r = low;
    } else if val > high {
        r = high;
    } else {
        r = val;
    }
    r
}

/// Computes the remainder of a/b. Provides the behavior that the
/// modulus of a negative number is always positive.
pub fn mod_t<T>(a: T, b: T) -> T
    where T: num::Zero + Copy + PartialOrd +
    Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + Div<T, Output = T>
{
    let result: T = a - (a / b) * b;
    if result < num::Zero::zero() {
        result + b
    } else {
        result
    }
}

/// Convert from angles expressed in degrees to radians.
pub fn radians(deg: Float) -> Float {
    (PI / 180.0) * deg
}

/// Convert from angles expressed in radians to degrees.
pub fn degrees(rad: Float) -> Float {
    (180.0 / PI) * rad
}

/// Determine if a given integer is an exact power of 2.
pub fn is_power_of_2<T>(v: T) -> bool
    where T: num::Zero + num::One + Copy + PartialOrd + BitAnd<T, Output = T> + Sub<T, Output = T>
{
    // https://doc.rust-lang.org/std/primitive.u32.html#method.is_power_of_two
    (v > num::Zero::zero()) && !((v & (v - num::One::one())) > num::Zero::zero())
}

/// Round an integer up to the next higher (or equal) power of 2.
pub fn round_up_pow2_32(v: i32) -> i32 {
    let mut ret: i32 = v; // copy value
    ret -= 1_i32;
    ret |= ret >> 1;
    ret |= ret >> 2;
    ret |= ret >> 4;
    ret |= ret >> 8;
    ret |= ret >> 16;
    ret + 1
}

/// Round an integer up to the next higher (or equal) power of 2.
pub fn round_up_pow2_64(v: i64) -> i64 {
    let mut ret: i64 = v; // copy value
    ret -= 1_i64;
    ret |= ret >> 1;
    ret |= ret >> 2;
    ret |= ret >> 4;
    ret |= ret >> 8;
    ret |= ret >> 16;
    ret + 1
}

/// Interpolate linearly between two provided values.
pub fn lerp(t: Float, v1: Float, v2: Float) -> Float {
    (1.0 as Float - t) * v1 + t * v2
}

/// Find solution(s) of the quadratic equation at<sup>2</sup> + bt + c = 0.
pub fn quadratic(a: Float, b: Float, c: Float, t0: &mut Float, t1: &mut Float) -> bool {
    // find quadratic discriminant
    let discrim: f64 = (b as f64) * (b as f64) - 4.0 * (a as f64) * (c as f64);
    if discrim < 0.0 {
        false
    } else {
        let root_discrim: f64 = discrim.sqrt();
        // compute quadratic _t_ values
        let q: f64;
        if b < 0.0 {
            q = -0.5 * (b as f64 - root_discrim);
        } else {
            q = -0.5 * (b as f64 + root_discrim);
        }
        *t0 = q as Float / a;
        *t1 = c / q as Float;
        if *t0 > *t1 {
            // std::swap(*t0, *t1);
            let swap = *t0;
            *t0 = *t1;
            *t1 = swap;
        }
        true
    }
}

// see efloat.h

/// Find solution(s) of the quadratic equation at<sup>2</sup> + bt + c = 0 using
/// *EFloat* instead of *Float* for error bounds.
pub fn quadratic_efloat(a: EFloat, b: EFloat, c: EFloat, t0: &mut EFloat, t1: &mut EFloat) -> bool {
    let discrim: f64 = b.v as f64 * b.v as f64 - 4.0f64 * a.v as f64 * c.v as f64;
    if discrim < 0.0 {
        false
    } else {
        let root_discrim: f64 = discrim.sqrt();
        let float_root_discrim: EFloat = EFloat::new(root_discrim as f32,
                                                     MACHINE_EPSILON as f32 * root_discrim as f32);
        // compute quadratic _t_ values
        let q: EFloat;
        if b.v < 0.0f32 {
            q = (b - float_root_discrim) * -0.5f32;
        } else {
            q = (b + float_root_discrim) * -0.5f32;
        }
        *t0 = q / a;
        *t1 = c / q;
        if (*t0).v > (*t1).v {
            let swap: EFloat = *t0;
            *t0 = *t1;
            *t1 = swap;
        }
        true
    }
}

#[derive(Debug,Default,Copy,Clone)]
pub struct EFloat {
    v: f32,
    low: f32,
    high: f32,
}

impl EFloat {
    pub fn new(v: f32, err: f32) -> Self {
        if err == 0.0 {
            EFloat {
                v: v,
                low: v,
                high: v,
            }
        } else {
            EFloat {
                v: v,
                low: next_float_down(v - err),
                high: next_float_up(v + err),
            }
        }
    }
    pub fn lower_bound(&self) -> f32 {
        self.low
    }
    pub fn upper_bound(&self) -> f32 {
        self.high
    }
}

impl PartialEq for EFloat {
    fn eq(&self, rhs: &EFloat) -> bool {
        self.v == rhs.v
    }
}

impl Add for EFloat {
    type Output = EFloat;
    fn add(self, rhs: EFloat) -> EFloat {
        // TODO: r.Check();
        EFloat {
            v: self.v + rhs.v,
            low: next_float_down(self.lower_bound() + rhs.lower_bound()),
            high: next_float_up(self.upper_bound() + rhs.upper_bound()),
        }
    }
}

impl Sub for EFloat {
    type Output = EFloat;
    fn sub(self, rhs: EFloat) -> EFloat {
        // TODO: r.Check();
        EFloat {
            v: self.v - rhs.v,
            low: next_float_down(self.lower_bound() - rhs.upper_bound()),
            high: next_float_up(self.upper_bound() - rhs.lower_bound()),
        }
    }
}

impl Mul for EFloat {
    type Output = EFloat;
    fn mul(self, rhs: EFloat) -> EFloat {
        let prod: [f32; 4] = [self.lower_bound() * rhs.lower_bound(),
                              self.upper_bound() * rhs.lower_bound(),
                              self.lower_bound() * rhs.upper_bound(),
                              self.upper_bound() * rhs.upper_bound()];
        // TODO: r.Check();
        EFloat {
            v: self.v * rhs.v,
            low: next_float_down(prod[0].min(prod[1]).min(prod[2].min(prod[3]))),
            high: next_float_up(prod[0].max(prod[1]).max(prod[2].max(prod[3]))),
        }
    }
}

impl Mul<f32> for EFloat {
    type Output = EFloat;
    fn mul(self, rhs: f32) -> EFloat {
        EFloat::new(rhs, 0.0) * self
    }
}

impl Div for EFloat {
    type Output = EFloat;
    fn div(self, rhs: EFloat) -> EFloat {
        let div: [f32; 4] = [self.lower_bound() / rhs.lower_bound(),
                             self.upper_bound() / rhs.lower_bound(),
                             self.lower_bound() / rhs.upper_bound(),
                             self.upper_bound() / rhs.upper_bound()];
        // TODO: r.Check();
        if rhs.low < 0.0 && rhs.high > 0.0 {
            // the interval we're dividing by straddles zero, so just
            // return an interval of everything
            EFloat {
                v: self.v / rhs.v,
                low: -std::f32::INFINITY,
                high: std::f32::INFINITY,
            }
        } else {
            EFloat {
                v: self.v / rhs.v,
                low: next_float_down(div[0].min(div[1]).min(div[2].min(div[3]))),
                high: next_float_up(div[0].max(div[1]).max(div[2].max(div[3]))),
            }
        }
    }
}

// parallel.h

#[derive(Debug,Default,Copy,Clone)]
pub struct AtomicFloat {
    bits: u32,
}

impl From<AtomicFloat> for Float {
    fn from(a: AtomicFloat) -> Float {
        bits_to_float(a.bits) as Float
    }
}

// see geometry.h

pub type Point2f = Point2<Float>;
pub type Point2i = Point2<i32>;
pub type Point3f = Point3<Float>;
pub type Point3i = Point3<i32>;
pub type Vector2f = Vector2<Float>;
pub type Vector2i = Vector2<i32>;
pub type Vector3f = Vector3<Float>;
pub type Vector3i = Vector3<i32>;
pub type Normal3f = Normal3<Float>;

#[derive(Debug,Default,Copy,Clone)]
pub struct Vector2<T> {
    pub x: T,
    pub y: T,
}

impl<T> Vector2<T> {
    pub fn length_squared(&self) -> T
        where T: Copy + Add<T, Output = T> + Mul<T, Output = T>
    {
        self.x * self.x + self.y * self.y
    }
    pub fn length(&self) -> T
        where T: num::Float
    {
        self.length_squared().sqrt()
    }
}

impl<T> MulAssign<T> for Vector2<T>
    where T: Copy + MulAssign
{
    fn mul_assign(&mut self, rhs: T) {
        self.x *= rhs;
        self.y *= rhs;
    }
}

/// Product of the Euclidean magnitudes of the two vectors and the
/// cosine of the angle between them. A return value of zero means
/// both vectors are orthogonal, a value if one means they are
/// codirectional.
pub fn vec2_dot<T>(v1: Vector2<T>, v2: Vector2<T>) -> T
    where T: Copy + Add<T, Output = T> + Mul<T, Output = T>
{
    v1.x * v2.x + v1.y * v2.y
}

#[derive(Debug,Default,Copy,Clone,PartialEq)]
pub struct Vector3<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T> Add for Vector3<T>
    where T: Copy + Add<T, Output = T>
{
    type Output = Vector3<T>;
    fn add(self, rhs: Vector3<T>) -> Vector3<T> {
        Vector3::<T> {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl<T> Sub for Vector3<T>
    where T: Copy + Sub<T, Output = T>
{
    type Output = Vector3<T>;
    fn sub(self, rhs: Vector3<T>) -> Vector3<T> {
        Vector3::<T> {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl<T> Mul<T> for Vector3<T>
    where T: Copy + Mul<T, Output = T>
{
    type Output = Vector3<T>;
    fn mul(self, rhs: T) -> Vector3<T>
        where T: Copy + Mul<T, Output = T>
    {
        Vector3::<T> {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

// work around bug
// https://github.com/rust-lang/rust/issues/40395
impl Div<Float> for Vector3<f32>
{
    type Output = Vector3<f32>;
    fn div(self, rhs: Float) -> Vector3<f32>
    {
        assert_ne!(rhs, 0.0 as Float);
        let inv: Float = 1.0 as Float / rhs;
        Vector3::<f32> {
            x: self.x * inv,
            y: self.y * inv,
            z: self.z * inv,
        }
    }
}

impl<T> Neg for Vector3<T>
    where T: Copy + Neg<Output = T>
{
    type Output = Vector3<T>;
    fn neg(self) -> Vector3<T> {
        Vector3::<T> {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl<T> Vector3<T> {
    pub fn abs(&self) -> Vector3<T>
        where T: num::Float
    {
        Vector3::<T> {
            x: self.x.abs(),
            y: self.y.abs(),
            z: self.z.abs(),
        }
    }
    pub fn length_squared(&self) -> T
        where T: Copy + Add<T, Output = T> + Mul<T, Output = T>
    {
        self.x * self.x + self.y * self.y + self.z * self.z
    }
    pub fn length(&self) -> T
        where T: num::Float
    {
        self.length_squared().sqrt()
    }
}

impl<T> Index<u8> for Vector3<T> {
    type Output = T;
    fn index(&self, index: u8) -> &T {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Check failed: i >= 0 && i <= 2"),
        }
    }
}

impl<T> IndexMut<u8> for Vector3<T> {
    fn index_mut(&mut self, index: u8) -> &mut T {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => panic!("Check failed: i >= 0 && i <= 2"),
        }
    }
}

impl<T> From<Point3<T>> for Vector3<T> {
    fn from(p: Point3<T>) -> Self {
        Vector3::<T> {
            x: p.x,
            y: p.y,
            z: p.z,
        }
    }
}

impl<T> From<Normal3<T>> for Vector3<T> {
    fn from(n: Normal3<T>) -> Self {
        Vector3::<T> {
            x: n.x,
            y: n.y,
            z: n.z,
        }
    }
}

/// Product of the Euclidean magnitudes of the two vectors and the
/// cosine of the angle between them. A return value of zero means
/// both vectors are orthogonal, a value if one means they are
/// codirectional.
pub fn vec3_dot_vec3<T>(v1: Vector3<T>, v2: Vector3<T>) -> T
    where T: Copy + Add<T, Output = T> + Mul<T, Output = T>
{
    v1.x * v2.x + v1.y * v2.y + v1.z * v2.z
}

/// Product of the Euclidean magnitudes of a vector (and a normal) and
/// the cosine of the angle between them. A return value of zero means
/// both are orthogonal, a value if one means they are codirectional.
pub fn vec3_dot_nrm<T>(v1: Vector3<T>, n2: Normal3<T>) -> T
    where T: Copy + Add<T, Output = T> + Mul<T, Output = T>
{
    // DCHECK(!v1.HasNaNs() && !n2.HasNaNs());
    return v1.x * n2.x + v1.y * n2.y + v1.z * n2.z;
}


/// Computes the absolute value of the dot product.
pub fn vec3_abs_dot_vec3<T>(v1: Vector3<T>, v2: Vector3<T>) -> T
    where T: num::Float
{
    vec3_dot_vec3(v1, v2).abs()
}

/// Computes the absolute value of the dot product.
pub fn vec3_abs_dot_nrm<T>(v1: Vector3<T>, n2: Normal3<T>) -> T
    where T: num::Float
{
    vec3_dot_nrm(v1, n2).abs()
}

/// Given two vectors in 3D, the cross product is a vector that is
/// perpendicular to both of them.
pub fn vec3_cross_vec3(v1: Vector3f, v2: Vector3f) -> Vector3f
{
    let v1x: f64 = v1.x as f64;
    let v1y: f64 = v1.y as f64;
    let v1z: f64 = v1.z as f64;
    let v2x: f64 = v2.x as f64;
    let v2y: f64 = v2.y as f64;
    let v2z: f64 = v2.z as f64;
    Vector3f {
        x: ((v1y * v2z) - (v1z * v2y)) as Float,
        y: ((v1z * v2x) - (v1x * v2z)) as Float,
        z: ((v1x * v2y) - (v1y * v2x)) as Float,
    }
}

/// Given a vectors and a normal in 3D, the cross product is a vector
/// that is perpendicular to both of them.
pub fn vec3_cross_nrm(v1: Vector3f, v2: Normal3f) -> Vector3f
{
    let v1x: f64 = v1.x as f64;
    let v1y: f64 = v1.y as f64;
    let v1z: f64 = v1.z as f64;
    let v2x: f64 = v2.x as f64;
    let v2y: f64 = v2.y as f64;
    let v2z: f64 = v2.z as f64;
    Vector3f {
        x: ((v1y * v2z) - (v1z * v2y)) as Float,
        y: ((v1z * v2x) - (v1x * v2z)) as Float,
        z: ((v1x * v2y) - (v1y * v2x)) as Float,
    }
}

/// Compute a new vector pointing in the same direction but with unit
/// length.
pub fn vec3_normalize(v: Vector3f) -> Vector3f
{
    v / v.length()
}

/// Return the largest coordinate value.
pub fn vec3_max_component<T>(v: Vector3<T>) -> T
    where T: num::Float
{
    v.x.max(v.y.max(v.z))
}

/// Return the index of the component with the largest value.
pub fn vec3_max_dimension<T>(v: Vector3<T>) -> usize
    where T: std::cmp::PartialOrd
{
    if v.x > v.y {
        if v.x > v.z {
            0_usize
        } else {
            2_usize
        }
    } else {
        if v.y > v.z {
            1_usize
        } else {
            2_usize
        }
    }
}

/// Permute the coordinate values according to the povided
/// permutation.
pub fn vec3_permute<T>(v: Vector3<T>, x: usize, y: usize, z: usize) -> Vector3<T>
    where T: Copy
{
    let v3: Vec<T> = vec![v.x, v.y, v.z];
    let xp: T = v3[x];
    let yp: T = v3[y];
    let zp: T = v3[z];
    Vector3::<T> {
        x: xp,
        y: yp,
        z: zp,
    }
}

/// Construct a local coordinate system given only a single 3D vector.
pub fn vec3_coordinate_system(v1: &Vector3f, v2: &mut Vector3f, v3: &mut Vector3f)
{
    if v1.x.abs() > v1.y.abs() {
        *v2 = Vector3f {
            x: -v1.z,
            y: 0.0 as Float,
            z: v1.x,
        } / (v1.x * v1.x + v1.z * v1.z).sqrt();
    } else {
        *v2 = Vector3f {
            x: 0.0 as Float,
            y: v1.z,
            z: -v1.y,
        } / (v1.y * v1.y + v1.z * v1.z).sqrt();
    }
    *v3 = vec3_cross_vec3(*v1, *v2);
}

#[derive(Debug,Default,Copy,Clone)]
pub struct Point2<T> {
    pub x: T,
    pub y: T,
}

impl<T> Add<Point2<T>> for Point2<T>
    where T: Add<T, Output = T>
{
    type Output = Point2<T>;
    fn add(self, rhs: Point2<T>) -> Point2<T> {
        // TODO: DCHECK(!v.HasNaNs());
        Point2::<T> {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl<T> Add<Vector2<T>> for Point2<T>
    where T: Add<T, Output = T>
{
    type Output = Point2<T>;
    fn add(self, rhs: Vector2<T>) -> Point2<T> {
        // TODO: DCHECK(!v.HasNaNs());
        Point2::<T> {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl<T> Sub<Point2<T>> for Point2<T>
    where T: Sub<T, Output = T>
{
    type Output = Vector2<T>;
    fn sub(self, rhs: Point2<T>) -> Vector2<T> {
        Vector2::<T> {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

impl<T> Sub<Vector2<T>> for Point2<T>
    where T: Sub<T, Output = T>
{
    type Output = Point2<T>;
    fn sub(self, rhs: Vector2<T>) -> Point2<T> {
        Point2::<T> {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

impl<T> Mul<T> for Point2<T>
    where T: Copy + Mul<T, Output = T>
{
    type Output = Point2<T>;
    fn mul(self, rhs: T) -> Point2<T>
        where T: Copy + Mul<T, Output = T>
    {
        Point2::<T> {
            x: self.x * rhs,
            y: self.y * rhs,
        }
    }
}

impl<T> Index<u8> for Point2<T> {
    type Output = T;
    fn index(&self, index: u8) -> &T {
        match index {
            0 => &self.x,
            1 => &self.y,
            _ => panic!("Check failed: i >= 0 && i <= 1"),
        }
    }
}

/// Apply floor operation component-wise.
pub fn pnt2_floor<T>(p: Point2<T>) -> Point2<T>
    where T: num::Float
{
    Point2 {
        x: p.x.floor(),
        y: p.y.floor(),
    }
}

/// Apply ceil operation component-wise.
pub fn pnt2_ceil<T>(p: Point2<T>) -> Point2<T>
    where T: num::Float
{
    Point2 {
        x: p.x.ceil(),
        y: p.y.ceil(),
    }
}

/// Apply std::cmp::min operation component-wise.
pub fn pnt2_min_pnt2<T>(pa: Point2<T>, pb: Point2<T>) -> Point2<T>
    where T: Ord
{
    Point2 {
        x: std::cmp::min(pa.x, pb.x),
        y: std::cmp::min(pa.y, pb.y),
    }
}

/// Apply std::cmp::max operation component-wise.
pub fn pnt2_max_pnt2<T>(pa: Point2<T>, pb: Point2<T>) -> Point2<T>
    where T: Ord
{
    Point2 {
        x: std::cmp::max(pa.x, pb.x),
        y: std::cmp::max(pa.y, pb.y),
    }
}

/// Is a 2D point inside a 2D bound?
pub fn pnt2_inside_exclusive<T>(pt: Point2<T>, b: Bounds2<T>) -> bool
    where T: PartialOrd
{
    pt.x >= b.p_min.x && pt.x < b.p_max.x && pt.y >= b.p_min.y && pt.y < b.p_max.y
}

#[derive(Debug,Default,Copy,Clone)]
pub struct Point3<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T> Add<Point3<T>> for Point3<T>
    where T: Add<T, Output = T>
{
    type Output = Point3<T>;
    fn add(self, rhs: Point3<T>) -> Point3<T> {
        // TODO: DCHECK(!v.HasNaNs());
        Point3::<T> {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl<T> Add<Vector3<T>> for Point3<T>
    where T: Add<T, Output = T>
{
    type Output = Point3<T>;
    fn add(self, rhs: Vector3<T>) -> Point3<T> {
        // TODO: DCHECK(!v.HasNaNs());
        Point3::<T> {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl<T> AddAssign<Vector3<T>> for Point3<T>
    where T: AddAssign
{
    fn add_assign(&mut self, rhs: Vector3<T>) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl<T> Sub<Point3<T>> for Point3<T>
    where T: Sub<T, Output = T>
{
    type Output = Vector3<T>;
    fn sub(self, rhs: Point3<T>) -> Vector3<T> {
        Vector3::<T> {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl<T> Sub<Vector3<T>> for Point3<T>
    where T: Sub<T, Output = T>
{
    type Output = Point3<T>;
    fn sub(self, rhs: Vector3<T>) -> Point3<T> {
        Point3::<T> {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl<T> Mul<T> for Point3<T>
    where T: Copy + Mul<T, Output = T>
{
    type Output = Point3<T>;
    fn mul(self, rhs: T) -> Point3<T>
        where T: Copy + Mul<T, Output = T>
    {
        Point3::<T> {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl<T> MulAssign<T> for Point3<T>
    where T: Copy + MulAssign
{
    fn mul_assign(&mut self, rhs: T) {
        self.x *= rhs;
        self.y *= rhs;
        self.z *= rhs;
    }
}

// work around bug
// https://github.com/rust-lang/rust/issues/40395
impl Div<Float> for Point3<f32>
{
    type Output = Point3<f32>;
    fn div(self, rhs: Float) -> Point3<f32>
    {
        assert_ne!(rhs, 0.0 as Float);
        let inv: Float = 1.0 as Float / rhs;
        Point3::<f32> {
            x: self.x * inv,
            y: self.y * inv,
            z: self.z * inv,
        }
    }
}

// work around bug
// https://github.com/rust-lang/rust/issues/40395
impl DivAssign<Float> for Point3<f32>
{
    fn div_assign(&mut self, rhs: Float) {
        assert_ne!(rhs, 0.0 as Float);
        let inv: Float = 1.0 as Float / rhs;
        self.x *= inv;
        self.y *= inv;
        self.z *= inv;
    }
}

impl<T> Index<u8> for Point3<T> {
    type Output = T;
    fn index(&self, index: u8) -> &T {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Check failed: i >= 0 && i <= 2"),
        }
    }
}

impl<T> IndexMut<u8> for Point3<T> {
    fn index_mut(&mut self, index: u8) -> &mut T {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => panic!("Check failed: i >= 0 && i <= 2"),
        }
    }
}

/// Permute the coordinate values according to the povided
/// permutation.
pub fn pnt3_permute<T>(v: Point3<T>, x: usize, y: usize, z: usize) -> Point3<T>
    where T: Copy
{
    let v3: Vec<T> = vec![v.x, v.y, v.z];
    let xp: T = v3[x];
    let yp: T = v3[y];
    let zp: T = v3[z];
    Point3::<T> {
        x: xp,
        y: yp,
        z: zp,
    }
}

/// Apply floor operation component-wise.
pub fn pnt3_floor<T>(p: Point3<T>) -> Point3<T>
    where T: num::Float
{
    Point3 {
        x: p.x.floor(),
        y: p.y.floor(),
        z: p.z.floor(),
    }
}

/// Apply ceil operation component-wise.
pub fn pnt3_ceil<T>(p: Point3<T>) -> Point3<T>
    where T: num::Float
{
    Point3 {
        x: p.x.ceil(),
        y: p.y.ceil(),
        z: p.z.ceil(),
    }
}

/// Apply abs operation component-wise.
pub fn pnt3_abs<T>(p: Point3<T>) -> Point3<T>
    where T: num::Float
{
    Point3 {
        x: p.x.abs(),
        y: p.y.abs(),
        z: p.z.abs(),
    }
}

/// The distance between two points is the length of the vector
/// between them.
pub fn pnt3_distance<T>(p1: Point3<T>, p2: Point3<T>) -> T
    where T: num::Float + Sub<T, Output = T>
{
    (p1 - p2).length()
}

/// The distance squared between two points is the length of the
/// vector between them squared.
pub fn pnt3_distance_squared<T>(p1: Point3<T>, p2: Point3<T>) -> T
    where T: num::Float + Sub<T, Output = T>
{
    (p1 - p2).length_squared()
}

/// When tracing spawned rays leaving the intersection point p, we
/// offset their origins enough to ensure that they are past the
/// boundary of the error box and thus won't incorrectly re-intersect
/// the surface.
pub fn pnt3_offset_ray_origin(p: Point3f, p_error: Vector3f, n: Normal3f, w: Vector3f) -> Point3f {
    //     Float d = Dot(Abs(n), pError);
    let d: Float = nrm_dot_vec3(nrm_abs(n), p_error);
    // #ifdef PBRT_FLOAT_AS_DOUBLE
    //     // We have tons of precision; for now bump up the offset a bunch just
    //     // to be extra sure that we start on the right side of the surface
    //     // (In case of any bugs in the epsilons code...)
    //     d *= 1024.;
    // #endif
    let mut offset: Vector3f = Vector3f::from(n) * d;
    if vec3_dot_nrm(w, n) < 0.0 as Float {
        offset = -offset;
    }
    let mut po: Point3f = p + offset;
    // round offset point _po_ away from _p_
    for i in 0..3 {
        if offset[i] > 0.0 as Float {
            po[i] = next_float_up(po[i]);
        } else {
            if offset[i] < 0.0 as Float {
                po[i] = next_float_down(po[i]);
            }
        }
    }
    po
}

/// Calculate appropriate direction vector from two angles.
pub fn spherical_direction(sin_theta: Float, cos_theta: Float, phi: Float) -> Vector3f {
    Vector3f {
        x: sin_theta * phi.cos(),
        y: sin_theta * phi.sin(),
        z: cos_theta,
    }
}

/// Take three basis vectors representing the x, y, and z axes and
/// return the appropriate direction vector with respect to the
/// coordinate frame defined by them.
pub fn spherical_direction_vec3(sin_theta: Float, cos_theta: Float, phi: Float,
                                x: &Vector3f, y: &Vector3f, z: &Vector3f) -> Vector3f {
    *x * sin_theta * phi.cos() + *y * sin_theta * phi.sin() + *z * cos_theta
}

/// Conversion of a direction to spherical angles. Note that
/// **spherical_theta()** assumes that the vector **v** has been
/// normalized before being passed in.
pub fn spherical_theta(v: &Vector3f) -> Float {
    clamp_t(v.z, -1.0 as Float, 1.0 as Float).acos()
}

/// Conversion of a direction to spherical angles.
pub fn spherical_phi(v: &Vector3f) -> Float {
    let p: Float = v.y.atan2(v.x);
    if p < 0.0 as Float {
        p + 2.0 as Float * PI
    } else {
        p
    }
}

#[derive(Debug,Default,Copy,Clone)]
pub struct Normal3<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T> Add for Normal3<T>
    where T: Copy + Add<T, Output = T>
{
    type Output = Normal3<T>;
    fn add(self, rhs: Normal3<T>) -> Normal3<T> {
        Normal3::<T> {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl<T> Sub for Normal3<T>
    where T: Copy + Sub<T, Output = T>
{
    type Output = Normal3<T>;
    fn sub(self, rhs: Normal3<T>) -> Normal3<T> {
        Normal3::<T> {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl<T> Mul<T> for Normal3<T>
    where T: Copy + Mul<T, Output = T>
{
    type Output = Normal3<T>;
    fn mul(self, rhs: T) -> Normal3<T>
        where T: Copy + Mul<T, Output = T>
    {
        Normal3::<T> {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl<T> MulAssign<T> for Normal3<T>
    where T: Copy + MulAssign
{
    fn mul_assign(&mut self, rhs: T) {
        self.x *= rhs;
        self.y *= rhs;
        self.z *= rhs;
    }
}

impl<T> Neg for Normal3<T>
    where T: Copy + Neg<Output = T>
{
    type Output = Normal3<T>;
    fn neg(self) -> Normal3<T> {
        Normal3::<T> {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl<T> Normal3<T> {
    pub fn length_squared(&self) -> T
        where T: Copy + Add<T, Output = T> + Mul<T, Output = T>
    {
        self.x * self.x + self.y * self.y + self.z * self.z
    }
    pub fn length(&self) -> T
        where T: num::Float
    {
        self.length_squared().sqrt()
    }
}

impl<T> PartialEq for Normal3<T>
    where T: std::cmp::PartialOrd
{
    fn eq(&self, rhs: &Normal3<T>) -> bool {
        if self.x == rhs.x && self.y == rhs.y && self.z == rhs.z {
            true
        } else {
            false
        }
    }
    fn ne(&self, rhs: &Normal3<T>) -> bool {
        !self.eq(rhs)
    }
}

// work around bug
// https://github.com/rust-lang/rust/issues/40395
impl Div<Float> for Normal3<f32>
{
    type Output = Normal3<f32>;
    fn div(self, rhs: Float) -> Normal3<f32>
    {
        assert_ne!(rhs, 0.0 as Float);
        let inv: Float = 1.0 as Float / rhs;
        Normal3::<f32> {
            x: self.x * inv,
            y: self.y * inv,
            z: self.z * inv,
        }
    }
}

impl<T> From<Vector3<T>> for Normal3<T> {
    fn from(v: Vector3<T>) -> Self {
        // TODO: DCHECK(!v.HasNaNs());
        Normal3::<T> {
            x: v.x,
            y: v.y,
            z: v.z,
        }
    }
}

/// Given a normal and a vector in 3D, the cross product is a vector
/// that is perpendicular to both of them.
pub fn nrm_cross_vec3(n1: Normal3f, v2: Vector3f) -> Vector3f
{
    let n1x: f64 = n1.x as f64;
    let n1y: f64 = n1.y as f64;
    let n1z: f64 = n1.z as f64;
    let v2x: f64 = v2.x as f64;
    let v2y: f64 = v2.y as f64;
    let v2z: f64 = v2.z as f64;
    Vector3f {
        x: ((n1y * v2z) - (n1z * v2y)) as Float,
        y: ((n1z * v2x) - (n1x * v2z)) as Float,
        z: ((n1x * v2y) - (n1y * v2x)) as Float,
    }
}

/// Compute a new normal pointing in the same direction but with unit
/// length.
pub fn nrm_normalize(n: Normal3f) -> Normal3f
{
    n / n.length()
}

/// Product of the Euclidean magnitudes of a normal (and another
/// normal) and the cosine of the angle between them. A return value
/// of zero means both are orthogonal, a value if one means they are
/// codirectional.
pub fn nrm_dot_nrm<T>(n1: Normal3<T>, n2: Normal3<T>) -> T
    where T: Copy + Add<T, Output = T> + Mul<T, Output = T>
{
    // TODO: DCHECK(!n1.HasNaNs() && !n2.HasNaNs());
    n1.x * n2.x + n1.y * n2.y + n1.z * n2.z
}

/// Product of the Euclidean magnitudes of a normal (and a vector) and
/// the cosine of the angle between them. A return value of zero means
/// both are orthogonal, a value if one means they are codirectional.
pub fn nrm_dot_vec3<T>(n1: Normal3<T>, v2: Vector3<T>) -> T
    where T: Copy + Add<T, Output = T> + Mul<T, Output = T>
{
    // TODO: DCHECK(!n1.HasNaNs() && !v2.HasNaNs());
    n1.x * v2.x + n1.y * v2.y + n1.z * v2.z
}

/// Computes the absolute value of the dot product.
pub fn nrm_abs_dot_vec3<T>(n1: Normal3<T>, v2: Vector3<T>) -> T
    where T: num::Float
{
    nrm_dot_vec3(n1, v2).abs()
}

/// Return normal with the absolute value of each coordinate.
pub fn nrm_abs<T>(n: Normal3<T>) -> Normal3<T>
    where T: num::Float
{
    Normal3::<T> {
        x: n.x.abs(),
        y: n.y.abs(),
        z: n.z.abs(),
    }
}

/// Flip a surface normal so that it lies in the same hemisphere as a
/// given vector.
pub fn nrm_faceforward_vec3(n: Normal3f, v: Vector3f) -> Normal3f {
    if nrm_dot_vec3(n, v) < 0.0 as Float {
        -n
    } else {
        n
    }
}

/// Flip a surface normal so that it lies in the same hemisphere as a
/// given normal.
pub fn nrm_faceforward_nrm(n: Normal3f, n2: Normal3f) -> Normal3f {
    if nrm_dot_nrm(n, n2) < 0.0 as Float {
        -n
    } else {
        n
    }
}

pub type Bounds2f = Bounds2<Float>;
pub type Bounds2i = Bounds2<i32>;
pub type Bounds3f = Bounds3<Float>;
pub type Bounds3i = Bounds3<i32>;

#[derive(Debug,Default,Copy,Clone)]
pub struct Bounds2<T> {
    pub p_min: Point2<T>,
    pub p_max: Point2<T>,
}

impl<T> Bounds2<T> {
    pub fn new(p1: Point2<T>, p2: Point2<T>) -> Self
        where T: Copy + Ord
    {
        let p_min: Point2<T> = Point2::<T> {
            x: std::cmp::min(p1.x, p2.x),
            y: std::cmp::min(p1.y, p2.y),
        };
        let p_max: Point2<T> = Point2::<T> {
            x: std::cmp::max(p1.x, p2.x),
            y: std::cmp::max(p1.y, p2.y),
        };
        Bounds2::<T> {
            p_min: p_min,
            p_max: p_max,
        }
    }
    pub fn diagonal(&self) -> Vector2<T>
        where T: Copy + Sub<T, Output = T>
    {
        self.p_max - self.p_min
    }
    pub fn area(&self) -> T
        where T: Copy + Sub<T, Output = T> + Mul<T, Output = T>
    {
        let d: Vector2<T> = self.p_max - self.p_min;
        d.x * d.y
    }
}

pub struct Bounds2Iterator<'a> {
    p: Point2i,
    bounds: &'a Bounds2i,
}

impl<'a> Iterator for Bounds2Iterator<'a> {
    type Item = Point2i;

    fn next(&mut self) -> Option<Point2i> {
        self.p.x += 1;
        if self.p.x == self.bounds.p_max.x {
            self.p.x = self.bounds.p_min.x;
            self.p.y += 1;
        }
        if self.p.y == self.bounds.p_max.y {
            None
        } else {
            Some(self.p)
        }
    }
}

impl<'a> IntoIterator for &'a Bounds2i {
    type Item = Point2i;
    type IntoIter = Bounds2Iterator<'a>;

    fn into_iter(self) -> Self::IntoIter {
        Bounds2Iterator {
            // need to start 1 before p_min.x as next() will be called
            // to get the first element
            p: Point2i {
                x: self.p_min.x - 1,
                y: self.p_min.y,
            },
            bounds: self,
        }
    }
}

/// The intersection of two bounding boxes can be found by computing
/// the maximum of their two respective minimum coordinates and the
/// minimum of their maximum coordinates.
pub fn bnd2_intersect_bnd2<T>(b1: Bounds2<T>, b2: Bounds2<T>) -> Bounds2<T>
    where T: Copy + Ord
{
    Bounds2::<T> {
        p_min: Point2::<T> {
            x: std::cmp::max(b1.p_min.x, b2.p_min.x),
            y: std::cmp::max(b1.p_min.y, b2.p_min.y),
        },
        p_max: Point2::<T> {
            x: std::cmp::min(b1.p_max.x, b2.p_max.x),
            y: std::cmp::min(b1.p_max.y, b2.p_max.y),
        },
    }
}

#[derive(Debug,Copy,Clone)]
pub struct Bounds3<T> {
    pub p_min: Point3<T>,
    pub p_max: Point3<T>,
}

// work around bug
// https://github.com/rust-lang/rust/issues/40395
impl Default for Bounds3<f32> {
    fn default() -> Bounds3<f32> {
        let min_num: Float = std::f32::MIN;
        let max_num: Float = std::f32::MAX;
        // Bounds3f
        Bounds3::<f32> {
            p_min: Point3f {
                x: max_num,
                y: max_num,
                z: max_num,
            },
            p_max: Point3f {
                x: min_num,
                y: min_num,
                z: min_num,
            },
        }
    }
}

impl<T> Bounds3<T> {
    pub fn new(p1: Point3<T>, p2: Point3<T>) -> Self
        where T: num::Float
    {
        let p_min: Point3<T> = Point3::<T> {
            x: p1.x.min(p2.x),
            y: p1.y.min(p2.y),
            z: p1.z.min(p2.z),
        };
        let p_max: Point3<T> = Point3::<T> {
            x: p1.x.max(p2.x),
            y: p1.y.max(p2.y),
            z: p1.z.max(p2.z),
        };
        Bounds3::<T> {
            p_min: p_min,
            p_max: p_max,
        }
    }
    pub fn diagonal(&self) -> Vector3<T>
        where T: Copy + Sub<T, Output = T>
    {
        self.p_max - self.p_min
    }
    pub fn surface_area(&self) -> T
        where T: Copy + Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T>
    {
        let d: Vector3<T> = self.diagonal();
        // 2 * (d.x * d.y + d.x * d.z + d.y * d.z)
        let r: T = d.x * d.y + d.x * d.z + d.y * d.z;
        r + r // avoid '2 *'
    }
    pub fn maximum_extent(&self) -> u8
        where T: Copy + std::cmp::PartialOrd + Sub<T, Output = T>
    {
        let d: Vector3<T> = self.diagonal();
        if d.x > d.y && d.x > d.z {
            0_u8
        } else if d.y > d.z {
            1_u8
        } else {
            2_u8
        }
    }
    pub fn offset(&self, p: Point3<T>) -> Vector3<T>
        where T: Copy + std::cmp::PartialOrd + Sub<T, Output = T> + DivAssign<T>
    {
        let mut o: Vector3<T> = p - self.p_min;
        if self.p_max.x > self.p_min.x {
            o.x /= self.p_max.x - self.p_min.x;
        }
        if self.p_max.y > self.p_min.y {
            o.y /= self.p_max.y - self.p_min.y;
        }
        if self.p_max.z > self.p_min.z {
            o.z /= self.p_max.z - self.p_min.z;
        }
        o
    }
    pub fn bounding_sphere(b: &Bounds3f, center: &mut Point3f, radius: &mut Float) {
        let p_min: Point3f = b.p_min as Point3f;
        let p_max: Point3f = b.p_max as Point3f;
        let sum: Point3f = p_min + p_max;
        *center = sum / 2.0;
        let center_copy: Point3f = *center as Point3f;
        let is_inside: bool = pnt3_inside_bnd3(center_copy, *b);
        if is_inside {
            *radius = pnt3_distance(center_copy, p_max);
        } else {
            *radius = 0.0;
        }
    }
}

impl Bounds3<Float> {
    pub fn lerp(&self, t: Point3f) -> Point3f {
        Point3f {
            x: lerp(t.x, self.p_min.x as Float, self.p_max.x as Float),
            y: lerp(t.y, self.p_min.y as Float, self.p_max.y as Float),
            z: lerp(t.z, self.p_min.z as Float, self.p_max.z as Float),
        }
    }
    pub fn intersect_p(&self, ray: &Ray, inv_dir: &Vector3f, dir_is_neg: [u8; 3]) -> bool {
        // check for ray intersection against $x$ and $y$ slabs
        let mut t_min: Float = (self[dir_is_neg[0]].x - ray.o.x) * inv_dir.x;
        let mut t_max: Float = (self[1_u8 - dir_is_neg[0]].x - ray.o.x) * inv_dir.x;
        let ty_min: Float = (self[dir_is_neg[1]].y - ray.o.y) * inv_dir.y;
        let mut ty_max: Float = (self[1_u8 - dir_is_neg[1]].y - ray.o.y) * inv_dir.y;
        // update _t_max_ and _ty_max_ to ensure robust bounds intersection
        t_max *= 1.0 + 2.0 * gamma(3_i32);
        ty_max *= 1.0 + 2.0 * gamma(3_i32);
        if t_min > ty_max || ty_min > t_max {
            return false;
        }
        if ty_min > t_min {
            t_min = ty_min;
        }
        if ty_max < t_max {
            t_max = ty_max;
        }
        // check for ray intersection against $z$ slab
        let tz_min: Float = (self[dir_is_neg[2]].z - ray.o.z) * inv_dir.z;
        let mut tz_max: Float = (self[1_u8 - dir_is_neg[2]].z - ray.o.z) * inv_dir.z;
        // update _tz_max_ to ensure robust bounds intersection
        tz_max *= 1.0 + 2.0 * gamma(3_i32);
        if t_min > tz_max || tz_min > t_max {
            return false;
        }
        if tz_min > t_min {
            t_min = tz_min;
        }
        if tz_max < t_max {
            t_max = tz_max;
        }
        (t_min < ray.t_max) && (t_max > 0.0)
    }
}

impl<T> Index<u8> for Bounds3<T> {
    type Output = Point3<T>;

    fn index(&self, i: u8) -> &Point3<T> {
        match i {
            0 => &self.p_min,
            1 => &self.p_max,
            _ => panic!("Invalid index!"),
        }
    }
}

// /// Minimum squared distance from point to box; returns zero if point
// /// is inside.
// pub fn pnt3_distance_squared_bnd3(p: Point3f, b: Bounds3f) -> Float {
//     let dx: Float = (b.p_min.x - p.x).max(num::Zero::zero()).max(p.x - b.p_max.x);
//     let dy: Float = (b.p_min.y - p.y).max(num::Zero::zero()).max(p.y - b.p_max.y);
//     let dz: Float = (b.p_min.z - p.z).max(num::Zero::zero()).max(p.z - b.p_max.z);
//     dx * dx + dy * dy + dz * dz
// }

// /// Minimum distance from point to box; returns zero if point is
// /// inside.
// pub fn pnt3_distance_bnd3(p: Point3f, b: Bounds3f) -> Float {
//     pnt3_distance_squared_bnd3(p, b).sqrt()
// }

/// Given a bounding box and a point, the **bnd3_union_pnt3()**
/// function returns a new bounding box that encompasses that point as
/// well as the original box.
pub fn bnd3_union_pnt3<T>(b: Bounds3<T>, p: Point3<T>) -> Bounds3<T>
    where T: num::Float
{
    let p_min: Point3<T> = Point3::<T> {
        x: b.p_min.x.min(p.x),
        y: b.p_min.y.min(p.y),
        z: b.p_min.z.min(p.z),
    };
    let p_max: Point3<T> = Point3::<T> {
        x: b.p_max.x.max(p.x),
        y: b.p_max.y.max(p.y),
        z: b.p_max.z.max(p.z),
    };
    Bounds3::new(p_min, p_max)
}

/// Construct a new box that bounds the space encompassed by two other
/// bounding boxes.
pub fn bnd3_union_bnd3<T>(b1: Bounds3<T>, b2: Bounds3<T>) -> Bounds3<T>
    where T: num::Float
{
    let p_min: Point3<T> = Point3::<T> {
        x: b1.p_min.x.min(b2.p_min.x),
        y: b1.p_min.y.min(b2.p_min.y),
        z: b1.p_min.z.min(b2.p_min.z),
    };
    let p_max: Point3<T> = Point3::<T> {
        x: b1.p_max.x.max(b2.p_max.x),
        y: b1.p_max.y.max(b2.p_max.y),
        z: b1.p_max.z.max(b2.p_max.z),
    };
    Bounds3::new(p_min, p_max)
}

/// Determine if a given point is inside the bounding box.
pub fn pnt3_inside_bnd3(p: Point3f, b: Bounds3f) -> bool {
    p.x >= b.p_min.x && p.x <= b.p_max.x && p.y >= b.p_min.y && p.y <= b.p_max.y &&
        p.z >= b.p_min.z && p.z <= b.p_max.z
}

#[derive(Debug,Default,Copy,Clone)]
pub struct Ray {
    /// origin
    pub o: Point3f,
    /// direction
    pub d: Vector3f,
    /// limits the ray to a segment along its infinite extent
    pub t_max: Float,
    /// used for animations
    pub time: Float,
    /// in C++: 'class RayDifferential : public Ray'
    pub differential: Option<RayDifferential>,
}

impl Ray {
    // Point3f operator()(Float t) const { return o + d * t; }
    fn position(&self, t: Float) -> Point3f {
        self.o + self.d * t
    }
    // from class RayDifferential
    pub fn scale_differentials(&mut self, s: Float) {
        if let Some(d) = self.differential.iter_mut().next() {
            d.rx_origin = self.o + (d.rx_origin - self.o) * s;
            d.ry_origin = self.o + (d.ry_origin - self.o) * s;
            d.rx_direction = self.d + (d.rx_direction - self.d) * s;
            d.ry_direction = self.d + (d.ry_direction - self.d) * s;
        }
    }
}

#[derive(Debug,Default,Copy,Clone)]
pub struct RayDifferential {
    pub rx_origin: Point3f,
    pub ry_origin: Point3f,
    pub rx_direction: Vector3f,
    pub ry_direction: Vector3f,
}

// see transform.h

#[derive(Debug,Copy,Clone)]
pub struct Matrix4x4 {
    pub m: [[Float; 4]; 4],
}

impl Default for Matrix4x4 {
    fn default() -> Self {
        Matrix4x4 {
            m: [[1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0]],
        }
    }
}

impl Matrix4x4 {
    pub fn new(t00: Float,
               t01: Float,
               t02: Float,
               t03: Float,
               t10: Float,
               t11: Float,
               t12: Float,
               t13: Float,
               t20: Float,
               t21: Float,
               t22: Float,
               t23: Float,
               t30: Float,
               t31: Float,
               t32: Float,
               t33: Float)
               -> Self {
        Matrix4x4 {
            m: [[t00, t01, t02, t03],
                [t10, t11, t12, t13],
                [t20, t21, t22, t23],
                [t30, t31, t32, t33]],
        }
    }
    pub fn transpose(m: Matrix4x4) -> Matrix4x4 {
        Matrix4x4 {
            m: [[m.m[0][0], m.m[1][0], m.m[2][0], m.m[3][0]],
                [m.m[0][1], m.m[1][1], m.m[2][1], m.m[3][1]],
                [m.m[0][2], m.m[1][2], m.m[2][2], m.m[3][2]],
                [m.m[0][3], m.m[1][3], m.m[2][3], m.m[3][3]]],
        }
    }
    pub fn inverse(m: Matrix4x4) -> Matrix4x4 {
        let mut indxc = vec![0; 4];
        let mut indxr = vec![0; 4];
        let mut ipiv = vec![0; 4];
        let mut minv: Matrix4x4 = Matrix4x4::new(m.m[0][0],
                                                 m.m[0][1],
                                                 m.m[0][2],
                                                 m.m[0][3],
                                                 m.m[1][0],
                                                 m.m[1][1],
                                                 m.m[1][2],
                                                 m.m[1][3],
                                                 m.m[2][0],
                                                 m.m[2][1],
                                                 m.m[2][2],
                                                 m.m[2][3],
                                                 m.m[3][0],
                                                 m.m[3][1],
                                                 m.m[3][2],
                                                 m.m[3][3]);
        for i in 0..4 {
            let mut irow = 0;
            let mut icol = 0;
            let mut big: Float = 0.0;
            // choose pivot
            for j in 0..4 {
                if ipiv[j] != 1 {
                    for k in 0..4 {
                        if ipiv[k] == 0 {
                            let abs: Float = (minv.m[j][k]).abs();
                            if abs >= big {
                                big = abs;
                                irow = j;
                                icol = k;
                            }
                        } else {
                            if ipiv[k] > 1 {
                                println!("Singular matrix in MatrixInvert");
                            }
                        }
                    }
                }
            }
            ipiv[icol] = ipiv[icol] + 1;
            // swap rows _irow_ and _icol_ for pivot
            if irow != icol {
                for k in 0..4 {
                    // C++: std::swap(minv[irow][k], minv[icol][k]);
                    let swap = minv.m[irow][k];
                    minv.m[irow][k] = minv.m[icol][k];
                    minv.m[icol][k] = swap;
                }
            }
            indxr[i] = irow;
            indxc[i] = icol;
            if minv.m[icol][icol] == 0.0 {
                println!("Singular matrix in MatrixInvert");
            }
            // set $m[icol][icol]$ to one by scaling row _icol_ appropriately
            let pivinv: Float = 1.0 / minv.m[icol][icol];
            minv.m[icol][icol] = 1.0;
            for j in 0..4 {
                minv.m[icol][j] = minv.m[icol][j] * pivinv;
            }
            // subtract this row from others to zero out their columns
            for j in 0..4 {
                if j != icol {
                    let save: Float = minv.m[j][icol];
                    minv.m[j][icol] = 0.0;
                    for k in 0..4 {
                        minv.m[j][k] = minv.m[j][k] - (minv.m[icol][k] * save);
                    }
                }
            }
        }
        // swap columns to reflect permutation
        for i in 0..4 {
            let j = 3 - i;
            if indxr[j] != indxc[j] {
                for k in 0..4 {
                    // C++: std::swap(minv[k][indxr[j]], minv[k][indxc[j]]);
                    let swap = minv.m[k][indxr[j]];
                    minv.m[k][indxr[j]] = minv.m[k][indxc[j]];
                    minv.m[k][indxc[j]] = swap;
                }
            }
        }
        minv
    }
}

// see transform.cpp

/// Finds the closed-form solution of a 2x2 linear system.
pub fn solve_linear_system_2x2(a: [[Float; 2]; 2],
                               b: [Float; 2],
                               x0: &mut Float,
                               x1: &mut Float)
                               -> bool {
    let det: Float = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    if det.abs() < 1e-10 as Float {
        return false;
    }
    *x0 = (a[1][1] * b[0] - a[0][1] * b[1]) / det;
    *x1 = (a[0][0] * b[1] - a[1][0] * b[0]) / det;
    if (*x0).is_nan() || (*x1).is_nan() {
        return false;
    }
    true
}

/// The product of two matrices.
pub fn mtx_mul(m1: Matrix4x4, m2: Matrix4x4) -> Matrix4x4 {
    let mut r: Matrix4x4 = Matrix4x4::default();
    for i in 0..4 {
        for j in 0..4 {
            r.m[i][j] = m1.m[i][0] * m2.m[0][j] + m1.m[i][1] * m2.m[1][j] +
                m1.m[i][2] * m2.m[2][j] + m1.m[i][3] * m2.m[3][j];
        }
    }
    r
}

#[derive(Debug,Copy,Clone)]
pub struct Transform {
    pub m: Matrix4x4,
    pub m_inv: Matrix4x4,
}

impl Default for Transform {
    fn default() -> Self {
        Transform {
            m: Matrix4x4::default(),
            m_inv: Matrix4x4::default(),
        }
    }
}

impl Mul for Transform {
    type Output = Transform;
    fn mul(self, rhs: Transform) -> Transform {
        Transform {
            m: mtx_mul(self.m, rhs.m),
            m_inv: mtx_mul(rhs.m_inv, self.m_inv),
        }
    }
}

impl Transform {
    pub fn new(t00: Float,
               t01: Float,
               t02: Float,
               t03: Float,
               t10: Float,
               t11: Float,
               t12: Float,
               t13: Float,
               t20: Float,
               t21: Float,
               t22: Float,
               t23: Float,
               t30: Float,
               t31: Float,
               t32: Float,
               t33: Float)
               -> Self {
        Transform {
            m: Matrix4x4::new(t00,
                              t01,
                              t02,
                              t03,
                              t10,
                              t11,
                              t12,
                              t13,
                              t20,
                              t21,
                              t22,
                              t23,
                              t30,
                              t31,
                              t32,
                              t33),
            m_inv: Matrix4x4::inverse(Matrix4x4::new(t00,
                                                     t01,
                                                     t02,
                                                     t03,
                                                     t10,
                                                     t11,
                                                     t12,
                                                     t13,
                                                     t20,
                                                     t21,
                                                     t22,
                                                     t23,
                                                     t30,
                                                     t31,
                                                     t32,
                                                     t33)),
        }
    }
    pub fn inverse(t: Transform) -> Transform {
        Transform {
            m: t.m_inv,
            m_inv: t.m,
        }
    }
    pub fn translate(delta: Vector3f) -> Transform {
        Transform {
            m: Matrix4x4::new(1.0,
                              0.0,
                              0.0,
                              delta.x,
                              0.0,
                              1.0,
                              0.0,
                              delta.y,
                              0.0,
                              0.0,
                              1.0,
                              delta.z,
                              0.0,
                              0.0,
                              0.0,
                              1.0),
            m_inv: Matrix4x4::new(1.0,
                                  0.0,
                                  0.0,
                                  -delta.x,
                                  0.0,
                                  1.0,
                                  0.0,
                                  -delta.y,
                                  0.0,
                                  0.0,
                                  1.0,
                                  -delta.z,
                                  0.0,
                                  0.0,
                                  0.0,
                                  1.0),
        }
    }
    pub fn scale(x: Float, y: Float, z: Float) -> Transform {
        Transform {
            m: Matrix4x4::new(x,
                              0.0,
                              0.0,
                              0.0,
                              0.0,
                              y,
                              0.0,
                              0.0,
                              0.0,
                              0.0,
                              z,
                              0.0,
                              0.0,
                              0.0,
                              0.0,
                              1.0),
            m_inv: Matrix4x4::new(1.0 / x,
                                  0.0,
                                  0.0,
                                  0.0,
                                  0.0,
                                  1.0 / y,
                                  0.0,
                                  0.0,
                                  0.0,
                                  0.0,
                                  1.0 / z,
                                  0.0,
                                  0.0,
                                  0.0,
                                  0.0,
                                  1.0),
        }
    }
    pub fn rotate_x(theta: Float) -> Transform {
        let sin_theta: Float = radians(theta).sin();
        let cos_theta: Float = radians(theta).cos();
        let m = Matrix4x4::new(1.0,
                               0.0,
                               0.0,
                               0.0,
                               0.0,
                               cos_theta,
                               -sin_theta,
                               0.0,
                               0.0,
                               sin_theta,
                               cos_theta,
                               0.0,
                               0.0,
                               0.0,
                               0.0,
                               1.0);
        Transform {
            m: m,
            m_inv: Matrix4x4::transpose(m),
        }
    }
    pub fn rotate_y(theta: Float) -> Transform {
        let sin_theta: Float = radians(theta).sin();
        let cos_theta: Float = radians(theta).cos();
        let m = Matrix4x4::new(cos_theta,
                               0.0,
                               sin_theta,
                               0.0,
                               0.0,
                               1.0,
                               0.0,
                               0.0,
                               -sin_theta,
                               0.0,
                               cos_theta,
                               0.0,
                               0.0,
                               0.0,
                               0.0,
                               1.0);
        Transform {
            m: m,
            m_inv: Matrix4x4::transpose(m),
        }
    }
    pub fn rotate_z(theta: Float) -> Transform {
        let sin_theta: Float = radians(theta).sin();
        let cos_theta: Float = radians(theta).cos();
        let m = Matrix4x4::new(cos_theta,
                               -sin_theta,
                               0.0,
                               0.0,
                               sin_theta,
                               cos_theta,
                               0.0,
                               0.0,
                               0.0,
                               0.0,
                               1.0,
                               0.0,
                               0.0,
                               0.0,
                               0.0,
                               1.0);
        Transform {
            m: m,
            m_inv: Matrix4x4::transpose(m),
        }
    }
    pub fn rotate(theta: Float, axis: Vector3f) -> Transform {
        let a: Vector3f = vec3_normalize(axis);
        let sin_theta: Float = radians(theta).sin();
        let cos_theta: Float = radians(theta).cos();
        let mut m = Matrix4x4::default();
        // compute rotation of first basis vector
        m.m[0][0] = a.x * a.x + (1.0 - a.x * a.x) * cos_theta;
        m.m[0][1] = a.x * a.y * (1.0 - cos_theta) - a.z * sin_theta;
        m.m[0][2] = a.x * a.z * (1.0 - cos_theta) + a.y * sin_theta;
        m.m[0][3] = 0.0;
        // compute rotations of second basis vectors
        m.m[1][0] = a.x * a.y * (1.0 - cos_theta) + a.z * sin_theta;
        m.m[1][1] = a.y * a.y + (1.0 - a.y * a.y) * cos_theta;
        m.m[1][2] = a.y * a.z * (1.0 - cos_theta) - a.x * sin_theta;
        m.m[1][3] = 0.0;
        // compute rotations of third basis vectors
        m.m[2][0] = a.x * a.z * (1.0 - cos_theta) - a.y * sin_theta;
        m.m[2][1] = a.y * a.z * (1.0 - cos_theta) + a.x * sin_theta;
        m.m[2][2] = a.z * a.z + (1.0 - a.z * a.z) * cos_theta;
        m.m[2][3] = 0.0;
        Transform {
            m: m,
            m_inv: Matrix4x4::transpose(m),
        }
    }
    pub fn look_at(pos: Point3f, look: Point3f, up: Vector3f) -> Transform {
        let mut camera_to_world = Matrix4x4::default();
        // initialize fourth column of viewing matrix
        camera_to_world.m[0][3] = pos.x;
        camera_to_world.m[1][3] = pos.y;
        camera_to_world.m[2][3] = pos.z;
        camera_to_world.m[3][3] = 1.0;
        // initialize first three columns of viewing matrix
        let dir: Vector3f = vec3_normalize(look - pos);
        if vec3_cross_vec3(vec3_normalize(up), dir).length() == 0.0 {
            println!("\"up\" vector ({}, {}, {}) and viewing direction ({}, {}, {}) passed to \
                      LookAt are pointing in the same direction.  Using the identity \
                      transformation.",
                     up.x,
                     up.y,
                     up.z,
                     dir.x,
                     dir.y,
                     dir.z);
            Transform::default()
        } else {
            let left: Vector3f = vec3_normalize(vec3_cross_vec3(vec3_normalize(up), dir));
            let new_up: Vector3f = vec3_cross_vec3(dir, left);
            camera_to_world.m[0][0] = left.x;
            camera_to_world.m[1][0] = left.y;
            camera_to_world.m[2][0] = left.z;
            camera_to_world.m[3][0] = 0.0;
            camera_to_world.m[0][1] = new_up.x;
            camera_to_world.m[1][1] = new_up.y;
            camera_to_world.m[2][1] = new_up.z;
            camera_to_world.m[3][1] = 0.0;
            camera_to_world.m[0][2] = dir.x;
            camera_to_world.m[1][2] = dir.y;
            camera_to_world.m[2][2] = dir.z;
            camera_to_world.m[3][2] = 0.0;
            Transform {
                m: Matrix4x4::inverse(camera_to_world),
                m_inv: camera_to_world,
            }
        }
    }
    pub fn perspective(fov: Float, n: Float, f: Float) -> Transform {
        // perform projective divide for perspective projection
        let persp = Matrix4x4::new(1.0,
                                   0.0,
                                   0.0,
                                   0.0,
                                   0.0,
                                   1.0,
                                   0.0,
                                   0.0,
                                   0.0,
                                   0.0,
                                   f / (f - n),
                                   -f * n / (f - n),
                                   0.0,
                                   0.0,
                                   1.0,
                                   0.0);
        // scale canonical perspective view to specified field of view
        let inv_tan_ang: Float = 1.0 / (radians(fov) / 2.0).tan();
        let scale: Transform = Transform::scale(inv_tan_ang, inv_tan_ang, 1.0);
        let persp_trans: Transform = Transform {
            m: persp,
            m_inv: Matrix4x4::inverse(persp),
        };
        scale * persp_trans
    }
    pub fn transform_point(&self, p: Point3<Float>) -> Point3<Float> {
        let x: Float = p.x;
        let y: Float = p.y;
        let z: Float = p.z;
        let xp: Float = self.m.m[0][0] * x + self.m.m[0][1] * y + self.m.m[0][2] * z +
            self.m.m[0][3];
        let yp: Float = self.m.m[1][0] * x + self.m.m[1][1] * y + self.m.m[1][2] * z +
            self.m.m[1][3];
        let zp: Float = self.m.m[2][0] * x + self.m.m[2][1] * y + self.m.m[2][2] * z +
            self.m.m[2][3];
        let wp: Float = self.m.m[3][0] * x + self.m.m[3][1] * y + self.m.m[3][2] * z +
            self.m.m[3][3];
        assert!(wp != 0.0, "wp = {:?} != 0.0", wp);
        if wp == 1. {
            Point3::<Float> {
                x: xp,
                y: yp,
                z: zp,
            }
        } else {
            let inv: Float = 1.0 as Float / wp;
            Point3::<Float> {
                x: inv * xp,
                y: inv * yp,
                z: inv * zp,
            }
        }
    }
    pub fn transform_vector(&self, v: Vector3<Float>) -> Vector3<Float> {
        let x: Float = v.x;
        let y: Float = v.y;
        let z: Float = v.z;
        Vector3::<Float> {
            x: self.m.m[0][0] * x + self.m.m[0][1] * y + self.m.m[0][2] * z,
            y: self.m.m[1][0] * x + self.m.m[1][1] * y + self.m.m[1][2] * z,
            z: self.m.m[2][0] * x + self.m.m[2][1] * y + self.m.m[2][2] * z,
        }
    }
    pub fn transform_normal(&self, n: Normal3<Float>) -> Normal3<Float> {
        let x: Float = n.x;
        let y: Float = n.y;
        let z: Float = n.z;
        Normal3::<Float> {
            x: self.m_inv.m[0][0] * x + self.m_inv.m[1][0] * y + self.m_inv.m[2][0] * z,
            y: self.m_inv.m[0][1] * x + self.m_inv.m[1][1] * y + self.m_inv.m[2][1] * z,
            z: self.m_inv.m[0][2] * x + self.m_inv.m[1][2] * y + self.m_inv.m[2][2] * z,
        }
    }
    pub fn transform_ray(&self, r: &mut Ray) {
        // Ray tr = (*this)(Ray(r));
        let mut o_error: Vector3f = Vector3f::default();
        let mut o: Point3f = self.transform_point_with_error(r.o, &mut o_error);
        let d: Vector3f = self.transform_vector(r.d);
        let length_squared: Float = d.length_squared();
        let mut t_max: Float = r.t_max;
        if length_squared > 0.0 as Float {
            let dt: Float = vec3_dot_vec3(d.abs(), o_error) / length_squared;
            o += d * dt;
            t_max -= dt;
        }
        r.o = o;
        r.d = d;
        r.t_max = t_max;
        if let Some(d) = r.differential {
            let diff: RayDifferential = RayDifferential {
                rx_origin: self.transform_point(d.rx_origin),
                ry_origin: self.transform_point(d.ry_origin),
                rx_direction: self.transform_vector(d.rx_direction),
                ry_direction: self.transform_vector(d.ry_direction),
            };
            r.differential = Some(diff);
        }
    }
    pub fn transform_bounds(&self, b: Bounds3f) -> Bounds3f {
        let m: Transform = *self;
        let p: Point3f = self.transform_point(Point3f {
            x: b.p_min.x,
            y: b.p_min.y,
            z: b.p_min.z,
        });
        let mut ret: Bounds3f = Bounds3f {
            p_min: p,
            p_max: p,
        };
        ret = bnd3_union_pnt3(ret,
                              m.transform_point(Point3f {
                                  x: b.p_max.x,
                                  y: b.p_min.y,
                                  z: b.p_min.z,
                              }));
        ret = bnd3_union_pnt3(ret,
                              m.transform_point(Point3f {
                                  x: b.p_min.x,
                                  y: b.p_max.y,
                                  z: b.p_min.z,
                              }));
        ret = bnd3_union_pnt3(ret,
                              m.transform_point(Point3f {
                                  x: b.p_min.x,
                                  y: b.p_min.y,
                                  z: b.p_max.z,
                              }));
        ret = bnd3_union_pnt3(ret,
                              m.transform_point(Point3f {
                                  x: b.p_min.x,
                                  y: b.p_max.y,
                                  z: b.p_max.z,
                              }));
        ret = bnd3_union_pnt3(ret,
                              m.transform_point(Point3f {
                                  x: b.p_max.x,
                                  y: b.p_max.y,
                                  z: b.p_min.z,
                              }));
        ret = bnd3_union_pnt3(ret,
                              m.transform_point(Point3f {
                                  x: b.p_max.x,
                                  y: b.p_min.y,
                                  z: b.p_max.z,
                              }));
        ret = bnd3_union_pnt3(ret,
                              m.transform_point(Point3f {
                                  x: b.p_max.x,
                                  y: b.p_max.y,
                                  z: b.p_max.z,
                              }));
        ret
    }
    pub fn transform_point_with_error(&self,
                                      p: Point3<Float>,
                                      p_error: &mut Vector3<Float>)
                                      -> Point3<Float> {
        let x: Float = p.x;
        let y: Float = p.y;
        let z: Float = p.z;
        // compute transformed coordinates from point _pt_
        let xp: Float = self.m.m[0][0] * x + self.m.m[0][1] * y + self.m.m[0][2] * z +
            self.m.m[0][3];
        let yp: Float = self.m.m[1][0] * x + self.m.m[1][1] * y + self.m.m[1][2] * z +
            self.m.m[1][3];
        let zp: Float = self.m.m[2][0] * x + self.m.m[2][1] * y + self.m.m[2][2] * z +
            self.m.m[2][3];
        let wp: Float = self.m.m[3][0] * x + self.m.m[3][1] * y + self.m.m[3][2] * z +
            self.m.m[3][3];
        // compute absolute error for transformed point
        let x_abs_sum: Float = (self.m.m[0][0] * x).abs() + (self.m.m[0][1] * y).abs() +
            (self.m.m[0][2] * z).abs() +
            self.m.m[0][3].abs();
        let y_abs_sum: Float = (self.m.m[1][0] * x).abs() + (self.m.m[1][1] * y).abs() +
            (self.m.m[1][2] * z).abs() +
            self.m.m[1][3].abs();
        let z_abs_sum: Float = (self.m.m[2][0] * x).abs() + (self.m.m[2][1] * y).abs() +
            (self.m.m[2][2] * z).abs() +
            self.m.m[2][3].abs();
        *p_error = Vector3::<Float> {
            x: x_abs_sum,
            y: y_abs_sum,
            z: z_abs_sum,
        } * gamma(3i32);
        assert!(wp != 0.0, "wp = {:?} != 0.0", wp);
        if wp == 1. {
            Point3::<Float> {
                x: xp,
                y: yp,
                z: zp,
            }
        } else {
            let inv: Float = 1.0 as Float / wp;
            Point3::<Float> {
                x: inv * xp,
                y: inv * yp,
                z: inv * zp,
            }
        }
    }
    pub fn transform_point_with_abs_error(&self,
                                          pt: Point3<Float>,
                                          pt_error: &Vector3<Float>,
                                          abs_error: &mut Vector3<Float>)
                                          -> Point3<Float> {
        let x: Float = pt.x;
        let y: Float = pt.y;
        let z: Float = pt.z;
        // compute transformed coordinates from point _pt_
        let xp: Float = self.m.m[0][0] * x + self.m.m[0][1] * y + self.m.m[0][2] * z +
            self.m.m[0][3];
        let yp: Float = self.m.m[1][0] * x + self.m.m[1][1] * y + self.m.m[1][2] * z +
            self.m.m[1][3];
        let zp: Float = self.m.m[2][0] * x + self.m.m[2][1] * y + self.m.m[2][2] * z +
            self.m.m[2][3];
        let wp: Float = self.m.m[3][0] * x + self.m.m[3][1] * y + self.m.m[3][2] * z +
            self.m.m[3][3];
        abs_error.x = (gamma(3i32) + 1.0 as Float) *
            (self.m.m[0][0].abs() * pt_error.x + self.m.m[0][1].abs() * pt_error.y +
             self.m.m[0][2].abs() * pt_error.z) +
            gamma(3i32) * ((self.m.m[0][0] * x).abs() +
                           (self.m.m[0][1] * y).abs() +
                           (self.m.m[0][2] * z).abs() + self.m.m[0][3].abs());
        abs_error.y = (gamma(3i32) + 1.0 as Float) *
            (self.m.m[1][0].abs() * pt_error.x + self.m.m[1][1].abs() * pt_error.y +
             self.m.m[1][2].abs() * pt_error.z) +
            gamma(3i32) * ((self.m.m[1][0] * x).abs() +
                           (self.m.m[1][1] * y).abs() +
                           (self.m.m[1][2] * z).abs() + self.m.m[1][3].abs());
        abs_error.z = (gamma(3i32) + 1.0 as Float) *
            (self.m.m[2][0].abs() * pt_error.x + self.m.m[2][1].abs() * pt_error.y +
             self.m.m[2][2].abs() * pt_error.z) +
            gamma(3i32) * ((self.m.m[2][0] * x).abs() +
                           (self.m.m[2][1] * y).abs() +
                           (self.m.m[2][2] * z).abs() + self.m.m[2][3].abs());
        assert!(wp != 0.0, "wp = {:?} != 0.0", wp);
        if wp == 1. {
            Point3::<Float> {
                x: xp,
                y: yp,
                z: zp,
            }
        } else {
            let inv: Float = 1.0 as Float / wp;
            Point3::<Float> {
                x: inv * xp,
                y: inv * yp,
                z: inv * zp,
            }
        }
    }
    pub fn transform_vector_with_error(&self,
                                       v: Vector3<Float>,
                                       abs_error: &mut Vector3<Float>)
                                       -> Vector3<Float> {
        let x: Float = v.x;
        let y: Float = v.y;
        let z: Float = v.z;
        let gamma: Float = gamma(3i32);
        abs_error.x = gamma *
            ((self.m.m[0][0] * v.x).abs() + (self.m.m[0][1] * v.y).abs() +
             (self.m.m[0][2] * v.z).abs());
        abs_error.y = gamma *
            ((self.m.m[1][0] * v.x).abs() + (self.m.m[1][1] * v.y).abs() +
             (self.m.m[1][2] * v.z).abs());
        abs_error.z = gamma *
            ((self.m.m[2][0] * v.x).abs() + (self.m.m[2][1] * v.y).abs() +
             (self.m.m[2][2] * v.z).abs());
        Vector3::<Float> {
            x: self.m.m[0][0] * x + self.m.m[0][1] * y + self.m.m[0][2] * z,
            y: self.m.m[1][0] * x + self.m.m[1][1] * y + self.m.m[1][2] * z,
            z: self.m.m[2][0] * x + self.m.m[2][1] * y + self.m.m[2][2] * z,
        }
    }
    pub fn transform_ray_with_error(&self,
                                    r: &Ray,
                                    o_error: &mut Vector3<Float>,
                                    d_error: &mut Vector3<Float>)
                                    -> Ray {
        let mut o: Point3f = self.transform_point_with_error(r.o, o_error);
        let d: Vector3f = self.transform_vector_with_error(r.d, d_error);
        let length_squared: Float = d.length_squared();
        if length_squared > 0.0 {
            let dt: Float = vec3_dot_vec3(d.abs(), *o_error) / length_squared;
            o += d * dt;
        }
        Ray {
            o: o,
            d: d,
            t_max: r.t_max,
            time: r.time,
            differential: None,
        }
    }
    pub fn transform_surface_interaction(&self, si: &SurfaceInteraction) -> SurfaceInteraction {
        let mut ret: SurfaceInteraction = SurfaceInteraction::default();
        // transform _p_ and _pError_ in _SurfaceInteraction_
        ret.p = self.transform_point_with_abs_error(si.p, &si.p_error, &mut ret.p_error);
        // transform remaining members of _SurfaceInteraction_
        ret.n = nrm_normalize(self.transform_normal(si.n));
        ret.wo = vec3_normalize(self.transform_vector(si.wo));
        ret.time = si.time;
        ret.uv = si.uv;
        ret.shape = None; // TODO? si.shape;
        ret.dpdu = self.transform_vector(si.dpdu);
        ret.dpdv = self.transform_vector(si.dpdv);
        ret.dndu = self.transform_normal(si.dndu);
        ret.dndv = self.transform_normal(si.dndv);
        ret.shading.n = nrm_normalize(self.transform_normal(si.shading.n));
        ret.shading.dpdu = self.transform_vector(si.shading.dpdu);
        ret.shading.dpdv = self.transform_vector(si.shading.dpdv);
        ret.shading.dndu = self.transform_normal(si.shading.dndu);
        ret.shading.dndv = self.transform_normal(si.shading.dndv);
        ret.dudx = si.dudx;
        ret.dvdx = si.dvdx;
        ret.dudy = si.dudy;
        ret.dvdy = si.dvdy;
        ret.dpdx = self.transform_vector(si.dpdx);
        ret.dpdy = self.transform_vector(si.dpdy);
        ret.bsdf = si.bsdf.clone();
        ret.primitive = None; // TODO? si.primitive;
        ret.shading.n = nrm_faceforward_nrm(ret.shading.n, ret.n);
        ret
    }
}

#[derive(Debug,Default,Copy,Clone)]
pub struct DerivativeTerm {
    kc: Float,
    kx: Float,
    ky: Float,
    kz: Float,
}

#[derive(Debug,Default,Copy,Clone)]
pub struct AnimatedTransform {
    start_transform: Transform,
    end_transform: Transform,
    start_time: Float,
    end_time: Float,
    actually_animated: bool,
    t: [Vector3f; 2],
    r: [Quaternion; 2],
    s: [Matrix4x4; 2],
    has_rotation: bool,
    c1: [DerivativeTerm; 3],
    c2: [DerivativeTerm; 3],
    c3: [DerivativeTerm; 3],
    c4: [DerivativeTerm; 3],
    c5: [DerivativeTerm; 3],
}

impl AnimatedTransform {
    pub fn new(start_transform: &Transform,
               start_time: Float,
               end_transform: &Transform,
               end_time: Float)
               -> Self {
        let mut at: AnimatedTransform = AnimatedTransform::default();
        at.start_transform = start_transform.clone();
        at.end_transform = end_transform.clone();
        at.start_time = start_time;
        at.end_time = end_time;
        AnimatedTransform::decompose(&start_transform.m, &mut at.t[0], &mut at.r[0], &mut at.s[0]);
        AnimatedTransform::decompose(&end_transform.m, &mut at.t[1], &mut at.r[1], &mut at.s[1]);
        // flip _r[1]_ if needed to select shortest path
        if quat_dot(at.r[0], at.r[1]) < 0.0 {
            at.r[1] = -at.r[1];
        }
        at.has_rotation = quat_dot(at.r[0], at.r[1]) < 0.9995;
        // compute terms of motion derivative function
        if at.has_rotation {
            let cos_theta: Float = quat_dot(at.r[0], at.r[1]);
            let theta: Float = (clamp_t(cos_theta, -1.0, 1.0)).acos();
            let qperp: Quaternion = quat_normalize(at.r[1] - at.r[0] * cos_theta);
            let t0x: Float = at.t[0].x;
            let t0y: Float = at.t[0].y;
            let t0z: Float = at.t[0].z;
            let t1x: Float = at.t[1].x;
            let t1y: Float = at.t[1].y;
            let t1z: Float = at.t[1].z;
            let q1x: Float = at.r[0].v.x;
            let q1y: Float = at.r[0].v.y;
            let q1z: Float = at.r[0].v.z;
            let q1w: Float = at.r[0].w;
            let qperpx: Float = qperp.v.x;
            let qperpy: Float = qperp.v.y;
            let qperpz: Float = qperp.v.z;
            let qperpw: Float = qperp.w;
            let s000: Float = at.s[0].m[0][0];
            let s001: Float = at.s[0].m[0][1];
            let s002: Float = at.s[0].m[0][2];
            let s010: Float = at.s[0].m[1][0];
            let s011: Float = at.s[0].m[1][1];
            let s012: Float = at.s[0].m[1][2];
            let s020: Float = at.s[0].m[2][0];
            let s021: Float = at.s[0].m[2][1];
            let s022: Float = at.s[0].m[2][2];
            let s100: Float = at.s[1].m[0][0];
            let s101: Float = at.s[1].m[0][1];
            let s102: Float = at.s[1].m[0][2];
            let s110: Float = at.s[1].m[1][0];
            let s111: Float = at.s[1].m[1][1];
            let s112: Float = at.s[1].m[1][2];
            let s120: Float = at.s[1].m[2][0];
            let s121: Float = at.s[1].m[2][1];
            let s122: Float = at.s[1].m[2][2];
            at.c1[0] = DerivativeTerm {
                kc: -t0x + t1x,
                kx: (-1.0 + q1y * q1y + q1z * q1z + qperpy * qperpy + qperpz * qperpz) * s000 +
                    q1w * q1z * s010 - qperpx * qperpy * s010 +
                    qperpw * qperpz * s010 - q1w * q1y * s020 -
                    qperpw * qperpy * s020 -
                    qperpx * qperpz * s020 + s100 -
                    q1y * q1y * s100 - q1z * q1z * s100 -
                    qperpy * qperpy * s100 -
                    qperpz * qperpz * s100 - q1w * q1z * s110 +
                    qperpx * qperpy * s110 -
                    qperpw * qperpz * s110 + q1w * q1y * s120 +
                    qperpw * qperpy * s120 +
                    qperpx * qperpz * s120 +
                    q1x * (-(q1y * s010) - q1z * s020 + q1y * s110 + q1z * s120),
                ky: (-1.0 + q1y * q1y + q1z * q1z + qperpy * qperpy + qperpz * qperpz) * s001 +
                    q1w * q1z * s011 - qperpx * qperpy * s011 +
                    qperpw * qperpz * s011 - q1w * q1y * s021 -
                    qperpw * qperpy * s021 -
                    qperpx * qperpz * s021 + s101 -
                    q1y * q1y * s101 - q1z * q1z * s101 -
                    qperpy * qperpy * s101 -
                    qperpz * qperpz * s101 - q1w * q1z * s111 +
                    qperpx * qperpy * s111 -
                    qperpw * qperpz * s111 + q1w * q1y * s121 +
                    qperpw * qperpy * s121 +
                    qperpx * qperpz * s121 +
                    q1x * (-(q1y * s011) - q1z * s021 + q1y * s111 + q1z * s121),
                kz: (-1.0 + q1y * q1y + q1z * q1z + qperpy * qperpy + qperpz * qperpz) * s002 +
                    q1w * q1z * s012 - qperpx * qperpy * s012 +
                    qperpw * qperpz * s012 - q1w * q1y * s022 -
                    qperpw * qperpy * s022 -
                    qperpx * qperpz * s022 + s102 -
                    q1y * q1y * s102 - q1z * q1z * s102 -
                    qperpy * qperpy * s102 -
                    qperpz * qperpz * s102 - q1w * q1z * s112 +
                    qperpx * qperpy * s112 -
                    qperpw * qperpz * s112 + q1w * q1y * s122 +
                    qperpw * qperpy * s122 +
                    qperpx * qperpz * s122 +
                    q1x * (-(q1y * s012) - q1z * s022 + q1y * s112 + q1z * s122),
            };
            at.c2[0] = DerivativeTerm {
                kc: 0.0,
                kx: -(qperpy * qperpy * s000) - qperpz * qperpz * s000 + qperpx * qperpy * s010 -
                    qperpw * qperpz * s010 +
                    qperpw * qperpy * s020 + qperpx * qperpz * s020 +
                    q1y * q1y * (s000 - s100) + q1z * q1z * (s000 - s100) +
                    qperpy * qperpy * s100 +
                    qperpz * qperpz * s100 - qperpx * qperpy * s110 +
                    qperpw * qperpz * s110 -
                    qperpw * qperpy * s120 - qperpx * qperpz * s120 +
                    2.0 * q1x * qperpy * s010 * theta -
                    2.0 * q1w * qperpz * s010 * theta +
                    2.0 * q1w * qperpy * s020 * theta +
                    2.0 * q1x * qperpz * s020 * theta +
                    q1y *
                    (q1x * (-s010 + s110) + q1w * (-s020 + s120) +
                     2.0 * (-2.0 * qperpy * s000 + qperpx * s010 + qperpw * s020) * theta) +
                    q1z *
                    (q1w * (s010 - s110) + q1x * (-s020 + s120) -
                     2.0 * (2.0 * qperpz * s000 + qperpw * s010 - qperpx * s020) * theta),
                ky: -(qperpy * qperpy * s001) - qperpz * qperpz * s001 + qperpx * qperpy * s011 -
                    qperpw * qperpz * s011 +
                    qperpw * qperpy * s021 + qperpx * qperpz * s021 +
                    q1y * q1y * (s001 - s101) + q1z * q1z * (s001 - s101) +
                    qperpy * qperpy * s101 +
                    qperpz * qperpz * s101 - qperpx * qperpy * s111 +
                    qperpw * qperpz * s111 -
                    qperpw * qperpy * s121 - qperpx * qperpz * s121 +
                    2.0 * q1x * qperpy * s011 * theta -
                    2.0 * q1w * qperpz * s011 * theta +
                    2.0 * q1w * qperpy * s021 * theta +
                    2.0 * q1x * qperpz * s021 * theta +
                    q1y *
                    (q1x * (-s011 + s111) + q1w * (-s021 + s121) +
                     2.0 * (-2.0 * qperpy * s001 + qperpx * s011 + qperpw * s021) * theta) +
                    q1z *
                    (q1w * (s011 - s111) + q1x * (-s021 + s121) -
                     2.0 * (2.0 * qperpz * s001 + qperpw * s011 - qperpx * s021) * theta),
                kz: -(qperpy * qperpy * s002) - qperpz * qperpz * s002 + qperpx * qperpy * s012 -
                    qperpw * qperpz * s012 +
                    qperpw * qperpy * s022 + qperpx * qperpz * s022 +
                    q1y * q1y * (s002 - s102) + q1z * q1z * (s002 - s102) +
                    qperpy * qperpy * s102 +
                    qperpz * qperpz * s102 - qperpx * qperpy * s112 +
                    qperpw * qperpz * s112 -
                    qperpw * qperpy * s122 - qperpx * qperpz * s122 +
                    2.0 * q1x * qperpy * s012 * theta -
                    2.0 * q1w * qperpz * s012 * theta +
                    2.0 * q1w * qperpy * s022 * theta +
                    2.0 * q1x * qperpz * s022 * theta +
                    q1y *
                    (q1x * (-s012 + s112) + q1w * (-s022 + s122) +
                     2.0 * (-2.0 * qperpy * s002 + qperpx * s012 + qperpw * s022) * theta) +
                    q1z *
                    (q1w * (s012 - s112) + q1x * (-s022 + s122) -
                     2.0 * (2.0 * qperpz * s002 + qperpw * s012 - qperpx * s022) * theta),
            };
            at.c3[0] = DerivativeTerm {
                kc: 0.0,
                kx: -2.0 *
                    (q1x * qperpy * s010 - q1w * qperpz * s010 + q1w * qperpy * s020 +
                     q1x * qperpz * s020 - q1x * qperpy * s110 +
                     q1w * qperpz * s110 - q1w * qperpy * s120 -
                     q1x * qperpz * s120 +
                     q1y *
                     (-2.0 * qperpy * s000 + qperpx * s010 + qperpw * s020 +
                      2.0 * qperpy * s100 - qperpx * s110 - qperpw * s120) +
                     q1z *
                     (-2.0 * qperpz * s000 - qperpw * s010 + qperpx * s020 +
                      2.0 * qperpz * s100 + qperpw * s110 - qperpx * s120)) *
                    theta,
                ky: -2.0 *
                    (q1x * qperpy * s011 - q1w * qperpz * s011 + q1w * qperpy * s021 +
                     q1x * qperpz * s021 - q1x * qperpy * s111 +
                     q1w * qperpz * s111 - q1w * qperpy * s121 -
                     q1x * qperpz * s121 +
                     q1y *
                     (-2.0 * qperpy * s001 + qperpx * s011 + qperpw * s021 +
                      2.0 * qperpy * s101 - qperpx * s111 - qperpw * s121) +
                     q1z *
                     (-2.0 * qperpz * s001 - qperpw * s011 + qperpx * s021 +
                      2.0 * qperpz * s101 + qperpw * s111 - qperpx * s121)) *
                    theta,
                kz: -2.0 *
                    (q1x * qperpy * s012 - q1w * qperpz * s012 + q1w * qperpy * s022 +
                     q1x * qperpz * s022 - q1x * qperpy * s112 +
                     q1w * qperpz * s112 - q1w * qperpy * s122 -
                     q1x * qperpz * s122 +
                     q1y *
                     (-2.0 * qperpy * s002 + qperpx * s012 + qperpw * s022 +
                      2.0 * qperpy * s102 - qperpx * s112 - qperpw * s122) +
                     q1z *
                     (-2.0 * qperpz * s002 - qperpw * s012 + qperpx * s022 +
                      2.0 * qperpz * s102 + qperpw * s112 - qperpx * s122)) *
                    theta,
            };
            at.c4[0] = DerivativeTerm {
                kc: 0.0,
                kx: -(q1x * qperpy * s010) + q1w * qperpz * s010 - q1w * qperpy * s020 -
                    q1x * qperpz * s020 + q1x * qperpy * s110 -
                    q1w * qperpz * s110 + q1w * qperpy * s120 +
                    q1x * qperpz * s120 + 2.0 * q1y * q1y * s000 * theta +
                    2.0 * q1z * q1z * s000 * theta -
                    2.0 * qperpy * qperpy * s000 * theta -
                    2.0 * qperpz * qperpz * s000 * theta +
                    2.0 * qperpx * qperpy * s010 * theta -
                    2.0 * qperpw * qperpz * s010 * theta +
                    2.0 * qperpw * qperpy * s020 * theta +
                    2.0 * qperpx * qperpz * s020 * theta +
                    q1y *
                    (-(qperpx * s010) - qperpw * s020 + 2.0 * qperpy * (s000 - s100) +
                     qperpx * s110 + qperpw * s120 -
                     2.0 * q1x * s010 * theta - 2.0 * q1w * s020 * theta) +
                    q1z *
                    (2.0 * qperpz * s000 + qperpw * s010 - qperpx * s020 - 2.0 * qperpz * s100 -
                     qperpw * s110 + qperpx * s120 +
                     2.0 * q1w * s010 * theta - 2.0 * q1x * s020 * theta),
                ky: -(q1x * qperpy * s011) + q1w * qperpz * s011 - q1w * qperpy * s021 -
                    q1x * qperpz * s021 + q1x * qperpy * s111 -
                    q1w * qperpz * s111 + q1w * qperpy * s121 +
                    q1x * qperpz * s121 + 2.0 * q1y * q1y * s001 * theta +
                    2.0 * q1z * q1z * s001 * theta -
                    2.0 * qperpy * qperpy * s001 * theta -
                    2.0 * qperpz * qperpz * s001 * theta +
                    2.0 * qperpx * qperpy * s011 * theta -
                    2.0 * qperpw * qperpz * s011 * theta +
                    2.0 * qperpw * qperpy * s021 * theta +
                    2.0 * qperpx * qperpz * s021 * theta +
                    q1y *
                    (-(qperpx * s011) - qperpw * s021 + 2.0 * qperpy * (s001 - s101) +
                     qperpx * s111 + qperpw * s121 -
                     2.0 * q1x * s011 * theta - 2.0 * q1w * s021 * theta) +
                    q1z *
                    (2.0 * qperpz * s001 + qperpw * s011 - qperpx * s021 - 2.0 * qperpz * s101 -
                     qperpw * s111 + qperpx * s121 +
                     2.0 * q1w * s011 * theta - 2.0 * q1x * s021 * theta),
                kz: -(q1x * qperpy * s012) + q1w * qperpz * s012 - q1w * qperpy * s022 -
                    q1x * qperpz * s022 + q1x * qperpy * s112 -
                    q1w * qperpz * s112 + q1w * qperpy * s122 +
                    q1x * qperpz * s122 + 2.0 * q1y * q1y * s002 * theta +
                    2.0 * q1z * q1z * s002 * theta -
                    2.0 * qperpy * qperpy * s002 * theta -
                    2.0 * qperpz * qperpz * s002 * theta +
                    2.0 * qperpx * qperpy * s012 * theta -
                    2.0 * qperpw * qperpz * s012 * theta +
                    2.0 * qperpw * qperpy * s022 * theta +
                    2.0 * qperpx * qperpz * s022 * theta +
                    q1y *
                    (-(qperpx * s012) - qperpw * s022 + 2.0 * qperpy * (s002 - s102) +
                     qperpx * s112 + qperpw * s122 -
                     2.0 * q1x * s012 * theta - 2.0 * q1w * s022 * theta) +
                    q1z *
                    (2.0 * qperpz * s002 + qperpw * s012 - qperpx * s022 - 2.0 * qperpz * s102 -
                     qperpw * s112 + qperpx * s122 +
                     2.0 * q1w * s012 * theta - 2.0 * q1x * s022 * theta),
            };
            at.c5[0] = DerivativeTerm {
                kc: 0.0,
                kx: 2.0 *
                    (qperpy * qperpy * s000 + qperpz * qperpz * s000 - qperpx * qperpy * s010 +
                     qperpw * qperpz * s010 -
                     qperpw * qperpy * s020 - qperpx * qperpz * s020 -
                     qperpy * qperpy * s100 -
                     qperpz * qperpz * s100 + q1y * q1y * (-s000 + s100) +
                     q1z * q1z * (-s000 + s100) + qperpx * qperpy * s110 -
                     qperpw * qperpz * s110 +
                     q1y * (q1x * (s010 - s110) + q1w * (s020 - s120)) +
                     qperpw * qperpy * s120 + qperpx * qperpz * s120 +
                     q1z * (-(q1w * s010) + q1x * s020 + q1w * s110 - q1x * s120)) *
                    theta,
                ky: 2.0 *
                    (qperpy * qperpy * s001 + qperpz * qperpz * s001 - qperpx * qperpy * s011 +
                     qperpw * qperpz * s011 -
                     qperpw * qperpy * s021 - qperpx * qperpz * s021 -
                     qperpy * qperpy * s101 -
                     qperpz * qperpz * s101 + q1y * q1y * (-s001 + s101) +
                     q1z * q1z * (-s001 + s101) + qperpx * qperpy * s111 -
                     qperpw * qperpz * s111 +
                     q1y * (q1x * (s011 - s111) + q1w * (s021 - s121)) +
                     qperpw * qperpy * s121 + qperpx * qperpz * s121 +
                     q1z * (-(q1w * s011) + q1x * s021 + q1w * s111 - q1x * s121)) *
                    theta,
                kz: 2.0 *
                    (qperpy * qperpy * s002 + qperpz * qperpz * s002 - qperpx * qperpy * s012 +
                     qperpw * qperpz * s012 -
                     qperpw * qperpy * s022 - qperpx * qperpz * s022 -
                     qperpy * qperpy * s102 -
                     qperpz * qperpz * s102 + q1y * q1y * (-s002 + s102) +
                     q1z * q1z * (-s002 + s102) + qperpx * qperpy * s112 -
                     qperpw * qperpz * s112 +
                     q1y * (q1x * (s012 - s112) + q1w * (s022 - s122)) +
                     qperpw * qperpy * s122 + qperpx * qperpz * s122 +
                     q1z * (-(q1w * s012) + q1x * s022 + q1w * s112 - q1x * s122)) *
                    theta,
            };
            at.c1[1] = DerivativeTerm {
                kc: -t0y + t1y,
                kx: -(qperpx * qperpy * s000) - qperpw * qperpz * s000 - s010 + q1z * q1z * s010 +
                    qperpx * qperpx * s010 +
                    qperpz * qperpz * s010 - q1y * q1z * s020 +
                    qperpw * qperpx * s020 - qperpy * qperpz * s020 +
                    qperpx * qperpy * s100 + qperpw * qperpz * s100 +
                    q1w * q1z * (-s000 + s100) +
                    q1x * q1x * (s010 - s110) + s110 - q1z * q1z * s110 -
                    qperpx * qperpx * s110 - qperpz * qperpz * s110 +
                    q1x * (q1y * (-s000 + s100) + q1w * (s020 - s120)) +
                    q1y * q1z * s120 - qperpw * qperpx * s120 +
                    qperpy * qperpz * s120,
                ky: -(qperpx * qperpy * s001) - qperpw * qperpz * s001 - s011 + q1z * q1z * s011 +
                    qperpx * qperpx * s011 +
                    qperpz * qperpz * s011 - q1y * q1z * s021 +
                    qperpw * qperpx * s021 - qperpy * qperpz * s021 +
                    qperpx * qperpy * s101 + qperpw * qperpz * s101 +
                    q1w * q1z * (-s001 + s101) +
                    q1x * q1x * (s011 - s111) + s111 - q1z * q1z * s111 -
                    qperpx * qperpx * s111 - qperpz * qperpz * s111 +
                    q1x * (q1y * (-s001 + s101) + q1w * (s021 - s121)) +
                    q1y * q1z * s121 - qperpw * qperpx * s121 +
                    qperpy * qperpz * s121,
                kz: -(qperpx * qperpy * s002) - qperpw * qperpz * s002 - s012 + q1z * q1z * s012 +
                    qperpx * qperpx * s012 +
                    qperpz * qperpz * s012 - q1y * q1z * s022 +
                    qperpw * qperpx * s022 - qperpy * qperpz * s022 +
                    qperpx * qperpy * s102 + qperpw * qperpz * s102 +
                    q1w * q1z * (-s002 + s102) +
                    q1x * q1x * (s012 - s112) + s112 - q1z * q1z * s112 -
                    qperpx * qperpx * s112 - qperpz * qperpz * s112 +
                    q1x * (q1y * (-s002 + s102) + q1w * (s022 - s122)) +
                    q1y * q1z * s122 - qperpw * qperpx * s122 +
                    qperpy * qperpz * s122,
            };
            at.c2[1] = DerivativeTerm {
                kc: 0.0,
                kx: qperpx * qperpy * s000 + qperpw * qperpz * s000 + q1z * q1z * s010 -
                    qperpx * qperpx * s010 - qperpz * qperpz * s010 -
                    q1y * q1z * s020 - qperpw * qperpx * s020 +
                    qperpy * qperpz * s020 -
                    qperpx * qperpy * s100 - qperpw * qperpz * s100 +
                    q1x * q1x * (s010 - s110) - q1z * q1z * s110 +
                    qperpx * qperpx * s110 +
                    qperpz * qperpz * s110 +
                    q1y * q1z * s120 + qperpw * qperpx * s120 -
                    qperpy * qperpz * s120 + 2.0 * q1z * qperpw * s000 * theta +
                    2.0 * q1y * qperpx * s000 * theta -
                    4.0 * q1z * qperpz * s010 * theta +
                    2.0 * q1z * qperpy * s020 * theta +
                    2.0 * q1y * qperpz * s020 * theta +
                    q1x *
                    (q1w * s020 + q1y * (-s000 + s100) - q1w * s120 + 2.0 * qperpy * s000 * theta -
                     4.0 * qperpx * s010 * theta -
                     2.0 * qperpw * s020 * theta) +
                    q1w *
                    (-(q1z * s000) + q1z * s100 + 2.0 * qperpz * s000 * theta -
                     2.0 * qperpx * s020 * theta),
                ky: qperpx * qperpy * s001 + qperpw * qperpz * s001 + q1z * q1z * s011 -
                    qperpx * qperpx * s011 - qperpz * qperpz * s011 -
                    q1y * q1z * s021 - qperpw * qperpx * s021 +
                    qperpy * qperpz * s021 -
                    qperpx * qperpy * s101 - qperpw * qperpz * s101 +
                    q1x * q1x * (s011 - s111) - q1z * q1z * s111 +
                    qperpx * qperpx * s111 +
                    qperpz * qperpz * s111 +
                    q1y * q1z * s121 + qperpw * qperpx * s121 -
                    qperpy * qperpz * s121 + 2.0 * q1z * qperpw * s001 * theta +
                    2.0 * q1y * qperpx * s001 * theta -
                    4.0 * q1z * qperpz * s011 * theta +
                    2.0 * q1z * qperpy * s021 * theta +
                    2.0 * q1y * qperpz * s021 * theta +
                    q1x *
                    (q1w * s021 + q1y * (-s001 + s101) - q1w * s121 + 2.0 * qperpy * s001 * theta -
                     4.0 * qperpx * s011 * theta -
                     2.0 * qperpw * s021 * theta) +
                    q1w *
                    (-(q1z * s001) + q1z * s101 + 2.0 * qperpz * s001 * theta -
                     2.0 * qperpx * s021 * theta),
                kz: qperpx * qperpy * s002 + qperpw * qperpz * s002 + q1z * q1z * s012 -
                    qperpx * qperpx * s012 - qperpz * qperpz * s012 -
                    q1y * q1z * s022 - qperpw * qperpx * s022 +
                    qperpy * qperpz * s022 -
                    qperpx * qperpy * s102 - qperpw * qperpz * s102 +
                    q1x * q1x * (s012 - s112) - q1z * q1z * s112 +
                    qperpx * qperpx * s112 +
                    qperpz * qperpz * s112 +
                    q1y * q1z * s122 + qperpw * qperpx * s122 -
                    qperpy * qperpz * s122 + 2.0 * q1z * qperpw * s002 * theta +
                    2.0 * q1y * qperpx * s002 * theta -
                    4.0 * q1z * qperpz * s012 * theta +
                    2.0 * q1z * qperpy * s022 * theta +
                    2.0 * q1y * qperpz * s022 * theta +
                    q1x *
                    (q1w * s022 + q1y * (-s002 + s102) - q1w * s122 + 2.0 * qperpy * s002 * theta -
                     4.0 * qperpx * s012 * theta -
                     2.0 * qperpw * s022 * theta) +
                    q1w *
                    (-(q1z * s002) + q1z * s102 + 2.0 * qperpz * s002 * theta -
                     2.0 * qperpx * s022 * theta),
            };
            at.c3[1] =
                DerivativeTerm {
                    kc: 0.0,
                    kx: 2.0 *
                        (-(q1x * qperpy * s000) - q1w * qperpz * s000 + 2.0 * q1x * qperpx * s010 +
                         q1x * qperpw * s020 + q1w * qperpx * s020 +
                         q1x * qperpy * s100 + q1w * qperpz * s100 -
                         2.0 * q1x * qperpx * s110 -
                         q1x * qperpw * s120 - q1w * qperpx * s120 +
                         q1z *
                         (2.0 * qperpz * s010 - qperpy * s020 + qperpw * (-s000 + s100) -
                          2.0 * qperpz * s110 + qperpy * s120) +
                         q1y * (-(qperpx * s000) - qperpz * s020 + qperpx * s100 + qperpz * s120)) *
                        theta,
                    ky: 2.0 *
                        (-(q1x * qperpy * s001) - q1w * qperpz * s001 + 2.0 * q1x * qperpx * s011 +
                         q1x * qperpw * s021 + q1w * qperpx * s021 +
                         q1x * qperpy * s101 + q1w * qperpz * s101 -
                         2.0 * q1x * qperpx * s111 -
                         q1x * qperpw * s121 - q1w * qperpx * s121 +
                         q1z *
                         (2.0 * qperpz * s011 - qperpy * s021 + qperpw * (-s001 + s101) -
                          2.0 * qperpz * s111 + qperpy * s121) +
                         q1y * (-(qperpx * s001) - qperpz * s021 + qperpx * s101 + qperpz * s121)) *
                        theta,
                    kz: 2.0 *
                        (-(q1x * qperpy * s002) - q1w * qperpz * s002 + 2.0 * q1x * qperpx * s012 +
                         q1x * qperpw * s022 + q1w * qperpx * s022 +
                         q1x * qperpy * s102 + q1w * qperpz * s102 -
                         2.0 * q1x * qperpx * s112 -
                         q1x * qperpw * s122 - q1w * qperpx * s122 +
                         q1z *
                         (2.0 * qperpz * s012 - qperpy * s022 + qperpw * (-s002 + s102) -
                          2.0 * qperpz * s112 + qperpy * s122) +
                         q1y * (-(qperpx * s002) - qperpz * s022 + qperpx * s102 + qperpz * s122)) *
                        theta,
                };
            at.c4[1] = DerivativeTerm {
                kc: 0.0,
                kx: -(q1x * qperpy * s000) - q1w * qperpz * s000 + 2.0 * q1x * qperpx * s010 +
                    q1x * qperpw * s020 + q1w * qperpx * s020 +
                    q1x * qperpy * s100 +
                    q1w * qperpz * s100 - 2.0 * q1x * qperpx * s110 -
                    q1x * qperpw * s120 -
                    q1w * qperpx * s120 + 2.0 * qperpx * qperpy * s000 * theta +
                    2.0 * qperpw * qperpz * s000 * theta +
                    2.0 * q1x * q1x * s010 * theta +
                    2.0 * q1z * q1z * s010 * theta -
                    2.0 * qperpx * qperpx * s010 * theta -
                    2.0 * qperpz * qperpz * s010 * theta +
                    2.0 * q1w * q1x * s020 * theta -
                    2.0 * qperpw * qperpx * s020 * theta +
                    2.0 * qperpy * qperpz * s020 * theta +
                    q1y *
                    (-(qperpx * s000) - qperpz * s020 + qperpx * s100 + qperpz * s120 -
                     2.0 * q1x * s000 * theta) +
                    q1z *
                    (2.0 * qperpz * s010 - qperpy * s020 + qperpw * (-s000 + s100) -
                     2.0 * qperpz * s110 + qperpy * s120 -
                     2.0 * q1w * s000 * theta - 2.0 * q1y * s020 * theta),
                ky: -(q1x * qperpy * s001) - q1w * qperpz * s001 + 2.0 * q1x * qperpx * s011 +
                    q1x * qperpw * s021 + q1w * qperpx * s021 +
                    q1x * qperpy * s101 +
                    q1w * qperpz * s101 - 2.0 * q1x * qperpx * s111 -
                    q1x * qperpw * s121 -
                    q1w * qperpx * s121 + 2.0 * qperpx * qperpy * s001 * theta +
                    2.0 * qperpw * qperpz * s001 * theta +
                    2.0 * q1x * q1x * s011 * theta +
                    2.0 * q1z * q1z * s011 * theta -
                    2.0 * qperpx * qperpx * s011 * theta -
                    2.0 * qperpz * qperpz * s011 * theta +
                    2.0 * q1w * q1x * s021 * theta -
                    2.0 * qperpw * qperpx * s021 * theta +
                    2.0 * qperpy * qperpz * s021 * theta +
                    q1y *
                    (-(qperpx * s001) - qperpz * s021 + qperpx * s101 + qperpz * s121 -
                     2.0 * q1x * s001 * theta) +
                    q1z *
                    (2.0 * qperpz * s011 - qperpy * s021 + qperpw * (-s001 + s101) -
                     2.0 * qperpz * s111 + qperpy * s121 -
                     2.0 * q1w * s001 * theta - 2.0 * q1y * s021 * theta),
                kz: -(q1x * qperpy * s002) - q1w * qperpz * s002 + 2.0 * q1x * qperpx * s012 +
                    q1x * qperpw * s022 + q1w * qperpx * s022 +
                    q1x * qperpy * s102 +
                    q1w * qperpz * s102 - 2.0 * q1x * qperpx * s112 -
                    q1x * qperpw * s122 -
                    q1w * qperpx * s122 + 2.0 * qperpx * qperpy * s002 * theta +
                    2.0 * qperpw * qperpz * s002 * theta +
                    2.0 * q1x * q1x * s012 * theta +
                    2.0 * q1z * q1z * s012 * theta -
                    2.0 * qperpx * qperpx * s012 * theta -
                    2.0 * qperpz * qperpz * s012 * theta +
                    2.0 * q1w * q1x * s022 * theta -
                    2.0 * qperpw * qperpx * s022 * theta +
                    2.0 * qperpy * qperpz * s022 * theta +
                    q1y *
                    (-(qperpx * s002) - qperpz * s022 + qperpx * s102 + qperpz * s122 -
                     2.0 * q1x * s002 * theta) +
                    q1z *
                    (2.0 * qperpz * s012 - qperpy * s022 + qperpw * (-s002 + s102) -
                     2.0 * qperpz * s112 + qperpy * s122 -
                     2.0 * q1w * s002 * theta - 2.0 * q1y * s022 * theta),
            };
            at.c5[1] = DerivativeTerm {
                kc: 0.,
                kx: -2.0 *
                    (qperpx * qperpy * s000 + qperpw * qperpz * s000 + q1z * q1z * s010 -
                     qperpx * qperpx * s010 - qperpz * qperpz * s010 -
                     q1y * q1z * s020 - qperpw * qperpx * s020 +
                     qperpy * qperpz * s020 -
                     qperpx * qperpy * s100 - qperpw * qperpz * s100 +
                     q1w * q1z * (-s000 + s100) +
                     q1x * q1x * (s010 - s110) - q1z * q1z * s110 +
                     qperpx * qperpx * s110 + qperpz * qperpz * s110 +
                     q1x * (q1y * (-s000 + s100) + q1w * (s020 - s120)) +
                     q1y * q1z * s120 + qperpw * qperpx * s120 -
                     qperpy * qperpz * s120) * theta,
                ky: -2.0 *
                    (qperpx * qperpy * s001 + qperpw * qperpz * s001 + q1z * q1z * s011 -
                     qperpx * qperpx * s011 - qperpz * qperpz * s011 -
                     q1y * q1z * s021 - qperpw * qperpx * s021 +
                     qperpy * qperpz * s021 -
                     qperpx * qperpy * s101 - qperpw * qperpz * s101 +
                     q1w * q1z * (-s001 + s101) +
                     q1x * q1x * (s011 - s111) - q1z * q1z * s111 +
                     qperpx * qperpx * s111 + qperpz * qperpz * s111 +
                     q1x * (q1y * (-s001 + s101) + q1w * (s021 - s121)) +
                     q1y * q1z * s121 + qperpw * qperpx * s121 -
                     qperpy * qperpz * s121) * theta,
                kz: -2.0 *
                    (qperpx * qperpy * s002 + qperpw * qperpz * s002 + q1z * q1z * s012 -
                     qperpx * qperpx * s012 - qperpz * qperpz * s012 -
                     q1y * q1z * s022 - qperpw * qperpx * s022 +
                     qperpy * qperpz * s022 -
                     qperpx * qperpy * s102 - qperpw * qperpz * s102 +
                     q1w * q1z * (-s002 + s102) +
                     q1x * q1x * (s012 - s112) - q1z * q1z * s112 +
                     qperpx * qperpx * s112 + qperpz * qperpz * s112 +
                     q1x * (q1y * (-s002 + s102) + q1w * (s022 - s122)) +
                     q1y * q1z * s122 + qperpw * qperpx * s122 -
                     qperpy * qperpz * s122) * theta,
            };
            at.c1[2] = DerivativeTerm {
                kc: -t0z + t1z,
                kx: (qperpw * qperpy * s000 - qperpx * qperpz * s000 - q1y * q1z * s010 -
                     qperpw * qperpx * s010 - qperpy * qperpz * s010 - s020 +
                     q1y * q1y * s020 +
                     qperpx * qperpx * s020 + qperpy * qperpy * s020 -
                     qperpw * qperpy * s100 +
                     qperpx * qperpz * s100 + q1x * q1z * (-s000 + s100) +
                     q1y * q1z * s110 +
                     qperpw * qperpx * s110 + qperpy * qperpz * s110 +
                     q1w * (q1y * (s000 - s100) + q1x * (-s010 + s110)) +
                     q1x * q1x * (s020 - s120) + s120 - q1y * q1y * s120 -
                     qperpx * qperpx * s120 - qperpy * qperpy * s120),
                ky: (qperpw * qperpy * s001 - qperpx * qperpz * s001 - q1y * q1z * s011 -
                     qperpw * qperpx * s011 - qperpy * qperpz * s011 - s021 +
                     q1y * q1y * s021 +
                     qperpx * qperpx * s021 + qperpy * qperpy * s021 -
                     qperpw * qperpy * s101 +
                     qperpx * qperpz * s101 + q1x * q1z * (-s001 + s101) +
                     q1y * q1z * s111 +
                     qperpw * qperpx * s111 + qperpy * qperpz * s111 +
                     q1w * (q1y * (s001 - s101) + q1x * (-s011 + s111)) +
                     q1x * q1x * (s021 - s121) + s121 - q1y * q1y * s121 -
                     qperpx * qperpx * s121 - qperpy * qperpy * s121),
                kz: (qperpw * qperpy * s002 - qperpx * qperpz * s002 - q1y * q1z * s012 -
                     qperpw * qperpx * s012 - qperpy * qperpz * s012 - s022 +
                     q1y * q1y * s022 +
                     qperpx * qperpx * s022 + qperpy * qperpy * s022 -
                     qperpw * qperpy * s102 +
                     qperpx * qperpz * s102 + q1x * q1z * (-s002 + s102) +
                     q1y * q1z * s112 +
                     qperpw * qperpx * s112 + qperpy * qperpz * s112 +
                     q1w * (q1y * (s002 - s102) + q1x * (-s012 + s112)) +
                     q1x * q1x * (s022 - s122) + s122 - q1y * q1y * s122 -
                     qperpx * qperpx * s122 - qperpy * qperpy * s122),
            };
            at.c2[2] = DerivativeTerm {
                kc: 0.0,
                kx: (q1w * q1y * s000 - q1x * q1z * s000 - qperpw * qperpy * s000 +
                     qperpx * qperpz * s000 - q1w * q1x * s010 -
                     q1y * q1z * s010 + qperpw * qperpx * s010 +
                     qperpy * qperpz * s010 + q1x * q1x * s020 +
                     q1y * q1y * s020 - qperpx * qperpx * s020 -
                     qperpy * qperpy * s020 - q1w * q1y * s100 +
                     q1x * q1z * s100 + qperpw * qperpy * s100 -
                     qperpx * qperpz * s100 + q1w * q1x * s110 +
                     q1y * q1z * s110 - qperpw * qperpx * s110 -
                     qperpy * qperpz * s110 - q1x * q1x * s120 -
                     q1y * q1y * s120 + qperpx * qperpx * s120 +
                     qperpy * qperpy * s120 -
                     2.0 * q1y * qperpw * s000 * theta +
                     2.0 * q1z * qperpx * s000 * theta -
                     2.0 * q1w * qperpy * s000 * theta +
                     2.0 * q1x * qperpz * s000 * theta +
                     2.0 * q1x * qperpw * s010 * theta +
                     2.0 * q1w * qperpx * s010 * theta +
                     2.0 * q1z * qperpy * s010 * theta +
                     2.0 * q1y * qperpz * s010 * theta -
                     4.0 * q1x * qperpx * s020 * theta -
                     4.0 * q1y * qperpy * s020 * theta),
                ky: (q1w * q1y * s001 - q1x * q1z * s001 - qperpw * qperpy * s001 +
                     qperpx * qperpz * s001 - q1w * q1x * s011 -
                     q1y * q1z * s011 + qperpw * qperpx * s011 +
                     qperpy * qperpz * s011 + q1x * q1x * s021 +
                     q1y * q1y * s021 - qperpx * qperpx * s021 -
                     qperpy * qperpy * s021 - q1w * q1y * s101 +
                     q1x * q1z * s101 + qperpw * qperpy * s101 -
                     qperpx * qperpz * s101 + q1w * q1x * s111 +
                     q1y * q1z * s111 - qperpw * qperpx * s111 -
                     qperpy * qperpz * s111 - q1x * q1x * s121 -
                     q1y * q1y * s121 + qperpx * qperpx * s121 +
                     qperpy * qperpy * s121 -
                     2.0 * q1y * qperpw * s001 * theta +
                     2.0 * q1z * qperpx * s001 * theta -
                     2.0 * q1w * qperpy * s001 * theta +
                     2.0 * q1x * qperpz * s001 * theta +
                     2.0 * q1x * qperpw * s011 * theta +
                     2.0 * q1w * qperpx * s011 * theta +
                     2.0 * q1z * qperpy * s011 * theta +
                     2.0 * q1y * qperpz * s011 * theta -
                     4.0 * q1x * qperpx * s021 * theta -
                     4.0 * q1y * qperpy * s021 * theta),
                kz: (q1w * q1y * s002 - q1x * q1z * s002 - qperpw * qperpy * s002 +
                     qperpx * qperpz * s002 - q1w * q1x * s012 -
                     q1y * q1z * s012 + qperpw * qperpx * s012 +
                     qperpy * qperpz * s012 + q1x * q1x * s022 +
                     q1y * q1y * s022 - qperpx * qperpx * s022 -
                     qperpy * qperpy * s022 - q1w * q1y * s102 +
                     q1x * q1z * s102 + qperpw * qperpy * s102 -
                     qperpx * qperpz * s102 + q1w * q1x * s112 +
                     q1y * q1z * s112 - qperpw * qperpx * s112 -
                     qperpy * qperpz * s112 - q1x * q1x * s122 -
                     q1y * q1y * s122 + qperpx * qperpx * s122 +
                     qperpy * qperpy * s122 -
                     2.0 * q1y * qperpw * s002 * theta +
                     2.0 * q1z * qperpx * s002 * theta -
                     2.0 * q1w * qperpy * s002 * theta +
                     2.0 * q1x * qperpz * s002 * theta +
                     2.0 * q1x * qperpw * s012 * theta +
                     2.0 * q1w * qperpx * s012 * theta +
                     2.0 * q1z * qperpy * s012 * theta +
                     2.0 * q1y * qperpz * s012 * theta -
                     4.0 * q1x * qperpx * s022 * theta -
                     4.0 * q1y * qperpy * s022 * theta),
            };
            at.c3[2] =
                DerivativeTerm {
                    kc: 0.0,
                    kx: -2.0 *
                        (-(q1w * qperpy * s000) + q1x * qperpz * s000 + q1x * qperpw * s010 +
                         q1w * qperpx * s010 - 2.0 * q1x * qperpx * s020 +
                         q1w * qperpy * s100 - q1x * qperpz * s100 -
                         q1x * qperpw * s110 - q1w * qperpx * s110 +
                         q1z * (qperpx * s000 + qperpy * s010 - qperpx * s100 - qperpy * s110) +
                         2.0 * q1x * qperpx * s120 +
                         q1y *
                         (qperpz * s010 - 2.0 * qperpy * s020 + qperpw * (-s000 + s100) -
                          qperpz * s110 + 2.0 * qperpy * s120)) * theta,
                    ky: -2.0 *
                        (-(q1w * qperpy * s001) + q1x * qperpz * s001 + q1x * qperpw * s011 +
                         q1w * qperpx * s011 - 2.0 * q1x * qperpx * s021 +
                         q1w * qperpy * s101 - q1x * qperpz * s101 -
                         q1x * qperpw * s111 - q1w * qperpx * s111 +
                         q1z * (qperpx * s001 + qperpy * s011 - qperpx * s101 - qperpy * s111) +
                         2.0 * q1x * qperpx * s121 +
                         q1y *
                         (qperpz * s011 - 2.0 * qperpy * s021 + qperpw * (-s001 + s101) -
                          qperpz * s111 + 2.0 * qperpy * s121)) * theta,
                    kz: -2.0 *
                        (-(q1w * qperpy * s002) + q1x * qperpz * s002 + q1x * qperpw * s012 +
                         q1w * qperpx * s012 - 2.0 * q1x * qperpx * s022 +
                         q1w * qperpy * s102 - q1x * qperpz * s102 -
                         q1x * qperpw * s112 - q1w * qperpx * s112 +
                         q1z * (qperpx * s002 + qperpy * s012 - qperpx * s102 - qperpy * s112) +
                         2.0 * q1x * qperpx * s122 +
                         q1y *
                         (qperpz * s012 - 2.0 * qperpy * s022 + qperpw * (-s002 + s102) -
                          qperpz * s112 + 2.0 * qperpy * s122)) * theta,
                };
            at.c4[2] = DerivativeTerm {
                kc: 0.0,
                kx: q1w * qperpy * s000 - q1x * qperpz * s000 - q1x * qperpw * s010 -
                    q1w * qperpx * s010 + 2.0 * q1x * qperpx * s020 -
                    q1w * qperpy * s100 + q1x * qperpz * s100 +
                    q1x * qperpw * s110 + q1w * qperpx * s110 -
                    2.0 * q1x * qperpx * s120 -
                    2.0 * qperpw * qperpy * s000 * theta +
                    2.0 * qperpx * qperpz * s000 * theta -
                    2.0 * q1w * q1x * s010 * theta +
                    2.0 * qperpw * qperpx * s010 * theta +
                    2.0 * qperpy * qperpz * s010 * theta +
                    2.0 * q1x * q1x * s020 * theta +
                    2.0 * q1y * q1y * s020 * theta -
                    2.0 * qperpx * qperpx * s020 * theta -
                    2.0 * qperpy * qperpy * s020 * theta +
                    q1z *
                    (-(qperpx * s000) - qperpy * s010 + qperpx * s100 + qperpy * s110 -
                     2.0 * q1x * s000 * theta) +
                    q1y *
                    (-(qperpz * s010) + 2.0 * qperpy * s020 + qperpw * (s000 - s100) +
                     qperpz * s110 - 2.0 * qperpy * s120 +
                     2.0 * q1w * s000 * theta - 2.0 * q1z * s010 * theta),
                ky: q1w * qperpy * s001 - q1x * qperpz * s001 - q1x * qperpw * s011 -
                    q1w * qperpx * s011 + 2.0 * q1x * qperpx * s021 -
                    q1w * qperpy * s101 + q1x * qperpz * s101 +
                    q1x * qperpw * s111 + q1w * qperpx * s111 -
                    2.0 * q1x * qperpx * s121 -
                    2.0 * qperpw * qperpy * s001 * theta +
                    2.0 * qperpx * qperpz * s001 * theta -
                    2.0 * q1w * q1x * s011 * theta +
                    2.0 * qperpw * qperpx * s011 * theta +
                    2.0 * qperpy * qperpz * s011 * theta +
                    2.0 * q1x * q1x * s021 * theta +
                    2.0 * q1y * q1y * s021 * theta -
                    2.0 * qperpx * qperpx * s021 * theta -
                    2.0 * qperpy * qperpy * s021 * theta +
                    q1z *
                    (-(qperpx * s001) - qperpy * s011 + qperpx * s101 + qperpy * s111 -
                     2.0 * q1x * s001 * theta) +
                    q1y *
                    (-(qperpz * s011) + 2.0 * qperpy * s021 + qperpw * (s001 - s101) +
                     qperpz * s111 - 2.0 * qperpy * s121 +
                     2.0 * q1w * s001 * theta - 2.0 * q1z * s011 * theta),
                kz: q1w * qperpy * s002 - q1x * qperpz * s002 - q1x * qperpw * s012 -
                    q1w * qperpx * s012 + 2.0 * q1x * qperpx * s022 -
                    q1w * qperpy * s102 + q1x * qperpz * s102 +
                    q1x * qperpw * s112 + q1w * qperpx * s112 -
                    2.0 * q1x * qperpx * s122 -
                    2.0 * qperpw * qperpy * s002 * theta +
                    2.0 * qperpx * qperpz * s002 * theta -
                    2.0 * q1w * q1x * s012 * theta +
                    2.0 * qperpw * qperpx * s012 * theta +
                    2.0 * qperpy * qperpz * s012 * theta +
                    2.0 * q1x * q1x * s022 * theta +
                    2.0 * q1y * q1y * s022 * theta -
                    2.0 * qperpx * qperpx * s022 * theta -
                    2.0 * qperpy * qperpy * s022 * theta +
                    q1z *
                    (-(qperpx * s002) - qperpy * s012 + qperpx * s102 + qperpy * s112 -
                     2.0 * q1x * s002 * theta) +
                    q1y *
                    (-(qperpz * s012) + 2.0 * qperpy * s022 + qperpw * (s002 - s102) +
                     qperpz * s112 - 2.0 * qperpy * s122 +
                     2.0 * q1w * s002 * theta - 2.0 * q1z * s012 * theta),
            };
            at.c5[2] = DerivativeTerm {
                kc: 0.,
                kx: 2.0 *
                    (qperpw * qperpy * s000 - qperpx * qperpz * s000 + q1y * q1z * s010 -
                     qperpw * qperpx * s010 - qperpy * qperpz * s010 -
                     q1y * q1y * s020 + qperpx * qperpx * s020 +
                     qperpy * qperpy * s020 +
                     q1x * q1z * (s000 - s100) - qperpw * qperpy * s100 +
                     qperpx * qperpz * s100 +
                     q1w * (q1y * (-s000 + s100) + q1x * (s010 - s110)) -
                     q1y * q1z * s110 + qperpw * qperpx * s110 +
                     qperpy * qperpz * s110 + q1y * q1y * s120 -
                     qperpx * qperpx * s120 - qperpy * qperpy * s120 +
                     q1x * q1x * (-s020 + s120)) * theta,
                ky: 2.0 *
                    (qperpw * qperpy * s001 - qperpx * qperpz * s001 + q1y * q1z * s011 -
                     qperpw * qperpx * s011 - qperpy * qperpz * s011 -
                     q1y * q1y * s021 + qperpx * qperpx * s021 +
                     qperpy * qperpy * s021 +
                     q1x * q1z * (s001 - s101) - qperpw * qperpy * s101 +
                     qperpx * qperpz * s101 +
                     q1w * (q1y * (-s001 + s101) + q1x * (s011 - s111)) -
                     q1y * q1z * s111 + qperpw * qperpx * s111 +
                     qperpy * qperpz * s111 + q1y * q1y * s121 -
                     qperpx * qperpx * s121 - qperpy * qperpy * s121 +
                     q1x * q1x * (-s021 + s121)) * theta,
                kz: 2.0 *
                    (qperpw * qperpy * s002 - qperpx * qperpz * s002 + q1y * q1z * s012 -
                     qperpw * qperpx * s012 - qperpy * qperpz * s012 -
                     q1y * q1y * s022 + qperpx * qperpx * s022 +
                     qperpy * qperpy * s022 +
                     q1x * q1z * (s002 - s102) - qperpw * qperpy * s102 +
                     qperpx * qperpz * s102 +
                     q1w * (q1y * (-s002 + s102) + q1x * (s012 - s112)) -
                     q1y * q1z * s112 + qperpw * qperpx * s112 +
                     qperpy * qperpz * s112 + q1y * q1y * s122 -
                     qperpx * qperpx * s122 - qperpy * qperpy * s122 +
                     q1x * q1x * (-s022 + s122)) * theta,
            };
        }
        at
    }
    pub fn decompose(m: &Matrix4x4, t: &mut Vector3f, rquat: &mut Quaternion, s: &mut Matrix4x4) {
        // extract translation from transformation matrix
        t.x = m.m[0][3];
        t.y = m.m[1][3];
        t.z = m.m[2][3];
        // compute new transformation matrix _m_ without translation
        let mut matrix: Matrix4x4 = m.clone();
        for i in 0..3 {
            matrix.m[i][3] = 0.0;
            matrix.m[3][i] = 0.0;
        }
        matrix.m[3][3] = 1.0;
        // extract rotation _r_ from transformation matrix
        let mut norm: Float;
        let mut count: u8 = 0;
        let mut r: Matrix4x4 = matrix.clone();
        loop {
            // compute next matrix _rnext_ in series
            let mut rnext: Matrix4x4 = Matrix4x4::default();
            let rit: Matrix4x4 = Matrix4x4::inverse(Matrix4x4::transpose(r));
            for i in 0..4 {
                for j in 0..4 {
                    rnext.m[i][j] = 0.5 * (r.m[i][j] + rit.m[i][j]);
                }
            }
            // compute norm of difference between _r_ and _rnext_
            norm = 0.0;
            for i in 0..3 {
                let n: Float = (r.m[i][0] - rnext.m[i][0]).abs() +
                               (r.m[i][1] - rnext.m[i][1]).abs() +
                               (r.m[i][2] - rnext.m[i][2]).abs();
                norm = norm.max(n);
            }
            r = rnext.clone();
            count += 1;
            if count >= 100 || norm <= 0.0001 {
                break;
            }
        }
        // XXX TODO FIXME deal with flip...
        let transform: Transform = Transform {
            m: r.clone(),
            m_inv: Matrix4x4::inverse(r.clone()),
        };
        *rquat = Quaternion::new(transform);

        // compute scale _S_ using rotation and original matrix
        *s = mtx_mul(Matrix4x4::inverse(r), *m);
    }
    pub fn transform_ray(&self, r: &mut Ray) {
        // if !self.actually_animated || self.r.time <= self.start_time {
        // } else if ...
        // TODO: above
        self.start_transform.transform_ray(r)
    }
}

// see quaternion.h

#[derive(Debug,Copy,Clone)]
pub struct Quaternion {
    pub v: Vector3f,
    pub w: Float,
}

impl Default for Quaternion {
    fn default() -> Self {
        Quaternion {
            v: Vector3f::default(),
            w: 1.0,
        }
    }
}

impl Sub for Quaternion {
    type Output = Quaternion;
    fn sub(self, rhs: Quaternion) -> Quaternion {
        Quaternion {
            v: self.v - rhs.v,
            w: self.w - rhs.w,
        }
    }
}

impl Mul<Float> for Quaternion {
    type Output = Quaternion;
    fn mul(self, rhs: Float) -> Quaternion {
        Quaternion {
            v: self.v * rhs,
            w: self.w * rhs,
        }
    }
}

impl Div<Float> for Quaternion {
    type Output = Quaternion;
    fn div(self, rhs: Float) -> Quaternion {
        Quaternion {
            v: self.v / rhs,
            w: self.w / rhs,
        }
    }
}

impl Neg for Quaternion {
    type Output = Quaternion;
    fn neg(self) -> Quaternion {
        Quaternion {
            v: -self.v,
            w: -self.w,
        }
    }
}

impl Quaternion {
    pub fn new(t: Transform) -> Self {
        let m: Matrix4x4 = t.m.clone();
        let trace: Float = m.m[0][0] + m.m[1][1] + m.m[2][2];
        // if (trace > 0.f) {
        if trace > 0.0 {
            // compute w from matrix trace, then xyz
            // 4w^2 = m[0][0] + m[1][1] + m[2][2] + m[3][3] (but m[3][3] == 1)
            let mut s: Float = (trace + 1.0).sqrt();
            let w: Float = s / 2.0;
            s = 0.5 / s;
            Quaternion {
                v: Vector3f {
                    x: (m.m[2][1] - m.m[1][2]) * s,
                    y: (m.m[0][2] - m.m[2][0]) * s,
                    z: (m.m[1][0] - m.m[0][1]) * s,
                },
                w: w,
            }
        } else {
            // compute largest of $x$, $y$, or $z$, then remaining components
            let nxt: [usize; 3] = [1, 2, 0];
            let mut q: [Float; 3] = [0.0; 3];
            let mut i: usize = 0;
            if m.m[1][1] > m.m[0][0] {
                i = 1;
            }
            if m.m[2][2] > m.m[i][i] {
                i = 2;
            }
            let j = nxt[i];
            let k = nxt[j];
            let mut s: Float = ((m.m[i][i] - (m.m[j][j] + m.m[k][k])) + 1.0).sqrt();
            q[i] = s * 0.5;
            if s != 0.0 {
                s = 0.5 / s;
            }
            let w: Float = (m.m[k][j] - m.m[j][k]) * s;
            q[j] = (m.m[j][i] + m.m[i][j]) * s;
            q[k] = (m.m[k][i] + m.m[i][k]) * s;
            Quaternion {
                v: Vector3f {
                    x: q[0],
                    y: q[1],
                    z: q[2],
                },
                w: w,
            }
        }
    }
}

/// The inner product of two quaterions.
pub fn quat_dot(q1: Quaternion, q2: Quaternion) -> Float {
    vec3_dot_vec3(q1.v, q2.v) + q1.w * q2.w
}

/// A quaternion can be normalized by dividing by its length.
pub fn quat_normalize(q: Quaternion) -> Quaternion {
    q / quat_dot(q, q)
}

// see interaction.h

pub trait Interaction {
    fn is_surface_interaction(&self) -> bool;
    fn is_medium_interaction(&self) -> bool;
    fn spawn_ray(&self, d: Vector3f) -> Ray;
    fn get_p(&self) -> Point3f;
    fn get_time(&self) -> Float;
    fn get_p_error(&self) -> Vector3f;
    fn get_wo(&self) -> Vector3f;
    fn get_n(&self) -> Normal3f;
}

#[derive(Debug,Default,Copy,Clone)]
pub struct InteractionCommon {
    // Interaction Public Data
    pub p: Point3f,
    pub time: Float,
    pub p_error: Vector3f,
    pub wo: Vector3f,
    pub n: Normal3f,
}

impl InteractionCommon {
    pub fn spawn_ray_to(&self, it: InteractionCommon) -> Ray {
        let origin: Point3f = pnt3_offset_ray_origin(self.p, self.p_error, self.n, it.p - self.p);
        let target: Point3f = pnt3_offset_ray_origin(it.p, it.p_error, it.n, origin - it.p);
        let d: Vector3f = target - origin;
        Ray {
            o: origin,
            d: d,
            t_max: 1.0 - SHADOW_EPSILON,
            time: self.time,
            differential: None,
        }
    }
}

#[derive(Debug,Default,Copy,Clone)]
pub struct Shading {
    pub n: Normal3f,
    pub dpdu: Vector3f,
    pub dpdv: Vector3f,
    pub dndu: Normal3f,
    pub dndv: Normal3f,
}

#[derive(Default,Clone)]
pub struct SurfaceInteraction<'a, 'b> {
    // Interaction Public Data
    pub p: Point3f,
    pub time: Float,
    pub p_error: Vector3f,
    pub wo: Vector3f,
    pub n: Normal3f,
    // TODO: MediumInterface mediumInterface;
    // SurfaceInteraction Public Data
    pub uv: Point2f,
    pub dpdu: Vector3f,
    pub dpdv: Vector3f,
    pub dndu: Normal3f,
    pub dndv: Normal3f,
    pub dpdx: Vector3f,
    pub dpdy: Vector3f,
    pub dudx: Float,
    pub dvdx: Float,
    pub dudy: Float,
    pub dvdy: Float,
    pub primitive: Option<&'a GeometricPrimitive>,
    pub shading: Shading,
    pub bsdf: Option<Arc<Bsdf>>,
    pub shape: Option<&'b Shape>,
}

impl<'a, 'b> SurfaceInteraction<'a, 'b> {
    pub fn new(p: Point3f,
               p_error: Vector3f,
               uv: Point2f,
               wo: Vector3f,
               dpdu: Vector3f,
               dpdv: Vector3f,
               dndu: Normal3f,
               dndv: Normal3f,
               time: Float,
               sh: Option<&'b Shape>)
               -> Self {
        let nv: Vector3f = vec3_normalize(vec3_cross_vec3(dpdu, dpdv));
        // TODO: Adjust normal based on orientation and handedness
        let n: Normal3f = Normal3f {
            x: nv.x,
            y: nv.y,
            z: nv.z,
        };
        // initialize shading geometry from true geometry
        let shading: Shading = Shading {
            n: n,
            dpdu: Vector3f::from(dpdu),
            dpdv: Vector3f::from(dpdv),
            dndu: dndu,
            dndv: dndv,
        };
        SurfaceInteraction {
            p: p,
            time: time,
            p_error: p_error,
            wo: vec3_normalize(wo),
            n: n,
            uv: uv,
            dpdu: dpdu,
            dpdv: dpdv,
            dndu: dndu,
            dndv: dndv,
            dpdx: Vector3f::default(),
            dpdy: Vector3f::default(),
            dudx: 0.0 as Float,
            dvdx: 0.0 as Float,
            dudy: 0.0 as Float,
            dvdy: 0.0 as Float,
            primitive: None,
            shading: shading,
            bsdf: None,
            shape: sh,
        }
    }
    pub fn set_shading_geometry(&mut self,
                                dpdus: Vector3f,
                                dpdvs: Vector3f,
                                dndus: Normal3f,
                                dndvs: Normal3f,
                                orientation_is_authoritative: bool) {
        // compute _shading.n_ for _SurfaceInteraction_
        self.shading.n = nrm_normalize(Normal3f::from(vec3_cross_vec3(dpdus, dpdvs)));
        if let Some(shape) = self.shape {
            if shape.get_reverse_orientation() ^ shape.get_transform_swaps_handedness() {
                self.shading.n = -self.shading.n;
            }
        }
        if orientation_is_authoritative {
            self.n = nrm_faceforward_nrm(self.n, self.shading.n);
        } else {
            self.shading.n = nrm_faceforward_nrm(self.shading.n, self.n);
        }
        // initialize _shading_ partial derivative values
        self.shading.dpdu = dpdus;
        self.shading.dpdv = dpdvs;
        self.shading.dndu = dndus;
        self.shading.dndv = dndvs;
    }
    pub fn compute_scattering_functions(&mut self,
                                        ray: &Ray,
                                        // arena: &mut Arena,
                                        allow_multiple_lobes: bool,
                                        mode: TransportMode) {
        self.compute_differentials(ray);
        if let Some(primitive) = self.primitive {
            primitive.compute_scattering_functions(self, // arena,
                                                   mode,
                                                   allow_multiple_lobes);
        }
    }
    pub fn compute_differentials(&mut self, ray: &Ray) {
        if let Some(ref diff) = ray.differential {
            // estimate screen space change in $\pt{}$ and $(u,v)$

            // compute auxiliary intersection points with plane
            let d: Float = vec3_dot_vec3(Vector3f::from(self.n),
                                         Vector3f {
                                             x: self.p.x,
                                             y: self.p.y,
                                             z: self.p.z,
                                         });
            let tx: Float =
                -(vec3_dot_vec3(Vector3f::from(self.n), Vector3f::from(diff.rx_origin)) - d) /
                vec3_dot_vec3(Vector3f::from(self.n), diff.rx_direction);
            if tx.is_nan() {
                self.dudx = 0.0 as Float;
                self.dvdx = 0.0 as Float;
                self.dudy = 0.0 as Float;
                self.dvdy = 0.0 as Float;
                self.dpdx = Vector3f::default();
                self.dpdy = Vector3f::default();
            } else {
                let px: Point3f = diff.rx_origin + diff.rx_direction * tx;
                let ty: Float =
                    -(vec3_dot_vec3(Vector3f::from(self.n), Vector3f::from(diff.ry_origin)) - d) /
                    vec3_dot_vec3(Vector3f::from(self.n), diff.ry_direction);
                if ty.is_nan() {
                    self.dudx = 0.0 as Float;
                    self.dvdx = 0.0 as Float;
                    self.dudy = 0.0 as Float;
                    self.dvdy = 0.0 as Float;
                    self.dpdx = Vector3f::default();
                    self.dpdy = Vector3f::default();
                } else {
                    let py: Point3f = diff.ry_origin + diff.ry_direction * ty;
                    self.dpdx = px - self.p;
                    self.dpdy = py - self.p;

                    // compute $(u,v)$ offsets at auxiliary points

                    // choose two dimensions to use for ray offset computation
                    let mut dim: [u8; 2] = [0_u8; 2];
                    if self.n.x.abs() > self.n.y.abs() && self.n.x.abs() > self.n.z.abs() {
                        dim[0] = 1;
                        dim[1] = 2;
                    } else if self.n.y.abs() > self.n.z.abs() {
                        dim[0] = 0;
                        dim[1] = 2;
                    } else {
                        dim[0] = 0;
                        dim[1] = 1;
                    }

                    // initialize _a_, _bx_, and _by_ matrices for offset computation
                    let a: [[Float; 2]; 2] = [[self.dpdu[dim[0]], self.dpdv[dim[0]]],
                                              [self.dpdu[dim[1]], self.dpdv[dim[1]]]];
                    let bx: [Float; 2] = [px[dim[0]] - self.p[dim[0]], px[dim[1]] - self.p[dim[1]]];
                    let by: [Float; 2] = [py[dim[0]] - self.p[dim[0]], py[dim[1]] - self.p[dim[1]]];
                    if !solve_linear_system_2x2(a, bx, &mut self.dudx, &mut self.dvdx) {
                        self.dudx = 0.0 as Float;
                        self.dvdx = 0.0 as Float;
                    }
                    if !solve_linear_system_2x2(a, by, &mut self.dudy, &mut self.dvdy) {
                        self.dudy = 0.0 as Float;
                        self.dvdy = 0.0 as Float;
                    }
                }
            }
        } else {
            self.dudx = 0.0 as Float;
            self.dvdx = 0.0 as Float;
            self.dudy = 0.0 as Float;
            self.dvdy = 0.0 as Float;
            self.dpdx = Vector3f::default();
            self.dpdy = Vector3f::default();
        }
    }
    pub fn le(&self, w: Vector3f) -> Spectrum {
        if let Some(primitive) = self.primitive {
            if let Some(area_light) = primitive.get_area_light() {
                // create InteractionCommon from self
                let interaction: InteractionCommon = InteractionCommon {
                    p: self.p,
                    time: self.time,
                    p_error: self.p_error,
                    wo: self.wo,
                    n: self.n,
                };
                return area_light.l(&interaction, w);
            }
        }
        Spectrum::default()
    }
}

impl<'a, 'b> Interaction for SurfaceInteraction<'a, 'b> {
    fn is_surface_interaction(&self) -> bool {
        self.n != Normal3f::default()
    }
    fn is_medium_interaction(&self) -> bool {
        !self.is_surface_interaction()
    }
    fn spawn_ray(&self, d: Vector3f) -> Ray {
        let o: Point3f = pnt3_offset_ray_origin(self.p, self.p_error, self.n, d);
        Ray {
            o: o,
            d: d,
            t_max: std::f32::INFINITY,
            time: self.time,
            differential: None,
        }
    }
    fn get_p(&self) -> Point3f {
        self.p.clone()
    }
    fn get_time(&self) -> Float {
        self.time
    }
    fn get_p_error(&self) -> Vector3f {
        self.p_error.clone()
    }
    fn get_wo(&self) -> Vector3f {
        self.wo.clone()
    }
    fn get_n(&self) -> Normal3f {
        self.n.clone()
    }
}

// see shape.h

pub trait Shape {
    fn object_bound(&self) -> Bounds3f;
    fn world_bound(&self) -> Bounds3f;
    fn intersect(&self, r: &Ray) -> Option<(SurfaceInteraction, Float)>;
    fn intersect_p(&self, r: &Ray) -> bool;
    fn get_reverse_orientation(&self) -> bool;
    fn get_transform_swaps_handedness(&self) -> bool;
    fn area(&self) -> Float;
    fn sample(&self, u: Point2f, pdf: &mut Float) -> InteractionCommon;
    fn sample_with_ref_point(&self, iref: &InteractionCommon, u: Point2f, pdf: &mut Float) -> InteractionCommon;
    fn pdf(&self, iref: &Interaction, wi: Vector3f) -> Float;
}

// see primitive.h

pub trait Primitive {
    fn world_bound(&self) -> Bounds3f;
    fn intersect(&self, ray: &mut Ray) -> Option<SurfaceInteraction>;
    fn intersect_p(&self, r: &Ray) -> bool;
    fn get_area_light(&self) -> Option<Arc<AreaLight + Send + Sync>>;
    fn get_material(&self) -> Option<Arc<Material + Send + Sync>>;
    fn compute_scattering_functions(&self,
                                    isect: &mut SurfaceInteraction,
                                    // arena: &mut Arena,
                                    mode: TransportMode,
                                    allow_multiple_lobes: bool) {
        if let Some(ref material) = self.get_material() {
            material.compute_scattering_functions(isect, // arena,
                                                  mode,
                                                  allow_multiple_lobes);
        }
        assert!(nrm_dot_nrm(isect.n, isect.shading.n) > 0.0);
    }
}

pub struct GeometricPrimitive {
    pub shape: Arc<Shape + Send + Sync>,
    pub material: Option<Arc<Material + Send + Sync>>,
    pub area_light: Option<Arc<AreaLight + Send + Sync>>,
    // TODO: MediumInterface mediumInterface;
}

impl GeometricPrimitive {
    pub fn new(shape: Arc<Shape + Send + Sync>,
               material: Arc<Material + Send + Sync>,
               area_light: Option<Arc<AreaLight + Send + Sync>>) -> Self {
        if let Some(area_light) = area_light {
            GeometricPrimitive {
                shape: shape,
                material: Some(material),
                area_light: Some(area_light),
            }
        } else {
            GeometricPrimitive {
                shape: shape,
                material: Some(material),
                area_light: None,
            }
        }
    }
}

impl Primitive for GeometricPrimitive {
    fn world_bound(&self) -> Bounds3f {
        self.shape.world_bound()
    }
    fn intersect(&self, ray: &mut Ray) -> Option<SurfaceInteraction> {
        self.shape.intersect(ray).map(|(mut isect, t_hit)| {
            isect.primitive = Some(self.clone());
            ray.t_max = t_hit;
            isect
        })
    }
    fn intersect_p(&self, r: &Ray) -> bool {
        self.shape.intersect_p(r)
    }
    fn get_material(&self) -> Option<Arc<Material + Send + Sync>> {
        if let Some(ref material) = self.material {
            Some(material.clone())
        } else {
            None
        }
    }
    fn get_area_light(&self) -> Option<Arc<AreaLight + Send + Sync>> {
        if let Some(ref area_light) = self.area_light {
            Some(area_light.clone())
        } else {
            None
        }
    }
}

// see cylinder.h

#[derive(Clone)]
pub struct Cylinder {
    pub radius: Float,
    pub z_min: Float,
    pub z_max: Float,
    pub phi_max: Float,
    // inherited from class Shape (see shape.h)
    object_to_world: Transform,
    world_to_object: Transform,
    reverse_orientation: bool,
    transform_swaps_handedness: bool,
    pub material: Option<Arc<Material + Send + Sync>>,
}

impl Default for Cylinder {
    fn default() -> Self {
        Cylinder {
            // Shape
            object_to_world: Transform::default(),
            world_to_object: Transform::default(),
            reverse_orientation: false,
            transform_swaps_handedness: false,
            // Cylinder
            radius: 1.0,
            z_min: -1.0,
            z_max: 1.0,
            phi_max: radians(360.0),
            material: None,
        }
    }
}

impl Cylinder {
    pub fn new(object_to_world: Transform,
               world_to_object: Transform,
               reverse_orientation: bool,
               radius: Float,
               z_min: Float,
               z_max: Float,
               phi_max: Float)
               -> Self {
        Cylinder {
            // Shape
            object_to_world: object_to_world,
            world_to_object: world_to_object,
            reverse_orientation: reverse_orientation,
            transform_swaps_handedness: false,
            // Cylinder
            radius: radius,
            z_min: z_min.min(z_max),
            z_max: z_min.max(z_max),
            phi_max: radians(clamp_t(phi_max, 0.0, 360.0)),
            material: None,
        }
    }
}

impl Shape for Cylinder {
    fn object_bound(&self) -> Bounds3f {
        Bounds3f {
            p_min: Point3f {
                x: -self.radius,
                y: -self.radius,
                z: self.z_min,
            },
            p_max: Point3f {
                x: self.radius,
                y: self.radius,
                z: self.z_max,
            },
        }
    }
    fn world_bound(&self) -> Bounds3f {
        // in C++: Bounds3f Shape::WorldBound() const { return (*ObjectToWorld)(ObjectBound()); }
        self.object_to_world.transform_bounds(self.object_bound())
    }
    fn intersect(&self, r: &Ray) -> Option<(SurfaceInteraction, Float)> {
        // TODO: ProfilePhase p(Prof::ShapeIntersect);
        // transform _Ray_ to object space
        let mut o_err: Vector3f = Vector3f::default();
        let mut d_err: Vector3f = Vector3f::default();
        let ray: Ray = self.world_to_object.transform_ray_with_error(r, &mut o_err, &mut d_err);

        // compute quadratic cylinder coefficients

        // initialize _EFloat_ ray coordinate values
        let ox = EFloat::new(ray.o.x as f32, o_err.x as f32);
        let oy = EFloat::new(ray.o.y as f32, o_err.y as f32);
        // let oz = EFloat::new(ray.o.z as f32, o_err.z as f32);
        let dx = EFloat::new(ray.d.x as f32, d_err.x as f32);
        let dy = EFloat::new(ray.d.y as f32, d_err.y as f32);
        // let dz = EFloat::new(ray.d.z as f32, d_err.z as f32);
        let a: EFloat = dx * dx + dy * dy;
        let b: EFloat = (dx * ox + dy * oy) * 2.0f32;
        let c: EFloat = ox * ox + oy * oy -
                        EFloat::new(self.radius as f32, 0.0) * EFloat::new(self.radius as f32, 0.0);

        // Solve quadratic equation for _t_ values
        let mut t0: EFloat = EFloat::default();
        let mut t1: EFloat = EFloat::default();
        if !quadratic_efloat(a, b, c, &mut t0, &mut t1) {
            return None;
        }
        // check quadric shape _t0_ and _t1_ for nearest intersection
        if t0.upper_bound() > ray.t_max as f32 || t1.lower_bound() <= 0.0f32 {
            return None;
        }
        let mut t_shape_hit: EFloat = t0;
        if t_shape_hit.lower_bound() <= 0.0f32 {
            t_shape_hit = t1;
            if t_shape_hit.upper_bound() > ray.t_max as f32 {
                return None;
            }
        }
        // compute cylinder hit point and $\phi$
        let mut p_hit: Point3f = ray.position(t_shape_hit.v);
        // refine cylinder intersection point
        let hit_rad: Float = (p_hit.x * p_hit.x + p_hit.y * p_hit.y).sqrt();
        p_hit.x *= self.radius / hit_rad;
        p_hit.y *= self.radius / hit_rad;
        let mut phi: Float = p_hit.y.atan2(p_hit.x);
        if phi < 0.0 as Float {
            phi += 2.0 as Float * PI;
        }
        // test cylinder intersection against clipping parameters
        if p_hit.z < self.z_min || p_hit.z > self.z_max || phi > self.phi_max {
            if t_shape_hit == t1 {
                return None;
            }
            t_shape_hit = t1;
            if t1.upper_bound() > ray.t_max {
                return None;
            }
            // compute cylinder hit point and $\phi$
            p_hit = ray.position(t_shape_hit.v);

            // refine cylinder intersection point
            let hit_rad: Float = (p_hit.x * p_hit.x + p_hit.y * p_hit.y).sqrt();
            p_hit.x *= self.radius / hit_rad;
            p_hit.y *= self.radius / hit_rad;
            phi = p_hit.y.atan2(p_hit.x);
            if phi < 0.0 as Float {
                phi += 2.0 as Float * PI;
            }
            if p_hit.z < self.z_min || p_hit.z > self.z_max || phi > self.phi_max {
                return None;
            }
        }
        // find parametric representation of cylinder hit
        let u: Float = phi / self.phi_max;
        let v: Float = (p_hit.z - self.z_min) / (self.z_max - self.z_min);
        // Compute cylinder $\dpdu$ and $\dpdv$
        let dpdu: Vector3f = Vector3f {
            x: -self.phi_max * p_hit.y,
            y: self.phi_max * p_hit.x,
            z: 0.0,
        };
        let dpdv: Vector3f = Vector3f {
            x: 0.0,
            y: 0.0,
            z: self.z_max - self.z_min,
        };
        // compute cylinder $\dndu$ and $\dndv$
        let d2_p_duu: Vector3f = Vector3f {
            x: p_hit.x,
            y: p_hit.y,
            z: 0.0,
        } * -self.phi_max * self.phi_max;
        let d2_p_duv: Vector3f = Vector3f {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };
        let d2_p_dvv: Vector3f = Vector3f {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };
        // compute coefficients for fundamental forms
        let ec: Float = vec3_dot_vec3(dpdu, dpdu);
        let fc: Float = vec3_dot_vec3(dpdu, dpdv);
        let gc: Float = vec3_dot_vec3(dpdv, dpdv);
        let nc: Vector3f = vec3_normalize(vec3_cross_vec3(dpdu, dpdv));
        let el: Float = vec3_dot_vec3(nc, d2_p_duu);
        let fl: Float = vec3_dot_vec3(nc, d2_p_duv);
        let gl: Float = vec3_dot_vec3(nc, d2_p_dvv);
        // compute $\dndu$ and $\dndv$ from fundamental form coefficients
        let inv_egf2: Float = 1.0 / (ec * gc - fc * fc);
        let dndu = dpdu * (fl * fc - el * gc) * inv_egf2 + dpdv * (el * fc - fl * ec) * inv_egf2;
        let dndu = Normal3f {
            x: dndu.x,
            y: dndu.y,
            z: dndu.z,
        };
        let dndv = dpdu * (gl * fc - fl * gc) * inv_egf2 + dpdv * (fl * fc - gl * ec) * inv_egf2;
        let dndv = Normal3f {
            x: dndv.x,
            y: dndv.y,
            z: dndv.z,
        };
        // compute error bounds for cylinder intersection
        let p_error: Vector3f = Vector3f {
                x: p_hit.x,
                y: p_hit.y,
                z: 0.0,
            }
            .abs() * gamma(3_i32);
        // initialize _SurfaceInteraction_ from parametric information
        let uv_hit: Point2f = Point2f { x: u, y: v };
        let wo: Vector3f = -ray.d;
        let si: SurfaceInteraction =
            SurfaceInteraction::new(p_hit, p_error, uv_hit, wo, dpdu, dpdv, dndu, dndv, ray.time, None);
        let mut isect: SurfaceInteraction = self.object_to_world.transform_surface_interaction(&si);
        if let Some(_shape) = si.shape {
            isect.shape = si.shape;
        }
        if let Some(_primitive) = si.primitive {
            isect.primitive = si.primitive;
        }
        Some((isect, t_shape_hit.v as Float))
    }
    fn intersect_p(&self, r: &Ray) -> bool {
        // TODO: ProfilePhase p(Prof::ShapeIntersect);
        // transform _Ray_ to object space
        let mut o_err: Vector3f = Vector3f::default();
        let mut d_err: Vector3f = Vector3f::default();
        let ray: Ray = self.world_to_object.transform_ray_with_error(r, &mut o_err, &mut d_err);

        // compute quadratic cylinder coefficients

        // initialize _EFloat_ ray coordinate values
        let ox = EFloat::new(ray.o.x as f32, o_err.x as f32);
        let oy = EFloat::new(ray.o.y as f32, o_err.y as f32);
        // let oz = EFloat::new(ray.o.z as f32, o_err.z as f32);
        let dx = EFloat::new(ray.d.x as f32, d_err.x as f32);
        let dy = EFloat::new(ray.d.y as f32, d_err.y as f32);
        // let dz = EFloat::new(ray.d.z as f32, d_err.z as f32);
        let a: EFloat = dx * dx + dy * dy;
        let b: EFloat = (dx * ox + dy * oy) * 2.0f32;
        let c: EFloat = ox * ox + oy * oy -
                        EFloat::new(self.radius as f32, 0.0) * EFloat::new(self.radius as f32, 0.0);

        // Solve quadratic equation for _t_ values
        let mut t0: EFloat = EFloat::default();
        let mut t1: EFloat = EFloat::default();
        if !quadratic_efloat(a, b, c, &mut t0, &mut t1) {
            return false;
        }
        // check quadric shape _t0_ and _t1_ for nearest intersection
        if t0.upper_bound() > ray.t_max as f32 || t1.lower_bound() <= 0.0f32 {
            return false;
        }
        let mut t_shape_hit: EFloat = t0;
        if t_shape_hit.lower_bound() <= 0.0f32 {
            t_shape_hit = t1;
            if t_shape_hit.upper_bound() > ray.t_max as f32 {
                return false;
            }
        }
        // compute cylinder hit point and $\phi$
        let mut p_hit: Point3f = ray.position(t_shape_hit.v);
        // refine cylinder intersection point
        let hit_rad: Float = (p_hit.x * p_hit.x + p_hit.y * p_hit.y).sqrt();
        p_hit.x *= self.radius / hit_rad;
        p_hit.y *= self.radius / hit_rad;
        let mut phi: Float = p_hit.y.atan2(p_hit.x);
        if phi < 0.0 as Float {
            phi += 2.0 as Float * PI;
        }
        // test cylinder intersection against clipping parameters
        if p_hit.z < self.z_min || p_hit.z > self.z_max || phi > self.phi_max {
            if t_shape_hit == t1 {
                return false;
            }
            t_shape_hit = t1;
            if t1.upper_bound() > ray.t_max {
                return false;
            }
            // compute cylinder hit point and $\phi$
            p_hit = ray.position(t_shape_hit.v);

            // refine cylinder intersection point
            let hit_rad: Float = (p_hit.x * p_hit.x + p_hit.y * p_hit.y).sqrt();
            p_hit.x *= self.radius / hit_rad;
            p_hit.y *= self.radius / hit_rad;
            phi = p_hit.y.atan2(p_hit.x);
            if phi < 0.0 as Float {
                phi += 2.0 as Float * PI;
            }
            if p_hit.z < self.z_min || p_hit.z > self.z_max || phi > self.phi_max {
                return false;
            }
        }
        true
    }
    fn get_reverse_orientation(&self) -> bool {
        self.reverse_orientation
    }
    fn get_transform_swaps_handedness(&self) -> bool {
        self.transform_swaps_handedness
    }
    fn area(&self) -> Float {
        (self.z_max - self.z_min) * self.radius * self.phi_max
    }
    fn sample(&self, u: Point2f, pdf: &mut Float) -> InteractionCommon {
        let z: Float = lerp(u[0], self.z_min, self.z_max);
        let phi: Float = u[1] * self.phi_max;
        let mut p_obj: Point3f = Point3f {
            x: self.radius * phi.cos(),
            y: self.radius * phi.sin(),
            z: z,
        };
        let mut it: InteractionCommon = InteractionCommon::default();
        it.n = nrm_normalize(self.object_to_world.transform_normal(Normal3f { x: p_obj.x,
                                                                              y: p_obj.y,
                                                                              z: 0.0,
        }));
        if self.reverse_orientation {
            it.n *= -1.0 as Float;
        }
        // reproject _p_obj_ to cylinder surface and compute _p_obj_error_
        let hit_rad: Float = (p_obj.x * p_obj.x + p_obj.y * p_obj.y).sqrt();
        p_obj.x *= self.radius / hit_rad;
        p_obj.y *= self.radius / hit_rad;
        let p_obj_error: Vector3f = Vector3f{
            x: p_obj.x,
            y: p_obj.y,
            z: 0.0,
        }.abs() * gamma(3_i32);
        it.p = self.object_to_world.transform_point_with_abs_error(p_obj, &p_obj_error, &mut it.p_error);
        *pdf = 1.0 as Float / self.area();
        it
    }
    fn sample_with_ref_point(&self, iref: &InteractionCommon, u: Point2f, pdf: &mut Float) -> InteractionCommon {
        let intr: InteractionCommon = self.sample(u, pdf);
        let mut wi: Vector3f = intr.p - iref.p;
        if wi.length_squared() == 0.0 as Float {
            *pdf = 0.0 as Float;
        } else {
            wi = vec3_normalize(wi);
            // convert from area measure, as returned by the Sample()
            // call above, to solid angle measure.
            *pdf *= pnt3_distance_squared(iref.p, intr.p) / nrm_abs_dot_vec3(intr.n, -wi);
            if (*pdf).is_infinite() {
                *pdf = 0.0 as Float;
            }
        }
        intr
    }
    fn pdf(&self, iref: &Interaction, wi: Vector3f) -> Float {
        // intersect sample ray with area light geometry
        let ray: Ray = iref.spawn_ray(wi);
        // ignore any alpha textures used for trimming the shape when
        // performing this intersection. Hack for the "San Miguel"
        // scene, where this is used to make an invisible area light.
        if let Some((isect_light, _t_hit)) = self.intersect(&ray) {
            // convert light sample weight to solid angle measure
            let mut pdf: Float = pnt3_distance_squared(iref.get_p(), isect_light.p) /
                (nrm_abs_dot_vec3(isect_light.n, -wi) * self.area());
            if pdf.is_infinite() {
                pdf = 0.0 as Float;
            }
            pdf
        } else {
            0.0 as Float
        }
    }
}

// see disk.h

#[derive(Clone)]
pub struct Disk {
    pub height: Float,
    pub radius: Float,
    pub inner_radius: Float,
    pub phi_max: Float,
    // inherited from class Shape (see shape.h)
    object_to_world: Transform,
    world_to_object: Transform,
    reverse_orientation: bool,
    transform_swaps_handedness: bool,
    pub material: Option<Arc<Material + Send + Sync>>,
}

impl Default for Disk {
    fn default() -> Self {
        Disk {
            // Shape
            object_to_world: Transform::default(),
            world_to_object: Transform::default(),
            reverse_orientation: false,
            transform_swaps_handedness: false,
            // Disk
            height: 0.0,
            radius: 1.0,
            inner_radius: 0.0,
            phi_max: radians(360.0),
            material: None,
        }
    }
}

impl Disk {
    pub fn new(object_to_world: Transform,
               world_to_object: Transform,
               reverse_orientation: bool,
               transform_swaps_handedness: bool,
               height: Float,
               radius: Float,
               inner_radius: Float,
               phi_max: Float)
               -> Self {
        Disk {
            // Shape
            object_to_world: object_to_world,
            world_to_object: world_to_object,
            reverse_orientation: reverse_orientation,
            transform_swaps_handedness: transform_swaps_handedness,
            // Disk
            height: height,
            radius: radius,
            inner_radius: inner_radius,
            phi_max: radians(clamp_t(phi_max, 0.0, 360.0)),
            material: None,
        }
    }
}

impl Shape for Disk {
    fn object_bound(&self) -> Bounds3f {
        Bounds3f {
            p_min: Point3f {
                x: -self.radius,
                y: -self.radius,
                z: self.height,
            },
            p_max: Point3f {
                x: self.radius,
                y: self.radius,
                z: self.height,
            },
        }
    }
    fn world_bound(&self) -> Bounds3f {
        // in C++: Bounds3f Shape::WorldBound() const { return (*ObjectToWorld)(ObjectBound()); }
        self.object_to_world.transform_bounds(self.object_bound())
    }
    fn intersect(&self, r: &Ray) -> Option<(SurfaceInteraction, Float)> {
        // TODO: ProfilePhase p(Prof::ShapeIntersect);
        // transform _Ray_ to object space
        let mut o_err: Vector3f = Vector3f::default();
        let mut d_err: Vector3f = Vector3f::default();
        let ray: Ray = self.world_to_object
            .transform_ray_with_error(r, &mut o_err, &mut d_err);

        // compute plane intersection for disk

        // reject disk intersections for rays parallel to the disk's plane
        if ray.d.z == 0.0 {
            return None;
        }
        let t_shape_hit: Float = (self.height - ray.o.z) / ray.d.z;
        if t_shape_hit <= 0.0 || t_shape_hit >= ray.t_max {
            return None;
        }
        // see if hit point is inside disk radii and $\phimax$
        let mut p_hit: Point3f = ray.position(t_shape_hit);
        let dist2: Float = p_hit.x * p_hit.x + p_hit.y * p_hit.y;
        if dist2 > self.radius * self.radius || dist2 < self.inner_radius * self.inner_radius {
            return None;
        }
        // test disk $\phi$ value against $\phimax$
        let mut phi: Float = p_hit.y.atan2(p_hit.x);
        if phi < 0.0 {
            phi += 2.0_f32 * PI;
        }
        if phi > self.phi_max {
            return None;
        }
        // find parametric representation of disk hit
        let u: Float = phi / self.phi_max;
        let r_hit: Float = dist2.sqrt();
        let one_minus_v: Float = (r_hit - self.inner_radius) / (self.radius - self.inner_radius);
        let v: Float = 1.0 - one_minus_v;
        let dpdu: Vector3f = Vector3f {
            x: -self.phi_max * p_hit.y,
            y: self.phi_max * p_hit.x,
            z: 0.0,
        };
        let dpdv: Vector3f = Vector3f {
            x: p_hit.x,
            y: p_hit.y,
            z: 0.0,
        } * (self.inner_radius - self.radius) / r_hit;
        let dndu: Normal3f = Normal3f::default();
        let dndv: Normal3f = Normal3f::default();
        // refine disk intersection point
        p_hit.z = self.height;
        // compute error bounds for disk intersection
        let p_error: Vector3f = Vector3f::default();
        // initialize _SurfaceInteraction_ from parametric information
        let uv_hit: Point2f = Point2f { x: u, y: v };
        let wo: Vector3f = -ray.d;
        let si: SurfaceInteraction = SurfaceInteraction::new(p_hit,
                                                             p_error,
                                                             uv_hit,
                                                             wo,
                                                             dpdu,
                                                             dpdv,
                                                             dndu,
                                                             dndv,
                                                             ray.time,
                                                             None);
        let mut isect: SurfaceInteraction = self.object_to_world.transform_surface_interaction(&si);
        if let Some(_shape) = si.shape {
            isect.shape = si.shape;
        }
        if let Some(_primitive) = si.primitive {
            isect.primitive = si.primitive;
        }
        Some((isect, t_shape_hit))
    }
    fn intersect_p(&self, r: &Ray) -> bool {
        // TODO: ProfilePhase p(Prof::ShapeIntersectP);
        // transform _Ray_ to object space
        let mut o_err: Vector3f = Vector3f::default();
        let mut d_err: Vector3f = Vector3f::default();
        let ray: Ray = self.world_to_object
            .transform_ray_with_error(r, &mut o_err, &mut d_err);

        // compute plane intersection for disk

        // reject disk intersections for rays parallel to the disk's plane
        if ray.d.z == 0.0 {
            return false;
        }
        let t_shape_hit: Float = (self.height - ray.o.z) / ray.d.z;
        if t_shape_hit <= 0.0 || t_shape_hit >= ray.t_max {
            return false;
        }
        // see if hit point is inside disk radii and $\phimax$
        let p_hit: Point3f = ray.position(t_shape_hit);
        let dist2: Float = p_hit.x * p_hit.x + p_hit.y * p_hit.y;
        if dist2 > self.radius * self.radius || dist2 < self.inner_radius * self.inner_radius {
            return false;
        }
        // test disk $\phi$ value against $\phimax$
        let mut phi: Float = p_hit.y.atan2(p_hit.x);
        if phi < 0.0 {
            phi += 2.0_f32 * PI;
        }
        if phi > self.phi_max {
            return false;
        }
        true
    }
    fn get_reverse_orientation(&self) -> bool {
        self.reverse_orientation
    }
    fn get_transform_swaps_handedness(&self) -> bool {
        self.transform_swaps_handedness
    }
    fn area(&self) -> Float {
        self.phi_max * 0.5 as Float *
        (self.radius * self.radius - self.inner_radius * self.inner_radius)
    }
    fn sample(&self, u: Point2f, pdf: &mut Float) -> InteractionCommon {
        let pd: Point2f = concentric_sample_disk(u);
        let p_obj: Point3f = Point3f {
            x: pd.x * self.radius,
            y: pd.y * self.radius,
            z: self.height,
        };
        let mut it: InteractionCommon = InteractionCommon::default();
        it.n = nrm_normalize(self.object_to_world.transform_normal(Normal3f {
            x: 0.0 as Float,
            y: 0.0 as Float,
            z: 1.0 as Float,
        }));
        if self.reverse_orientation {
            it.n *= -1.0 as Float;
        }
        let pt_error: Vector3f = Vector3f::default();
        it.p = self.object_to_world.transform_point_with_abs_error(p_obj,
                                                                   &pt_error,
                                                                   &mut it.p_error);
        *pdf = 1.0 as Float / self.area();
        it
    }
    fn sample_with_ref_point(&self, iref: &InteractionCommon, u: Point2f, pdf: &mut Float) -> InteractionCommon {
        let intr: InteractionCommon = self.sample(u, pdf);
        let mut wi: Vector3f = intr.p - iref.p;
        if wi.length_squared() == 0.0 as Float {
            *pdf = 0.0 as Float;
        } else {
            wi = vec3_normalize(wi);
            // convert from area measure, as returned by the Sample()
            // call above, to solid angle measure.
            *pdf *= pnt3_distance_squared(iref.p, intr.p) / nrm_abs_dot_vec3(intr.n, -wi);
            if (*pdf).is_infinite() {
                *pdf = 0.0 as Float;
            }
        }
        intr
    }
    fn pdf(&self, iref: &Interaction, wi: Vector3f) -> Float {
        // intersect sample ray with area light geometry
        let ray: Ray = iref.spawn_ray(wi);
        // ignore any alpha textures used for trimming the shape when
        // performing this intersection. Hack for the "San Miguel"
        // scene, where this is used to make an invisible area light.
        if let Some((isect_light, _t_hit)) = self.intersect(&ray) {
            // convert light sample weight to solid angle measure
            let mut pdf: Float = pnt3_distance_squared(iref.get_p(), isect_light.p) /
                (nrm_abs_dot_vec3(isect_light.n, -wi) * self.area());
            if pdf.is_infinite() {
                pdf = 0.0 as Float;
            }
            pdf
        } else {
            0.0 as Float
        }
    }
}

// see sphere.h

#[derive(Clone)]
pub struct Sphere {
    pub radius: Float,
    pub z_min: Float,
    pub z_max: Float,
    pub theta_min: Float,
    pub theta_max: Float,
    pub phi_max: Float,
    // inherited from class Shape (see shape.h)
    object_to_world: Transform,
    world_to_object: Transform,
    reverse_orientation: bool,
    transform_swaps_handedness: bool,
    pub material: Option<Arc<Material + Send + Sync>>,
}

impl Default for Sphere {
    fn default() -> Self {
        Sphere {
            // Shape
            object_to_world: Transform::default(),
            world_to_object: Transform::default(),
            reverse_orientation: false,
            transform_swaps_handedness: false,
            // Sphere
            radius: 1.0,
            z_min: -1.0,
            z_max: 1.0,
            theta_min: (-1.0 as Float).acos(),
            theta_max: (1.0 as Float).acos(),
            phi_max: radians(360.0),
            material: None,
        }
    }
}

impl Sphere {
    pub fn new(object_to_world: Transform,
               world_to_object: Transform,
               reverse_orientation: bool,
               transform_swaps_handedness: bool,
               radius: Float,
               z_min: Float,
               z_max: Float,
               phi_max: Float)
               -> Self {
        Sphere {
            // Shape
            object_to_world: object_to_world,
            world_to_object: world_to_object,
            reverse_orientation: reverse_orientation,
            transform_swaps_handedness: transform_swaps_handedness,
            // Sphere
            radius: radius,
            z_min: clamp_t(z_min.min(z_max), -radius, radius),
            z_max: clamp_t(z_min.max(z_max), -radius, radius),
            theta_min: clamp_t(z_min.min(z_max) / radius, -1.0, 1.0).acos(),
            theta_max: clamp_t(z_min.max(z_max) / radius, -1.0, 1.0).acos(),
            phi_max: radians(clamp_t(phi_max, 0.0, 360.0)),
            material: None,
        }
    }
}

impl Shape for Sphere {
    fn object_bound(&self) -> Bounds3f {
        Bounds3f {
            p_min: Point3f {
                x: -self.radius,
                y: -self.radius,
                z: self.z_min,
            },
            p_max: Point3f {
                x: self.radius,
                y: self.radius,
                z: self.z_max,
            },
        }
    }
    fn world_bound(&self) -> Bounds3f {
        // in C++: Bounds3f Shape::WorldBound() const { return (*ObjectToWorld)(ObjectBound()); }
        self.object_to_world.transform_bounds(self.object_bound())
    }
    fn intersect(&self, r: &Ray) -> Option<(SurfaceInteraction, Float)> {
        // transform _Ray_ to object space
        let mut o_err: Vector3f = Vector3f::default();
        let mut d_err: Vector3f = Vector3f::default();
        let ray: Ray = self.world_to_object.transform_ray_with_error(r, &mut o_err, &mut d_err);

        // compute quadratic sphere coefficients

        // initialize _EFloat_ ray coordinate values
        let ox = EFloat::new(ray.o.x as f32, o_err.x as f32);
        let oy = EFloat::new(ray.o.y as f32, o_err.y as f32);
        let oz = EFloat::new(ray.o.z as f32, o_err.z as f32);
        let dx = EFloat::new(ray.d.x as f32, d_err.x as f32);
        let dy = EFloat::new(ray.d.y as f32, d_err.y as f32);
        let dz = EFloat::new(ray.d.z as f32, d_err.z as f32);
        let a: EFloat = dx * dx + dy * dy + dz * dz;
        let b: EFloat = (dx * ox + dy * oy + dz * oz) * 2.0f32;
        let c: EFloat = ox * ox + oy * oy + oz * oz -
                        EFloat::new(self.radius as f32, 0.0) * EFloat::new(self.radius as f32, 0.0);

        // solve quadratic equation for _t_ values
        let mut t0: EFloat = EFloat::default();
        let mut t1: EFloat = EFloat::default();
        if !quadratic_efloat(a, b, c, &mut t0, &mut t1) {
            return None;
        }
        // check quadric shape _t0_ and _t1_ for nearest intersection
        if t0.upper_bound() > ray.t_max as f32 || t1.lower_bound() <= 0.0f32 {
            return None;
        }
        let mut t_shape_hit: EFloat = t0;
        if t_shape_hit.lower_bound() <= 0.0f32 {
            t_shape_hit = t1;
            if t_shape_hit.upper_bound() > ray.t_max as f32 {
                return None;
            }
        }
        // compute sphere hit position and $\phi$
        let mut p_hit: Point3f = ray.position(t_shape_hit.v);
        // refine sphere intersection point
        p_hit *= self.radius / pnt3_distance(p_hit, Point3f::default());
        if p_hit.x == 0.0 && p_hit.y == 0.0 {
            p_hit.x = 1e-5_f32 * self.radius;
        }
        let mut phi: Float = p_hit.y.atan2(p_hit.x);
        if phi < 0.0 {
            phi += 2.0_f32 * PI;
        }
        // test sphere intersection against clipping parameters
        if (self.z_min > -self.radius && p_hit.z < self.z_min) ||
           (self.z_max < self.radius && p_hit.z > self.z_max) || phi > self.phi_max {
            if t_shape_hit == t1 {
                return None;
            }
            if t1.upper_bound() > ray.t_max as f32 {
                return None;
            }
            t_shape_hit = t1;
            // compute sphere hit position and $\phi$
            p_hit = ray.position(t_shape_hit.v);

            // refine sphere intersection point
            p_hit *= self.radius / pnt3_distance(p_hit, Point3f::default());
            if p_hit.x == 0.0 && p_hit.y == 0.0 {
                p_hit.x = 1e-5_f32 * self.radius;
            }
            phi = p_hit.y.atan2(p_hit.x);
            if phi < 0.0 {
                phi += 2.0_f32 * PI;
            }
            if (self.z_min > -self.radius && p_hit.z < self.z_min) ||
               (self.z_max < self.radius && p_hit.z > self.z_max) ||
               phi > self.phi_max {
                return None;
            }
        }
        // find parametric representation of sphere hit
        let u: Float = phi / self.phi_max;
        let theta: Float = clamp_t(p_hit.z / self.radius, -1.0, 1.0).acos();
        let v: Float = (theta - self.theta_min) / (self.theta_max - self.theta_min);
        // compute sphere $\dpdu$ and $\dpdv$
        let z_radius: Float = (p_hit.x * p_hit.x + p_hit.y * p_hit.y).sqrt();
        let inv_z_radius: Float = 1.0 / z_radius;
        let cos_phi: Float = p_hit.x * inv_z_radius;
        let sin_phi: Float = p_hit.y * inv_z_radius;
        let dpdu: Vector3f = Vector3f {
            x: -self.phi_max * p_hit.y,
            y: self.phi_max * p_hit.x,
            z: 0.0,
        };
        let dpdv: Vector3f = Vector3f {
            x: p_hit.z * cos_phi,
            y: p_hit.z * sin_phi,
            z: -self.radius * theta.sin(),
        } * (self.theta_max - self.theta_min);
        // compute sphere $\dndu$ and $\dndv$
        let d2_p_duu: Vector3f = Vector3f {
            x: p_hit.x,
            y: p_hit.y,
            z: 0.0,
        } * -self.phi_max * self.phi_max;
        let d2_p_duv: Vector3f = Vector3f {
            x: -sin_phi,
            y: cos_phi,
            z: 0.0,
        } * (self.theta_max - self.theta_min) * p_hit.z *
                                 self.phi_max;
        let d2_p_dvv: Vector3f = Vector3f {
            x: p_hit.x,
            y: p_hit.y,
            z: p_hit.z,
        } * -(self.theta_max - self.theta_min) *
                                 (self.theta_max - self.theta_min);
        // compute coefficients for fundamental forms
        let ec: Float = vec3_dot_vec3(dpdu, dpdu);
        let fc: Float = vec3_dot_vec3(dpdu, dpdv);
        let gc: Float = vec3_dot_vec3(dpdv, dpdv);
        let nc: Vector3f = vec3_normalize(vec3_cross_vec3(dpdu, dpdv));
        let el: Float = vec3_dot_vec3(nc, d2_p_duu);
        let fl: Float = vec3_dot_vec3(nc, d2_p_duv);
        let gl: Float = vec3_dot_vec3(nc, d2_p_dvv);
        // compute $\dndu$ and $\dndv$ from fundamental form coefficients
        let inv_egf2: Float = 1.0 / (ec * gc - fc * fc);
        let dndu = dpdu * (fl * fc - el * gc) * inv_egf2 + dpdv * (el * fc - fl * ec) * inv_egf2;
        let dndu = Normal3f {
            x: dndu.x,
            y: dndu.y,
            z: dndu.z,
        };
        let dndv = dpdu * (gl * fc - fl * gc) * inv_egf2 + dpdv * (fl * fc - gl * ec) * inv_egf2;
        let dndv = Normal3f {
            x: dndv.x,
            y: dndv.y,
            z: dndv.z,
        };
        // compute error bounds for sphere intersection
        let p_error: Vector3f = Vector3f {
                x: p_hit.x,
                y: p_hit.y,
                z: p_hit.z,
            }
            .abs() * gamma(5_i32);
        // initialize _SurfaceInteraction_ from parametric information
        let uv_hit: Point2f = Point2f { x: u, y: v };
        let wo: Vector3f = -ray.d;
        let si: SurfaceInteraction =
            SurfaceInteraction::new(p_hit, p_error, uv_hit, wo, dpdu, dpdv, dndu, dndv, ray.time, None);
        let mut isect: SurfaceInteraction = self.object_to_world.transform_surface_interaction(&si);
        if let Some(_shape) = si.shape {
            isect.shape = si.shape;
        }
        if let Some(_primitive) = si.primitive {
            isect.primitive = si.primitive;
        }
        Some((isect, t_shape_hit.v as Float))
    }
    fn intersect_p(&self, r: &Ray) -> bool {
        // transform _Ray_ to object space
        let mut o_err: Vector3f = Vector3f::default();
        let mut d_err: Vector3f = Vector3f::default();
        let ray: Ray = self.world_to_object.transform_ray_with_error(r, &mut o_err, &mut d_err);

        // compute quadratic sphere coefficients

        // initialize _EFloat_ ray coordinate values
        let ox = EFloat::new(ray.o.x as f32, o_err.x as f32);
        let oy = EFloat::new(ray.o.y as f32, o_err.y as f32);
        let oz = EFloat::new(ray.o.z as f32, o_err.z as f32);
        let dx = EFloat::new(ray.d.x as f32, d_err.x as f32);
        let dy = EFloat::new(ray.d.y as f32, d_err.y as f32);
        let dz = EFloat::new(ray.d.z as f32, d_err.z as f32);
        let a: EFloat = dx * dx + dy * dy + dz * dz;
        let b: EFloat = (dx * ox + dy * oy + dz * oz) * 2.0f32;
        let c: EFloat = ox * ox + oy * oy + oz * oz -
                        EFloat::new(self.radius as f32, 0.0) * EFloat::new(self.radius as f32, 0.0);

        // solve quadratic equation for _t_ values
        let mut t0: EFloat = EFloat::default();
        let mut t1: EFloat = EFloat::default();
        if !quadratic_efloat(a, b, c, &mut t0, &mut t1) {
            return false;
        }
        // check quadric shape _t0_ and _t1_ for nearest intersection
        if t0.upper_bound() > ray.t_max as f32 || t1.lower_bound() <= 0.0f32 {
            return false;
        }
        let mut t_shape_hit: EFloat = t0;
        if t_shape_hit.lower_bound() <= 0.0f32 {
            t_shape_hit = t1;
            if t_shape_hit.upper_bound() > ray.t_max as f32 {
                return false;
            }
        }
        // compute sphere hit position and $\phi$
        let mut p_hit: Point3f = ray.position(t_shape_hit.v);
        // refine sphere intersection point
        p_hit *= self.radius / pnt3_distance(p_hit, Point3f::default());
        if p_hit.x == 0.0 && p_hit.y == 0.0 {
            p_hit.x = 1e-5_f32 * self.radius;
        }
        let mut phi: Float = p_hit.y.atan2(p_hit.x);
        if phi < 0.0 {
            phi += 2.0_f32 * PI;
        }
        // test sphere intersection against clipping parameters
        if (self.z_min > -self.radius && p_hit.z < self.z_min) ||
           (self.z_max < self.radius && p_hit.z > self.z_max) || phi > self.phi_max {
            if t_shape_hit == t1 {
                return false;
            }
            if t1.upper_bound() > ray.t_max as f32 {
                return false;
            }
            t_shape_hit = t1;
            // compute sphere hit position and $\phi$
            p_hit = ray.position(t_shape_hit.v);

            // refine sphere intersection point
            p_hit *= self.radius / pnt3_distance(p_hit, Point3f::default());
            if p_hit.x == 0.0 && p_hit.y == 0.0 {
                p_hit.x = 1e-5_f32 * self.radius;
            }
            phi = p_hit.y.atan2(p_hit.x);
            if phi < 0.0 {
                phi += 2.0_f32 * PI;
            }
            if (self.z_min > -self.radius && p_hit.z < self.z_min) ||
               (self.z_max < self.radius && p_hit.z > self.z_max) ||
               phi > self.phi_max {
                return false;
            }
        }
        true
    }
    fn get_reverse_orientation(&self) -> bool {
        self.reverse_orientation
    }
    fn get_transform_swaps_handedness(&self) -> bool {
        self.transform_swaps_handedness
    }
    fn area(&self) -> Float {
        self.phi_max * self.radius * (self.z_max - self.z_min)
    }
    fn sample(&self, u: Point2f, pdf: &mut Float) -> InteractionCommon {
        let mut p_obj: Point3f = Point3f::default() + uniform_sample_sphere(u) * self.radius;
        let mut it: InteractionCommon = InteractionCommon::default();
        it.n = nrm_normalize(self.object_to_world.transform_normal(Normal3f { x: p_obj.x,
                                                                              y: p_obj.y,
                                                                              z: p_obj.z,
        }));
        if self.reverse_orientation {
            it.n *= -1.0 as Float;
        }
        // reproject _p_obj_ to sphere surface and compute _p_obj_error_
        p_obj *= self.radius / pnt3_distance(p_obj, Point3f::default());
        let p_obj_error: Vector3f = Vector3f::from(p_obj).abs() * gamma(5_i32);
        it.p = self.object_to_world.transform_point_with_abs_error(p_obj, &p_obj_error, &mut it.p_error);
        *pdf = 1.0 as Float / self.area();
        it
    }
    fn sample_with_ref_point(&self, iref: &InteractionCommon, u: Point2f, pdf: &mut Float) -> InteractionCommon {
        let p_center: Point3f = self.object_to_world.transform_point(Point3f::default());
        // sample uniformly on sphere if $\pt{}$ is inside it
        let p_origin: Point3f = pnt3_offset_ray_origin(iref.p, iref.p_error, iref.n, p_center - iref.p);
        if pnt3_distance_squared(p_origin, p_center) <= self.radius * self.radius {
            let intr: InteractionCommon = self.sample(u, pdf);
            let mut wi: Vector3f = intr.p - iref.p;
            if wi.length_squared() == 0.0 as Float {
                *pdf = 0.0 as Float;
            } else {
                // convert from area measure returned by Sample() call
                // above to solid angle measure.
                wi = vec3_normalize(wi);
                *pdf *= pnt3_distance_squared(iref.p, intr.p) / nrm_abs_dot_vec3(intr.n, -wi);
            }
            if (*pdf).is_infinite() {
                *pdf = 0.0 as Float;
            }
            return intr;
        }

        // compute coordinate system for sphere sampling
        let wc: Vector3f= vec3_normalize(p_center - iref.p);
        let mut wc_x: Vector3f = Vector3f::default();
        let mut wc_y: Vector3f = Vector3f::default();
        vec3_coordinate_system(&wc, &mut wc_x, &mut wc_y);
        // sample sphere uniformly inside subtended cone

        // compute $\theta$ and $\phi$ values for sample in cone
        let sin_theta_max2: Float = self.radius * self.radius / pnt3_distance_squared(iref.p, p_center);
        let cos_theta_max: Float = (0.0 as Float).max(1.0 as Float - sin_theta_max2).sqrt();
        let cos_theta: Float = (1.0 as Float - u[0]) + u[0] * cos_theta_max;
        let sin_theta: Float = (0.0 as Float).max(1.0 as Float - cos_theta * cos_theta).sqrt();
        let phi: Float = u[1] * 2.0 as Float * PI;
        // compute angle $\alpha$ from center of sphere to sampled point on surface
        let dc: Float = pnt3_distance(iref.p, p_center);
        let ds: Float = dc * cos_theta -
            (0.0 as Float).max(self.radius * self.radius - dc * dc * sin_theta * sin_theta).sqrt();
        let cos_alpha: Float = (dc * dc + self.radius * self.radius - ds * ds) / (2.0 as Float * dc * self.radius);
        let sin_alpha: Float = (0.0 as Float).max(1.0 as Float - cos_alpha * cos_alpha).sqrt();
        // compute surface normal and sampled point on sphere
        let n_world: Vector3f = spherical_direction_vec3(sin_alpha, cos_alpha, phi,
                                                         &(-wc_x), &(-wc_y), &(-wc));
        let p_world: Point3f = p_center + Point3f {
            x: n_world.x,
            y: n_world.y,
            z: n_world.z,
        } * self.radius;
        // return _Interaction_ for sampled point on sphere
        let mut it: InteractionCommon = InteractionCommon::default();
        it.p = p_world;
        it.p_error = Vector3f::from(p_world).abs() * gamma(5_i32);
        it.n = Normal3f::from(n_world);
        if self.reverse_orientation {
            it.n *= -1.0 as Float;
        }
        // uniform cone PDF.
        *pdf = 1.0 as Float / (2.0 as Float * PI * (1.0 as Float - cos_theta_max));
        it
    }
    fn pdf(&self, iref: &Interaction, wi: Vector3f) -> Float {
        let p_center: Point3f = self.object_to_world.transform_point(Point3f::default());
        // return uniform PDF if point is inside sphere
        let p_origin: Point3f = pnt3_offset_ray_origin(iref.get_p(),
                                                       iref.get_p_error(),
                                                       iref.get_n(),
                                                       p_center - iref.get_p());
        if pnt3_distance_squared(p_origin, p_center) <= self.radius * self.radius {
            // return Shape::Pdf(ref, wi);

            // intersect sample ray with area light geometry
            let ray: Ray = iref.spawn_ray(wi);
            // ignore any alpha textures used for trimming the shape when
            // performing this intersection. Hack for the "San Miguel"
            // scene, where this is used to make an invisible area light.
            if let Some((isect_light, _t_hit)) = self.intersect(&ray) {
                // convert light sample weight to solid angle measure
                let mut pdf: Float = pnt3_distance_squared(iref.get_p(), isect_light.p) /
                    (nrm_abs_dot_vec3(isect_light.n, -wi) * self.area());
                if pdf.is_infinite() {
                    pdf = 0.0 as Float;
                }
                return pdf;
            } else {
                return 0.0 as Float;
            }
        }
        // compute general sphere PDF
        let sin_theta_max2: Float = self.radius * self.radius / pnt3_distance_squared(iref.get_p(), p_center);
        let cos_theta_max: Float = (0.0 as Float).max(1.0 as Float - sin_theta_max2).sqrt();
        return uniform_cone_pdf(cos_theta_max);

    }
}

// see triangle.h

#[derive(Debug,Clone)]
pub struct TriangleMesh {
    /// the total number of triangles in the mesh
    pub n_triangles: usize,
    /// vector of vertex indices
    pub vertex_indices: Vec<usize>,
    /// the total number of vertices in the mesh
    pub n_vertices: usize,
    /// vector of *n_vertices* vertex positions
    pub p: Vec<Point3f>,
    /// an optional vector of normal vectors (can be empty)
    pub n: Vec<Normal3f>,
    /// an optional vector of tangent vectors (can be empty)
    pub s: Vec<Vector3f>,
    /// an optional vector of paramtric (u, v) values (texture coordinates)
    pub uv: Vec<Point2f>,
    // TODO: std::shared_ptr<Texture<Float>> alphaMask, shadowAlphaMask;
    // inherited from class Shape (see shape.h)
    pub object_to_world: Transform, // TODO: not pub?
    pub world_to_object: Transform, // TODO: not pub?
    reverse_orientation: bool,
    pub transform_swaps_handedness: bool, // TODO: not pub?
}

impl TriangleMesh {
    pub fn new(object_to_world: Transform,
               world_to_object: Transform,
               reverse_orientation: bool,
               transform_swaps_handedness: bool,
               n_triangles: usize,
               vertex_indices: Vec<usize>,
               n_vertices: usize,
               p: Vec<Point3f>,
               s: Vec<Vector3f>,
               n: Vec<Normal3f>,
               uv: Vec<Point2f>)
               -> Self {
        TriangleMesh {
            // Shape
            object_to_world: object_to_world,
            world_to_object: world_to_object,
            reverse_orientation: reverse_orientation,
            transform_swaps_handedness: transform_swaps_handedness,
            // TriangleMesh
            n_triangles: n_triangles,
            vertex_indices: vertex_indices,
            n_vertices: n_vertices,
            p: p,
            n: n,
            s: s,
            uv: uv,
        }
    }
}

#[derive(Clone)]
pub struct Triangle {
    mesh: Arc<TriangleMesh>,
    pub id: usize,
    // inherited from class Shape (see shape.h)
    object_to_world: Transform,
    world_to_object: Transform,
    reverse_orientation: bool,
    transform_swaps_handedness: bool,
    pub material: Option<Arc<Material + Send + Sync>>,
}

impl Triangle {
    pub fn new(object_to_world: Transform,
               world_to_object: Transform,
               reverse_orientation: bool,
               mesh: Arc<TriangleMesh>,
               tri_number: usize)
               -> Self {
        Triangle {
            mesh: mesh,
            id: tri_number,
            object_to_world: object_to_world,
            world_to_object: world_to_object,
            reverse_orientation: reverse_orientation,
            transform_swaps_handedness: false,
            material: None,
        }
    }
    pub fn get_uvs(&self) -> [Point2f; 3] {
        if self.mesh.uv.is_empty() {
            [Point2f { x: 0.0, y: 0.0 }, Point2f { x: 1.0, y: 0.0 }, Point2f { x: 1.0, y: 1.0 }]
        } else {
            [self.mesh.uv[self.mesh.vertex_indices[self.id * 3 + 0]],
             self.mesh.uv[self.mesh.vertex_indices[self.id * 3 + 1]],
             self.mesh.uv[self.mesh.vertex_indices[self.id * 3 + 2]]]
        }
    }
}

impl Shape for Triangle {
    fn object_bound(&self) -> Bounds3f {
        let p0: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 0]];
        let p1: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 1]];
        let p2: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 2]];
        bnd3_union_pnt3(Bounds3f::new(self.world_to_object.transform_point(p0),
                                      self.world_to_object.transform_point(p1)),
                        self.world_to_object.transform_point(p2))
    }
    fn world_bound(&self) -> Bounds3f {
        let p0: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 0]];
        let p1: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 1]];
        let p2: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 2]];
        bnd3_union_pnt3(Bounds3f::new(p0, p1), p2)
    }
    fn intersect(&self, ray: &Ray) -> Option<(SurfaceInteraction, Float)> {
        // get triangle vertices in _p0_, _p1_, and _p2_
        let p0: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 0]];
        let p1: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 1]];
        let p2: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 2]];
        // translate vertices based on ray origin
        let mut p0t: Point3f = p0 -
                               Vector3f {
            x: ray.o.x,
            y: ray.o.y,
            z: ray.o.z,
        };
        let mut p1t: Point3f = p1 -
                               Vector3f {
            x: ray.o.x,
            y: ray.o.y,
            z: ray.o.z,
        };
        let mut p2t: Point3f = p2 -
                               Vector3f {
            x: ray.o.x,
            y: ray.o.y,
            z: ray.o.z,
        };
        // permute components of triangle vertices and ray direction
        let kz: usize = vec3_max_dimension(ray.d.abs());
        let mut kx: usize = kz + 1;
        if kx == 3 {
            kx = 0;
        }
        let mut ky: usize = kx + 1;
        if ky == 3 {
            ky = 0;
        }
        let d: Vector3f = vec3_permute(ray.d, kx, ky, kz);
        p0t = pnt3_permute(p0t, kx, ky, kz);
        p1t = pnt3_permute(p1t, kx, ky, kz);
        p2t = pnt3_permute(p2t, kx, ky, kz);
        // apply shear transformation to translated vertex positions
        let sx: Float = -d.x / d.z;
        let sy: Float = -d.y / d.z;
        let sz: Float = 1.0 / d.z;
        p0t.x += sx * p0t.z;
        p0t.y += sy * p0t.z;
        p1t.x += sx * p1t.z;
        p1t.y += sy * p1t.z;
        p2t.x += sx * p2t.z;
        p2t.y += sy * p2t.z;
        // compute edge function coefficients _e0_, _e1_, and _e2_
        let mut e0: Float = p1t.x * p2t.y - p1t.y * p2t.x;
        let mut e1: Float = p2t.x * p0t.y - p2t.y * p0t.x;
        let mut e2: Float = p0t.x * p1t.y - p0t.y * p1t.x;
        // fall back to double precision test at triangle edges
        if mem::size_of::<Float>() == mem::size_of::<f32>() && (e0 == 0.0 || e1 == 0.0 || e2 == 0.0) {
            let p2txp1ty: f64 = p2t.x as f64 * p1t.y as f64;
            let p2typ1tx: f64 = p2t.y as f64 * p1t.x as f64;
            e0 = (p2typ1tx - p2txp1ty) as Float;
            let p0txp2ty = p0t.x as f64 * p2t.y as f64;
            let p0typ2tx = p0t.y as f64 * p2t.x as f64;
            e1 = (p0typ2tx - p0txp2ty) as Float;
            let p1txp0ty = p1t.x as f64 * p0t.y as f64;
            let p1typ0tx = p1t.y as f64 * p0t.x as f64;
            e2 = (p1typ0tx - p1txp0ty) as Float;
        }
        // perform triangle edge and determinant tests
        if (e0 < 0.0 || e1 < 0.0 || e2 < 0.0) && (e0 > 0.0 || e1 > 0.0 || e2 > 0.0) {
            return None;
        }
        let det: Float = e0 + e1 + e2;
        if det == 0.0 {
            return None;
        }
        // compute scaled hit distance to triangle and test against ray $t$ range
        p0t.z *= sz;
        p1t.z *= sz;
        p2t.z *= sz;
        let t_scaled: Float = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
        if det < 0.0 && (t_scaled >= 0.0 || t_scaled < ray.t_max * det) {
            return None;
        } else if det > 0.0 && (t_scaled <= 0.0 || t_scaled > ray.t_max * det) {
            return None;
        }
        // compute barycentric coordinates and $t$ value for triangle intersection
        let inv_det: Float = 1.0 / det;
        let b0: Float = e0 * inv_det;
        let b1: Float = e1 * inv_det;
        let b2: Float = e2 * inv_det;
        let t: Float = t_scaled * inv_det;

        // ensure that computed triangle $t$ is conservatively greater than zero

        // compute $\delta_z$ term for triangle $t$ error bounds
        let max_zt: Float = vec3_max_component(Vector3f {
                x: p0t.z,
                y: p1t.z,
                z: p2t.z,
            }
            .abs());
        let delta_z: Float = gamma(3_i32) * max_zt;
        // compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
        let max_xt: Float = vec3_max_component(Vector3f {
                x: p0t.x,
                y: p1t.x,
                z: p2t.x,
            }
            .abs());
        let max_yt: Float = vec3_max_component(Vector3f {
                x: p0t.y,
                y: p1t.y,
                z: p2t.y,
            }
            .abs());
        let delta_x: Float = gamma(5) * (max_xt + max_zt);
        let delta_y: Float = gamma(5) * (max_yt + max_zt);
        // compute $\delta_e$ term for triangle $t$ error bounds
        let delta_e: Float = 2.0 *
                             (gamma(2) * max_xt * max_yt + delta_y * max_xt + delta_x * max_yt);
        // compute $\delta_t$ term for triangle $t$ error bounds and check _t_
        let max_e: Float = vec3_max_component(Vector3f {
                x: e0,
                y: e1,
                z: e2,
            }
            .abs());
        let delta_t: Float =
            3.0 * (gamma(3) * max_e * max_zt + delta_e * max_zt + delta_z * max_e) * inv_det.abs();
        if t <= delta_t {
            return None;
        }
        // compute triangle partial derivatives
        let uv: [Point2f; 3] = self.get_uvs();
        // compute deltas for triangle partial derivatives
        let duv02: Vector2f = uv[0] - uv[2];
        let duv12: Vector2f = uv[1] - uv[2];
        let dp02: Vector3f = p0 - p2;
        let dp12: Vector3f = p1 - p2;
        let determinant: Float = duv02.x * duv12.y - duv02.y * duv12.x;
        let degenerate_uv: bool = determinant.abs() < 1e-8 as Float;
        // Vector3f dpdu, dpdv;
        let mut dpdu: Vector3f = Vector3f::default();
        let mut dpdv: Vector3f = Vector3f::default();
        if !degenerate_uv {
            let invdet: Float = 1.0 / determinant;
            dpdu = (dp02 * duv12.y - dp12 * duv02.y) * invdet;
            dpdv = (dp02 * -duv12.x + dp12 * duv02.x) * invdet;
        }
        if degenerate_uv || vec3_cross_vec3(dpdu, dpdv).length_squared() == 0.0 {
            // handle zero determinant for triangle partial derivative matrix
            vec3_coordinate_system(&vec3_normalize(vec3_cross_vec3(p2 - p0, p1 - p0)),
                                   &mut dpdu,
                                   &mut dpdv);
        }
        // compute error bounds for triangle intersection
        let x_abs_sum: Float = (b0 * p0.x).abs() + (b1 * p1.x).abs() + (b2 * p2.x).abs();
        let y_abs_sum: Float = (b0 * p0.y).abs() + (b1 * p1.y).abs() + (b2 * p2.y).abs();
        let z_abs_sum: Float = (b0 * p0.z).abs() + (b1 * p1.z).abs() + (b2 * p2.z).abs();
        let p_error: Vector3f = Vector3f {
            x: x_abs_sum,
            y: y_abs_sum,
            z: z_abs_sum,
        } * gamma(7);
        // interpolate $(u,v)$ parametric coordinates and hit point
        let p_hit: Point3f = p0 * b0 + p1 * b1 + p2 * b2;
        let uv_hit: Point2f = uv[0] * b0 + uv[1] * b1 + uv[2] * b2;
        // TODO: test intersection against alpha texture, if present
        // if (testAlphaTexture && mesh->alphaMask) {
        //     SurfaceInteraction isectLocal(p_hit, Vector3f(0, 0, 0), uv_hit, -ray.d,
        //                                   dpdu, dpdv, Normal3f(0, 0, 0),
        //                                   Normal3f(0, 0, 0), ray.time, this);
        //     if (mesh->alphaMask->Evaluate(isectLocal) == 0) return false;
        // }
        // fill in _SurfaceInteraction_ from triangle hit
        let dndu: Normal3f = Normal3f::default();
        let dndv: Normal3f = Normal3f::default();
        let wo: Vector3f = -ray.d;
        let mut si: SurfaceInteraction =
            SurfaceInteraction::new(p_hit, p_error, uv_hit, wo, dpdu, dpdv, dndu, dndv, ray.time, Some(self));
        // override surface normal in _isect_ for triangle
        let surface_normal: Normal3f = Normal3f::from(vec3_normalize(vec3_cross_vec3(dp02, dp12)));
        si.n = surface_normal;
        si.shading.n = surface_normal;
        if !self.mesh.n.is_empty() || !self.mesh.s.is_empty() {
            // initialize _Triangle_ shading geometry

            // compute shading normal _ns_ for triangle
            let mut ns: Normal3f;
            if !self.mesh.n.is_empty() {
                let n0 = self.mesh.n[self.mesh.vertex_indices[self.id * 3 + 0]];
                let n1 = self.mesh.n[self.mesh.vertex_indices[self.id * 3 + 1]];
                let n2 = self.mesh.n[self.mesh.vertex_indices[self.id * 3 + 2]];
                ns = Normal3::from(n0) * b0 + Normal3::from(n1) * b1 + Normal3::from(n2) * b2;
                if ns.length_squared() > 0.0 {
                    ns = nrm_normalize(ns);
                } else {
                    ns = si.n;
                }
            } else {
                ns = si.n;
            }
            // compute shading tangent _ss_ for triangle
            let mut ss: Vector3f;
            if !self.mesh.s.is_empty() {
                let s0 = self.mesh.s[self.mesh.vertex_indices[self.id * 3 + 0]];
                let s1 = self.mesh.s[self.mesh.vertex_indices[self.id * 3 + 1]];
                let s2 = self.mesh.s[self.mesh.vertex_indices[self.id * 3 + 2]];
                ss = s0 * b0 + s1 * b1 + s2 * b2;
                if ss.length_squared() > 0.0 {
                    ss = vec3_normalize(ss);
                } else {
                    ss = vec3_normalize(si.dpdu);
                }
            } else {
                ss = vec3_normalize(si.dpdu);
            }
            // compute shading bitangent _ts_ for triangle and adjust _ss_
            let mut ts: Vector3f = vec3_cross_nrm(ss, ns);
            if ts.length_squared() > 0.0 {
                ts = vec3_normalize(ts);
                ss = vec3_cross_nrm(ts, ns);
            } else {
                vec3_coordinate_system(&Vector3f::from(ns), &mut ss, &mut ts);
            }
            // compute $\dndu$ and $\dndv$ for triangle shading geometry
            let dndu: Normal3f;
            let dndv: Normal3f;
            if !self.mesh.n.is_empty() {
                // compute deltas for triangle partial derivatives of normal
                let duv02: Vector2f = uv[0] - uv[2];
                let duv12: Vector2f = uv[1] - uv[2];
                let dn1: Normal3f =
                    Normal3::from(self.mesh.n[self.mesh.vertex_indices[self.id * 3 + 0]]) -
                    Normal3::from(self.mesh.n[self.mesh.vertex_indices[self.id * 3 + 2]]);
                let dn2: Normal3f =
                    Normal3::from(self.mesh.n[self.mesh.vertex_indices[self.id * 3 + 1]]) -
                    Normal3::from(self.mesh.n[self.mesh.vertex_indices[self.id * 3 + 2]]);
                let determinant: Float = duv02.x * duv12.y - duv02.y * duv12.x;
                let degenerate_uv: bool = determinant.abs() < 1e-8;
                if degenerate_uv {
                    dndu = Normal3f::default();
                    dndv = Normal3f::default();
                } else {
                    let inv_det: Float = 1.0 / determinant;
                    dndu = (dn1 * duv12.y - dn2 * duv02.y) * inv_det;
                    dndv = (dn1 * -duv12.x + dn2 * duv02.x) * inv_det;
                }
            } else {
                dndu = Normal3f::default();
                dndv = Normal3f::default();
            }
            si.set_shading_geometry(ss, ts, dndu, dndv, true);
        }
        // ensure correct orientation of the geometric normal
        if !self.mesh.n.is_empty() {
            si.n = nrm_faceforward_nrm(si.n, si.shading.n);
        } else if self.reverse_orientation ^ self.transform_swaps_handedness {
            si.n = -si.n;
            si.shading.n = -si.n;
        }
        Some((si, t as Float))
    }
    fn intersect_p(&self, ray: &Ray) -> bool {
        // TODO: ProfilePhase p(Prof::TriIntersectP);
        // TODO: ++nTests;
        // get triangle vertices in _p0_, _p1_, and _p2_
        let p0: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 0]];
        let p1: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 1]];
        let p2: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 2]];
        // translate vertices based on ray origin
        let mut p0t: Point3f = p0 -
                               Vector3f {
            x: ray.o.x,
            y: ray.o.y,
            z: ray.o.z,
        };
        let mut p1t: Point3f = p1 -
                               Vector3f {
            x: ray.o.x,
            y: ray.o.y,
            z: ray.o.z,
        };
        let mut p2t: Point3f = p2 -
                               Vector3f {
            x: ray.o.x,
            y: ray.o.y,
            z: ray.o.z,
        };
        // permute components of triangle vertices and ray direction
        let kz: usize = vec3_max_dimension(ray.d.abs());
        let mut kx: usize = kz + 1;
        if kx == 3 {
            kx = 0;
        }
        let mut ky: usize = kx + 1;
        if ky == 3 {
            ky = 0;
        }
        let d: Vector3f = vec3_permute(ray.d, kx, ky, kz);
        p0t = pnt3_permute(p0t, kx, ky, kz);
        p1t = pnt3_permute(p1t, kx, ky, kz);
        p2t = pnt3_permute(p2t, kx, ky, kz);
        // apply shear transformation to translated vertex positions
        let sx: Float = -d.x / d.z;
        let sy: Float = -d.y / d.z;
        let sz: Float = 1.0 / d.z;
        p0t.x += sx * p0t.z;
        p0t.y += sy * p0t.z;
        p1t.x += sx * p1t.z;
        p1t.y += sy * p1t.z;
        p2t.x += sx * p2t.z;
        p2t.y += sy * p2t.z;
        // compute edge function coefficients _e0_, _e1_, and _e2_
        let mut e0: Float = p1t.x * p2t.y - p1t.y * p2t.x;
        let mut e1: Float = p2t.x * p0t.y - p2t.y * p0t.x;
        let mut e2: Float = p0t.x * p1t.y - p0t.y * p1t.x;
        // fall back to double precision test at triangle edges
        if mem::size_of::<Float>() == mem::size_of::<f32>() && (e0 == 0.0 || e1 == 0.0 || e2 == 0.0) {
            let p2txp1ty: f64 = p2t.x as f64 * p1t.y as f64;
            let p2typ1tx: f64 = p2t.y as f64 * p1t.x as f64;
            e0 = (p2typ1tx - p2txp1ty) as Float;
            let p0txp2ty = p0t.x as f64 * p2t.y as f64;
            let p0typ2tx = p0t.y as f64 * p2t.x as f64;
            e1 = (p0typ2tx - p0txp2ty) as Float;
            let p1txp0ty = p1t.x as f64 * p0t.y as f64;
            let p1typ0tx = p1t.y as f64 * p0t.x as f64;
            e2 = (p1typ0tx - p1txp0ty) as Float;
        }
        // perform triangle edge and determinant tests
        if (e0 < 0.0 || e1 < 0.0 || e2 < 0.0) && (e0 > 0.0 || e1 > 0.0 || e2 > 0.0) {
            return false;
        }
        let det: Float = e0 + e1 + e2;
        if det == 0.0 {
            return false;
        }
        // compute scaled hit distance to triangle and test against ray $t$ range
        p0t.z *= sz;
        p1t.z *= sz;
        p2t.z *= sz;
        let t_scaled: Float = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
        if det < 0.0 && (t_scaled >= 0.0 || t_scaled < ray.t_max * det) {
            return false;
        } else if det > 0.0 && (t_scaled <= 0.0 || t_scaled > ray.t_max * det) {
            return false;
        }
        // compute barycentric coordinates and $t$ value for triangle intersection
        let inv_det: Float = 1.0 / det;
        // let b0: Float = e0 * inv_det;
        // let b1: Float = e1 * inv_det;
        // let b2: Float = e2 * inv_det;
        let t: Float = t_scaled * inv_det;

        // ensure that computed triangle $t$ is conservatively greater than zero

        // compute $\delta_z$ term for triangle $t$ error bounds
        let max_zt: Float = vec3_max_component(Vector3f {
                x: p0t.z,
                y: p1t.z,
                z: p2t.z,
            }
            .abs());
        let delta_z: Float = gamma(3_i32) * max_zt;
        // compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
        let max_xt: Float = vec3_max_component(Vector3f {
                x: p0t.x,
                y: p1t.x,
                z: p2t.x,
            }
            .abs());
        let max_yt: Float = vec3_max_component(Vector3f {
                x: p0t.y,
                y: p1t.y,
                z: p2t.y,
            }
            .abs());
        let delta_x: Float = gamma(5) * (max_xt + max_zt);
        let delta_y: Float = gamma(5) * (max_yt + max_zt);
        // compute $\delta_e$ term for triangle $t$ error bounds
        let delta_e: Float = 2.0 *
                             (gamma(2) * max_xt * max_yt + delta_y * max_xt + delta_x * max_yt);
        // compute $\delta_t$ term for triangle $t$ error bounds and check _t_
        let max_e: Float = vec3_max_component(Vector3f {
                x: e0,
                y: e1,
                z: e2,
            }
            .abs());
        let delta_t: Float =
            3.0 * (gamma(3) * max_e * max_zt + delta_e * max_zt + delta_z * max_e) * inv_det.abs();
        if t <= delta_t {
            return false;
        }
        // TODO: if (testAlphaTexture && (mesh->alphaMask || mesh->shadowAlphaMask)) { ... }
        // TODO: ++nHits;
        true
    }
    fn get_reverse_orientation(&self) -> bool {
        self.reverse_orientation
    }
    fn get_transform_swaps_handedness(&self) -> bool {
        self.transform_swaps_handedness
    }
    fn area(&self) -> Float {
        // get triangle vertices in _p0_, _p1_, and _p2_
        let p0: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 0]];
        let p1: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 1]];
        let p2: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 2]];
        0.5 as Float * vec3_cross_vec3(p1 - p0, p2 - p0).length()
    }
    fn sample(&self, u: Point2f, pdf: &mut Float) -> InteractionCommon {
        let b: Point2f = uniform_sample_triangle(u);
        // get triangle vertices in _p0_, _p1_, and _p2_
        let p0: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 0]];
        let p1: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 1]];
        let p2: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 2]];
        let mut it: InteractionCommon = InteractionCommon::default();
        it.p = p0 * b[0] + p1 * b[1] + p2* (1.0 as Float - b[0] - b[1]);
        // compute surface normal for sampled point on triangle
        it.n = nrm_normalize(Normal3f::from(vec3_cross_vec3(p1 - p0, p2 - p0)));
        // ensure correct orientation of the geometric normal; follow
        // the same approach as was used in Triangle::Intersect().
        if !self.mesh.n.is_empty() {
            let ns: Normal3f =
                Normal3f::from(self.mesh.n[self.mesh.vertex_indices[self.id * 3 + 0]] * b[0] +
                               self.mesh.n[self.mesh.vertex_indices[self.id * 3 + 1]] * b[1] +
                               self.mesh.n[self.mesh.vertex_indices[self.id * 3 + 2]] * (1.0 as Float - b[0] - b[1]));
            it.n = nrm_faceforward_nrm(it.n, ns);
        } else if self.reverse_orientation ^ self.transform_swaps_handedness {
            it.n *= -1.0 as Float;
        }
        // compute error bounds for sampled point on triangle
        let p_abs_sum: Point3f =
            pnt3_abs(p0 * b[0]) + pnt3_abs(p1 * b[1]) + pnt3_abs(p2 * (1.0 as Float - b[0] - b[1]));
        it.p_error = Vector3f{ x: p_abs_sum.x,
                               y: p_abs_sum.y,
                               z: p_abs_sum.z, } * gamma(6);
        *pdf = 1.0 as Float / self.area();
        it
    }
    fn sample_with_ref_point(&self, iref: &InteractionCommon, u: Point2f, pdf: &mut Float) -> InteractionCommon {
        let intr: InteractionCommon = self.sample(u, pdf);
        let mut wi: Vector3f = intr.p - iref.p;
        if wi.length_squared() == 0.0 as Float {
            *pdf = 0.0 as Float;
        } else {
            wi = vec3_normalize(wi);
            // convert from area measure, as returned by the Sample()
            // call above, to solid angle measure.
            *pdf *= pnt3_distance_squared(iref.p, intr.p) / nrm_abs_dot_vec3(intr.n, -wi);
            if (*pdf).is_infinite() {
                *pdf = 0.0 as Float;
            }
        }
        intr
    }
    fn pdf(&self, iref: &Interaction, wi: Vector3f) -> Float {
        // intersect sample ray with area light geometry
        let ray: Ray = iref.spawn_ray(wi);
        // ignore any alpha textures used for trimming the shape when
        // performing this intersection. Hack for the "San Miguel"
        // scene, where this is used to make an invisible area light.
        if let Some((isect_light, _t_hit)) = self.intersect(&ray) {
            // convert light sample weight to solid angle measure
            let mut pdf: Float = pnt3_distance_squared(iref.get_p(), isect_light.p) /
                (nrm_abs_dot_vec3(isect_light.n, -wi) * self.area());
            if pdf.is_infinite() {
                pdf = 0.0 as Float;
            }
            pdf
        } else {
            0.0 as Float
        }
    }
}

// see spectrum.h

// #[derive(Debug,Clone)]
// pub enum SpectrumType {
//     Reflectance,
//     Illuminant,
// }

#[derive(Debug,Default,Copy,Clone)]
pub struct RGBSpectrum {
    pub c: [Float; 3],
}

impl RGBSpectrum {
    pub fn new(v: Float) -> Self {
        // let n_spectrum_samples = 3; // RGB
        RGBSpectrum { c: [v, v, v] }
        // TODO: DCHECK(!HasNaNs());
    }
    pub fn rgb(r: Float, g: Float, b: Float) -> RGBSpectrum {
        RGBSpectrum { c: [r, g, b] }
    }
    pub fn from_rgb(rgb: &[Float; 3]) -> Spectrum {
        let mut s: Spectrum = Spectrum::new(0.0 as Float);
        s.c[0] = rgb[0];
        s.c[1] = rgb[1];
        s.c[2] = rgb[2];
        // TODO: DCHECK(!s.HasNaNs());
        s
    }
    pub fn to_rgb(&self, rgb: &mut [Float; 3]) {
        rgb[0] = self.c[0];
        rgb[1] = self.c[1];
        rgb[2] = self.c[2];
    }
    pub fn to_xyz(&self, xyz: &mut [Float; 3]) {
        rgb_to_xyz(&self.c, xyz);
    }
    pub fn from_srgb(rgb: &[u8; 3]) -> RGBSpectrum {
        fn convert(v: u8) -> Float {
            let value = v as Float / 255.0;
            // see InverseGammaCorrect(Float value) in pbrt.h
            if value <= 0.04045 {
                value / 12.92
            } else {
                ((value + 0.055) * 1.0 / 1.055).powf(2.4)
            }
        }
        RGBSpectrum::rgb(convert(rgb[0]), convert(rgb[1]), convert(rgb[2]))
    }
    pub fn y(&self) -> Float {
        let y_weight: [Float; 3] = [0.212671, 0.715160, 0.072169];
        y_weight[0] * self.c[0] + y_weight[1] * self.c[1] + y_weight[2] * self.c[2]
    }
    // from CoefficientSpectrum
    pub fn is_black(&self) -> bool {
        for i in 0..3 {
            if self.c[i] != 0.0 as Float {
                return false;
            }
        }
        true
    }
    pub fn clamp(&self, low: Float, high: Float) -> RGBSpectrum {
        let mut ret: RGBSpectrum = RGBSpectrum::default();
        let n_spectrum_samples: usize = 3; // RGB
        for i in 0..n_spectrum_samples {
            ret.c[i] = clamp_t(self.c[i], low, high);
        }
        assert!(!ret.has_nans());
        ret
    }
    pub fn max_component_value(&self) -> Float {
        let mut m: Float = self.c[0];
        let n_spectrum_samples: usize = 3; // RGB
        for i in 1..n_spectrum_samples {
            m = m.max(self.c[i]);
        }
        m
    }
    pub fn has_nans(&self) -> bool {
        for i in 0..3 {
            if self.c[i].is_nan() {
                return true;
            }
        }
        false
    }
}

impl PartialEq for RGBSpectrum {
    fn eq(&self, rhs: &RGBSpectrum) -> bool {
        for i in 0..3 {
            if self.c[i] != rhs.c[i] {
                return false;
            }
        }
        true
    }
}

impl Add for RGBSpectrum {
    type Output = RGBSpectrum;
    fn add(self, rhs: RGBSpectrum) -> RGBSpectrum {
        RGBSpectrum { c: [self.c[0] + rhs.c[0], self.c[1] + rhs.c[1], self.c[2] + rhs.c[2]] }
    }
}

impl AddAssign for RGBSpectrum {
    fn add_assign(&mut self, rhs: RGBSpectrum) {
        // TODO: DCHECK(!s2.HasNaNs());
        self.c[0] += rhs.c[0];
        self.c[1] += rhs.c[1];
        self.c[2] += rhs.c[2];
    }
}

impl Mul for RGBSpectrum {
    type Output = RGBSpectrum;
    fn mul(self, rhs: RGBSpectrum) -> RGBSpectrum {
        RGBSpectrum { c: [self.c[0] * rhs.c[0], self.c[1] * rhs.c[1], self.c[2] * rhs.c[2]] }
    }
}

impl Mul<Float> for RGBSpectrum {
    type Output = RGBSpectrum;
    fn mul(self, rhs: Float) -> RGBSpectrum {
        RGBSpectrum { c: [self.c[0] * rhs, self.c[1] * rhs, self.c[2] * rhs] }
    }
}

impl MulAssign for RGBSpectrum {
    fn mul_assign(&mut self, rhs: RGBSpectrum) {
        // TODO: DCHECK(!HasNaNs());
        self.c[0] *= rhs.c[0];
        self.c[1] *= rhs.c[1];
        self.c[2] *= rhs.c[2];
    }
}

impl Sub for RGBSpectrum {
    type Output = RGBSpectrum;
    fn sub(self, rhs: RGBSpectrum) -> RGBSpectrum {
        RGBSpectrum { c: [self.c[0] - rhs.c[0], self.c[1] - rhs.c[1], self.c[2] - rhs.c[2]] }
    }
}

impl Div for RGBSpectrum {
    type Output = RGBSpectrum;
    fn div(self, rhs: RGBSpectrum) -> RGBSpectrum {
        RGBSpectrum { c: [self.c[0] / rhs.c[0], self.c[1] / rhs.c[1], self.c[2] / rhs.c[2]] }
    }
}

impl Div<Float> for RGBSpectrum {
    type Output = RGBSpectrum;
    fn div(self, rhs: Float) -> RGBSpectrum {
        assert_ne!(rhs, 0.0 as Float);
        assert!(!rhs.is_nan(), "rhs is NaN");
        let ret: RGBSpectrum =
            RGBSpectrum { c: [self.c[0] / rhs, self.c[1] / rhs, self.c[2] / rhs] };
        assert!(!ret.has_nans());
        ret
    }
}

impl Zero for RGBSpectrum {
    fn zero() -> RGBSpectrum {
        RGBSpectrum::new(0.0 as Float)
    }

    fn is_zero(&self) -> bool {
        self.is_black()
    }
}

impl Index<usize> for RGBSpectrum {
    type Output = Float;
    fn index(&self, index: usize) -> &Float {
        match index {
            0 => &self.c[0],
            1 => &self.c[1],
            2 => &self.c[2],
            _ => panic!("Check failed: i >= 0 && i <= 2"),
        }
    }
}

impl IndexMut<usize> for RGBSpectrum {
    fn index_mut(&mut self, index: usize) -> &mut Float {
        match index {
            0 => &mut self.c[0],
            1 => &mut self.c[1],
            2 => &mut self.c[2],
            _ => panic!("Check failed: i >= 0 && i <= 2"),
        }
    }
}

impl From<Float> for RGBSpectrum {
    fn from(f: Float) -> Self {
        RGBSpectrum::new(f)
    }
}

/// Calculate RGB coefficients from a XYZ representation.
pub fn xyz_to_rgb(xyz: &[Float; 3], rgb: &mut[Float; 3]) {
    rgb[0] = 3.240479 * xyz[0] - 1.537150 * xyz[1] - 0.498535 * xyz[2];
    rgb[1] = -0.969256 * xyz[0] + 1.875991 * xyz[1] + 0.041556 * xyz[2];
    rgb[2] = 0.055648 * xyz[0] - 0.204043 * xyz[1] + 1.057311 * xyz[2];
}

/// Calculate XYZ representation from RGB coefficients.
pub fn rgb_to_xyz(rgb: &[Float; 3], xyz: &mut[Float; 3]) {
    xyz[0] = 0.412453 * rgb[0] + 0.357580 * rgb[1] + 0.180423 * rgb[2];
    xyz[1] = 0.212671 * rgb[0] + 0.715160 * rgb[1] + 0.072169 * rgb[2];
    xyz[2] = 0.019334 * rgb[0] + 0.119193 * rgb[1] + 0.950227 * rgb[2];
}

/// Interpolate linearly between two provided values.
pub fn lerp_rgb(t: Float, s1: Spectrum, s2: Spectrum) -> Spectrum {
    s1 * (1.0 as Float - t) + s2 * t
}

// see bvh.h

#[derive(Debug,Clone)]
pub enum SplitMethod {
    SAH,
    HLBVH,
    Middle,
    EqualCounts,
}

#[derive(Debug,Default,Copy,Clone)]
pub struct BVHPrimitiveInfo {
    primitive_number: usize,
    bounds: Bounds3f,
    centroid: Point3f,
}

impl BVHPrimitiveInfo {
    pub fn new(primitive_number: usize, bounds: Bounds3f) -> Self {
        BVHPrimitiveInfo {
            primitive_number: primitive_number,
            bounds: bounds,
            centroid: bounds.p_min * 0.5 + bounds.p_max * 0.5,
        }
    }
}

#[derive(Debug)]
pub struct BVHBuildNode<'a> {
    pub bounds: Bounds3f,
    pub child1: Option<&'a mut BVHBuildNode<'a>>,
    pub child2: Option<&'a mut BVHBuildNode<'a>>,
    pub split_axis: u8,
    pub first_prim_offset: usize,
    pub n_primitives: usize,
}

impl<'a> Default for BVHBuildNode<'a> {
    fn default() -> Self {
        BVHBuildNode {
            bounds: Bounds3f::default(),
            child1: None,
            child2: None,
            split_axis: 0_u8,
            first_prim_offset: 0_usize,
            n_primitives: 0_usize,
        }
    }
}

impl<'a> BVHBuildNode<'a> {
    pub fn init_leaf(&mut self, first: usize, n: usize, b: &Bounds3f) {
        self.first_prim_offset = first;
        self.n_primitives = n;
        self.bounds = *b;
        self.child1 = None;
        self.child2 = None;
    }
    pub fn init_interior(&mut self, axis: u8, c0: &'a mut BVHBuildNode<'a>, c1: &'a mut BVHBuildNode<'a>) {
        self.n_primitives = 0;
        self.bounds = bnd3_union_bnd3(c0.bounds, c1.bounds);
        self.child1 = Some(c0);
        self.child2 = Some(c1);
        self.split_axis = axis;
    }
}

#[derive(Debug,Copy,Clone)]
struct BucketInfo {
    count: usize,
    bounds: Bounds3f,
}

impl Default for BucketInfo {
    fn default() -> Self {
        BucketInfo {
            count: 0_usize,
            bounds: Bounds3f::default(),
        }
    }
}

#[derive(Debug,Default,Copy,Clone)]
pub struct LinearBVHNode {
    bounds: Bounds3f,
    // in C++ a union { int primitivesOffset;     // leaf
    //                  int secondChildOffset; }; // interior
    offset: usize,
    n_primitives: usize,
    axis: u8, // TODO? pad
}

// BVHAccel -> Aggregate -> Primitive
pub struct BVHAccel {
    max_prims_in_node: usize,
    split_method: SplitMethod,
    pub primitives: Vec<Arc<Primitive + Sync + Send>>,
    pub nodes: Vec<LinearBVHNode>,
}

impl BVHAccel {
    pub fn new(p: Vec<Arc<Primitive + Sync + Send>>,
               max_prims_in_node: usize,
               split_method: SplitMethod)
               -> Self {
        let bvh = Arc::new(BVHAccel {
            max_prims_in_node: std::cmp::min(max_prims_in_node, 255),
            split_method: split_method.clone(),
            primitives: p,
            nodes: Vec::new(),
        });
        let num_prims = bvh.primitives.len();
        let mut primitive_info = vec![BVHPrimitiveInfo::default(); num_prims];
        for i in 0..num_prims {
            let world_bound = bvh.primitives[i].world_bound();
            primitive_info[i] = BVHPrimitiveInfo::new(i, world_bound);
        }
        // TODO: if (splitMethod == SplitMethod::HLBVH)
        let mut arena: Arena<BVHBuildNode> = Arena::with_capacity(1024 * 1024);
        let mut total_nodes: usize = 0;
        let mut ordered_prims: Vec<Arc<Primitive + Sync + Send>> = Vec::with_capacity(num_prims);
        println!("BVHAccel::recursive_build(...)");
        let start = PreciseTime::now();
        let root = BVHAccel::recursive_build(bvh.clone(), // instead of self
                                             &mut arena,
                                             &mut primitive_info,
                                             0,
                                             num_prims,
                                             &mut total_nodes,
                                             &mut ordered_prims);
        let end = PreciseTime::now();
        println!("{} seconds for building BVH ...", start.to(end));
        // flatten first
        let mut nodes = vec![LinearBVHNode::default(); total_nodes];
        let mut offset: usize = 0;
        println!("BVHAccel::flatten_bvh_tree(...)");
        let start = PreciseTime::now();
        BVHAccel::flatten_bvh_tree(root, &mut nodes, &mut offset);
        let end = PreciseTime::now();
        println!("{} seconds for flattening BVH ...", start.to(end));
        assert!(nodes.len() == total_nodes);
        // primitives.swap(orderedPrims);
        let bvh_ordered_prims = Arc::new(BVHAccel {
            max_prims_in_node: std::cmp::min(max_prims_in_node, 255),
            split_method: split_method.clone(),
            primitives: ordered_prims,
            nodes: nodes,
        });
        let unwrapped = Arc::try_unwrap(bvh_ordered_prims);
        unwrapped.ok().unwrap()
    }
    pub fn world_bound(&self) -> Bounds3f {
        if self.nodes.len() > 0 {
            self.nodes[0].bounds
        } else {
            Bounds3f::default()
        }
    }
    pub fn intersect(&self, ray: &mut Ray) -> Option<SurfaceInteraction> {
        if self.nodes.len() == 0 {
            return None;
        }
        // TODO: ProfilePhase p(Prof::AccelIntersect);
        let mut hit: bool = false;
        let inv_dir: Vector3f = Vector3f {
            x: 1.0 / ray.d.x,
            y: 1.0 / ray.d.y,
            z: 1.0 / ray.d.z,
        };
        let dir_is_neg: [u8; 3] =
            [(inv_dir.x < 0.0) as u8, (inv_dir.y < 0.0) as u8, (inv_dir.z < 0.0) as u8];
        // follow ray through BVH nodes to find primitive intersections
        let mut to_visit_offset: u32 = 0;
        let mut current_node_index: u32 = 0;
        let mut nodes_to_visit: [u32; 64] = [0_u32; 64];
        let mut si: SurfaceInteraction = SurfaceInteraction::default();
        loop {
            let node: LinearBVHNode = self.nodes[current_node_index as usize];
            // check ray against BVH node
            let intersects: bool = node.bounds.intersect_p(ray, &inv_dir, dir_is_neg);
            if intersects {
                if node.n_primitives > 0 {
                    // intersect ray with primitives in leaf BVH node
                    for i in 0..node.n_primitives {
                        // see primitive.h GeometricPrimitive::Intersect() ...
                        if let Some(isect) = self.primitives[node.offset + i].intersect(ray) {
                            // TODO: CHECK_GE(...)
                            si = isect;
                            hit = true;
                        }
                    }
                    if to_visit_offset == 0_u32 {
                        break;
                    }
                    to_visit_offset -= 1_u32;
                    current_node_index = nodes_to_visit[to_visit_offset as usize];
                } else {
                    // put far BVH node on _nodesToVisit_ stack,
                    // advance to near node
                    if dir_is_neg[node.axis as usize] == 1_u8 {
                        nodes_to_visit[to_visit_offset as usize] = current_node_index + 1_u32;
                        to_visit_offset += 1_u32;
                        current_node_index = node.offset as u32;
                    } else {
                        nodes_to_visit[to_visit_offset as usize] = node.offset as u32;
                        to_visit_offset += 1_u32;
                        current_node_index += 1_u32;
                    }
                }
            } else {
                if to_visit_offset == 0_u32 {
                    break;
                }
                to_visit_offset -= 1_u32;
                current_node_index = nodes_to_visit[to_visit_offset as usize];
            }
        }
        if hit {
            Some(si)
        } else {
            None
        }
    }
    pub fn intersect_p(&self, ray: &mut Ray) -> bool {
        if self.nodes.len() == 0 {
            return false;
        }
        // TODO: ProfilePhase p(Prof::AccelIntersectP);
        let inv_dir: Vector3f = Vector3f {
            x: 1.0 / ray.d.x,
            y: 1.0 / ray.d.y,
            z: 1.0 / ray.d.z,
        };
        let dir_is_neg: [u8; 3] =
            [(inv_dir.x < 0.0) as u8, (inv_dir.y < 0.0) as u8, (inv_dir.z < 0.0) as u8];
        let mut to_visit_offset: u32 = 0;
        let mut current_node_index: u32 = 0;
        let mut nodes_to_visit: [u32; 64] = [0_u32; 64];
        loop {
            let node: LinearBVHNode = self.nodes[current_node_index as usize];
            let intersects: bool = node.bounds.intersect_p(ray, &inv_dir, dir_is_neg);
            if intersects {
                // process BVH node _node_ for traversal
                if node.n_primitives > 0 {
                    for i in 0..node.n_primitives {
                        if self.primitives[node.offset + i].intersect_p(ray) {
                            return true;
                        }
                    }
                    if to_visit_offset == 0_u32 {
                        break;
                    }
                    to_visit_offset -= 1_u32;
                    current_node_index = nodes_to_visit[to_visit_offset as usize];
                } else {
                    if dir_is_neg[node.axis as usize] == 1_u8 {
                        nodes_to_visit[to_visit_offset as usize] = current_node_index + 1_u32;
                        to_visit_offset += 1_u32;
                        current_node_index = node.offset as u32;
                    } else {
                        nodes_to_visit[to_visit_offset as usize] = node.offset as u32;
                        to_visit_offset += 1_u32;
                        current_node_index += 1_u32;
                    }
                }
            } else {
                if to_visit_offset == 0_u32 {
                    break;
                }
                to_visit_offset -= 1_u32;
                current_node_index = nodes_to_visit[to_visit_offset as usize];
            }
        }
        false
    }
    pub fn recursive_build<'a>(bvh: Arc<BVHAccel>,
                               arena: &'a Arena<BVHBuildNode<'a>>,
                               primitive_info: &mut Vec<BVHPrimitiveInfo>,
                               start: usize,
                               end: usize,
                               total_nodes: &mut usize,
                               ordered_prims: &mut Vec<Arc<Primitive + Sync + Send>>)
                               -> &'a mut BVHBuildNode<'a> {
        assert_ne!(start, end);
        let node: &mut BVHBuildNode<'a> = arena.alloc(BVHBuildNode::default());
        *total_nodes += 1_usize;
        // compute bounds of all primitives in BVH node
        let mut bounds: Bounds3f = Bounds3f::default();
        for i in start..end {
            bounds = bnd3_union_bnd3(bounds, primitive_info[i].bounds);
        }
        let n_primitives: usize = end - start;
        if n_primitives == 1 {
            // create leaf _BVHBuildNode_
            let first_prim_offset: usize = ordered_prims.len();
            for i in start..end {
                let prim_num: usize = primitive_info[i].primitive_number;
                ordered_prims.push(bvh.primitives[prim_num].clone());
            }
            node.init_leaf(first_prim_offset, n_primitives, &bounds);
            return node;
        } else {
            // compute bound of primitive centroids, choose split dimension _dim_
            let mut centroid_bounds: Bounds3f = Bounds3f::default();
            for i in start..end {
                centroid_bounds = bnd3_union_pnt3(centroid_bounds, primitive_info[i].centroid);
            }
            let dim: u8 = centroid_bounds.maximum_extent();
            // partition primitives into two sets and build children
            let mut mid: usize = (start + end) / 2_usize;
            if centroid_bounds.p_max[dim] == centroid_bounds.p_min[dim] {
                // create leaf _BVHBuildNode_
                let first_prim_offset: usize = ordered_prims.len();
                for i in start..end {
                    let prim_num: usize = primitive_info[i].primitive_number;
                    ordered_prims.push(bvh.primitives[prim_num].clone());
                }
                node.init_leaf(first_prim_offset, n_primitives, &bounds);
                return node;
            } else {
                // partition primitives based on _splitMethod_
                match bvh.split_method {
                    SplitMethod::Middle => {
                        // TODO
                    }
                    SplitMethod::EqualCounts => {
                        // TODO
                    }
                    SplitMethod::SAH | SplitMethod::HLBVH => {
                        if n_primitives <= 2 {
                            mid = (start + end) / 2;
                            if start != end - 1 {
                                if primitive_info[end - 1].centroid[dim] <
                                   primitive_info[start].centroid[dim] {
                                    primitive_info.swap(start, end - 1);
                                }
                            }
                        } else {
                            // allocate _BucketInfo_ for SAH partition buckets
                            let n_buckets: usize = 12;
                            let mut buckets: [BucketInfo; 12] = [BucketInfo::default(); 12];
                            // initialize _BucketInfo_ for SAH partition buckets
                            for i in start..end {
                                let mut b: usize =
                                    (n_buckets as Float *
                                     centroid_bounds.offset(primitive_info[i].centroid)[dim])
                                    as usize;
                                if b == n_buckets {
                                    b = n_buckets - 1;
                                }
                                // assert!(b >= 0_usize, "b >= 0");
                                assert!(b < n_buckets, "b < {}", n_buckets);
                                buckets[b].count += 1;
                                buckets[b].bounds = bnd3_union_bnd3(buckets[b].bounds,
                                                                    primitive_info[i].bounds);
                            }
                            // compute costs for splitting after each bucket
                            let mut cost: [Float; 11] = [0.0; 11];
                            for i in 0..(n_buckets - 1) {
                                let mut b0: Bounds3f = Bounds3f::default();
                                let mut b1: Bounds3f = Bounds3f::default();
                                let mut count0: usize = 0;
                                let mut count1: usize = 0;
                                for j in 0..(i + 1) {
                                    b0 = bnd3_union_bnd3(b0, buckets[j].bounds);
                                    count0 += buckets[j].count;
                                }
                                for j in (i + 1)..n_buckets {
                                    b1 = bnd3_union_bnd3(b1, buckets[j].bounds);
                                    count1 += buckets[j].count;
                                }
                                cost[i] = 1.0 +
                                          (count0 as Float * b0.surface_area() +
                                           count1 as Float * b1.surface_area()) /
                                          bounds.surface_area();
                            }
                            // find bucket to split at that minimizes SAH metric
                            let mut min_cost: Float = cost[0];
                            let mut min_cost_split_bucket: usize = 0;
                            for i in 0..(n_buckets - 1) {
                                if cost[i] < min_cost {
                                    min_cost = cost[i];
                                    min_cost_split_bucket = i;
                                }
                            }
                            // either create leaf or split primitives
                            // at selected SAH bucket
                            let leaf_cost: Float = n_primitives as Float;
                            if n_primitives > bvh.max_prims_in_node || min_cost < leaf_cost {
                                let (mut left, mut right): (Vec<BVHPrimitiveInfo>, Vec<BVHPrimitiveInfo>) =
                                    primitive_info[start..end].into_iter()
                                    .partition(|&pi| {
                                        let mut b: usize =
                                            (n_buckets as Float *
                                             centroid_bounds.offset(pi.centroid)[dim]) as usize;
                                        if b == n_buckets {b = n_buckets - 1;}
                                        // assert!(b >= 0_usize, "b >= 0");
                                        assert!(b < n_buckets, "b < {}", n_buckets);
                                        b <= min_cost_split_bucket
                                    });
                                mid = start + left.len();
                                let combined_len = left.len() + right.len();
                                if combined_len == primitive_info.len() {
                                    primitive_info.clear();
                                    primitive_info.append(&mut left);
                                    primitive_info.append(&mut right);
                                } else {
                                    primitive_info.splice(start..mid, left.iter().cloned());
                                    primitive_info.splice(mid..end, right.iter().cloned());
                                }
                            } else {
                                // create leaf _BVHBuildNode_
                                let first_prim_offset: usize = ordered_prims.len();
                                for i in start..end {
                                    let prim_num: usize = primitive_info[i].primitive_number;
                                    ordered_prims.push(bvh.primitives[prim_num].clone());
                                }
                                node.init_leaf(first_prim_offset, n_primitives, &bounds);
                                return node;
                            }
                        }
                    }
                }
                // make sure we get result for c1 before c0
                let c1 = BVHAccel::recursive_build(bvh.clone(),
                                                   arena,
                                                   primitive_info,
                                                   mid,
                                                   end,
                                                   total_nodes,
                                                   ordered_prims);
                let c0 = BVHAccel::recursive_build(bvh.clone(),
                                                   arena,
                                                   primitive_info,
                                                   start,
                                                   mid,
                                                   total_nodes,
                                                   ordered_prims);
                node.init_interior(dim, c0, c1);
            }
        }
        return node;
    }
    fn flatten_bvh_tree<'a>(node: &mut BVHBuildNode<'a>,
                            nodes: &mut Vec<LinearBVHNode>,
                            offset: &mut usize)
                            -> usize {
        let my_offset: usize = *offset;
        *offset += 1;
        if node.n_primitives > 0 {
            // leaf
            let linear_node = LinearBVHNode {
                bounds: node.bounds,
                offset: node.first_prim_offset,
                n_primitives: node.n_primitives,
                axis: 0_u8,
            };
            nodes[my_offset] = linear_node;
        } else {
            // interior
            if let Some(ref mut child1) = node.child1 {
                BVHAccel::flatten_bvh_tree(child1, nodes, offset);
            }
            if let Some(ref mut child2) = node.child2 {
                let linear_node = LinearBVHNode {
                    bounds: node.bounds,
                    offset: BVHAccel::flatten_bvh_tree(child2, nodes, offset),
                    n_primitives: 0_usize,
                    axis: node.split_axis,
                };
                nodes[my_offset] = linear_node;
            }
        }
        my_offset
    }
}

// see sampler.h

pub trait Sampler {
    fn get_1d(&mut self) -> Float;
    fn get_2d(&mut self) -> Point2f;
}

// see zerotwosequence.h

#[derive(Debug,Clone)]
pub struct ZeroTwoSequenceSampler {
    pub samples_per_pixel: i64,
    pub n_sampled_dimensions: i64,
    // inherited from class PixelSampler (see sampler.h)
    pub samples_1d: Vec<Vec<Float>>,
    pub samples_2d: Vec<Vec<Point2f>>,
    pub current_1d_dimension: i32,
    pub current_2d_dimension: i32,
    pub rng: Rng,
    // inherited from class Sampler (see sampler.h)
    pub current_pixel: Point2i,
    pub current_pixel_sample_index: i64,
    pub samples_1d_array_sizes: Vec<i32>,
    pub samples_2d_array_sizes: Vec<i32>,
    pub samples_1d_array: Vec<Vec<Float>>,
    pub samples_2d_array: Vec<Vec<Point2f>>,
    pub array_1d_offset: usize,
    pub array_2d_offset: usize,
}

impl Default for ZeroTwoSequenceSampler {
    fn default() -> Self {
        let mut lds: ZeroTwoSequenceSampler = ZeroTwoSequenceSampler {
            samples_per_pixel: 1_i64,
            n_sampled_dimensions: 4_i64,
            samples_1d: Vec::new(),
            samples_2d: Vec::new(),
            current_1d_dimension: 0_i32,
            current_2d_dimension: 0_i32,
            rng: Rng::default(),
            current_pixel: Point2i::default(),
            current_pixel_sample_index: 0_i64,
            samples_1d_array_sizes: Vec::new(),
            samples_2d_array_sizes: Vec::new(),
            samples_1d_array: Vec::new(),
            samples_2d_array: Vec::new(),
            array_1d_offset: 0_usize,
            array_2d_offset: 0_usize,
        };
        for _i in 0..lds.n_sampled_dimensions {
            let additional_1d: Vec<Float> = vec![0.0; lds.samples_per_pixel as usize];
            let additional_2d: Vec<Point2f> =
                vec![Point2f::default(); lds.samples_per_pixel as usize];
            lds.samples_1d.push(additional_1d);
            lds.samples_2d.push(additional_2d);
        }
        lds
    }
}

impl ZeroTwoSequenceSampler {
    pub fn new(samples_per_pixel: i64, n_sampled_dimensions: i64) -> Self {
        let mut lds: ZeroTwoSequenceSampler = ZeroTwoSequenceSampler {
            samples_per_pixel: samples_per_pixel,
            n_sampled_dimensions: n_sampled_dimensions,
            samples_1d: Vec::new(),
            samples_2d: Vec::new(),
            current_1d_dimension: 0_i32,
            current_2d_dimension: 0_i32,
            rng: Rng::default(),
            current_pixel: Point2i::default(),
            current_pixel_sample_index: 0_i64,
            samples_1d_array_sizes: Vec::new(),
            samples_2d_array_sizes: Vec::new(),
            samples_1d_array: Vec::new(),
            samples_2d_array: Vec::new(),
            array_1d_offset: 0_usize,
            array_2d_offset: 0_usize,
        };
        for _i in 0..lds.n_sampled_dimensions {
            let additional_1d: Vec<Float> = vec![0.0; lds.samples_per_pixel as usize];
            let additional_2d: Vec<Point2f> =
                vec![Point2f::default(); lds.samples_per_pixel as usize];
            lds.samples_1d.push(additional_1d);
            lds.samples_2d.push(additional_2d);
        }
        lds
    }
    pub fn get_camera_sample(&mut self, p_raster: Point2i) -> CameraSample {
        let mut cs: CameraSample = CameraSample::default();
        cs.p_film = Point2f {
            x: p_raster.x as Float,
            y: p_raster.y as Float,
        } + self.get_2d();
        cs.time = self.get_1d();
        cs.p_lens = self.get_2d();
        cs
    }
    pub fn start_pixel(&mut self, p: Point2i) {
        // TODO: ProfilePhase _(Prof::StartPixel);
        // generate 1D and 2D pixel sample components using $(0,2)$-sequence
        for samples in &mut self.samples_1d {
            van_der_corput(1, self.samples_per_pixel as i32, samples, &mut self.rng);
        }
        for samples in &mut self.samples_2d {
            sobol_2d(1, self.samples_per_pixel as i32, samples, &mut self.rng);
        }
        // generate 1D and 2D array samples using $(0,2)$-sequence
        for i in 0..self.samples_1d_array_sizes.len() {
            let samples: &mut [Float] = self.samples_1d_array[i].as_mut_slice();
            van_der_corput(self.samples_1d_array_sizes[i],
                           self.samples_per_pixel as i32,
                           samples,
                           &mut self.rng);
        }
        for i in 0..self.samples_2d_array_sizes.len() {
            let samples: &mut [Point2f] = self.samples_2d_array[i].as_mut_slice();
            sobol_2d(self.samples_2d_array_sizes[i],
                     self.samples_per_pixel as i32,
                     samples,
                     &mut self.rng);
        }
        // PixelSampler::StartPixel(p);
        self.current_pixel = p;
        self.current_pixel_sample_index = 0_i64;
        // reset array offsets for next pixel sample
        self.array_1d_offset = 0_usize;
        self.array_2d_offset = 0_usize;
    }
    pub fn round_count(&self, count: i32) -> i32 {
        round_up_pow2_32(count)
    }
    // PixelSampler public methods
    pub fn start_next_sample(&mut self) -> bool {
        self.current_1d_dimension = 0_i32;
        self.current_2d_dimension = 0_i32;
        // Sampler::StartNextSample()
        // reset array offsets for next pixel sample
        self.array_1d_offset = 0_usize;
        self.array_2d_offset = 0_usize;
        self.current_pixel_sample_index += 1_i64;
        self.current_pixel_sample_index < self.samples_per_pixel
    }
    // inherited from class Sampler (see sampler.h)
    pub fn request_2d_array(&mut self, n: i32) {
        assert_eq!(self.round_count(n), n);
        self.samples_2d_array_sizes.push(n);
        let size: usize = (n * self.samples_per_pixel as i32) as usize;
        let additional_points: Vec<Point2f> = vec![Point2f::default(); size];
        self.samples_2d_array.push(additional_points);
    }
    pub fn get_2d_array(&mut self, n: i32) -> Vec<Point2f> {
        let mut samples: Vec<Point2f> = Vec::new();
        if self.array_2d_offset == self.samples_2d_array.len() {
            return samples;
        }
        assert_eq!(self.samples_2d_array_sizes[self.array_2d_offset], n);
        assert!(self.current_pixel_sample_index < self.samples_per_pixel,
                "self.current_pixel_sample_index ({}) < self.samples_per_pixel ({})",
                self.current_pixel_sample_index,
                self.samples_per_pixel);
        let start: usize = (self.current_pixel_sample_index * n as i64) as usize;
        let end: usize = start + n as usize;
        samples = self.samples_2d_array[self.array_2d_offset][start..end].to_vec();
        self.array_2d_offset += 1;
        samples
    }
    pub fn get_current_sample_number(&self) -> i64 {
        self.current_pixel_sample_index
    }
}

impl Sampler for ZeroTwoSequenceSampler {
    fn get_1d(&mut self) -> Float {
        // TODO: ProfilePhase _(Prof::GetSample);
        assert!(self.current_pixel_sample_index < self.samples_per_pixel,
                "current_pixel_sample_index = {}, samples_per_pixel = {}",
                self.current_pixel_sample_index,
                self.samples_per_pixel);
        if self.current_1d_dimension < self.samples_1d.len() as i32 {
            let sample: Float = self.samples_1d[self.current_1d_dimension
                                                as usize][self.current_pixel_sample_index
                                                          as usize];
            self.current_1d_dimension += 1;
            sample
        } else {
            self.rng.uniform_float()
        }
    }
    fn get_2d(&mut self) -> Point2f {
        // TODO: ProfilePhase _(Prof::GetSample);
        assert!(self.current_pixel_sample_index < self.samples_per_pixel,
                "current_pixel_sample_index = {}, samples_per_pixel = {}",
                self.current_pixel_sample_index,
                self.samples_per_pixel);
        if self.current_2d_dimension < self.samples_2d.len() as i32 {
            let sample: Point2f = self.samples_2d[self.current_2d_dimension
                                                  as usize][self.current_pixel_sample_index
                                                            as usize];
            self.current_2d_dimension += 1;
            sample
        } else {
            // C++ call order for Point2f(rng.UniformFloat(), rng.UniformFloat());
            let y = self.rng.uniform_float();
            let x = self.rng.uniform_float();
            Point2f {
                x: x,
                y: y,
            }
        }
    }
}

// see lowdiscrepancy.h

/// The bits of an integer quantity can be efficiently reversed with a
/// series of logical bit operations.
pub fn reverse_bits_32(n: u32) -> u32 {
    let mut n = (n << 16) | (n >> 16);
    n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8);
    n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4);
    n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2);
    n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1);
    n
}

/// The bits of a 64-bit value can be reversed by reversing the two
/// 32-bit components individually and then interchanging them.
pub fn reverse_bits_64(n: u64) -> u64 {
    let n0: u64 = reverse_bits_32(n as u32) as u64;
    let n1: u64 = reverse_bits_32((n >> 32) as u32) as u64;
    (n0 << 32) | n1
}

/// Takes a generator matrix *c*, a number of 1D samples to generate
/// *n*, and stores the corresponding samples in memory at the
/// location pointed to by *p*.
pub fn gray_code_sample_1d(c: [u32; 32], n: u32, scramble: u32, p: &mut [Float]) {
    let mut v: u32 = scramble;
    for i in 0..n as usize {
        // 1/2^32
        p[i] = (v as Float * 2.3283064365386963e-10 as Float).min(ONE_MINUS_EPSILON);
        v ^= c[(i + 1).trailing_zeros() as usize];
    }
}

/// Takes two generator matrices *c0* and *c1*, a number of 2D samples
/// to generate *n*, and stores the corresponding samples in memory at
/// the location pointed to by *p*.
pub fn gray_code_sample_2d(c0: &[u32], c1: &[u32], n: u32, scramble: &Point2i, p: &mut [Point2f]) {
    let mut v: [u32; 2] = [scramble.x as u32, scramble.y as u32];
    for i in 0..n as usize {
        p[i].x = (v[0] as Float * 2.3283064365386963e-10 as Float).min(ONE_MINUS_EPSILON);
        p[i].y = (v[1] as Float * 2.3283064365386963e-10 as Float).min(ONE_MINUS_EPSILON);
        v[0] ^= c0[(i + 1).trailing_zeros() as usize];
        v[1] ^= c1[(i + 1).trailing_zeros() as usize];
    }
}

/// Generates a number of scrambled 1D sample values using the Gray
/// code-based sampling machinery.
pub fn van_der_corput(n_samples_per_pixel_sample: i32,
                      n_pixel_samples: i32,
                      samples: &mut [Float],
                      rng: &mut Rng) {
    let scramble: u32 = rng.uniform_uint32();
    let c_van_der_corput: [u32; 32] = [0x80000000, 0x40000000, 0x20000000, 0x10000000, 0x8000000,
                                       0x4000000, 0x2000000, 0x1000000, 0x800000, 0x400000,
                                       0x200000, 0x100000, 0x80000, 0x40000, 0x20000, 0x10000,
                                       0x8000, 0x4000, 0x2000, 0x1000, 0x800, 0x400, 0x200, 0x100,
                                       0x80, 0x40, 0x20, 0x10, 0x8, 0x4, 0x2, 0x1];
    let total_samples: i32 = n_samples_per_pixel_sample * n_pixel_samples;
    gray_code_sample_1d(c_van_der_corput, total_samples as u32, scramble, samples);
    // randomly shuffle 1D sample points
    for i in 0..n_pixel_samples as usize {
        shuffle(&mut samples[(i * n_samples_per_pixel_sample as usize)..],
                n_samples_per_pixel_sample,
                1,
                rng);
    }
    shuffle(&mut samples[..],
            n_pixel_samples,
            n_samples_per_pixel_sample,
            rng);
}

/// Similar to *van_der_corput()*, but uses two generator matrices to
/// generate the first two dimensions of Sobol' points.
pub fn sobol_2d(n_samples_per_pixel_sample: i32,
                n_pixel_samples: i32,
                samples: &mut [Point2f],
                rng: &mut Rng) {
    let x: i32 = rng.uniform_uint32() as i32;
    let y: i32 = rng.uniform_uint32() as i32;
    let scramble: Point2i = Point2i { x: x, y: y };
    // define 2D Sobol$'$ generator matrices _c_sobol[2]_
    let c_sobol: [[u32; 32]; 2] =
        [[0x80000000_u32,
          0x40000000,
          0x20000000,
          0x10000000,
          0x8000000,
          0x4000000,
          0x2000000,
          0x1000000,
          0x800000,
          0x400000,
          0x200000,
          0x100000,
          0x80000,
          0x40000,
          0x20000,
          0x10000,
          0x8000,
          0x4000,
          0x2000,
          0x1000,
          0x800,
          0x400,
          0x200,
          0x100,
          0x80,
          0x40,
          0x20,
          0x10,
          0x8,
          0x4,
          0x2,
          0x1],
         [0x80000000, 0xc0000000, 0xa0000000, 0xf0000000, 0x88000000, 0xcc000000, 0xaa000000,
          0xff000000, 0x80800000, 0xc0c00000, 0xa0a00000, 0xf0f00000, 0x88880000, 0xcccc0000,
          0xaaaa0000, 0xffff0000, 0x80008000, 0xc000c000, 0xa000a000, 0xf000f000, 0x88008800,
          0xcc00cc00, 0xaa00aa00, 0xff00ff00, 0x80808080, 0xc0c0c0c0, 0xa0a0a0a0, 0xf0f0f0f0,
          0x88888888, 0xcccccccc, 0xaaaaaaaa, 0xffffffff]];
    gray_code_sample_2d(&c_sobol[0],
                        &c_sobol[1],
                        (n_samples_per_pixel_sample * n_pixel_samples) as u32,
                        &scramble,
                        &mut samples[..]);
    for i in 0..n_pixel_samples as usize {
        shuffle(&mut samples[(i * n_samples_per_pixel_sample as usize)..],
                n_samples_per_pixel_sample,
                1,
                rng);
    }
    shuffle(&mut samples[..],
            n_pixel_samples,
            n_samples_per_pixel_sample,
            rng);
}

// see lowdiscrepancy.cpp

/// Once we have an appropriate prime number, use it to compute the
/// radical inverse.
pub fn radical_inverse_specialized(base: u16, a: u64) -> Float {
    let inv_base: Float = 1.0 as Float / base as Float;
    let mut reversed_digits: u64 = 0_u64;
    let mut inv_base_n: Float = 1.0 as Float;
    let mut a: u64 = a; // shadowing input parameter
    while a != 0_u64 {
        let next: u64 = a / base as u64;
        let digit: u64 = a - next * base as u64;
        reversed_digits = reversed_digits * base as u64 + digit;
        inv_base_n *= inv_base;
        a = next;
    }
    assert!(reversed_digits as Float * inv_base_n < 1.00001 as Float);
    (reversed_digits as Float * inv_base_n).min(ONE_MINUS_EPSILON)
}

/// Map to an appropriate prime number and delegate to another
/// function to compute the radical inverse.
pub fn radical_inverse(base_index: u16, a: u64) -> Float {
    match base_index {
        0 => {
            // TODO: #ifndef PBRT_HAVE_HEX_FP_CONSTANTS
            // 0x1p-64 = (2.0 as Float).powi(-64 as i32)
            return reverse_bits_64(a) as Float * (2.0 as Float).powi(-64 as i32);
        },
        1 => {
            return radical_inverse_specialized(3_u16, a);
        },
        2 => {
            return radical_inverse_specialized(5_u16, a);
        },
        3 => {
            return radical_inverse_specialized(7_u16, a);
        },
        4 => {
            return radical_inverse_specialized(11_u16, a);
        },
        5 => {
            return radical_inverse_specialized(13_u16, a);
        },
        6 => {
            return radical_inverse_specialized(17_u16, a);
        },
        // WORK
        _ => {
            panic!("TODO: radical_inverse({:?}, {:?})", base_index, a);
        },
    };
}

// see filter.h

pub trait Filter {
    fn evaluate(&self, p: Point2f) -> Float;
    fn get_radius(&self) -> Vector2f;
}

// see box.h

#[derive(Debug,Default,Copy,Clone)]
pub struct BoxFilter {
    // inherited from Filter (see filter.h)
    pub radius: Vector2f,
    pub inv_radius: Vector2f,
}

impl Filter for BoxFilter {
    fn evaluate(&self, _p: Point2f) -> Float {
        1.0
    }
    fn get_radius(&self) -> Vector2f {
        Vector2f {
            x: self.radius.x,
            y: self.radius.y,
        }
    }
}

// see gaussian.h

#[derive(Debug,Default,Copy,Clone)]
pub struct GaussianFilter {
    pub alpha: Float,
    pub exp_x: Float,
    pub exp_y: Float,
    // inherited from Filter (see filter.h)
    pub radius: Vector2f,
    pub inv_radius: Vector2f,
}

impl GaussianFilter {
    pub fn gaussian(&self, d: Float, expv: Float) -> Float {
        (0.0 as Float).max((-self.alpha * d * d).exp() - expv)
    }
}

impl Filter for GaussianFilter {
    fn evaluate(&self, p: Point2f) -> Float {
        self.gaussian(p.x, self.exp_x) * self.gaussian(p.y, self.exp_y)
    }
    fn get_radius(&self) -> Vector2f {
        Vector2f {
            x: self.radius.x,
            y: self.radius.y,
        }
    }
}

// see film.h

const FILTER_TABLE_WIDTH: usize = 16;

#[derive(Debug,Default,Copy,Clone)]
pub struct Pixel {
    xyz: [Float; 3],
    filter_weight_sum: Float,
    splat_xyz: [AtomicFloat; 3],
    pad: Float,
}

#[derive(Debug,Default,Copy,Clone)]
pub struct FilmTilePixel {
    contrib_sum: Spectrum,
    filter_weight_sum: Float,
}

pub struct FilmTile<'a> {
    pub pixel_bounds: Bounds2i,
    filter_radius: Vector2f,
    inv_filter_radius: Vector2f,
    filter_table: &'a [Float; FILTER_TABLE_WIDTH * FILTER_TABLE_WIDTH],
    filter_table_size: usize,
    pixels: Vec<FilmTilePixel>,
    max_sample_luminance: Float,
}

impl<'a> FilmTile<'a> {
    pub fn new(pixel_bounds: Bounds2i,
               filter_radius: Vector2f,
               filter_table: &'a [Float; FILTER_TABLE_WIDTH * FILTER_TABLE_WIDTH],
               filter_table_size: usize,
               max_sample_luminance: Float)
               -> Self {
        FilmTile {
            pixel_bounds: pixel_bounds,
            filter_radius: filter_radius,
            inv_filter_radius: Vector2f {
                x: 1.0 / filter_radius.x,
                y: 1.0 / filter_radius.y,
            },
            filter_table: filter_table,
            filter_table_size: filter_table_size,
            // TODO: pixels = std::vector<FilmTilePixel>(std::max(0, pixelBounds.Area()));
            pixels: vec![FilmTilePixel::default(); pixel_bounds.area() as usize],
            max_sample_luminance: max_sample_luminance,
        }
    }
    pub fn add_sample(&mut self, p_film: Point2f, l: &mut Spectrum, sample_weight: Float) {
        // TODO: ProfilePhase _(Prof::AddFilmSample);
        if l.y() > self.max_sample_luminance {
            *l *= Spectrum::new(self.max_sample_luminance / l.y());
        }
        // compute sample's raster bounds
        let p_film_discrete: Point2f = p_film - Vector2f{ x: 0.5, y: 0.5, };
        let p0f: Point2f = pnt2_ceil(p_film_discrete - self.filter_radius);
        let mut p0: Point2i = Point2i { x: p0f.x as i32, y: p0f.y as i32, };
        let p1f: Point2f = pnt2_floor(p_film_discrete + self.filter_radius);
        let mut p1: Point2i = Point2i { x: p1f.x as i32 + 1, y: p1f.y as i32 + 1, };
        p0 = pnt2_max_pnt2(p0, self.pixel_bounds.p_min);
        p1 = pnt2_min_pnt2(p1, self.pixel_bounds.p_max);

        // loop over filter support and add sample to pixel arrays

        // precompute $x$ and $y$ filter table offsets
        let mut ifx: Vec<usize> = Vec::with_capacity(p1.x as usize - p0.x as usize);
        for x in p0.x..p1.x {
            let fx: Float = ((x as Float - p_film_discrete.x) * self.inv_filter_radius.x *
                      self.filter_table_size as Float)
                .abs();
            ifx.push(fx.floor().min(self.filter_table_size as Float - 1.0) as usize);
        }
        let mut ify: Vec<usize> = Vec::with_capacity(p1.y as usize - p0.y as usize);
        for y in p0.y..p1.y {
            let fy: Float = ((y as Float - p_film_discrete.y) * self.inv_filter_radius.y *
                             self.filter_table_size as Float)
                .abs();
            ify.push(fy.floor().min(self.filter_table_size as Float - 1.0) as usize);
        }
        for y in p0.y..p1.y {
            for x in p0.x..p1.x {
                // evaluate filter value at $(x,y)$ pixel
                let offset: usize = ify[(y - p0.y) as usize] * self.filter_table_size + ifx[(x - p0.x) as usize];
                let filter_weight: Float = self.filter_table[offset];
                // update pixel values with filtered sample contribution
                let idx = self.get_pixel_index(x, y);
                let ref mut pixel = self.pixels[idx];
                pixel.contrib_sum += *l * Spectrum::new(sample_weight) * Spectrum::new(filter_weight);
                pixel.filter_weight_sum += filter_weight;
            }
        }
    }
    fn get_pixel_index(&self, x: i32, y: i32) -> usize {
        let width: i32 = self.pixel_bounds.p_max.x - self.pixel_bounds.p_min.x;
        let pidx = (y - self.pixel_bounds.p_min.y) * width + (x - self.pixel_bounds.p_min.x);
        pidx as usize
    }
}

pub struct Film {
    // Film Public Data
    /// The overall resolution of the image in pixels
    pub full_resolution: Point2i,
    /// The length of the diagonal of the film's physical area (specified in mm, stored in meters)
    pub diagonal: Float,
    /// A filter function
    pub filter: Arc<Filter + Sync + Send>,
    /// The filename of the output image
    pub filename: String,
    /// A crop window that may specify a subset of the image to render
    pub cropped_pixel_bounds: Bounds2i,

    // Film Private Data
    pub pixels: RwLock<Vec<Pixel>>,
    filter_table: [Float; FILTER_TABLE_WIDTH * FILTER_TABLE_WIDTH],
    scale: Float,
    max_sample_luminance: Float,
}

impl Film {
    pub fn new(resolution: Point2i,
               crop_window: Bounds2f,
               filter: Arc<Filter + Sync + Send>,
               diagonal: Float,
               filename: String,
               scale: Float,
               max_sample_luminance: Float)
               -> Self {
        let cropped_pixel_bounds: Bounds2i = Bounds2i {
            p_min: Point2i {
                x: (resolution.x as Float * crop_window.p_min.x).ceil() as i32,
                y: (resolution.y as Float * crop_window.p_min.y).ceil() as i32,
            },
            p_max: Point2i {
                x: (resolution.x as Float * crop_window.p_max.x).ceil() as i32,
                y: (resolution.y as Float * crop_window.p_max.y).ceil() as i32,
            },
        };
        // allocate film image storage
        // let pixels: Vec<Pixel> = vec![Pixel::default(); cropped_pixel_bounds.area() as usize];
        // precompute filter weight table
        let mut filter_table: [Float; FILTER_TABLE_WIDTH * FILTER_TABLE_WIDTH] =
            [0.0; FILTER_TABLE_WIDTH * FILTER_TABLE_WIDTH];
        let mut offset: usize = 0;
        let filter_radius: Vector2f = filter.get_radius();
        for y in 0..FILTER_TABLE_WIDTH {
            for x in 0..FILTER_TABLE_WIDTH {
                let p: Point2f = Point2f {
                    x: (x as Float + 0.5) * filter_radius.x / FILTER_TABLE_WIDTH as Float,
                    y: (y as Float + 0.5) * filter_radius.y / FILTER_TABLE_WIDTH as Float,
                };
                filter_table[offset] = filter.evaluate(p);
                offset += 1;
            }
        }
        Film {
            full_resolution: resolution,
            diagonal: diagonal * 0.001,
            filter: filter,
            filename: filename,
            cropped_pixel_bounds: cropped_pixel_bounds,
            pixels: RwLock::new(vec![Pixel::default(); cropped_pixel_bounds.area() as usize]),
            filter_table: filter_table,
            scale: scale,
            max_sample_luminance: max_sample_luminance,
        }
    }
    pub fn get_sample_bounds(&self) -> Bounds2i {
        let f: Point2f = pnt2_floor(Point2f {
            x: self.cropped_pixel_bounds.p_min.x as Float,
            y: self.cropped_pixel_bounds.p_min.y as Float,
        } + Vector2f { x: 0.5, y: 0.5 } - self.filter.get_radius());
        let c: Point2f = pnt2_ceil(Point2f {
            x: self.cropped_pixel_bounds.p_max.x as Float,
            y: self.cropped_pixel_bounds.p_max.y as Float,
        } - Vector2f { x: 0.5, y: 0.5 } + self.filter.get_radius());
        let float_bounds: Bounds2f = Bounds2f {
            p_min: f,
            p_max: c,
        };
        Bounds2i {
            p_min: Point2i {
                x: float_bounds.p_min.x as i32,
                y: float_bounds.p_min.y as i32,
            },
            p_max: Point2i {
                x: float_bounds.p_max.x as i32,
                y: float_bounds.p_max.y as i32,
            },
        }
    }
    pub fn get_film_tile(&self, sample_bounds: Bounds2i) -> FilmTile {
        // bound image pixels that samples in _sample_bounds_ contribute to
        let half_pixel: Vector2f = Vector2f { x: 0.5, y: 0.5 };
        let float_bounds: Bounds2f = Bounds2f {
            p_min: Point2f {
                x: sample_bounds.p_min.x as Float,
                y: sample_bounds.p_min.y as Float,
            },
            p_max: Point2f {
                x: sample_bounds.p_max.x as Float,
                y: sample_bounds.p_max.y as Float,
            },
        };
        let p_min: Point2f = float_bounds.p_min - half_pixel - self.filter.get_radius();
        let p0: Point2i = Point2i {
            x: p_min.x.ceil() as i32,
            y: p_min.y.ceil() as i32,
        };
        let p_max: Point2f = float_bounds.p_max - half_pixel + self.filter.get_radius();
        let p1: Point2i = Point2i {
            x: p_max.x.floor() as i32,
            y: p_max.y.floor() as i32,
        } + Point2i { x: 1, y: 1 };
        let tile_pixel_bounds: Bounds2i = bnd2_intersect_bnd2(Bounds2i {
                                                                  p_min: p0,
                                                                  p_max: p1,
                                                              },
                                                              self.cropped_pixel_bounds);
        FilmTile::new(tile_pixel_bounds,
                      self.filter.get_radius(),
                      &self.filter_table,
                      FILTER_TABLE_WIDTH,
                      self.max_sample_luminance)
    }
    pub fn merge_film_tile(&self, tile: &FilmTile) {
        // TODO: ProfilePhase p(Prof::MergeFilmTile);
        // println!("Merging film tile {:?}", tile.pixel_bounds);
        // TODO: std::lock_guard<std::mutex> lock(mutex);
        for pixel in &tile.pixel_bounds {
            // merge _pixel_ into _Film::pixels_
            let idx = tile.get_pixel_index(pixel.x, pixel.y);
            let ref tile_pixel = tile.pixels[idx];
            // START let mut merge_pixel: &mut Pixel = self.get_pixel_mut(pixel);
            assert!(pnt2_inside_exclusive(pixel, self.cropped_pixel_bounds));
            let width: i32 = self.cropped_pixel_bounds.p_max.x - self.cropped_pixel_bounds.p_min.x;
            let offset: i32 = (pixel.x - self.cropped_pixel_bounds.p_min.x) +
                (pixel.y - self.cropped_pixel_bounds.p_min.y) * width;
            let mut pixels_write = self.pixels.write().unwrap();
            let mut merge_pixel = pixels_write[offset as usize];
            // END let mut merge_pixel: &mut Pixel = self.get_pixel_mut(pixel);
            let mut xyz: [Float; 3] = [0.0; 3];
            tile_pixel.contrib_sum.to_xyz(&mut xyz);
            for i in 0..3 {
                merge_pixel.xyz[i] += xyz[i];
            }
            merge_pixel.filter_weight_sum += tile_pixel.filter_weight_sum;
            // write pixel back
            pixels_write[offset as usize] = merge_pixel;
        }
    }
    pub fn write_image(&self, splat_scale: Float) {
        println!("Converting image to RGB and computing final weighted pixel values");
        let mut rgb: Vec<Float> =
            vec![0.0 as Float; (3 * self.cropped_pixel_bounds.area()) as usize];
        let mut exr: Vec<(Float, Float, Float)> = // copy data for OpenEXR image
            vec![(0.0_f32, 0.0_f32, 0.0_f32); self.cropped_pixel_bounds.area() as usize];
        let mut offset: usize = 0;
        for p in &self.cropped_pixel_bounds {
            // convert pixel XYZ color to RGB
            let pixel: Pixel = self.get_pixel(p);
            let start = 3 * offset;
            let mut rgb_array: [Float; 3] = [0.0 as Float; 3];
            xyz_to_rgb(&pixel.xyz, &mut rgb_array); // TODO: Use 'rgb' directly.
            rgb[start + 0] = rgb_array[0];
            rgb[start + 1] = rgb_array[1];
            rgb[start + 2] = rgb_array[2];
            // normalize pixel with weight sum
            let filter_weight_sum: Float = pixel.filter_weight_sum;
            if filter_weight_sum != 0.0 as Float {
                let inv_wt: Float = 1.0 as Float / filter_weight_sum;
                rgb[start + 0] = (rgb[start + 0] * inv_wt).max(0.0 as Float);
                rgb[start + 1] = (rgb[start + 1] * inv_wt).max(0.0 as Float);
                rgb[start + 2] = (rgb[start + 2] * inv_wt).max(0.0 as Float);
            }
            // add splat value at pixel
            let mut splat_rgb: [Float; 3] = [0.0 as Float; 3];
            let splat_xyz: [Float; 3] = [Float::from(pixel.splat_xyz[0]),
                                         Float::from(pixel.splat_xyz[1]),
                                         Float::from(pixel.splat_xyz[2])];
            xyz_to_rgb(&splat_xyz, &mut splat_rgb);
            rgb[start + 0] += splat_scale * splat_rgb[0];
            rgb[start + 1] += splat_scale * splat_rgb[1];
            rgb[start + 2] += splat_scale * splat_rgb[2];
            // scale pixel value by _scale_
            rgb[start + 0] *= self.scale;
            rgb[start + 1] *= self.scale;
            rgb[start + 2] *= self.scale;
            // copy data for OpenEXR image
            exr[offset].0 = rgb[start + 0];
            exr[offset].1 = rgb[start + 1];
            exr[offset].2 = rgb[start + 2];
            offset += 1;
        }
        let filename = "pbrt.png";
        println!("Writing image {:?} with bounds {:?}",
                 filename, // TODO: self.filename,
                 self.cropped_pixel_bounds);
        // TODO: pbrt::WriteImage(filename, &rgb[0], croppedPixelBounds, fullResolution);
        let mut buffer: Vec<u8> = vec![0.0 as u8; (3 * self.cropped_pixel_bounds.area()) as usize];
        // 8-bit format; apply gamma (see WriteImage(...) in imageio.cpp)
        let width: u32 = (self.cropped_pixel_bounds.p_max.x -
                          self.cropped_pixel_bounds.p_min.x) as u32;
        let height: u32 = (self.cropped_pixel_bounds.p_max.y -
                           self.cropped_pixel_bounds.p_min.y) as u32;
        // OpenEXR
        let filename = "pbrt_rust.exr";
        println!("Writing image {:?} with bounds {:?}",
                 filename, // TODO: self.filename,
                 self.cropped_pixel_bounds);
        let mut file = std::fs::File::create("pbrt_rust.exr").unwrap();
        let mut output_file = ScanlineOutputFile::new(&mut file,
                                                      Header::new()
                                                      .set_resolution(width, height)
                                                      .add_channel("R", PixelType::FLOAT)
                                                      .add_channel("G", PixelType::FLOAT)
                                                      .add_channel("B", PixelType::FLOAT)).unwrap();
        let mut fb = FrameBuffer::new(width as u32, height as u32);
        fb.insert_channels(&["R", "G", "B"], &exr);
        output_file.write_pixels(&fb).unwrap();

        // OpenEXR
        for y in 0..height {
            for x in 0..width {
                // red
                let index: usize = (3 * (y * width + x) + 0) as usize;
                buffer[index] = clamp_t(255.0 as Float * gamma_correct(rgb[index]) + 0.5,
                                        0.0 as Float,
                                        255.0 as Float) as u8;
                // green
                let index: usize = (3 * (y * width + x) + 1) as usize;
                buffer[index] = clamp_t(255.0 as Float * gamma_correct(rgb[index]) + 0.5,
                                        0.0 as Float,
                                        255.0 as Float) as u8;
                // blue
                let index: usize = (3 * (y * width + x) + 2) as usize;
                buffer[index] = clamp_t(255.0 as Float * gamma_correct(rgb[index]) + 0.5,
                                        0.0 as Float,
                                        255.0 as Float) as u8;
            }
        }
        // write "pbrt.png" to disk
        image::save_buffer(&Path::new("pbrt.png"),
                           &buffer,
                           width,
                           height,
                           image::RGB(8))
            .unwrap();
    }
    pub fn get_pixel(&self, p: Point2i) -> Pixel {
        assert!(pnt2_inside_exclusive(p, self.cropped_pixel_bounds));
        let width: i32 = self.cropped_pixel_bounds.p_max.x - self.cropped_pixel_bounds.p_min.x;
        let offset: i32 = (p.x - self.cropped_pixel_bounds.p_min.x) +
                          (p.y - self.cropped_pixel_bounds.p_min.y) * width;
        self.pixels.read().unwrap()[offset as usize]
    }
}

// see camera.h

pub trait Camera {
    fn generate_ray_differential(&self, sample: &CameraSample, ray: &mut Ray) -> Float;
    fn get_film(&self) -> Arc<Film>;
}

#[derive(Debug,Default,Copy,Clone)]
pub struct CameraSample {
    pub p_film: Point2f,
    pub p_lens: Point2f,
    pub time: Float,
}

// see perspective.h

pub struct PerspectiveCamera {
    // inherited from Camera (see camera.h)
    pub camera_to_world: AnimatedTransform,
    pub shutter_open: Float,
    pub shutter_close: Float,
    pub film: Arc<Film>,
    // TODO: const Medium *medium;
    // inherited from ProjectiveCamera (see camera.h)
    // camera_to_screen: Transform,
    raster_to_camera: Transform,
    // screen_to_raster: Transform,
    // raster_to_screen: Transform,
    // lens_radius: Float,
    // focal_distance: Float,
    // private data (see perspective.h)
    dx_camera: Vector3f,
    dy_camera: Vector3f,
    // a: Float,
}

impl PerspectiveCamera {
    pub fn new(camera_to_world: AnimatedTransform,
               screen_window: Bounds2f,
               shutter_open: Float,
               shutter_close: Float,
               _lens_radius: Float,
               _focal_distance: Float,
               fov: Float,
               film: Arc<Film>,
               /* const Medium *medium */)
               -> Self {
        // see perspective.cpp
        let camera_to_screen: Transform = Transform::perspective(fov, 1e-2, 1000.0);
        // see camera.h
        // compute projective camera screen transformations
        let scale1 = Transform::scale(film.full_resolution.x as Float,
                                      film.full_resolution.y as Float,
                                      1.0);
        let scale2 = Transform::scale(1.0 / (screen_window.p_max.x - screen_window.p_min.x),
                                      1.0 / (screen_window.p_min.y - screen_window.p_max.y),
                                      1.0);
        let translate = Transform::translate(Vector3f {
            x: -screen_window.p_min.x,
            y: -screen_window.p_max.y,
            z: 0.0,
        });
        let screen_to_raster = scale1 * scale2 * translate;
        let raster_to_screen = Transform::inverse(screen_to_raster);
        let raster_to_camera = Transform::inverse(camera_to_screen) * raster_to_screen;
        // see perspective.cpp
        // compute differential changes in origin for perspective camera rays
        let dx_camera: Vector3f = raster_to_camera.transform_point(Point3f {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        }) -
                                  raster_to_camera.transform_point(Point3f {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        });
        let dy_camera: Vector3f = raster_to_camera.transform_point(Point3f {
            x: 0.0,
            y: 1.0,
            z: 0.0,
        }) -
                                  raster_to_camera.transform_point(Point3f {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        });
        // compute image plane bounds at $z=1$ for _PerspectiveCamera_
        let res: Point2i = film.full_resolution;
        let mut p_min: Point3f = raster_to_camera.transform_point(Point3f {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        });
        // Point3f p_max = RasterToCamera(Point3f(res.x, res.y, 0));
        let mut p_max: Point3f = raster_to_camera.transform_point(Point3f {
            x: res.x as Float,
            y: res.y as Float,
            z: 0.0,
        });
        p_min /= p_min.z;
        p_max /= p_max.z;
        // let a: Float = ((p_max.x - p_min.x) * (p_max.y - p_min.y)).abs();

        PerspectiveCamera {
            camera_to_world: camera_to_world,
            shutter_open: shutter_open,
            shutter_close: shutter_close,
            film: film,
            // camera_to_screen: camera_to_screen,
            raster_to_camera: raster_to_camera,
            // screen_to_raster: screen_to_raster,
            // raster_to_screen: raster_to_screen,
            // lens_radius: lens_radius,
            // focal_distance: focal_distance,
            dx_camera: dx_camera,
            dy_camera: dy_camera,
            // a: a,
        }
    }
}

impl Camera for PerspectiveCamera {
    fn generate_ray_differential(&self, sample: &CameraSample, ray: &mut Ray) -> Float {
        // TODO: ProfilePhase prof(Prof::GenerateCameraRay);
        // compute raster and camera sample positions
        let p_film: Point3f = Point3f {
            x: sample.p_film.x,
            y: sample.p_film.y,
            z: 0.0,
        };
        let p_camera: Point3f = self.raster_to_camera.transform_point(p_film);
        let dir: Vector3f = vec3_normalize(Vector3f {
            x: p_camera.x,
            y: p_camera.y,
            z: p_camera.z,
        });
        // *ray = RayDifferential(Point3f(0, 0, 0), dir);
        ray.o = Point3f::default();
        ray.d = dir;
        ray.t_max = std::f32::INFINITY;
        ray.time = lerp(sample.time, self.shutter_open, self.shutter_close);
        // TODO: modify ray for depth of field
        // TODO: if (lensRadius > 0) { ... } else {
        let diff: RayDifferential = RayDifferential {
            rx_origin: ray.o,
            ry_origin: ray.o,
            rx_direction: vec3_normalize(Vector3f {
                x: p_camera.x,
                y: p_camera.y,
                z: p_camera.z,
            } + self.dx_camera),
            ry_direction: vec3_normalize(Vector3f {
                x: p_camera.x,
                y: p_camera.y,
                z: p_camera.z,
            } + self.dy_camera),
        };
        // TODO: ray->medium = medium;
        // TODO: *ray = CameraToWorld(*ray);
        // ray->hasDifferentials = true;
        ray.differential = Some(diff);
        self.camera_to_world.transform_ray(ray);
        1.0
    }
    fn get_film(&self) -> Arc<Film> {
        self.film.clone()
    }
}

// see reflection.h

pub struct Bsdf {
    pub eta: Float,
    /// shading normal
    pub ns: Normal3f,
    /// geometric normal
    pub ng: Normal3f,
    pub ss: Vector3f,
    pub ts: Vector3f,
    pub bxdfs: Vec<Box<Bxdf + Sync + Send>>,
}

impl Bsdf {
    pub fn new(si: &SurfaceInteraction, eta: Float, bxdfs: Vec<Box<Bxdf + Sync + Send>>) -> Self {
        let ss = vec3_normalize(si.shading.dpdu);
        Bsdf {
            eta: eta,
            ns: si.shading.n,
            ng: si.n,
            ss: ss,
            ts: nrm_cross_vec3(si.shading.n, ss),
            bxdfs: bxdfs,
        }
    }
    pub fn num_components(&self, flags: u8) -> u8 {
        let mut num: u8 = 0;
        let n_bxdfs: usize = self.bxdfs.len();
        for i in 0..n_bxdfs {
            if self.bxdfs[i].matches_flags(flags) {
                num += 1;
            }
        }
        num
    }
    pub fn world_to_local(&self, v: Vector3f) -> Vector3f {
        Vector3f {
            x: vec3_dot_vec3(v, self.ss),
            y: vec3_dot_vec3(v, self.ts),
            z: vec3_dot_vec3(v, Vector3f::from(self.ns)),
        }
    }
    pub fn local_to_world(&self, v: Vector3f) -> Vector3f {
        Vector3f {
            x: self.ss.x * v.x + self.ts.x * v.y + self.ns.x * v.z,
            y: self.ss.y * v.x + self.ts.y * v.y + self.ns.y * v.z,
            z: self.ss.z * v.x + self.ts.z * v.y + self.ns.z * v.z,
        }
    }
    pub fn f(&self, wo_w: Vector3f, wi_w: Vector3f, flags: u8) -> Spectrum {
        // TODO: ProfilePhase pp(Prof::BSDFEvaluation);
        let wi: Vector3f = self.world_to_local(wi_w);
        let wo: Vector3f = self.world_to_local(wo_w);
        if wo.z == 0.0 as Float {
            return Spectrum::new(0.0 as Float);
        }
        let reflect: bool = (vec3_dot_vec3(wi_w, Vector3f::from(self.ng)) *
                             vec3_dot_vec3(wo_w, Vector3f::from(self.ng))) >
                            0.0 as Float;
        let mut f: Spectrum = Spectrum::new(0.0 as Float);
        let n_bxdfs: usize = self.bxdfs.len();
        for i in 0..n_bxdfs {
            if self.bxdfs[i].matches_flags(flags) &&
               ((reflect && (self.bxdfs[i].get_type() & BxdfType::BsdfReflection as u8 > 0_u8)) ||
                (!reflect &&
                 (self.bxdfs[i].get_type() & BxdfType::BsdfTransmission as u8 > 0_u8)))
            {
                f += self.bxdfs[i].f(wo, wi);
            }
        }
        f
    }
    /// Calls the individual Bxdf::sample_f() methods to generate samples.
    pub fn sample_f(&self,
                    wo_world: Vector3f,
                    wi_world: &mut Vector3f,
                    u: Point2f,
                    pdf: &mut Float,
                    bsdf_flags: u8,
                    sampled_type: &mut u8)
                    -> Spectrum {
        // TODO: ProfilePhase pp(Prof::BSDFSampling);
        // choose which _BxDF_ to sample
        let matching_comps: u8 = self.num_components(bsdf_flags);
        if matching_comps == 0 {
            *pdf = 0.0 as Float;
            *sampled_type = 0_u8;
            return Spectrum::default();
        }
        let comp: u8 = std::cmp::min((u[0] * matching_comps as Float).floor() as u8,
                                     matching_comps - 1_u8);
        // get _BxDF_ pointer for chosen component
        let mut bxdf: Option<&Box<Bxdf + Sync + Send>> = None;
        let mut count: i8 = comp as i8;
        let n_bxdfs: usize = self.bxdfs.len();
        let mut bxdf_index: usize = 0_usize;
        for i in 0..n_bxdfs {
            let matches: bool = self.bxdfs[i].matches_flags(bsdf_flags);
            if  matches && count == 0 {
                count -= 1_i8;
                bxdf = self.bxdfs.get(i);
                bxdf_index = i;
                break;
            } else {
                // fix count
                if matches {
                    // C++ version does this in a single line:
                    // if (bxdfs[i]->MatchesFlags(type) && count-- == 0)
                    count -= 1_i8;
                }
            }
        }
        if bxdf.is_some() {
            let bxdf = bxdf.unwrap();
            // TODO: println!("BSDF::Sample_f chose comp = {:?} /
            // matching = {:?}, bxdf: {:?}", comp, matching_comps,
            // bxdf);

            // remap _BxDF_ sample _u_ to $[0,1)^2$
            let u_remapped: Point2f = Point2f {
                x: (u[0] * matching_comps as Float - comp as Float).min(ONE_MINUS_EPSILON),
                y: u[1],
            };
            // sample chosen _BxDF_
            let mut wi: Vector3f = Vector3f::default();
            let wo: Vector3f = self.world_to_local(wo_world);
            if wo.z == 0.0 as Float {
                return Spectrum::default();
            }
            *pdf = 0.0 as Float;
            if *sampled_type != 0_u8 {
                *sampled_type = bxdf.get_type();
            }
            let mut f: Spectrum = bxdf.sample_f(wo, &mut wi, u_remapped, pdf, sampled_type);
            // let mut ratio: Spectrum = Spectrum::default();
            // if *pdf > 0.0 as Float {
            //     ratio = f / *pdf;
            // }
            // println!("For wo = {:?}, sampled f = {:?}, pdf = {:?}, ratio = {:?}, wi = {:?}",
            //          wo,
            //          f,
            //          *pdf,
            //          ratio,
            //          wi);
            if *pdf == 0.0 as Float {
                if *sampled_type != 0_u8 {
                    *sampled_type = 0_u8;
                }
                return Spectrum::default();
            }
            *wi_world = self.local_to_world(wi);
            // compute overall PDF with all matching _BxDF_s
            if (bxdf.get_type() & BxdfType::BsdfSpecular as u8 == 0_u8) && matching_comps > 1_u8 {
                for i in 0..n_bxdfs {
                    // instead of self.bxdfs[i] != bxdf we compare stored index
                    if bxdf_index != i && self.bxdfs[i].matches_flags(bsdf_flags) {
                        *pdf += self.bxdfs[i].pdf(wo, wi);
                    }
                }
            }
            if matching_comps > 1_u8 {
                *pdf /= matching_comps as Float;
            }
            // compute value of BSDF for sampled direction
            if (bxdf.get_type() & BxdfType::BsdfSpecular as u8 == 0_u8) && matching_comps > 1_u8 {
                let reflect: bool = vec3_dot_nrm(*wi_world, self.ng) *
                                    vec3_dot_nrm(wo_world, self.ng) >
                                    0.0 as Float;
                f = Spectrum::default();
                for i in 0..n_bxdfs {
                    if self.bxdfs[i].matches_flags(bsdf_flags) &&
                       ((reflect && (bxdf.get_type() & BxdfType::BsdfReflection as u8) != 0_u8) ||
                        (reflect && (bxdf.get_type() & BxdfType::BsdfTransmission as u8) == 0_u8)) {
                        f += self.bxdfs[i].f(wo, wi);
                    }
                }
            }
            // let mut ratio: Spectrum = Spectrum::default();
            // if *pdf > 0.0 as Float {
            //     ratio = f / *pdf;
            // }
            // println!("Overall f = {:?}, pdf = {:?}, ratio = {:?}", f, *pdf, ratio);
            return f;
        } else {
            panic!("CHECK_NOTNULL(bxdf)");
        }
    }
    pub fn pdf(&self, wo_world: Vector3f, wi_world: Vector3f, bsdf_flags: u8) -> Float {
        // TODO: ProfilePhase pp(Prof::BSDFPdf);
        let n_bxdfs: usize = self.bxdfs.len();
        if n_bxdfs == 0 {
            return 0.0 as Float;
        }
        let wo: Vector3f = self.world_to_local(wo_world);
        let wi: Vector3f = self.world_to_local(wi_world);
        if wo.z == 0.0 as Float {
            return 0.0 as Float;
        }
        let mut pdf: Float = 0.0 as Float;
        let mut matching_comps: u8 = 0;
        for i in 0..n_bxdfs {
            if self.bxdfs[i].matches_flags(bsdf_flags) {
                matching_comps += 1;
                pdf += self.bxdfs[i].pdf(wo, wi);
            }
        }
        let mut v: Float = 0.0 as Float;
        if matching_comps > 0 {
            v = pdf / matching_comps as Float;
        }
        v
    }
}

#[repr(u8)]
pub enum BxdfType {
    BsdfReflection = 1,
    BsdfTransmission = 2,
    BsdfDiffuse = 4,
    BsdfGlossy = 8,
    BsdfSpecular = 16,
    BsdfAll = 31,
}

pub trait Bxdf {
    fn matches_flags(&self, t: u8) -> bool {
        self.get_type() & t == self.get_type()
    }
    fn f(&self, wo: Vector3f, wi: Vector3f) -> Spectrum;
    fn sample_f(&self,
                wo: Vector3f,
                wi: &mut Vector3f,
                u: Point2f,
                pdf: &mut Float,
                sampled_type: &mut u8)
                -> Spectrum;
    fn pdf(&self, wo: Vector3f, wi: Vector3f) -> Float {
        if vec3_same_hemisphere_vec3(wo, wi) {
            abs_cos_theta(wi) * INV_PI
        } else {
            0.0 as Float
        }
    }
    fn get_type(&self) -> u8;
}

pub trait Fresnel {
    fn evaluate(&self, cos_theta_i: &mut Float) -> Spectrum;
}

#[derive(Debug,Default,Copy,Clone)]
pub struct FresnelDielectric {
    pub eta_i: Float,
    pub eta_t: Float,
}

impl Fresnel for FresnelDielectric {
    fn evaluate(&self, cos_theta_i: &mut Float) -> Spectrum {
        Spectrum::new(fr_dielectric(cos_theta_i, self.eta_i, self.eta_t))
    }
}

#[derive(Debug,Default,Copy,Clone)]
pub struct FresnelNoOp {
}

impl Fresnel for FresnelNoOp {
    fn evaluate(&self, _cos_theta_i: &mut Float) -> Spectrum {
        Spectrum::new(1.0 as Float)
    }
}

#[derive(Clone)]
pub struct SpecularReflection {
    pub r: Spectrum,
    pub fresnel: Arc<Fresnel + Send + Sync>,
}

impl SpecularReflection {
    pub fn new(r: Spectrum, fresnel: Arc<Fresnel + Send + Sync>) -> Self {
        SpecularReflection {
            r: r,
            fresnel: fresnel,
        }
    }
}

impl Bxdf for SpecularReflection {
    fn f(&self, _wo: Vector3f, _wi: Vector3f) -> Spectrum {
        Spectrum::new(0.0 as Float)
    }
    fn sample_f(&self,
                wo: Vector3f,
                wi: &mut Vector3f,
                _sample: Point2f,
                pdf: &mut Float,
                _sampled_type: &mut u8)
                -> Spectrum {
        // compute perfect specular reflection direction
        *wi = Vector3f { x: -wo.x, y: -wo.y, z: wo.z, };
        *pdf = 1.0 as Float;
        let mut cos_theta_i: Float = cos_theta(*wi);
        self.fresnel.evaluate(&mut cos_theta_i) * self.r / abs_cos_theta(*wi)
    }
    fn get_type(&self) -> u8 {
        BxdfType::BsdfReflection as u8 | BxdfType::BsdfSpecular as u8
    }
}

pub struct SpecularTransmission {
    pub t: Spectrum,
    pub eta_a: Float,
    pub eta_b: Float,
    pub fresnel: FresnelDielectric,
    pub mode: TransportMode,
}

impl SpecularTransmission {
    pub fn new(t: Spectrum, eta_a: Float, eta_b: Float, mode: TransportMode) -> Self {
        SpecularTransmission {
            t: t,
            eta_a: eta_a,
            eta_b: eta_b,
            fresnel: FresnelDielectric { eta_i: eta_a, eta_t: eta_b, },
            mode: mode,
        }
    }
}

impl Bxdf for SpecularTransmission {
    fn f(&self, _wo: Vector3f, _wi: Vector3f) -> Spectrum {
        Spectrum::new(0.0 as Float)
    }
    fn sample_f(&self,
                wo: Vector3f,
                wi: &mut Vector3f,
                _sample: Point2f,
                pdf: &mut Float,
                _sampled_type: &mut u8)
                -> Spectrum {
        // figure out which $\eta$ is incident and which is transmitted
        let entering: bool = cos_theta(wo) > 0.0;
        let mut eta_i: Float = self.eta_b;
        if entering {
            eta_i = self.eta_a;
        }
        let mut eta_t: Float = self.eta_a;
        if entering {
            eta_t = self.eta_b;
        }
        // compute ray direction for specular transmission
        if !refract(wo,
                    nrm_faceforward_vec3(Normal3f { x: 0.0, y: 0.0, z: 1.0 }, wo),
                    eta_i / eta_t, wi) {
            return Spectrum::default();
        }
        *pdf = 1.0;
        let mut ft: Spectrum = self.t * (Spectrum::new(1.0 as Float) - self.fresnel.evaluate(&mut cos_theta(*wi)));
        // account for non-symmetry with transmission to different medium
        if self.mode == TransportMode::Radiance {
            ft *= Spectrum::new((eta_i * eta_i) / (eta_t * eta_t));
        }
        ft / abs_cos_theta(*wi)
    }
    fn get_type(&self) -> u8 {
        BxdfType::BsdfTransmission as u8 | BxdfType::BsdfSpecular as u8
    }
}

pub struct FresnelSpecular {
    pub r: Spectrum,
    pub t: Spectrum,
    pub eta_a: Float,
    pub eta_b: Float,
    pub mode: TransportMode,
}

impl FresnelSpecular {
    pub fn new(r: Spectrum, t: Spectrum, eta_a: Float, eta_b: Float, mode: TransportMode) -> Self {
        FresnelSpecular {
            r: r,
            t: t,
            eta_a: eta_a,
            eta_b: eta_b,
            mode: mode,
        }
    }
}

impl Bxdf for FresnelSpecular {
    fn f(&self, _wo: Vector3f, _wi: Vector3f) -> Spectrum {
        Spectrum::new(0.0 as Float)
    }
    fn sample_f(&self,
                wo: Vector3f,
                wi: &mut Vector3f,
                sample: Point2f,
                pdf: &mut Float,
                sampled_type: &mut u8)
                -> Spectrum {
        let mut ct: Float = cos_theta(wo);
        let f: Float = fr_dielectric(&mut ct, self.eta_a, self.eta_b);
        if sample[0] < f {
            // compute specular reflection for _FresnelSpecular_

            // compute perfect specular reflection direction
            *wi = Vector3f {
                x: -wo.x,
                y: -wo.y,
                z: wo.z,
            };
            if *sampled_type != 0_u8 {
                *sampled_type = self.get_type();
            }
            *pdf = f;
            return self.r * f / abs_cos_theta(*wi);
        } else {
            // compute specular transmission for _FresnelSpecular_

            // figure out which $\eta$ is incident and which is transmitted
            let entering: bool = cos_theta(wo) > 0.0 as Float;
            let eta_i: Float;
            if entering {
                eta_i = self.eta_a;
            } else {
                eta_i = self.eta_b;
            }
            let eta_t: Float;
            if entering {
                eta_t = self.eta_b;
            } else {
                eta_t = self.eta_a;
            }
            // compute ray direction for specular transmission
            if !refract(wo, nrm_faceforward_vec3(Normal3f { x: 0.0, y: 0.0, z: 1.0, }, wo), eta_i / eta_t, wi) {
                return Spectrum::default();
            }
            let mut ft: Spectrum = self.t * (1.0 as Float - f);
            // account for non-symmetry with transmission to different medium
            if self.mode == TransportMode::Radiance {
                ft *= Spectrum::new((eta_i * eta_i) / (eta_t * eta_t));
            }
            if *sampled_type != 0_u8 {
                *sampled_type = self.get_type();
            }
            *pdf = 1.0 as Float - f;
            return ft / abs_cos_theta(*wi);
        }
    }
    fn get_type(&self) -> u8 {
        BxdfType::BsdfReflection as u8 | BxdfType::BsdfTransmission as u8 | BxdfType::BsdfSpecular as u8
    }
}

#[derive(Debug,Default,Copy,Clone)]
pub struct LambertianReflection {
    pub r: Spectrum,
}

impl LambertianReflection {
    pub fn new(r: Spectrum) -> Self {
        LambertianReflection { r: r }
    }
}

impl Bxdf for LambertianReflection {
    fn f(&self, _wo: Vector3f, _wi: Vector3f) -> Spectrum {
        self.r * Spectrum::new(INV_PI)
    }
    fn sample_f(&self,
                wo: Vector3f,
                wi: &mut Vector3f,
                u: Point2f,
                pdf: &mut Float,
                _sampled_type: &mut u8)
                -> Spectrum {
        *wi = cosine_sample_hemisphere(u);
        if wo.z < 0.0 as Float {
            wi.z *= -1.0 as Float;
        }
        *pdf = self.pdf(wo, *wi);
        self.f(wo, *wi)
    }
    fn get_type(&self) -> u8 {
        BxdfType::BsdfDiffuse as u8 | BxdfType::BsdfReflection as u8
    }
}

pub struct OrenNayar {
    pub r: Spectrum,
    pub a: Float,
    pub b: Float,
}

impl OrenNayar {
    pub fn new(r: Spectrum, sigma: Float) -> Self {
        let sigma = radians(sigma);
        let sigma2: Float = sigma * sigma;
        OrenNayar {
            r: r,
            a: 1.0 - (sigma2 / (2.0 * (sigma2 + 0.33))),
            b: 0.45 * sigma2 / (sigma2 + 0.09),
        }
    }
}

impl Bxdf for OrenNayar {
    fn f(&self, wo: Vector3f, wi: Vector3f) -> Spectrum {
        let sin_theta_i: Float = sin_theta(wi);
        let sin_theta_o: Float = sin_theta(wo);
        // compute cosine term of Oren-Nayar model
        let mut max_cos: Float = 0.0 as Float;
        if sin_theta_i > 1.0e-4 && sin_theta_o > 1.0e-4 {
            let sin_phi_i: Float = sin_phi(wi);
            let cos_phi_i: Float = cos_phi(wi);
            let sin_phi_o: Float = sin_phi(wo);
            let cos_phi_o: Float = cos_phi(wo);
            let d_cos: Float = cos_phi_i * cos_phi_o + sin_phi_i * sin_phi_o;
            max_cos = d_cos.max(0.0 as Float);
        }
        // compute sine and tangent terms of Oren-Nayar model
        let sin_alpha: Float;
        let tan_beta: Float;
        if abs_cos_theta(wi) > abs_cos_theta(wo) {
            sin_alpha = sin_theta_o;
            tan_beta = sin_theta_i / abs_cos_theta(wi);
        } else {
            sin_alpha = sin_theta_i;
            tan_beta = sin_theta_o / abs_cos_theta(wo);
        }
        self.r * Spectrum::new(INV_PI * (self.a + self.b * max_cos * sin_alpha * tan_beta))
    }
    fn sample_f(&self,
                wo: Vector3f,
                wi: &mut Vector3f,
                u: Point2f,
                pdf: &mut Float,
                _sampled_type: &mut u8)
                -> Spectrum {
        *wi = cosine_sample_hemisphere(u);
        if wo.z > 0.0 as Float {
            wi.z *= -1.0 as Float;
        }
        *pdf = self.pdf(wo, *wi);
        self.f(wo, *wi)
    }
    fn get_type(&self) -> u8 {
        BxdfType::BsdfDiffuse as u8 | BxdfType::BsdfReflection as u8
    }
}

pub struct MicrofacetReflection {
    pub r: Spectrum,
    pub distribution: Option<TrowbridgeReitzDistribution>, // TODO: MicrofacetDistribution,
    pub fresnel: Arc<Fresnel + Send + Sync>,
}

impl MicrofacetReflection {
    pub fn new(r: Spectrum,
               distribution: Option<TrowbridgeReitzDistribution>,
               fresnel: Arc<Fresnel + Send + Sync>) -> Self {
        MicrofacetReflection {
            r: r,
            distribution: distribution,
            fresnel: fresnel,
        }
    }
}

impl Bxdf for MicrofacetReflection {
    fn f(&self, wo: Vector3f, wi: Vector3f) -> Spectrum {
        let cos_theta_o: Float = abs_cos_theta(wo);
        let cos_theta_i: Float = abs_cos_theta(wi);
        let mut wh: Vector3f = wi + wo;
        // handle degenerate cases for microfacet reflection
        if cos_theta_i == 0.0 || cos_theta_o == 0.0 {
            return Spectrum::new(0.0);
        }
        if wh.x == 0.0 && wh.y == 0.0 && wh.z == 0.0 {
            return Spectrum::new(0.0);
        }
        wh = vec3_normalize(wh);
        let mut dot: Float = vec3_dot_vec3(wi, wh);
        let f: Spectrum = self.fresnel.evaluate(&mut dot);
        if let Some(ref distribution) = self.distribution {
            return self.r * distribution.d(wh) * distribution.g(wo, wi) * f /
                (4.0 as Float * cos_theta_i * cos_theta_o);
        } else {
            panic!("MicrofacetReflection::f() needs self.distribution");
        }
    }
    fn sample_f(&self,
                wo: Vector3f,
                wi: &mut Vector3f,
                u: Point2f,
                pdf: &mut Float,
                _sampled_type: &mut u8)
                -> Spectrum {
        // sample microfacet orientation $\wh$ and reflected direction $\wi$
        if wo.z == 0.0 as Float {
            return Spectrum::default();
        }
        if let Some(ref distribution) = self.distribution {
            let wh: Vector3f = distribution.sample_wh(wo, u);
            *wi = reflect(wo, wh);
            if !vec3_same_hemisphere_vec3(wo, *wi) {
                return Spectrum::default();
            }
            // compute PDF of _wi_ for microfacet reflection
            *pdf = distribution.pdf(wo, wh) / (4.0 * vec3_dot_vec3(wo, wh));
            return self.f(wo, *wi);
        } else {
            panic!("MicrofacetReflection::f() needs self.distribution");
        }
    }
    fn pdf(&self, wo: Vector3f, wi: Vector3f) -> Float {
        if !vec3_same_hemisphere_vec3(wo, wi) {
            return 0.0 as Float;
        }
        let wh: Vector3f = vec3_normalize(wo + wi);
        if let Some(ref distribution) = self.distribution {
            return distribution.pdf(wo, wh) / (4.0 * vec3_dot_vec3(wo, wh));
        } else {
            panic!("MicrofacetReflection::f() needs self.distribution");
        }
    }
    fn get_type(&self) -> u8 {
        BxdfType::BsdfReflection as u8 | BxdfType::BsdfGlossy as u8
    }
}

/// Utility function to calculate cosine via spherical coordinates.
pub fn cos_theta(w: Vector3f) -> Float {
    w.z
}

/// Utility function to calculate the square cosine via spherical
/// coordinates.
pub fn cos_2_theta(w: Vector3f) -> Float {
    w.z * w.z
}

/// Utility function to calculate the absolute value of the cosine via
/// spherical coordinates.
pub fn abs_cos_theta(w: Vector3f) -> Float {
    w.z.abs()
}

/// Utility function to calculate the square sine via spherical
/// coordinates.
pub fn sin_2_theta(w: Vector3f) -> Float {
    (0.0 as Float).max(1.0 as Float - cos_2_theta(w))
}

/// Utility function to calculate sine via spherical coordinates.
pub fn sin_theta(w: Vector3f) -> Float {
    sin_2_theta(w).sqrt()
}

/// Utility function to calculate the tangent via spherical
/// coordinates.
pub fn tan_theta(w: Vector3f) -> Float {
    sin_theta(w) / cos_theta(w)
}

/// Utility function to calculate the square tangent via spherical
/// coordinates.
pub fn tan_2_theta(w: Vector3f) -> Float {
    sin_2_theta(w) / cos_2_theta(w)
}

/// Utility function to calculate cosine via spherical coordinates.
pub fn cos_phi(w: Vector3f) -> Float {
    let sin_theta: Float = sin_theta(w);
    if sin_theta == 0.0 as Float {
        1.0 as Float
    } else {
        clamp_t(w.x / sin_theta, -1.0, 1.0)
    }
}

/// Utility function to calculate sine via spherical coordinates.
pub fn sin_phi(w: Vector3f) -> Float {
    let sin_theta: Float = sin_theta(w);
    if sin_theta == 0.0 as Float {
        0.0 as Float
    } else {
        clamp_t(w.y / sin_theta, -1.0, 1.0)
    }
}

/// Utility function to calculate square cosine via spherical coordinates.
pub fn cos_2_phi(w: Vector3f) -> Float {
    cos_phi(w) * cos_phi(w)
}

/// Utility function to calculate square sine via spherical coordinates.
pub fn sin_2_phi(w: Vector3f) -> Float {
    sin_phi(w) * sin_phi(w)
}

/// Computes the reflection direction given an incident direction and
/// a surface normal.
pub fn reflect(wo: Vector3f, n: Vector3f) -> Vector3f {
    -wo + n * 2.0 as Float * vec3_dot_vec3(wo, n)
}

/// Computes the refraction direction given an incident direction, a
/// surface normal, and the ratio of indices of refraction (incident
/// and transmitted).
pub fn refract(wi: Vector3f, n: Normal3f, eta: Float, wt: &mut Vector3f) -> bool {
    // compute $\cos \theta_\roman{t}$ using Snell's law
    let cos_theta_i: Float = nrm_dot_vec3(n, wi);
    let sin2_theta_i: Float = (0.0 as Float).max(1.0 as Float - cos_theta_i * cos_theta_i);
    let sin2_theta_t: Float = eta * eta * sin2_theta_i;
    // handle total internal reflection for transmission
    if sin2_theta_t >= 1.0 as Float {
        return false;
    }
    let cos_theta_t: Float = (1.0 as Float - sin2_theta_t).sqrt();
    *wt = -wi * eta + Vector3f::from(n) * (eta * cos_theta_i - cos_theta_t);
    true
}

/// Check that two vectors lie on the same side of of the surface.
pub fn vec3_same_hemisphere_vec3(w: Vector3f, wp: Vector3f) -> bool {
    w.z * wp.z > 0.0 as Float
}

// see reflection.cpp

/// Computes the Fresnel reflection formula for dielectric materials
/// and unpolarized light.
pub fn fr_dielectric(cos_theta_i: &mut Float, eta_i: Float, eta_t: Float) -> Float {
    let not_clamped: Float = *cos_theta_i;
    *cos_theta_i = clamp_t(not_clamped, -1.0, 1.0);
    // potentially swap indices of refraction
    let entering: bool = *cos_theta_i > 0.0;
    // use local copies because of potential swap (otherwise eta_i and
    // eta_t would have to be mutable)
    let mut local_eta_i = eta_i;
    let mut local_eta_t = eta_t;
    if !entering {
        std::mem::swap(&mut local_eta_i, &mut local_eta_t);
        *cos_theta_i = (*cos_theta_i).abs();
    }
    // compute _cos_theta_t_ using Snell's law
    let sin_theta_i: Float = (0.0 as Float).max(1.0 as Float - *cos_theta_i * *cos_theta_i).sqrt();
    let sin_theta_t: Float = local_eta_i / local_eta_t * sin_theta_i;
    // handle total internal reflection
    if sin_theta_t >= 1.0 as Float {
        return 1.0 as Float;
    }
    let cos_theta_t: Float = (0.0 as Float).max(1.0 as Float - sin_theta_t * sin_theta_t).sqrt();
    let r_parl: Float = ((local_eta_t * *cos_theta_i) - (local_eta_i * cos_theta_t)) /
                        ((local_eta_t * *cos_theta_i) + (local_eta_i * cos_theta_t));
    let r_perp: Float = ((local_eta_i * *cos_theta_i) - (local_eta_t * cos_theta_t)) /
                        ((local_eta_i * *cos_theta_i) + (local_eta_t * cos_theta_t));
    (r_parl * r_parl + r_perp * r_perp) / 2.0
}

// see microfacet.h

pub trait MicrofacetDistribution {
    fn d(&self, wh: Vector3f) -> Float;
    fn lambda(&self, w: Vector3f) -> Float;
    fn g1(&self, w: Vector3f) -> Float {
        1.0 as Float / (1.0 as Float + self.lambda(w))
    }
    fn g(&self, wo: Vector3f, wi: Vector3f) -> Float {
        1.0 as Float / (1.0 as Float + self.lambda(wo) + self.lambda(wi))
    }
    fn pdf(&self, wo: Vector3f, wh: Vector3f) -> Float {
        if self.get_sample_visible_area() {
            self.d(wh) * self.g1(wo) * vec3_abs_dot_vec3(wo, wh) / abs_cos_theta(wo)
        } else {
            self.d(wh) * abs_cos_theta(wh)
        }
    }
    fn get_sample_visible_area(&self) -> bool;
}

pub struct TrowbridgeReitzDistribution {
    pub alpha_x: Float,
    pub alpha_y: Float,
    // inherited from class MicrofacetDistribution (see microfacet.h)
    pub sample_visible_area: bool,
}

impl TrowbridgeReitzDistribution {
    pub fn new(alpha_x: Float, alpha_y: Float, sample_visible_area: bool) -> Self {
        TrowbridgeReitzDistribution {
            alpha_x: alpha_x,
            alpha_y: alpha_y,
            sample_visible_area: sample_visible_area,
        }
    }
    /// Microfacet distribution function: In comparison to the
    /// Beckmann-Spizzichino model, Trowbridge-Reitz has higher tails - it
    /// falls off to zero more slowly for directions far from the surface
    /// normal.
    pub fn roughness_to_alpha(roughness: Float) -> Float {
        let mut roughness = roughness;
        let limit: Float = 1e-3 as Float;
        if limit > roughness {
            roughness = limit;
        }
        let x: Float = roughness.ln(); // natural (base e) logarithm
        1.62142 + 0.819955 * x + 0.1734 * x * x + 0.0171201 * x * x * x +
        0.000640711 * x * x * x * x
    }
    pub fn sample_wh(&self, wo: Vector3f, u: Point2f) -> Vector3f {
        let mut wh: Vector3f;
        if !self.sample_visible_area {
            let cos_theta;
            let mut phi: Float = (2.0 * PI) * u[1];
            if self.alpha_x == self.alpha_y {
                let tan_theta2: Float = self.alpha_x * self.alpha_x * u[0] / (1.0 - u[0]);
                cos_theta = 1.0 / (1.0 + tan_theta2).sqrt();
            } else {
                phi = (self.alpha_y / self.alpha_x * (2.0 * PI * u[1] + 0.5 * PI).tan()).atan();
                if u[1] > 0.5 {
                    phi += PI;
                }
                let sin_phi: Float = phi.sin();
                let cos_phi: Float = phi.cos();
                let alphax2: Float = self.alpha_x * self.alpha_x;
                let alphay2: Float = self.alpha_y * self.alpha_y;
                let alpha2: Float = 1.0 / (cos_phi * cos_phi / alphax2 + sin_phi * sin_phi / alphay2);
                let tan_theta2: Float = alpha2 * u[0] / (1.0 - u[0]);
                cos_theta = 1.0 / (1.0 + tan_theta2).sqrt();
            }
            let sin_theta: Float = (0.0 as Float).max(1.0 - cos_theta * cos_theta).sqrt();
            wh = spherical_direction(sin_theta, cos_theta, phi);
            if !vec3_same_hemisphere_vec3(wo, wh) {
                wh = -wh;
            }
        } else {
            let flip: bool = wo.z < 0.0;
            if flip {
                wh = trowbridge_reitz_sample(-wo, self.alpha_x, self.alpha_y, u[0], u[1]);
                wh = -wh;
            } else {
                wh = trowbridge_reitz_sample(wo, self.alpha_x, self.alpha_y, u[0], u[1]);
            }
        }
        wh
    }
}

impl MicrofacetDistribution for TrowbridgeReitzDistribution {
    fn d(&self, wh: Vector3f) -> Float {
        let tan_2_theta: Float = tan_2_theta(wh);
        if tan_2_theta.is_infinite() {
            return 0.0 as Float;
        }
        let cos_4_theta: Float = cos_2_theta(wh) * cos_2_theta(wh);
        let e: Float = (cos_2_phi(wh) / (self.alpha_x * self.alpha_x) +
                        sin_2_phi(wh) / (self.alpha_y * self.alpha_y)) *
                       tan_2_theta;
        1.0 as Float /
        (PI * self.alpha_x * self.alpha_y * cos_4_theta * (1.0 as Float + e) * (1.0 as Float + e))
    }
    fn lambda(&self, w: Vector3f) -> Float {
        let abs_tan_theta: Float = tan_theta(w).abs();
        if abs_tan_theta.is_infinite() {
            return 0.0;
        }
        // compute _alpha_ for direction _w_
        let alpha: Float = (cos_2_phi(w) * self.alpha_x * self.alpha_x +
                            sin_2_phi(w) * self.alpha_y * self.alpha_y)
                .sqrt();
        let alpha_2_tan_2_theta: Float = (alpha * abs_tan_theta) * (alpha * abs_tan_theta);
        (-1.0 as Float + (1.0 as Float + alpha_2_tan_2_theta).sqrt()) / 2.0 as Float
    }
    fn get_sample_visible_area(&self) -> bool {
        self.sample_visible_area
    }
}

fn trowbridge_reitz_sample_11(cos_theta: Float,
                              u1: Float,
                              u2: Float,
                              slope_x: &mut Float,
                              slope_y: &mut Float) {
    // special case (normal incidence)
    if cos_theta > 0.9999 {
        let r: Float = (u1 / (1.0 - u1)).sqrt();
        let phi: Float = 6.28318530718 * u2;
        *slope_x = r * phi.cos();
        *slope_y = r * phi.sin();
        return;
    }

    let sin_theta: Float = (0.0 as Float).max(1.0 as Float - cos_theta * cos_theta).sqrt();
    let tan_theta: Float = sin_theta / cos_theta;
    let a: Float = 1.0 / tan_theta;
    let g1: Float = 2.0 / (1.0 + (1.0 + 1.0 / (a * a)).sqrt());

    // sample slope_x
    let a: Float = 2.0 * u1 / g1 - 1.0;
    let mut tmp: Float = 1.0 / (a * a - 1.0);
    if tmp > 1e10
    {
        tmp = 1e10;
    }
    let b: Float = tan_theta;
    let d: Float = (b * b * tmp * tmp - (a * a - b * b) * tmp).max(0.0 as Float).sqrt();
    let slope_x_1: Float = b * tmp - d;
    let slope_x_2: Float = b * tmp + d;
    if a < 0.0 || slope_x_2 > 1.0 / tan_theta {
        *slope_x = slope_x_1;
    } else {
        *slope_x = slope_x_2;
    }

    // sample slope_y
    let s: Float;
    let new_u2: Float;
    if u2 > 0.5 {
        s = 1.0;
        new_u2 = 2.0 * (u2 - 0.5);
    } else {
        s = -1.0;
        new_u2 = 2.0 * (0.5 - u2);
    }
    let z: Float =
        (new_u2 * (new_u2 * (new_u2 * 0.27385 - 0.73369) + 0.46341)) /
        (new_u2 * (new_u2 * (new_u2 * 0.093073 + 0.309420) - 1.0) + 0.597999);
    *slope_y = s * z * (1.0 + *slope_x * *slope_x).sqrt();

    assert!(!(*slope_y).is_infinite());
    assert!(!(*slope_y).is_nan());
}

fn trowbridge_reitz_sample(wi: Vector3f,
                           alpha_x: Float,
                           alpha_y: Float,
                           u1: Float,
                           u2: Float) -> Vector3f {
    // 1. stretch wi
    let wi_stretched: Vector3f = vec3_normalize(Vector3f {
        x: alpha_x * wi.x,
        y: alpha_y * wi.y,
        z: wi.z,
    });

    // 2. simulate P22_{wi}(x_slope, y_slope, 1, 1)
    let mut slope_x: Float = 0.0;
    let mut slope_y: Float = 0.0;
    trowbridge_reitz_sample_11(cos_theta(wi_stretched), u1, u2, &mut slope_x, &mut slope_y);

    // 3. rotate
    let tmp: Float = cos_phi(wi_stretched) * slope_x - sin_phi(wi_stretched) * slope_y;
    slope_y = sin_phi(wi_stretched) * slope_x + cos_phi(wi_stretched) * slope_y;
    slope_x = tmp;

    // 4. unstretch
    slope_x = alpha_x * slope_x;
    slope_y = alpha_y * slope_y;

    // 5. compute normal
    vec3_normalize(Vector3f {
        x: -slope_x,
        y: -slope_y,
        z:1.0,
    })
}

// see material.h

/// **Material** defines the interface that material implementations
/// must provide.
pub trait Material {
    /// The method is given a **SurfaceInteraction** object that
    /// contains geometric properties at an intersection point on the
    /// surface of a shape and is responsible for determining the
    /// reflective properties at the point and initializing some
    /// member variables.
    fn compute_scattering_functions(&self,
                                    si: &mut SurfaceInteraction,
                                    // arena: &mut Arena,
                                    mode: TransportMode,
                                    allow_multiple_lobes: bool);
}

// see matte.h

/// Describes a purely diffuse surface.
pub struct MatteMaterial {
    pub kd: Arc<Texture<Spectrum> + Sync + Send>, // default: 0.5
    pub sigma: Float, // default: 0.0
    // TODO: bump_map
}

impl MatteMaterial {
    pub fn new(kd: Arc<Texture<Spectrum> + Send + Sync>, sigma: Float) -> Self {
        MatteMaterial {
            kd: kd,
            sigma: sigma,
        }
    }
    pub fn bsdf(&self, si: &SurfaceInteraction) -> Bsdf {
        let mut bxdfs: Vec<Box<Bxdf + Send + Sync>> = Vec::new();
        let r: Spectrum = self.kd.evaluate(si).clamp(0.0 as Float, std::f32::INFINITY as Float);
        if !r.is_black() {
            if self.sigma == 0.0 {
                bxdfs.push(Box::new(LambertianReflection::new(r)));
            } else {
                bxdfs.push(Box::new(OrenNayar::new(r, self.sigma)));
            }
        }
        Bsdf::new(si, 1.5, bxdfs)
    }
}

impl Material for MatteMaterial {
    fn compute_scattering_functions(&self,
                                    si: &mut SurfaceInteraction,
                                    // arena: &mut Arena,
                                    _mode: TransportMode,
                                    _allow_multiple_lobes: bool) {
        si.bsdf = Some(Arc::new(self.bsdf(si)));
    }
}

// see plastic.h

/// Plastic can be modeled as a mixture of a diffuse and glossy
/// scattering function.
pub struct PlasticMaterial {
    pub kd: Arc<Texture<Spectrum> + Sync + Send>, // default: 0.25
    pub ks: Arc<Texture<Spectrum> + Sync + Send>, // default: 0.25
    pub roughness: Arc<Texture<Float> + Sync + Send>, // default: 0.1
    // TODO: bump_map
    pub remap_roughness: bool,
}

impl PlasticMaterial {
    pub fn new(kd: Arc<Texture<Spectrum> + Send + Sync>,
               ks: Arc<Texture<Spectrum> + Send + Sync>,
               roughness: Arc<Texture<Float> + Sync + Send>,
               remap_roughness: bool) -> Self {
        PlasticMaterial {
            kd: kd,
            ks: ks,
            roughness: roughness,
            remap_roughness: remap_roughness,
        }
    }
    pub fn bsdf(&self, si: &SurfaceInteraction) -> Bsdf {
        let mut bxdfs: Vec<Box<Bxdf + Send + Sync>> = Vec::new();
        // initialize diffuse component of plastic material
        let kd: Spectrum = self.kd.evaluate(si).clamp(0.0 as Float, std::f32::INFINITY as Float);
        if !kd.is_black() {
            bxdfs.push(Box::new(LambertianReflection::new(kd)));
        }
        // initialize specular component of plastic material
        let ks: Spectrum = self.ks.evaluate(si).clamp(0.0 as Float, std::f32::INFINITY as Float);
        if !ks.is_black() {
            let fresnel = Arc::new(FresnelDielectric {
                eta_i: 1.5 as Float,
                eta_t: 1.0 as Float,
            });
            // create microfacet distribution _distrib_ for plastic material
            let mut rough: Float = self.roughness.evaluate(si);
            if self.remap_roughness {
                rough = TrowbridgeReitzDistribution::roughness_to_alpha(rough);
            }
            let distrib: TrowbridgeReitzDistribution = TrowbridgeReitzDistribution {
                alpha_x: rough,
                alpha_y: rough,
                sample_visible_area: true,
            };
            bxdfs.push(Box::new(MicrofacetReflection::new(ks, Some(distrib), fresnel)));
        }
        Bsdf::new(si, 1.0, bxdfs)
    }
}

impl Material for PlasticMaterial {
    fn compute_scattering_functions(&self,
                                    si: &mut SurfaceInteraction,
                                    // arena: &mut Arena,
                                    _mode: TransportMode,
                                    _allow_multiple_lobes: bool) {
        si.bsdf = Some(Arc::new(self.bsdf(si)));
    }
}

// see glass.h

/// Perfect or glossy specular reflection and transmission, weighted
/// by Fresnel terms for accurate angular-dependent variation.
pub struct GlassMaterial {
    pub kr: Arc<Texture<Spectrum> + Sync + Send>, // default: 1.0
    pub kt: Arc<Texture<Spectrum> + Sync + Send>, // default: 1.0
    pub u_roughness: Arc<Texture<Float> + Sync + Send>, // default: 0.0
    pub v_roughness: Arc<Texture<Float> + Sync + Send>, // default: 0.0
    pub index: Arc<Texture<Float> + Sync + Send>, // TODO: bump_map
    pub remap_roughness: bool,
}

impl GlassMaterial {
    pub fn bsdf(&self, si: &SurfaceInteraction, mode: TransportMode, allow_multiple_lobes: bool) -> Bsdf {
        let mut bxdfs: Vec<Box<Bxdf + Send + Sync>> = Vec::new();
        let eta: Float = self.index.evaluate(si);
        let mut urough: Float = self.u_roughness.evaluate(si);
        let mut vrough: Float = self.v_roughness.evaluate(si);
        let r: Spectrum = self.kr.evaluate(si).clamp(0.0 as Float, std::f32::INFINITY as Float);
        let t: Spectrum = self.kt.evaluate(si).clamp(0.0 as Float, std::f32::INFINITY as Float);
        let is_specular: bool = urough == 0.0 as Float && vrough == 0.0 as Float;
        if is_specular && allow_multiple_lobes {
            bxdfs.push(Box::new(FresnelSpecular::new(r, t, 1.0 as Float, eta, mode)));
        } else {
            if self.remap_roughness {
                urough = TrowbridgeReitzDistribution::roughness_to_alpha(urough);
                vrough = TrowbridgeReitzDistribution::roughness_to_alpha(vrough);
            }
            let distrib: Option<TrowbridgeReitzDistribution> = match is_specular {
                true => None,
                false => Some(TrowbridgeReitzDistribution::new(urough, vrough, true)),
            };
            if !r.is_black() {
                let fresnel = Arc::new(FresnelDielectric {
                    eta_i: 1.0 as Float,
                    eta_t: eta,
                });
                if is_specular {
                    bxdfs.push(Box::new(SpecularReflection::new(r, fresnel)));
                } else {
                    bxdfs.push(Box::new(MicrofacetReflection::new(r, distrib, fresnel)));
                }
            }
            if !t.is_black() {
                if is_specular {
                    bxdfs.push(Box::new(SpecularTransmission::new(t, 1.0, eta, mode)));
                } else {
                    // TODO: si->bsdf->Add(ARENA_ALLOC(arena, MicrofacetTransmission)(
                    // T, distrib, 1.f, eta, mode));
                }
            }
        }
        Bsdf::new(si, eta, bxdfs)
    }
}

impl Material for GlassMaterial {
    fn compute_scattering_functions(&self,
                                    si: &mut SurfaceInteraction,
                                    // arena: &mut Arena,
                                    mode: TransportMode,
                                    allow_multiple_lobes: bool) {
        si.bsdf = Some(Arc::new(self.bsdf(si, mode, allow_multiple_lobes)));
    }
}

// see mirror.h

/// A simple mirror, modeled with perfect specular reflection.
pub struct MirrorMaterial {
    pub kr: Arc<Texture<Spectrum> + Sync + Send>, // default: 0.9
    // TODO: bump_map
}

impl MirrorMaterial {
    pub fn new(kr: Arc<Texture<Spectrum> + Send + Sync>) -> Self {
        MirrorMaterial {
            kr: kr,
        }
    }
    pub fn bsdf(&self, si: &SurfaceInteraction) -> Bsdf {
        let mut bxdfs: Vec<Box<Bxdf + Send + Sync>> = Vec::new();
        let r: Spectrum = self.kr.evaluate(si).clamp(0.0 as Float, std::f32::INFINITY as Float);
        let fresnel = Arc::new(FresnelNoOp{});
        bxdfs.push(Box::new(SpecularReflection::new(r, fresnel)));
        Bsdf::new(si, 1.5, bxdfs)
    }
}

impl Material for MirrorMaterial {
    fn compute_scattering_functions(&self,
                                    si: &mut SurfaceInteraction,
                                    // arena: &mut Arena,
                                    _mode: TransportMode,
                                    _allow_multiple_lobes: bool) {
        si.bsdf = Some(Arc::new(self.bsdf(si)));
    }
}

// see texture.h

pub trait TextureMapping2D {
    fn map(&self, si: &SurfaceInteraction, dstdx: &mut Vector2f, dstdy: &mut Vector2f) -> Point2f;
}

#[derive(Debug,Default,Copy,Clone)]
pub struct UVMapping2D {
    pub su: Float,
    pub sv: Float,
    pub du: Float,
    pub dv: Float,
}

impl TextureMapping2D for UVMapping2D {
    fn map(&self, si: &SurfaceInteraction, dstdx: &mut Vector2f, dstdy: &mut Vector2f) -> Point2f {
        // compute texture differentials for 2D identity mapping
        *dstdx = Vector2f {
            x: si.dudx * self.su,
            y: si.dvdx * self.sv,
        };
        *dstdy = Vector2f {
            x: si.dudy * self.su,
            y: si.dvdy * self.sv,
        };
        Point2f {
            x: si.uv[0] * self.su + self.du,
            y: si.uv[1] * self.sv + self.dv
        }
    }
}

#[derive(Debug,Default,Copy,Clone)]
pub struct PlanarMapping2D {
    pub vs: Vector3f,
    pub vt: Vector3f,
    pub ds: Float,
    pub dt: Float,
}

impl TextureMapping2D for PlanarMapping2D {
    fn map(&self, si: &SurfaceInteraction, dstdx: &mut Vector2f, dstdy: &mut Vector2f) -> Point2f {
        let vec: Vector3f = Vector3f {
            x: si.p.x,
            y: si.p.y,
            z: si.p.z,
        };
        *dstdx = Vector2f {
            x: vec3_dot_vec3(si.dpdx, self.vs),
            y: vec3_dot_vec3(si.dpdx, self.vt),
        };
        *dstdy = Vector2f {
            x: vec3_dot_vec3(si.dpdy, self.vs),
            y: vec3_dot_vec3(si.dpdy, self.vt),
        };
        Point2f {
            x: self.ds + vec3_dot_vec3(vec, self.vs),
            y: self.dt + vec3_dot_vec3(vec, self.vt),
        }
    }
}

pub trait Texture<T> {
    fn evaluate(&self, si: &SurfaceInteraction) -> T;
}

pub fn lanczos(x: Float, tau: Float) -> Float {
    let mut x: Float = x;
    x = x.abs();
    if x < 1e-5 as Float {
        return 1.0 as Float;
    }
    if x > 1.0 as Float {
        return 0.0 as Float;
    }
    x *= PI;
    let s: Float = (x * tau).sin() / (x * tau);
    let lanczos: Float = x.sin() / x;
    s * lanczos
}

// see constant.h

pub struct ConstantTexture<T> {
    pub value: T,
}

impl<T: Copy> ConstantTexture<T> {
    pub fn new(value: T) -> Self {
        ConstantTexture { value: value }
    }
}

impl<T: Copy> Texture<T> for ConstantTexture<T> {
    fn evaluate(&self, _si: &SurfaceInteraction) -> T {
        self.value
    }
}

// checkerboard.h

pub struct Checkerboard2DTexture<T> {
    pub tex1: Arc<Texture<T> + Send + Sync>,
    pub tex2: Arc<Texture<T> + Send + Sync>,
    pub mapping: Box<TextureMapping2D + Send + Sync>,
    // TODO: const AAMethod aaMethod;
}

impl<T: Copy> Checkerboard2DTexture<T> {
    pub fn new(mapping: Box<TextureMapping2D + Send + Sync>,
               tex1: Arc<Texture<T> + Send + Sync>,
               tex2: Arc<Texture<T> + Send + Sync>// , TODO: aaMethod
    ) -> Self {
        Checkerboard2DTexture {
            tex1: tex1,
            tex2: tex2,
            mapping: mapping,
        }
    }
}

impl<T: Copy> Texture<T> for Checkerboard2DTexture<T> {
    fn evaluate(&self, si: &SurfaceInteraction) -> T {
        let mut dstdx: Vector2f = Vector2f::default();
        let mut dstdy: Vector2f = Vector2f::default();
        let st: Point2f = self.mapping.map(si, &mut dstdx, &mut dstdy);
        // TODO: if (aaMethod == AAMethod::None) {
        if (st.x.floor() as u32 + st.y.floor() as u32) % 2 == 0 {
            self.tex1.evaluate(si)
        } else {
            self.tex2.evaluate(si)
        }
    }
}

// see mipmap.h

const WEIGHT_LUT_SIZE: usize = 128;

#[derive(Debug,Clone)]
pub enum ImageWrap {
    Repeat,
    Black,
    Clamp,
}

#[derive(Debug,Default,Copy,Clone)]
pub struct ResampleWeight {
    pub first_texel: i32,
    pub weight: [Float; 4],
}

pub struct MipMap {
    // MIPMap Private Data
    pub do_trilinear: bool,
    pub max_anisotropy: Float,
    pub wrap_mode: ImageWrap,
    pub resolution: Point2i,
    pub pyramid: Vec<BlockedArray<Spectrum>>,
    // TODO: static Float weightLut[WeightLUTSize];
    pub weight_lut: [Float; WEIGHT_LUT_SIZE],
}

impl MipMap {
    pub fn new(res: &Point2i,
               img: &[Spectrum],
               do_trilinear: bool,
               max_anisotropy: Float,
               wrap_mode: ImageWrap) -> Self {
        let mut resolution = *res;
        let mut resampled_image: Vec<Spectrum> = Vec::new();
        if !is_power_of_2(resolution.x) || !is_power_of_2(resolution.y) {
            // resample image to power-of-two resolution
            let res_pow_2: Point2i = Point2i {
                x: round_up_pow2_32(resolution.x),
                y: round_up_pow2_32(resolution.y),
            };
            println!("Resampling MIPMap from {:?} to {:?}. Ratio= {:?}",
                     resolution, res_pow_2,
                     (res_pow_2.x * res_pow_2.y) as Float /
                     (resolution.x * resolution.y) as Float);
            // resample image in $s$ direction
            let s_weights: Vec<ResampleWeight> = MipMap::resample_weights(resolution.x, res_pow_2.x);
            // TODO: resampled_image.reset(new T[resPow2[0] * resPow2[1]]);
            resampled_image = vec![Spectrum::default(); (res_pow_2.x * res_pow_2.y) as usize];
            // apply _s_weights_ to zoom in $s$ direction
            // TODO: ParallelFor([&](int t) {
            for t in 0..resolution.y { // chunk size 16
                for s in 0..res_pow_2.x {
                    // compute texel $(s,t)$ in $s$-zoomed image
                    resampled_image[(t * res_pow_2.x + s) as usize] = Spectrum::new(0.0 as Float);
                    for j in 0..4 {
                        let mut orig_s: i32 = s_weights[s as usize].first_texel + j as i32;
                        orig_s = match wrap_mode {
                            ImageWrap::Repeat => mod_t(orig_s, resolution.x),
                            ImageWrap::Clamp => clamp_t(orig_s, 0_i32, resolution.x - 1_i32),
                            _ => orig_s,
                        };
                        if orig_s >= 0_i32 && orig_s < resolution.x {
                            resampled_image[(t * res_pow_2.x + s) as usize] +=
                                img[(t * resolution.x + orig_s) as usize] *
                                s_weights[s as usize].weight[j];

                        }
                    }
                }
            }
            // TODO: }, resolution[1], 16);
            // resample image in $t$ direction
            let t_weights: Vec<ResampleWeight> = MipMap::resample_weights(resolution.y, res_pow_2.y);
            // std::vector<T *> resample_bufs;
            // int nThreads = MaxThreadIndex();
            // for (int i = 0; i < nThreads; ++i)
            //     resample_bufs.push_back(new T[resPow2[1]]);
            // let resampled_bufs: Vec<Spectrum> = vec![Spectrum::default(); res_pow_2.y as usize]; // single-threaded
            let mut work_data: Vec<Spectrum> = vec![Spectrum::default(); res_pow_2.y as usize]; // single-threaded
            // TODO: ParallelFor([&](int s) {
            // T *work_data = resample_bufs[ThreadIndex];
            for s in 0..res_pow_2.x { // chunk size 32
                for t in 0..res_pow_2.y {
                    work_data[t as usize] = Spectrum::new(0.0 as Float);
                    for j in 0..4 {
                        let mut offset: i32 = t_weights[t as usize].first_texel + j as i32;
                        offset = match wrap_mode {
                            ImageWrap::Repeat => mod_t(offset, resolution.y),
                            ImageWrap::Clamp => clamp_t(offset, 0_i32, resolution.y - 1_i32),
                            _ => offset,
                        };
                        if offset >= 0_i32 && offset < resolution.y {
                            work_data[t as usize] += 
                                resampled_image[(offset * res_pow_2.x + s) as usize] *
                                t_weights[t as usize].weight[j];
                        }
                    }
                }
                for t in 0..res_pow_2.y {
                    resampled_image[(t * res_pow_2.x + s) as usize] = MipMap::clamp(work_data[t as usize]);
                }
            }
            // TODO: }, resPow2[0], 32);
            // for (auto ptr : resample_bufs) delete[] ptr;
            resolution = res_pow_2;
        }
        let mut mipmap = MipMap {
            do_trilinear: do_trilinear,
            max_anisotropy: max_anisotropy,
            wrap_mode: wrap_mode,
            resolution: resolution,
            pyramid: Vec::new(),
            weight_lut: [0.0 as Float; WEIGHT_LUT_SIZE],
        };
        // initialize levels of MipMap for image
        let n_levels = 1 + (std::cmp::max(resolution.x, resolution.y) as Float).log2() as usize;
        // initialize most detailed level of MipMap
        let img_data: &[Spectrum] = if resampled_image.is_empty() {
            img
        } else {
            &resampled_image[..]
        };
        mipmap.pyramid.push(BlockedArray::new_from(resolution.x as usize,
                                                   resolution.y as usize,
                                                   img_data));
        for i in 1..n_levels {
            // initialize $i$th MipMap level from $i-1$st level
            let s_res = std::cmp::max(1, mipmap.pyramid[i - 1].u_size() / 2);
            let t_res = std::cmp::max(1, mipmap.pyramid[i - 1].v_size() / 2);
            let mut ba = BlockedArray::new(s_res, t_res);
            // filter 4 texels from finer level of pyramid
            for t in 0..t_res {
                for s in 0..s_res {
                    let (si, ti) = (s as isize, t as isize);
                    ba[(s, t)] = (*mipmap.texel(i - 1, 2 * si, 2 * ti) +
                                  *mipmap.texel(i - 1, 2 * si + 1, 2 * ti) +
                                  *mipmap.texel(i - 1, 2 * si, 2 * ti + 1) +
                                  *mipmap.texel(i - 1, 2 * si + 1, 2 * ti + 1)) * 0.25;
                }
            }
            mipmap.pyramid.push(ba);
        }
        // initialize EWA filter weights if needed
        if mipmap.weight_lut[0] == 0.0 as Float {
            for i in 0..WEIGHT_LUT_SIZE {
                let alpha: Float = 2.0 as Float;
                let r2: Float = i as Float / (WEIGHT_LUT_SIZE - 1) as Float;
                mipmap.weight_lut[i] = (-alpha * r2).exp() - (-alpha).exp();
            }
        }
        // TODO: mipMapMemory += (4 * resolution[0] * resolution[1] * sizeof(T)) / 3;
        mipmap
    }
    pub fn width(&self) -> i32 {
        self.resolution.x
    }
    pub fn height(&self) -> i32 {
        self.resolution.y
    }
    pub fn levels(&self) -> usize {
        self.pyramid.len()
    }
    pub fn texel(&self, level: usize, s: isize, t: isize) -> &Spectrum {
        let l = &self.pyramid[level];
        let (u_size, v_size) = (l.u_size() as isize, l.v_size() as isize);
        let (ss, tt): (usize, usize) = match self.wrap_mode {
            ImageWrap::Repeat => (mod_t(s as usize, u_size as usize),
                                  mod_t(t as usize, v_size as usize)),
            ImageWrap::Clamp => {
                (clamp_t(s, 0, u_size - 1) as usize, clamp_t(t, 0, v_size - 1) as usize)
            }
            ImageWrap::Black => {
                // TODO: let black: T = num::Zero::zero();
                if s < 0 || s >= u_size || t < 0 || t >= v_size {
                    // TODO: return &black;
                    (clamp_t(s, 0, u_size - 1) as usize, clamp_t(t, 0, v_size - 1) as usize) // TMP
                } else {
                    (s as usize, t as usize)
                }
            }
        };
        &l[(ss, tt)]
    }
    pub fn lookup_pnt_flt(&self, st: &Point2f, width: Float) -> Spectrum {
        // TODO: ++nTrilerpLookups;
        // TODO: ProfilePhase p(Prof::TexFiltTrilerp);
        // compute MIPMap level for trilinear filtering
        let level: Float = self.levels() as Float - 1.0 as Float + width.max(1e-8 as Float).log2();
        // perform trilinear interpolation at appropriate MIPMap level
        if level < 0.0 as Float {
            return self.triangle(0_usize, st);
        } else if level >= self.levels() as Float - 1 as Float {
            return *self.texel(self.levels(), 0_isize, 0_isize);
        } else {
            let i_level: usize = level.floor() as usize;
            let delta: Float = level - i_level as Float;
            return lerp_rgb(delta, self.triangle(i_level, st), self.triangle(i_level + 1_usize, st));
        }
    }
    pub fn lookup_pnt_vec_vec(&self, st: &Point2f, dst0: &mut Vector2f, dst1: &mut Vector2f) -> Spectrum {
        if self.do_trilinear {
            let width: Float = dst0.x.abs().max(dst0.y.abs()).max(dst1.x.abs().max(dst1.y.abs()));
            println!("TODO: Lookup(st, 2 * width) = Lookup({:?}, {:?});", st, 2.0 as Float * width)
        }
        // TODO: ++nEWALookups;
        // TODO: ProfilePhase p(Prof::TexFiltEWA);
        // compute ellipse minor and major axes
        if dst0.length_squared() < dst1.length_squared() {
            // std::swap(dst0, dst1);
            let swap: Vector2f = Vector2f { x: dst0.x, y: dst0.y, };
            // dst0 = dst1
            dst0.x = dst1.x;
            dst0.y = dst1.y;
            // dst1 = dst0
            dst1.x = swap.x;
            dst1.y = swap.y;
        }
        let major_length: Float = dst0.length();
        let mut minor_length: Float = dst1.length();
        // clamp ellipse eccentricity if too large
        if minor_length * self.max_anisotropy < major_length && minor_length > 0.0 as Float {
            let scale: Float = major_length / (minor_length * self.max_anisotropy);
            *dst1 *= scale;
            minor_length *= scale;
        }
        if minor_length == 0.0 as Float
        {
            return self.triangle(0, st);
        }
        // choose level of detail for EWA lookup and perform EWA filtering
        let lod: Float = (0.0 as Float).max(self.levels() as Float - 1.0 as Float + minor_length.log2() as Float);
        let ilod: usize = lod.floor() as usize;
        let col2: Spectrum = self.ewa(ilod + 1, st.clone(), dst0.clone(), dst1.clone());
        let col1: Spectrum = self.ewa(ilod, st.clone(), dst0.clone(), dst1.clone());
        let ret: Spectrum = lerp_rgb(lod - ilod as Float, col1, col2);
        ret
    }
    fn resample_weights(old_res: i32, new_res: i32) -> Vec<ResampleWeight> {
        assert!(new_res >= old_res);
        let mut wt: Vec<ResampleWeight> = Vec::with_capacity(new_res as usize);
        let filterwidth: Float = 2.0 as Float;
        for i in 0..new_res {
            // compute image resampling weights for _i_th texel
            let center: Float = (i as Float + 0.5 as Float) * old_res as Float / new_res as Float;
            let mut rw: ResampleWeight = ResampleWeight::default();
            rw.first_texel = ((center - filterwidth) + 0.5 as Float).floor() as i32;
            for j in 0..4 {
                let pos: Float = rw.first_texel as Float + j as Float + 0.5 as Float;
                rw.weight[j] = lanczos((pos - center) / filterwidth, 2.0 as Float);
            }
            // normalize filter weights for texel resampling
            let inv_sum_wts: Float =
                1.0 as Float / (rw.weight[0] + rw.weight[1] + rw.weight[2] + rw.weight[3]);
            for j in 0..4 {
                rw.weight[j] *= inv_sum_wts;
            }
            wt.push(rw); // add to vector
        }
        wt
    }
    fn clamp(v: Spectrum) -> Spectrum {
        v.clamp(0.0 as Float, std::f32::INFINITY as Float)
    }
    fn triangle(&self, level: usize, st: &Point2f) -> Spectrum {
        let level: usize = clamp_t(level, 0_usize, self.levels() - 1_usize);
        let s: Float = st.x * self.pyramid[level].u_size() as Float - 0.5;
        let t: Float = st.y * self.pyramid[level].v_size() as Float - 0.5;
        let s0: isize = s.floor() as isize;
        let t0: isize = t.floor() as isize;
        let ds: Float = s - s0 as Float;
        let dt: Float = t - t0 as Float;
        *self.texel(level, s0, t0) * (1.0 - ds) * (1.0 - dt) +
        *self.texel(level, s0, t0 + 1) * (1.0 - ds) * dt +
        *self.texel(level, s0 + 1, t0) * ds * (1.0 - dt) +
        *self.texel(level, s0 + 1, t0 + 1) * ds * dt
    }
    fn ewa(&self, level: usize, st: Point2f, dst0: Vector2f, dst1: Vector2f) -> Spectrum {
        if level >= self.levels() {
            return *self.texel(self.levels() - 1, 0, 0);
        }
        // convert EWA coordinates to appropriate scale for level
        let mut new_st: Vector2f = Vector2f { x: st.x, y: st.y, };
        new_st.x = new_st.x * self.pyramid[level].u_size() as Float - 0.5 as Float;
        new_st.y = new_st.y * self.pyramid[level].v_size() as Float - 0.5 as Float;
        let mut new_dst0: Vector2f = Vector2f { x: dst0.x, y: dst0.y, };
        let mut new_dst1: Vector2f = Vector2f { x: dst1.x, y: dst1.y, };
        new_dst0.x *= self.pyramid[level].u_size() as Float;
        new_dst0.y *= self.pyramid[level].v_size() as Float;
        new_dst1.x *= self.pyramid[level].u_size() as Float;
        new_dst1.y *= self.pyramid[level].v_size() as Float;
        // compute ellipse coefficients to bound EWA filter region
        let mut a: Float = new_dst0.y * new_dst0.y + new_dst1.y * new_dst1.y + 1.0 as Float;
        let mut b: Float = -2.0 as Float * (new_dst0.x * new_dst0.y + new_dst1.x * new_dst1.y);
        let mut c: Float = new_dst0.x * new_dst0.x + new_dst1.x * new_dst1.x + 1.0 as Float;
        let inv_f: Float = 1.0 as Float / (a * c - b * b * 0.25 as Float);
        a *= inv_f;
        b *= inv_f;
        c *= inv_f;
        // compute the ellipse's $(s,t)$ bounding box in texture space
        let det: Float = -b * b + 4.0 as Float * a * c;
        let inv_det: Float = 1.0 as Float / det;
        let u_sqrt: Float = (det * c).sqrt();
        let v_sqrt: Float = (a * det).sqrt();
        let s0: isize = (new_st.x - 2.0 as Float * inv_det * u_sqrt).ceil() as isize;
        let s1: isize  = (new_st.x + 2.0 as Float * inv_det * u_sqrt).floor() as isize;
        let t0: isize  = (new_st.y - 2.0 as Float * inv_det * v_sqrt).ceil() as isize;
        let t1: isize  = (new_st.y + 2.0 as Float * inv_det * v_sqrt).floor() as isize;
        // scan over ellipse bound and compute quadratic equation
        let mut sum: Spectrum = Spectrum::new(0.0 as Float);
        let mut sum_wts: Float = 0.0;
        for it in t0..(t1 + 1) {
            let tt: Float = it as Float - new_st.y;
            for is in s0..(s1 + 1) {
                let ss: Float = is as Float - new_st.x;
                // compute squared radius and filter texel if inside ellipse
                let r2: Float = a * ss * ss + b * ss * tt + c * tt * tt;
                if r2 < 1.0 as Float {
                    let index: usize = std::cmp::min((r2 * WEIGHT_LUT_SIZE as Float) as usize, WEIGHT_LUT_SIZE - 1);
                    let weight: Float = self.weight_lut[index];
                    sum += *self.texel(level, is as isize, it as isize) * weight;
                    sum_wts += weight;
                }
            }
        }
        sum / sum_wts
    }
}

// see imagemap.h

pub struct ImageTexture {
    pub mapping: Box<TextureMapping2D + Send + Sync>,
    pub mipmap: Arc<MipMap>,
}

impl ImageTexture {
    pub fn new(mapping: Box<TextureMapping2D + Send + Sync>,
               filename: String,
               do_trilinear: bool,
               max_aniso: Float,
               wrap_mode: ImageWrap,
               _scale: Float,
               _gamma: bool)
               -> Self {
        let path = Path::new(&filename);
        let img_result: ImageResult<DynamicImage> = image::open(path);
        if !img_result.is_ok() {
            panic!("Error reading \"{}\"", filename);
        }
        let buf = img_result.unwrap();
        let rgb = buf.to_rgb();
        let res = Point2i {
            x: rgb.width() as i32,
            y: rgb.height() as i32,
        };
        // instead of convertIn(texels[i], &convertedTexels[i], scale, gamma);
        let mut texels: Vec<Spectrum> = rgb.pixels()
            .map(|p| Spectrum::from_srgb(&p.data))
            .collect();
        // flip image in y; texture coordinate space has (0,0) at the
        // lower left corner.
        for y in 0..res.y / 2 {
            for x in 0..res.x {
                let o1 = (y * res.x + x) as usize;
                let o2 = ((res.y - 1 - y) * res.x + x) as usize;
                texels.swap(o1, o2);
            }
        }
        // create _MipMap_ from converted texels (see above)
        let mipmap = Arc::new(MipMap::new(&res, &texels[..], do_trilinear, max_aniso, wrap_mode));
        ImageTexture {
            mapping: mapping,
            mipmap: mipmap,
        }
    }
    pub fn convert_out(from: &Spectrum, to: &mut Spectrum) {
        let mut rgb: [Float; 3] = [0.0 as Float; 3];
        from.to_rgb(&mut rgb);
        *to = Spectrum::from_rgb(&rgb);
    }
}

impl Texture<Spectrum> for ImageTexture {
    fn evaluate(&self, si: &SurfaceInteraction) -> Spectrum {
        // Vector2f dstdx, dstdy;
        // Point2f st = mapping->Map(si, &dstdx, &dstdy);
        // Tmemory mem = mipmap->Lookup(st, dstdx, dstdy);
        // Treturn ret;
        // convertOut(mem, &ret);
        // return ret;
        let mut dstdx: Vector2f = Vector2f::default();
        let mut dstdy: Vector2f = Vector2f::default();
        let st: Point2f = self.mapping.map(si, &mut dstdx, &mut dstdy);
        let mem: Spectrum = self.mipmap.lookup_pnt_vec_vec(&st, &mut dstdx, &mut dstdy);
        let mut ret: Spectrum = Spectrum::new(0.0);
        ImageTexture::convert_out(&mem, &mut ret);
        ret
    }
}

// see light.h

#[repr(u8)]
pub enum LightFlags {
    DeltaPosition = 1,
    DeltaDirection = 2,
    Area = 4,
    Infinite = 8,
}

pub trait Light {
    /// Returns the radiance arriving at a point at a certain time due
    /// to the light, assuming there are no occluding objects between
    /// them.
    fn sample_li(&self,
                 iref: &InteractionCommon,
                 u: Point2f,
                 wi: &mut Vector3f,
                 pdf: &mut Float,
                 vis: &mut VisibilityTester)
                 -> Spectrum;
    fn power(&self) -> Spectrum;
    fn preprocess(&self, scene: &Scene);
    fn le(&self, _ray: &mut Ray) -> Spectrum;
    fn pdf_li(&self, iref: &Interaction, wi: Vector3f) -> Float;
    fn get_flags(&self) -> u8;
    fn get_n_samples(&self) -> i32;
}

/// Check if LightFlags::DeltaPosition or LightFlags::DeltaDirection
/// is set.
pub fn is_delta_light(flags: u8) -> bool {
    let mut pos: bool = false;
    let mut dir: bool = false;
    if (flags & LightFlags::DeltaPosition as u8) > 0 {
        pos = true;
    }
    if (flags & LightFlags::DeltaDirection as u8) > 0 {
        dir = true;
    }
    pos || dir
}

/// A closure - an object that encapsulates a small amount of data and
/// some computation that is yet to be done.
#[derive(Debug,Default,Copy,Clone)]
pub struct VisibilityTester {
    pub p0: InteractionCommon, // TODO: private
    pub p1: InteractionCommon, // TODO: private
}

impl VisibilityTester {
    pub fn unoccluded(&self, scene: &Scene) -> bool {
        !scene.intersect_p(&mut self.p0.spawn_ray_to(self.p1))
    }
}

// see point.h

#[derive(Debug,Copy,Clone)]
pub struct PointLight {
    // private data (see point.h)
    pub p_light: Point3f,
    pub i: Spectrum,
    // inherited from class Light (see light.h)
    flags: u8,
    n_samples: i32,
}

impl PointLight {
    pub fn new(light_to_world: &Transform, i: &Spectrum) -> Self {
        PointLight {
            p_light: light_to_world.transform_point(Point3f::default()),
            i: *i,
            flags: LightFlags::DeltaPosition as u8,
            n_samples: 1_i32,
        }
    }
}

impl Light for PointLight {
    fn sample_li(&self,
                 iref: &InteractionCommon,
                 _u: Point2f,
                 wi: &mut Vector3f,
                 pdf: &mut Float,
                 vis: &mut VisibilityTester)
                 -> Spectrum {
        // TODO: ProfilePhase _(Prof::LightSample);
        *wi = vec3_normalize(self.p_light - iref.p);
        *pdf = 1.0 as Float;
        *vis = VisibilityTester {
            p0: InteractionCommon {
                p: iref.p,
                time: iref.time,
                p_error: iref.p_error,
                wo: iref.wo,
                n: iref.n,
            },
            p1: InteractionCommon {
                p: self.p_light,
                time: iref.time,
                p_error: Vector3f::default(),
                wo: Vector3f::default(),
                n: Normal3f::default(),
            },
        };
        self.i / pnt3_distance_squared(self.p_light, iref.p)
    }
    fn power(&self) -> Spectrum {
        Spectrum::default()
    }
    fn preprocess(&self, _scene: &Scene) {
    }
    /// Default implementation returns no emitted radiance for a ray
    /// that escapes the scene bounds.
    fn le(&self, _ray: &mut Ray) -> Spectrum {
        Spectrum::new(0.0 as Float)
    }
    fn pdf_li(&self, _iref: &Interaction, _wi: Vector3f) -> Float {
        0.0 as Float
    }
    fn get_flags(&self) -> u8 {
        self.flags
    }
    fn get_n_samples(&self) -> i32 {
        self.n_samples
    }
}

// see distant.h

#[derive(Debug)]
pub struct DistantLight {
    // private data (see distant.h)
    pub l: Spectrum,
    pub w_light: Vector3f,
    pub world_center: RwLock<Point3f>,
    pub world_radius: RwLock<Float>,
    // inherited from class Light (see light.h)
    flags: u8,
    n_samples: i32,
    // TODO: const MediumInterface mediumInterface;
    light_to_world: Transform,
    world_to_light: Transform,
}

impl DistantLight {
    pub fn new(light_to_world: &Transform, l: &Spectrum, w_light: &Vector3f) -> Self {
        DistantLight {
            l: *l,
            w_light: vec3_normalize(light_to_world.transform_vector(*w_light)),
            world_center: RwLock::new(Point3f::default()),
            world_radius: RwLock::new(0.0),
            flags: LightFlags::DeltaDirection as u8,
            n_samples: 1_i32,
            light_to_world: Transform::default(),
            world_to_light: Transform::default(),
        }
    }
}

impl Light for DistantLight {
    fn sample_li(&self,
                     iref: &InteractionCommon,
                     _u: Point2f,
                     wi: &mut Vector3f,
                     pdf: &mut Float,
                     vis: &mut VisibilityTester)
                     -> Spectrum {
        // TODO: ProfilePhase _(Prof::LightSample);
        *wi = self.w_light;
        *pdf = 1.0 as Float;
        let p_outside: Point3f = iref.p + self.w_light * (2.0 as Float * *self.world_radius.read().unwrap());
        *vis = VisibilityTester {
            p0: InteractionCommon {
                p: iref.p,
                time: iref.time,
                p_error: iref.p_error,
                wo: iref.wo,
                n: iref.n,
            },
            p1: InteractionCommon {
                p: p_outside,
                time: iref.time,
                p_error: Vector3f::default(),
                wo: Vector3f::default(),
                n: Normal3f::default(),
            },
        };
        self.l
    }
    fn power(&self) -> Spectrum {
        Spectrum::default()
    }
    /// Some of the **DistanceLight** methods need to know the bounds
    /// of the scene. Because lights are created before the scene
    /// geometry, these bounds aren't available when the
    /// **DistanceLight** constructor runs. Therefore,
    /// **DistanceLight** implements the optional *preprocess()*
    /// method to get the bound. This method is called at the end of
    /// the **Scene** constructor.
    fn preprocess(&self, scene: &Scene) {
        let mut world_center_ref = self.world_center.write().unwrap();
        let mut world_radius_ref = self.world_radius.write().unwrap();
        Bounds3f::bounding_sphere(&scene.world_bound(),
                                  &mut world_center_ref,
                                  &mut world_radius_ref);
    }
    /// Default implementation returns no emitted radiance for a ray
    /// that escapes the scene bounds.
    fn le(&self, _ray: &mut Ray) -> Spectrum {
        Spectrum::new(0.0 as Float)
    }
    fn pdf_li(&self, _iref: &Interaction, _wi: Vector3f) -> Float {
        0.0 as Float
    }
    fn get_flags(&self) -> u8 {
        self.flags
    }
    fn get_n_samples(&self) -> i32 {
        self.n_samples
    }
}

// see infinte.h

pub struct InfiniteAreaLight {
    // private data (see infinte.h)
    pub lmap: Arc<MipMap>,
    pub world_center: RwLock<Point3f>,
    pub world_radius: RwLock<Float>,
    pub distribution: Arc<Distribution2D>,
    // inherited from class Light (see light.h)
    flags: u8,
    n_samples: i32,
    // TODO: const MediumInterface mediumInterface;
    light_to_world: Transform,
    world_to_light: Transform,
}

impl InfiniteAreaLight {
    pub fn new(light_to_world: &Transform, l: &Spectrum, n_samples: i32, texmap: String) -> Self {
        // read texel data from _texmap_ and initialize _Lmap_
        if texmap != String::from("") {
            // https://cessen.github.io/openexr-rs/openexr/index.html
            let mut resolution: Point2i = Point2i::default();
            let mut names_and_fills: Vec<(&str, f64)> = Vec::new();
            // header
            let mut file = std::fs::File::open(texmap.clone()).unwrap();
            let input_file = InputFile::new(&mut file).unwrap();
            // get resolution
            let (width, height) = input_file.header().data_dimensions();
            resolution.x = width as i32;
            resolution.y = height as i32;
            println!("resolution = {:?}", resolution);
            // make sure the image properties are the same (see incremental_io.rs in github/openexr-rs)
            for channel_name in ["R", "G", "B"].iter() {
                let channel = input_file
                    .header()
                    .get_channel(channel_name)
                    .expect(&format!("Didn't find channel {}.", channel_name));
                assert!(channel.pixel_type == PixelType::HALF);
                names_and_fills.push((channel_name, 0.0_f64));
            }
            let mut pixel_data = vec![(f16::from_f32(0.0), f16::from_f32(0.0), f16::from_f32(0.0)); (resolution.x*resolution.y) as usize];
            {
                // read pixels
                let mut file = std::fs::File::open(texmap.clone()).unwrap();
                let mut input_file = InputFile::new(&mut file).unwrap();
                let mut fb = FrameBufferMut::new(resolution.x as u32, resolution.y as u32);
                fb.insert_channels(&names_and_fills[..], &mut pixel_data);
                input_file.read_pixels(&mut fb).unwrap();
            }
            // convert pixel data into Vec<Spectrum> (and on the way multiply by _l_)
            let mut texels: Vec<Spectrum> = Vec::new();
            for i in 0..(resolution.x*resolution.y) {
                let (r, g, b) = pixel_data[i as usize];
                texels.push(Spectrum::rgb(decode_f16(r.as_bits()),
                                          decode_f16(g.as_bits()),
                                          decode_f16(b.as_bits())) * *l);
            }
            // create _MipMap_ from converted texels (see above)
            let do_trilinear: bool = false;
            let max_aniso: Float = 8.0 as Float;
            let wrap_mode: ImageWrap = ImageWrap::Repeat;
            let lmap = Arc::new(MipMap::new(&resolution, &texels[..], do_trilinear, max_aniso, wrap_mode));

            // initialize sampling PDFs for infinite area light

            // compute scalar-valued image _img_ from environment map
            let width: i32 = 2_i32 * lmap.width();
            let height: i32 = 2_i32 * lmap.height();
            let mut img: Vec<Float> = Vec::new();
            let fwidth: Float = 0.5 as Float / (width as Float).min(height as Float);
            // TODO: ParallelFor(...) {...}
            for v in 0..height {
                let vp: Float = (v as Float + 0.5 as Float) / height as Float;
                let sin_theta: Float = (PI * (v as Float + 0.5 as Float) / height as Float).sin();
                for u in 0..width {
                    let up: Float = (u as Float + 0.5 as Float) / width as Float;
                    let st: Point2f = Point2f { x: up, y: vp, };
                    img.push(lmap.lookup_pnt_flt(&st,
                                                 fwidth).y() * sin_theta);
                }
            }
            let distribution: Arc<Distribution2D> = Arc::new(Distribution2D::new(img, width, height));
            InfiniteAreaLight {
                lmap: lmap,
                world_center: RwLock::new(Point3f::default()),
                world_radius: RwLock::new(0.0),
                distribution: distribution,
                flags: LightFlags::DeltaDirection as u8,
                n_samples: std::cmp::max(1_i32, n_samples),
                light_to_world: *light_to_world,
                world_to_light: Transform::inverse(*light_to_world),
            }
        } else {
            let resolution: Point2i = Point2i {
                x: 1_i32,
                y: 1_i32,
            };
            let texels: Vec<Spectrum> = vec![Spectrum::new(1.0 as Float)];
            let do_trilinear: bool = false;
            let max_aniso: Float = 8.0 as Float;
            let wrap_mode: ImageWrap = ImageWrap::Repeat;
            let lmap = Arc::new(MipMap::new(&resolution, &texels[..], do_trilinear, max_aniso, wrap_mode));

            // initialize sampling PDFs for infinite area light

            // compute scalar-valued image _img_ from environment map
            let width: i32 = 2_i32 * lmap.width();
            let height: i32 = 2_i32 * lmap.height();
            let mut img: Vec<Float> = Vec::new();
            let fwidth: Float = 0.5 as Float / (width as Float).min(height as Float);
            // TODO: ParallelFor(...) {...}
            for v in 0..height {
                let vp: Float = (v as Float + 0.5 as Float) / height as Float;
                let sin_theta: Float = (PI * (v as Float + 0.5 as Float) / height as Float).sin();
                for u in 0..width {
                    let up: Float = (u as Float + 0.5 as Float) / width as Float;
                    let st: Point2f = Point2f { x: up, y: vp, };
                    img.push(lmap.lookup_pnt_flt(&st,
                                                 fwidth).y() * sin_theta);
                }
            }
            let distribution: Arc<Distribution2D> = Arc::new(Distribution2D::new(img, width, height));
            InfiniteAreaLight {
                lmap: lmap,
                world_center: RwLock::new(Point3f::default()),
                world_radius: RwLock::new(0.0),
                distribution: distribution,
                flags: LightFlags::DeltaDirection as u8,
                n_samples: std::cmp::max(1_i32, n_samples),
                light_to_world: Transform::default(),
                world_to_light: Transform::default(),
            }
        }
    }
}

impl Light for InfiniteAreaLight {
    fn sample_li(&self,
                 iref: &InteractionCommon,
                 u: Point2f,
                 wi: &mut Vector3f,
                 pdf: &mut Float,
                 vis: &mut VisibilityTester)
                 -> Spectrum {
        // TODO: ProfilePhase _(Prof::LightSample);
        // find $(u,v)$ sample coordinates in infinite light texture
        let mut map_pdf: Float = 0.0 as Float;
        let uv: Point2f = self.distribution.sample_continuous(&u, &mut map_pdf);
        if map_pdf == 0 as Float {
            return Spectrum::default();
        }
        // convert infinite light sample point to direction
        let theta: Float = uv[1] * PI;
        let phi: Float = uv[0] * 2.0 as Float * PI;
        let cos_theta: Float = theta.cos();
        let sin_theta: Float = theta.sin();
        let sin_phi: Float = phi.sin();
        let cos_phi: Float = phi.cos();
        let vec: Vector3f = Vector3f {
            x: sin_theta * cos_phi,
            y: sin_theta * sin_phi,
            z: cos_theta,
        };
        *wi = self.light_to_world.transform_vector(vec);
        // compute PDF for sampled infinite light direction
        *pdf = map_pdf / (2.0 as Float * PI * PI * sin_theta);
        if sin_theta == 0.0 as Float {
            *pdf = 0.0 as Float;
        }
        // return radiance value for infinite light direction
        let world_radius: Float = *self.world_radius.read().unwrap();
        *vis = VisibilityTester {
            p0: InteractionCommon {
                p: iref.p,
                time: iref.time,
                p_error: iref.p_error,
                wo: iref.wo,
                n: iref.n,
            },
            p1: InteractionCommon {
                p: iref.p + *wi * (2.0 as Float * world_radius),
                time: iref.time,
                p_error: Vector3f::default(),
                wo: Vector3f::default(),
                n: Normal3f::default(),
            },
        };
        // TODO: SpectrumType::Illuminant
        self.lmap.lookup_pnt_flt(&uv, 0.0 as Float)
    }
    /// Like directional lights, the total power from the infinite
    /// area light is related to the surface area of the scene. Like
    /// many other lights the power computed here is approximate.
    fn power(&self) -> Spectrum {
        let p: Point2f = Point2f { x: 0.5, y: 0.5, };
        let world_radius: Float = *self.world_radius.read().unwrap();
        // TODO: SpectrumType::Illuminant
        self.lmap.lookup_pnt_flt(&p, 0.5 as Float) * Spectrum::new(PI * world_radius * world_radius)
    }
    /// Like **DistanceLights**, **InfiniteAreaLights** also need the
    /// scene bounds; here again, the **preprocess()** method finds
    /// the scene bounds after all of the scene geometry has been
    /// created.
    fn preprocess(&self, scene: &Scene) {
        let mut world_center_ref = self.world_center.write().unwrap();
        let mut world_radius_ref = self.world_radius.write().unwrap();
        Bounds3f::bounding_sphere(&scene.world_bound(),
                                  &mut world_center_ref,
                                  &mut world_radius_ref);
    }
    /// Because infinte area lights need to be able to contribute
    /// radiance to rays that don't hit any geometry in the scene,
    /// we'll add a method to the base **Light** class that returns
    /// emitted radiance due to that light along a ray that escapes
    /// the scene bounds. It's the responsibility of the integrators
    /// to call this method for these rays.
    fn le(&self, ray: &mut Ray) -> Spectrum {
        let w: Vector3f = vec3_normalize(self.world_to_light.transform_vector(ray.d));
        let st: Point2f = Point2f {
            x: spherical_phi(&w) * INV_2_PI,
            y: spherical_theta(&w) * INV_PI,
        };
        // TODO: SpectrumType::Illuminant
        self.lmap.lookup_pnt_flt(&st, 0.0 as Float)
    }
    fn pdf_li(&self, _iref: &Interaction, w: Vector3f) -> Float {
        // TODO: ProfilePhase _(Prof::LightPdf);
        let wi: Vector3f = self.world_to_light.transform_vector(w);
        let theta: Float = spherical_theta(&wi);
        let phi: Float = spherical_phi(&wi);
        let sin_theta: Float = theta.sin();
        if sin_theta == 0 as Float {
            return 0 as Float;
        }
        let p: Point2f = Point2f {
            x: phi * INV_2_PI,
            y: theta * INV_PI,
        };
        self.distribution.pdf(&p) / (2.0 as Float * PI * PI * sin_theta)
    }
    fn get_flags(&self) -> u8 {
        self.flags
    }
    fn get_n_samples(&self) -> i32 {
        self.n_samples
    }
}

// see light.h

/// Area lights are light sources defined by one or more **Shapes**
/// that emit light from their surface, with some directional
/// distribution of radiance at each point on the surface.
pub trait AreaLight: Light {
    fn l(&self, intr: &InteractionCommon, w: Vector3f) -> Spectrum;
}

// see diffuse.h

pub struct DiffuseAreaLight {
    pub l_emit: Spectrum,
    pub shape: Arc<Shape + Send + Sync>,
    pub two_sided: bool,
    pub area: Float,
    // inherited from class Light (see light.h)
    flags: u8,
    n_samples: i32,
    // TODO: const MediumInterface mediumInterface;
    // light_to_world: Transform,
    // world_to_light: Transform,
}

impl DiffuseAreaLight {
    pub fn new(_light_to_world: &Transform,
               l_emit: &Spectrum,
               n_samples: i32,
               shape: Arc<Shape + Send + Sync>,
               two_sided: bool) -> Self {
        let area: Float = shape.area();
        DiffuseAreaLight {
            l_emit: *l_emit,
            shape: shape,
            two_sided: two_sided,
            area: area,
            // inherited from class Light (see light.h)
            flags: LightFlags::Area as u8,
            n_samples: std::cmp::max(1_i32, n_samples),
            // TODO: const MediumInterface mediumInterface;
            // light_to_world: *light_to_world,
            // world_to_light: Transform::inverse(*light_to_world),
        }
    }
}

impl Light for DiffuseAreaLight {
    fn sample_li(&self,
                 iref: &InteractionCommon,
                 u: Point2f,
                 wi: &mut Vector3f,
                 pdf: &mut Float,
                 vis: &mut VisibilityTester)
                 -> Spectrum {
        // TODO: ProfilePhase _(Prof::LightSample);
        let p_shape: InteractionCommon = self.shape.sample_with_ref_point(&iref, u, pdf);
        // TODO: iref.mediumInterface = mediumInterface;
        if *pdf == 0.0 as Float ||
            (p_shape.p - iref.p).length_squared() == 0.0 as Float {
                *pdf = 0.0 as Float;
                return Spectrum::default();
            }
        let new_wi: Vector3f = vec3_normalize(p_shape.p - iref.p);
        *wi = new_wi;
        vis.p0 = InteractionCommon {
            p: iref.p,
            time: iref.time,
            p_error: iref.p_error,
            wo: iref.wo,
            n: iref.n,
        };
        vis.p1 = InteractionCommon {
            p: p_shape.p,
            time: p_shape.time,
            p_error: p_shape.p_error,
            wo: p_shape.wo,
            n: p_shape.n,
        };
        self.l(&p_shape, -new_wi)
    }
    fn power(&self) -> Spectrum {
        Spectrum::default()
    }
    fn preprocess(&self, _scene: &Scene) {
        // TODO?
    }
    fn le(&self, _ray: &mut Ray) -> Spectrum {
        Spectrum::default()
    }
    fn pdf_li(&self, iref: &Interaction, wi: Vector3f) -> Float {
        // TODO: ProfilePhase _(Prof::LightPdf);
        self.shape.pdf(iref, wi)
    }
    fn get_flags(&self) -> u8 {
        self.flags
    }
    fn get_n_samples(&self) -> i32 {
        self.n_samples
    }
}

impl AreaLight for DiffuseAreaLight {
    fn l(&self, intr: &InteractionCommon, w: Vector3f) -> Spectrum {
        if self.two_sided || nrm_dot_vec3(intr.n, w) > 0.0 as Float {
            self.l_emit
        } else {
            Spectrum::new(0.0 as Float)
        }
    }
}

// see integrator.h

pub trait SamplerIntegrator {
    // TODO: use Sampler trait
    fn preprocess(&mut self, scene: &Scene, sampler: &mut ZeroTwoSequenceSampler);
    fn li(&self,
          ray: &mut Ray,
          scene: &Scene,
          sampler: &mut ZeroTwoSequenceSampler,
          // arena: &mut Arena,
          depth: i32)
          -> Spectrum;
    fn get_pixel_bounds(&self) -> Bounds2i;
}

// see integrator.cpp

/// Most basic direct lighting strategy.
pub fn uniform_sample_all_lights(it: &SurfaceInteraction,
                                 scene: &Scene,
                                 sampler: &mut ZeroTwoSequenceSampler,
                                 n_light_samples: &Vec<i32>,
                                 handle_media: bool)
                                 -> Spectrum {
    // TODO: ProfilePhase p(Prof::DirectLighting);
    let mut l: Spectrum = Spectrum::new(0.0);
    for j in 0..scene.lights.len() {
        // accumulate contribution of _j_th light to _L_
        let ref light = scene.lights[j];
        let n_samples = n_light_samples[j];
        let u_light_array: Vec<Point2f> = sampler.get_2d_array(n_samples);
        let u_scattering_array: Vec<Point2f> = sampler.get_2d_array(n_samples);
        if u_light_array.is_empty() || u_scattering_array.is_empty() {
            // use a single sample for illumination from _light_
            let u_light: Point2f = sampler.get_2d();
            let u_scattering: Point2f = sampler.get_2d();
            l += estimate_direct(it,
                                 u_scattering,
                                 light.clone(),
                                 u_light,
                                 scene,
                                 sampler, // arena,
                                 handle_media,
                                 false);
        } else {
            // estimate direct lighting using sample arrays
            let mut ld: Spectrum = Spectrum::new(0.0);
            for k in 0..n_samples {
                ld += estimate_direct(it,
                                      u_scattering_array[k as usize],
                                      light.clone(),
                                      u_light_array[k as usize],
                                      scene,
                                      sampler, // arena,
                                      handle_media,
                                      false);
            }
            l += ld / n_samples as Float;
        }
    }
    l
}

/// Estimate direct lighting for only one randomly chosen light and
/// multiply the result by the number of lights to compensate.
pub fn uniform_sample_one_light(it: &SurfaceInteraction,
                                scene: &Scene,
                                sampler: &mut ZeroTwoSequenceSampler,
                                handle_media: bool,
                                light_distrib: Option<&Distribution1D>)
                                -> Spectrum {
    // TODO: ProfilePhase p(Prof::DirectLighting);

    // randomly choose a single light to sample, _light_
    let n_lights: usize = scene.lights.len();
    if n_lights == 0_usize {
        return Spectrum::default();
    }
    let light_num: usize;
    let mut light_pdf: Option<Float> = Some(0.0 as Float);
    let pdf: Float;
    if let Some(light_distribution) = light_distrib {
        // if !light_distrib.is_null() {
        light_num = light_distribution.sample_discrete(sampler.get_1d(), light_pdf.as_mut());
        pdf = light_pdf.unwrap();
        if pdf == 0.0 as Float {
            return Spectrum::default();
        }
    } else {
        light_num = std::cmp::min((sampler.get_1d() * n_lights as Float) as usize, n_lights - 1);
        pdf = 1.0 as Float / n_lights as Float;
    }
    let light = &scene.lights[light_num];
    let u_light: Point2f = sampler.get_2d();
    let u_scattering: Point2f = sampler.get_2d();
    estimate_direct(it, u_scattering, light.clone(), u_light,
                    scene, sampler, handle_media, false) / pdf
}

/// Computes a direct lighting estimate for a single light source sample.
pub fn estimate_direct(it: &SurfaceInteraction,
                       u_scattering: Point2f,
                       light: Arc<Light + Send + Sync>,
                       u_light: Point2f,
                       scene: &Scene,
                       _sampler: &mut ZeroTwoSequenceSampler,
                       // TODO: arena
                       handle_media: bool,
                       specular: bool)
                       -> Spectrum {
    let mut bsdf_flags: u8 = BxdfType::BsdfAll as u8;
    if !specular {
        // bitwise not in Rust is ! (not the ~ operator like in C)
        bsdf_flags = BxdfType::BsdfAll as u8 & !(BxdfType::BsdfSpecular as u8);
    }
    let mut ld: Spectrum = Spectrum::new(0.0);
    // sample light source with multiple importance sampling
    let mut wi: Vector3f = Vector3f::default();
    let mut light_pdf: Float = 0.0 as Float;
    let mut scattering_pdf: Float = 0.0 as Float;
    let mut visibility: VisibilityTester = VisibilityTester::default();
    let it_common: InteractionCommon = InteractionCommon {
        p: it.get_p(),
        time: it.get_time(),
        p_error: it.get_p_error(),
        wo: it.get_wo(),
        n: it.get_n(),
    };
    let mut li: Spectrum = light.sample_li(&it_common, u_light, &mut wi, &mut light_pdf, &mut visibility);
    // TODO: println!("EstimateDirect uLight: {:?} -> Li: {:?}, wi:
    // {:?}, pdf: {:?}", u_light, li, wi, light_pdf);
    if light_pdf > 0.0 as Float && !li.is_black() {
        // compute BSDF or phase function's value for light sample
        let mut f: Spectrum = Spectrum::new(0.0);
        if it.is_surface_interaction() {
            // evaluate BSDF for light sampling strategy
            if let Some(ref bsdf) = it.bsdf {
                f = bsdf.f(it.get_wo(), wi, bsdf_flags) *
                    Spectrum::new(vec3_abs_dot_nrm(wi, it.shading.n));
                scattering_pdf = bsdf.pdf(it.get_wo(), wi, bsdf_flags);
                // TODO: println!("  surf f*dot :{:?}, scatteringPdf: {:?}", f, scattering_pdf);
            }
        } else {
            // evaluate phase function for light sampling strategy
            // TODO
            println!("TODO: evaluate phase function for light sampling strategy");
        }
        if !f.is_black() {
            // compute effect of visibility for light source sample
            if handle_media {
                // TODO: li *= tr(scene, sampler);
                // TODO: VLOG(2) << "  after Tr, Li: " << Li;
            } else {
                if !visibility.unoccluded(scene) {
                    // TODO: println!("  shadow ray blocked");
                    li = Spectrum::new(0.0 as Float);
                } else {
                    // TODO: println!("  shadow ray unoccluded");
                }
            }
            // add light's contribution to reflected radiance
            if !li.is_black() {
                if is_delta_light(light.get_flags()) {
                    ld += f * li / light_pdf;
                } else {
                    let weight: Float = power_heuristic(1_u8, light_pdf, 1_u8, scattering_pdf);
                    ld += f * li * Spectrum::new(weight) / light_pdf;
                }
            }
        }
    }
    // sample BSDF with multiple importance sampling
    if !is_delta_light(light.get_flags()) {
        let mut f: Spectrum = Spectrum::new(0.0);
        let mut sampled_specular: bool = false;
        if it.is_surface_interaction() {
            // sample scattered direction for surface interactions
            let mut sampled_type: u8 = 0_u8;
            if let Some(ref bsdf) = it.bsdf {
                f = bsdf.sample_f(it.get_wo(),
                                  &mut wi,
                                  u_scattering,
                                  &mut scattering_pdf,
                                  bsdf_flags,
                                  &mut sampled_type);
                f *= Spectrum::new(vec3_abs_dot_nrm(wi, it.shading.n));
                sampled_specular = (sampled_type & BxdfType::BsdfSpecular as u8) != 0_u8;
            } else {
                println!("TODO: if let Some(ref bsdf) = it.bsdf failed");
            }
        } else {
            // TODO
            println!("TODO: estimate_direct 1");
        }
        // TODO: println!("  BSDF / phase sampling f: {:?}, scatteringPdf: {:?}",
        //          f, scattering_pdf);
        if !f.is_black() && scattering_pdf > 0.0 {
            // account for light contributions along sampled direction _wi_
            let mut weight: Float = 1.0;
            if !sampled_specular {
                light_pdf = light.pdf_li(it, wi);
                if light_pdf == 0.0 {
                    return ld;
                }
                weight = power_heuristic(1, scattering_pdf, 1, light_pdf);
            }
            // find intersection and compute transmittance
            let mut ray: Ray = it.spawn_ray(wi);
            let tr: Spectrum = Spectrum::new(1.0 as Float);
            let mut found_surface_interaction: bool = false;
            // add light contribution from material sampling
            let mut li: Spectrum = Spectrum::default();
            if handle_media {
                // TODO: scene.IntersectTr(ray, sampler, &lightIsect, &Tr)
            } else {
                if let Some(light_isect) = scene.intersect(&mut ray) {
                    found_surface_interaction = true;
                    if let Some(primitive) = light_isect.primitive {
                        if let Some(area_light) = primitive.get_area_light() {
                            let pa = &*area_light as *const _ as *const usize;
                            let pl = &*light as *const _ as *const usize;
                            if pa == pl {
                                li = light_isect.le(-wi);
                            }
                        }
                    }
                }
            }
            if !found_surface_interaction {
                li = light.le(&mut ray);
            }
            if !li.is_black() {
                ld += f * li * tr * weight / scattering_pdf;
            }
        }
    }
    ld
}

// see sampling.h

#[derive(Debug,Default,Clone)]
pub struct Distribution1D {
    pub func: Vec<Float>,
    pub cdf: Vec<Float>,
    pub func_int: Float,
}

impl Distribution1D {
    pub fn new(f: Vec<Float>) -> Self {
        let n: usize = f.len();
        // compute integral of step function at $x_i$
        let mut cdf: Vec<Float> = Vec::new();
        cdf.push(0.0 as Float);
        for i in 1..(n + 1) {
            let previous: Float = cdf[i - 1];
            cdf.push(previous + f[i - 1] / n as Float);
        }
        // transform step function integral into CDF
        let func_int: Float = cdf[n];
        if func_int == 0.0 as Float {
            for i in 1..(n + 1) {
                cdf[i] = i as Float / n as Float;
            }
        } else {
            for i in 1..(n + 1) {
                cdf[i] /= func_int;
            }
        }
        Distribution1D {
            func: f,
            cdf: cdf,
            func_int: func_int,
        }
    }
    pub fn count(&self) -> usize {
        self.func.len()
    }
    pub fn sample_continuous(&self, u: Float, pdf: Option<&mut Float>, off: Option<&mut usize>) -> Float {
        // find surrounding CDF segments and _offset_
        // int offset = find_interval((int)cdf.size(),
        //                           [&](int index) { return cdf[index] <= u; });

        // see pbrt.h (int FindInterval(int size, const Predicate &pred) {...})
        let mut first: usize = 0;
        let mut len: usize = self.cdf.len();
        while len > 0 as usize {
            let half: usize = len >> 1;
            let middle: usize = first + half;
            // bisect range based on value of _pred_ at _middle_
            if self.cdf[middle] <= u {
                first = middle + 1;
                len -= half + 1;
            } else {
                len = half;
            }
        }
        let offset: usize = clamp_t(first as isize - 1_isize, 0 as isize, self.cdf.len() as isize - 2_isize) as usize;
        if let Some(off_ref) = off {
            *off_ref = offset;
        }
        // compute offset along CDF segment
        let mut du: Float = u - self.cdf[offset];
        if (self.cdf[offset + 1] - self.cdf[offset]) > 0.0 as Float {
            assert!(self.cdf[offset + 1] > self.cdf[offset]);
            du /= self.cdf[offset + 1] - self.cdf[offset];
        }
        assert!(!du.is_nan());
        // compute PDF for sampled offset
        if pdf.is_some() {
            if self.func_int > 0.0 as Float {
                *pdf.unwrap() = self.func[offset] / self.func_int;
            } else {
                *pdf.unwrap() = 0.0;
            }
        }
        // return $x\in{}[0,1)$ corresponding to sample
        (offset as Float + du) / self.count() as Float
    }
    pub fn sample_discrete(&self,
                           u: Float,
                           pdf: Option<&mut Float> /* TODO: Float *uRemapped = nullptr */
    ) -> usize {
        // find surrounding CDF segments and _offset_
        // let offset: usize = find_interval(cdf.size(),
        //                           [&](int index) { return cdf[index] <= u; });

        // see pbrt.h (int FindInterval(int size, const Predicate &pred) {...})
        let mut first: usize = 0;
        let mut len: usize = self.cdf.len();
        while len > 0 as usize {
            let half: usize = len >> 1;
            let middle: usize = first + half;
            // bisect range based on value of _pred_ at _middle_
            if self.cdf[middle] <= u {
                first = middle + 1;
                len -= half + 1;
            } else {
                len = half;
            }
        }
        let offset: usize = clamp_t(first as isize - 1_isize, 0 as isize, self.cdf.len() as isize - 2_isize) as usize;
        if pdf.is_some() {
            if self.func_int > 0.0 as Float {
                *pdf.unwrap() = self.func[offset] / (self.func_int * self.func.len() as Float);
            } else {
                *pdf.unwrap() = 0.0;
            }
        }
        // TODO: if (uRemapped)
        //     *uRemapped = (u - cdf[offset]) / (cdf[offset + 1] - cdf[offset]);
        // if (uRemapped) CHECK(*uRemapped >= 0.f && *uRemapped <= 1.f);
        offset
    }
}

#[derive(Debug,Default,Clone)]
pub struct Distribution2D {
    pub p_conditional_v: Vec<Arc<Distribution1D>>,
    pub p_marginal: Arc<Distribution1D>,
}

impl Distribution2D {
    pub fn new(func: Vec<Float>, nu: i32, nv: i32) -> Self {
        let mut p_conditional_v: Vec<Arc<Distribution1D>> = Vec::new();
        for v in 0..nv {
            // compute conditional sampling distribution for $\tilde{v}$
            let f: Vec<Float> = func[(v * nu) as usize..((v+1) * nu) as usize].to_vec();
            p_conditional_v.push(Arc::new(Distribution1D::new(f)));
        }
        // compute marginal sampling distribution $p[\tilde{v}]$
        let mut marginal_func: Vec<Float> = Vec::with_capacity(nv as usize);
        for v in 0..nv {
            marginal_func.push(p_conditional_v[v as usize].func_int);
        }
        let p_marginal: Arc<Distribution1D> = Arc::new(Distribution1D::new(marginal_func));
        Distribution2D {
            p_conditional_v: p_conditional_v,
            p_marginal: p_marginal,
        }
    }
    pub fn sample_continuous(&self, u: &Point2f, pdf: &mut Float) -> Point2f {
        let mut pdfs: [Float; 2] = [0.0 as Float; 2];
        let mut v: usize = 0_usize;
        let d1: Float = self.p_marginal.sample_continuous(u[1], Some(&mut (pdfs[1])), Some(&mut v));
        let d0: Float = self.p_conditional_v[v].sample_continuous(u[0], Some(&mut (pdfs[0])), None);
        *pdf = pdfs[0] * pdfs[1];
        Point2f {
            x: d0,
            y: d1,
        }
    }
    pub fn pdf(&self, p: &Point2f) -> Float {
        let iu: usize = clamp_t((p[0] as usize * self.p_conditional_v[0].count()) as usize,
                                0_usize,
                                self.p_conditional_v[0].count() - 1_usize);
        let iv: usize = clamp_t((p[1] as usize * self.p_marginal.count()) as usize,
                                0_usize,
                                self.p_marginal.count() - 1_usize);
        self.p_conditional_v[iv].func[iu] / self.p_marginal.func_int
    }
}

/// Randomly permute an array of *count* sample values, each of which
/// has *n_dimensions* dimensions.
pub fn shuffle<T>(samp: &mut [T], count: i32, n_dimensions: i32, rng: &mut Rng) {
    for i in 0..count {
        let other: i32 = i + rng.uniform_uint32_bounded((count - i) as u32) as i32;
        for j in 0..n_dimensions {
            samp.swap((n_dimensions * i + j) as usize,
                      (n_dimensions * other + j) as usize);
        }
    }
}

/// Cosine-weighted hemisphere sampling using Malley's method.
pub fn cosine_sample_hemisphere(u: Point2f) -> Vector3f {
    let d: Point2f = concentric_sample_disk(u);
    let z: Float = (0.0 as Float).max(1.0 as Float - d.x * d.x - d.y * d.y).sqrt();
    Vector3f {
        x: d.x,
        y: d.y,
        z: z,
    }
}

/// Returns a weight of cos_theta / PI.
pub fn cosine_hemisphere_pdf(cos_theta: Float) -> Float {
    cos_theta * INV_PI
}

/// Reducing the variance according to Veach's heuristic.
pub fn power_heuristic(nf: u8, f_pdf: Float, ng: u8, g_pdf: Float) -> Float {
    let f: Float = nf as Float * f_pdf;
    let g: Float = ng as Float  * g_pdf;
    (f * f) / (f * f + g * g)
}

// see sampling.cpp

/// Uniformly sample rays in a hemisphere. Choose a direction.
pub fn uniform_sample_hemisphere(u: Point2f) -> Vector3f {
    let z: Float = u[0_u8];
    let r: Float = (0.0 as Float).max(1.0 as Float - z * z).sqrt();
    let phi: Float = 2.0 as Float * PI * u[1_u8];
    Vector3f {
        x: r * phi.cos(),
        y: r * phi.sin(),
        z: z
    }
}

/// Uniformly sample rays in a hemisphere. Probability density
/// function (PDF).
pub fn uniform_hemisphere_pdf() -> Float {
    INV_2_PI
}

/// Uniformly sample rays in a full sphere. Choose a direction.
pub fn uniform_sample_sphere(u: Point2f) -> Vector3f {
    let z: Float = 1.0 as Float - 2.0 as Float * u[0];
    let r: Float = (0.0 as Float).max(1.0 as Float - z * z).sqrt();
    let phi: Float = 2.0 as Float * PI * u[1];
    Vector3f {
        x: r * phi.cos(),
        y: r * phi.sin(),
        z: z,
    }
}

/// Uniformly distribute samples over a unit disk.
pub fn concentric_sample_disk(u: Point2f) -> Point2f {
    // map uniform random numbers to $[-1,1]^2$
    let u_offset: Point2f = u * 2.0 as Float - Vector2f { x: 1.0 , y: 1.0 };
    // handle degeneracy at the origin
    if u_offset.x == 0.0 as Float && u_offset.y == 0.0 as Float{
        return Point2f::default();
    }
    // apply concentric mapping to point
    let theta: Float;
    let r: Float;
    if u_offset.x.abs() > u_offset.y.abs() {
        r = u_offset.x;
        theta = PI_OVER_4 * (u_offset.y / u_offset.x);
    } else {
        r = u_offset.y;
        theta = PI_OVER_2 - PI_OVER_4 * (u_offset.x / u_offset.y);
    }
    Point2f { x: theta.cos(), y: theta.sin(), } * r
}

/// Uniformly sample rays in a cone of directions. Probability density
/// function (PDF).
pub fn uniform_cone_pdf(cos_theta_max: Float) -> Float {
    1.0 as Float / (2.0 as Float * PI * (1.0 as Float - cos_theta_max))
}

/// Uniformly distributing samples over isosceles right triangles
/// actually works for any triangle.
pub fn uniform_sample_triangle(u: Point2f) -> Point2f {
    let su0: Float = u[0].sqrt();
    Point2f {
        x: 1.0 as Float - su0,
        y: u[1] * su0,
    }
}

// see lightdistrib.h

pub trait LightDistribution {
    /// Given a point |p| in space, this method returns a (hopefully
    /// effective) sampling distribution for light sources at that
    /// point.
    fn lookup(&self, p: Point3f) -> Arc<Distribution1D>;
}

#[derive(Debug,Default)]
struct HashEntry {
    packed_pos: AtomicU64,
    distribution: RwLock<Option<Arc<Distribution1D>>>,
}

pub struct UniformLightDistribution {
    pub distrib: Arc<Distribution1D>,
}

impl UniformLightDistribution {
    pub fn new(scene: &Scene) -> Self {
        let prob: Vec<Float> = vec![1.0 as Float; scene.lights.len()];
        UniformLightDistribution {
            distrib: Arc::new(Distribution1D::new(prob)),
        }
    }
}

impl LightDistribution for UniformLightDistribution {
    fn lookup(&self, _p: Point3f) -> Arc<Distribution1D> {
        self.distrib.clone()
    }
}

/// A spatially-varying light distribution that adjusts the
/// probability of sampling a light source based on an estimate of its
/// contribution to a region of space.  A fixed voxel grid is imposed
/// over the scene bounds and a sampling distribution is computed as
/// needed for each voxel.
pub struct SpatialLightDistribution {
    pub scene: Scene,
    pub n_voxels: [i32; 3],
    hash_table: Arc<Vec<HashEntry>>,
    pub hash_table_size: usize,
}

impl SpatialLightDistribution {
    pub fn new(scene: &Scene, max_voxels: u32) -> Self {
        // compute the number of voxels so that the widest scene
        // bounding box dimension has maxVoxels voxels and the other
        // dimensions have a number of voxels so that voxels are
        // roughly cube shaped.
        let b: Bounds3f = scene.world_bound();
        let diag: Vector3f = b.diagonal();
        let bmax: Float = diag[b.maximum_extent()];
        let mut n_voxels: [i32; 3] = [0_i32; 3];
        for i in 0..3 {
            n_voxels[i] = std::cmp::max(1 as i32, (diag[i as u8] /
                                                   bmax *
                                                   max_voxels as Float).round()
                                        as i32);
            // in the Lookup() method, we require that 20 or fewer
            // bits be sufficient to represent each coordinate
            // value. It's fairly hard to imagine that this would ever
            // be a problem.
            assert!(n_voxels[i] < (1 << 20));
        }
        let hash_table_size: usize =
            (4 as i32 * n_voxels[0] * n_voxels[1] * n_voxels[2]) as usize;
        let mut hash_table: Vec<HashEntry> = Vec::new();
        // let null: *mut Distribution1D = std::ptr::null_mut();
        for _i in 0..hash_table_size {
            let hash_entry: HashEntry = HashEntry {
                packed_pos: AtomicU64::new(INVALID_PACKED_POS),
                distribution: RwLock::new(None),
            };
            hash_table.push(hash_entry);
        }
        println!("SpatialLightDistribution: scene bounds {:?}, voxel res ({:?}, {:?}, {:?})",
                 b, n_voxels[0], n_voxels[1], n_voxels[2]);
        SpatialLightDistribution {
            scene: scene.clone(),
            n_voxels: n_voxels,
            hash_table: Arc::new(hash_table),
            hash_table_size: hash_table_size,
        }
    }
    /// Compute the sampling distribution for the voxel with integer
    /// coordiantes given by "pi".
    pub fn compute_distribution(&self, pi: Point3i) -> Distribution1D {
        // Compute the world-space bounding box of the voxel
        // corresponding to |pi|.
        let p0: Point3f = Point3f { x: pi[0] as Float / self.n_voxels[0] as Float,
                                    y: pi[1] as Float / self.n_voxels[1] as Float,
                                    z: pi[2] as Float / self.n_voxels[2] as Float,
        };
        let p1: Point3f = Point3f { x: (pi[0] + 1) as Float / self.n_voxels[0] as Float,
                                    y: (pi[1] + 1) as Float / self.n_voxels[1] as Float,
                                    z: (pi[2] + 1) as Float / self.n_voxels[2] as Float,
        };
        let voxel_bounds: Bounds3f = Bounds3f {
            p_min: self.scene.world_bound().lerp(p0),
            p_max: self.scene.world_bound().lerp(p1),
        };
        // Compute the sampling distribution. Sample a number of
        // points inside voxelBounds using a 3D Halton sequence; at
        // each one, sample each light source and compute a weight
        // based on Li/pdf for the light's sample (ignoring visibility
        // between the point in the voxel and the point on the light
        // source) as an approximation to how much the light is likely
        // to contribute to illumination in the voxel.
        let n_samples: usize = 128;
        let mut light_contrib: Vec<Float> = vec![0.0 as Float; self.scene.lights.len()];
        for i in 0..n_samples {
            let po: Point3f = voxel_bounds.lerp(Point3f {
                x: radical_inverse(0, i as u64),
                y: radical_inverse(1, i as u64),
                z: radical_inverse(2, i as u64),
            });
            let time: Float = 0.0;
            let intr: InteractionCommon = InteractionCommon {
                p: po,
                time: time,
                p_error: Vector3f::default(),
                wo: Vector3f {
                    x: 1.0,
                    y: 0.0,
                    z: 0.0,
                },
                n: Normal3f::default(),
            };
            // Use the next two Halton dimensions to sample a point on the
            // light source.
            let u: Point2f = Point2f {
                x: radical_inverse(3, i as u64),
                y: radical_inverse(4, i as u64),
            };
            for j in 0..self.scene.lights.len() {
                let mut pdf: Float = 0.0 as Float;
                let mut wi: Vector3f = Vector3f::default();
                let mut vis: VisibilityTester = VisibilityTester::default();
                let li: Spectrum = self.scene.lights[j].sample_li(&intr, u, &mut wi, &mut pdf, &mut vis);
                if pdf > 0.0 as Float {
                    // TODO: look at tracing shadow rays / computing
                    // beam transmittance. Probably shouldn't give
                    // those full weight but instead e.g. have an
                    // occluded shadow ray scale down the contribution
                    // by 10 or something.
                    light_contrib[j] += li.y() / pdf;
                }
            }
        }
        // We don't want to leave any lights with a zero probability;
        // it's possible that a light contributes to points in the
        // voxel even though we didn't find such a point when sampling
        // above. Therefore, compute a minimum (small) weight and
        // ensure that all lights are given at least the corresponding
        // probability.
        let sum_contrib: Float = light_contrib.iter().sum();
        let avg_contrib: Float = sum_contrib / (n_samples * light_contrib.len()) as Float;
        let min_contrib: Float;
        if avg_contrib > 0.0 as Float {
            min_contrib = 0.001 * avg_contrib;
        } else {
            min_contrib = 1.0 as Float;
        }
        for i in 0..light_contrib.len() {
            // println!("Voxel pi = {:?}, light {:?} contrib = {:?}",
            //          pi, i, light_contrib[i]);
            light_contrib[i] = light_contrib[i].max(min_contrib);
        }
        // println!("Initialized light distribution in voxel pi= {:?}, avg_contrib = {:?}",
        //          pi, avg_contrib);
        // Compute a sampling distribution from the accumulated contributions.
        Distribution1D::new(light_contrib)
    }
}

impl LightDistribution for SpatialLightDistribution {
    fn lookup(&self, p: Point3f) -> Arc<Distribution1D> {
        // TODO: ProfilePhase _(Prof::LightDistribLookup);
        // TODO: ++nLookups;

        // first, compute integer voxel coordinates for the given
        // point |p| with respect to the overall voxel grid.
        let offset: Vector3f = self.scene.world_bound().offset(p); // offset in [0,1].
        let mut pi: Point3i = Point3i::default();
        for i in 0..3 {
            // the clamp should almost never be necessary, but is
            // there to be robust to computed intersection points
            // being slightly outside the scene bounds due to
            // floating-point roundoff error.
            pi[i] = clamp_t((offset[i] * self.n_voxels[i as usize] as Float) as i32,
                            0_i32,
                            self.n_voxels[i as usize] - 1_i32);
        }
        // pack the 3D integer voxel coordinates into a single 64-bit value.
        let packed_pos: u64 = ((pi[0] as u64) << 40) | ((pi[1] as u64) << 20) | (pi[2] as u64);
        assert_ne!(packed_pos, INVALID_PACKED_POS);
        // Compute a hash value from the packed voxel coordinates. We
        // could just take packed_Pos mod the hash table size, but
        // since packed_Pos isn't necessarily well distributed on its
        // own, it's worthwhile to do a little work to make sure that
        // its bits values are individually fairly random. For details
        // of and motivation for the following, see:
        // http://zimbry.blogspot.ch/2011/09/better-bit-mixing-improving-on.html
        let mut hash: u64 = packed_pos;
        hash ^= hash >> 31;
        // hash *= 0x7fb5d329728ea185;
        let (mul, _overflow) = hash.overflowing_mul(0x7fb5d329728ea185);
        hash = mul;
        hash ^= hash >> 27;
        // hash *= 0x81dadef4bc2dd44d;
        let (mul, _overflow) = hash.overflowing_mul(0x81dadef4bc2dd44d);
        hash = mul;
        hash ^= hash >> 33;
        hash %= self.hash_table_size as u64;
        // // hash ^= hash >> 31;
        // let (shr, _overflow) = hash.overflowing_shr(31);
        // hash ^= shr;
        // // hash ^= hash >> 27;
        // let (shr, _overflow) = hash.overflowing_shr(27);
        // hash ^= shr;
        // // hash ^= hash >> 33;
        // let (shr, _overflow) = hash.overflowing_shr(33);
        // hash ^= shr;
        // // hash %= self.hash_table_size as u64;
        // let (rem, _overflow) = hash.overflowing_rem(self.hash_table_size as u64);
        // hash = rem;
        // BELOW: comparison is useless due to type limits
        // assert!(hash >= 0_u64, "hash needs to be greater or equal zero");
        // Now, see if the hash table already has an entry for the
        // voxel. We'll use quadratic probing when the hash table
        // entry is already used for another value; step stores the
        // square root of the probe step.
        let mut step: u64 = 1;
        // TODO: int nProbes = 0;
        loop {
            // TODO: ++nProbes;
            let entry: &HashEntry = &self.hash_table[hash as usize];
            // does the hash table entry at offset |hash| match the current point?
            let entry_packed_pos: u64 = entry.packed_pos.load(Ordering::Acquire);
            if entry_packed_pos == packed_pos {
                // Yes! Most of the time, there should already by a light
                // sampling distribution available.
                let option: &Option<Arc<Distribution1D>> = &*entry.distribution.read().unwrap();
                if let Some(ref dist) = *option {
                    // We have a valid sampling distribution.
                    return Arc::clone(dist);
                }
            } else if entry_packed_pos != INVALID_PACKED_POS {
                // The hash table entry we're checking has already
                // been allocated for another voxel. Advance to the
                // next entry with quadratic probing.
                hash += step * step;
                if hash >= self.hash_table_size as u64 {
                    hash %= self.hash_table_size as u64;
                }
                step += 1_u64;
            } else {
                // We have found an invalid entry. (Though this may
                // have changed since the load into entryPackedPos
                // above.)  Use an atomic compare/exchange to try to
                // claim this entry for the current position.
                let invalid: u64 = INVALID_PACKED_POS;
                let success = match entry.packed_pos.compare_exchange_weak(invalid,
                                                                           packed_pos,
                                                                           Ordering::SeqCst,
                                                                           Ordering::Relaxed) {
                    Ok(_) => true,
                    Err(_) => false,
                };
                if success {
                    // Success; we've claimed this position for this
                    // voxel's distribution. Now compute the sampling
                    // distribution and add it to the hash table.
                    let dist: Distribution1D = self.compute_distribution(pi);
                    let arc_dist: Arc<Distribution1D> = Arc::new(dist.clone());
                    let mut distribution = entry.distribution.write().unwrap();
                    *distribution = Some(arc_dist.clone());
                    return arc_dist;
                }
            }
        }
    }
}

// see lightdistrib.cpp

const INVALID_PACKED_POS: u64 = 0xffffffffffffffff;

/// Decides based on the name and the number of scene lights which
/// light distribution to return.
pub fn create_light_sample_distribution(name: String, scene: &Scene)
                                        -> Option<Arc<LightDistribution + Send + Sync>> {
    if name == String::from("uniform") || scene.lights.len() == 1 {
        return Some(Arc::new(UniformLightDistribution::new(scene)));
    } else if name == String::from("power") {
        println!("TODO: PowerLightDistribution");
        // return std::unique_ptr<LightDistribution>{
        //     new PowerLightDistribution(scene)};
    } else if name == String::from("spatial") {
        return Some(Arc::new(SpatialLightDistribution::new(scene, 64)));
    } else {
        println!("Light sample distribution type \"{:?}\" unknown. Using \"spatial\".",
                 name);
        // return std::unique_ptr<LightDistribution>{
        //     new SpatialLightDistribution(scene)};
    }
    None
}

// see path.h

/// Path Tracing (Global Illumination)
pub struct PathIntegrator {
    // inherited from SamplerIntegrator (see integrator.h)
    pixel_bounds: Bounds2i,
    // see path.h
    max_depth: u32,
    rr_threshold: Float, // 1.0
    light_sample_strategy: String, // "spatial"
    // TODO: std::unique_ptr<LightDistribution> lightDistribution;
    light_distribution: Option<Arc<LightDistribution + Send + Sync>>,
}

impl PathIntegrator {
    pub fn new(max_depth: u32,
               _perspective_camera: &PerspectiveCamera,
               _sampler: &ZeroTwoSequenceSampler,
               pixel_bounds: Bounds2i,
               rr_threshold: Float,
               light_sample_strategy: String)
               -> Self
    {
        PathIntegrator {
            pixel_bounds: pixel_bounds,
            max_depth: max_depth,
            rr_threshold: rr_threshold,
            light_sample_strategy: light_sample_strategy,
            light_distribution: None,
        }
    }
}

impl SamplerIntegrator for PathIntegrator {
    fn preprocess(&mut self, scene: &Scene, _sampler: &mut ZeroTwoSequenceSampler) {
        self.light_distribution =
            create_light_sample_distribution(self.light_sample_strategy.clone(),
                                             scene);
    }
    fn li(&self,
          r: &mut Ray,
          scene: &Scene,
          sampler: &mut ZeroTwoSequenceSampler,
          // arena: &mut Arena,
          _depth: i32)
          -> Spectrum {
        // TODO: ProfilePhase p(Prof::SamplerIntegratorLi);
        let mut l: Spectrum = Spectrum::default();
        let mut beta: Spectrum = Spectrum::new(1.0 as Float);
        let mut ray: Ray = Ray {
            o: r.o,
            d: r.d,
            t_max: r.t_max,
            time: r.time,
            differential: r.differential,
        };
        let mut specular_bounce: bool = false;
        let mut bounces: u32 = 0_u32;
        // Added after book publication: etaScale tracks the
        // accumulated effect of radiance scaling due to rays passing
        // through refractive boundaries (see the derivation on p. 527
        // of the third edition). We track this value in order to
        // remove it from beta when we apply Russian roulette; this is
        // worthwhile, since it lets us sometimes avoid terminating
        // refracted rays that are about to be refracted back out of a
        // medium and thus have their beta value increased.
        let mut eta_scale: Float = 1.0;
        loop {
            // find next path vertex and accumulate contribution
            // println!("Path tracer bounce {:?}, current L = {:?}, beta = {:?}",
            //          bounces, l, beta);
            // intersect _ray_ with scene and store intersection in _isect_
            if let Some(mut isect) = scene.intersect(&mut ray) {
                // possibly add emitted light at intersection
                if bounces == 0 || specular_bounce {
                    // add emitted light at path vertex
                    l += beta * isect.le(-ray.d);
                    // println!("Added Le -> L = {:?}", l);
                }
                // terminate path if _maxDepth_ was reached
                if bounces >= self.max_depth {
                    break;
                }
                // compute scattering functions and skip over medium boundaries
                let mode: TransportMode = TransportMode::Radiance;
                isect.compute_scattering_functions(&mut ray, true, mode);
                // if (!isect.bsdf) {
                //     VLOG(2) << "Skipping intersection due to null bsdf";
                //     ray = isect.SpawnRay(ray.d);
                //     bounces--;
                //     continue;
                // }
                if let Some(ref light_distribution) = self.light_distribution {
                    let distrib: Arc<Distribution1D> = light_distribution.lookup(isect.p);
                    // Sample illumination from lights to find path contribution.
                    // (But skip this for perfectly specular BSDFs.)
                    let bsdf_flags: u8 = BxdfType::BsdfAll as u8 & !(BxdfType::BsdfSpecular as u8);
                    if let Some(ref bsdf) = isect.bsdf {
                        if bsdf.num_components(bsdf_flags) > 0 {
                            // TODO: ++total_paths;
                            let ld: Spectrum = beta * uniform_sample_one_light(&isect,
                                                                               scene,
                                                                               sampler,
                                                                               false,
                                                                               Some(Arc::borrow(&distrib)));
                            // TODO: println!("Sampled direct lighting Ld = {:?}", ld);
                            // TODO: if ld.is_black() {
                            //     ++zero_radiance_paths;
                            // }
                            assert!(ld.y() >= 0.0 as Float, "ld = {:?}", ld);
                            l += ld;
                        }
                        // Sample BSDF to get new path direction
                        let wo: Vector3f = -ray.d;
                        let mut wi: Vector3f = Vector3f::default();
                        let mut pdf: Float = 0.0 as Float;
                        let bsdf_flags: u8 = BxdfType::BsdfAll as u8;
                        let mut sampled_type: u8 = u8::max_value(); // != 0
                        let f: Spectrum = bsdf.sample_f(wo,
                                                        &mut wi,
                                                        sampler.get_2d(),
                                                        &mut pdf,
                                                        bsdf_flags,
                                                        &mut sampled_type);

                        // println!("Sampled BSDF, f = {:?}, pdf = {:?}", f, pdf);
                        if f.is_black() || pdf == 0.0 as Float {
                            break;
                        }
                        beta *= (f * vec3_abs_dot_nrm(wi, isect.shading.n)) / pdf;
                        // println!("Updated beta = {:?}", beta);
                        assert!(beta.y() >= 0.0 as Float);
                        assert!(!(beta.y().is_infinite()));
                        specular_bounce = (sampled_type & BxdfType::BsdfSpecular as u8) != 0_u8;
                        if ((sampled_type & BxdfType::BsdfSpecular as u8) != 0_u8) &&
                            ((sampled_type & BxdfType::BsdfTransmission as u8) != 0_u8)
                        {
                            let eta: Float = bsdf.eta;
                            // Update the term that tracks radiance
                            // scaling for refraction depending on
                            // whether the ray is entering or leaving
                            // the medium.
                            if vec3_dot_nrm(wo, isect.n) > 0.0 as Float {
                                eta_scale *= eta * eta;
                            } else {
                                eta_scale *= 1.0 as Float / (eta * eta);
                            }
                        }
                        ray = isect.spawn_ray(wi);

                        // Account for subsurface scattering, if applicable
                        // TODO: if (isect.bssrdf && ((sampled_type & BxdfType::BsdfTransmission as u8) != 0_u8)) {
                        // Importance sample the BSSRDF
                        //     SurfaceInteraction pi;
                        //     Spectrum S = isect.bssrdf->Sample_S(
                        //         scene, sampler.Get1D(), sampler.Get2D(), arena, &pi, &pdf);
                        //     DCHECK(!std::isinf(beta.y()));
                        //     if (S.IsBlack() || pdf == 0) break;
                        //     beta *= S / pdf;
                        //     // Account for the direct subsurface scattering component
                        //     L += beta * UniformSampleOneLight(pi, scene, arena, sampler, false,
                        //                                       lightDistribution->Lookup(pi.p));
                        //     // Account for the indirect subsurface scattering component
                        //     Spectrum f = pi.bsdf->Sample_f(pi.wo, &wi, sampler.Get2D(), &pdf,
                        //                                    BSDF_ALL, &sampled_type);
                        //     if (f.IsBlack() || pdf == 0) break;
                        //     beta *= f * AbsDot(wi, pi.shading.n) / pdf;
                        //     DCHECK(!std::isinf(beta.y()));
                        //     specularBounce = (sampled_type & BSDF_SPECULAR) != 0;
                        //     ray = pi.SpawnRay(wi);
                        // }

                        // Possibly terminate the path with Russian roulette.
                        // Factor out radiance scaling due to refraction in rr_beta.
                        let rr_beta: Spectrum = beta * eta_scale;
                        if rr_beta.max_component_value() < self.rr_threshold && bounces > 3 {
                            let q: Float = (0.05 as Float).max(1.0 as Float - rr_beta.max_component_value());
                            if sampler.get_1d() < q {
                                break;
                            }
                            beta = beta / (1.0 as Float - q);
                            assert!(!(beta.y().is_infinite()));
                        }
                    } else {
                        println!("TODO: if let Some(ref bsdf) = isect.bsdf failed");
                    }
                }
            } else {
                // add emitted light from the environment
                if bounces == 0 || specular_bounce {
                    // for (const auto &light : scene.infiniteLights)
                    for light in &scene.infinite_lights {
                        l += beta * light.le(&mut ray);
                    }
                    // println!("Added infinite area lights -> L = {:?}", l);
                }
                // terminate path if ray escaped
                break;
            }
            bounces += 1_u32;
        }
        l
    }
    fn get_pixel_bounds(&self) -> Bounds2i {
        self.pixel_bounds
    }
}

// see ao.h

/// Ambient Occlusion
pub struct AOIntegrator {
    // inherited from SamplerIntegrator (see integrator.h)
    pixel_bounds: Bounds2i,
    // see ao.h
    cos_sample: bool,
    n_samples: i32,
}

impl AOIntegrator {
    pub fn new(cos_sample: bool,
               n_samples: i32,
               _perspective_camera: &PerspectiveCamera,
               _sampler: &ZeroTwoSequenceSampler,
               pixel_bounds: Bounds2i)
               -> Self {
        AOIntegrator {
            pixel_bounds: pixel_bounds,
            cos_sample: cos_sample,
            n_samples: n_samples,
        }
    }
}

impl SamplerIntegrator for AOIntegrator {
    fn preprocess(&mut self, _scene: &Scene, sampler: &mut ZeroTwoSequenceSampler) {
        sampler.request_2d_array(self.n_samples);
    }
    fn li(&self,
          r: &mut Ray,
          scene: &Scene,
          sampler: &mut ZeroTwoSequenceSampler,
          // arena: &mut Arena,
          _depth: i32)
          -> Spectrum {
        // TODO: ProfilePhase p(Prof::SamplerIntegratorLi);
        let mut l: Spectrum = Spectrum::default();
        let mut ray: Ray = Ray {
            o: r.o,
            d: r.d,
            t_max: r.t_max,
            time: r.time,
            differential: r.differential,
        };
        if let Some(mut isect) = scene.intersect(&mut ray) {
            let mode: TransportMode = TransportMode::Radiance;
            isect.compute_scattering_functions(&mut ray, true, mode);
            // if (!isect.bsdf) {
            //     VLOG(2) << "Skipping intersection due to null bsdf";
            //     ray = isect.SpawnRay(ray.d);
            //     goto retry;
            // }
            // compute coordinate frame based on true geometry, not
            // shading geometry.
            let n: Normal3f = nrm_faceforward_vec3(isect.n, -ray.d);
            let s: Vector3f = vec3_normalize(isect.dpdu);
            let t: Vector3f = nrm_cross_vec3(isect.n, s);
            let u: Vec<Point2f> = sampler.get_2d_array(self.n_samples);
            for i in 0..self.n_samples as usize {
                // Vector3f wi;
                let mut wi: Vector3f;
                let pdf: Float;
                if self.cos_sample {
                    wi = cosine_sample_hemisphere(u[i]);
                    pdf = cosine_hemisphere_pdf(wi.z.abs());
                } else {
                    wi = uniform_sample_hemisphere(u[i]);
                    pdf = uniform_hemisphere_pdf();
                }
                // transform wi from local frame to world space.
                wi = Vector3f {
                    x: s.x * wi.x + t.x * wi.y + n.x * wi.z,
                    y: s.y * wi.x + t.y * wi.y + n.y * wi.z,
                    z: s.z * wi.x + t.z * wi.y + n.z * wi.z,
                };
                let mut ray: Ray = isect.spawn_ray(wi);
                if !scene.intersect_p(&mut ray) {
                    l += Spectrum::new(vec3_dot_nrm(wi, n) / (pdf * self.n_samples as Float));
                }
            }
        }
        l
    }
    fn get_pixel_bounds(&self) -> Bounds2i {
        self.pixel_bounds
    }
}

// see directlighting.h

#[derive(Debug,Clone,PartialEq)]
pub enum LightStrategy {
    UniformSampleAll,
    UniformSampleOne,
}

/// Direct Lighting (no Global Illumination)
pub struct DirectLightingIntegrator {
    // inherited from SamplerIntegrator (see integrator.h)
    pixel_bounds: Bounds2i,
    // see directlighting.h
    strategy: LightStrategy,
    max_depth: i64,
    n_light_samples: Vec<i32>,
}

impl DirectLightingIntegrator {
    pub fn new(strategy: LightStrategy,
               max_depth: i64,
               pixel_bounds: Bounds2i)
               -> Self {
        DirectLightingIntegrator {
            pixel_bounds: pixel_bounds,
            strategy: strategy,
            max_depth: max_depth,
            n_light_samples: Vec::new(),
        }
    }
    pub fn specular_reflect(&self,
                            ray: &Ray,
                            isect: &SurfaceInteraction,
                            scene: &Scene,
                            sampler: &mut ZeroTwoSequenceSampler,
                            // arena: &mut Arena,
                            depth: i32)
                            -> Spectrum {
        // compute specular reflection direction _wi_ and BSDF value
        let wo: Vector3f = isect.wo;
        let mut wi: Vector3f = Vector3f::default();
        let mut pdf: Float = 0.0 as Float;
        let ns: Normal3f = isect.shading.n;
        let mut sampled_type: u8 = 0_u8;
        let bsdf_flags: u8 = BxdfType::BsdfReflection as u8 | BxdfType::BsdfSpecular as u8;
        let f: Spectrum;
        if let Some(ref bsdf) = isect.bsdf {
            f = bsdf.sample_f(wo,
                              &mut wi,
                              sampler.get_2d(),
                              &mut pdf,
                              bsdf_flags,
                              &mut sampled_type);
            if pdf > 0.0 as Float && !f.is_black() && vec3_abs_dot_nrm(wi, ns) != 0.0 as Float {
                // compute ray differential _rd_ for specular reflection
                let mut rd: Ray = isect.spawn_ray(wi);
                if let Some(d) = ray.differential.iter().next() {
                    let dndx: Normal3f = isect.shading.dndu * isect.dudx +
                                         isect.shading.dndv * isect.dvdx;
                    let dndy: Normal3f = isect.shading.dndu * isect.dudy +
                                         isect.shading.dndv * isect.dvdy;
                    let dwodx: Vector3f = -d.rx_direction - wo;
                    let dwody: Vector3f = -d.ry_direction - wo;
                    let ddndx: Float = vec3_dot_nrm(dwodx, ns) + vec3_dot_nrm(wo, dndx);
                    let ddndy: Float = vec3_dot_nrm(dwody, ns) + vec3_dot_nrm(wo, dndy);
                    // compute differential reflected directions
                    let diff: RayDifferential = RayDifferential {
                        rx_origin: isect.p + isect.dpdx,
                        ry_origin: isect.p + isect.dpdy,
                        rx_direction: wi - dwodx +
                                      Vector3f::from(dndx * vec3_dot_nrm(wo, ns) + ns * ddndx) *
                                      2.0 as Float,
                        ry_direction: wi - dwody +
                                      Vector3f::from(dndy * vec3_dot_nrm(wo, ns) + ns * ddndy) *
                                      2.0 as Float,
                    };
                    rd.differential = Some(diff);
                }
                return f * self.li(&mut rd, scene, sampler, depth + 1) *
                       Spectrum::new(vec3_abs_dot_nrm(wi, ns) / pdf);
            } else {
                Spectrum::new(0.0)
            }
        } else {
            Spectrum::new(0.0)
        }
    }
    pub fn specular_transmit(&self,
                             ray: &Ray,
                             isect: &SurfaceInteraction,
                             scene: &Scene,
                             sampler: &mut ZeroTwoSequenceSampler,
                             // arena: &mut Arena,
                             depth: i32)
                             -> Spectrum {
        let wo: Vector3f = isect.wo;
        let mut wi: Vector3f = Vector3f::default();
        let mut pdf: Float = 0.0 as Float;
        // let p: Point3f = isect.p;
        let ns: Normal3f = isect.shading.n;
        let mut sampled_type: u8 = 0_u8;
        let bsdf_flags: u8 = BxdfType::BsdfTransmission as u8 | BxdfType::BsdfSpecular as u8;
        let f: Spectrum;
        if let Some(ref bsdf) = isect.bsdf {
            f = bsdf.sample_f(wo,
                              &mut wi,
                              sampler.get_2d(),
                              &mut pdf,
                              bsdf_flags,
                              &mut sampled_type);
            if pdf > 0.0 as Float && !f.is_black() && vec3_abs_dot_nrm(wi, ns) != 0.0 as Float {
                // compute ray differential _rd_ for specular transmission
                let mut rd: Ray = isect.spawn_ray(wi);
                if let Some(d) = ray.differential.iter().next() {
                    let mut eta: Float = bsdf.eta;
                    let w: Vector3f = -wo;
                    if vec3_dot_nrm(wo, ns) < 0.0 as Float {
                        eta = 1.0 / eta;
                    }
                    let dndx: Normal3f = isect.shading.dndu * isect.dudx +
                                         isect.shading.dndv * isect.dvdx;
                    let dndy: Normal3f = isect.shading.dndu * isect.dudy +
                                         isect.shading.dndv * isect.dvdy;
                    let dwodx: Vector3f = -d.rx_direction - wo;
                    let dwody: Vector3f = -d.ry_direction - wo;
                    let ddndx: Float = vec3_dot_nrm(dwodx, ns) + vec3_dot_nrm(wo, dndx);
                    let ddndy: Float = vec3_dot_nrm(dwody, ns) + vec3_dot_nrm(wo, dndy);
                    let mu: Float = eta * vec3_dot_nrm(w, ns) - vec3_dot_nrm(wi, ns);
                    let dmudx: Float = (eta - (eta * eta * vec3_dot_nrm(w, ns)) / vec3_dot_nrm(wi, ns)) * ddndx;
                    let dmudy: Float = (eta - (eta * eta * vec3_dot_nrm(w, ns)) / vec3_dot_nrm(wi, ns)) * ddndy;
                    let diff: RayDifferential = RayDifferential {
                        rx_origin: isect.p + isect.dpdx,
                        ry_origin: isect.p + isect.dpdy,
                        rx_direction: wi + dwodx * eta - Vector3f::from(dndx * mu + ns * dmudx),
                        ry_direction: wi + dwody * eta - Vector3f::from(dndy * mu + ns * dmudy),
                    };
                    rd.differential = Some(diff);
                }
                return f * self.li(&mut rd, scene, sampler, depth + 1) *
                       Spectrum::new(vec3_abs_dot_nrm(wi, ns) / pdf);
            } else {
                Spectrum::new(0.0)
            }
        } else {
            Spectrum::new(0.0)
        }
    }
}

impl SamplerIntegrator for DirectLightingIntegrator {
    fn preprocess(&mut self, scene: &Scene, sampler: &mut ZeroTwoSequenceSampler) {
        if self.strategy == LightStrategy::UniformSampleAll {
            // compute number of samples to use for each light
            for li in 0..scene.lights.len() {
                let ref light = scene.lights[li];
                self.n_light_samples.push(sampler.round_count(light.get_n_samples()));
            }
            // request samples for sampling all lights
            for _i in 0..self.max_depth {
                for j in 0..scene.lights.len() {
                    sampler.request_2d_array(self.n_light_samples[j]);
                    sampler.request_2d_array(self.n_light_samples[j]);
                }
            }
        }
    }
    fn li(&self,
          ray: &mut Ray,
          scene: &Scene,
          sampler: &mut ZeroTwoSequenceSampler,
          // arena: &mut Arena,
          depth: i32)
          -> Spectrum {
        // TODO: ProfilePhase p(Prof::SamplerIntegratorLi);
        let mut l: Spectrum = Spectrum::new(0.0 as Float);
        // find closest ray intersection or return background radiance
        if let Some(mut isect) = scene.intersect(ray) {
            // compute scattering functions for surface interaction
            let mode: TransportMode = TransportMode::Radiance;
            isect.compute_scattering_functions(ray /* arena, */, false, mode);
            // if (!isect.bsdf)
            //     return Li(isect.SpawnRay(ray.d), scene, sampler, arena, depth);
            let wo: Vector3f = isect.wo;
            l += isect.le(wo);
            if scene.lights.len() > 0 {
                // compute direct lighting for _DirectLightingIntegrator_ integrator
                if self.strategy == LightStrategy::UniformSampleAll {
                    l += uniform_sample_all_lights(&isect,
                                                   scene,
                                                   sampler,
                                                   &self.n_light_samples,
                                                   false);
                } else {
                    l += uniform_sample_one_light(&isect,
                                                  scene,
                                                  sampler,
                                                  false,
                                                  None);
                }
            }
            if ((depth + 1_i32) as i64) < self.max_depth {
                // trace rays for specular reflection and refraction
                l += self.specular_reflect(ray, &isect, scene, sampler, // arena,
                                           depth);
                l += self.specular_transmit(ray, &isect, scene, sampler, // arena,
                                            depth);
            }
        } else {
            for light in &scene.lights {
                l += light.le(ray);
            }
        }
        l
    }
    fn get_pixel_bounds(&self) -> Bounds2i {
        self.pixel_bounds
    }
}

// see paramset.h

pub struct ParamSetItem<T> {
    pub name: String,
    pub values: Vec<T>,
    pub n_values: usize,
    pub looked_up: bool, // false
}

#[derive(Default)]
pub struct ParamSet {
    pub key_word: String,
    pub name: String,
    pub tex_type: String,
    pub tex_name: String,
    pub bools: Vec<ParamSetItem<bool>>,
    pub ints: Vec<ParamSetItem<i32>>,
    pub floats: Vec<ParamSetItem<Float>>,
    pub point2fs: Vec<ParamSetItem<Point2f>>,
    pub vector2fs: Vec<ParamSetItem<Vector2f>>,
    pub point3fs: Vec<ParamSetItem<Point3f>>,
    pub vector3fs: Vec<ParamSetItem<Vector3f>>,
    pub normals: Vec<ParamSetItem<Normal3f>>,
    pub spectra: Vec<ParamSetItem<Spectrum>>,
    pub strings: Vec<ParamSetItem<String>>,
    pub textures: Vec<ParamSetItem<String>>,
}

impl ParamSet {
    pub fn reset(&mut self, key_word: String, name: String, tex_type: String, tex_name: String) {
        self.key_word = key_word;
        self.name = name;
        self.tex_type = tex_type;
        self.tex_name = tex_name;
        self.bools.clear();
        self.ints.clear();
        self.floats.clear();
        self.point2fs.clear();
        self.vector2fs.clear();
        self.point3fs.clear();
        self.vector3fs.clear();
        self.normals.clear();
        self.spectra.clear();
        self.strings.clear();
        self.textures.clear();
    }
    pub fn add_float(&mut self, name: String, value: Float) {
        self.floats.push(ParamSetItem::<Float> {
            name: name,
            values: vec!(value),
            n_values: 1_usize,
            looked_up: false,
        });
    }
    pub fn add_floats(&mut self, name: String, values: Vec<Float>) {
        let n_values: usize = values.len();
        self.floats.push(ParamSetItem::<Float> {
            name: name,
            values: values,
            n_values: n_values,
            looked_up: false,
        });
    }
    pub fn add_int(&mut self, name: String, value: i32) {
        self.ints.push(ParamSetItem::<i32> {
            name: name,
            values: vec!(value),
            n_values: 1_usize,
            looked_up: false,
        });
    }
    pub fn add_ints(&mut self, name: String, values: Vec<i32>) {
        let n_values: usize = values.len();
        self.ints.push(ParamSetItem::<i32> {
            name: name,
            values: values,
            n_values: n_values,
            looked_up: false,
        });
    }
    pub fn add_bool(&mut self, name: String, value: bool) {
        self.bools.push(ParamSetItem::<bool> {
            name: name,
            values: vec!(value),
            n_values: 1_usize,
            looked_up: false,
        });
    }
    pub fn add_point3f(&mut self, name: String, value: Point3f) {
        self.point3fs.push(ParamSetItem::<Point3f> {
            name: name,
            values: vec!(value),
            n_values: 1_usize,
            looked_up: false,
        });
    }
    pub fn add_point3fs(&mut self, name: String, values: Vec<Float>) {
        let n_values: usize = values.len();
        let mut p_values: Vec<Point3f> = Vec::new();
        let n_points: usize = values.len() / 3_usize;
        assert!(n_values % 3 == 0, "point parameters need 3 coordinates");
        for i in 0..n_points {
            let x: Float = values[i*3+0];
            let y: Float = values[i*3+1];
            let z: Float = values[i*3+2];
            p_values.push(Point3f { x: x, y: y, z: z, });
        }
        self.point3fs.push(ParamSetItem::<Point3f> {
            name: name,
            values: p_values,
            n_values: n_points,
            looked_up: false,
        });
    }
    pub fn add_string(&mut self, name: String, value: String) {
        self.strings.push(ParamSetItem::<String> {
            name: name,
            values: vec!(value),
            n_values: 1_usize,
            looked_up: false,
        });
    }
    pub fn add_texture(&mut self, name: String, value: String) {
        self.textures.push(ParamSetItem::<String> {
            name: name,
            values: vec!(value),
            n_values: 1_usize,
            looked_up: false,
        });
    }
    pub fn add_vector3f(&mut self, name: String, value: Vector3f) {
        self.vector3fs.push(ParamSetItem::<Vector3f> {
            name: name,
            values: vec!(value),
            n_values: 1_usize,
            looked_up: false,
        });
    }
    pub fn add_normal3f(&mut self, name: String, value: Normal3f) {
        self.normals.push(ParamSetItem::<Normal3f> {
            name: name,
            values: vec!(value),
            n_values: 1_usize,
            looked_up: false,
        });
    }
    pub fn add_normal3fs(&mut self, name: String, values: Vec<Float>) {
        let n_values: usize = values.len();
        let mut p_values: Vec<Normal3f> = Vec::new();
        let n_normals: usize = values.len() / 3_usize;
        assert!(n_values % 3 == 0, "normal parameters need 3 coordinates");
        for i in 0..n_normals {
            let x: Float = values[i*3+0];
            let y: Float = values[i*3+1];
            let z: Float = values[i*3+2];
            p_values.push(Normal3f { x: x, y: y, z: z, });
        }
        self.normals.push(ParamSetItem::<Normal3f> {
            name: name,
            values: p_values,
            n_values: n_normals,
            looked_up: false,
        });
    }
    pub fn add_rgb_spectrum(&mut self, name: String, value: Spectrum) {
        self.spectra.push(ParamSetItem::<Spectrum> {
            name: name,
            values: vec!(value),
            n_values: 1_usize,
            looked_up: false,
        });
    }
    pub fn copy_from(&mut self, param_set: &ParamSet) {
        self.key_word = param_set.key_word.clone();
        // self.name = param_set.name.clone();
        self.bools.clear();
        self.ints.clear();
        for i in &param_set.ints {
            let mut values: Vec<i32> = Vec::new();
            for ix in 0..i.n_values {
                values.push(i.values[ix]);
            }
            self.ints.push(ParamSetItem::<i32> {
                name: i.name.clone(),
                values: values,
                n_values: i.n_values,
                looked_up: false,
            });
        }
        self.floats.clear();
        for f in &param_set.floats {
            let mut values: Vec<Float> = Vec::new();
            for ix in 0..f.n_values {
                values.push(f.values[ix]);
            }
            self.floats.push(ParamSetItem::<Float> {
                name: f.name.clone(),
                values: values,
                n_values: f.n_values,
                looked_up: false,
            });
        }
        self.point2fs.clear();
        self.vector2fs.clear();
        self.point3fs.clear();
        for p in &param_set.point3fs {
            let mut values: Vec<Point3f> = Vec::new();
            for ix in 0..p.n_values {
                values.push(p.values[ix].clone());
            }
            self.point3fs.push(ParamSetItem::<Point3f> {
                name: p.name.clone(),
                values: values,
                n_values: p.n_values,
                looked_up: false,
            });
        }
        self.vector3fs.clear();
        self.normals.clear();
        self.spectra.clear();
        for s in &param_set.spectra {
            let mut values: Vec<Spectrum> = Vec::new();
            for ix in 0..s.n_values {
                values.push(s.values[ix].clone());
            }
            self.spectra.push(ParamSetItem::<Spectrum> {
                name: s.name.clone(),
                values: values,
                n_values: s.n_values,
                looked_up: false,
            });
        }
        self.strings.clear();
        for s in &param_set.strings {
            let mut values: Vec<String> = Vec::new();
            for ix in 0..s.n_values {
                values.push(s.values[ix].clone());
            }
            self.strings.push(ParamSetItem::<String> {
                name: s.name.clone(),
                values: values,
                n_values: s.n_values,
                looked_up: false,
            });
        }
        self.textures.clear();
        for s in &param_set.textures {
            let mut values: Vec<String> = Vec::new();
            for ix in 0..s.n_values {
                values.push(s.values[ix].clone());
            }
            self.textures.push(ParamSetItem::<String> {
                name: s.name.clone(),
                values: values,
                n_values: s.n_values,
                looked_up: false,
            });
        }
    }
    pub fn find_one_float(&self, name: String, d: Float) -> Float {
        for v in &self.floats {
            if v.name == name && v.n_values == 1 {
                // v.looked_up = true;
                return v.values[0];
            }
        }
        d
    }
    pub fn find_one_int(&self, name: String, d: i32) -> i32 {
        for v in &self.ints {
            if v.name == name && v.n_values == 1 {
                // v.looked_up = true;
                return v.values[0];
            }
        }
        d
    }
    pub fn find_one_bool(&self, name: String, d: bool) -> bool {
        for v in &self.bools {
            if v.name == name && v.n_values == 1 {
                // v.looked_up = true;
                return v.values[0];
            }
        }
        d
    }
    pub fn find_one_point3f(&self, name: String, d: Point3f) -> Point3f {
        for v in &self.point3fs {
            if v.name == name && v.n_values == 1 {
                // v.looked_up = true;
                return v.values[0];
            }
        }
        d
    }
    pub fn find_one_vector3f(&self, name: String, d: Vector3f) -> Vector3f {
        for v in &self.vector3fs {
            if v.name == name && v.n_values == 1 {
                // v.looked_up = true;
                return v.values[0];
            }
        }
        d
    }
    pub fn find_one_spectrum(&self, name: String, d: Spectrum) -> Spectrum {
        for v in &self.spectra {
            if v.name == name && v.n_values == 1 {
                // v.looked_up = true;
                return v.values[0];
            }
        }
        d
    }
    pub fn find_one_string(&self, name: String, d: String) -> String {
        for v in &self.strings {
            if v.name == name && v.n_values == 1 {
                // v.looked_up = true;
                return v.values[0].clone();
            }
        }
        d
    }
    pub fn find_one_filename(&self, name: String, d: String) -> String {
        let filename: String = self.find_one_string(name, String::new());
        if filename == String::new() {
            return d;
        }
        // TODO: filename = AbsolutePath(ResolveFilename(filename));
        filename
    }
    pub fn find_texture(&self, name: String) -> String {
        let d: String = String::new();
        lookup_one(&self.textures, name, d)
    }
    pub fn find_int(&self, name: String) -> Vec<i32> {
        let mut values: Vec<i32> = Vec::new();
        for v in &self.ints {
            if v.name == name {
                let n_values = v.n_values;
                // v.looked_up = true;
                for i in 0..n_values {
                    values.push(v.values[i]);
                }
            }
        }
        values
    }
    pub fn find_float(&self, name: String) -> Vec<Float> {
        let mut values: Vec<Float> = Vec::new();
        for v in &self.floats {
            if v.name == name {
                let n_values = v.n_values;
                // v.looked_up = true;
                for i in 0..n_values {
                    values.push(v.values[i]);
                }
            }
        }
        values
    }
    pub fn find_point2f(&self, name: String) -> Vec<Point2f> {
        let mut values: Vec<Point2f> = Vec::new();
        for v in &self.point2fs {
            if v.name == name {
                let n_values = v.n_values;
                // v.looked_up = true;
                for i in 0..n_values {
                    values.push(v.values[i]);
                }
            }
        }
        values
    }
    pub fn find_vector2f(&self, name: String) -> Vec<Vector2f> {
        let mut values: Vec<Vector2f> = Vec::new();
        for v in &self.vector2fs {
            if v.name == name {
                let n_values = v.n_values;
                // v.looked_up = true;
                for i in 0..n_values {
                    values.push(v.values[i]);
                }
            }
        }
        values
    }
    pub fn find_point3f(&self, name: String) -> Vec<Point3f> {
        let mut values: Vec<Point3f> = Vec::new();
        for v in &self.point3fs {
            if v.name == name {
                let n_values = v.n_values;
                // v.looked_up = true;
                for i in 0..n_values {
                    values.push(v.values[i]);
                }
            }
        }
        values
    }
    pub fn find_vector3f(&self, name: String) -> Vec<Vector3f> {
        let mut values: Vec<Vector3f> = Vec::new();
        for v in &self.vector3fs {
            if v.name == name {
                let n_values = v.n_values;
                // v.looked_up = true;
                for i in 0..n_values {
                    values.push(v.values[i]);
                }
            }
        }
        values
    }
    pub fn find_normal3f(&self, name: String) -> Vec<Normal3f> {
        let mut values: Vec<Normal3f> = Vec::new();
        for v in &self.normals {
            if v.name == name {
                let n_values = v.n_values;
                // v.looked_up = true;
                for i in 0..n_values {
                    values.push(v.values[i]);
                }
            }
        }
        values
    }
}

#[derive(Default)]
pub struct TextureParams {
    pub float_textures: HashMap<String, Arc<Texture<Float> + Send + Sync>>,
    pub spectrum_textures: HashMap<String, Arc<Texture<Spectrum> + Send + Sync>>,
    pub geom_params: ParamSet,
    pub material_params: ParamSet,
}

impl TextureParams {
    pub fn get_spectrum_texture(&mut self, n: String, def: Spectrum) -> Arc<Texture<Spectrum> + Send + Sync>
    {
        let mut name: String = self.geom_params.find_texture(n.clone());
        if name == String::new() {
            name = self.material_params.find_texture(n.clone());
        }
        if name != String::new() {
            match self.spectrum_textures.get(name.as_str()) {
                Some(spectrum_texture) => {
                    return spectrum_texture.clone();
                },
                None => {
                    panic!("Couldn't find spectrum texture named \"{}\" for parameter \"{}\"",
                           name, n);
                },
            }
        }
        let mut val: Spectrum = self.material_params.find_one_spectrum(n.clone(), def);
        val = self.geom_params.find_one_spectrum(n.clone(), val);
        Arc::new(ConstantTexture { value: val })
    }
    pub fn get_float_texture(&mut self, n: String, def: Float) -> Arc<Texture<Float> + Send + Sync>
    {
        let mut name: String = self.geom_params.find_texture(n.clone());
        if name == String::new() {
            name = self.material_params.find_texture(n.clone());
        }
        if name != String::new() {
            match self.float_textures.get(name.as_str()) {
                Some(_float_texture) => {
                },
                None => {
                    panic!("Couldn't find float texture named \"{}\" for parameter \"{}\"",
                           name, n);
                },
            }
        }
        let mut val: Float = self.material_params.find_one_float(n.clone(), def);
        val = self.geom_params.find_one_float(n.clone(), val);
        Arc::new(ConstantTexture { value: val })
    }
    pub fn get_float_texture_or_null(&mut self, n: String) -> Option<Arc<Texture<Float> + Send + Sync>>
    {
        let mut name: String = self.geom_params.find_texture(n.clone());
        if name == String::new() {
            name = self.material_params.find_texture(n.clone());
        }
        if name != String::new() {
            match self.float_textures.get(name.as_str()) {
                Some(_float_texture) => {
                },
                None => {
                    println!("Couldn't find float texture named \"{}\" for parameter \"{}\"",
                             name, n);
                    return None;
                },
            }
        }
        let mut val: Vec<Float> = self.material_params.find_float(n.clone());
        if val.len() == 0_usize {
            val = self.geom_params.find_float(n.clone());
        }
        if val.len() == 0_usize {
            None
        } else {
            Some(Arc::new(ConstantTexture { value: val[0] }))
        }
    }
    pub fn find_float(&mut self, name: String, d: Float) -> Float {
        self.geom_params.find_one_float(name.clone(),
                                        self.material_params.find_one_float(name.clone(),
                                                                              d))
    }
    pub fn find_string(&mut self, name: String, d: String) -> String {
        self.geom_params.find_one_string(name.clone(),
                                         self.material_params.find_one_string(name.clone(),
                                                                              d))
    }
    pub fn find_filename(&mut self, name: String, d: String) -> String {
        self.geom_params.find_one_filename(name.clone(),
                                           self.material_params.find_one_filename(name.clone(),
                                                                                  d))
    }
    pub fn find_int(&mut self, name: String, d: i32) -> i32 {
        self.geom_params.find_one_int(name.clone(),
                                      self.material_params.find_one_int(name.clone(),
                                                                        d))
    }
    pub fn find_bool(&mut self, name: String, d: bool) -> bool {
        self.geom_params.find_one_bool(name.clone(),
                                       self.material_params.find_one_bool(name.clone(),
                                                                          d))
    }
    pub fn find_vector3f(&mut self, name: String, d: Vector3f) -> Vector3f {
        self.geom_params.find_one_vector3f(name.clone(),
                                           self.material_params.find_one_vector3f(name.clone(),
                                                                                  d))
    }
}

/// Replaces a macro on the C++ side.
pub fn lookup_one<T>(vec: &Vec<ParamSetItem<T>>, name: String, d: T) -> T
    where T: Clone
{
    for v in vec {
        if v.name == name && v.n_values == 1_usize {
            // v.looked_up = true;
            return v.values[0].clone();
        }
    }
    d
}

// see api.cpp

#[derive(Debug,Default,Copy,Clone)]
pub struct TransformSet {
    pub t: [Transform; 2],
}

pub struct RenderOptions {
    pub transform_start_time: Float,
    pub transform_end_time: Float,
    pub filter_name: String, // "box"
    pub filter_params: ParamSet,
    pub film_name: String, // "image"
    pub film_params: ParamSet,
    pub sampler_name: String, // "halton";
    pub sampler_params: ParamSet,
    pub accelerator_name: String, // "bvh";
    pub accelerator_params: ParamSet,
    pub integrator_name: String, // "path";
    pub integrator_params: ParamSet,
    pub camera_name: String, // "perspective";
    pub camera_params: ParamSet,
    pub camera_to_world: TransformSet,
    // TODO: std::map<std::string, std::shared_ptr<Medium>> namedMedia;
    pub lights: Vec<Arc<Light + Sync + Send>>,
    pub primitives: Vec<Arc<Primitive + Sync + Send>>,
    // TODO: std::map<std::string, std::vector<std::shared_ptr<Primitive>>> instances;
    // TODO: std::vector<std::shared_ptr<Primitive>> *currentInstance = nullptr;
    pub have_scattering_media: bool, // false
}

impl Default for RenderOptions {
    fn default() -> RenderOptions {
        RenderOptions {
            transform_start_time: 0.0 as Float,
            transform_end_time: 1.0 as Float,
            filter_name: String::from("box"),
            filter_params: ParamSet::default(),
            film_name: String::from("image"),
            film_params: ParamSet::default(),
            sampler_name: String::from("halton"),
            sampler_params: ParamSet::default(),
            accelerator_name: String::from("bvh"),
            accelerator_params: ParamSet::default(),
            integrator_name: String::from("image"),
            integrator_params: ParamSet::default(),
            camera_name: String::from("perspective"),
            camera_params: ParamSet::default(),
            camera_to_world: TransformSet {
                t: [Transform {
                    m: Matrix4x4 {
                        m: [[1.0, 0.0, 0.0, 0.0],
                            [0.0, 1.0, 0.0, 0.0],
                            [0.0, 0.0, 1.0, 0.0],
                            [0.0, 0.0, 0.0, 1.0]],
                    },
                    m_inv: Matrix4x4 {
                        m: [[1.0, 0.0, 0.0, 0.0],
                            [0.0, 1.0, 0.0, 0.0],
                            [0.0, 0.0, 1.0, 0.0],
                            [0.0, 0.0, 0.0, 1.0]],
                    },
                }; 2],
            },
            lights: Vec::new(),
            primitives: Vec::new(),
            have_scattering_media: false,
        }
    }
}

impl RenderOptions {
    // pub fn make_integrator(&self) -> Integrator {
    // }
    // pub fn make_camera(&self) -> Camera {
    // }
}

#[derive(Default)]
pub struct GraphicsState {
    // std::string currentInsideMedium, currentOutsideMedium;
    pub float_textures: HashMap<String, Arc<Texture<Float> + Send + Sync>>,
    pub spectrum_textures: HashMap<String, Arc<Texture<Spectrum> + Send + Sync>>,
    pub material_params: ParamSet,
    pub material: String,
    pub named_materials: HashMap<String, Arc<Material + Send + Sync>>,
    pub current_named_material: String,
    pub area_light_params: ParamSet,
    pub area_light: String,
    // bool reverseOrientation = false;
}

// see rng.h

const ONE_MINUS_EPSILON: Float = 0.99999994;
const PCG32_DEFAULT_STATE: u64 = 0x853c49e6748fea9b;
const PCG32_DEFAULT_STREAM: u64 = 0xda3e39cb94b95bdb;
const PCG32_MULT: u64 = 0x5851f42d4c957f2d;

/// Random number generator
#[derive(Debug,Default,Copy,Clone)]
pub struct Rng {
    state: u64,
    inc: u64,
}

impl Rng {
    pub fn new() -> Self {
        Rng {
            state: PCG32_DEFAULT_STATE,
            inc: PCG32_DEFAULT_STREAM,
        }
    }
    pub fn set_sequence(&mut self, initseq: u64) {
        self.state = 0_u64;
        let (shl, _overflow) = initseq.overflowing_shl(1);
        self.inc = shl | 1;
        self.uniform_uint32();
        let (add, _overflow) = self.state.overflowing_add(PCG32_DEFAULT_STATE);
        self.state = add;
        self.uniform_uint32();
    }
    pub fn uniform_uint32(&mut self) -> u32 {
        let oldstate: u64 = self.state;
        let (mul, _overflow) = oldstate.overflowing_mul(PCG32_MULT);
        let (add, _overflow) = mul.overflowing_add(self.inc);
        self.state = add;
        let (shr, _overflow) = oldstate.overflowing_shr(18);
        let combine = shr ^ oldstate;
        let (shr, _overflow) = combine.overflowing_shr(27);
        let xorshifted: u32 = shr as u32;
        let (shr, _overflow) = oldstate.overflowing_shr(59);
        let rot: u32 = shr as u32;
        // bitwise not in Rust is ! (not the ~ operator like in C)
        let (shr, _overflow) = xorshifted.overflowing_shr(rot);
        let (neg, _overflow) = rot.overflowing_neg();
        let (shl, _overflow) = xorshifted.overflowing_shl(neg & 31);
        shr | shl
    }
    pub fn uniform_uint32_bounded(&mut self, b: u32) -> u32 {
        // bitwise not in Rust is ! (not the ~ operator like in C)
        let threshold = (!b + 1) & b;
        loop {
            let r = self.uniform_uint32();
            if r >= threshold {
                return r % b;
            }
        }
    }
    pub fn uniform_float(&mut self) -> Float {
        (self.uniform_uint32() as Float * 2.3283064365386963e-10 as Float).min(ONE_MINUS_EPSILON)
    }
}

// see material.h

/// Is used to inform non-symetric BSDFs about the transported
/// quantity so that they can correctly switch between the adjoint and
/// non-adjoint forms.
#[derive(PartialEq)]
pub enum TransportMode {
    Radiance,
    Importance,
}

// see memory.h

const LOG_BLOCK_SIZE: usize = 3;
const BLOCK_SIZE: usize = 1 << LOG_BLOCK_SIZE;

fn round_up(x: usize) -> usize {
    (x + BLOCK_SIZE - 1) & !(BLOCK_SIZE - 1)
}

#[derive(Debug,Clone,Default)]
pub struct BlockedArray<T> {
    pub data: Vec<T>,
    pub u_res: usize,
    pub v_res: usize,
    pub u_blocks: usize,
    log_block_size: usize,
    block_size: usize,
}

impl<T> BlockedArray<T>
    where T: num::Zero + Clone + Add<T, Output = T>
{
    pub fn new(u_res: usize, v_res: usize) -> BlockedArray<T> {
        let data = vec![num::Zero::zero(); round_up(u_res) * round_up(v_res)];
        BlockedArray {
            u_res: u_res,
            v_res: v_res,
            u_blocks: round_up(u_res) >> LOG_BLOCK_SIZE,
            log_block_size: LOG_BLOCK_SIZE,
            block_size: BLOCK_SIZE,
            data: data,
        }
    }
    pub fn new_from(u_res: usize, v_res: usize, d: &[T]) -> BlockedArray<T> {
        let mut ba = Self::new(u_res, v_res);
        for u in 0..u_res {
            for v in 0..v_res {
                ba[(u, v)] = d[v * u_res + u].clone();
            }
        }
        ba
    }
    pub fn u_size(&self) -> usize {
        self.u_res
    }
    pub fn v_size(&self) -> usize {
        self.v_res
    }
    pub fn block_size(&self) -> usize {
        self.block_size
    }
    pub fn block(&self, a: usize) -> usize {
        a >> self.log_block_size
    }
    pub fn offset(&self, a: usize) -> usize {
        a & (self.block_size() - 1)
    }
}

impl<T> Index<(usize, usize)> for BlockedArray<T>
    where T: num::Zero + std::clone::Clone + Add<T, Output = T>
{
    type Output = T;
    fn index(&self, i: (usize, usize)) -> &T {
        let (u, v) = i;
        let bu = self.block(u);
        let bv = self.block(v);
        let ou = self.offset(u);
        let ov = self.offset(v);
        let offset = self.block_size() * self.block_size() * (self.u_blocks * bv + bu) +
                     self.block_size() * ov + ou;
        &self.data[offset]
    }
}

impl<T> IndexMut<(usize, usize)> for BlockedArray<T>
    where T: num::Zero + std::clone::Clone + Add<T, Output = T>
{
    fn index_mut(&mut self, i: (usize, usize)) -> &mut T {
        let (u, v) = i;
        let bu = self.block(u);
        let bv = self.block(v);
        let ou = self.offset(u);
        let ov = self.offset(v);
        let offset = self.block_size() * self.block_size() * (self.u_blocks * bv + bu) +
                     self.block_size() * ov + ou;
        &mut self.data[offset]
    }
}

// see github/tray_rust/src/sampler/block_queue.rs

/// The queue of blocks to be worked on shared immutably between worker threads.
pub struct BlockQueue {
    /// The block indices of blocks to work on for the image
    blocks: Vec<(u32, u32)>,
    /// Get the dimensions of an individual block
    dimensions: (u32, u32),
    /// Index of the next block to be worked on
    next: AtomicUsize,
}

impl BlockQueue {
    /// Create a block queue for the image with dimensions `img`.
    /// Panics if the image is not evenly broken into blocks of dimension `dim`
    pub fn new(img: (u32, u32), dim: (u32, u32), select_blocks: (usize, usize)) -> BlockQueue {
        if img.0 % dim.0 != 0 || img.1 % dim.1 != 0 {
            panic!("Image with dimension {:?} not evenly divided by blocks of {:?}", img, dim);
        }
        let num_blocks = (img.0 / dim.0, img.1 / dim.1);
        // TODO: the .. operator precedence is very low so we need this paren here at the moment
        // once (hopefully) it's raised we can remove the parens
        let mut blocks: Vec<(u32, u32)> = (0..num_blocks.0 * num_blocks.1)
            .map(|i| (i % num_blocks.0, i / num_blocks.0)).collect();
        blocks.sort_by(|a, b| morton2(a).cmp(&morton2(b)));
        // If we're only rendering a subset of the blocks then filter our list down
        if select_blocks.1 > 0 {
            blocks = blocks.into_iter().skip(select_blocks.0).take(select_blocks.1).collect();
        }
        if blocks.is_empty() {
            println!("Warning: This block queue is empty!");
        }
        BlockQueue { blocks: blocks, dimensions: dim, next: AtomicUsize::new(0) }
    }
    /// Get the dimensions of an individual block in the queue
    pub fn block_dim(&self) -> (u32, u32) { self.dimensions }
    /// Get an iterator to work through the queue
    pub fn iter(&self) -> BlockQueueIterator { BlockQueueIterator { queue: self } }
    /// Get the next block in the queue or None if the queue is finished
    fn next(&self) -> Option<(u32, u32)> {
        let i = self.next.fetch_add(1, Ordering::AcqRel);
        if i >= self.blocks.len() {
            None
        } else {
            Some(self.blocks[i])
        }
    }
    /// Get the length of the queue
    pub fn len(&self) -> usize { self.blocks.len() }
    /// Check if the queue is empty
    pub fn is_empty(&self) -> bool {
        self.next.load(Ordering::AcqRel) >= self.blocks.len()
    }
}

/// Iterator to work through the queue safely
pub struct BlockQueueIterator<'a> {
    queue: &'a BlockQueue,
}

impl<'a> Iterator for BlockQueueIterator<'a> {
    type Item = (u32, u32);
    fn next(&mut self) -> Option<(u32, u32)> {
        self.queue.next()
    }
}

// see github/tray_rust/src/sampler/morton.rs

///! Provides utilities for 2D Morton code generation using Fabian Giesen's Morton
///! code decoding functions, see [his post on Morton codes](https://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/)

/// Insert a 0 bit between each of the low 16 bits of x
pub fn part1_by1(mut x: u32) -> u32 {
	// x = ---- ---- ---- ---- fedc ba98 7654 3210
	x &= 0x0000ffff;
	// x = ---- ---- fedc ba98 ---- ---- 7654 3210
	x = (x ^ (x << 8)) & 0x00ff00ff;
	// x = ---- fedc ---- ba98 ---- 7654 ---- 3210
	x = (x ^ (x << 4)) & 0x0f0f0f0f;
	// x = --fe --dc --ba --98 --76 --54 --32 --10
	x = (x ^ (x << 2)) & 0x33333333;
	// x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
	(x ^ (x << 1)) & 0x55555555
}
/// Compute the Morton code for the `(x, y)` position.
pub fn morton2(p: &(u32, u32)) -> u32 {
	(part1_by1(p.1) << 1) + part1_by1(p.0)
}

/// **Main function** to **render** a scene mutli-threaded (using all
/// available cores).
pub fn render(scene: &Scene,
              camera: Arc<Camera + Send + Sync>,
              mut sampler: &mut ZeroTwoSequenceSampler,
              mut integrator: &mut Arc<SamplerIntegrator + Send + Sync>) {
    // SamplerIntegrator::Render (integrator.cpp)
    let film = camera.get_film();
    let sample_bounds: Bounds2i = film.get_sample_bounds();
    println!("sample_bounds = {:?}", sample_bounds);
    // let mut integrator: DirectLightingIntegrator =
    //     DirectLightingIntegrator::new(LightStrategy::UniformSampleAll, 10, sample_bounds);
    // create and preprocess sampler
    let integrator_option = Arc::get_mut(&mut integrator);
    let integrator: &mut (SamplerIntegrator + Send + Sync) = integrator_option.unwrap();
    integrator.preprocess(scene, &mut sampler);
    // use camera below
    let sample_extent: Vector2i = sample_bounds.diagonal();
    println!("sample_extent = {:?}", sample_extent);
    let tile_size: i32 = 16;
    let x: i32 = (sample_extent.x + tile_size - 1) / tile_size;
    let y: i32 = (sample_extent.y + tile_size - 1) / tile_size;
    let n_tiles: Point2i = Point2i { x: x, y: y };
    println!("n_tiles = {:?}", n_tiles);
    // TODO: ProgressReporter reporter(nTiles.x * nTiles.y, "Rendering");
    println!("Rendering");
    let num_cores: usize = num_cpus::get();
    // DEBUG: let num_cores: usize = 1; // TMP
    {
        let block_queue = BlockQueue::new(((n_tiles.x * tile_size) as u32,
                                           (n_tiles.y * tile_size) as u32),
                                          (tile_size as u32, tile_size as u32),
                                          (0, 0));
        println!("block_queue.len() = {}", block_queue.len());
        let integrator = &integrator;
        let bq = &block_queue;
        let sampler = &sampler;
        let camera = &camera;
        let film = &film;
        let pixel_bounds = integrator.get_pixel_bounds().clone();
        crossbeam::scope(|scope| {
            let (pixel_tx, pixel_rx) = mpsc::channel();
            // spawn worker threads
            for _ in 0..num_cores {
                let pixel_tx = pixel_tx.clone();
                scope.spawn(move || {
                    while let Some((x, y)) = bq.next() {
                        let tile: Point2i = Point2i {
                            x: x as i32,
                            y: y as i32,
                        };
                        let seed: i32 = tile.y * n_tiles.x + tile.x;
                        // don't use ZeroTwoSequenceSampler::Clone(int seed)
                        let mut tile_sampler = ZeroTwoSequenceSampler::clone(sampler);
                        // adjust the seed here
                        tile_sampler.rng.set_sequence(seed as u64);
                        let x0: i32 = sample_bounds.p_min.x + tile.x * tile_size;
                        let x1: i32 = std::cmp::min(x0 + tile_size, sample_bounds.p_max.x);
                        let y0: i32 = sample_bounds.p_min.y + tile.y * tile_size;
                        let y1: i32 = std::cmp::min(y0 + tile_size, sample_bounds.p_max.y);
                        let tile_bounds: Bounds2i = Bounds2i::new(Point2i { x: x0, y: y0 },
                                                                  Point2i { x: x1, y: y1 });
                        // println!("Starting image tile {:?}", tile_bounds);
                        let mut film_tile = film.get_film_tile(tile_bounds);
                        for pixel in &tile_bounds {
                            tile_sampler.start_pixel(pixel);
                            if !pnt2_inside_exclusive(pixel, pixel_bounds) {
                                continue;
                            }
                            let mut done: bool = false;
                            while !done {
                                // let's use the copy_arena crate instead of pbrt's MemoryArena
                                // let mut arena: Arena = Arena::with_capacity(262144); // 256kB

                                // initialize _CameraSample_ for current sample
                                let camera_sample: CameraSample =
                                    tile_sampler.get_camera_sample(pixel);
                                // generate camera ray for current sample
                                let mut ray: Ray = Ray::default();
                                let ray_weight: Float = camera
                                        .generate_ray_differential(&camera_sample, &mut ray);
                                ray.scale_differentials(1.0 as Float /
                                                        (tile_sampler.samples_per_pixel as Float)
                                    .sqrt());
                                // TODO: ++nCameraRays;
                                // evaluate radiance along camera ray
                                let mut l: Spectrum = Spectrum::new(0.0 as Float);
                                let y: Float = l.y();
                                if ray_weight > 0.0 {
                                    l = integrator.li(&mut ray,
                                                      scene,
                                                      &mut tile_sampler, // &mut arena,
                                                      0_i32);
                                }
                                if l.has_nans() {
                                    println!("Not-a-number radiance value returned for pixel \
                                              ({:?}, {:?}), sample {:?}. Setting to black.",
                                             pixel.x,
                                             pixel.y,
                                             tile_sampler.get_current_sample_number());
                                    l = Spectrum::new(0.0);
                                } else if y < -10.0e-5 as Float {
                                    println!("Negative luminance value, {:?}, returned for pixel \
                                              ({:?}, {:?}), sample {:?}. Setting to black.",
                                             y,
                                             pixel.x,
                                             pixel.y,
                                             tile_sampler.get_current_sample_number());
                                    l = Spectrum::new(0.0);
                                } else if y.is_infinite() {
                                    println!("Infinite luminance value returned for pixel ({:?}, \
                                              {:?}), sample {:?}. Setting to black.",
                                             pixel.x,
                                             pixel.y,
                                             tile_sampler.get_current_sample_number());
                                    l = Spectrum::new(0.0);
                                }
                                // println!("Camera sample: {:?} -> ray: {:?} -> L = {:?}",
                                //          camera_sample, ray, l);
                                // add camera ray's contribution to image
                                film_tile.add_sample(camera_sample.p_film, &mut l, ray_weight);
                                done = !tile_sampler.start_next_sample();
                            } // arena is dropped here !
                        }
                        // send the tile through the channel to main thread
                        pixel_tx.send(film_tile).expect(&format!("Failed to send tile"));
                    }
                });
            }
            // spawn thread to collect pixels and render image to file
            scope.spawn(move || {
                for _ in pbr::PbIter::new(0..bq.len()) {
                    let film_tile = pixel_rx.recv().unwrap();
                    // merge image tile into _Film_
                    film.merge_film_tile(&film_tile);
                }
            });
        });
    }
    println!("Rendering finished");
    film.write_image(1.0 as Float);
}
