//! # Shapes
//!
//! Careful abstraction of geometric shapes in a ray tracer is a key
//! component of a clean system design, and shapes are the ideal
//! candidate for an object-oriented approach. All geometric
//! primitives implement a common interface, and the rest of the
//! renderer can use this interface without needing any details about
//! the underlying shape. This makes it possible to separate the
//! geometric and the shading subsystem of pbrt.
//!
//! ## Spheres
//!
//! Spheres are a special case of a general type of surfaces called
//! quadrics. They are the simplest type of curved surfaces that is
//! useful to a ray tracer and are a good starting point for general
//! ray intersection routines.
//!
//! ```rust
//! extern crate pbrt;
//! 
//! use pbrt::core::pbrt::Float;
//! use pbrt::core::transform::Transform;
//! use pbrt::geometry::Vector3f;
//! use pbrt::shapes::sphere::Sphere;
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
//! ## Triangle Meshes
//!
//! While a natural representation would be to have a **Triangle**
//! shape implementation where each triangle stored the positions of
//! its three vertices, a more memory-efficient representation is to
//! separately store entire triangle meshes with an array of vertex
//! positions where each individual triangle just stores three offsets
//! into this array for its three vertices.
//!
//! ```rust
//! extern crate pbrt;
//! 
//! use pbrt::core::transform::Transform;
//! use pbrt::geometry::{Normal3f, Point2f, Point3f, Vector3f};
//! use pbrt::shapes::triangle::{Triangle, TriangleMesh};
//! use std::sync::Arc;
//! 
//! fn main() {
//!     let translate: Transform = Transform::translate(Vector3f {
//!                                                         x: 0.25,
//!                                                         y: 0.0,
//!                                                         z: 0.0,
//!                                                     });
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
//! ## Cones
//!
//! TODO
//!
//! ## Disks
//!
//! The disk is an interesting quadric since it has a particularly
//! straightforward intersection routine that avoids solving the
//! quadric equation.
//!
//! ```rust
//! extern crate pbrt;
//! 
//! use pbrt::core::pbrt::Float;
//! use pbrt::core::transform::Transform;
//! use pbrt::geometry::Vector3f;
//! use pbrt::shapes::disk::Disk;
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
//! ## Cylinders
//!
//! Another useful quadric is the cylinder. Cylinder shapes are
//! centered around the z axis.
//!
//! ```rust
//! extern crate pbrt;
//! 
//! use pbrt::core::pbrt::Float;
//! use pbrt::core::transform::Transform;
//! use pbrt::geometry::Vector3f;
//! use pbrt::shapes::cylinder::Cylinder;
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
//! ## Hyperboloids
//!
//! TODO
//!
//! ## Paraboloids
//!
//! TODO
//!

// pbrt
use core::pbrt::Float;
use core::interaction::{Interaction, InteractionCommon, SurfaceInteraction};
use geometry::{Bounds3f, Point2f, Ray, Vector3f};

pub mod cylinder;
pub mod disk;
pub mod plymesh;
pub mod sphere;
pub mod triangle;

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
    fn sample_with_ref_point(&self,
                             iref: &InteractionCommon,
                             u: Point2f,
                             pdf: &mut Float)
                             -> InteractionCommon;
    fn pdf(&self, iref: &Interaction, wi: Vector3f) -> Float;
}
