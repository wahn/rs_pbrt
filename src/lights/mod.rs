//! In order for objects in a scene to be visible, there must be a
//! source of illumination so that some light is reflected from them
//! to the camera sensor.
//!
//! - DiffuseAreaLight
//! - DistantLight
//! - GonioPhotometricLight
//! - InfiniteAreaLight
//! - PointLight
//! - ProjectionLight
//! - SpotLight
//!
//! ## Diffuse Area Lights
//!
//! Area lights are light sources defined by one or more **Shapes**
//! that emit light from their surface, with some directional
//! distribution of radiance at each point on the surface.
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
//! use pbrt::core::pbrt::{Float, Spectrum};
//! use pbrt::core::shape::Shape;
//! use pbrt::core::transform::Transform;
//! use pbrt::core::geometry::Vector3f;
//! use pbrt::lights::diffuse::DiffuseAreaLight;
//! use pbrt::shapes::disk::Disk;
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
//! ## Distant Lights
//!
//! A distant light, also known as directional light, describes an
//! emitter that deposits illumination from the same direction at
//! every point in space.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::core::pbrt::Spectrum;
//! use pbrt::core::transform::Transform;
//! use pbrt::core::geometry::{Point3f, Vector3f};
//! use pbrt::lights::distant::DistantLight;
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
//! ## Goniophotometric Diagram Lights
//!
//! TODO
//!
//! ## Infinite Area Lights
//!
//! Area lights are light sources defined by one or more **Shapes**
//! that emit light from their surface, with some directional
//! distribution of radiance at each point on the surface.
//!
//! An infinite far away area light source that surrounds the entire
//! scene. One way to visualize this light is as an enormous sphere
//! that casts light into the scene from every direction.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::core::pbrt::Spectrum;
//! use pbrt::core::transform::Transform;
//! use pbrt::lights::infinite::InfiniteAreaLight;
//!
//! fn main() {
//!     let l: Spectrum = Spectrum::new(50.0);
//!     let light_to_world: Transform = Transform::default();
//!     let n_samples: i32 = 1;
//!     let texmap: String = String::from("");
//!     let infinite_light: InfiniteAreaLight =
//!         InfiniteAreaLight::new(&light_to_world, &l, n_samples, texmap);
//! }
//! ```
//!
//! ## Point Lights
//!
//! Isotropic point light source that emits the same amount of light
//! in all directions.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::core::pbrt::Spectrum;
//! use pbrt::core::transform::Transform;
//! use pbrt::lights::point::PointLight;
//!
//! fn main() {
//!     let i: Spectrum = Spectrum::new(50.0);
//!     let light_to_world: Transform = Transform::default();
//!     let point_light: PointLight = PointLight::new(&light_to_world, &i);
//!     println!("point_light = {:?}", point_light);
//! }
//! ```
//!
//! ## Texture Projection Lights
//!
//! TODO
//!
//! ## Spotlights
//!
//! TODO
//!

pub mod diffuse;
pub mod distant;
pub mod infinite;
pub mod point;
