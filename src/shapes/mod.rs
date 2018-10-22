//! Careful abstraction of geometric shapes in a ray tracer is a key
//! component of a clean system design, and shapes are the ideal
//! candidate for an object-oriented approach. All geometric
//! primitives implement a common interface, and the rest of the
//! renderer can use this interface without needing any details about
//! the underlying shape. This makes it possible to separate the
//! geometric and the shading subsystem of pbrt.
//!
//! - Cone
//! - Curve
//! - Cylinder
//! - Disk
//! - Hyperboloid
//! - Paraboloid
//! - Sphere
//! - Triangle
//!
//! ## Cones
//!
//! TODO
//!
//! ## Curves
//!
//! TODO
//!
//! ## Spheres
//!
//! Spheres are a special case of a general type of surfaces called
//! quadrics. They are the simplest type of curved surfaces that is
//! useful to a ray tracer and are a good starting point for general
//! ray intersection routines.
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
//! ## Disks
//!
//! The disk is an interesting quadric since it has a particularly
//! straightforward intersection routine that avoids solving the
//! quadric equation.
//!
//! ## Cylinders
//!
//! Another useful quadric is the cylinder. Cylinder shapes are
//! centered around the z axis.
//!
//! ## Hyperboloids
//!
//! TODO
//!
//! ## Paraboloids
//!
//! TODO
//!

pub mod curve;
pub mod cylinder;
pub mod disk;
pub mod loopsubdiv;
pub mod nurbs;
pub mod plymesh;
pub mod sphere;
pub mod triangle;
