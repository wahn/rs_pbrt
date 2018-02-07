//! Acceleration structures are one of the components at the heart of
//! any ray tracer. Without algorithms to reduce the number of
//! unnecessary ray intersection tests, tracing a single ray through a
//! scene would take time linear in the number of primitives in the
//! scene, since the ray would need to be tested against each
//! primitive in turn to find the closest intersection.
//!
//! - BVHAccel
//! - KdTreeAccel

pub mod bvh;
