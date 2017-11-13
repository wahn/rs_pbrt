//! # Ambient Occlusion (AO)
//!
//! TODO
//!
//! # Direct Lighting
//!
//! The **DirectLightingIntegrator** accounts only for direct lighting
//! &mdash; light that has traveled directly from a light source to the
//! point being shaded &mdash; and ignores indirect illumination from
//! objects that are not themselfes emissive, except for basic
//! specular reflection and transmission effects.
//!
//! # Path Tracing
//!
//! TODO

pub mod ao;
pub mod directlighting;
pub mod path;
