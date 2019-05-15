//! Implementations of the **Medium** base class provide various
//! representations of volumetric scattering properties in a region of
//! space.
//!
//! - GridDensityMedium
//! - HomogeneousMedium
//!
//! ## Grid Density Medium
//!
//! ![Smoke from a CFD Simulation](/doc/img/smoke_plume_pbrt_rust_volpath.png)
//!
//! ## Homogeneous Medium
//!
//! ![A Volumetric Caustic](/doc/img/volume_caustic_pbrt_rust_mlt.png)

pub mod grid;
pub mod homogeneous;
