//! The abstract **Material** class defines the interface that
//! material implementations must provide.
//!
//! - DisneyMaterial
//! - FourierMaterial
//! - GlassMaterial
//! - HairMaterial
//! - KdSubsurfaceMaterial
//! - MatteMaterial
//! - MetalMaterial
//! - MirrorMaterial
//! - MixMaterial
//! - PlasticMaterial
//! - SubstrateMaterial
//! - SubsurfaceMaterial
//! - TranslucentMaterial
//! - UberMaterial
//!
//! ## HairMaterial
//!
//! ![HairMaterial](/doc/img/hair_pbrt_rust.png)
//!
//! ## SubstrateMaterial
//!
//! ![SubstrateMaterial](/doc/img/ganesha_pbrt_rust.png)

pub mod disney;
pub mod fourier;
pub mod glass;
pub mod hair;
pub mod matte;
pub mod metal;
pub mod mirror;
pub mod mixmat;
pub mod plastic;
pub mod substrate;
// pub mod subsurface;
pub mod translucent;
pub mod uber;
