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
//! ![HairMaterial](https://www.janwalter.org/assets/hair_rust_pbrt.png)
//!
//! ## SubstrateMaterial
//!
//! ![SubstrateMaterial](https://www.janwalter.org/assets/ganesha.png)

pub mod fourier;
pub mod glass;
pub mod hair;
pub mod matte;
pub mod metal;
pub mod mirror;
pub mod mixmat;
pub mod plastic;
pub mod substrate;
pub mod uber;
