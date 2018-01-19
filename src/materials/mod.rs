//! # Materials
//!
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
//! ## SubstrateMaterial
//!
//! ![SubstrateMaterial](https://www.janwalter.org/assets/ganesha.png)

pub mod glass;
pub mod hair;
pub mod matte;
pub mod metal;
pub mod mirror;
pub mod mixmat;
pub mod plastic;
pub mod substrate;
pub mod uber;
