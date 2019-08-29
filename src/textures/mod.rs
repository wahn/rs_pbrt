//! **Texture** is a template class parameterized by return type of
//! its evaluation function. This design makes it possible to reuse
//! almost all of the code among textures that return different
//! types. PBRT currently uses only **Float** and **Spectrum**
//! textures.
//!
//! - BilerpTexture
//! - Checkerboard2DTexture
//! - ConstantTexture
//! - DotsTexture
//! - FBmTexture
//! - ImageTexture
//! - MarbleTexture
//! - MixTexture
//! - PtexTexture
//! - ScaleTexture
//! - UVTexture
//! - WindyTexture
//! - WrinkledTexture
//!
//! ## Checkerboard2DTexture
//!
//! ![Checkerboard2DTexture](/doc/img/checkerboard_pbrt_rust.png)
//!
//! ## DotsTexture
//!
//! ![DotsTexture](/doc/img/dots_pbrt_rust.png)
//!
//! ## FBmTexture
//!
//! ![FBmTexture](/doc/img/fbm_pbrt_rust.png)
//!
//! ## MarbleTexture
//!
//! ![MarbleTexture](/doc/img/marble_pbrt_rust.png)
//!
//! ## WindyTexture
//!
//! ![WindyTexture](/doc/img/windy_pbrt_rust.png)
//!
//! ## WrinkledTexture
//!
//! ![WrinkledTexture](/doc/img/wrinkled_pbrt_rust.png)

pub mod checkerboard;
pub mod constant;
pub mod dots;
pub mod fbm;
pub mod imagemap;
pub mod marble;
pub mod mix;
pub mod scale;
pub mod windy;
pub mod wrinkled;
