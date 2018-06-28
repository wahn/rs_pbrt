//! # pbrt
//!
//! [Rust][rust] crate to implement at least parts of the [PBRT
//! book][book]'s C++ code. You can find a copy of the current code
//! [here][repo].
//!
//! [rust]: https://www.rust-lang.org/en-US
//! [book]: http://www.pbrt.org
//! [repo]: https://github.com/wahn/rs_pbrt
//!

#[macro_use]
extern crate hexf;
extern crate atomic;
#[cfg(feature = "openexr")]
extern crate half;
extern crate image;
extern crate num;
#[cfg(feature = "openexr")]
extern crate openexr;
extern crate ply_rs;
extern crate rayon;
extern crate time;
extern crate typed_arena;

pub mod accelerators;
pub mod blockqueue;
pub mod cameras;
pub mod core;
pub mod filters;
pub mod integrators;
pub mod lights;
pub mod materials;
pub mod media;
pub mod samplers;
pub mod shapes;
pub mod textures;
