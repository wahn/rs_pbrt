//! # pbrt
//!
//! [Rust][rust] crate to implement at least parts of the [PBRT
//! book][book]'s C++ code. You can find a copy of the current code
//! [here][repo].
//!
//! The main render loop for integrators implementing the
//! `SamplerIntegrator` trait can be found [here].
//!
//! There are two more render loops:
//!
//! 1. [render_bdpt][render_bdpt] for bidirectional path tracing
//! 2. [render_mlt][render_mlt] for Metropolis Light Transport
//!
//! [rust]: https://www.rust-lang.org/en-US
//! [book]: http://www.pbrt.org
//! [repo]: https://github.com/wahn/rs_pbrt
//! [here]: https://www.rs-pbrt.org/doc/crates/pbrt/integrators/fn.render.html
//! [render_bdpt]: https://www.rs-pbrt.org/doc/crates/pbrt/integrators/bdpt/fn.render_bdpt.html
//! [render_mlt]: https://www.rs-pbrt.org/doc/crates/pbrt/integrators/mlt/fn.render_mlt.html

#[macro_use]
extern crate hexf;
extern crate atomic;
extern crate byteorder;
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
