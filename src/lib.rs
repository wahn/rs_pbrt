//! # pbrt
//!
//! [Rust][rust] crate to implement at least parts of the [PBRT
//! book][book]'s C++ code. You can find a copy of the current code
//! [here][repo].
//!
//! The main render loop for integrators implementing the
//! `SamplerIntegrator` can be found [here].
//!
//! There are three more render loops:
//!
//! 1. [render loop][render_bdpt] for bidirectional path tracing
//! 2. [render loop][render_mlt] for Metropolis Light Transport
//! 2. [render loop][render_sppm] for Stochastic Progressive Photon Mapping
//!
//! [rust]: https://www.rust-lang.org
//! [book]: http://www.pbrt.org
//! [repo]: https://github.com/wahn/rs_pbrt
//! [here]: core/integrator/enum.SamplerIntegrator.html#method.render
//! [render_bdpt]: integrators/bdpt/struct.BDPTIntegrator.html#method.render
//! [render_mlt]: integrators/mlt/struct.MLTIntegrator.html#method.render
//! [render_sppm]: integrators/sppm/struct.SPPMIntegrator.html#method.render

#[macro_use] extern crate impl_ops;

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
