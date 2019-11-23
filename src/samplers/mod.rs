//! The **Sampler** base class not only defines the interface to
//! samplers but also provides some common functionality for use by
//! **Sampler** implementations.
//!
//! - HaltonSampler
//! - MaxMinDistSampler
//! - RandomSampler
//! - SobolSampler
//! - StratifiedSampler
//! - ZeroTwoSequenceSampler
//!
//! ## Halton Sampler
//!
//! The Halton Sampler generates not only points that are guaranteed
//! to not clump too closely together, but also generates points that
//! are simultaneously well distributed over all the dimensions of the
//! sample vector.
//!
//! ![halton](/doc/img/cornell_box_pbrt_rust_halton.png)
//!
//! ## Random Sampler
//!
//! The Random Sampler is using the random number generetor class
//! **RNG** based on an unpublished manuscript by O'Neill: A family of
//! simple fast space-efficient statistically good algorithms for
//! random number generation.
//!
//! ![random](/doc/img/cornell_box_pbrt_rust_random.png)
//!
//! ## Sobol Sampler
//!
//! The Sobol Sampler is very efficient to implement while also being
//! extremly well distributed over all dimensions of the sample
//! vector. The weakness of the Sobol' points is that they are prone
//! to structural grid artefacts before convergence.
//!
//! ![sobol](/doc/img/cornell_box_pbrt_rust_sobol.png)
//!
//! ## (0,2)-Sequence Sampler
//!
//! Certain low-discrepancy sequences allow us to satisfy two
//! desirable properties of samples: they generate sample vectors for
//! a pixel's worth of image samples such that the sample values for
//! each pixel sample are well distributed with respect to each other,
//! and simultaneuously such that the aggregate collection of sample
//! values for all of the pixel samples in the pixel are collectively
//! well distributed.
//!
//! ![lowdiscrepancy](/doc/img/cornell_box_pbrt_rust_lowdiscrepancy.png)
//!

pub mod halton;
pub mod random;
// pub mod sobol;
// pub mod zerotwosequence;
