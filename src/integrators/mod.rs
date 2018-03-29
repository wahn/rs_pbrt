//! **Integrator** is an abstract base class that defines the
//! **render()** method that must be provided by all integrators.
//!
//! - AOIntegrator
//! - BDPTIntegrator
//! - DirectLightingIntegrator
//! - MLTIntegrator
//! - PathIntegrator
//! - SPPMIntegrator
//! - VolPathIntegrator
//! - WhittedIntegrator
//!
//! ## Ambient Occlusion (AO)
//!
//! Ambient Occlusion is most often calculated by casting rays in
//! every direction from the surface. Rays which reach the background
//! or sky increase the brightness of the surface, whereas a ray which
//! hits any other object contributes no illumination. As a result,
//! points surrounded by a large amount of geometry are rendered dark,
//! whereas points with little geometry on the visible hemisphere
//! appear light.
//!
//! ![Ambient Occlusion](https://www.janwalter.org/assets/cornell_box_pbrt_rust_ao.png)
//!
//! ## Direct Lighting
//!
//! The **DirectLightingIntegrator** accounts only for direct lighting
//! &mdash; light that has traveled directly from a light source to the
//! point being shaded &mdash; and ignores indirect illumination from
//! objects that are not themselfes emissive, except for basic
//! specular reflection and transmission effects.
//!
//! ![Direct Lighting](https://www.janwalter.org/assets/cornell_box_v0.1.13.png)
//!
//! ## Path Tracing
//!
//! Path tracing incrementally generates paths of scattering events
//! starting at the camera and ending at light sources in the scene.
//!
//! ![Path Tracing](https://www.janwalter.org/assets/cornell_box_v0.2.0_high.png)

//! ## Bidirectional Path Tracing (BDPT)
//!
//! Bidirectional path tracing is a generalization of the standard
//! pathtracing algorithm that can be much more efficient. It
//! constructs paths that start from the camera on one end, from the
//! light on the other end, and connects in the middle with a
//! visibility ray.
//!
//! ![Bidirectional Path
//! Tracing](https://www.janwalter.org/assets/art_gallery_pbrt_rust_07n.png)

pub mod ao;
pub mod bdpt;
pub mod directlighting;
pub mod mlt;
pub mod path;
