//! In order for objects in a scene to be visible, there must be a
//! source of illumination so that some light is reflected from them
//! to the camera sensor.
//!
//! - DiffuseAreaLight
//! - DistantLight
//! - GonioPhotometricLight
//! - InfiniteAreaLight
//! - PointLight
//! - ProjectionLight
//! - SpotLight
//!
//! ## Diffuse Area Lights
//!
//! Area lights are light sources defined by one or more **Shapes**
//! that emit light from their surface, with some directional
//! distribution of radiance at each point on the surface.
//!
//! **DiffuseAreaLight** implements a basic area light source with a
//! uniform spatial and directional radiance distribution. The surface
//! it emits from is defined by a **Shape**. It only emits light on
//! the side of the surface with outward-facing surface normal; there
//! is no emission from the other side.
//!
//! ## Distant Lights
//!
//! A distant light, also known as directional light, describes an
//! emitter that deposits illumination from the same direction at
//! every point in space.
//!
//! ## Goniophotometric Diagram Lights
//!
//! TODO
//!
//! ## Infinite Area Lights
//!
//! Area lights are light sources defined by one or more **Shapes**
//! that emit light from their surface, with some directional
//! distribution of radiance at each point on the surface.
//!
//! An infinite far away area light source that surrounds the entire
//! scene. One way to visualize this light is as an enormous sphere
//! that casts light into the scene from every direction.
//!
//! ## Point Lights
//!
//! Isotropic point light source that emits the same amount of light
//! in all directions.
//!
//! ## Texture Projection Lights
//!
//! TODO
//!
//! ## Spotlights
//!
//! TODO
//!

pub mod diffuse;
pub mod distant;
pub mod infinite;
pub mod point;
