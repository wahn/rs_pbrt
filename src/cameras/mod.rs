//! The abstract **Camera** base class holds generic camera options
//! and defines the interface that all camera implementations must
//! provide.
//!
//! - EnvironmentCamera
//! - OrthographicCamera
//! - PerspectiveCamera
//! - RealisticCamera
//!
//! ## Perspective Camera
//!
//! The two most characteristic features of perspective are that
//! objects are smaller as their distance from the observer increases;
//! and that they are subject to **foreshortening**, meaning that an
//! object's dimensions along the line of sight are shorter than its
//! dimensions across the line of sight.
//!
//! ![Perspective Camera](/doc/img/wikipedia_perspective_camera.png)

pub mod perspective;
