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
//!
//! ## Orthographic Camera
//!
//! The orthographic camera is based on the orthographic projection
//! transformation. The orthographic transformation takes a
//! rectangular region of the scene and projects it onto the front
//! face of the box that defines the region. It doesn’t give the
//! effect of foreshortening - objects becoming smaller on the image
//! plane as they get farther away - but it does leave parallel lines
//! parallel, and it preserves relative distance between objects.
//!
//! ![Orthographic Camera](/doc/img/wikipedia_orthographic_camera.png)
//!
//! ## Environment Camera
//!
//! A camera model that traces rays in all directions around a point
//! in the scene, giving a 2D view of everything that is visible from
//! that point. One important use of this image representation is
//! environment lighting - a rendering technique that uses image-based
//! representations of light in a scene.

pub mod environment;
pub mod orthographic;
pub mod perspective;
