// pbrt
use core::interaction::{Interaction, InteractionCommon};
use core::pbrt::{Float, Spectrum};
use core::scene::Scene;
use geometry::{Point2f, Ray, Vector3f};

// see light.h

#[repr(u8)]
pub enum LightFlags {
    DeltaPosition = 1,
    DeltaDirection = 2,
    Area = 4,
    Infinite = 8,
}

pub trait Light {
    /// Returns the radiance arriving at a point at a certain time due
    /// to the light, assuming there are no occluding objects between
    /// them.
    fn sample_li(&self,
                 iref: &InteractionCommon,
                 u: Point2f,
                 wi: &mut Vector3f,
                 pdf: &mut Float,
                 vis: &mut VisibilityTester)
                 -> Spectrum;
    fn power(&self) -> Spectrum;
    fn preprocess(&self, scene: &Scene);
    fn le(&self, _ray: &mut Ray) -> Spectrum;
    fn pdf_li(&self, iref: &Interaction, wi: Vector3f) -> Float;
    fn get_flags(&self) -> u8;
    fn get_n_samples(&self) -> i32;
}

/// Check if LightFlags::DeltaPosition or LightFlags::DeltaDirection
/// is set.
pub fn is_delta_light(flags: u8) -> bool {
    let mut pos: bool = false;
    let mut dir: bool = false;
    if (flags & LightFlags::DeltaPosition as u8) > 0 {
        pos = true;
    }
    if (flags & LightFlags::DeltaDirection as u8) > 0 {
        dir = true;
    }
    pos || dir
}

/// A closure - an object that encapsulates a small amount of data and
/// some computation that is yet to be done.
#[derive(Debug,Default,Copy,Clone)]
pub struct VisibilityTester {
    pub p0: InteractionCommon, // TODO: private
    pub p1: InteractionCommon, // TODO: private
}

impl VisibilityTester {
    pub fn unoccluded(&self, scene: &Scene) -> bool {
        !scene.intersect_p(&mut self.p0.spawn_ray_to(self.p1))
    }
}

/// Area lights are light sources defined by one or more **Shapes**
/// that emit light from their surface, with some directional
/// distribution of radiance at each point on the surface.
pub trait AreaLight: Light {
    fn l(&self, intr: &InteractionCommon, w: Vector3f) -> Spectrum;
}
