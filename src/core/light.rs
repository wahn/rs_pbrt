//! In order for objects in a scene to be visible, there must be a
//! source of illumination so that some light is reflected from them
//! to the camera sensor.

// std
use std::sync::Arc;
// pbrt
use core::geometry::{Normal3f, Point2f, Ray, Vector3f};
use core::interaction::{Interaction, InteractionCommon};
use core::medium::MediumInterface;
use core::pbrt::{Float, Spectrum};
use core::primitive::Primitive;
use core::sampler::Sampler;
use core::scene::Scene;

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
    fn sample_li(
        &self,
        iref: &InteractionCommon,
        u: &Point2f,
        wi: &mut Vector3f,
        pdf: &mut Float,
        vis: &mut VisibilityTester,
    ) -> Spectrum;
    fn power(&self) -> Spectrum;
    fn preprocess(&self, scene: &Scene);
    fn le(&self, _ray: &mut Ray) -> Spectrum;
    fn pdf_li(&self, iref: &Interaction, wi: Vector3f) -> Float;
    fn sample_le(
        &self,
        u1: &Point2f,
        u2: &Point2f,
        time: Float,
        ray: &mut Ray,
        n_light: &mut Normal3f,
        pdf_pos: &mut Float,
        pdf_dir: &mut Float,
    ) -> Spectrum;
    fn pdf_le(&self, ray: &Ray, n_light: &Normal3f, pdf_pos: &mut Float, pdf_dir: &mut Float);
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
#[derive(Default, Clone)]
pub struct VisibilityTester {
    pub p0: InteractionCommon, // TODO: private
    pub p1: InteractionCommon, // TODO: private
}

impl VisibilityTester {
    pub fn unoccluded(&self, scene: &Scene) -> bool {
        !scene.intersect_p(&mut self.p0.spawn_ray_to(&self.p1))
    }
    pub fn tr(&self, scene: &Scene, sampler: &mut Box<Sampler + Send + Sync>) -> Spectrum {
        let mut ray: Ray = self.p0.spawn_ray_to(&self.p1);
        let mut tr: Spectrum = Spectrum::new(1.0 as Float);
        loop {
            let mut it: InteractionCommon = InteractionCommon::default();
            let mut medium_interface: Option<Arc<MediumInterface>> = None;
            if let Some(isect) = scene.intersect(&mut ray) {
                // handle opaque surface along ray's path
                if let Some(primitive) = isect.primitive {
                    if let Some(_material) = primitive.get_material() {
                        return Spectrum::default();
                    } else {
                        // update transmittance for current ray segment
                        if let Some(ref medium_arc) = ray.medium {
                            tr *= medium_arc.tr(&ray, sampler);
                        }
                    }
                }
                if let Some(mi_arc) = isect.medium_interface {
                    medium_interface = Some(mi_arc.clone());
                }
                it.p = isect.p;
                it.time = isect.time;
                it.p_error = isect.p_error;
                it.wo = isect.wo;
                it.n = isect.n;
                it.medium_interface = medium_interface;
            } else {
                // update transmittance for current ray segment
                if let Some(ref medium_arc) = ray.medium {
                    tr *= medium_arc.tr(&ray, sampler);
                }
                break;
            }
            ray = it.spawn_ray_to(&self.p1);
        }
        tr
    }
}

/// Area lights are light sources defined by one or more **Shapes**
/// that emit light from their surface, with some directional
/// distribution of radiance at each point on the surface.
pub trait AreaLight: Light {
    fn l(&self, intr: &InteractionCommon, w: &Vector3f) -> Spectrum;
}
