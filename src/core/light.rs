//! In order for objects in a scene to be visible, there must be a
//! source of illumination so that some light is reflected from them
//! to the camera sensor.

// std
use std::sync::Arc;
// pbrt
use crate::core::geometry::{Normal3f, Point2f, Ray, Vector3f};
use crate::core::interaction::{Interaction, InteractionCommon, SurfaceInteraction};
use crate::core::medium::MediumInterface;
use crate::core::pbrt::{Float, Spectrum};
use crate::core::sampler::Sampler;
use crate::core::scene::Scene;
use crate::lights::diffuse::DiffuseAreaLight;
use crate::lights::distant::DistantLight;
use crate::lights::goniometric::GonioPhotometricLight;
use crate::lights::infinite::InfiniteAreaLight;
use crate::lights::point::PointLight;
// use crate::lights::projection::ProjectionLight;
// use crate::lights::spot::SpotLight;

// see light.h

#[repr(u8)]
pub enum LightFlags {
    DeltaPosition = 1,
    DeltaDirection = 2,
    Area = 4,
    Infinite = 8,
}

pub enum Light {
    DiffuseArea(Box<DiffuseAreaLight>),
    Distant(Box<DistantLight>),
    GonioPhotometric(Box<GonioPhotometricLight>),
    InfiniteArea(Box<InfiniteAreaLight>),
    Point(Box<PointLight>),
    // Projection(Box<ProjectionLight>),
    // Spot(Box<SpotLight>),
}

impl Light {
    /// Returns the radiance arriving at a point at a certain time due
    /// to the light, assuming there are no occluding objects between
    /// them.
    pub fn sample_li(
        &self,
        iref: &InteractionCommon,
        u: Point2f,
        wi: &mut Vector3f,
        pdf: &mut Float,
    ) -> (Spectrum, Option<VisibilityTester>) {
        match self {
            Light::DiffuseArea(light) => light.sample_li(iref, u, wi, pdf),
            Light::Distant(light) => light.sample_li(iref, u, wi, pdf),
            Light::GonioPhotometric(light) => light.sample_li(iref, u, wi, pdf),
            Light::InfiniteArea(light) => light.sample_li(iref, u, wi, pdf),
            Light::Point(light) => light.sample_li(iref, u, wi, pdf),
            // Light::Projection(light) => light.sample_li(iref, u, wi, pdf),
            // Light::Spot(light) => light.sample_li(iref, u, wi, pdf),
        }
    }
    pub fn power(&self) -> Spectrum {
        match self {
            Light::DiffuseArea(light) => light.power(),
            Light::Distant(light) => light.power(),
            Light::GonioPhotometric(light) => light.power(),
            Light::InfiniteArea(light) => light.power(),
            Light::Point(light) => light.power(),
            // Light::Projection(light) => light.power(),
            // Light::Spot(light) => light.power(),
        }
    }
    pub fn preprocess(&self, scene: &Scene) {
        match self {
            Light::DiffuseArea(light) => light.preprocess(scene),
            Light::Distant(light) => light.preprocess(scene),
            Light::GonioPhotometric(light) => light.preprocess(scene),
            Light::InfiniteArea(light) => light.preprocess(scene),
            Light::Point(light) => light.preprocess(scene),
            // Light::Projection(light) => light.preprocess(scene),
            // Light::Spot(light) => light.preprocess(scene),
        }
    }
    pub fn le(&self, ray: &mut Ray) -> Spectrum {
        match self {
            Light::DiffuseArea(light) => light.le(ray),
            Light::Distant(light) => light.le(ray),
            Light::GonioPhotometric(light) => light.le(ray),
            Light::InfiniteArea(light) => light.le(ray),
            Light::Point(light) => light.le(ray),
            // Light::Projection(light) => light.le(ray),
            // Light::Spot(light) => light.le(ray),
        }
    }
    pub fn pdf_li(&self, iref: &dyn Interaction, wi: Vector3f) -> Float {
        match self {
            Light::DiffuseArea(light) => light.pdf_li(iref, wi),
            Light::Distant(light) => light.pdf_li(iref, wi),
            Light::GonioPhotometric(light) => light.pdf_li(iref, wi),
            Light::InfiniteArea(light) => light.pdf_li(iref, wi),
            Light::Point(light) => light.pdf_li(iref, wi),
            // Light::Projection(light) => light.pdf_li(iref, wi),
            // Light::Spot(light) => light.pdf_li(iref, wi),
        }
    }
    pub fn sample_le(
        &self,
        u1: Point2f,
        u2: Point2f,
        time: Float,
        ray: &mut Ray,
        n_light: &mut Normal3f,
        pdf_pos: &mut Float,
        pdf_dir: &mut Float,
    ) -> Spectrum {
        match self {
            Light::DiffuseArea(light) => {
                light.sample_le(u1, u2, time, ray, n_light, pdf_pos, pdf_dir)
            }
            Light::Distant(light) => light.sample_le(u1, u2, time, ray, n_light, pdf_pos, pdf_dir),
            Light::GonioPhotometric(light) => {
                light.sample_le(u1, u2, time, ray, n_light, pdf_pos, pdf_dir)
            }
            Light::InfiniteArea(light) => {
                light.sample_le(u1, u2, time, ray, n_light, pdf_pos, pdf_dir)
            }
            Light::Point(light) => light.sample_le(u1, u2, time, ray, n_light, pdf_pos, pdf_dir),
            // Light::Projection(light) => {
            //     light.sample_le(u1, u2, time, ray, n_light, pdf_pos, pdf_dir)
            // }
            // Light::Spot(light) => light.sample_le(u1, u2, time, ray, n_light, pdf_pos, pdf_dir),
        }
    }
    pub fn pdf_le(&self, ray: &Ray, n_light: &Normal3f, pdf_pos: &mut Float, pdf_dir: &mut Float) {
        match self {
            Light::DiffuseArea(light) => light.pdf_le(ray, n_light, pdf_pos, pdf_dir),
            Light::Distant(light) => light.pdf_le(ray, n_light, pdf_pos, pdf_dir),
            Light::GonioPhotometric(light) => light.pdf_le(ray, n_light, pdf_pos, pdf_dir),
            Light::InfiniteArea(light) => light.pdf_le(ray, n_light, pdf_pos, pdf_dir),
            Light::Point(light) => light.pdf_le(ray, n_light, pdf_pos, pdf_dir),
            // Light::Projection(light) => light.pdf_le(ray, n_light, pdf_pos, pdf_dir),
            // Light::Spot(light) => light.pdf_le(ray, n_light, pdf_pos, pdf_dir),
        }
    }
    pub fn get_flags(&self) -> u8 {
        match self {
            Light::DiffuseArea(light) => light.get_flags(),
            Light::Distant(light) => light.get_flags(),
            Light::GonioPhotometric(light) => light.get_flags(),
            Light::InfiniteArea(light) => light.get_flags(),
            Light::Point(light) => light.get_flags(),
            // Light::Projection(light) => light.get_flags(),
            // Light::Spot(light) => light.get_flags(),
        }
    }
    pub fn get_n_samples(&self) -> i32 {
        match self {
            Light::DiffuseArea(light) => light.get_n_samples(),
            Light::Distant(light) => light.get_n_samples(),
            Light::GonioPhotometric(light) => light.get_n_samples(),
            Light::InfiniteArea(light) => light.get_n_samples(),
            Light::Point(light) => light.get_n_samples(),
            // Light::Projection(light) => light.get_n_samples(),
            // Light::Spot(light) => light.get_n_samples(),
        }
    }
    // AreaLight
    pub fn l(&self, intr: &InteractionCommon, w: &Vector3f) -> Spectrum {
        match self {
            Light::DiffuseArea(light) => light.l(intr, w),
            _ => panic!("Not an area light"),
        }
    }
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

/// VisibilityTesters are created by providing two Interaction
/// objects, one for each end point of the shadow ray to be traced.
#[derive(Default, Clone)]
pub struct VisibilityTester {
    pub p0: InteractionCommon, // TODO: private
    pub p1: InteractionCommon, // TODO: private
}

impl VisibilityTester {
    pub fn unoccluded(&self, scene: &Scene) -> bool {
        !scene.intersect_p(&mut self.p0.spawn_ray_to(&self.p1))
    }
    pub fn tr(&self, scene: &Scene, sampler: &mut Sampler) -> Spectrum {
        let mut ray: Ray = self.p0.spawn_ray_to(&self.p1);
        let mut tr: Spectrum = Spectrum::new(1.0 as Float);
        loop {
            let mut it: InteractionCommon = InteractionCommon::default();
            let mut medium_interface: Option<Arc<MediumInterface>> = None;
            let mut isect: SurfaceInteraction = SurfaceInteraction::default();
            if scene.intersect(&mut ray, &mut isect) {
                // handle opaque surface along ray's path
                if let Some(primitive_raw) = isect.primitive {
                    let primitive = unsafe { &*primitive_raw };
                    if let Some(_material) = primitive.get_material() {
                        return Spectrum::default();
                    } else {
                        // update transmittance for current ray segment
                        if let Some(ref medium_arc) = ray.medium {
                            tr *= medium_arc.tr(&ray, sampler);
                        }
                    }
                }
                if let Some(mi_arc) = &isect.common.medium_interface {
                    medium_interface = Some(mi_arc.clone());
                }
                it.p = isect.common.p;
                it.time = isect.common.time;
                it.p_error = isect.common.p_error;
                it.wo = isect.common.wo;
                it.n = isect.common.n;
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

// Area lights are light sources defined by one or more **Shapes**
// that emit light from their surface, with some directional
// distribution of radiance at each point on the surface.

// pub trait AreaLight: Light {
//     fn l(&self, intr: &InteractionCommon, w: &Vector3f) -> Spectrum;
// }
