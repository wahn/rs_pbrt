// std
use std;
use std::f32::consts::PI;
use std::sync::Arc;
// pbrt
use crate::core::geometry::pnt3_distance_squared;
use crate::core::geometry::{Bounds2f, Normal3f, Point2f, Point3f, Ray, Vector3f};
use crate::core::interaction::{Interaction, InteractionCommon};
use crate::core::light::{Light, LightFlags, VisibilityTester};
use crate::core::medium::{Medium, MediumInterface};
use crate::core::pbrt::{Float, Spectrum};
use crate::core::sampling::{uniform_sample_sphere, uniform_sphere_pdf};
use crate::core::scene::Scene;
use crate::core::transform::Transform;

// see projection.h

pub struct ProjectionLight {
    // private data (see projection.h)
    pub p_light: Point3f,
    pub i: Spectrum,
    pub light_projection: Transform,
    pub hither: Float,
    pub yon: Float,
    pub screen_bounds: Bounds2f,
    pub cos_total_width: Float,
    // inherited from class Light (see light.h)
    pub flags: u8,
    pub n_samples: i32,
    pub medium_interface: MediumInterface,
}

impl ProjectionLight {}

impl Light for ProjectionLight {
    fn sample_li(
        &self,
        iref: &InteractionCommon,
        _u: &Point2f,
        wi: &mut Vector3f,
        pdf: &mut Float,
        vis: &mut VisibilityTester,
    ) -> Spectrum {
        // WORK
        Spectrum::default()
    }
    fn power(&self) -> Spectrum {
        // WORK
        Spectrum::default()
    }
    fn preprocess(&self, _scene: &Scene) {}
    /// Default implementation returns no emitted radiance for a ray
    /// that escapes the scene bounds.
    fn le(&self, _ray: &mut Ray) -> Spectrum {
        Spectrum::new(0.0 as Float)
    }
    fn pdf_li(&self, _iref: &dyn Interaction, _wi: Vector3f) -> Float {
        0.0 as Float
    }
    fn sample_le(
        &self,
        u1: &Point2f,
        _u2: &Point2f,
        time: Float,
        ray: &mut Ray,
        n_light: &mut Normal3f,
        pdf_pos: &mut Float,
        pdf_dir: &mut Float,
    ) -> Spectrum {
        // WORK
        Spectrum::default()
    }
    fn get_flags(&self) -> u8 {
        self.flags
    }
    fn get_n_samples(&self) -> i32 {
        self.n_samples
    }
    fn pdf_le(&self, _ray: &Ray, _n_light: &Normal3f, pdf_pos: &mut Float, pdf_dir: &mut Float) {
        *pdf_pos = 0.0 as Float;
        *pdf_dir = uniform_sphere_pdf();
    }
}
