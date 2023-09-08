// std
use std::cell::Cell;
use std::f32::consts::PI;
use std::sync::Arc;
// pbrt
use crate::core::geometry::pnt3_distance_squaredf;
use crate::core::geometry::{Normal3f, Point2f, Point3f, Ray, Vector3f};
use crate::core::interaction::{Interaction, InteractionCommon};
use crate::core::light::{LightFlags, VisibilityTester};
use crate::core::medium::{Medium, MediumInterface};
use crate::core::pbrt::{Float, Spectrum};
use crate::core::sampling::{uniform_sample_sphere, uniform_sphere_pdf};
use crate::core::scene::Scene;
use crate::core::transform::Transform;

// see point.h

#[derive(Clone)]
pub struct PointLight {
    // private data (see point.h)
    pub p_light: Point3f,
    pub i: Spectrum,
    // inherited from class Light (see light.h)
    pub flags: u8,
    pub n_samples: i32,
    pub medium_interface: MediumInterface,
}

impl PointLight {
    pub fn new(
        light_to_world: &Transform,
        medium_interface: &MediumInterface,
        i: &Spectrum,
    ) -> Self {
        let mut inside: Option<Arc<Medium>> = None;
        let mut outside: Option<Arc<Medium>> = None;
        if let Some(ref mi_inside) = medium_interface.inside {
            inside = Some(mi_inside.clone());
        }
        if let Some(ref mi_outside) = medium_interface.outside {
            outside = Some(mi_outside.clone());
        }
        PointLight {
            p_light: light_to_world.transform_point(&Point3f::default()),
            i: *i,
            flags: LightFlags::DeltaPosition as u8,
            n_samples: 1_i32,
            medium_interface: MediumInterface { inside, outside },
        }
    }
    // Light
    pub fn sample_li<'a, 'b>(
        &'b self,
        iref: &'a InteractionCommon,
        light_intr: &'b mut InteractionCommon,
        _u: Point2f,
        wi: &mut Vector3f,
        pdf: &mut Float,
        vis: &mut VisibilityTester<'a, 'b>,
    ) -> Spectrum {
        // TODO: ProfilePhase _(Prof::LightSample);
        *wi = (self.p_light - iref.p).normalize();
        *pdf = 1.0 as Float;
        light_intr.p = self.p_light;
        light_intr.time = iref.time;
        vis.p0 = Some(iref);
        vis.p1 = Some(light_intr);
        self.i / pnt3_distance_squaredf(&self.p_light, &iref.p)
    }
    pub fn power(&self) -> Spectrum {
        self.i * (4.0 as Float * PI)
    }
    pub fn preprocess(&self, _scene: &Scene) {}
    /// Default implementation returns no emitted radiance for a ray
    /// that escapes the scene bounds.
    pub fn le(&self, _ray: &Ray) -> Spectrum {
        Spectrum::new(0.0 as Float)
    }
    pub fn pdf_li(&self, _iref: &dyn Interaction, _wi: &Vector3f) -> Float {
        0.0 as Float
    }
    pub fn sample_le(
        &self,
        u1: Point2f,
        _u2: Point2f,
        time: Float,
        ray: &mut Ray,
        n_light: &mut Normal3f,
        pdf_pos: &mut Float,
        pdf_dir: &mut Float,
    ) -> Spectrum {
        // TODO: ProfilePhase _(Prof::LightSample);
        *ray = Ray {
            o: self.p_light,
            d: uniform_sample_sphere(u1),
            t_max: Cell::new(std::f32::INFINITY),
            time,
            differential: None,
            medium: None,
        };
        *n_light = Normal3f::from(ray.d);
        *pdf_pos = 1.0 as Float;
        *pdf_dir = uniform_sphere_pdf();
        self.i
    }
    pub fn get_flags(&self) -> u8 {
        self.flags
    }
    pub fn get_n_samples(&self) -> i32 {
        self.n_samples
    }
    pub fn pdf_le(
        &self,
        _ray: &Ray,
        _n_light: &Normal3f,
        pdf_pos: &mut Float,
        pdf_dir: &mut Float,
    ) {
        *pdf_pos = 0.0 as Float;
        *pdf_dir = uniform_sphere_pdf();
    }
}
