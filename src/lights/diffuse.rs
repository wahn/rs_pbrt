// std
use std::f32::consts::PI;
use std::rc::Rc;
use std::sync::Arc;
// pbrt
use crate::core::geometry::{nrm_abs_dot_vec3, nrm_dot_vec3, vec3_coordinate_system};
use crate::core::geometry::{Normal3f, Point2f, Ray, Vector3f, XYEnum};
use crate::core::interaction::{Interaction, InteractionCommon};
use crate::core::light::{LightFlags, VisibilityTester};
use crate::core::medium::{Medium, MediumInterface};
use crate::core::pbrt::{Float, Spectrum};
use crate::core::rng::FLOAT_ONE_MINUS_EPSILON;
use crate::core::sampling::{cosine_hemisphere_pdf, cosine_sample_hemisphere};
use crate::core::scene::Scene;
use crate::core::shape::Shape;
use crate::core::transform::Transform;

// see diffuse.h

pub struct DiffuseAreaLight {
    pub l_emit: Spectrum,
    pub shape: Arc<Shape>,
    pub two_sided: bool,
    pub area: Float,
    // inherited from class Light (see light.h)
    pub flags: u8,
    pub n_samples: i32,
    pub medium_interface: MediumInterface,
    // light_to_world: Transform,
    // world_to_light: Transform,
}

impl DiffuseAreaLight {
    pub fn new(
        _light_to_world: &Transform,
        medium_interface: &MediumInterface,
        l_emit: &Spectrum,
        n_samples: i32,
        shape: Arc<Shape>,
        two_sided: bool,
    ) -> Self {
        let area: Float = shape.area();
        let mut inside: Option<Arc<Medium>> = None;
        let mut outside: Option<Arc<Medium>> = None;
        if let Some(ref mi_inside) = medium_interface.inside {
            inside = Some(mi_inside.clone());
        }
        if let Some(ref mi_outside) = medium_interface.outside {
            outside = Some(mi_outside.clone());
        }
        DiffuseAreaLight {
            l_emit: *l_emit,
            shape,
            two_sided,
            area,
            // inherited from class Light (see light.h)
            flags: LightFlags::Area as u8,
            n_samples: std::cmp::max(1_i32, n_samples),
            medium_interface: MediumInterface { inside, outside },
            // light_to_world: *light_to_world,
            // world_to_light: Transform::inverse(*light_to_world),
        }
    }
    // Light
    pub fn sample_li(
        &self,
        iref: Rc<InteractionCommon>,
        u: Point2f,
        wi: &mut Vector3f,
        pdf: &mut Float,
    ) -> (Spectrum, Option<VisibilityTester>) {
        // TODO: ProfilePhase _(Prof::LightSample);
        let p_shape: Rc<InteractionCommon> = self.shape.sample_with_ref_point(iref.clone(), u, pdf);
        // TODO: iref.mediumInterface = mediumInterface;
        if *pdf == 0.0 as Float || (p_shape.p - iref.p).length_squared() == 0.0 as Float {
            *pdf = 0.0 as Float;
            return (Spectrum::default(), None);
        }
        let new_wi: Vector3f = (p_shape.p - iref.p).normalize();
        *wi = new_wi;
        (
            self.l(&p_shape, &-new_wi),
            Some(VisibilityTester {
                p0: Some(Rc::new(InteractionCommon {
                    p: iref.p,
                    time: iref.time,
                    p_error: iref.p_error,
                    wo: iref.wo,
                    n: iref.n,
                    medium_interface: None,
                })),
                p1: Some(p_shape),
            }),
        )
    }
    pub fn power(&self) -> Spectrum {
        // return (twoSided ? 2 : 1) * Lemit * area * Pi;
        let factor = if self.two_sided {
            2.0 as Float
        } else {
            1.0 as Float
        };
        self.l_emit * factor * self.area * PI
    }
    pub fn preprocess(&self, _scene: &Scene) {
        // TODO?
    }
    pub fn le(&self, _ray: &mut Ray) -> Spectrum {
        Spectrum::default()
    }
    pub fn pdf_li(&self, iref: &dyn Interaction, wi: Vector3f) -> Float {
        // TODO: ProfilePhase _(Prof::LightPdf);
        self.shape.pdf_with_ref_point(iref, &wi)
    }
    pub fn sample_le(
        &self,
        u1: Point2f,
        u2: Point2f,
        _time: Float,
        ray: &mut Ray,
        n_light: &mut Normal3f,
        pdf_pos: &mut Float,
        pdf_dir: &mut Float,
    ) -> Spectrum {
        // TODO: ProfilePhase _(Prof::LightSample);

        // sample a point on the area light's _Shape_, _p_shape_
        let ic: Rc<InteractionCommon> = self.shape.sample(u1, pdf_pos);
        // TODO: p_shape.mediumInterface = mediumInterface;
        *n_light = ic.n;
        // sample a cosine-weighted outgoing direction _w_ for area light
        let mut w: Vector3f;
        if self.two_sided {
            let mut u: Point2f = Point2f { x: u2.x, y: u2.y };
            // choose a side to sample and then remap u[0] to [0,1]
            // before applying cosine-weighted hemisphere sampling for
            // the chosen side.
            if u[XYEnum::X] < 0.5 as Float {
                u[XYEnum::X] = (u[XYEnum::X] * 2.0 as Float).min(FLOAT_ONE_MINUS_EPSILON);
                w = cosine_sample_hemisphere(u);
            } else {
                u[XYEnum::X] =
                    ((u[XYEnum::X] - 0.5 as Float) * 2.0 as Float).min(FLOAT_ONE_MINUS_EPSILON);
                w = cosine_sample_hemisphere(u);
                w.z *= -1.0 as Float;
            }
            *pdf_dir = 0.5 as Float * cosine_hemisphere_pdf(w.z.abs());
        } else {
            w = cosine_sample_hemisphere(u2);
            *pdf_dir = cosine_hemisphere_pdf(w.z);
        }
        let n: Vector3f = Vector3f::from(ic.n);
        let mut v1: Vector3f = Vector3f::default();
        let mut v2: Vector3f = Vector3f::default();
        vec3_coordinate_system(&n, &mut v1, &mut v2);
        w = v1 * w.x + v2 * w.y + n * w.z;
        *ray = ic.spawn_ray(&w);
        self.l(&ic, &w)
    }
    pub fn pdf_le(&self, ray: &Ray, n: &Normal3f, pdf_pos: &mut Float, pdf_dir: &mut Float) {
        *pdf_pos = self.shape.pdf(Rc::new(InteractionCommon::default()));
        if self.two_sided {
            *pdf_dir = 0.5 as Float * cosine_hemisphere_pdf(nrm_abs_dot_vec3(&n, &ray.d));
        } else {
            *pdf_dir = cosine_hemisphere_pdf(nrm_dot_vec3(&n, &ray.d));
        }
    }
    pub fn get_flags(&self) -> u8 {
        self.flags
    }
    pub fn get_n_samples(&self) -> i32 {
        self.n_samples
    }
    // AreaLight
    pub fn l(&self, intr: &InteractionCommon, w: &Vector3f) -> Spectrum {
        if self.two_sided || nrm_dot_vec3(&intr.n, &w) > 0.0 as Float {
            self.l_emit
        } else {
            Spectrum::new(0.0 as Float)
        }
    }
}
