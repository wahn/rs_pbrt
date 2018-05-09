// std
use std;
use std::f32::consts::PI;
use std::sync::Arc;
// pbrt
use core::geometry::{Normal3f, Point2f, Point3f, Ray, Vector3f};
use core::geometry::{pnt3_distance_squared, vec3_normalize};
use core::interaction::{Interaction, InteractionCommon};
use core::light::{Light, LightFlags, VisibilityTester};
use core::medium::{Medium, MediumInterface};
use core::pbrt::radians;
use core::pbrt::{Float, Spectrum};
use core::reflection::cos_theta;
use core::sampling::{uniform_cone_pdf, uniform_sample_cone};
use core::scene::Scene;
use core::transform::Transform;

// see spot.h

#[derive(Clone)]
pub struct SpotLight {
    // private data (see spot.h)
    pub p_light: Point3f,
    pub i: Spectrum,
    pub cos_total_width: Float,
    pub cos_falloff_start: Float,
    // inherited from class Light (see light.h)
    pub flags: u8,
    pub n_samples: i32,
    pub medium_interface: MediumInterface,
    pub light_to_world: Transform,
    pub world_to_light: Transform,
}

impl SpotLight {
    pub fn new(
        light_to_world: &Transform,
        medium_interface: &MediumInterface,
        i: &Spectrum,
        total_width: Float,
        falloff_start: Float,
    ) -> Self {
        let mut inside: Option<Arc<Medium + Send + Sync>> = None;
        let mut outside: Option<Arc<Medium + Send + Sync>> = None;
        if let Some(ref mi_outside) = medium_interface.outside {
            // in C++: MediumInterface(const Medium *medium) : inside(medium), outside(medium)
            inside = Some(mi_outside.clone());
            outside = Some(mi_outside.clone());
        }
        SpotLight {
            p_light: light_to_world.transform_point(&Point3f::default()),
            i: *i,
            cos_total_width: radians(total_width).cos(),
            cos_falloff_start: radians(falloff_start).cos(),
            flags: LightFlags::DeltaPosition as u8,
            n_samples: 1_i32,
            medium_interface: MediumInterface {
                inside: inside,
                outside: outside,
            },
            light_to_world: *light_to_world,
            world_to_light: Transform::inverse(light_to_world),
        }
    }
    pub fn falloff(&self, w: &Vector3f) -> Float {
        let wl: Vector3f = vec3_normalize(&self.world_to_light.transform_vector(w));
        let cos_theta: Float = wl.z;
        if cos_theta < self.cos_total_width {
            return 0.0 as Float;
        }
        if cos_theta >= self.cos_falloff_start {
            return 1.0 as Float;
        }
        // compute falloff inside spotlight cone
        let delta: Float =
            (cos_theta - self.cos_total_width) / (self.cos_falloff_start - self.cos_total_width);
        (delta * delta) * (delta * delta)
    }
}

impl Light for SpotLight {
    fn sample_li(
        &self,
        iref: &InteractionCommon,
        _u: &Point2f,
        wi: &mut Vector3f,
        pdf: &mut Float,
        vis: &mut VisibilityTester,
    ) -> Spectrum {
        // TODO: ProfilePhase _(Prof::LightSample);
        *wi = vec3_normalize(&(self.p_light - iref.p));
        *pdf = 1.0 as Float;
        let mut inside: Option<Arc<Medium + Send + Sync>> = None;
        let mut outside: Option<Arc<Medium + Send + Sync>> = None;
        if let Some(ref mi_inside_arc) = self.medium_interface.inside {
            inside = Some(mi_inside_arc.clone());
        }
        if let Some(ref mi_outside_arc) = self.medium_interface.outside {
            outside = Some(mi_outside_arc.clone());
        }
        let medium_interface_arc: Arc<MediumInterface> = Arc::new(MediumInterface::new(inside, outside));
        *vis = VisibilityTester {
            p0: InteractionCommon {
                p: iref.p,
                time: iref.time,
                p_error: iref.p_error,
                wo: iref.wo,
                n: iref.n,
                medium_interface: Some(medium_interface_arc.clone()),
            },
            p1: InteractionCommon {
                p: self.p_light,
                time: iref.time,
                p_error: Vector3f::default(),
                wo: Vector3f::default(),
                n: Normal3f::default(),
                medium_interface: Some(medium_interface_arc.clone()),
            },
        };
        self.i * self.falloff(&-*wi) / pnt3_distance_squared(&self.p_light, &iref.p)
    }
    fn power(&self) -> Spectrum {
        self.i * 2.0 as Float * PI
            * (1.0 as Float - 0.5 as Float * (self.cos_falloff_start + self.cos_total_width))
    }
    fn preprocess(&self, _scene: &Scene) {}
    /// Default implementation returns no emitted radiance for a ray
    /// that escapes the scene bounds.
    fn le(&self, _ray: &mut Ray) -> Spectrum {
        Spectrum::new(0.0 as Float)
    }
    fn pdf_li(&self, _iref: &Interaction, _wi: Vector3f) -> Float {
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
        // TODO: ProfilePhase _(Prof::LightSample);
        let w: Vector3f = uniform_sample_cone(u1, self.cos_total_width);
        let mut inside: Option<Arc<Medium + Send + Sync>> = None;
        if let Some(ref mi_inside) = self.medium_interface.inside {
            inside = Some(mi_inside.clone());
        }
        *ray = Ray {
            o: self.p_light,
            d: self.light_to_world.transform_vector(&w),
            t_max: std::f32::INFINITY,
            time: time,
            differential: None,
            medium: inside,
        };
        *n_light = Normal3f::from(ray.d);
        *pdf_pos = 1.0 as Float;
        *pdf_dir = uniform_cone_pdf(self.cos_total_width);
        self.i * self.falloff(&ray.d)
    }
    fn get_flags(&self) -> u8 {
        self.flags
    }
    fn get_n_samples(&self) -> i32 {
        self.n_samples
    }
    fn pdf_le(&self, ray: &Ray, _n_light: &Normal3f, pdf_pos: &mut Float, pdf_dir: &mut Float) {
        *pdf_pos = 0.0 as Float;
        if cos_theta(&self.world_to_light.transform_vector(&ray.d)) > self.cos_total_width {
            *pdf_dir = uniform_cone_pdf(self.cos_total_width);
        } else {
            *pdf_dir = 0.0 as Float;
        }
    }
}
