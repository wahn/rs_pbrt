// pbrt
use core::interaction::{Interaction, InteractionCommon};
use core::light::{Light, LightFlags, VisibilityTester};
use core::pbrt::{Float, Spectrum};
use core::scene::Scene;
use core::transform::Transform;
use geometry::{Normal3f, Point2f, Point3f, Ray, Vector3f};
use geometry::{pnt3_distance_squared, vec3_normalize};

// see point.h

#[derive(Debug,Copy,Clone)]
pub struct PointLight {
    // private data (see point.h)
    pub p_light: Point3f,
    pub i: Spectrum,
    // inherited from class Light (see light.h)
    flags: u8,
    n_samples: i32,
}

impl PointLight {
    pub fn new(light_to_world: &Transform, i: &Spectrum) -> Self {
        PointLight {
            p_light: light_to_world.transform_point(Point3f::default()),
            i: *i,
            flags: LightFlags::DeltaPosition as u8,
            n_samples: 1_i32,
        }
    }
}

impl Light for PointLight {
    fn sample_li(&self,
                 iref: &InteractionCommon,
                 _u: Point2f,
                 wi: &mut Vector3f,
                 pdf: &mut Float,
                 vis: &mut VisibilityTester)
                 -> Spectrum {
        // TODO: ProfilePhase _(Prof::LightSample);
        *wi = vec3_normalize(self.p_light - iref.p);
        *pdf = 1.0 as Float;
        *vis = VisibilityTester {
            p0: InteractionCommon {
                p: iref.p,
                time: iref.time,
                p_error: iref.p_error,
                wo: iref.wo,
                n: iref.n,
            },
            p1: InteractionCommon {
                p: self.p_light,
                time: iref.time,
                p_error: Vector3f::default(),
                wo: Vector3f::default(),
                n: Normal3f::default(),
            },
        };
        self.i / pnt3_distance_squared(self.p_light, iref.p)
    }
    fn power(&self) -> Spectrum {
        Spectrum::default()
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
    fn get_flags(&self) -> u8 {
        self.flags
    }
    fn get_n_samples(&self) -> i32 {
        self.n_samples
    }
}
