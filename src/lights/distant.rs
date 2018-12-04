// std
use std;
use std::f32::consts::PI;
use std::sync::RwLock;
// pbrt
use core::geometry::{vec3_coordinate_system, vec3_normalize};
use core::geometry::{Bounds3f, Normal3f, Point2f, Point3f, Ray, Vector3f};
use core::interaction::{Interaction, InteractionCommon};
use core::light::{Light, LightFlags, VisibilityTester};
use core::medium::MediumInterface;
use core::pbrt::{Float, Spectrum};
use core::sampling::concentric_sample_disk;
use core::scene::Scene;
use core::transform::Transform;

// see distant.h

pub struct DistantLight {
    // private data (see distant.h)
    pub l: Spectrum,
    pub w_light: Vector3f,
    pub world_center: RwLock<Point3f>,
    pub world_radius: RwLock<Float>,
    // inherited from class Light (see light.h)
    pub flags: u8,
    pub n_samples: i32,
    pub medium_interface: MediumInterface,
    pub light_to_world: Transform,
    pub world_to_light: Transform,
}

impl DistantLight {
    pub fn new(light_to_world: &Transform, l: &Spectrum, w_light: &Vector3f) -> Self {
        DistantLight {
            l: *l,
            w_light: vec3_normalize(&light_to_world.transform_vector(&*w_light)),
            world_center: RwLock::new(Point3f::default()),
            world_radius: RwLock::new(0.0),
            flags: LightFlags::DeltaDirection as u8,
            n_samples: 1_i32,
            medium_interface: MediumInterface::default(),
            light_to_world: Transform::default(),
            world_to_light: Transform::default(),
        }
    }
}

impl Light for DistantLight {
    fn sample_li(
        &self,
        iref: &InteractionCommon,
        _u: &Point2f,
        wi: &mut Vector3f,
        pdf: &mut Float,
        vis: &mut VisibilityTester,
    ) -> Spectrum {
        // TODO: ProfilePhase _(Prof::LightSample);
        *wi = self.w_light;
        *pdf = 1.0 as Float;
        let p_outside: Point3f =
            iref.p + self.w_light * (2.0 as Float * *self.world_radius.read().unwrap());
        *vis = VisibilityTester {
            p0: InteractionCommon {
                p: iref.p,
                time: iref.time,
                p_error: iref.p_error,
                wo: iref.wo,
                n: iref.n,
                medium_interface: None,
            },
            p1: InteractionCommon {
                p: p_outside,
                time: iref.time,
                p_error: Vector3f::default(),
                wo: Vector3f::default(),
                n: Normal3f::default(),
                medium_interface: None,
            },
        };
        self.l
    }
    fn power(&self) -> Spectrum {
        let world_radius: Float = *self.world_radius.read().unwrap();
        self.l * PI * world_radius * world_radius
    }
    /// Some of the **DistanceLight** methods need to know the bounds
    /// of the scene. Because lights are created before the scene
    /// geometry, these bounds aren't available when the
    /// **DistanceLight** constructor runs. Therefore,
    /// **DistanceLight** implements the optional *preprocess()*
    /// method to get the bound. This method is called at the end of
    /// the **Scene** constructor.
    fn preprocess(&self, scene: &Scene) {
        let mut world_center_ref = self.world_center.write().unwrap();
        let mut world_radius_ref = self.world_radius.write().unwrap();
        Bounds3f::bounding_sphere(
            &scene.world_bound(),
            &mut world_center_ref,
            &mut world_radius_ref,
        );
    }
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

        // choose point on disk oriented toward infinite light direction
        let mut v1: Vector3f = Vector3f::default();
        let mut v2: Vector3f = Vector3f::default();
        vec3_coordinate_system(&self.w_light, &mut v1, &mut v2);
        let cd: Point2f = concentric_sample_disk(u1);
        let world_center: Point3f = *self.world_center.read().unwrap();
        let world_radius: Float = *self.world_radius.read().unwrap();
        let p_disk: Point3f = world_center + (v1 * cd.x + v2 * cd.y) * world_radius;
        // set ray origin and direction for infinite light ray
        *ray = Ray {
            o: p_disk + self.w_light * world_radius,
            d: -self.w_light,
            t_max: std::f32::INFINITY,
            time: time,
            differential: None,
            medium: None,
        };
        *n_light = Normal3f::from(ray.d);
        *pdf_pos = 1.0 as Float / (PI * world_radius * world_radius);
        *pdf_dir = 1.0 as Float;
        self.l
    }
    fn pdf_le(&self, _ray: &Ray, _n_light: &Normal3f, pdf_pos: &mut Float, pdf_dir: &mut Float) {
        let world_radius: Float = *self.world_radius.read().unwrap();
        *pdf_pos = 1.0 as Float / (PI * world_radius * world_radius);
        *pdf_dir = 0.0 as Float;
    }
    fn get_flags(&self) -> u8 {
        self.flags
    }
    fn get_n_samples(&self) -> i32 {
        self.n_samples
    }
}
