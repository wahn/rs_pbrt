// std
use std;
use std::f32::consts::PI;
use std::io::BufReader;
use std::sync::Arc;
// pbrt
use crate::core::geometry::{pnt3_distance_squared, spherical_phi, spherical_theta};
use crate::core::geometry::{Normal3f, Point2f, Point2i, Point3f, Ray, Vector3f};
use crate::core::interaction::{Interaction, InteractionCommon};
use crate::core::light::{Light, LightFlags, VisibilityTester};
use crate::core::medium::MediumInterface;
use crate::core::mipmap::{ImageWrap, MipMap};
use crate::core::pbrt::{Float, Spectrum};
use crate::core::pbrt::{INV_2_PI, INV_PI};
use crate::core::sampling::{uniform_sample_sphere, uniform_sphere_pdf};
use crate::core::scene::Scene;
use crate::core::transform::Transform;

// see goniometric.h

#[derive(Clone)]
pub struct GonioPhotometricLight {
    pub p_light: Point3f,
    pub i: Spectrum,
    pub mipmap: Option<Arc<MipMap<Spectrum>>>,
    // inherited from class Light (see light.h)
    pub flags: u8,
    pub n_samples: i32,
    pub medium_interface: MediumInterface,
    pub light_to_world: Transform,
    pub world_to_light: Transform,
}

impl GonioPhotometricLight {
    #[cfg(not(feature = "openexr"))]
    pub fn new(
        light_to_world: &Transform,
        medium_interface: &MediumInterface,
        i: &Spectrum,
        texname: String,
    ) -> Self {
        GonioPhotometricLight::new_hdr(light_to_world, medium_interface, i, texname)
    }
    pub fn new_hdr(
        light_to_world: &Transform,
        _medium_interface: &MediumInterface,
        i: &Spectrum,
        texname: String,
    ) -> Self {
        if texname != String::from("") {
            let file = std::fs::File::open(texname.clone()).unwrap();
            let reader = BufReader::new(file);
            let img_result = image::hdr::HDRDecoder::with_strictness(reader, false);
            if img_result.is_ok() {
                if let Some(hdr) = img_result.ok() {
                    let meta = hdr.metadata();
                    let resolution: Point2i = Point2i {
                        x: meta.width as i32,
                        y: meta.height as i32,
                    };
                    let img_result = hdr.read_image_transform(|p| {
                        let rgb = p.to_hdr();
                        Spectrum::rgb(rgb[0], rgb[1], rgb[2])
                    });
                    if img_result.is_ok() {
                        let texels = img_result.ok().unwrap();
                        // create _MipMap_ from converted texels (see above)
                        let do_trilinear: bool = false;
                        let max_aniso: Float = 8.0 as Float;
                        let wrap_mode: ImageWrap = ImageWrap::Repeat;
                        let projection_map = Arc::new(MipMap::new(
                            &resolution,
                            &texels[..],
                            do_trilinear,
                            max_aniso,
                            wrap_mode,
                        ));
                        let p_light: Point3f = light_to_world.transform_point(&Point3f::default());
                        return GonioPhotometricLight {
                            p_light: p_light,
                            i: *i,
                            mipmap: Some(projection_map),
                            flags: LightFlags::DeltaPosition as u8,
                            n_samples: 1_i32,
                            medium_interface: MediumInterface::default(),
                            light_to_world: *light_to_world,
                            world_to_light: Transform::inverse(&*light_to_world),
                        };
                    }
                }
            } else {
                println!("WARNING: ProjectionLight::new() ... no OpenEXR support !!!");
            }
        }
        GonioPhotometricLight {
            p_light: light_to_world.transform_point(&Point3f::default()),
            i: *i,
            mipmap: None,
            flags: LightFlags::DeltaPosition as u8,
            n_samples: 1_i32,
            medium_interface: MediumInterface::default(),
            light_to_world: Transform::default(),
            world_to_light: Transform::default(),
        }
    }
    pub fn scale(&self, w: &Vector3f) -> Spectrum {
        let mut wp: Vector3f = self.world_to_light.transform_vector(w).normalize();
        std::mem::swap(&mut wp.y, &mut wp.z);
        let theta: Float = spherical_theta(&wp);
        let phi: Float = spherical_phi(&wp);
        if let Some(mipmap) = &self.mipmap {
            let st: Point2f = Point2f {
                x: phi * INV_2_PI,
                y: theta * INV_PI,
            };
            mipmap.lookup_pnt_flt(&st, 0.0 as Float)
        } else {
            Spectrum::new(1.0 as Float)
        }
    }
}

impl Light for GonioPhotometricLight {
    fn sample_li(
        &self,
        iref: &InteractionCommon,
        _u: &Point2f,
        wi: &mut Vector3f,
        pdf: &mut Float,
        vis: &mut VisibilityTester,
    ) -> Spectrum {
        *wi = (self.p_light - iref.p).normalize();
        *pdf = 1.0 as Float;
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
                p: self.p_light,
                time: iref.time,
                p_error: Vector3f::default(),
                wo: Vector3f::default(),
                n: Normal3f::default(),
                medium_interface: None,
            },
        };
        self.i * self.scale(&-*wi) / pnt3_distance_squared(&self.p_light, &iref.p)
    }
    fn power(&self) -> Spectrum {
        if let Some(mipmap) = &self.mipmap {
            mipmap.lookup_pnt_flt(
                &Point2f {
                    x: 0.5 as Float,
                    y: 0.5 as Float,
                },
                0.5 as Float,
            ) * self.i
                * 4.0 as Float
                * PI
        } else {
            Spectrum::new(1.0 as Float) * self.i * 4.0 as Float * PI
        }
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
        *ray = Ray {
            o: self.p_light,
            d: uniform_sample_sphere(u1),
            t_max: std::f32::INFINITY,
            time,
            differential: None,
            medium: None,
        };
        *n_light = Normal3f::from(ray.d);
        *pdf_pos = 1.0 as Float;
        *pdf_dir = uniform_sphere_pdf();
        self.i * self.scale(&ray.d)
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
