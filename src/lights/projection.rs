// std
use std;
use std::f32::consts::PI;
use std::io::BufReader;
use std::sync::Arc;
// pbrt
use crate::core::geometry::{pnt2_inside_bnd2, pnt3_distance_squared};
use crate::core::geometry::{Bounds2f, Normal3f, Point2f, Point2i, Point3f, Ray, Vector3f};
use crate::core::interaction::{Interaction, InteractionCommon};
use crate::core::light::{Light, LightFlags, VisibilityTester};
use crate::core::medium::{Medium, MediumInterface};
use crate::core::mipmap::{ImageWrap, MipMap};
use crate::core::pbrt::{Float, Spectrum};
use crate::core::sampling::{uniform_sample_sphere, uniform_sphere_pdf};
use crate::core::scene::Scene;
use crate::core::transform::Transform;

// see projection.h

pub struct ProjectionLight {
    // private data (see projection.h)
    pub projection_map: Option<Arc<MipMap<Spectrum>>>,
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
    pub light_to_world: Transform,
    pub world_to_light: Transform,
}

impl ProjectionLight {
    #[cfg(not(feature = "openexr"))]
    pub fn new(
        light_to_world: &Transform,
        medium_interface: &MediumInterface,
        i: &Spectrum,
        texname: String,
        fov: Float,
    ) -> Self {
        ProjectionLight::new_hdr(light_to_world, medium_interface, i, texname, fov)
    }
    pub fn new_hdr(
        light_to_world: &Transform,
        medium_interface: &MediumInterface,
        i: &Spectrum,
        texname: String,
        fov: Float,
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
                        Spectrum::rgb(rgb[0], rgb[1], rgb[2]) * *i
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
                        let aspect: Float = resolution.x as Float / resolution.y as Float;
                        let screen_bounds: Bounds2f;
                        if aspect > 1.0 as Float {
                            screen_bounds = Bounds2f {
                                p_min: Point2f {
                                    x: -aspect,
                                    y: -1.0 as Float,
                                },
                                p_max: Point2f {
                                    x: aspect,
                                    y: 1.0 as Float,
                                },
                            };
                        } else {
                            screen_bounds = Bounds2f {
                                p_min: Point2f {
                                    x: -1.0 as Float,
                                    y: -1.0 as Float / aspect,
                                },
                                p_max: Point2f {
                                    x: 1.0 as Float,
                                    y: 1.0 as Float / aspect,
                                },
                            };
                        }
                        let hither: Float = 1e-3 as Float;
                        let yon: Float = 1e30 as Float;
                        let light_projection: Transform = Transform::perspective(fov, hither, yon);
                        let screen_to_light: Transform = Transform::inverse(&light_projection);
                        let p_corner: Point3f = Point3f {
                            x: screen_bounds.p_max.x,
                            y: screen_bounds.p_max.y,
                            z: 0.0 as Float,
                        };
                        let w_corner: Vector3f =
                            Vector3f::from(screen_to_light.transform_point(&p_corner)).normalize();
                        let cos_total_width: Float = w_corner.z;
                        return ProjectionLight {
                            projection_map: Some(projection_map),
                            p_light: p_light,
                            i: *i,
                            light_projection: light_projection,
                            hither: hither,
                            yon: yon,
                            screen_bounds: screen_bounds,
                            cos_total_width: cos_total_width,
                            flags: LightFlags::Infinite as u8,
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
        ProjectionLight {
            projection_map: None,
            p_light: Point3f::default(),
            i: Spectrum::default(),
            light_projection: Transform::default(),
            hither: 0.0 as Float,
            yon: 0.0 as Float,
            screen_bounds: Bounds2f::default(),
            cos_total_width: 0.0 as Float,
            flags: LightFlags::Infinite as u8,
            n_samples: 1_i32,
            medium_interface: MediumInterface::default(),
            light_to_world: Transform::default(),
            world_to_light: Transform::default(),
        }
    }
    pub fn projection(&self, w: &Vector3f) -> Spectrum {
        // Vector3f wl = world_to_light(w);
        let wl: Vector3f = self.world_to_light.transform_vector(w);
        // discard directions behind projection light
        if wl.z < self.hither {
            return Spectrum::default();
        }
        // project point onto projection plane and compute light
        let p: Point3f = self.light_projection.transform_point(&Point3f {
            x: wl.x,
            y: wl.y,
            z: wl.z,
        });
        if !pnt2_inside_bnd2(&Point2f { x: p.x, y: p.y }, &self.screen_bounds) {
            return Spectrum::default();
        }
        if let Some(projection_map) = &self.projection_map {
            let st: Point2f = Point2f::from(self.screen_bounds.offset(&Point2f { x: p.x, y: p.y }));
            projection_map.lookup_pnt_flt(&st, 0.0 as Float)
        } else {
            Spectrum::new(1.0 as Float)
        }
    }
}

impl Light for ProjectionLight {
    fn sample_li(
        &self,
        iref: &InteractionCommon,
        _u: &Point2f,
        wi: &mut Vector3f,
        pdf: &mut Float,
        vis: &mut VisibilityTester,
    ) -> Spectrum {
        // TODO: ProfilePhase _(Prof::LightSample);
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
        self.i * self.projection(&-*wi) / pnt3_distance_squared(&self.p_light, &iref.p)
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
