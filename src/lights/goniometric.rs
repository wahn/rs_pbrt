// std
use std;
use std::f32::consts::PI;
use std::io::BufReader;
use std::sync::Arc;
// pbrt
use crate::core::geometry::{Bounds2f, Normal3f, Point2f, Point2i, Point3f, Ray, Vector3f};
use crate::core::interaction::{Interaction, InteractionCommon};
use crate::core::light::{Light, LightFlags, VisibilityTester};
use crate::core::medium::{Medium, MediumInterface};
use crate::core::mipmap::{ImageWrap, MipMap};
use crate::core::pbrt::{Float, Spectrum};
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
                        return GonioPhotometricLight {
                            p_light: p_light,
                            i: *i,
                            mipmap: Some(projection_map),
                            flags: LightFlags::DeltaPosition as u8,
                            n_samples: 1_i32,
                            medium_interface: MediumInterface::default(),
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
        // WORK
        Spectrum::default()
    }
    fn power(&self) -> Spectrum {
        // WORK
        Spectrum::default()
    }
    fn preprocess(&self, _scene: &Scene) {}
    fn le(&self, _ray: &mut Ray) -> Spectrum {
        // WORK
        Spectrum::default()
    }
    fn pdf_li(&self, _iref: &dyn Interaction, _wi: Vector3f) -> Float {
        // WORK
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
        // WORK
    }
}
