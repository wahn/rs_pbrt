// std
use std;
use std::f32::consts::PI;
use std::io::BufReader;
use std::sync::Arc;
// others
#[cfg(feature = "openexr")]
use half::f16;
#[cfg(feature = "openexr")]
use openexr::{FrameBufferMut, InputFile, PixelType};
// pbrt
use crate::core::geometry::{pnt2_inside_bnd2, pnt3_distance_squared};
use crate::core::geometry::{Bounds2f, Normal3f, Point2f, Point2i, Point3f, Ray, Vector3f};
use crate::core::interaction::{Interaction, InteractionCommon};
use crate::core::light::{LightFlags, VisibilityTester};
use crate::core::medium::{Medium, MediumInterface};
use crate::core::mipmap::{ImageWrap, MipMap};
use crate::core::pbrt::{Float, Spectrum};
use crate::core::reflection::cos_theta;
use crate::core::sampling::{uniform_cone_pdf, uniform_sample_cone};
use crate::core::scene::Scene;
use crate::core::transform::Transform;

// see https://stackoverflow.com/questions/36008434/how-can-i-decode-f16-to-f32-using-only-the-stable-standard-library
#[inline]
#[cfg(feature = "openexr")]
fn decode_f16(half: u16) -> f32 {
    let exp: u16 = half >> 10 & 0x1f;
    let mant: u16 = half & 0x3ff;
    let val: f32 = if exp == 0 {
        (mant as f32) * (2.0f32).powi(-24)
    } else if exp != 31 {
        (mant as f32 + 1024f32) * (2.0f32).powi(exp as i32 - 25)
    } else if mant == 0 {
        ::std::f32::INFINITY
    } else {
        ::std::f32::NAN
    };
    if half & 0x8000 != 0 {
        -val
    } else {
        val
    }
}

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
    #[cfg(feature = "openexr")]
    pub fn new(
        light_to_world: &Transform,
        medium_interface: &MediumInterface,
        i: &Spectrum,
        texname: String,
        fov: Float,
    ) -> Self {
        if texname != String::from("") {
            // https://cessen.github.io/openexr-rs/openexr/index.html
            let mut resolution: Point2i = Point2i::default();
            let mut names_and_fills: Vec<(&str, f64)> = Vec::new();
            // header
            let file_result = std::fs::File::open(texname.clone());
            if file_result.is_ok() {
                let mut file = file_result.unwrap();
                let input_file_result = InputFile::new(&mut file);
                if input_file_result.is_ok() {
                    let input_file = input_file_result.unwrap();
                    // get resolution
                    let (width, height) = input_file.header().data_dimensions();
                    resolution.x = width as i32;
                    resolution.y = height as i32;
                    // make sure the image properties are the same (see incremental_io.rs in github/openexr-rs)
                    for channel_name in ["R", "G", "B"].iter() {
                        let channel = input_file
                            .header()
                            .get_channel(channel_name)
                            .expect(&format!("Didn't find channel {}.", channel_name));
                        assert!(channel.pixel_type == PixelType::HALF);
                        names_and_fills.push((channel_name, 0.0_f64));
                    }
                    let mut pixel_data =
                        vec![
                            (f16::from_f32(0.0), f16::from_f32(0.0), f16::from_f32(0.0));
                            (resolution.x * resolution.y) as usize
                        ];
                    {
                        // read pixels
                        let mut file = std::fs::File::open(texname.clone()).unwrap();
                        let mut input_file = InputFile::new(&mut file).unwrap();
                        let mut fb = FrameBufferMut::new(resolution.x as u32, resolution.y as u32);
                        fb.insert_channels(&names_and_fills[..], &mut pixel_data);
                        input_file.read_pixels(&mut fb).unwrap();
                    }
                    // convert pixel data into Vec<Spectrum>
                    let mut texels: Vec<Spectrum> = Vec::new();
                    for idx in 0..(resolution.x * resolution.y) {
                        let (r, g, b) = pixel_data[idx as usize];
                        texels.push(Spectrum::rgb(
                            decode_f16(r.to_bits()),
                            decode_f16(g.to_bits()),
                            decode_f16(b.to_bits()),
                        ));
                    }
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
                        flags: LightFlags::DeltaPosition as u8,
                        n_samples: 1_i32,
                        medium_interface: MediumInterface::default(),
                        light_to_world: *light_to_world,
                        world_to_light: Transform::inverse(&*light_to_world),
                    };
                } else {
                    // try to open an HDR image instead (TODO: check extension upfront)
                    return ProjectionLight::new_hdr(
                        light_to_world,
                        medium_interface,
                        i,
                        texname,
                        fov,
                    );
                }
            } else {
                // try to open an HDR image instead (TODO: check extension upfront)
                return ProjectionLight::new_hdr(light_to_world, medium_interface, i, texname, fov);
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
            flags: LightFlags::DeltaPosition as u8,
            n_samples: 1_i32,
            medium_interface: MediumInterface::default(),
            light_to_world: Transform::default(),
            world_to_light: Transform::default(),
        }
    }
    pub fn new_hdr(
        light_to_world: &Transform,
        _medium_interface: &MediumInterface,
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
                    let mut texels: Vec<Spectrum> =
                        vec![Spectrum::default(); (resolution.x * resolution.y) as usize];
                    let img_result = hdr.read_image_transform(
                        |p| {
                            let rgb = p.to_hdr();
                            Spectrum::rgb(rgb[0], rgb[1], rgb[2])
                        },
                        &mut texels,
                    );
                    if img_result.is_ok() {
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
        ProjectionLight {
            projection_map: None,
            p_light: Point3f::default(),
            i: Spectrum::default(),
            light_projection: Transform::default(),
            hither: 0.0 as Float,
            yon: 0.0 as Float,
            screen_bounds: Bounds2f::default(),
            cos_total_width: 0.0 as Float,
            flags: LightFlags::DeltaPosition as u8,
            n_samples: 1_i32,
            medium_interface: MediumInterface::default(),
            light_to_world: Transform::default(),
            world_to_light: Transform::default(),
        }
    }
    pub fn projection(&self, w: &Vector3f) -> Spectrum {
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
    // Light
    pub fn sample_li(
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
        self.i * self.projection(&-*wi) / pnt3_distance_squared(&self.p_light, &iref.p)
    }
    pub fn power(&self) -> Spectrum {
        if let Some(projection_map) = &self.projection_map {
            projection_map.lookup_pnt_flt(
                &Point2f {
                    x: 0.5 as Float,
                    y: 0.5 as Float,
                },
                0.5 as Float,
            ) * self.i
                * 2.0 as Float
                * PI
                * (1.0 as Float - self.cos_total_width)
        } else {
            Spectrum::new(1.0 as Float)
                * self.i
                * 2.0 as Float
                * PI
                * (1.0 as Float - self.cos_total_width)
        }
    }
    pub fn preprocess(&self, _scene: &Scene) {}
    /// Default implementation returns no emitted radiance for a ray
    /// that escapes the scene bounds.
    pub fn le(&self, _ray: &mut Ray) -> Spectrum {
        Spectrum::new(0.0 as Float)
    }
    pub fn pdf_li(&self, _iref: &dyn Interaction, _wi: Vector3f) -> Float {
        0.0 as Float
    }
    pub fn sample_le(
        &self,
        u1: &Point2f,
        _u2: &Point2f,
        time: Float,
        ray: &mut Ray,
        n_light: &mut Normal3f,
        pdf_pos: &mut Float,
        pdf_dir: &mut Float,
    ) -> Spectrum {
        let v: Vector3f = uniform_sample_cone(u1, self.cos_total_width);
        let mut inside: Option<Arc<Medium>> = None;
        if let Some(ref mi_inside) = self.medium_interface.inside {
            inside = Some(mi_inside.clone());
        }
        *ray = Ray {
            o: self.p_light,
            d: self.light_to_world.transform_vector(&v),
            t_max: std::f32::INFINITY,
            time,
            differential: None,
            medium: inside,
        };
        *n_light = Normal3f::from(ray.d);
        *pdf_pos = 1.0 as Float;
        *pdf_dir = uniform_cone_pdf(self.cos_total_width);
        self.i * self.projection(&ray.d)
    }
    pub fn get_flags(&self) -> u8 {
        self.flags
    }
    pub fn get_n_samples(&self) -> i32 {
        self.n_samples
    }
    pub fn pdf_le(&self, ray: &Ray, _n_light: &Normal3f, pdf_pos: &mut Float, pdf_dir: &mut Float) {
        *pdf_pos = 0.0 as Float;
        if cos_theta(&self.world_to_light.transform_vector(&ray.d)) >= self.cos_total_width {
            *pdf_dir = uniform_cone_pdf(self.cos_total_width);
        } else {
            *pdf_dir = 0.0 as Float;
        }
    }
}
