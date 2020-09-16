// std
use std::cell::Cell;
use std::f32::consts::PI;
use std::io::BufReader;
use std::sync::Arc;
// others
#[cfg(feature = "openexr")]
use half::f16;
#[cfg(feature = "openexr")]
use openexr::{FrameBufferMut, InputFile, PixelType};
// pbrt
use crate::core::geometry::{pnt3_distance_squaredf, spherical_phi, spherical_theta};
use crate::core::geometry::{Normal3f, Point2f, Point2i, Point3f, Ray, Vector3f};
use crate::core::interaction::{Interaction, InteractionCommon};
use crate::core::light::{LightFlags, VisibilityTester};
use crate::core::medium::MediumInterface;
use crate::core::mipmap::{ImageWrap, MipMap};
use crate::core::pbrt::{Float, Spectrum};
use crate::core::pbrt::{INV_2_PI, INV_PI};
use crate::core::sampling::{uniform_sample_sphere, uniform_sphere_pdf};
use crate::core::scene::Scene;
use crate::core::transform::Transform;

// see https://stackoverflow.com/questions/36008434/how-can-i-decode-f16-to-f32-using-only-the-stable-standard-library
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
    #[cfg(feature = "openexr")]
    pub fn new(
        light_to_world: &Transform,
        medium_interface: &MediumInterface,
        i: &Spectrum,
        texname: String,
    ) -> Self {
        // read texel data from _texname_ and initialize _mipmap_
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
                        resolution,
                        &texels[..],
                        do_trilinear,
                        max_aniso,
                        wrap_mode,
                    ));
                    GonioPhotometricLight {
                        p_light: light_to_world.transform_point(&Point3f::default()),
                        i: *i,
                        mipmap: Some(projection_map),
                        flags: LightFlags::DeltaPosition as u8,
                        n_samples: 1_i32,
                        medium_interface: MediumInterface::default(),
                        light_to_world: Transform::default(),
                        world_to_light: Transform::default(),
                    }
                } else {
                    // try to open an HDR image instead (TODO: check extension upfront)
                    GonioPhotometricLight::new_hdr(light_to_world, medium_interface, i, texname)
                }
            } else {
                // try to open an HDR image instead (TODO: check extension upfront)
                GonioPhotometricLight::new_hdr(light_to_world, medium_interface, i, texname)
            }
        } else {
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
    }
    pub fn new_hdr(
        light_to_world: &Transform,
        _medium_interface: &MediumInterface,
        i: &Spectrum,
        texname: String,
    ) -> Self {
        if texname != "" {
            let file = std::fs::File::open(texname).unwrap();
            let reader = BufReader::new(file);
            let img_result = image::hdr::HdrDecoder::with_strictness(reader, false);
            if img_result.is_ok() {
                if let Ok(hdr) = img_result {
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
                            resolution,
                            &texels[..],
                            do_trilinear,
                            max_aniso,
                            wrap_mode,
                        ));
                        let p_light: Point3f = light_to_world.transform_point(&Point3f::default());
                        return GonioPhotometricLight {
                            p_light,
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
            mipmap.lookup_pnt_flt(st, 0.0 as Float)
        } else {
            Spectrum::new(1.0 as Float)
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
        *wi = (self.p_light - iref.p).normalize();
        *pdf = 1.0 as Float;
        light_intr.p = self.p_light;
        light_intr.time = iref.time;
        vis.p0 = Some(&iref);
        vis.p1 = Some(light_intr);
        self.i * self.scale(&-*wi) / pnt3_distance_squaredf(&self.p_light, &iref.p)
    }
    pub fn power(&self) -> Spectrum {
        if let Some(mipmap) = &self.mipmap {
            mipmap.lookup_pnt_flt(
                Point2f {
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
        self.i * self.scale(&ray.d)
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
