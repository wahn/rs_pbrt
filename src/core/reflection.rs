//! When light is incident on the surface, the surface scatters the
//! light, reflecting some of it back into the environment. There are
//! two main effects that need to be described to model this
//! reflection: the spectral distribution of the reflected light and
//! its directional distribution.

// std
use std;
use std::f32::consts::PI;
use std::fs::File;
use std::io::Read;
use std::path::Path;
use std::sync::Arc;
// others
use byteorder::{LittleEndian, ReadBytesExt};
use num::Zero;
// pbrt
use crate::core::bssrdf::SeparableBssrdfAdapter;
use crate::core::geometry::{
    nrm_cross_vec3, nrm_dot_vec3, nrm_faceforward_vec3, vec3_abs_dot_vec3, vec3_dot_nrm,
    vec3_dot_vec3,
};
use crate::core::geometry::{Normal3f, Point2f, Vector3f};
use crate::core::interaction::SurfaceInteraction;
use crate::core::interpolation::{
    catmull_rom_weights, fourier, sample_catmull_rom_2d, sample_fourier,
};
use crate::core::material::TransportMode;
use crate::core::microfacet::{MicrofacetDistribution, TrowbridgeReitzDistribution};
use crate::core::pbrt::INV_PI;
use crate::core::pbrt::{clamp_t, lerp, radians};
use crate::core::pbrt::{Float, Spectrum};
use crate::core::rng::FLOAT_ONE_MINUS_EPSILON;
use crate::core::sampling::cosine_sample_hemisphere;
use crate::materials::disney::{
    DisneyClearCoat, DisneyDiffuse, DisneyFakeSS, DisneyRetro, DisneySheen,
};
use crate::materials::hair::HairBSDF;

/// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
///
/// The Schlick Fresnel approximation is:
///
/// R = R(0) + (1 - R(0)) (1 - cos theta)^5,
///
/// where R(0) is the reflectance at normal indicence.
#[inline]
fn schlick_weight(cos_theta: Float) -> Float {
    let m = clamp_t(1.0 - cos_theta, 0.0, 1.0);
    (m * m) * (m * m) * m
}

#[inline]
pub fn fr_schlick(r0: Float, cos_theta: Float) -> Float {
    lerp(schlick_weight(cos_theta), r0, 1.0)
}

#[inline]
fn fr_schlick_spectrum(r0: Spectrum, cos_theta: Float) -> Spectrum {
    lerp(schlick_weight(cos_theta), r0, Spectrum::from(1.0))
}

// see reflection.h

#[derive(Default)]
pub struct FourierBSDFTable {
    pub eta: Float,
    pub m_max: i32,
    pub n_channels: i32,
    pub n_mu: i32,
    pub mu: Vec<Float>,
    pub m: Vec<i32>,
    pub a_offset: Vec<i32>,
    pub a: Vec<Float>,
    pub a0: Vec<Float>,
    pub cdf: Vec<Float>,
    pub recip: Vec<Float>,
}

impl FourierBSDFTable {
    pub fn read(&mut self, filename: &String) -> bool {
        let path = Path::new(&filename);
        let result = File::open(path);
        if !result.is_ok() {
            println!("ERROR: Unable to open tabulated BSDF file {:?}", filename);
            return false;
        }
        // header
        let mut file = result.unwrap();
        let mut buffer = [0; 8];
        let io_result = file.read_exact(&mut buffer);
        if io_result.is_ok() {
            let header_exp: [u8; 8] = [b'S', b'C', b'A', b'T', b'F', b'U', b'N', 0x01_u8];
            if buffer == header_exp {
                let mut buffer: [i32; 9] = [0; 9]; // 9 32-bit (signed) integers (the last 3 are unused)
                let io_result = file.read_i32_into::<LittleEndian>(&mut buffer);
                if io_result.is_ok() {
                    let flags: i32 = buffer[0];
                    self.n_mu = buffer[1];
                    let n_coeffs: i32 = buffer[2];
                    self.m_max = buffer[3];
                    self.n_channels = buffer[4];
                    let n_bases: i32 = buffer[5];
                    let mut buffer: [f32; 1] = [0_f32; 1]; // 1 32-bit float
                    let io_result = file.read_f32_into::<LittleEndian>(&mut buffer);
                    if io_result.is_ok() {
                        self.eta = buffer[0];
                        let mut buffer: [i32; 4] = [0; 4]; // 4 32-bit (signed) integers are unused
                        let io_result = file.read_i32_into::<LittleEndian>(&mut buffer);
                        if io_result.is_ok() {
                            // only a subset of BSDF files are
                            // supported for simplicity, in
                            // particular: monochromatic and RGB files
                            // with uniform (i.e. non-textured)
                            // material properties
                            if flags != 1_i32
                                || (self.n_channels != 1_i32 && self.n_channels != 3_i32)
                                || n_bases != 1_i32
                            {
                                panic!(
                                    "ERROR: Tabulated BSDF file {:?} has an incompatible file format or version."
                                );
                            }
                            // self.mu
                            self.mu.reserve_exact(self.n_mu as usize);
                            for _ in 0..self.n_mu as usize {
                                let f: f32 = file.read_f32::<LittleEndian>().unwrap();
                                self.mu.push(f as Float);
                            }
                            // self.cdf
                            self.cdf
                                .reserve_exact(self.n_mu as usize * self.n_mu as usize);
                            for _ in 0..(self.n_mu as usize * self.n_mu as usize) {
                                let f: f32 = file.read_f32::<LittleEndian>().unwrap();
                                self.cdf.push(f as Float);
                            }
                            // self.a0
                            self.a0
                                .reserve_exact(self.n_mu as usize * self.n_mu as usize);
                            // offset_and_length
                            let mut offset_and_length: Vec<i32> = Vec::with_capacity(
                                self.n_mu as usize * self.n_mu as usize * 2_usize,
                            );
                            for _ in 0..(self.n_mu as usize * self.n_mu as usize * 2_usize) {
                                let i: i32 = file.read_i32::<LittleEndian>().unwrap();
                                offset_and_length.push(i);
                            }
                            // self.a_offset
                            self.a_offset
                                .reserve_exact(self.n_mu as usize * self.n_mu as usize);
                            // self.m
                            self.m
                                .reserve_exact(self.n_mu as usize * self.n_mu as usize);
                            // self.a
                            self.a.reserve_exact(n_coeffs as usize);
                            for _ in 0..n_coeffs as usize {
                                let f: f32 = file.read_f32::<LittleEndian>().unwrap();
                                self.a.push(f as Float);
                            }
                            // fill self.a_offset, self.m, and self.a0 vectors
                            for i in 0..(self.n_mu as usize * self.n_mu as usize) {
                                let offset: i32 = offset_and_length[(2 * i) as usize];
                                let length: i32 = offset_and_length[(2 * i + 1) as usize];
                                self.a_offset.push(offset);
                                self.m.push(length);
                                if length > 0 {
                                    self.a0.push(self.a[offset as usize]);
                                } else {
                                    self.a0.push(0.0 as Float);
                                }
                            }
                            // self.recip
                            self.recip.reserve_exact(self.m_max as usize);
                            for i in 0..self.m_max as usize {
                                self.recip.push(1.0 as Float / i as Float);
                            }
                        } else {
                            panic!(
                                "ERROR: Tabulated BSDF file {:?} has an incompatible file format or version."
                            );
                        }
                    } else {
                        panic!(
                            "ERROR: Tabulated BSDF file {:?} has an incompatible file format or version."
                        );
                    }
                } else {
                    panic!(
                        "ERROR: Tabulated BSDF file {:?} has an incompatible file format or version."
                    );
                }
            } else {
                panic!(
                    "ERROR: Tabulated BSDF file {:?} has an incompatible file format or version."
                );
            }
        }
        true
    }
    pub fn get_ak(&self, offset_i: i32, offset_o: i32, mptr: &mut i32) -> i32 {
        let idx: i32 = offset_o * self.n_mu + offset_i;
        assert!(
            idx >= 0,
            "get_ak({:?}, {:?}, ...) with idx = {:?}",
            offset_i,
            offset_o,
            idx
        );
        *mptr = self.m[idx as usize];
        self.a_offset[idx as usize]
    }
    pub fn get_weights_and_offset(
        &self,
        cos_theta: Float,
        offset: &mut i32,
        weights: &mut [Float; 4],
    ) -> bool {
        catmull_rom_weights(&self.mu, cos_theta, offset, weights)
    }
}

pub struct Bsdf {
    pub eta: Float,
    /// shading normal
    pub ns: Normal3f,
    /// geometric normal
    pub ng: Normal3f,
    pub ss: Vector3f,
    pub ts: Vector3f,
    pub bxdfs: Vec<Bxdf>,
}

impl Bsdf {
    pub fn new(si: &SurfaceInteraction, eta: Float, bxdfs: Vec<Bxdf>) -> Self {
        let ss = si.shading.dpdu.normalize();
        Bsdf {
            eta,
            ns: si.shading.n,
            ng: si.n,
            ss,
            ts: nrm_cross_vec3(&si.shading.n, &ss),
            bxdfs,
        }
    }
    pub fn num_components(&self, flags: u8) -> u8 {
        let mut num: u8 = 0;
        let n_bxdfs: usize = self.bxdfs.len();
        for i in 0..n_bxdfs {
            if self.bxdfs[i].matches_flags(flags) {
                num += 1;
            }
        }
        num
    }
    pub fn world_to_local(&self, v: &Vector3f) -> Vector3f {
        Vector3f {
            x: vec3_dot_vec3(v, &self.ss),
            y: vec3_dot_vec3(v, &self.ts),
            z: vec3_dot_vec3(v, &Vector3f::from(self.ns)),
        }
    }
    pub fn local_to_world(&self, v: &Vector3f) -> Vector3f {
        Vector3f {
            x: self.ss.x * v.x + self.ts.x * v.y + self.ns.x * v.z,
            y: self.ss.y * v.x + self.ts.y * v.y + self.ns.y * v.z,
            z: self.ss.z * v.x + self.ts.z * v.y + self.ns.z * v.z,
        }
    }
    pub fn f(&self, wo_w: &Vector3f, wi_w: &Vector3f, flags: u8) -> Spectrum {
        // TODO: ProfilePhase pp(Prof::BSDFEvaluation);
        let wi: Vector3f = self.world_to_local(wi_w);
        let wo: Vector3f = self.world_to_local(wo_w);
        if wo.z == 0.0 as Float {
            return Spectrum::new(0.0 as Float);
        }
        let reflect: bool = (vec3_dot_vec3(wi_w, &Vector3f::from(self.ng))
            * vec3_dot_vec3(wo_w, &Vector3f::from(self.ng)))
            > 0.0 as Float;
        let mut f: Spectrum = Spectrum::new(0.0 as Float);
        let n_bxdfs: usize = self.bxdfs.len();
        for i in 0..n_bxdfs {
            if self.bxdfs[i].matches_flags(flags)
                && ((reflect && (self.bxdfs[i].get_type() & BxdfType::BsdfReflection as u8 > 0_u8))
                    || (!reflect
                        && (self.bxdfs[i].get_type() & BxdfType::BsdfTransmission as u8 > 0_u8)))
            {
                f += self.bxdfs[i].f(&wo, &wi);
            }
        }
        f
    }
    /// Calls the individual Bxdf::sample_f() methods to generate samples.
    pub fn sample_f(
        &self,
        wo_world: &Vector3f,
        wi_world: &mut Vector3f,
        u: &Point2f,
        pdf: &mut Float,
        bsdf_flags: u8,
        sampled_type: &mut u8,
    ) -> Spectrum {
        // TODO: ProfilePhase pp(Prof::BSDFSampling);
        // choose which _BxDF_ to sample
        let matching_comps: u8 = self.num_components(bsdf_flags);
        if matching_comps == 0 {
            *pdf = 0.0 as Float;
            *sampled_type = 0_u8;
            return Spectrum::default();
        }
        let comp: u8 = std::cmp::min(
            (u[0] * matching_comps as Float).floor() as u8,
            matching_comps - 1_u8,
        );
        // get _BxDF_ pointer for chosen component
        let mut bxdf: Option<&Bxdf> = None;
        let mut count: i8 = comp as i8;
        let n_bxdfs: usize = self.bxdfs.len();
        let mut bxdf_index: usize = 0_usize;
        for i in 0..n_bxdfs {
            let matches: bool = self.bxdfs[i].matches_flags(bsdf_flags);
            if matches && count == 0 {
                count -= 1_i8;
                bxdf = self.bxdfs.get(i);
                bxdf_index = i;
                break;
            } else {
                // fix count
                if matches {
                    // C++ version does this in a single line:
                    // if (bxdfs[i]->MatchesFlags(type) && count-- == 0)
                    count -= 1_i8;
                }
            }
        }
        if bxdf.is_some() {
            let bxdf = bxdf.unwrap();
            // TODO: println!("BSDF::Sample_f chose comp = {:?} /
            // matching = {:?}, bxdf: {:?}", comp, matching_comps,
            // bxdf);

            // remap _BxDF_ sample _u_ to $[0,1)^2$
            let u_remapped: Point2f = Point2f {
                x: (u[0] * matching_comps as Float - comp as Float).min(FLOAT_ONE_MINUS_EPSILON),
                y: u[1],
            };
            // sample chosen _BxDF_
            let mut wi: Vector3f = Vector3f::default();
            let wo: Vector3f = self.world_to_local(wo_world);
            if wo.z == 0.0 as Float {
                return Spectrum::default();
            }
            *pdf = 0.0 as Float;
            if *sampled_type != 0_u8 {
                *sampled_type = bxdf.get_type();
            }
            let mut f: Spectrum = bxdf.sample_f(&wo, &mut wi, &u_remapped, pdf, sampled_type);
            // let mut ratio: Spectrum = Spectrum::default();
            // if *pdf > 0.0 as Float {
            //     ratio = f / *pdf;
            // }
            // println!("For wo = {:?}, sampled f = {:?}, pdf = {:?}, ratio = {:?}, wi = {:?}",
            //          wo,
            //          f,
            //          *pdf,
            //          ratio,
            //          wi);
            if *pdf == 0.0 as Float {
                if *sampled_type != 0_u8 {
                    *sampled_type = 0_u8;
                }
                return Spectrum::default();
            }
            *wi_world = self.local_to_world(&wi);
            // compute overall PDF with all matching _BxDF_s
            if (bxdf.get_type() & BxdfType::BsdfSpecular as u8 == 0_u8) && matching_comps > 1_u8 {
                for i in 0..n_bxdfs {
                    // instead of self.bxdfs[i] != bxdf we compare stored index
                    if bxdf_index != i && self.bxdfs[i].matches_flags(bsdf_flags) {
                        *pdf += self.bxdfs[i].pdf(&wo, &wi);
                    }
                }
            }
            if matching_comps > 1_u8 {
                *pdf /= matching_comps as Float;
            }
            // compute value of BSDF for sampled direction
            if bxdf.get_type() & BxdfType::BsdfSpecular as u8 == 0_u8 {
                let reflect: bool = vec3_dot_nrm(&*wi_world, &self.ng)
                    * vec3_dot_nrm(wo_world, &self.ng)
                    > 0.0 as Float;
                f = Spectrum::default();
                for i in 0..n_bxdfs {
                    if self.bxdfs[i].matches_flags(bsdf_flags)
                        && ((reflect
                            && ((self.bxdfs[i].get_type() & BxdfType::BsdfReflection as u8)
                                != 0_u8))
                            || (!reflect
                                && ((self.bxdfs[i].get_type() & BxdfType::BsdfTransmission as u8)
                                    != 0_u8)))
                    {
                        f += self.bxdfs[i].f(&wo, &wi);
                    }
                }
            }
            // let mut ratio: Spectrum = Spectrum::default();
            // if *pdf > 0.0 as Float {
            //     ratio = f / *pdf;
            // }
            // println!("Overall f = {:?}, pdf = {:?}, ratio = {:?}", f, *pdf, ratio);
            return f;
        } else {
            panic!("CHECK_NOTNULL(bxdf)");
        }
    }
    pub fn pdf(&self, wo_world: &Vector3f, wi_world: &Vector3f, bsdf_flags: u8) -> Float {
        // TODO: ProfilePhase pp(Prof::BSDFPdf);
        let n_bxdfs: usize = self.bxdfs.len();
        if n_bxdfs == 0 {
            return 0.0 as Float;
        }
        let wo: Vector3f = self.world_to_local(wo_world);
        let wi: Vector3f = self.world_to_local(wi_world);
        if wo.z == 0.0 as Float {
            return 0.0 as Float;
        }
        let mut pdf: Float = 0.0 as Float;
        let mut matching_comps: u8 = 0;
        for i in 0..n_bxdfs {
            if self.bxdfs[i].matches_flags(bsdf_flags) {
                matching_comps += 1;
                pdf += self.bxdfs[i].pdf(&wo, &wi);
            }
        }
        let mut v: Float = 0.0 as Float;
        if matching_comps > 0 {
            v = pdf / matching_comps as Float;
        }
        v
    }
}

#[repr(u8)]
pub enum BxdfType {
    BsdfReflection = 1,
    BsdfTransmission = 2,
    BsdfDiffuse = 4,
    BsdfGlossy = 8,
    BsdfSpecular = 16,
    BsdfAll = 31,
}

pub enum Bxdf {
    // Scaled(ScaledBxDF),
    SpecRefl(SpecularReflection),
    SpecTrans(SpecularTransmission),
    FresnelSpec(FresnelSpecular),
    LambertianRefl(LambertianReflection),
    LambertianTrans(LambertianTransmission),
    OrenNayarRefl(OrenNayar),
    MicrofacetRefl(MicrofacetReflection),
    MicrofacetTrans(MicrofacetTransmission),
    FresnelBlnd(FresnelBlend),
    Fourier(FourierBSDF),
    // bssrdf.rs
    Bssrdf(SeparableBssrdfAdapter),
    // disney.rs
    DisDiff(DisneyDiffuse),
    DisSS(DisneyFakeSS),
    DisRetro(DisneyRetro),
    DisSheen(DisneySheen),
    DisClearCoat(DisneyClearCoat),
    // hair.rs
    Hair(HairBSDF),
}

impl Bxdf {
    pub fn matches_flags(&self, t: u8) -> bool {
        match self {
            // Bxdf::Scaled(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
            Bxdf::SpecRefl(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
            Bxdf::SpecTrans(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
            Bxdf::FresnelSpec(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
            Bxdf::LambertianRefl(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
            Bxdf::LambertianTrans(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
            Bxdf::OrenNayarRefl(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
            Bxdf::MicrofacetRefl(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
            Bxdf::MicrofacetTrans(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
            Bxdf::FresnelBlnd(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
            Bxdf::Fourier(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
            Bxdf::Bssrdf(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
            Bxdf::DisDiff(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
            Bxdf::DisSS(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
            Bxdf::DisRetro(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
            Bxdf::DisSheen(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
            Bxdf::DisClearCoat(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
            Bxdf::Hair(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
        }
    }
    pub fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        match self {
            // Bxdf::Scaled(bxdf) => bxdf.f(wo, wi),
            Bxdf::SpecRefl(bxdf) => bxdf.f(wo, wi),
            Bxdf::SpecTrans(bxdf) => bxdf.f(wo, wi),
            Bxdf::FresnelSpec(bxdf) => bxdf.f(wo, wi),
            Bxdf::LambertianRefl(bxdf) => bxdf.f(wo, wi),
            Bxdf::LambertianTrans(bxdf) => bxdf.f(wo, wi),
            Bxdf::OrenNayarRefl(bxdf) => bxdf.f(wo, wi),
            Bxdf::MicrofacetRefl(bxdf) => bxdf.f(wo, wi),
            Bxdf::MicrofacetTrans(bxdf) => bxdf.f(wo, wi),
            Bxdf::FresnelBlnd(bxdf) => bxdf.f(wo, wi),
            Bxdf::Fourier(bxdf) => bxdf.f(wo, wi),
            Bxdf::Bssrdf(bxdf) => bxdf.f(wo, wi),
            Bxdf::DisDiff(bxdf) => bxdf.f(wo, wi),
            Bxdf::DisSS(bxdf) => bxdf.f(wo, wi),
            Bxdf::DisRetro(bxdf) => bxdf.f(wo, wi),
            Bxdf::DisSheen(bxdf) => bxdf.f(wo, wi),
            Bxdf::DisClearCoat(bxdf) => bxdf.f(wo, wi),
            Bxdf::Hair(bxdf) => bxdf.f(wo, wi),
        }
    }
    /// Sample the BxDF for the given outgoing direction, using the given pair of uniform samples.
    ///
    /// The default implementation uses importance sampling by using a cosine-weighted
    /// distribution.
    pub fn sample_f(
        &self,
        wo: &Vector3f,
        wi: &mut Vector3f,
        u: &Point2f,
        pdf: &mut Float,
        sampled_type: &mut u8,
    ) -> Spectrum {
        match self {
            // Bxdf::Scaled(bxdf) => bxdf.sample_f(wo, wi, u, pdf, sampled_type),
            Bxdf::SpecRefl(bxdf) => bxdf.sample_f(wo, wi, u, pdf, sampled_type),
            Bxdf::SpecTrans(bxdf) => bxdf.sample_f(wo, wi, u, pdf, sampled_type),
            Bxdf::FresnelSpec(bxdf) => bxdf.sample_f(wo, wi, u, pdf, sampled_type),
            Bxdf::LambertianRefl(bxdf) => bxdf.sample_f(wo, wi, u, pdf, sampled_type),
            Bxdf::LambertianTrans(bxdf) => bxdf.sample_f(wo, wi, u, pdf, sampled_type),
            Bxdf::OrenNayarRefl(bxdf) => bxdf.sample_f(wo, wi, u, pdf, sampled_type),
            Bxdf::MicrofacetRefl(bxdf) => bxdf.sample_f(wo, wi, u, pdf, sampled_type),
            Bxdf::MicrofacetTrans(bxdf) => bxdf.sample_f(wo, wi, u, pdf, sampled_type),
            Bxdf::FresnelBlnd(bxdf) => bxdf.sample_f(wo, wi, u, pdf, sampled_type),
            Bxdf::Fourier(bxdf) => bxdf.sample_f(wo, wi, u, pdf, sampled_type),
            Bxdf::Bssrdf(_bxdf) => self.sample_f(wo, wi, u, pdf, sampled_type),
            Bxdf::DisDiff(_bxdf) => self.sample_f(wo, wi, u, pdf, sampled_type),
            Bxdf::DisSS(_bxdf) => self.sample_f(wo, wi, u, pdf, sampled_type),
            Bxdf::DisRetro(_bxdf) => self.sample_f(wo, wi, u, pdf, sampled_type),
            Bxdf::DisSheen(_bxdf) => self.sample_f(wo, wi, u, pdf, sampled_type),
            Bxdf::DisClearCoat(bxdf) => bxdf.sample_f(wo, wi, u, pdf, sampled_type),
            Bxdf::Hair(bxdf) => bxdf.sample_f(wo, wi, u, pdf, sampled_type),
        }
    }
    /// Evaluate the PDF for the given outgoing and incoming directions.
    ///
    /// Note: this method needs to be consistent with ```Bxdf::sample_f()```.
    pub fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        match self {
            // Bxdf::Scaled(bxdf) => bxdf.pdf(wo, wi),
            Bxdf::SpecRefl(bxdf) => bxdf.pdf(wo, wi),
            Bxdf::SpecTrans(bxdf) => bxdf.pdf(wo, wi),
            Bxdf::FresnelSpec(bxdf) => bxdf.pdf(wo, wi),
            Bxdf::LambertianRefl(bxdf) => bxdf.pdf(wo, wi),
            Bxdf::LambertianTrans(bxdf) => bxdf.pdf(wo, wi),
            Bxdf::OrenNayarRefl(bxdf) => bxdf.pdf(wo, wi),
            Bxdf::MicrofacetRefl(bxdf) => bxdf.pdf(wo, wi),
            Bxdf::MicrofacetTrans(bxdf) => bxdf.pdf(wo, wi),
            Bxdf::FresnelBlnd(bxdf) => bxdf.pdf(wo, wi),
            Bxdf::Fourier(bxdf) => bxdf.pdf(wo, wi),
            Bxdf::Bssrdf(_bxdf) => self.pdf(wo, wi),
            Bxdf::DisDiff(_bxdf) => self.pdf(wo, wi),
            Bxdf::DisSS(_bxdf) => self.pdf(wo, wi),
            Bxdf::DisRetro(_bxdf) => self.pdf(wo, wi),
            Bxdf::DisSheen(_bxdf) => self.pdf(wo, wi),
            Bxdf::DisClearCoat(bxdf) => bxdf.pdf(wo, wi),
            Bxdf::Hair(bxdf) => bxdf.pdf(wo, wi),
        }
    }
    pub fn get_type(&self) -> u8 {
        match self {
            // Bxdf::Scaled(bxdf) => bxdf.get_type(),
            Bxdf::SpecRefl(bxdf) => bxdf.get_type(),
            Bxdf::SpecTrans(bxdf) => bxdf.get_type(),
            Bxdf::FresnelSpec(bxdf) => bxdf.get_type(),
            Bxdf::LambertianRefl(bxdf) => bxdf.get_type(),
            Bxdf::LambertianTrans(bxdf) => bxdf.get_type(),
            Bxdf::OrenNayarRefl(bxdf) => bxdf.get_type(),
            Bxdf::MicrofacetRefl(bxdf) => bxdf.get_type(),
            Bxdf::MicrofacetTrans(bxdf) => bxdf.get_type(),
            Bxdf::FresnelBlnd(bxdf) => bxdf.get_type(),
            Bxdf::Fourier(bxdf) => bxdf.get_type(),
            Bxdf::Bssrdf(bxdf) => bxdf.get_type(),
            Bxdf::DisDiff(bxdf) => bxdf.get_type(),
            Bxdf::DisSS(bxdf) => bxdf.get_type(),
            Bxdf::DisRetro(bxdf) => bxdf.get_type(),
            Bxdf::DisSheen(bxdf) => bxdf.get_type(),
            Bxdf::DisClearCoat(bxdf) => bxdf.get_type(),
            Bxdf::Hair(bxdf) => bxdf.get_type(),
        }
    }
}

// pub struct ScaledBxDF {
//     pub bxdf: Bxdf,
//     pub scale: Spectrum,
// }

// impl ScaledBxDF {
//     pub fn new(bxdf: Bxdf, scale: Spectrum) -> Self {
//         ScaledBxDF { bxdf, scale }
//     }
//     pub fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
//         self.scale * self.bxdf.f(wo, wi)
//     }
//     pub fn sample_f(
//         &self,
//         wo: &Vector3f,
//         wi: &mut Vector3f,
//         sample: &Point2f,
//         pdf: &mut Float,
//         sampled_type: &mut u8,
//     ) -> Spectrum {
//         let f: Spectrum = self.bxdf.sample_f(wo, wi, sample, pdf, sampled_type);
//         self.scale * f
//     }
//     pub fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
//         self.bxdf.pdf(wo, wi)
//     }
//     pub fn get_type(&self) -> u8 {
//         self.bxdf.get_type()
//     }
// }

pub enum Fresnel {
    NoOp(FresnelNoOp),
    Conductor(FresnelConductor),
    Dielectric(FresnelDielectric),
    Disney(DisneyFresnel),
}

impl Fresnel {
    pub fn evaluate(&self, cos_theta_i: Float) -> Spectrum {
        match self {
            Fresnel::NoOp(fresnel) => fresnel.evaluate(cos_theta_i),
            Fresnel::Conductor(fresnel) => fresnel.evaluate(cos_theta_i),
            Fresnel::Dielectric(fresnel) => fresnel.evaluate(cos_theta_i),
            Fresnel::Disney(fresnel) => fresnel.evaluate(cos_theta_i),
        }
    }
}

/// Specialized Fresnel function used for the specular component, based on
/// a mixture between dielectric and the Schlick Fresnel approximation.
#[derive(Debug, Clone, Copy)]
pub struct DisneyFresnel {
    r0: Spectrum,
    metallic: Float,
    eta: Float,
}

impl DisneyFresnel {
    pub fn new(r0: Spectrum, metallic: Float, eta: Float) -> DisneyFresnel {
        DisneyFresnel { r0, metallic, eta }
    }
    pub fn evaluate(&self, cos_i: Float) -> Spectrum {
        lerp(
            self.metallic,
            Spectrum::from(fr_dielectric(cos_i, 1.0, self.eta)),
            fr_schlick_spectrum(self.r0, cos_i),
        )
    }
}

#[derive(Debug, Default, Copy, Clone)]
pub struct FresnelConductor {
    pub eta_i: Spectrum,
    pub eta_t: Spectrum,
    pub k: Spectrum,
}

impl FresnelConductor {
    pub fn evaluate(&self, cos_theta_i: Float) -> Spectrum {
        fr_conductor(cos_theta_i, self.eta_i, self.eta_t, self.k)
    }
}

#[derive(Debug, Default, Copy, Clone)]
pub struct FresnelDielectric {
    pub eta_i: Float,
    pub eta_t: Float,
}

impl FresnelDielectric {
    pub fn evaluate(&self, cos_theta_i: Float) -> Spectrum {
        Spectrum::new(fr_dielectric(cos_theta_i, self.eta_i, self.eta_t))
    }
}

#[derive(Debug, Default, Copy, Clone)]
pub struct FresnelNoOp {}

impl FresnelNoOp {
    pub fn evaluate(&self, _cos_theta_i: Float) -> Spectrum {
        Spectrum::new(1.0 as Float)
    }
}

pub struct SpecularReflection {
    pub r: Spectrum,
    pub fresnel: Fresnel,
}

impl SpecularReflection {
    pub fn new(r: Spectrum, fresnel: Fresnel) -> Self {
        SpecularReflection { r, fresnel }
    }
    pub fn f(&self, _wo: &Vector3f, _wi: &Vector3f) -> Spectrum {
        Spectrum::new(0.0 as Float)
    }
    pub fn sample_f(
        &self,
        wo: &Vector3f,
        wi: &mut Vector3f,
        _sample: &Point2f,
        pdf: &mut Float,
        _sampled_type: &mut u8,
    ) -> Spectrum {
        // compute perfect specular reflection direction
        *wi = Vector3f {
            x: -wo.x,
            y: -wo.y,
            z: wo.z,
        };
        *pdf = 1.0 as Float;
        let cos_theta_i: Float = cos_theta(&*wi);
        self.fresnel.evaluate(cos_theta_i) * self.r / abs_cos_theta(&*wi)
    }
    pub fn pdf(&self, _wo: &Vector3f, _wi: &Vector3f) -> Float {
        0.0 as Float
    }
    pub fn get_type(&self) -> u8 {
        BxdfType::BsdfReflection as u8 | BxdfType::BsdfSpecular as u8
    }
}

pub struct SpecularTransmission {
    pub t: Spectrum,
    pub eta_a: Float,
    pub eta_b: Float,
    pub fresnel: FresnelDielectric,
    pub mode: TransportMode,
}

impl SpecularTransmission {
    pub fn new(t: Spectrum, eta_a: Float, eta_b: Float, mode: TransportMode) -> Self {
        SpecularTransmission {
            t,
            eta_a,
            eta_b,
            fresnel: FresnelDielectric {
                eta_i: eta_a,
                eta_t: eta_b,
            },
            mode,
        }
    }
    pub fn f(&self, _wo: &Vector3f, _wi: &Vector3f) -> Spectrum {
        Spectrum::new(0.0 as Float)
    }
    pub fn sample_f(
        &self,
        wo: &Vector3f,
        wi: &mut Vector3f,
        _sample: &Point2f,
        pdf: &mut Float,
        _sampled_type: &mut u8,
    ) -> Spectrum {
        // figure out which $\eta$ is incident and which is transmitted
        let entering: bool = cos_theta(wo) > 0.0;
        let mut eta_i: Float = self.eta_b;
        if entering {
            eta_i = self.eta_a;
        }
        let mut eta_t: Float = self.eta_a;
        if entering {
            eta_t = self.eta_b;
        }
        // compute ray direction for specular transmission
        if !refract(
            wo,
            &nrm_faceforward_vec3(
                &Normal3f {
                    x: 0.0,
                    y: 0.0,
                    z: 1.0,
                },
                wo,
            ),
            eta_i / eta_t,
            wi,
        ) {
            return Spectrum::default();
        }
        *pdf = 1.0;
        let mut ft: Spectrum =
            self.t * (Spectrum::new(1.0 as Float) - self.fresnel.evaluate(cos_theta(&*wi)));
        // account for non-symmetry with transmission to different medium
        if self.mode == TransportMode::Radiance {
            ft *= Spectrum::new((eta_i * eta_i) / (eta_t * eta_t));
        }
        ft / abs_cos_theta(&*wi)
    }
    pub fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        if vec3_same_hemisphere_vec3(wo, wi) {
            abs_cos_theta(wi) * INV_PI
        } else {
            0.0 as Float
        }
    }
    pub fn get_type(&self) -> u8 {
        BxdfType::BsdfTransmission as u8 | BxdfType::BsdfSpecular as u8
    }
}

pub struct FresnelSpecular {
    pub r: Spectrum,
    pub t: Spectrum,
    pub eta_a: Float,
    pub eta_b: Float,
    pub mode: TransportMode,
}

impl FresnelSpecular {
    pub fn new(r: Spectrum, t: Spectrum, eta_a: Float, eta_b: Float, mode: TransportMode) -> Self {
        FresnelSpecular {
            r,
            t,
            eta_a,
            eta_b,
            mode,
        }
    }
    pub fn f(&self, _wo: &Vector3f, _wi: &Vector3f) -> Spectrum {
        Spectrum::new(0.0 as Float)
    }
    pub fn sample_f(
        &self,
        wo: &Vector3f,
        wi: &mut Vector3f,
        sample: &Point2f,
        pdf: &mut Float,
        sampled_type: &mut u8,
    ) -> Spectrum {
        let ct: Float = cos_theta(wo);
        let f: Float = fr_dielectric(ct, self.eta_a, self.eta_b);
        if sample[0] < f {
            // compute specular reflection for _FresnelSpecular_

            // compute perfect specular reflection direction
            *wi = Vector3f {
                x: -wo.x,
                y: -wo.y,
                z: wo.z,
            };
            if *sampled_type != 0_u8 {
                *sampled_type = BxdfType::BsdfReflection as u8 | BxdfType::BsdfSpecular as u8
            }
            *pdf = f;
            return self.r * f / abs_cos_theta(&*wi);
        } else {
            // compute specular transmission for _FresnelSpecular_

            // figure out which $\eta$ is incident and which is transmitted
            let entering: bool = cos_theta(wo) > 0.0 as Float;
            let eta_i: Float;
            if entering {
                eta_i = self.eta_a;
            } else {
                eta_i = self.eta_b;
            }
            let eta_t: Float;
            if entering {
                eta_t = self.eta_b;
            } else {
                eta_t = self.eta_a;
            }
            // compute ray direction for specular transmission
            if !refract(
                wo,
                &nrm_faceforward_vec3(
                    &Normal3f {
                        x: 0.0,
                        y: 0.0,
                        z: 1.0,
                    },
                    wo,
                ),
                eta_i / eta_t,
                wi,
            ) {
                return Spectrum::default();
            }
            let mut ft: Spectrum = self.t * (1.0 as Float - f);
            // account for non-symmetry with transmission to different medium
            if self.mode == TransportMode::Radiance {
                ft *= Spectrum::new((eta_i * eta_i) / (eta_t * eta_t));
            }
            if *sampled_type != 0_u8 {
                *sampled_type = BxdfType::BsdfTransmission as u8 | BxdfType::BsdfSpecular as u8
            }
            *pdf = 1.0 as Float - f;
            return ft / abs_cos_theta(&*wi);
        }
    }
    pub fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        if vec3_same_hemisphere_vec3(wo, wi) {
            abs_cos_theta(wi) * INV_PI
        } else {
            0.0 as Float
        }
    }
    pub fn get_type(&self) -> u8 {
        BxdfType::BsdfReflection as u8
            | BxdfType::BsdfTransmission as u8
            | BxdfType::BsdfSpecular as u8
    }
}

#[derive(Debug, Default, Copy, Clone)]
pub struct LambertianReflection {
    pub r: Spectrum,
}

impl LambertianReflection {
    pub fn new(r: Spectrum) -> Self {
        LambertianReflection { r }
    }
    pub fn f(&self, _wo: &Vector3f, _wi: &Vector3f) -> Spectrum {
        self.r * Spectrum::new(INV_PI)
    }
    pub fn sample_f(
        &self,
        wo: &Vector3f,
        wi: &mut Vector3f,
        u: &Point2f,
        pdf: &mut Float,
        _sampled_type: &mut u8,
    ) -> Spectrum {
        *wi = cosine_sample_hemisphere(u);
        if wo.z < 0.0 as Float {
            wi.z *= -1.0 as Float;
        }
        *pdf = self.pdf(wo, &*wi);
        self.f(wo, &*wi)
    }
    pub fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        if vec3_same_hemisphere_vec3(wo, wi) {
            abs_cos_theta(wi) * INV_PI
        } else {
            0.0 as Float
        }
    }
    pub fn get_type(&self) -> u8 {
        BxdfType::BsdfDiffuse as u8 | BxdfType::BsdfReflection as u8
    }
}

#[derive(Copy, Clone, Debug)]
pub struct LambertianTransmission {
    t: Spectrum,
}

impl LambertianTransmission {
    pub fn new(t: Spectrum) -> LambertianTransmission {
        LambertianTransmission { t }
    }
    pub fn f(&self, _wo: &Vector3f, _wi: &Vector3f) -> Spectrum {
        self.t * INV_PI
    }
    pub fn sample_f(
        &self,
        wo: &Vector3f,
        wi: &mut Vector3f,
        u: &Point2f,
        pdf: &mut Float,
        _sampled_type: &mut u8,
    ) -> Spectrum {
        *wi = cosine_sample_hemisphere(u);
        if wo.z > 0.0 as Float {
            wi.z *= -1.0 as Float;
        }
        *pdf = self.pdf(wo, &*wi);
        self.f(wo, &*wi)
    }
    pub fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        if !vec3_same_hemisphere_vec3(wo, wi) {
            abs_cos_theta(wi) * INV_PI
        } else {
            0.0 as Float
        }
    }
    pub fn get_type(&self) -> u8 {
        BxdfType::BsdfDiffuse as u8 | BxdfType::BsdfTransmission as u8
    }
}

pub struct OrenNayar {
    pub r: Spectrum,
    pub a: Float,
    pub b: Float,
}

impl OrenNayar {
    pub fn new(r: Spectrum, sigma: Float) -> Self {
        let sigma = radians(sigma);
        let sigma2: Float = sigma * sigma;
        OrenNayar {
            r,
            a: 1.0 - (sigma2 / (2.0 * (sigma2 + 0.33))),
            b: 0.45 * sigma2 / (sigma2 + 0.09),
        }
    }
    pub fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        let sin_theta_i: Float = sin_theta(wi);
        let sin_theta_o: Float = sin_theta(wo);
        // compute cosine term of Oren-Nayar model
        let mut max_cos: Float = 0.0 as Float;
        if sin_theta_i > 1.0e-4 && sin_theta_o > 1.0e-4 {
            let sin_phi_i: Float = sin_phi(wi);
            let cos_phi_i: Float = cos_phi(wi);
            let sin_phi_o: Float = sin_phi(wo);
            let cos_phi_o: Float = cos_phi(wo);
            let d_cos: Float = cos_phi_i * cos_phi_o + sin_phi_i * sin_phi_o;
            max_cos = d_cos.max(0.0 as Float);
        }
        // compute sine and tangent terms of Oren-Nayar model
        let sin_alpha: Float;
        let tan_beta: Float;
        if abs_cos_theta(wi) > abs_cos_theta(wo) {
            sin_alpha = sin_theta_o;
            tan_beta = sin_theta_i / abs_cos_theta(wi);
        } else {
            sin_alpha = sin_theta_i;
            tan_beta = sin_theta_o / abs_cos_theta(wo);
        }
        self.r * Spectrum::new(INV_PI * (self.a + self.b * max_cos * sin_alpha * tan_beta))
    }
    pub fn sample_f(
        &self,
        wo: &Vector3f,
        wi: &mut Vector3f,
        u: &Point2f,
        pdf: &mut Float,
        _sampled_type: &mut u8,
    ) -> Spectrum {
        *wi = cosine_sample_hemisphere(u);
        if wo.z < 0.0 as Float {
            wi.z *= -1.0 as Float;
        }
        *pdf = self.pdf(wo, &*wi);
        self.f(wo, &*wi)
    }
    pub fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        if vec3_same_hemisphere_vec3(wo, wi) {
            abs_cos_theta(wi) * INV_PI
        } else {
            0.0 as Float
        }
    }
    pub fn get_type(&self) -> u8 {
        BxdfType::BsdfDiffuse as u8 | BxdfType::BsdfReflection as u8
    }
}

pub struct MicrofacetReflection {
    pub r: Spectrum,
    pub distribution: Arc<dyn MicrofacetDistribution + Send + Sync>,
    pub fresnel: Fresnel,
}

impl MicrofacetReflection {
    pub fn new(
        r: Spectrum,
        distribution: Arc<dyn MicrofacetDistribution + Send + Sync>,
        fresnel: Fresnel,
    ) -> Self {
        MicrofacetReflection {
            r,
            distribution,
            fresnel,
        }
    }
    pub fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        let cos_theta_o: Float = abs_cos_theta(wo);
        let cos_theta_i: Float = abs_cos_theta(wi);
        let mut wh: Vector3f = *wi + *wo;
        // handle degenerate cases for microfacet reflection
        if cos_theta_i == 0.0 || cos_theta_o == 0.0 {
            return Spectrum::new(0.0);
        }
        if wh.x == 0.0 && wh.y == 0.0 && wh.z == 0.0 {
            return Spectrum::new(0.0);
        }
        wh = wh.normalize();
        let dot: Float = vec3_dot_vec3(wi, &wh);
        let f: Spectrum = self.fresnel.evaluate(dot);
        self.r * self.distribution.d(&wh) * self.distribution.g(wo, wi) * f
            / (4.0 as Float * cos_theta_i * cos_theta_o)
    }

    pub fn sample_f(
        &self,
        wo: &Vector3f,
        wi: &mut Vector3f,
        u: &Point2f,
        pdf: &mut Float,
        _sampled_type: &mut u8,
    ) -> Spectrum {
        // sample microfacet orientation $\wh$ and reflected direction $\wi$
        if wo.z == 0.0 as Float {
            return Spectrum::default();
        }
        let wh: Vector3f = self.distribution.sample_wh(wo, u);
        *wi = reflect(wo, &wh);
        if !vec3_same_hemisphere_vec3(wo, &*wi) {
            return Spectrum::default();
        }
        // compute PDF of _wi_ for microfacet reflection
        *pdf = self.distribution.pdf(wo, &wh) / (4.0 * vec3_dot_vec3(wo, &wh));
        self.f(wo, &*wi)
    }

    pub fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        if !vec3_same_hemisphere_vec3(wo, wi) {
            return 0.0 as Float;
        }
        let wh: Vector3f = (*wo + *wi).normalize();
        self.distribution.pdf(wo, &wh) / (4.0 * vec3_dot_vec3(wo, &wh))
    }

    pub fn get_type(&self) -> u8 {
        BxdfType::BsdfReflection as u8 | BxdfType::BsdfGlossy as u8
    }
}

// MicrofacetTransmission
pub struct MicrofacetTransmission {
    t: Spectrum,
    distribution: Arc<dyn MicrofacetDistribution + Send + Sync>,
    eta_a: Float,
    eta_b: Float,
    fresnel: FresnelDielectric,
    mode: TransportMode,
}

impl MicrofacetTransmission {
    pub fn new(
        t: Spectrum,
        distribution: Arc<dyn MicrofacetDistribution + Send + Sync>,
        eta_a: Float,
        eta_b: Float,
        mode: TransportMode,
    ) -> MicrofacetTransmission {
        MicrofacetTransmission {
            t,
            distribution,
            eta_a,
            eta_b,
            fresnel: FresnelDielectric {
                eta_i: eta_a,
                eta_t: eta_b,
            },
            mode,
        }
    }
    pub fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        if vec3_same_hemisphere_vec3(wo, wi) {
            // transmission only
            return Spectrum::zero();
        }

        let cos_theta_o = cos_theta(wo);
        let cos_theta_i = cos_theta(wi);
        // Handle degenerate case for microfacet reflection
        if cos_theta_o == 0.0 || cos_theta_i == 0.0 {
            return Spectrum::zero();
        }

        let eta = if cos_theta_o > 0.0 {
            self.eta_b / self.eta_a
        } else {
            self.eta_a / self.eta_b
        };

        let mut wh: Vector3f = (*wo + *wi * eta).normalize();
        if wh.z < 0.0 {
            wh = -wh;
        }

        let f = self.fresnel.evaluate(vec3_dot_vec3(wo, &wh));

        let sqrt_denom = vec3_dot_vec3(wo, &wh) + eta * vec3_dot_vec3(wi, &wh);
        let factor = match self.mode {
            TransportMode::Radiance => 1.0 / eta,
            _ => 1.0,
        };

        (Spectrum::new(1.0) - f)
            * self.t
            * Float::abs(
                self.distribution.d(&wh)
                    * self.distribution.g(wo, wi)
                    * eta
                    * eta
                    * vec3_abs_dot_vec3(wi, &wh)
                    * vec3_abs_dot_vec3(wo, &wh)
                    * factor
                    * factor
                    / (cos_theta_i * cos_theta_o * sqrt_denom * sqrt_denom),
            )
    }
    pub fn get_type(&self) -> u8 {
        BxdfType::BsdfTransmission as u8 | BxdfType::BsdfGlossy as u8
    }
    /// Override sample_f() to use a better importance sampling method than weighted cosine based
    /// on the microface distribution
    pub fn sample_f(
        &self,
        wo: &Vector3f,
        wi: &mut Vector3f,
        u: &Point2f,
        pdf: &mut Float,
        _sampled_type: &mut u8,
    ) -> Spectrum {
        if wo.z == 0.0 {
            return Spectrum::zero();
        }

        let wh: Vector3f = self.distribution.sample_wh(wo, u);
        let eta = if cos_theta(wo) > 0.0 {
            self.eta_a / self.eta_b
        } else {
            self.eta_b / self.eta_a
        };

        if refract(wo, &wh.into(), eta, wi) {
            *pdf = self.pdf(wo, &wi);
            self.f(wo, wi)
        } else {
            Spectrum::zero()
        }
    }
    pub fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        if vec3_same_hemisphere_vec3(wo, wi) {
            return 0.0;
        }

        let eta = if cos_theta(wo) > 0.0 {
            self.eta_b / self.eta_a
        } else {
            self.eta_a / self.eta_b
        };
        let wh: Vector3f = (*wo + *wi * eta).normalize();

        let sqrt_denom = vec3_dot_vec3(wo, &wh) + eta * vec3_dot_vec3(wi, &wh);
        let dwh_dwi = ((eta * eta * vec3_dot_vec3(wi, &wh)) / (sqrt_denom * sqrt_denom)).abs();

        self.distribution.pdf(wo, &wh) * dwh_dwi
    }
}

pub struct FresnelBlend {
    pub rd: Spectrum,
    pub rs: Spectrum,
    pub distribution: Option<TrowbridgeReitzDistribution>, // TODO: MicrofacetDistribution,
}

impl FresnelBlend {
    pub fn new(
        rd: Spectrum,
        rs: Spectrum,
        distribution: Option<TrowbridgeReitzDistribution>,
    ) -> Self {
        FresnelBlend {
            rd,
            rs,
            distribution,
        }
    }
    pub fn schlick_fresnel(&self, cos_theta: Float) -> Spectrum {
        self.rs + (Spectrum::new(1.0) - self.rs) * pow5(1.0 - cos_theta)
    }
    pub fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        let diffuse: Spectrum = self.rd
            * (Spectrum::new(1.0 as Float) - self.rs)
            * (28.0 as Float / (23.0 as Float * PI))
            * (1.0 - pow5(1.0 - 0.5 * abs_cos_theta(wi)))
            * (1.0 - pow5(1.0 - 0.5 * abs_cos_theta(wo)));
        let mut wh: Vector3f = *wi + *wo;
        if wh.x == 0.0 && wh.y == 0.0 && wh.z == 0.0 {
            return Spectrum::new(0.0 as Float);
        }
        wh = wh.normalize();
        if let Some(ref distribution) = self.distribution {
            let schlick_fresnel: Spectrum = self.schlick_fresnel(vec3_dot_vec3(wi, &wh));
            assert!(schlick_fresnel.c[0] >= 0.0, "wi = {:?}; wh = {:?}", wi, wh);
            let specular: Spectrum = schlick_fresnel
                * (distribution.d(&wh)
                    / (4.0
                        * vec3_dot_vec3(wi, &wh).abs()
                        * f32::max(abs_cos_theta(wi), abs_cos_theta(wo))));
            diffuse + specular
        } else {
            diffuse
        }
    }
    pub fn sample_f(
        &self,
        wo: &Vector3f,
        wi: &mut Vector3f,
        sample: &Point2f,
        pdf: &mut Float,
        _sampled_type: &mut u8,
    ) -> Spectrum {
        let mut u: Point2f = *sample;
        if u[0] < 0.5 as Float {
            u[0] = Float::min(2.0 * u[0], FLOAT_ONE_MINUS_EPSILON);
            // cosine-sample the hemisphere, flipping the direction if necessary
            *wi = cosine_sample_hemisphere(&u);
            if wo.z < 0.0 as Float {
                wi.z *= -1.0 as Float;
            }
        } else {
            u[0] = Float::min(2.0 * (u[0] - 0.5 as Float), FLOAT_ONE_MINUS_EPSILON);
            // sample microfacet orientation $\wh$ and reflected direction $\wi$
            if let Some(ref distribution) = self.distribution {
                let wh: Vector3f = distribution.sample_wh(wo, &u);
                *wi = reflect(wo, &wh);
                if !vec3_same_hemisphere_vec3(wo, &*wi) {
                    return Spectrum::new(0.0);
                }
            }
        }
        *pdf = self.pdf(wo, &*wi);
        self.f(wo, &*wi)
    }
    pub fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        // if (!SameHemisphere(wo, wi)) return 0;
        if !vec3_same_hemisphere_vec3(wo, wi) {
            return 0.0 as Float;
        }
        let wh: Vector3f = (*wo + *wi).normalize();
        if let Some(ref distribution) = self.distribution {
            let pdf_wh: Float = distribution.pdf(wo, &wh);
            0.5 as Float * (abs_cos_theta(wi) * INV_PI + pdf_wh / (4.0 * vec3_dot_vec3(wo, &wh)))
        } else {
            0.0 as Float
        }
    }
    pub fn get_type(&self) -> u8 {
        BxdfType::BsdfReflection as u8 | BxdfType::BsdfGlossy as u8
    }
}

pub struct FourierBSDF {
    pub bsdf_table: Arc<FourierBSDFTable>,
    pub mode: TransportMode,
}

impl FourierBSDF {
    pub fn new(bsdf_table: Arc<FourierBSDFTable>, mode: TransportMode) -> Self {
        FourierBSDF { bsdf_table, mode }
    }
    pub fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        // find the zenith angle cosines and azimuth difference angle
        let mu_i: Float = cos_theta(&-(*wi));
        let mu_o: Float = cos_theta(wo);
        let cos_phi: Float = cos_d_phi(&-(*wi), wo);
        // compute Fourier coefficients

        // determine offsets and weights
        let mut offset_i: i32 = 0;
        let mut offset_o: i32 = 0;
        let mut weights_i: [Float; 4] = [0.0 as Float; 4];
        let mut weights_o: [Float; 4] = [0.0 as Float; 4];
        if !self
            .bsdf_table
            .get_weights_and_offset(mu_i, &mut offset_i, &mut weights_i)
            || !self
                .bsdf_table
                .get_weights_and_offset(mu_o, &mut offset_o, &mut weights_o)
        {
            return Spectrum::default();
        }
        // allocate storage to accumulate _ak_ coefficients
        let mut ak: Vec<Float> =
            Vec::with_capacity((self.bsdf_table.m_max * self.bsdf_table.n_channels) as usize);
        for _i in 0..(self.bsdf_table.m_max * self.bsdf_table.n_channels) as usize {
            ak.push(0.0 as Float); // initialize with 0
        }
        // accumulate weighted sums of nearby $a_k$ coefficients
        let mut m_max: i32 = 0;
        for b in 0..4 {
            for a in 0..4 {
                // add contribution of _(a, b)_ to $a_k$ values
                let weight: Float = weights_i[a] * weights_o[b];
                if weight != 0.0 as Float {
                    let mut m: i32 = 0;
                    let a_idx: i32 =
                        self.bsdf_table
                            .get_ak(offset_i + a as i32, offset_o + b as i32, &mut m);
                    m_max = std::cmp::max(m_max, m);
                    for c in 0..self.bsdf_table.n_channels as usize {
                        for k in 0..m as usize {
                            ak[c * self.bsdf_table.m_max as usize + k] += weight
                                * self.bsdf_table.a[(a_idx + c as i32 * m + k as i32) as usize];
                        }
                    }
                }
            }
        }
        // evaluate Fourier expansion for angle $\phi$
        let y: Float = (0.0 as Float).max(fourier(&ak, 0_usize, m_max, cos_phi as f64));
        let mut scale: Float = 0.0 as Float;
        if mu_i != 0.0 as Float {
            scale = 1.0 as Float / mu_i.abs();
        }
        // update _scale_ to account for adjoint light transport
        if self.mode == TransportMode::Radiance && (mu_i * mu_o) > 0.0 as Float {
            let eta: Float;
            if mu_i > 0.0 as Float {
                eta = 1.0 as Float / self.bsdf_table.eta;
            } else {
                eta = self.bsdf_table.eta;
            }
            scale *= eta * eta;
        }
        if self.bsdf_table.n_channels == 1_i32 {
            Spectrum::new(y * scale)
        } else {
            // compute and return RGB colors for tabulated BSDF
            let r: Float = fourier(
                &ak,
                (1_i32 * self.bsdf_table.m_max) as usize,
                m_max,
                cos_phi as f64,
            );
            let b: Float = fourier(
                &ak,
                (2_i32 * self.bsdf_table.m_max) as usize,
                m_max,
                cos_phi as f64,
            );
            let g: Float = 1.39829 as Float * y - 0.100913 as Float * b - 0.297375 as Float * r;
            let rgb: [Float; 3] = [r * scale, g * scale, b * scale];
            Spectrum::from_rgb(&rgb).clamp(0.0 as Float, std::f32::INFINITY as Float)
        }
    }
    pub fn sample_f(
        &self,
        wo: &Vector3f,
        wi: &mut Vector3f,
        sample: &Point2f,
        pdf: &mut Float,
        _sampled_type: &mut u8,
    ) -> Spectrum {
        // sample zenith angle component for _FourierBSDF_
        let mu_o: Float = cos_theta(wo);
        let mut pdf_mu: Float = 0.0;
        let mu_i: Float = sample_catmull_rom_2d(
            &self.bsdf_table.mu,
            &self.bsdf_table.mu,
            &self.bsdf_table.a0,
            &self.bsdf_table.cdf,
            mu_o,
            sample[1],
            None,
            Some(&mut pdf_mu),
        );
        // compute Fourier coefficients $a_k$ for $(\mui, \muo)$

        // determine offsets and weights for $\mui$ and $\muo$
        let mut offset_i: i32 = 0;
        let mut offset_o: i32 = 0;
        let mut weights_i: [Float; 4] = [0.0 as Float; 4];
        let mut weights_o: [Float; 4] = [0.0 as Float; 4];
        if !self
            .bsdf_table
            .get_weights_and_offset(mu_i, &mut offset_i, &mut weights_i)
            || !self
                .bsdf_table
                .get_weights_and_offset(mu_o, &mut offset_o, &mut weights_o)
        {
            return Spectrum::default();
        }
        // allocate storage to accumulate _ak_ coefficients
        let mut ak: Vec<Float> =
            Vec::with_capacity((self.bsdf_table.m_max * self.bsdf_table.n_channels) as usize);
        for _i in 0..(self.bsdf_table.m_max * self.bsdf_table.n_channels) as usize {
            ak.push(0.0 as Float); // initialize with 0
        }
        // accumulate weighted sums of nearby $a_k$ coefficients
        let mut m_max: i32 = 0;
        for b in 0..4 {
            for a in 0..4 {
                // add contribution of _(a, b)_ to $a_k$ values
                let weight: Float = weights_i[a] * weights_o[b];
                if weight != 0.0 as Float {
                    let mut m: i32 = 0;
                    let a_idx =
                        self.bsdf_table
                            .get_ak(offset_i + a as i32, offset_o + b as i32, &mut m);
                    m_max = std::cmp::max(m_max, m);
                    for c in 0..self.bsdf_table.n_channels as usize {
                        for k in 0..m as usize {
                            ak[c * self.bsdf_table.m_max as usize + k] += weight
                                * self.bsdf_table.a[(a_idx + c as i32 * m + k as i32) as usize];
                        }
                    }
                }
            }
        }
        // importance sample the luminance Fourier expansion
        let mut phi: Float = 0.0;
        let mut pdf_phi: Float = 0.0;
        let y: Float = sample_fourier(
            &ak,
            &self.bsdf_table.recip,
            m_max,
            sample[0],
            &mut pdf_phi,
            &mut phi,
        );
        *pdf = (0.0 as Float).max(pdf_phi * pdf_mu);
        // compute the scattered direction for _FourierBSDF_
        let sin_2_theta_i: Float = (0.0 as Float).max(1.0 as Float - mu_i * mu_i);
        let mut norm: Float = (sin_2_theta_i / sin_2_theta(wo)).sqrt();
        if norm.is_infinite() {
            norm = 0.0;
        }
        let sin_phi: Float = phi.sin();
        let cos_phi: Float = phi.cos();
        *wi = -Vector3f {
            x: norm * (cos_phi * wo.x - sin_phi * wo.y),
            y: norm * (sin_phi * wo.x + cos_phi * wo.y),
            z: mu_i,
        };
        // Mathematically, wi will be normalized (if wo was). However,
        // in practice, floating-point rounding error can cause some
        // error to accumulate in the computed value of wi here. This
        // can be catastrophic: if the ray intersects an object with
        // the FourierBSDF again and the wo (based on such a wi) is
        // nearly perpendicular to the surface, then the wi computed
        // at the next intersection can end up being substantially
        // (like 4x) longer than normalized, which leads to all sorts
        // of errors, including negative spectral values. Therefore,
        // we normalize again here.
        *wi = wi.normalize();
        // evaluate remaining Fourier expansions for angle $\phi$
        let mut scale: Float = 0.0 as Float;
        if mu_i != 0.0 as Float {
            scale = 1.0 as Float / mu_i.abs();
        }
        // update _scale_ to account for adjoint light transport
        if self.mode == TransportMode::Radiance && (mu_i * mu_o) > 0.0 as Float {
            let eta: Float;
            if mu_i > 0.0 as Float {
                eta = 1.0 as Float / self.bsdf_table.eta;
            } else {
                eta = self.bsdf_table.eta;
            }
            scale *= eta * eta;
        }
        if self.bsdf_table.n_channels == 1_i32 {
            Spectrum::new(y * scale)
        } else {
            // compute and return RGB colors for tabulated BSDF
            let r: Float = fourier(
                &ak,
                (1_i32 * self.bsdf_table.m_max) as usize,
                m_max,
                cos_phi as f64,
            );
            let b: Float = fourier(
                &ak,
                (2_i32 * self.bsdf_table.m_max) as usize,
                m_max,
                cos_phi as f64,
            );
            let g: Float = 1.39829 as Float * y - 0.100913 as Float * b - 0.297375 as Float * r;
            let rgb: [Float; 3] = [r * scale, g * scale, b * scale];
            Spectrum::from_rgb(&rgb).clamp(0.0 as Float, std::f32::INFINITY as Float)
        }
    }
    pub fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        // find the zenith angle cosines and azimuth difference angle
        let mu_i: Float = cos_theta(&-(*wi));
        let mu_o: Float = cos_theta(wo);
        let cos_phi: Float = cos_d_phi(&-(*wi), wo);
        // compute luminance Fourier coefficients
        let mut offset_i: i32 = 0;
        let mut offset_o: i32 = 0;
        let mut weights_i: [Float; 4] = [0.0 as Float; 4];
        let mut weights_o: [Float; 4] = [0.0 as Float; 4];
        if !self
            .bsdf_table
            .get_weights_and_offset(mu_i, &mut offset_i, &mut weights_i)
            || !self
                .bsdf_table
                .get_weights_and_offset(mu_o, &mut offset_o, &mut weights_o)
        {
            return 0.0 as Float;
        }
        let mut ak: Vec<Float> = Vec::with_capacity(self.bsdf_table.m_max as usize);
        for _i in 0..self.bsdf_table.m_max as usize {
            ak.push(0.0 as Float); // initialize with 0
        }
        let mut m_max: i32 = 0;
        for o in 0..4 {
            for i in 0..4 {
                let weight: Float = weights_i[i] * weights_o[o];
                if weight == 0.0 as Float {
                    continue;
                }
                let mut order: i32 = 0;
                let a_idx: i32 =
                    self.bsdf_table
                        .get_ak(offset_i + i as i32, offset_o + o as i32, &mut order);
                m_max = std::cmp::max(m_max, order);
                for k in 0..order as usize {
                    ak[k] += weight * self.bsdf_table.a[(a_idx + k as i32) as usize];
                }
            }
        }
        // evaluate probability of sampling _wi_
        let mut rho: Float = 0.0;
        for o in 0..4 {
            if weights_o[o] == 0.0 as Float {
                continue;
            }
            rho += weights_o[o]
                * self.bsdf_table.cdf[(offset_o as usize + o) * self.bsdf_table.n_mu as usize
                    + self.bsdf_table.n_mu as usize
                    - 1 as usize]
                * (2.0 as Float * PI);
        }
        let y: Float = (0.0 as Float).max(fourier(&ak, 0_usize, m_max, cos_phi as f64));
        if rho > 0.0 as Float && y > 0.0 as Float {
            y / rho
        } else {
            0.0 as Float
        }
    }
    pub fn get_type(&self) -> u8 {
        BxdfType::BsdfReflection as u8
            | BxdfType::BsdfTransmission as u8
            | BxdfType::BsdfGlossy as u8
    }
}

/// Utility function to calculate cosine via spherical coordinates.
pub fn cos_theta(w: &Vector3f) -> Float {
    w.z
}

/// Utility function to calculate the square cosine via spherical
/// coordinates.
pub fn cos_2_theta(w: &Vector3f) -> Float {
    w.z * w.z
}

/// Utility function to calculate the absolute value of the cosine via
/// spherical coordinates.
pub fn abs_cos_theta(w: &Vector3f) -> Float {
    w.z.abs()
}

/// Utility function to calculate the square sine via spherical
/// coordinates.
pub fn sin_2_theta(w: &Vector3f) -> Float {
    (0.0 as Float).max(1.0 as Float - cos_2_theta(w))
}

/// Utility function to calculate sine via spherical coordinates.
pub fn sin_theta(w: &Vector3f) -> Float {
    sin_2_theta(w).sqrt()
}

/// Utility function to calculate the tangent via spherical
/// coordinates.
pub fn tan_theta(w: &Vector3f) -> Float {
    sin_theta(w) / cos_theta(w)
}

/// Utility function to calculate the square tangent via spherical
/// coordinates.
pub fn tan_2_theta(w: &Vector3f) -> Float {
    sin_2_theta(w) / cos_2_theta(w)
}

/// Utility function to calculate cosine via spherical coordinates.
pub fn cos_phi(w: &Vector3f) -> Float {
    let sin_theta: Float = sin_theta(w);
    if sin_theta == 0.0 as Float {
        1.0 as Float
    } else {
        clamp_t(w.x / sin_theta, -1.0, 1.0)
    }
}

/// Utility function to calculate sine via spherical coordinates.
pub fn sin_phi(w: &Vector3f) -> Float {
    let sin_theta: Float = sin_theta(w);
    if sin_theta == 0.0 as Float {
        0.0 as Float
    } else {
        clamp_t(w.y / sin_theta, -1.0, 1.0)
    }
}

/// Utility function to calculate square cosine via spherical coordinates.
pub fn cos_2_phi(w: &Vector3f) -> Float {
    cos_phi(w) * cos_phi(w)
}

/// Utility function to calculate square sine via spherical coordinates.
pub fn sin_2_phi(w: &Vector3f) -> Float {
    sin_phi(w) * sin_phi(w)
}

/// Utility function to calculate the cosine of the angle between two
/// vectors in the shading coordinate system.
pub fn cos_d_phi(wa: &Vector3f, wb: &Vector3f) -> Float {
    clamp_t(
        (wa.x * wb.x + wa.y * wb.y)
            / ((wa.x * wa.x + wa.y * wa.y) * (wb.x * wb.x + wb.y * wb.y)).sqrt(),
        -1.0 as Float,
        1.0 as Float,
    )
}

/// Computes the reflection direction given an incident direction and
/// a surface normal.
pub fn reflect(wo: &Vector3f, n: &Vector3f) -> Vector3f {
    -(*wo) + *n * 2.0 as Float * vec3_dot_vec3(wo, n)
}

/// Computes the refraction direction given an incident direction, a
/// surface normal, and the ratio of indices of refraction (incident
/// and transmitted).
pub fn refract(wi: &Vector3f, n: &Normal3f, eta: Float, wt: &mut Vector3f) -> bool {
    // compute $\cos \theta_\roman{t}$ using Snell's law
    let cos_theta_i: Float = nrm_dot_vec3(n, wi);
    let sin2_theta_i: Float = (0.0 as Float).max(1.0 as Float - cos_theta_i * cos_theta_i);
    let sin2_theta_t: Float = eta * eta * sin2_theta_i;
    // handle total internal reflection for transmission
    if sin2_theta_t >= 1.0 as Float {
        return false;
    }
    let cos_theta_t: Float = (1.0 as Float - sin2_theta_t).sqrt();
    *wt = -(*wi) * eta + Vector3f::from(*n) * (eta * cos_theta_i - cos_theta_t);
    true
}

/// Check that two vectors lie on the same side of of the surface.
pub fn vec3_same_hemisphere_vec3(w: &Vector3f, wp: &Vector3f) -> bool {
    w.z * wp.z > 0.0 as Float
}

// see reflection.cpp

/// Computes the Fresnel reflection formula for dielectric materials
/// and unpolarized light.
pub fn fr_dielectric(cos_theta_i: Float, eta_i: Float, eta_t: Float) -> Float {
    let mut cos_theta_i = clamp_t(cos_theta_i, -1.0, 1.0);
    // potentially swap indices of refraction
    let entering: bool = cos_theta_i > 0.0;
    // use local copies because of potential swap (otherwise eta_i and
    // eta_t would have to be mutable)
    let mut local_eta_i = eta_i;
    let mut local_eta_t = eta_t;
    if !entering {
        std::mem::swap(&mut local_eta_i, &mut local_eta_t);
        cos_theta_i = cos_theta_i.abs();
    }
    // compute _cos_theta_t_ using Snell's law
    let sin_theta_i: Float = (0.0 as Float)
        .max(1.0 as Float - cos_theta_i * cos_theta_i)
        .sqrt();
    let sin_theta_t: Float = local_eta_i / local_eta_t * sin_theta_i;
    // handle total internal reflection
    if sin_theta_t >= 1.0 as Float {
        return 1.0 as Float;
    }
    let cos_theta_t: Float = (0.0 as Float)
        .max(1.0 as Float - sin_theta_t * sin_theta_t)
        .sqrt();
    let r_parl: Float = ((local_eta_t * cos_theta_i) - (local_eta_i * cos_theta_t))
        / ((local_eta_t * cos_theta_i) + (local_eta_i * cos_theta_t));
    let r_perp: Float = ((local_eta_i * cos_theta_i) - (local_eta_t * cos_theta_t))
        / ((local_eta_i * cos_theta_i) + (local_eta_t * cos_theta_t));
    (r_parl * r_parl + r_perp * r_perp) / 2.0
}

/// Computes the Fresnel reflectance at the boundary between a
/// conductor and a dielectric medium.
pub fn fr_conductor(cos_theta_i: Float, eta_i: Spectrum, eta_t: Spectrum, k: Spectrum) -> Spectrum {
    let not_clamped: Float = cos_theta_i;
    let cos_theta_i: Float = clamp_t(not_clamped, -1.0, 1.0);
    let eta: Spectrum = eta_t / eta_i;
    let eta_k: Spectrum = k / eta_i;
    let cos_theta_i2: Float = cos_theta_i * cos_theta_i;
    let sin_theta_i2: Float = 1.0 as Float - cos_theta_i2;
    let eta_2: Spectrum = eta * eta;
    let eta_k2: Spectrum = eta_k * eta_k;
    let t0: Spectrum = eta_2 - eta_k2 - Spectrum::new(sin_theta_i2);
    let a2_plus_b2: Spectrum = (t0 * t0 + eta_2 * eta_k2 * Spectrum::new(4 as Float)).sqrt();
    let t1: Spectrum = a2_plus_b2 + Spectrum::new(cos_theta_i2);
    let a: Spectrum = ((a2_plus_b2 + t0) * 0.5 as Float).sqrt();
    let t2: Spectrum = a * 2.0 as Float * cos_theta_i;
    let rs: Spectrum = (t1 - t2) / (t1 + t2);
    let t3: Spectrum = a2_plus_b2 * cos_theta_i2 + Spectrum::new(sin_theta_i2 * sin_theta_i2);
    let t4: Spectrum = t2 * sin_theta_i2;
    let rp: Spectrum = rs * (t3 - t4) / (t3 + t4);
    (rp + rs) * Spectrum::new(0.5 as Float)
}

#[inline]
fn pow5(v: Float) -> Float {
    (v * v) * (v * v) * v
}
