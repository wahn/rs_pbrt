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
// pbrt
use core::geometry::{
    nrm_cross_vec3, nrm_dot_vec3, nrm_faceforward_vec3, vec3_dot_nrm, vec3_dot_vec3, vec3_normalize,
};
use core::geometry::{Normal3f, Point2f, Vector3f};
use core::interaction::SurfaceInteraction;
use core::interpolation::{catmull_rom_weights, fourier};
use core::material::TransportMode;
use core::microfacet::{MicrofacetDistribution, TrowbridgeReitzDistribution};
use core::pbrt::INV_PI;
use core::pbrt::{clamp_t, radians};
use core::pbrt::{Float, Spectrum};
use core::rng::FLOAT_ONE_MINUS_EPSILON;
use core::sampling::cosine_sample_hemisphere;

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
                    println!("WORK: flags = {:?}", flags);
                    self.n_mu = buffer[1];
                    println!("WORK: n_mu = {:?}", self.n_mu);
                    let n_coeffs: i32 = buffer[2];
                    println!("WORK: n_coeffs = {:?}", n_coeffs);
                    self.m_max = buffer[3];
                    println!("WORK: m_max = {:?}", self.m_max);
                    self.n_channels = buffer[4];
                    println!("WORK: n_channels = {:?}", self.n_channels);
                    let n_bases: i32 = buffer[5];
                    println!("WORK: n_bases = {:?}", n_bases);
                    let mut buffer: [f32; 1] = [0_f32; 1]; // 1 32-bit float
                    let io_result = file.read_f32_into::<LittleEndian>(&mut buffer);
                    if io_result.is_ok() {
                        self.eta = buffer[0];
                        println!("WORK: eta = {:?}", self.eta);
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
                            println!("WORK: {} f32 values read for mu", self.mu.len());
                            // self.cdf
                            self.cdf
                                .reserve_exact(self.n_mu as usize * self.n_mu as usize);
                            for _ in 0..(self.n_mu as usize * self.n_mu as usize) {
                                let f: f32 = file.read_f32::<LittleEndian>().unwrap();
                                self.cdf.push(f as Float);
                            }
                            println!("WORK: {} f32 values read for cdf", self.cdf.len());
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
                            println!(
                                "WORK: {} f32 values read for offset_and_length",
                                offset_and_length.len()
                            );
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
                            println!("WORK: {} f32 values read for a", self.a.len());
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
    pub fn get_ak(&self, offset_i: i32, offset_o: i32, mptr: &mut i32) -> usize {
        *mptr = self.m[(offset_o * self.n_mu + offset_i) as usize];
        self.a_offset[(offset_o * self.n_mu + offset_i) as usize] as usize
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
    pub bxdfs: Vec<Arc<Bxdf + Sync + Send>>,
}

impl Bsdf {
    pub fn new(si: &SurfaceInteraction, eta: Float, bxdfs: Vec<Arc<Bxdf + Sync + Send>>) -> Self {
        let ss = vec3_normalize(&si.shading.dpdu);
        Bsdf {
            eta: eta,
            ns: si.shading.n,
            ng: si.n,
            ss: ss,
            ts: nrm_cross_vec3(&si.shading.n, &ss),
            bxdfs: bxdfs,
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
        let mut bxdf: Option<&Arc<Bxdf + Sync + Send>> = None;
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
            if (bxdf.get_type() & BxdfType::BsdfSpecular as u8 == 0_u8) && matching_comps > 1_u8 {
                let reflect: bool = vec3_dot_nrm(&*wi_world, &self.ng)
                    * vec3_dot_nrm(wo_world, &self.ng)
                    > 0.0 as Float;
                f = Spectrum::default();
                for i in 0..n_bxdfs {
                    if self.bxdfs[i].matches_flags(bsdf_flags)
                        && ((reflect && (bxdf.get_type() & BxdfType::BsdfReflection as u8) != 0_u8)
                            || (!reflect
                                && (bxdf.get_type() & BxdfType::BsdfTransmission as u8) != 0_u8))
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

pub trait Bxdf {
    fn matches_flags(&self, t: u8) -> bool {
        self.get_type() & t == self.get_type()
    }
    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum;
    fn sample_f(
        &self,
        wo: &Vector3f,
        wi: &mut Vector3f,
        u: &Point2f,
        pdf: &mut Float,
        sampled_type: &mut u8,
    ) -> Spectrum;
    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float;
    fn get_type(&self) -> u8;
}

pub struct ScaledBxDF {
    pub bxdf: Arc<Bxdf + Sync + Send>,
    pub scale: Spectrum,
}

impl ScaledBxDF {
    pub fn new(bxdf: Arc<Bxdf + Send + Sync>, scale: Spectrum) -> Self {
        ScaledBxDF {
            bxdf: bxdf,
            scale: scale,
        }
    }
}

impl Bxdf for ScaledBxDF {
    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        self.scale * self.bxdf.f(wo, wi)
    }
    fn sample_f(
        &self,
        wo: &Vector3f,
        wi: &mut Vector3f,
        sample: &Point2f,
        pdf: &mut Float,
        sampled_type: &mut u8,
    ) -> Spectrum {
        let f: Spectrum = self.bxdf.sample_f(wo, wi, sample, pdf, sampled_type);
        self.scale * f
    }
    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        self.bxdf.pdf(wo, wi)
    }
    fn get_type(&self) -> u8 {
        self.bxdf.get_type()
    }
}

pub trait Fresnel {
    fn evaluate(&self, cos_theta_i: &mut Float) -> Spectrum;
}

#[derive(Debug, Default, Copy, Clone)]
pub struct FresnelConductor {
    pub eta_i: Spectrum,
    pub eta_t: Spectrum,
    pub k: Spectrum,
}

impl Fresnel for FresnelConductor {
    fn evaluate(&self, cos_theta_i: &mut Float) -> Spectrum {
        fr_conductor(cos_theta_i, self.eta_i, self.eta_t, self.k)
    }
}

#[derive(Debug, Default, Copy, Clone)]
pub struct FresnelDielectric {
    pub eta_i: Float,
    pub eta_t: Float,
}

impl Fresnel for FresnelDielectric {
    fn evaluate(&self, cos_theta_i: &mut Float) -> Spectrum {
        Spectrum::new(fr_dielectric(cos_theta_i, self.eta_i, self.eta_t))
    }
}

#[derive(Debug, Default, Copy, Clone)]
pub struct FresnelNoOp {}

impl Fresnel for FresnelNoOp {
    fn evaluate(&self, _cos_theta_i: &mut Float) -> Spectrum {
        Spectrum::new(1.0 as Float)
    }
}

#[derive(Clone)]
pub struct SpecularReflection {
    pub r: Spectrum,
    pub fresnel: Arc<Fresnel + Send + Sync>,
}

impl SpecularReflection {
    pub fn new(r: Spectrum, fresnel: Arc<Fresnel + Send + Sync>) -> Self {
        SpecularReflection {
            r: r,
            fresnel: fresnel,
        }
    }
}

impl Bxdf for SpecularReflection {
    fn f(&self, _wo: &Vector3f, _wi: &Vector3f) -> Spectrum {
        Spectrum::new(0.0 as Float)
    }
    fn sample_f(
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
        let mut cos_theta_i: Float = cos_theta(&*wi);
        self.fresnel.evaluate(&mut cos_theta_i) * self.r / abs_cos_theta(&*wi)
    }
    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        if vec3_same_hemisphere_vec3(wo, wi) {
            abs_cos_theta(wi) * INV_PI
        } else {
            0.0 as Float
        }
    }
    fn get_type(&self) -> u8 {
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
            t: t,
            eta_a: eta_a,
            eta_b: eta_b,
            fresnel: FresnelDielectric {
                eta_i: eta_a,
                eta_t: eta_b,
            },
            mode: mode,
        }
    }
}

impl Bxdf for SpecularTransmission {
    fn f(&self, _wo: &Vector3f, _wi: &Vector3f) -> Spectrum {
        Spectrum::new(0.0 as Float)
    }
    fn sample_f(
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
            self.t * (Spectrum::new(1.0 as Float) - self.fresnel.evaluate(&mut cos_theta(&*wi)));
        // account for non-symmetry with transmission to different medium
        if self.mode == TransportMode::Radiance {
            ft *= Spectrum::new((eta_i * eta_i) / (eta_t * eta_t));
        }
        ft / abs_cos_theta(&*wi)
    }
    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        if vec3_same_hemisphere_vec3(wo, wi) {
            abs_cos_theta(wi) * INV_PI
        } else {
            0.0 as Float
        }
    }
    fn get_type(&self) -> u8 {
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
            r: r,
            t: t,
            eta_a: eta_a,
            eta_b: eta_b,
            mode: mode,
        }
    }
}

impl Bxdf for FresnelSpecular {
    fn f(&self, _wo: &Vector3f, _wi: &Vector3f) -> Spectrum {
        Spectrum::new(0.0 as Float)
    }
    fn sample_f(
        &self,
        wo: &Vector3f,
        wi: &mut Vector3f,
        sample: &Point2f,
        pdf: &mut Float,
        sampled_type: &mut u8,
    ) -> Spectrum {
        let mut ct: Float = cos_theta(wo);
        let f: Float = fr_dielectric(&mut ct, self.eta_a, self.eta_b);
        if sample[0] < f {
            // compute specular reflection for _FresnelSpecular_

            // compute perfect specular reflection direction
            *wi = Vector3f {
                x: -wo.x,
                y: -wo.y,
                z: wo.z,
            };
            if *sampled_type != 0_u8 {
                *sampled_type = self.get_type();
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
                *sampled_type = self.get_type();
            }
            *pdf = 1.0 as Float - f;
            return ft / abs_cos_theta(&*wi);
        }
    }
    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        if vec3_same_hemisphere_vec3(wo, wi) {
            abs_cos_theta(wi) * INV_PI
        } else {
            0.0 as Float
        }
    }
    fn get_type(&self) -> u8 {
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
        LambertianReflection { r: r }
    }
}

impl Bxdf for LambertianReflection {
    fn f(&self, _wo: &Vector3f, _wi: &Vector3f) -> Spectrum {
        self.r * Spectrum::new(INV_PI)
    }
    fn sample_f(
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
    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        if vec3_same_hemisphere_vec3(wo, wi) {
            abs_cos_theta(wi) * INV_PI
        } else {
            0.0 as Float
        }
    }
    fn get_type(&self) -> u8 {
        BxdfType::BsdfDiffuse as u8 | BxdfType::BsdfReflection as u8
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
            r: r,
            a: 1.0 - (sigma2 / (2.0 * (sigma2 + 0.33))),
            b: 0.45 * sigma2 / (sigma2 + 0.09),
        }
    }
}

impl Bxdf for OrenNayar {
    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
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
    fn sample_f(
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
    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        if vec3_same_hemisphere_vec3(wo, wi) {
            abs_cos_theta(wi) * INV_PI
        } else {
            0.0 as Float
        }
    }
    fn get_type(&self) -> u8 {
        BxdfType::BsdfDiffuse as u8 | BxdfType::BsdfReflection as u8
    }
}

pub struct MicrofacetReflection {
    pub r: Spectrum,
    pub distribution: Option<TrowbridgeReitzDistribution>, // TODO: MicrofacetDistribution,
    pub fresnel: Arc<Fresnel + Send + Sync>,
}

impl MicrofacetReflection {
    pub fn new(
        r: Spectrum,
        distribution: Option<TrowbridgeReitzDistribution>,
        fresnel: Arc<Fresnel + Send + Sync>,
    ) -> Self {
        MicrofacetReflection {
            r: r,
            distribution: distribution,
            fresnel: fresnel,
        }
    }
}

impl Bxdf for MicrofacetReflection {
    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
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
        wh = vec3_normalize(&wh);
        let mut dot: Float = vec3_dot_vec3(wi, &wh);
        let f: Spectrum = self.fresnel.evaluate(&mut dot);
        if let Some(ref distribution) = self.distribution {
            return self.r * distribution.d(&wh) * distribution.g(wo, wi) * f
                / (4.0 as Float * cos_theta_i * cos_theta_o);
        } else {
            panic!("MicrofacetReflection::f() needs self.distribution");
        }
    }
    fn sample_f(
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
        if let Some(ref distribution) = self.distribution {
            let wh: Vector3f = distribution.sample_wh(wo, u);
            *wi = reflect(wo, &wh);
            if !vec3_same_hemisphere_vec3(wo, &*wi) {
                return Spectrum::default();
            }
            // compute PDF of _wi_ for microfacet reflection
            *pdf = distribution.pdf(wo, &wh) / (4.0 * vec3_dot_vec3(wo, &wh));
            return self.f(wo, &*wi);
        } else {
            panic!("MicrofacetReflection::f() needs self.distribution");
        }
    }
    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        if !vec3_same_hemisphere_vec3(wo, wi) {
            return 0.0 as Float;
        }
        let wh: Vector3f = vec3_normalize(&(*wo + *wi));
        if let Some(ref distribution) = self.distribution {
            return distribution.pdf(wo, &wh) / (4.0 * vec3_dot_vec3(wo, &wh));
        } else {
            panic!("MicrofacetReflection::f() needs self.distribution");
        }
    }
    fn get_type(&self) -> u8 {
        BxdfType::BsdfReflection as u8 | BxdfType::BsdfGlossy as u8
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
            rd: rd,
            rs: rs,
            distribution: distribution,
        }
    }
    pub fn schlick_fresnel(&self, cos_theta: Float) -> Spectrum {
        self.rs + (Spectrum::new(1.0) - self.rs) * pow5(1.0 - cos_theta)
    }
}

impl Bxdf for FresnelBlend {
    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        let diffuse: Spectrum = self.rd
            * (Spectrum::new(1.0 as Float) - self.rs)
            * (28.0 as Float / (23.0 as Float * PI))
            * (1.0 - pow5(1.0 - 0.5 * abs_cos_theta(wi)))
            * (1.0 - pow5(1.0 - 0.5 * abs_cos_theta(wo)));
        let mut wh: Vector3f = *wi + *wo;
        if wh.x == 0.0 && wh.y == 0.0 && wh.z == 0.0 {
            return Spectrum::new(0.0 as Float);
        }
        wh = vec3_normalize(&wh);
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
    fn sample_f(
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
    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        // if (!SameHemisphere(wo, wi)) return 0;
        if !vec3_same_hemisphere_vec3(wo, wi) {
            return 0.0 as Float;
        }
        let wh: Vector3f = vec3_normalize(&(*wo + *wi));
        if let Some(ref distribution) = self.distribution {
            let pdf_wh: Float = distribution.pdf(wo, &wh);
            0.5 as Float * (abs_cos_theta(wi) * INV_PI + pdf_wh / (4.0 * vec3_dot_vec3(wo, &wh)))
        } else {
            0.0 as Float
        }
    }
    fn get_type(&self) -> u8 {
        BxdfType::BsdfReflection as u8 | BxdfType::BsdfGlossy as u8
    }
}

pub struct FourierBSDF {
    pub bsdf_table: Arc<FourierBSDFTable>,
    pub mode: TransportMode,
}

impl FourierBSDF {
    pub fn new(bsdf_table: Arc<FourierBSDFTable>, mode: TransportMode) -> Self {
        FourierBSDF {
            bsdf_table: bsdf_table,
            mode: mode,
        }
    }
}

impl Bxdf for FourierBSDF {
    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
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
                    let a_idx = self.bsdf_table.get_ak(offset_i, offset_o, &mut m);
                    m_max = std::cmp::max(m_max, m);
                    for c in 0..self.bsdf_table.n_channels as usize {
                        for k in 0..m as usize {
                            ak[c * self.bsdf_table.m_max as usize + k] +=
                                weight * self.bsdf_table.a[a_idx + c * m as usize + k];
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
            let r: Float = fourier(&ak, (1_i32 * self.bsdf_table.m_max) as usize, m_max, cos_phi as f64);
            let b: Float = fourier(&ak, (2_i32 * self.bsdf_table.m_max) as usize, m_max, cos_phi as f64);
            let g: Float = 1.39829 as Float * y - 0.100913 as Float * b - 0.297375 as Float * r;
            let mut rgb: [Float; 3] = [r * scale, g * scale, b * scale];
            Spectrum::from_rgb(&rgb).clamp(0.0 as Float, std::f32::INFINITY as Float)
        }
    }
    fn sample_f(
        &self,
        wo: &Vector3f,
        wi: &mut Vector3f,
        sample: &Point2f,
        pdf: &mut Float,
        _sampled_type: &mut u8,
    ) -> Spectrum {
        // WORK
        Spectrum::default()
    }
    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        // WORK
        0.0 as Float
    }
    fn get_type(&self) -> u8 {
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
        ((wa.x * wb.x + wa.y * wb.y) / ((wa.x * wa.x + wa.y * wa.y) * (wb.x * wb.x + wb.y * wb.y)))
            .sqrt(),
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
pub fn fr_dielectric(cos_theta_i: &mut Float, eta_i: Float, eta_t: Float) -> Float {
    let not_clamped: Float = *cos_theta_i;
    *cos_theta_i = clamp_t(not_clamped, -1.0, 1.0);
    // potentially swap indices of refraction
    let entering: bool = *cos_theta_i > 0.0;
    // use local copies because of potential swap (otherwise eta_i and
    // eta_t would have to be mutable)
    let mut local_eta_i = eta_i;
    let mut local_eta_t = eta_t;
    if !entering {
        std::mem::swap(&mut local_eta_i, &mut local_eta_t);
        *cos_theta_i = (*cos_theta_i).abs();
    }
    // compute _cos_theta_t_ using Snell's law
    let sin_theta_i: Float = (0.0 as Float)
        .max(1.0 as Float - *cos_theta_i * *cos_theta_i)
        .sqrt();
    let sin_theta_t: Float = local_eta_i / local_eta_t * sin_theta_i;
    // handle total internal reflection
    if sin_theta_t >= 1.0 as Float {
        return 1.0 as Float;
    }
    let cos_theta_t: Float = (0.0 as Float)
        .max(1.0 as Float - sin_theta_t * sin_theta_t)
        .sqrt();
    let r_parl: Float = ((local_eta_t * *cos_theta_i) - (local_eta_i * cos_theta_t))
        / ((local_eta_t * *cos_theta_i) + (local_eta_i * cos_theta_t));
    let r_perp: Float = ((local_eta_i * *cos_theta_i) - (local_eta_t * cos_theta_t))
        / ((local_eta_i * *cos_theta_i) + (local_eta_t * cos_theta_t));
    (r_parl * r_parl + r_perp * r_perp) / 2.0
}

/// Computes the Fresnel reflectance at the boundary between a
/// conductor and a dielectric medium.
pub fn fr_conductor(
    cos_theta_i: &mut Float,
    eta_i: Spectrum,
    eta_t: Spectrum,
    k: Spectrum,
) -> Spectrum {
    let not_clamped: Float = *cos_theta_i;
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
