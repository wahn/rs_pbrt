// std
use std;
use std::sync::Arc;
// pbrt
use core::interaction::SurfaceInteraction;
use core::material::TransportMode;
use core::microfacet::{MicrofacetDistribution, TrowbridgeReitzDistribution};
use core::pbrt::INV_PI;
use core::pbrt::{Float, Spectrum};
use core::pbrt::{clamp_t, radians};
use core::rng::FLOAT_ONE_MINUS_EPSILON;
use core::sampling::cosine_sample_hemisphere;
use geometry::{Normal3f, Point2f, Vector3f};
use geometry::{nrm_cross_vec3, nrm_dot_vec3, nrm_faceforward_vec3, vec3_dot_nrm, vec3_dot_vec3,
               vec3_normalize};

// see reflection.h

pub struct Bsdf {
    pub eta: Float,
    /// shading normal
    pub ns: Normal3f,
    /// geometric normal
    pub ng: Normal3f,
    pub ss: Vector3f,
    pub ts: Vector3f,
    pub bxdfs: Vec<Box<Bxdf + Sync + Send>>,
}

impl Bsdf {
    pub fn new(si: &SurfaceInteraction, eta: Float, bxdfs: Vec<Box<Bxdf + Sync + Send>>) -> Self {
        let ss = vec3_normalize(si.shading.dpdu);
        Bsdf {
            eta: eta,
            ns: si.shading.n,
            ng: si.n,
            ss: ss,
            ts: nrm_cross_vec3(si.shading.n, ss),
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
    pub fn world_to_local(&self, v: Vector3f) -> Vector3f {
        Vector3f {
            x: vec3_dot_vec3(v, self.ss),
            y: vec3_dot_vec3(v, self.ts),
            z: vec3_dot_vec3(v, Vector3f::from(self.ns)),
        }
    }
    pub fn local_to_world(&self, v: Vector3f) -> Vector3f {
        Vector3f {
            x: self.ss.x * v.x + self.ts.x * v.y + self.ns.x * v.z,
            y: self.ss.y * v.x + self.ts.y * v.y + self.ns.y * v.z,
            z: self.ss.z * v.x + self.ts.z * v.y + self.ns.z * v.z,
        }
    }
    pub fn f(&self, wo_w: Vector3f, wi_w: Vector3f, flags: u8) -> Spectrum {
        // TODO: ProfilePhase pp(Prof::BSDFEvaluation);
        let wi: Vector3f = self.world_to_local(wi_w);
        let wo: Vector3f = self.world_to_local(wo_w);
        if wo.z == 0.0 as Float {
            return Spectrum::new(0.0 as Float);
        }
        let reflect: bool = (vec3_dot_vec3(wi_w, Vector3f::from(self.ng)) *
                             vec3_dot_vec3(wo_w, Vector3f::from(self.ng))) >
                            0.0 as Float;
        let mut f: Spectrum = Spectrum::new(0.0 as Float);
        let n_bxdfs: usize = self.bxdfs.len();
        for i in 0..n_bxdfs {
            if self.bxdfs[i].matches_flags(flags) &&
               ((reflect && (self.bxdfs[i].get_type() & BxdfType::BsdfReflection as u8 > 0_u8)) ||
                (!reflect &&
                 (self.bxdfs[i].get_type() & BxdfType::BsdfTransmission as u8 > 0_u8))) {
                f += self.bxdfs[i].f(wo, wi);
            }
        }
        f
    }
    /// Calls the individual Bxdf::sample_f() methods to generate samples.
    pub fn sample_f(&self,
                    wo_world: Vector3f,
                    wi_world: &mut Vector3f,
                    u: Point2f,
                    pdf: &mut Float,
                    bsdf_flags: u8,
                    sampled_type: &mut u8)
                    -> Spectrum {
        // TODO: ProfilePhase pp(Prof::BSDFSampling);
        // choose which _BxDF_ to sample
        let matching_comps: u8 = self.num_components(bsdf_flags);
        if matching_comps == 0 {
            *pdf = 0.0 as Float;
            *sampled_type = 0_u8;
            return Spectrum::default();
        }
        let comp: u8 = std::cmp::min((u[0] * matching_comps as Float).floor() as u8,
                                     matching_comps - 1_u8);
        // get _BxDF_ pointer for chosen component
        let mut bxdf: Option<&Box<Bxdf + Sync + Send>> = None;
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
            let mut f: Spectrum = bxdf.sample_f(wo, &mut wi, u_remapped, pdf, sampled_type);
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
            *wi_world = self.local_to_world(wi);
            // compute overall PDF with all matching _BxDF_s
            if (bxdf.get_type() & BxdfType::BsdfSpecular as u8 == 0_u8) && matching_comps > 1_u8 {
                for i in 0..n_bxdfs {
                    // instead of self.bxdfs[i] != bxdf we compare stored index
                    if bxdf_index != i && self.bxdfs[i].matches_flags(bsdf_flags) {
                        *pdf += self.bxdfs[i].pdf(wo, wi);
                    }
                }
            }
            if matching_comps > 1_u8 {
                *pdf /= matching_comps as Float;
            }
            // compute value of BSDF for sampled direction
            if (bxdf.get_type() & BxdfType::BsdfSpecular as u8 == 0_u8) && matching_comps > 1_u8 {
                let reflect: bool = vec3_dot_nrm(*wi_world, self.ng) *
                                    vec3_dot_nrm(wo_world, self.ng) >
                                    0.0 as Float;
                f = Spectrum::default();
                for i in 0..n_bxdfs {
                    if self.bxdfs[i].matches_flags(bsdf_flags) &&
                       ((reflect && (bxdf.get_type() & BxdfType::BsdfReflection as u8) != 0_u8) ||
                        (reflect && (bxdf.get_type() & BxdfType::BsdfTransmission as u8) == 0_u8)) {
                        f += self.bxdfs[i].f(wo, wi);
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
    pub fn pdf(&self, wo_world: Vector3f, wi_world: Vector3f, bsdf_flags: u8) -> Float {
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
                pdf += self.bxdfs[i].pdf(wo, wi);
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
    fn f(&self, wo: Vector3f, wi: Vector3f) -> Spectrum;
    fn sample_f(&self,
                wo: Vector3f,
                wi: &mut Vector3f,
                u: Point2f,
                pdf: &mut Float,
                sampled_type: &mut u8)
                -> Spectrum;
    fn pdf(&self, wo: Vector3f, wi: Vector3f) -> Float {
        if vec3_same_hemisphere_vec3(wo, wi) {
            abs_cos_theta(wi) * INV_PI
        } else {
            0.0 as Float
        }
    }
    fn get_type(&self) -> u8;
}

pub trait Fresnel {
    fn evaluate(&self, cos_theta_i: &mut Float) -> Spectrum;
}

#[derive(Debug,Default,Copy,Clone)]
pub struct FresnelDielectric {
    pub eta_i: Float,
    pub eta_t: Float,
}

impl Fresnel for FresnelDielectric {
    fn evaluate(&self, cos_theta_i: &mut Float) -> Spectrum {
        Spectrum::new(fr_dielectric(cos_theta_i, self.eta_i, self.eta_t))
    }
}

#[derive(Debug,Default,Copy,Clone)]
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
    fn f(&self, _wo: Vector3f, _wi: Vector3f) -> Spectrum {
        Spectrum::new(0.0 as Float)
    }
    fn sample_f(&self,
                wo: Vector3f,
                wi: &mut Vector3f,
                _sample: Point2f,
                pdf: &mut Float,
                _sampled_type: &mut u8)
                -> Spectrum {
        // compute perfect specular reflection direction
        *wi = Vector3f {
            x: -wo.x,
            y: -wo.y,
            z: wo.z,
        };
        *pdf = 1.0 as Float;
        let mut cos_theta_i: Float = cos_theta(*wi);
        self.fresnel.evaluate(&mut cos_theta_i) * self.r / abs_cos_theta(*wi)
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
    fn f(&self, _wo: Vector3f, _wi: Vector3f) -> Spectrum {
        Spectrum::new(0.0 as Float)
    }
    fn sample_f(&self,
                wo: Vector3f,
                wi: &mut Vector3f,
                _sample: Point2f,
                pdf: &mut Float,
                _sampled_type: &mut u8)
                -> Spectrum {
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
        if !refract(wo,
                    nrm_faceforward_vec3(Normal3f {
                                             x: 0.0,
                                             y: 0.0,
                                             z: 1.0,
                                         },
                                         wo),
                    eta_i / eta_t,
                    wi) {
            return Spectrum::default();
        }
        *pdf = 1.0;
        let mut ft: Spectrum =
            self.t * (Spectrum::new(1.0 as Float) - self.fresnel.evaluate(&mut cos_theta(*wi)));
        // account for non-symmetry with transmission to different medium
        if self.mode == TransportMode::Radiance {
            ft *= Spectrum::new((eta_i * eta_i) / (eta_t * eta_t));
        }
        ft / abs_cos_theta(*wi)
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
    fn f(&self, _wo: Vector3f, _wi: Vector3f) -> Spectrum {
        Spectrum::new(0.0 as Float)
    }
    fn sample_f(&self,
                wo: Vector3f,
                wi: &mut Vector3f,
                sample: Point2f,
                pdf: &mut Float,
                sampled_type: &mut u8)
                -> Spectrum {
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
            return self.r * f / abs_cos_theta(*wi);
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
            if !refract(wo,
                        nrm_faceforward_vec3(Normal3f {
                                                 x: 0.0,
                                                 y: 0.0,
                                                 z: 1.0,
                                             },
                                             wo),
                        eta_i / eta_t,
                        wi) {
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
            return ft / abs_cos_theta(*wi);
        }
    }
    fn get_type(&self) -> u8 {
        BxdfType::BsdfReflection as u8 | BxdfType::BsdfTransmission as u8 |
        BxdfType::BsdfSpecular as u8
    }
}

#[derive(Debug,Default,Copy,Clone)]
pub struct LambertianReflection {
    pub r: Spectrum,
}

impl LambertianReflection {
    pub fn new(r: Spectrum) -> Self {
        LambertianReflection { r: r }
    }
}

impl Bxdf for LambertianReflection {
    fn f(&self, _wo: Vector3f, _wi: Vector3f) -> Spectrum {
        self.r * Spectrum::new(INV_PI)
    }
    fn sample_f(&self,
                wo: Vector3f,
                wi: &mut Vector3f,
                u: Point2f,
                pdf: &mut Float,
                _sampled_type: &mut u8)
                -> Spectrum {
        *wi = cosine_sample_hemisphere(u);
        if wo.z < 0.0 as Float {
            wi.z *= -1.0 as Float;
        }
        *pdf = self.pdf(wo, *wi);
        self.f(wo, *wi)
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
    fn f(&self, wo: Vector3f, wi: Vector3f) -> Spectrum {
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
    fn sample_f(&self,
                wo: Vector3f,
                wi: &mut Vector3f,
                u: Point2f,
                pdf: &mut Float,
                _sampled_type: &mut u8)
                -> Spectrum {
        *wi = cosine_sample_hemisphere(u);
        if wo.z > 0.0 as Float {
            wi.z *= -1.0 as Float;
        }
        *pdf = self.pdf(wo, *wi);
        self.f(wo, *wi)
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
    pub fn new(r: Spectrum,
               distribution: Option<TrowbridgeReitzDistribution>,
               fresnel: Arc<Fresnel + Send + Sync>)
               -> Self {
        MicrofacetReflection {
            r: r,
            distribution: distribution,
            fresnel: fresnel,
        }
    }
}

impl Bxdf for MicrofacetReflection {
    fn f(&self, wo: Vector3f, wi: Vector3f) -> Spectrum {
        let cos_theta_o: Float = abs_cos_theta(wo);
        let cos_theta_i: Float = abs_cos_theta(wi);
        let mut wh: Vector3f = wi + wo;
        // handle degenerate cases for microfacet reflection
        if cos_theta_i == 0.0 || cos_theta_o == 0.0 {
            return Spectrum::new(0.0);
        }
        if wh.x == 0.0 && wh.y == 0.0 && wh.z == 0.0 {
            return Spectrum::new(0.0);
        }
        wh = vec3_normalize(wh);
        let mut dot: Float = vec3_dot_vec3(wi, wh);
        let f: Spectrum = self.fresnel.evaluate(&mut dot);
        if let Some(ref distribution) = self.distribution {
            return self.r * distribution.d(wh) * distribution.g(wo, wi) * f /
                   (4.0 as Float * cos_theta_i * cos_theta_o);
        } else {
            panic!("MicrofacetReflection::f() needs self.distribution");
        }
    }
    fn sample_f(&self,
                wo: Vector3f,
                wi: &mut Vector3f,
                u: Point2f,
                pdf: &mut Float,
                _sampled_type: &mut u8)
                -> Spectrum {
        // sample microfacet orientation $\wh$ and reflected direction $\wi$
        if wo.z == 0.0 as Float {
            return Spectrum::default();
        }
        if let Some(ref distribution) = self.distribution {
            let wh: Vector3f = distribution.sample_wh(wo, u);
            *wi = reflect(wo, wh);
            if !vec3_same_hemisphere_vec3(wo, *wi) {
                return Spectrum::default();
            }
            // compute PDF of _wi_ for microfacet reflection
            *pdf = distribution.pdf(wo, wh) / (4.0 * vec3_dot_vec3(wo, wh));
            return self.f(wo, *wi);
        } else {
            panic!("MicrofacetReflection::f() needs self.distribution");
        }
    }
    fn pdf(&self, wo: Vector3f, wi: Vector3f) -> Float {
        if !vec3_same_hemisphere_vec3(wo, wi) {
            return 0.0 as Float;
        }
        let wh: Vector3f = vec3_normalize(wo + wi);
        if let Some(ref distribution) = self.distribution {
            return distribution.pdf(wo, wh) / (4.0 * vec3_dot_vec3(wo, wh));
        } else {
            panic!("MicrofacetReflection::f() needs self.distribution");
        }
    }
    fn get_type(&self) -> u8 {
        BxdfType::BsdfReflection as u8 | BxdfType::BsdfGlossy as u8
    }
}

/// Utility function to calculate cosine via spherical coordinates.
pub fn cos_theta(w: Vector3f) -> Float {
    w.z
}

/// Utility function to calculate the square cosine via spherical
/// coordinates.
pub fn cos_2_theta(w: Vector3f) -> Float {
    w.z * w.z
}

/// Utility function to calculate the absolute value of the cosine via
/// spherical coordinates.
pub fn abs_cos_theta(w: Vector3f) -> Float {
    w.z.abs()
}

/// Utility function to calculate the square sine via spherical
/// coordinates.
pub fn sin_2_theta(w: Vector3f) -> Float {
    (0.0 as Float).max(1.0 as Float - cos_2_theta(w))
}

/// Utility function to calculate sine via spherical coordinates.
pub fn sin_theta(w: Vector3f) -> Float {
    sin_2_theta(w).sqrt()
}

/// Utility function to calculate the tangent via spherical
/// coordinates.
pub fn tan_theta(w: Vector3f) -> Float {
    sin_theta(w) / cos_theta(w)
}

/// Utility function to calculate the square tangent via spherical
/// coordinates.
pub fn tan_2_theta(w: Vector3f) -> Float {
    sin_2_theta(w) / cos_2_theta(w)
}

/// Utility function to calculate cosine via spherical coordinates.
pub fn cos_phi(w: Vector3f) -> Float {
    let sin_theta: Float = sin_theta(w);
    if sin_theta == 0.0 as Float {
        1.0 as Float
    } else {
        clamp_t(w.x / sin_theta, -1.0, 1.0)
    }
}

/// Utility function to calculate sine via spherical coordinates.
pub fn sin_phi(w: Vector3f) -> Float {
    let sin_theta: Float = sin_theta(w);
    if sin_theta == 0.0 as Float {
        0.0 as Float
    } else {
        clamp_t(w.y / sin_theta, -1.0, 1.0)
    }
}

/// Utility function to calculate square cosine via spherical coordinates.
pub fn cos_2_phi(w: Vector3f) -> Float {
    cos_phi(w) * cos_phi(w)
}

/// Utility function to calculate square sine via spherical coordinates.
pub fn sin_2_phi(w: Vector3f) -> Float {
    sin_phi(w) * sin_phi(w)
}

/// Computes the reflection direction given an incident direction and
/// a surface normal.
pub fn reflect(wo: Vector3f, n: Vector3f) -> Vector3f {
    -wo + n * 2.0 as Float * vec3_dot_vec3(wo, n)
}

/// Computes the refraction direction given an incident direction, a
/// surface normal, and the ratio of indices of refraction (incident
/// and transmitted).
pub fn refract(wi: Vector3f, n: Normal3f, eta: Float, wt: &mut Vector3f) -> bool {
    // compute $\cos \theta_\roman{t}$ using Snell's law
    let cos_theta_i: Float = nrm_dot_vec3(n, wi);
    let sin2_theta_i: Float = (0.0 as Float).max(1.0 as Float - cos_theta_i * cos_theta_i);
    let sin2_theta_t: Float = eta * eta * sin2_theta_i;
    // handle total internal reflection for transmission
    if sin2_theta_t >= 1.0 as Float {
        return false;
    }
    let cos_theta_t: Float = (1.0 as Float - sin2_theta_t).sqrt();
    *wt = -wi * eta + Vector3f::from(n) * (eta * cos_theta_i - cos_theta_t);
    true
}

/// Check that two vectors lie on the same side of of the surface.
pub fn vec3_same_hemisphere_vec3(w: Vector3f, wp: Vector3f) -> bool {
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
    let r_parl: Float = ((local_eta_t * *cos_theta_i) - (local_eta_i * cos_theta_t)) /
                        ((local_eta_t * *cos_theta_i) + (local_eta_i * cos_theta_t));
    let r_perp: Float = ((local_eta_i * *cos_theta_i) - (local_eta_t * cos_theta_t)) /
                        ((local_eta_i * *cos_theta_i) + (local_eta_t * cos_theta_t));
    (r_parl * r_parl + r_perp * r_perp) / 2.0
}
