//std
use std;
use std::f32::consts::PI;
use std::sync::Arc;
// pbrt
use crate::core::geometry::{Point2f, Vector3f};
use crate::core::interaction::SurfaceInteraction;
use crate::core::material::{Material, TransportMode};
use crate::core::paramset::TextureParams;
use crate::core::pbrt::{clamp_t, radians};
use crate::core::pbrt::{Float, Spectrum};
use crate::core::reflection::{abs_cos_theta, fr_dielectric};
use crate::core::reflection::{Bsdf, Bxdf, BxdfType};
use crate::core::texture::Texture;
use crate::textures::constant::ConstantTexture;

pub struct HairMaterial {
    pub sigma_a: Option<Arc<dyn Texture<Spectrum> + Sync + Send>>,
    pub color: Option<Arc<dyn Texture<Spectrum> + Sync + Send>>,
    pub eumelanin: Option<Arc<dyn Texture<Float> + Sync + Send>>,
    pub pheomelanin: Option<Arc<dyn Texture<Float> + Sync + Send>>,
    pub eta: Arc<dyn Texture<Float> + Sync + Send>, // default: 1.55
    pub beta_m: Arc<dyn Texture<Float> + Sync + Send>, // default: 0.3
    pub beta_n: Arc<dyn Texture<Float> + Sync + Send>, // default: 0.3
    pub alpha: Arc<dyn Texture<Float> + Sync + Send>, // default: 2.0
}

impl HairMaterial {
    pub fn new(
        sigma_a: Option<Arc<dyn Texture<Spectrum> + Send + Sync>>,
        color: Option<Arc<dyn Texture<Spectrum> + Send + Sync>>,
        eumelanin: Option<Arc<dyn Texture<Float> + Send + Sync>>,
        pheomelanin: Option<Arc<dyn Texture<Float> + Send + Sync>>,
        eta: Arc<dyn Texture<Float> + Send + Sync>,
        beta_m: Arc<dyn Texture<Float> + Send + Sync>,
        beta_n: Arc<dyn Texture<Float> + Send + Sync>,
        alpha: Arc<dyn Texture<Float> + Send + Sync>,
    ) -> Self {
        HairMaterial {
            sigma_a,
            color,
            eumelanin,
            pheomelanin,
            eta,
            beta_m,
            beta_n,
            alpha,
        }
    }
    pub fn create(mp: &mut TextureParams) -> Arc<dyn Material + Send + Sync> {
        let mut sigma_a: Option<Arc<dyn Texture<Spectrum> + Send + Sync>> =
            mp.get_spectrum_texture_or_null("sigma_a");
        let color: Option<Arc<dyn Texture<Spectrum> + Send + Sync>> =
            mp.get_spectrum_texture_or_null("color");
        let eumelanin: Option<Arc<dyn Texture<Float> + Send + Sync>> =
            mp.get_float_texture_or_null("eumelanin");
        let pheomelanin: Option<Arc<dyn Texture<Float> + Send + Sync>> =
            mp.get_float_texture_or_null("pheomelanin");
        if let Some(_sigma_a) = sigma_a.clone() {
            if let Some(_color) = color.clone() {
                println!("WARNING: Ignoring \"color\" parameter since \"sigma_a\" was provided.");
            }
            if let Some(_eumelanin) = eumelanin.clone() {
                println!(
                    "WARNING: Ignoring \"eumelanin\" parameter since \"sigma_a\" was provided."
                );
            }
            if let Some(_pheomelanin) = pheomelanin.clone() {
                println!(
                    "WARNING: Ignoring \"pheomelanin\" parameter since \"sigma_a\" was provided."
                );
            }
        } else if let Some(_color) = color.clone() {
            if let Some(_sigma_a) = sigma_a.clone() {
                println!("WARNING: Ignoring \"sigma_a\" parameter since \"color\" was provided.");
            }
            if let Some(_eumelanin) = eumelanin.clone() {
                println!("WARNING: Ignoring \"eumelanin\" parameter since \"color\" was provided.");
            }
            if let Some(_pheomelanin) = pheomelanin.clone() {
                println!(
                    "WARNING: Ignoring \"pheomelanin\" parameter since \"color\" was provided."
                );
            }
        } else if let Some(_eumelanin) = eumelanin.clone() {
            if let Some(_sigma_a) = sigma_a.clone() {
                println!(
                    "WARNING: Ignoring \"sigma_a\" parameter since \"eumelanin\" was provided."
                );
            }
            if let Some(_color) = color.clone() {
                println!("WARNING: Ignoring \"color\" parameter since \"eumelanin\" was provided.");
            }
        } else if let Some(_pheomelanin) = pheomelanin.clone() {
            if let Some(_sigma_a) = sigma_a.clone() {
                println!(
                    "WARNING: Ignoring \"sigma_a\" parameter since \"pheomelanin\" was provided."
                );
            }
            if let Some(_color) = color.clone() {
                println!(
                    "WARNING: Ignoring \"color\" parameter since \"pheomelanin\" was provided."
                );
            }
        } else {
            // default: brown-ish hair.
            sigma_a = Some(Arc::new(ConstantTexture::new(
                HairBSDF::sigma_a_from_concentration(1.3 as Float, 0.0 as Float),
            )));
        }
        let eta = mp.get_float_texture("eta", 1.55);
        let beta_m = mp.get_float_texture("beta_m", 0.3);
        let beta_n = mp.get_float_texture("beta_n", 0.3);
        let alpha = mp.get_float_texture("alpha", 2.0);
        Arc::new(HairMaterial::new(
            sigma_a,
            color,
            eumelanin,
            pheomelanin,
            eta,
            beta_m,
            beta_n,
            alpha,
        ))
    }
}

impl Material for HairMaterial {
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        // arena: &mut Arena,
        _mode: TransportMode,
        _allow_multiple_lobes: bool,
        _material: Option<Arc<dyn Material + Send + Sync>>,
        scale_opt: Option<Spectrum>,
    ) -> Vec<Bxdf> {
        let mut use_scale: bool = false;
        let mut sc: Spectrum = Spectrum::default();
        if let Some(scale) = scale_opt {
            use_scale = true;
            sc = scale;
        }
        let mut bxdfs: Vec<Bxdf> = Vec::new();
        let bm: Float = self.beta_m.evaluate(si);
        let bn: Float = self.beta_n.evaluate(si);
        let a: Float = radians(self.alpha.evaluate(si));
        let e: Float = self.eta.evaluate(si);
        let sig_a: Spectrum;
        if let Some(ref sigma_a) = self.sigma_a {
            sig_a = sigma_a.evaluate(si);
        } else if let Some(ref color) = self.color {
            let c: Spectrum = color
                .evaluate(si)
                .clamp(0.0 as Float, std::f32::INFINITY as Float);
            sig_a = HairBSDF::sigma_a_from_reflectance(c, bn);
        } else {
            let mut ce: Float = 0.0 as Float;
            let mut cp: Float = 0.0 as Float;
            if let Some(ref eumelanin) = self.eumelanin {
                ce = (0.0 as Float).max(eumelanin.evaluate(si));
                if let Some(ref pheomelanin) = self.pheomelanin {
                    cp = (0.0 as Float).max(pheomelanin.evaluate(si));
                }
            } else {
                if let Some(ref pheomelanin) = self.pheomelanin {
                    cp = (0.0 as Float).max(pheomelanin.evaluate(si));
                }
            }
            sig_a = HairBSDF::sigma_a_from_concentration(ce, cp);
        }
        let h: Float = -1.0 as Float + 2.0 as Float * si.uv[1];
        if use_scale {
            bxdfs.push(Bxdf::Hair(HairBSDF::new(h, e, sig_a, bm, bn, a, Some(sc))));
        } else {
            bxdfs.push(Bxdf::Hair(HairBSDF::new(h, e, sig_a, bm, bn, a, None)));
        }
        si.bsdf = Some(Arc::new(Bsdf::new(si, 1.0, Vec::new())));
        bxdfs
    }
}

pub const P_MAX: u8 = 3_u8;
pub const SQRT_PI_OVER_8: Float = 0.626657069 as Float;

pub struct HairBSDF {
    pub h: Float,
    pub gamma_o: Float,
    pub eta: Float,
    pub sigma_a: Spectrum,
    pub beta_m: Float,
    pub beta_n: Float,
    pub v: [Float; (P_MAX + 1) as usize],
    pub s: Float,
    pub sin_2k_alpha: [Float; 3],
    pub cos_2k_alpha: [Float; 3],
    sc_opt: Option<Spectrum>,
}

impl HairBSDF {
    pub fn new(
        h: Float,
        eta: Float,
        sigma_a: Spectrum,
        beta_m: Float,
        beta_n: Float,
        alpha: Float,
        sc_opt: Option<Spectrum>,
    ) -> Self {
        assert!(h >= -1.0 as Float && h <= 1.0 as Float);
        assert!(beta_m >= 0.0 as Float && beta_m <= 1.0 as Float);
        assert!(beta_n >= 0.0 as Float && beta_n <= 1.0 as Float);
        // gamma_o(SafeASin(h))
        assert!(h >= -1.0001 as Float && h <= 1.0001 as Float);
        let gamma_o: Float = clamp_t(h, -1.0 as Float, 1.0 as Float).asin();
        // compute longitudinal variance from $\beta_m$
        assert!(
            P_MAX >= 3_u8,
            "Longitudinal variance code must be updated to handle low P_MAX"
        );
        let mut v: [Float; (P_MAX + 1) as usize] = [0.0 as Float; (P_MAX + 1) as usize];
        let beta_m_2: Float = beta_m * beta_m;
        let beta_m_4: Float = beta_m_2 * beta_m_2;
        let beta_m_8: Float = beta_m_4 * beta_m_4;
        let beta_m_16: Float = beta_m_8 * beta_m_8;
        let beta_m_20: Float = beta_m_16 * beta_m_4;
        let f: Float =
            0.726 as Float * beta_m + 0.812 as Float * beta_m_2 + 3.7 as Float * beta_m_20;
        v[0] = f * f;
        v[1] = 0.25 as Float * v[0];
        v[2] = 4.0 as Float * v[0];
        for p in 3..(P_MAX + 1) {
            // TODO: is there anything better here?
            v[p as usize] = v[2];
        }
        // compute azimuthal logistic scale factor from $\beta_n$
        let beta_n_2: Float = beta_n * beta_n;
        let beta_n_4: Float = beta_n_2 * beta_n_2;
        let beta_n_8: Float = beta_n_4 * beta_n_4;
        let beta_n_16: Float = beta_n_8 * beta_n_8;
        let beta_n_22: Float = beta_n_16 * beta_n_4 * beta_n_2;
        let s: Float = SQRT_PI_OVER_8
            * (0.265 as Float * beta_n + 1.194 as Float * beta_n_2 + 5.372 as Float * beta_n_22);
        assert!(!s.is_nan());
        // compute $\alpha$ terms for hair scales
        let mut sin_2k_alpha: [Float; 3] = [0.0 as Float; 3];
        sin_2k_alpha[0] = radians(alpha).sin();
        // cos_2k_alpha[0] = SafeSqrt(1 - Sqr(sin_2k_alpha[0]));
        let mut cos_2k_alpha: [Float; 3] = [0.0 as Float; 3];
        let sqr: Float = sin_2k_alpha[0] * sin_2k_alpha[0];
        let one_minus_sqr: Float = 1.0 as Float - sqr;
        assert!(one_minus_sqr >= -1e-4);
        cos_2k_alpha[0] = (0.0 as Float).max(one_minus_sqr).sqrt();
        for i in 1..3 {
            sin_2k_alpha[i] = 2.0 as Float * cos_2k_alpha[i - 1] * sin_2k_alpha[i - 1];
            cos_2k_alpha[i] = (cos_2k_alpha[i - 1] * cos_2k_alpha[i - 1])
                - (sin_2k_alpha[i - 1] * sin_2k_alpha[i - 1]);
        }
        HairBSDF {
            h,
            gamma_o,
            eta,
            sigma_a,
            beta_m,
            beta_n,
            v,
            s,
            sin_2k_alpha,
            cos_2k_alpha,
            sc_opt,
        }
    }
    pub fn compute_ap_pdf(&self, cos_theta_o: Float) -> [Float; (P_MAX + 1) as usize] {
        // compute array of $A_p$ values for _cos_theta_o_
        let x: Float = 1.0 as Float - (cos_theta_o * cos_theta_o);
        assert!(x >= -1e-4);
        let sin_theta_o: Float = (0.0 as Float).max(x).sqrt();
        // compute $\cos \thetat$ for refracted ray
        let sin_theta_t: Float = sin_theta_o / self.eta;
        let x: Float = 1.0 as Float - (sin_theta_t * sin_theta_t);
        assert!(x >= -1e-4);
        let cos_theta_t: Float = (0.0 as Float).max(x).sqrt();
        // compute $\gammat$ for refracted ray
        let etap: Float = (self.eta * self.eta - (sin_theta_o * sin_theta_o)).sqrt() / cos_theta_o;
        let sin_gamma_t: Float = self.h / etap;
        let x: Float = 1.0 as Float - (sin_gamma_t * sin_gamma_t);
        assert!(x >= -1e-4);
        let cos_gamma_t: Float = (0.0 as Float).max(x).sqrt();
        // compute the transmittance _t_ of a single path through the cylinder
        let t: Spectrum = (-self.sigma_a * (2.0 as Float * cos_gamma_t / cos_theta_t)).exp();
        let ap: [Spectrum; (P_MAX + 1) as usize] = ap(cos_theta_o, self.eta, self.h, t);
        // compute $A_p$ PDF from individual $A_p$ terms
        let mut ap_pdf: [Float; (P_MAX + 1) as usize] = [0.0 as Float; (P_MAX + 1) as usize];
        let mut sum_y: Float = 0.0 as Float;
        for i in 0..(P_MAX + 1) {
            sum_y += ap[i as usize].y();
        }
        for i in 0..(P_MAX + 1) {
            ap_pdf[i as usize] = ap[i as usize].y() / sum_y;
        }
        ap_pdf
    }
    pub fn sigma_a_from_concentration(ce: Float, cp: Float) -> Spectrum {
        let mut sigma_a: [Float; 3] = [0.0 as Float; 3];
        let eumelanin_sigma_a: [Float; 3] = [0.419 as Float, 0.697 as Float, 1.37 as Float];
        let pheomelanin_sigma_a: [Float; 3] = [0.187 as Float, 0.4 as Float, 1.05 as Float];
        for i in 0..3 {
            sigma_a[i] = ce * eumelanin_sigma_a[i] + cp * pheomelanin_sigma_a[i];
        }
        Spectrum::from_rgb(&sigma_a)
    }
    pub fn sigma_a_from_reflectance(c: Spectrum, beta_n: Float) -> Spectrum {
        let mut sigma_a: Spectrum = Spectrum::default();
        for i in 0..3 {
            let sqr: Float = beta_n * beta_n;
            let pow3: Float = sqr * beta_n;
            let pow4: Float = pow3 * beta_n;
            let pow5: Float = pow4 * beta_n;
            let f: Float = c.c[i].ln()
                / (5.969 as Float - 0.215 as Float * beta_n + 2.532 as Float * sqr
                    - 10.73 as Float * pow3
                    + 5.574 as Float * pow4
                    + 0.245 as Float * pow5);
            sigma_a.c[i] = f * f;
        }
        sigma_a
    }
    // Bxdf
    pub fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        // compute hair coordinate system terms related to _wo_
        let sin_theta_o: Float = wo.x;
        // Float cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
        let x: Float = 1.0 as Float - (sin_theta_o * sin_theta_o);
        assert!(x >= -1e-4);
        let cos_theta_o: Float = (0.0 as Float).max(x).sqrt();
        let phi_o: Float = wo.z.atan2(wo.y);
        // compute hair coordinate system terms related to _wi_
        let sin_theta_i: Float = wi.x;
        // Float cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));
        let x: Float = 1.0 as Float - (sin_theta_i * sin_theta_i);
        assert!(x >= -1e-4);
        let cos_theta_i: Float = (0.0 as Float).max(x).sqrt();
        let phi_i: Float = wi.z.atan2(wi.y);
        // compute $\cos \thetat$ for refracted ray
        let sin_theta_t: Float = sin_theta_o / self.eta;
        // Float cosThetaT = SafeSqrt(1 - Sqr(sinThetaT));
        let x: Float = 1.0 as Float - (sin_theta_t * sin_theta_t);
        assert!(x >= -1e-4);
        let cos_theta_t: Float = (0.0 as Float).max(x).sqrt();
        // compute $\gammat$ for refracted ray
        let etap: Float = (self.eta * self.eta - (sin_theta_o * sin_theta_o)).sqrt() / cos_theta_o;
        let sin_gamma_t: Float = self.h / etap;
        // Float cosGammaT = SafeSqrt(1 - Sqr(sinGammaT));
        let x: Float = 1.0 as Float - (sin_gamma_t * sin_gamma_t);
        assert!(x >= -1e-4);
        let cos_gamma_t: Float = (0.0 as Float).max(x).sqrt();
        // Float gammaT = SafeASin(sinGammaT);
        let x: Float = sin_gamma_t;
        assert!(x >= -1.0001 && x <= 1.0001);
        let gamma_t: Float = clamp_t(x, -1.0 as Float, 1.0 as Float).asin();
        // compute the transmittance _t_ of a single path through the cylinder
        let t: Spectrum = (-self.sigma_a * (2.0 as Float * cos_gamma_t / cos_theta_t)).exp();
        // evaluate hair BSDF
        let phi: Float = phi_i - phi_o;
        let ap: [Spectrum; (P_MAX + 1) as usize] = ap(cos_theta_o, self.eta, self.h, t);
        let mut fsum: Spectrum = Spectrum::default();
        for p in 0..P_MAX {
            // compute $\sin \thetai$ and $\cos \thetai$ terms accounting for scales
            let sin_theta_ip: Float;
            let mut cos_theta_ip: Float;
            if p == 0 {
                sin_theta_ip =
                    sin_theta_i * self.cos_2k_alpha[1] + cos_theta_i * self.sin_2k_alpha[1];
                cos_theta_ip =
                    cos_theta_i * self.cos_2k_alpha[1] - sin_theta_i * self.sin_2k_alpha[1];
            } else if p == 1 {
                sin_theta_ip =
                    sin_theta_i * self.cos_2k_alpha[0] - cos_theta_i * self.sin_2k_alpha[0];
                cos_theta_ip =
                    cos_theta_i * self.cos_2k_alpha[0] + sin_theta_i * self.sin_2k_alpha[0];
            } else if p == 2 {
                sin_theta_ip =
                    sin_theta_i * self.cos_2k_alpha[2] - cos_theta_i * self.sin_2k_alpha[2];
                cos_theta_ip =
                    cos_theta_i * self.cos_2k_alpha[2] + sin_theta_i * self.sin_2k_alpha[2];
            } else {
                sin_theta_ip = sin_theta_i;
                cos_theta_ip = cos_theta_i;
            }
            // handle out-of-range $\cos \thetai$ from scale adjustment
            cos_theta_ip = cos_theta_ip.abs();
            fsum += ap[p as usize]
                * mp(
                    cos_theta_ip,
                    cos_theta_o,
                    sin_theta_ip,
                    sin_theta_o,
                    self.v[p as usize],
                )
                * np(phi, p as i32, self.s, self.gamma_o, gamma_t);
        }
        // compute contribution of remaining terms after _pMax_
        fsum += ap[P_MAX as usize]
            * mp(
                cos_theta_i,
                cos_theta_o,
                sin_theta_i,
                sin_theta_o,
                self.v[P_MAX as usize],
            )
            / (2.0 as Float * PI);
        if abs_cos_theta(wi) > 0.0 as Float {
            fsum = fsum / abs_cos_theta(wi);
        }
        assert!(!fsum.y().is_infinite() && !fsum.y().is_nan());
        if let Some(sc) = self.sc_opt {
            sc * fsum
        } else {
            fsum
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
        // compute hair coordinate system terms related to _wo_
        let sin_theta_o: Float = wo.x;
        let x: Float = 1.0 as Float - (sin_theta_o * sin_theta_o);
        assert!(x >= -1e-4);
        let cos_theta_o: Float = (0.0 as Float).max(x).sqrt();
        let phi_o: Float = wo.z.atan2(wo.y);
        // derive four random samples from _sample_
        let mut u: [Point2f; 2] = [demux_float(sample[0]), demux_float(sample[1])];
        // determine which term $p$ to sample for hair scattering
        let ap_pdf: [Float; (P_MAX + 1) as usize] = self.compute_ap_pdf(cos_theta_o);
        let mut p: usize = 0;
        for i in 0..P_MAX {
            p = i as usize; // store index in p for later
            if u[0][0] < ap_pdf[p] {
                break;
            }
            u[0][0] -= ap_pdf[p];
        }
        // sample $M_p$ to compute $\thetai$
        u[1][0] = u[1][0].max(1e-5 as Float);
        let cos_theta: Float = 1.0 as Float
            + self.v[p]
                * (u[1][0] + (1.0 as Float - u[1][0]) * (-2.0 as Float / self.v[p]).exp()).ln();
        let x: Float = 1.0 as Float - (cos_theta * cos_theta);
        assert!(x >= -1e-4);
        let sin_theta: Float = (0.0 as Float).max(x).sqrt();
        let cos_phi: Float = (2.0 as Float * PI * u[1][1]).cos();
        let mut sin_theta_i: Float = -cos_theta * sin_theta_o + sin_theta * cos_phi * cos_theta_o;
        let x: Float = 1.0 as Float - (sin_theta_i * sin_theta_i);
        assert!(x >= -1e-4);
        let mut cos_theta_i: Float = (0.0 as Float).max(x).sqrt();
        // update sampled $\sin \thetai$ and $\cos \thetai$ to account for scales
        let mut sin_theta_ip: Float = sin_theta_i;
        let mut cos_theta_ip: Float = cos_theta_i;
        if p == 0_usize {
            sin_theta_ip = sin_theta_i * self.cos_2k_alpha[1] - cos_theta_i * self.sin_2k_alpha[1];
            cos_theta_ip = cos_theta_i * self.cos_2k_alpha[1] + sin_theta_i * self.sin_2k_alpha[1];
        } else if p == 1_usize {
            sin_theta_ip = sin_theta_i * self.cos_2k_alpha[0] + cos_theta_i * self.sin_2k_alpha[0];
            cos_theta_ip = cos_theta_i * self.cos_2k_alpha[0] - sin_theta_i * self.sin_2k_alpha[0];
        } else if p == 2_usize {
            sin_theta_ip = sin_theta_i * self.cos_2k_alpha[2] + cos_theta_i * self.sin_2k_alpha[2];
            cos_theta_ip = cos_theta_i * self.cos_2k_alpha[2] - sin_theta_i * self.sin_2k_alpha[2];
        }
        sin_theta_i = sin_theta_ip;
        cos_theta_i = cos_theta_ip;

        // sample $N_p$ to compute $\Delta\phi$

        // compute $\gammat$ for refracted ray
        let etap: Float = (self.eta * self.eta - sin_theta_o * sin_theta_o).sqrt() / cos_theta_o;
        let sin_gamma_t: Float = self.h / etap;
        assert!(sin_gamma_t >= -1.0001 as Float && sin_gamma_t <= 1.0001 as Float);
        let gamma_t: Float = clamp_t(sin_gamma_t, -1.0 as Float, 1.0 as Float).asin();
        let dphi: Float;
        if p < P_MAX as usize {
            dphi = phi_fn(p as i32, self.gamma_o, gamma_t)
                + sample_trimmed_logistic(u[0][1], self.s, -PI, PI);
        } else {
            dphi = 2.0 as Float * PI * u[0][1];
        }
        // compute _wi_ from sampled hair scattering angles
        let phi_i: Float = phi_o + dphi;
        *wi = Vector3f {
            x: sin_theta_i,
            y: cos_theta_i * phi_i.cos(),
            z: cos_theta_i * phi_i.sin(),
        };

        // compute PDF for sampled hair scattering direction _wi_
        *pdf = 0.0 as Float;
        for p in 0..P_MAX {
            // compute $\sin \thetai$ and $\cos \thetai$ terms accounting for scales
            //     Float sin_theta_ip, cos_theta_ip;
            if p == 0 {
                sin_theta_ip =
                    sin_theta_i * self.cos_2k_alpha[1] + cos_theta_i * self.sin_2k_alpha[1];
                cos_theta_ip =
                    cos_theta_i * self.cos_2k_alpha[1] - sin_theta_i * self.sin_2k_alpha[1];
            } else if p == 1 {
                sin_theta_ip =
                    sin_theta_i * self.cos_2k_alpha[0] - cos_theta_i * self.sin_2k_alpha[0];
                cos_theta_ip =
                    cos_theta_i * self.cos_2k_alpha[0] + sin_theta_i * self.sin_2k_alpha[0];
            } else if p == 2 {
                sin_theta_ip =
                    sin_theta_i * self.cos_2k_alpha[2] - cos_theta_i * self.sin_2k_alpha[2];
                cos_theta_ip =
                    cos_theta_i * self.cos_2k_alpha[2] + sin_theta_i * self.sin_2k_alpha[2];
            } else {
                sin_theta_ip = sin_theta_i;
                cos_theta_ip = cos_theta_i;
            }

            // handle out-of-range $\cos \thetai$ from scale adjustment
            cos_theta_ip = cos_theta_ip.abs();
            *pdf += ap_pdf[p as usize]
                * mp(
                    cos_theta_ip,
                    cos_theta_o,
                    sin_theta_ip,
                    sin_theta_o,
                    self.v[p as usize],
                )
                * np(dphi, p as i32, self.s, self.gamma_o, gamma_t);
        }
        *pdf += ap_pdf[P_MAX as usize]
            * mp(
                cos_theta_i,
                cos_theta_o,
                sin_theta_i,
                sin_theta_o,
                self.v[P_MAX as usize],
            )
            * (1.0 as Float / (2.0 as Float * PI));
        if let Some(sc) = self.sc_opt {
            sc * self.f(wo, &*wi)
        } else {
            self.f(wo, &*wi)
        }
    }
    pub fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        // compute hair coordinate system terms related to _wo_
        let sin_theta_o: Float = wo.x;
        let x: Float = 1.0 as Float - (sin_theta_o * sin_theta_o);
        assert!(x >= -1e-4);
        let cos_theta_o: Float = (0.0 as Float).max(x).sqrt();
        let phi_o: Float = wo.z.atan2(wo.y);
        // compute hair coordinate system terms related to _wi_
        let sin_theta_i: Float = wi.x;
        let x: Float = 1.0 as Float - (sin_theta_i * sin_theta_i);
        assert!(x >= -1e-4);
        let cos_theta_i: Float = (0.0 as Float).max(x).sqrt();
        let phi_i: Float = wi.z.atan2(wi.y);
        // compute $\gammat$ for refracted ray
        let etap: Float = (self.eta * self.eta - (sin_theta_o * sin_theta_o)).sqrt() / cos_theta_o;
        let sin_gamma_t: Float = self.h / etap;
        let x: Float = sin_gamma_t;
        assert!(x >= -1.0001 && x <= 1.0001);
        let gamma_t: Float = clamp_t(x, -1.0 as Float, 1.0 as Float).asin();
        // compute PDF for $A_p$ terms
        let ap_pdf: [Float; (P_MAX + 1) as usize] = self.compute_ap_pdf(cos_theta_o);
        // compute PDF sum for hair scattering events
        let phi: Float = phi_i - phi_o;
        let mut pdf: Float = 0.0 as Float;
        for p in 0..P_MAX {
            // compute $\sin \thetai$ and $\cos \thetai$ terms accounting for scales
            let sin_theta_ip: Float;
            let mut cos_theta_ip: Float;
            if p == 0 {
                sin_theta_ip =
                    sin_theta_i * self.cos_2k_alpha[1] + cos_theta_i * self.sin_2k_alpha[1];
                cos_theta_ip =
                    cos_theta_i * self.cos_2k_alpha[1] - sin_theta_i * self.sin_2k_alpha[1];
            } else if p == 1 {
                sin_theta_ip =
                    sin_theta_i * self.cos_2k_alpha[0] - cos_theta_i * self.sin_2k_alpha[0];
                cos_theta_ip =
                    cos_theta_i * self.cos_2k_alpha[0] + sin_theta_i * self.sin_2k_alpha[0];
            } else if p == 2 {
                sin_theta_ip =
                    sin_theta_i * self.cos_2k_alpha[2] - cos_theta_i * self.sin_2k_alpha[2];
                cos_theta_ip =
                    cos_theta_i * self.cos_2k_alpha[2] + sin_theta_i * self.sin_2k_alpha[2];
            } else {
                sin_theta_ip = sin_theta_i;
                cos_theta_ip = cos_theta_i;
            }
            // handle out-of-range $\cos \thetai$ from scale adjustment
            cos_theta_ip = cos_theta_ip.abs();
            pdf += ap_pdf[p as usize]
                * mp(
                    cos_theta_ip,
                    cos_theta_o,
                    sin_theta_ip,
                    sin_theta_o,
                    self.v[p as usize],
                )
                * np(phi, p as i32, self.s, self.gamma_o, gamma_t);
        }
        pdf += ap_pdf[P_MAX as usize]
            * mp(
                cos_theta_i,
                cos_theta_o,
                sin_theta_i,
                sin_theta_o,
                self.v[P_MAX as usize],
            )
            * (1.0 as Float / (2.0 as Float * PI));
        pdf
    }
    pub fn get_type(&self) -> u8 {
        BxdfType::BsdfGlossy as u8
            | BxdfType::BsdfReflection as u8
            | BxdfType::BsdfTransmission as u8
    }
}

// https://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/
fn compact_1_by_1(mut x: u32) -> u32 {
    // TODO: as of Haswell, the PEXT instruction could do all this in a
    // single instruction.
    // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    x &= 0x55555555;
    // x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x ^ (x >> 1)) & 0x33333333;
    // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x ^ (x >> 2)) & 0x0f0f0f0f;
    // x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x ^ (x >> 4)) & 0x00ff00ff;
    // x = ---- ---- ---- ---- fedc ba98 7654 3210
    x = (x ^ (x >> 8)) & 0x0000ffff;
    x
}

fn demux_float(f: Float) -> Point2f {
    assert!(f >= 0.0 as Float && f < 1.0 as Float);
    let v: u64 = (f * (1u64 << 32) as Float) as u64;
    assert!(v < 0x100000000);
    let bits: [u32; 2] = [compact_1_by_1(v as u32), compact_1_by_1((v >> 1) as u32)];
    Point2f {
        x: bits[0] as Float / (1 << 16) as Float,
        y: bits[1] as Float / (1 << 16) as Float,
    }
}

// Hair Local Functions

fn mp(
    cos_theta_i: Float,
    cos_theta_o: Float,
    sin_theta_i: Float,
    sin_theta_o: Float,
    v: Float,
) -> Float {
    let a: Float = cos_theta_i * cos_theta_o / v;
    let b: Float = sin_theta_i * sin_theta_o / v;
    let mp: Float;
    if v <= 0.1 as Float {
        mp = (log_i0(a) - b - 1.0 as Float / v
            + 0.6931 as Float
            + (1.0 as Float / (2. as Float * v)).ln())
        .exp();
    } else {
        mp = ((-b).exp() * i0(a)) / ((1.0 as Float / v).sinh() * 2.0 as Float * v);
    }
    mp
}

#[inline]
fn i0(x: Float) -> Float {
    let mut val: Float = 0.0 as Float;
    let mut x2i: Float = 1.0 as Float;
    let mut ifact: i32 = 1_i32;
    let mut i4: i32 = 1_i32;
    // i0(x) \approx Sum_i x^(2i) / (4^i (i!)^2)
    for i in 0..10 {
        if i as i32 > 1_i32 {
            ifact *= i as i32;
        }
        val += x2i / (i4 as Float * (ifact as Float * ifact as Float));
        x2i *= x * x;
        i4 *= 4_i32;
    }
    val
}

#[inline]
fn log_i0(x: Float) -> Float {
    if x > 12.0 as Float {
        x + 0.5
            * (-((2.0 as Float * PI).ln())
                + (1.0 as Float / x).ln()
                + 1.0 as Float / (8.0 as Float * x))
    } else {
        i0(x).ln()
    }
}

fn ap(cos_theta_o: Float, eta: Float, h: Float, t: Spectrum) -> [Spectrum; (P_MAX + 1) as usize] {
    let mut ap: [Spectrum; (P_MAX + 1) as usize] = [Spectrum::default(); (P_MAX + 1) as usize];
    // compute $p=0$ attenuation at initial cylinder intersection
    // Float cosGammaO = SafeSqrt(1 - h * h);
    let x: Float = 1.0 as Float - (h * h);
    assert!(x >= -1e-4);
    let cos_gamma_o: Float = (0.0 as Float).max(x).sqrt();
    let cos_theta: Float = cos_theta_o * cos_gamma_o;
    let f: Float = fr_dielectric(cos_theta, 1.0 as Float, eta);
    ap[0] = Spectrum::new(f);
    // compute $p=1$ attenuation term
    let one_minus_f: Float = 1.0 as Float - f;
    ap[1] = t * (one_minus_f * one_minus_f);
    // compute attenuation terms up to $p=_P_MAX_$
    for p in 2..P_MAX {
        // TODO: is there anything better here?
        ap[p as usize] = ap[(p - 1) as usize] * t * f;
    }
    // compute attenuation term accounting for remaining orders of scattering
    ap[P_MAX as usize] = ap[(P_MAX - 1) as usize] * t * f / (Spectrum::new(1.0 as Float) - t * f);
    ap
}

#[inline]
fn phi_fn(p: i32, gamma_o: Float, gamma_t: Float) -> Float {
    2.0 as Float * p as Float * gamma_t - 2.0 as Float * gamma_o + p as Float * PI
}

#[inline]
fn logistic(x: Float, s: Float) -> Float {
    let x: Float = x.abs();
    let e: Float = (-x / s).exp();
    let one_plus_e: Float = 1.0 as Float + e;
    e / (s * (one_plus_e * one_plus_e))
}

#[inline]
fn logistic_cdf(x: Float, s: Float) -> Float {
    let e: Float = (-x / s).exp();
    let one_plus_e: Float = 1.0 as Float + e;
    1.0 as Float / one_plus_e
}

#[inline]
fn trimmed_logistic(x: Float, s: Float, a: Float, b: Float) -> Float {
    assert!(a < b);
    logistic(x, s) / (logistic_cdf(b, s) - logistic_cdf(a, s))
}

#[inline]
fn np(phi: Float, p: i32, s: Float, gamma_o: Float, gamma_t: Float) -> Float {
    let mut dphi: Float = phi - phi_fn(p, gamma_o, gamma_t);
    // remap _dphi_ to $[-\pi,\pi]$
    while dphi > PI {
        dphi -= 2.0 as Float * PI;
    }
    while dphi < -PI {
        dphi += 2.0 as Float * PI;
    }
    return trimmed_logistic(dphi, s, -PI, PI);
}

fn sample_trimmed_logistic(u: Float, s: Float, a: Float, b: Float) -> Float {
    assert!(a < b);
    let k: Float = logistic_cdf(b, s) - logistic_cdf(a, s);
    let x: Float = -s * (1.0 as Float / (u * k + logistic_cdf(a, s)) - 1.0 as Float).ln();
    assert!(!x.is_nan());
    clamp_t(x, a, b)
}
