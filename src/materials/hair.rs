//std
use std;
use std::sync::Arc;
// pbrt
use core::geometry::{Point2f, Vector3f};
use core::interaction::SurfaceInteraction;
use core::material::{Material, TransportMode};
use core::microfacet::TrowbridgeReitzDistribution;
use core::paramset::TextureParams;
use core::pbrt::{Float, Spectrum};
use core::pbrt::{clamp_t, radians};
use core::reflection::{Bsdf, Bxdf, BxdfType, FresnelConductor, MicrofacetReflection};
use core::texture::Texture;
use textures::constant::ConstantTexture;

pub struct HairMaterial {
    pub sigma_a: Option<Arc<Texture<Spectrum> + Sync + Send>>,
    pub color: Option<Arc<Texture<Spectrum> + Sync + Send>>,
    pub eumelanin: Option<Arc<Texture<Float> + Sync + Send>>,
    pub pheomelanin: Option<Arc<Texture<Float> + Sync + Send>>,
    pub eta: Arc<Texture<Float> + Sync + Send>, // default: 1.55
    pub beta_m: Arc<Texture<Float> + Sync + Send>, // default: 0.3
    pub beta_n: Arc<Texture<Float> + Sync + Send>, // default: 0.3
    pub alpha: Arc<Texture<Float> + Sync + Send>, // default: 2.0
}

impl HairMaterial {
    pub fn new(
        sigma_a: Option<Arc<Texture<Spectrum> + Send + Sync>>,
        color: Option<Arc<Texture<Spectrum> + Send + Sync>>,
        eumelanin: Option<Arc<Texture<Float> + Send + Sync>>,
        pheomelanin: Option<Arc<Texture<Float> + Send + Sync>>,
        eta: Arc<Texture<Float> + Send + Sync>,
        beta_m: Arc<Texture<Float> + Send + Sync>,
        beta_n: Arc<Texture<Float> + Send + Sync>,
        alpha: Arc<Texture<Float> + Send + Sync>,
    ) -> Self {
        HairMaterial {
            sigma_a: sigma_a,
            color: color,
            eumelanin: eumelanin,
            pheomelanin: pheomelanin,
            eta: eta,
            beta_m: beta_m,
            beta_n: beta_n,
            alpha: alpha,
        }
    }
    pub fn create(mp: &mut TextureParams) -> Arc<Material + Send + Sync> {
        let mut sigma_a: Option<Arc<Texture<Spectrum> + Send + Sync>> =
            mp.get_spectrum_texture_or_null(String::from("sigma_a"));
        let color: Option<Arc<Texture<Spectrum> + Send + Sync>> =
            mp.get_spectrum_texture_or_null(String::from("color"));
        let eumelanin: Option<Arc<Texture<Float> + Send + Sync>> =
            mp.get_float_texture_or_null(String::from("eumelanin"));
        let pheomelanin: Option<Arc<Texture<Float> + Send + Sync>> =
            mp.get_float_texture_or_null(String::from("pheomelanin"));
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
        let eta = mp.get_float_texture(String::from("eta"), 1.55);
        let beta_m = mp.get_float_texture(String::from("beta_m"), 0.3);
        let beta_n = mp.get_float_texture(String::from("beta_n"), 0.3);
        let alpha = mp.get_float_texture(String::from("alpha"), 2.0);
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
    pub fn bsdf(&self, si: &SurfaceInteraction) -> Bsdf {
        let mut bxdfs: Vec<Arc<Bxdf + Send + Sync>> = Vec::new();
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
        // TODO: bxdfs.push(Arc::new(HairBSDF::new(h, e, sig_a, bm, bn, a)));
        Bsdf::new(si, 1.0, bxdfs)
    }
}

impl Material for HairMaterial {
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        // arena: &mut Arena,
        _mode: TransportMode,
        _allow_multiple_lobes: bool,
    ) {
        si.bsdf = Some(Arc::new(self.bsdf(si)));
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
}

impl HairBSDF {
    pub fn new(
        h: Float,
        eta: Float,
        sigma_a: Spectrum,
        beta_m: Float,
        beta_n: Float,
        alpha: Float,
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
        for p in 3..P_MAX {
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
            h: h,
            gamma_o: gamma_o,
            eta: eta,
            sigma_a: sigma_a,
            beta_m: beta_m,
            beta_n: beta_n,
            v: v,
            s: s,
            sin_2k_alpha: sin_2k_alpha,
            cos_2k_alpha: cos_2k_alpha,
        }
    }
    pub fn sigma_a_from_concentration(ce: Float, cp: Float) -> Spectrum {
        let mut sigma_a: [Float; 3] = [0.0 as Float; 3];
        let eumelanin_sigma_a: [Float; 3] = [0.419 as Float, 0.697 as Float, 1.37 as Float];
        let pheomelanin_sigma_a: [Float; 3] = [0.187 as Float, 0.4 as Float, 1.05 as Float];
        for i in 0..3 {
            sigma_a[i] = (ce * eumelanin_sigma_a[i] + cp * pheomelanin_sigma_a[i]);
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
                    - 10.73 as Float * pow3 + 5.574 as Float * pow4
                    + 0.245 as Float * pow5);
            sigma_a.c[i] = f * f;
        }
        sigma_a
    }
}

impl Bxdf for HairBSDF {
    fn f(&self, _wo: Vector3f, _wi: Vector3f) -> Spectrum {
        // TODO
        Spectrum::new(0.0 as Float)
    }
    fn sample_f(
        &self,
        wo: Vector3f,
        wi: &mut Vector3f,
        _sample: Point2f,
        pdf: &mut Float,
        _sampled_type: &mut u8,
    ) -> Spectrum {
        // TODO
        Spectrum::new(0.0 as Float)
    }
    fn pdf(&self, wo: Vector3f, wi: Vector3f) -> Float {
        // TODO
        0.0 as Float
    }
    fn get_type(&self) -> u8 {
        BxdfType::BsdfGlossy as u8 | BxdfType::BsdfReflection as u8
            | BxdfType::BsdfTransmission as u8
    }
}
