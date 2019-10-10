//std
use std::sync::Arc;
// pbrt
use crate::core::interaction::SurfaceInteraction;
use crate::core::material::{Material, TransportMode};
use crate::core::microfacet::TrowbridgeReitzDistribution;
use crate::core::paramset::TextureParams;
use crate::core::pbrt::{Float, Spectrum};
use crate::core::reflection::{
    Bsdf, Bxdf, Fresnel, FresnelDielectric, LambertianReflection, LambertianTransmission,
    MicrofacetReflection, MicrofacetTransmission,
};
use crate::core::texture::Texture;

pub struct TranslucentMaterial {
    pub kd: Arc<dyn Texture<Spectrum> + Sync + Send>, // default: 0.25
    pub ks: Arc<dyn Texture<Spectrum> + Sync + Send>, // default: 0.25
    pub roughness: Arc<dyn Texture<Float> + Sync + Send>, // default: 0.1
    pub reflect: Arc<dyn Texture<Spectrum> + Sync + Send>, // default: 0.5
    pub transmit: Arc<dyn Texture<Spectrum> + Sync + Send>, // default: 0.5
    pub bump_map: Option<Arc<dyn Texture<Float> + Send + Sync>>,
    pub remap_roughness: bool, // default: true
}

impl TranslucentMaterial {
    pub fn new(
        kd: Arc<dyn Texture<Spectrum> + Send + Sync>,
        ks: Arc<dyn Texture<Spectrum> + Send + Sync>,
        roughness: Arc<dyn Texture<Float> + Sync + Send>,
        reflect: Arc<dyn Texture<Spectrum> + Send + Sync>,
        transmit: Arc<dyn Texture<Spectrum> + Send + Sync>,
        bump_map: Option<Arc<dyn Texture<Float> + Sync + Send>>,
        remap_roughness: bool,
    ) -> Self {
        TranslucentMaterial {
            kd,
            ks,
            roughness,
            reflect,
            transmit,
            bump_map,
            remap_roughness,
        }
    }
    pub fn create(mp: &mut TextureParams) -> Arc<dyn Material + Send + Sync> {
        let kd = mp.get_spectrum_texture("Kd", Spectrum::new(0.25 as Float));
        let ks = mp.get_spectrum_texture("Ks", Spectrum::new(0.25 as Float));
        let reflect = mp.get_spectrum_texture("reflect", Spectrum::new(0.5 as Float));
        let transmit = mp.get_spectrum_texture("transmit", Spectrum::new(0.5 as Float));
        let roughness = mp.get_float_texture("roughness", 0.1 as Float);
        let bump_map = mp.get_float_texture_or_null("bumpmap");
        let remap_roughness: bool = mp.find_bool("remaproughness", true);
        Arc::new(TranslucentMaterial::new(
            kd,
            ks,
            roughness,
            reflect,
            transmit,
            bump_map,
            remap_roughness,
        ))
    }
}

impl Material for TranslucentMaterial {
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        // arena: &mut Arena,
        mode: TransportMode,
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
        if let Some(ref bump) = self.bump_map {
            Self::bump(bump, si);
        }
        let mut bxdfs: Vec<Bxdf> = Vec::new();
        let eta: Float = 1.5;
        let r: Spectrum = self
            .reflect
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        let t: Spectrum = self
            .transmit
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        if r.is_black() && t.is_black() {
            si.bsdf = Some(Arc::new(Bsdf::new(si, eta, Vec::new())));
            return bxdfs;
        }
        let kd: Spectrum = self
            .kd
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        if !kd.is_black() {
            if !r.is_black() {
                if use_scale {
                    bxdfs.push(Bxdf::LambertianRefl(LambertianReflection::new(
                        r * kd,
                        Some(sc),
                    )));
                } else {
                    bxdfs.push(Bxdf::LambertianRefl(LambertianReflection::new(
                        r * kd,
                        None,
                    )));
                }
            }
            if !t.is_black() {
                if use_scale {
                    bxdfs.push(Bxdf::LambertianTrans(LambertianTransmission::new(
                        t * kd,
                        Some(sc),
                    )));
                } else {
                    bxdfs.push(Bxdf::LambertianTrans(LambertianTransmission::new(
                        t * kd,
                        None,
                    )));
                }
            }
        }
        let ks: Spectrum = self
            .ks
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        if !ks.is_black() && (!r.is_black() || !t.is_black()) {
            let mut rough: Float = self.roughness.evaluate(si);
            if self.remap_roughness {
                rough = TrowbridgeReitzDistribution::roughness_to_alpha(rough);
            }
            let distrib = Arc::new(TrowbridgeReitzDistribution::new(rough, rough, true));
            if !r.is_black() {
                let fresnel = Fresnel::Dielectric(FresnelDielectric {
                    eta_i: 1.0 as Float,
                    eta_t: eta,
                });
                if use_scale {
                    bxdfs.push(Bxdf::MicrofacetRefl(MicrofacetReflection::new(
                        r * ks,
                        distrib.clone(),
                        fresnel,
                        Some(sc),
                    )));
                } else {
                    bxdfs.push(Bxdf::MicrofacetRefl(MicrofacetReflection::new(
                        r * ks,
                        distrib.clone(),
                        fresnel,
                        None,
                    )));
                }
            }
            if !t.is_black() {
                if use_scale {
                    bxdfs.push(Bxdf::MicrofacetTrans(MicrofacetTransmission::new(
                        t * ks,
                        distrib.clone(),
                        1.0,
                        eta,
                        mode,
                        Some(sc),
                    )));
                } else {
                    bxdfs.push(Bxdf::MicrofacetTrans(MicrofacetTransmission::new(
                        t * ks,
                        distrib.clone(),
                        1.0,
                        eta,
                        mode,
                        None,
                    )));
                }
            }
        }
        si.bsdf = Some(Arc::new(Bsdf::new(si, eta, Vec::new())));
        bxdfs
    }
}
