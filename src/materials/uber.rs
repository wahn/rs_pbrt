//std
use std;
use std::sync::Arc;
// pbrt
use crate::core::interaction::SurfaceInteraction;
use crate::core::material::{Material, TransportMode};
use crate::core::microfacet::TrowbridgeReitzDistribution;
use crate::core::paramset::TextureParams;
use crate::core::pbrt::{Float, Spectrum};
use crate::core::reflection::{
    Bsdf, Bxdf, FresnelDielectric, LambertianReflection, MicrofacetReflection, SpecularReflection,
    SpecularTransmission,
};
use crate::core::texture::Texture;

// see uber.h

pub struct UberMaterial {
    pub kd: Arc<dyn Texture<Spectrum> + Sync + Send>, // default: 0.25
    pub ks: Arc<dyn Texture<Spectrum> + Sync + Send>, // default: 0.25
    pub kr: Arc<dyn Texture<Spectrum> + Sync + Send>, // default: 0.0
    pub kt: Arc<dyn Texture<Spectrum> + Sync + Send>, // default: 0.0
    pub opacity: Arc<dyn Texture<Spectrum> + Sync + Send>, // default: 1.0
    pub roughness: Arc<dyn Texture<Float> + Sync + Send>, // default: 0.1
    pub u_roughness: Option<Arc<dyn Texture<Float> + Sync + Send>>,
    pub v_roughness: Option<Arc<dyn Texture<Float> + Sync + Send>>,
    pub eta: Arc<dyn Texture<Float> + Sync + Send>, // default: 1.5
    pub bump_map: Option<Arc<dyn Texture<Float> + Sync + Send>>,
    pub remap_roughness: bool,
}

impl UberMaterial {
    pub fn new(
        kd: Arc<dyn Texture<Spectrum> + Sync + Send>,
        ks: Arc<dyn Texture<Spectrum> + Sync + Send>,
        kr: Arc<dyn Texture<Spectrum> + Sync + Send>,
        kt: Arc<dyn Texture<Spectrum> + Sync + Send>,
        roughness: Arc<dyn Texture<Float> + Sync + Send>,
        u_roughness: Option<Arc<dyn Texture<Float> + Sync + Send>>,
        v_roughness: Option<Arc<dyn Texture<Float> + Sync + Send>>,
        opacity: Arc<dyn Texture<Spectrum> + Sync + Send>,
        eta: Arc<dyn Texture<Float> + Send + Sync>,
        bump_map: Option<Arc<dyn Texture<Float> + Sync + Send>>,
        remap_roughness: bool,
    ) -> Self {
        UberMaterial {
            kd,
            ks,
            kr,
            kt,
            opacity,
            roughness,
            u_roughness,
            v_roughness,
            eta,
            bump_map,
            remap_roughness,
        }
    }
    pub fn create(mp: &mut TextureParams) -> Arc<dyn Material + Send + Sync> {
        let kd: Arc<dyn Texture<Spectrum> + Sync + Send> =
            mp.get_spectrum_texture("Kd", Spectrum::new(0.25));
        let ks: Arc<dyn Texture<Spectrum> + Sync + Send> =
            mp.get_spectrum_texture("Ks", Spectrum::new(0.25));
        let kr: Arc<dyn Texture<Spectrum> + Sync + Send> =
            mp.get_spectrum_texture("Kr", Spectrum::new(0.0));
        let kt: Arc<dyn Texture<Spectrum> + Sync + Send> =
            mp.get_spectrum_texture("Kt", Spectrum::new(0.0));
        let roughness: Arc<dyn Texture<Float> + Send + Sync> =
            mp.get_float_texture("roughness", 0.1 as Float);
        let u_roughness: Option<Arc<dyn Texture<Float> + Send + Sync>> =
            mp.get_float_texture_or_null("uroughness");
        let v_roughness: Option<Arc<dyn Texture<Float> + Send + Sync>> =
            mp.get_float_texture_or_null("vroughness");
        let opacity: Arc<dyn Texture<Spectrum> + Send + Sync> =
            mp.get_spectrum_texture("opacity", Spectrum::new(1.0));
        let bump_map: Option<Arc<dyn Texture<Float> + Send + Sync>> =
            mp.get_float_texture_or_null("bumpmap");
        let remap_roughness: bool = mp.find_bool("remaproughness", true);
        let eta_option: Option<Arc<dyn Texture<Float> + Send + Sync>> =
            mp.get_float_texture_or_null("eta");
        if let Some(ref eta) = eta_option {
            Arc::new(UberMaterial::new(
                kd,
                ks,
                kr,
                kt,
                roughness,
                u_roughness,
                v_roughness,
                opacity,
                eta.clone(),
                bump_map,
                remap_roughness,
            ))
        } else {
            let eta: Arc<dyn Texture<Float> + Send + Sync> = mp.get_float_texture("index", 1.5 as Float);
            Arc::new(UberMaterial::new(
                kd,
                ks,
                kr,
                kt,
                roughness,
                u_roughness,
                v_roughness,
                opacity,
                eta,
                bump_map,
                remap_roughness,
            ))
        }
    }
}

impl Material for UberMaterial {
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        // arena: &mut Arena,
        mode: TransportMode,
        _allow_multiple_lobes: bool,
        _material: Option<Arc<dyn Material + Send + Sync>>,
    ) {
        if let Some(ref bump_map) = self.bump_map {
            Self::bump(bump_map, si);
        }
        let mut bxdfs: Vec<Arc<dyn Bxdf + Send + Sync>> = Vec::new();
        let e: Float = self.eta.evaluate(si);
        let op: Spectrum = self
            .opacity
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        let t: Spectrum =
            (Spectrum::new(1.0) - op).clamp(0.0 as Float, std::f32::INFINITY as Float);
        if !t.is_black() {
            bxdfs.push(Arc::new(SpecularTransmission::new(
                t,
                1.0,
                1.0,
                mode.clone(),
            )));
        }
        let kd: Spectrum = op * self
            .kd
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        if !kd.is_black() {
            bxdfs.push(Arc::new(LambertianReflection::new(kd)));
        }
        let ks: Spectrum = op * self
            .ks
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        if !ks.is_black() {
            let fresnel = Arc::new(FresnelDielectric {
                eta_i: 1.0,
                eta_t: e,
            });
            let mut u_rough: Float;
            if let Some(ref u_roughness) = self.u_roughness {
                u_rough = u_roughness.evaluate(si);
            } else {
                u_rough = self.roughness.evaluate(si);
            }
            let mut v_rough: Float;
            if let Some(ref v_roughness) = self.v_roughness {
                v_rough = v_roughness.evaluate(si);
            } else {
                v_rough = self.roughness.evaluate(si);
            }
            if self.remap_roughness {
                u_rough = TrowbridgeReitzDistribution::roughness_to_alpha(u_rough);
                v_rough = TrowbridgeReitzDistribution::roughness_to_alpha(v_rough);
            }
            let distrib = Arc::new(TrowbridgeReitzDistribution::new(u_rough, v_rough, true));
            bxdfs.push(Arc::new(MicrofacetReflection::new(ks, distrib, fresnel)));
        }
        let kr: Spectrum = op * self
            .kr
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        if !kr.is_black() {
            let fresnel = Arc::new(FresnelDielectric {
                eta_i: 1.0,
                eta_t: e,
            });
            bxdfs.push(Arc::new(SpecularReflection::new(kr, fresnel)));
        }
        let kt: Spectrum = op * self
            .kt
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        if !kt.is_black() {
            bxdfs.push(Arc::new(SpecularTransmission::new(
                kt,
                1.0,
                e,
                mode.clone(),
            )));
        }
        if !t.is_black() {
            si.bsdf = Some(Arc::new(Bsdf::new(si, 1.0, bxdfs)));
        } else {
            si.bsdf = Some(Arc::new(Bsdf::new(si, e, bxdfs)));
        }
    }
}
