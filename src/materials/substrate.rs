//std
use std;
use std::sync::Arc;
// pbrt
use crate::core::interaction::SurfaceInteraction;
use crate::core::material::{Material, TransportMode};
use crate::core::microfacet::TrowbridgeReitzDistribution;
use crate::core::paramset::TextureParams;
use crate::core::pbrt::{Float, Spectrum};
use crate::core::reflection::{Bsdf, Bxdf, FresnelBlend};
use crate::core::texture::Texture;

// see substrate.h

pub struct SubstrateMaterial {
    pub kd: Arc<dyn Texture<Spectrum> + Sync + Send>, // default: 0.5
    pub ks: Arc<dyn Texture<Spectrum> + Sync + Send>, // default: 0.5
    pub nu: Arc<dyn Texture<Float> + Sync + Send>,    // default: 0.1
    pub nv: Arc<dyn Texture<Float> + Sync + Send>,    // default: 0.1
    pub bump_map: Option<Arc<dyn Texture<Float> + Send + Sync>>,
    pub remap_roughness: bool,
}

impl SubstrateMaterial {
    pub fn new(
        kd: Arc<dyn Texture<Spectrum> + Send + Sync>,
        ks: Arc<dyn Texture<Spectrum> + Send + Sync>,
        nu: Arc<dyn Texture<Float> + Sync + Send>,
        nv: Arc<dyn Texture<Float> + Sync + Send>,
        bump_map: Option<Arc<dyn Texture<Float> + Sync + Send>>,
        remap_roughness: bool,
    ) -> Self {
        SubstrateMaterial {
            kd,
            ks,
            nu,
            nv,
            bump_map,
            remap_roughness,
        }
    }
    pub fn create(mp: &mut TextureParams) -> Arc<dyn Material + Send + Sync> {
        let kd: Arc<dyn Texture<Spectrum> + Sync + Send> =
            mp.get_spectrum_texture("Kd", Spectrum::new(0.5));
        let ks: Arc<dyn Texture<Spectrum> + Sync + Send> =
            mp.get_spectrum_texture("Ks", Spectrum::new(0.5));
        let uroughness: Arc<dyn Texture<Float> + Sync + Send> =
            mp.get_float_texture("uroughness", 0.1);
        let vroughness: Arc<dyn Texture<Float> + Sync + Send> =
            mp.get_float_texture("vroughness", 0.1);
        let bump_map = mp.get_float_texture_or_null("bumpmap");
        let remap_roughness: bool = mp.find_bool("remaproughness", true);
        Arc::new(SubstrateMaterial::new(
            kd,
            ks,
            uroughness,
            vroughness,
            bump_map,
            remap_roughness,
        ))
    }
}

impl Material for SubstrateMaterial {
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        // arena: &mut Arena,
        _mode: TransportMode,
        _allow_multiple_lobes: bool,
        _material: Option<Arc<dyn Material + Send + Sync>>,
    ) {
        if let Some(ref bump) = self.bump_map {
            Self::bump(bump, si);
        }
        let mut bxdfs: Vec<Bxdf> = Vec::new();
        let d: Spectrum = self
            .kd
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        let s: Spectrum = self
            .ks
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        let mut roughu: Float = self.nu.evaluate(si);
        let mut roughv: Float = self.nv.evaluate(si);
        if !d.is_black() || !s.is_black() {
            if self.remap_roughness {
                roughu = TrowbridgeReitzDistribution::roughness_to_alpha(roughu);
                roughv = TrowbridgeReitzDistribution::roughness_to_alpha(roughv);
            }
            let distrib: Option<TrowbridgeReitzDistribution> =
                Some(TrowbridgeReitzDistribution::new(roughu, roughv, true));
            bxdfs.push(Bxdf::FresnelBlnd(FresnelBlend::new(d, s, distrib)));
        }
        si.bsdf = Some(Arc::new(Bsdf::new(si, 1.0, bxdfs)));
    }
}
