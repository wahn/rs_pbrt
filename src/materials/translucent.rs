//std
use std::sync::Arc;
// pbrt
use crate::core::interaction::SurfaceInteraction;
use crate::core::material::{Material, TransportMode};
use crate::core::microfacet::TrowbridgeReitzDistribution;
use crate::core::paramset::TextureParams;
use crate::core::pbrt::{Float, Spectrum};
use crate::core::reflection::{
    Bsdf, Bxdf, FresnelDielectric, LambertianReflection, LambertianTransmission,
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
    ) {
        if let Some(ref bump) = self.bump_map {
            Self::bump(bump, si);
        }
        let mut bxdfs: Vec<Arc<dyn Bxdf + Send + Sync>> = Vec::new();
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
            si.bsdf = Some(Arc::new(Bsdf::new(si, eta, bxdfs)));
            return;
        }
        let kd: Spectrum = self
            .kd
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        if !kd.is_black() {
            if !r.is_black() {
                bxdfs.push(Arc::new(LambertianReflection::new(r * kd)));
            }
            if !t.is_black() {
                bxdfs.push(Arc::new(LambertianTransmission::new(t * kd)));
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
                let fresnel = Arc::new(FresnelDielectric {
                    eta_i: 1.0 as Float,
                    eta_t: eta,
                });
                bxdfs.push(Arc::new(MicrofacetReflection::new(
                    r * ks,
                    distrib.clone(),
                    fresnel,
                )));
            }
            if !t.is_black() {
                bxdfs.push(Arc::new(MicrofacetTransmission::new(
                    t * ks,
                    distrib.clone(),
                    1.0,
                    eta,
                    mode,
                )));
            }
        }
        si.bsdf = Some(Arc::new(Bsdf::new(si, eta, bxdfs)));
    }
}
