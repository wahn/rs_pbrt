//std
use std;
use std::sync::Arc;
// pbrt
use core::interaction::SurfaceInteraction;
use core::material::{Material, TransportMode};
use core::microfacet::TrowbridgeReitzDistribution;
use core::paramset::TextureParams;
use core::pbrt::{Float, Spectrum};
use core::reflection::{
    Bsdf, Bxdf, FresnelDielectric, FresnelSpecular, MicrofacetReflection, MicrofacetTransmission,
    SpecularReflection, SpecularTransmission,
};
use core::texture::Texture;

// see glass.h

/// Perfect or glossy specular reflection and transmission, weighted
/// by Fresnel terms for accurate angular-dependent variation.
pub struct GlassMaterial {
    pub kr: Arc<Texture<Spectrum> + Sync + Send>, // default: 1.0
    pub kt: Arc<Texture<Spectrum> + Sync + Send>, // default: 1.0
    pub u_roughness: Arc<Texture<Float> + Sync + Send>, // default: 0.0
    pub v_roughness: Arc<Texture<Float> + Sync + Send>, // default: 0.0
    pub index: Arc<Texture<Float> + Sync + Send>,
    pub bump_map: Option<Arc<Texture<Float> + Send + Sync>>,
    pub remap_roughness: bool,
}

impl GlassMaterial {
    pub fn new(
        kr: Arc<Texture<Spectrum> + Sync + Send>,
        kt: Arc<Texture<Spectrum> + Sync + Send>,
        u_roughness: Arc<Texture<Float> + Sync + Send>,
        v_roughness: Arc<Texture<Float> + Sync + Send>,
        index: Arc<Texture<Float> + Send + Sync>,
        bump_map: Option<Arc<Texture<Float> + Sync + Send>>,
        remap_roughness: bool,
    ) -> Self {
        GlassMaterial {
            kr: kr,
            kt: kt,
            u_roughness: u_roughness,
            v_roughness: v_roughness,
            index: index,
            bump_map: bump_map,
            remap_roughness: remap_roughness,
        }
    }
    pub fn create(mp: &mut TextureParams) -> Arc<Material + Send + Sync> {
        let kr = mp.get_spectrum_texture("Kr", Spectrum::new(1.0 as Float));
        let kt = mp.get_spectrum_texture("Kt", Spectrum::new(1.0 as Float));
        let roughu = mp.get_float_texture("uroughness", 0.0 as Float);
        let roughv = mp.get_float_texture("vroughness", 0.0 as Float);
        let bump_map = mp.get_float_texture_or_null("bumpmap");
        let remap_roughness: bool = mp.find_bool("remaproughness", true);
        let eta_option: Option<Arc<Texture<Float> + Send + Sync>> =
            mp.get_float_texture_or_null("eta");
        if let Some(ref eta) = eta_option {
            Arc::new(GlassMaterial::new(
                kr,
                kt,
                roughu,
                roughv,
                eta.clone(),
                bump_map,
                remap_roughness,
            ))
        } else {
            let eta: Arc<Texture<Float> + Send + Sync> = mp.get_float_texture("eta", 1.5 as Float);
            Arc::new(GlassMaterial::new(
                kr,
                kt,
                roughu,
                roughv,
                eta,
                bump_map,
                remap_roughness,
            ))
        }
    }
}

impl Material for GlassMaterial {
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        // arena: &mut Arena,
        mode: TransportMode,
        allow_multiple_lobes: bool,
    ) {
        if let Some(ref bump_map) = self.bump_map {
            Self::bump(bump_map, si);
        }
        let mut bxdfs: Vec<Arc<Bxdf + Send + Sync>> = Vec::new();
        let eta: Float = self.index.evaluate(si);
        let mut urough: Float = self.u_roughness.evaluate(si);
        let mut vrough: Float = self.v_roughness.evaluate(si);
        let r: Spectrum = self
            .kr
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        let t: Spectrum = self
            .kt
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        let is_specular: bool = urough == 0.0 as Float && vrough == 0.0 as Float;
        if is_specular && allow_multiple_lobes {
            bxdfs.push(Arc::new(FresnelSpecular::new(
                r,
                t,
                1.0 as Float,
                eta,
                mode,
            )));
        } else {
            if self.remap_roughness {
                urough = TrowbridgeReitzDistribution::roughness_to_alpha(urough);
                vrough = TrowbridgeReitzDistribution::roughness_to_alpha(vrough);
            }
            if !r.is_black() {
                let fresnel = Arc::new(FresnelDielectric {
                    eta_i: 1.0 as Float,
                    eta_t: eta,
                });
                if is_specular {
                    bxdfs.push(Arc::new(SpecularReflection::new(r, fresnel)));
                } else {
                    let distrib = Arc::new(TrowbridgeReitzDistribution::new(urough, vrough, true));
                    bxdfs.push(Arc::new(MicrofacetReflection::new(r, distrib, fresnel)));
                }
            }
            if !t.is_black() {
                if is_specular {
                    bxdfs.push(Arc::new(SpecularTransmission::new(t, 1.0, eta, mode)));
                } else {
                    let distrib = Arc::new(TrowbridgeReitzDistribution::new(urough, vrough, true));
                    bxdfs.push(Arc::new(MicrofacetTransmission::new(
                        t, distrib, 1.0, eta, mode,
                    )));
                }
            }
        }
        si.bsdf = Some(Arc::new(Bsdf::new(si, eta, bxdfs)));
    }
}
