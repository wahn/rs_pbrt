//std
use std::sync::Arc;
// pbrt
use crate::core::interaction::SurfaceInteraction;
use crate::core::material::{Material, TransportMode};
use crate::core::microfacet::{MicrofacetDistribution, TrowbridgeReitzDistribution};
use crate::core::paramset::TextureParams;
use crate::core::pbrt::{Float, Spectrum};
use crate::core::reflection::{
    Bsdf, Bxdf, Fresnel, FresnelDielectric, FresnelSpecular, MicrofacetReflection,
    MicrofacetTransmission, SpecularReflection, SpecularTransmission,
};
use crate::core::texture::Texture;

// see glass.h

/// Perfect or glossy specular reflection and transmission, weighted
/// by Fresnel terms for accurate angular-dependent variation.
pub struct GlassMaterial {
    pub kr: Arc<dyn Texture<Spectrum> + Sync + Send>, // default: 1.0
    pub kt: Arc<dyn Texture<Spectrum> + Sync + Send>, // default: 1.0
    pub u_roughness: Arc<dyn Texture<Float> + Sync + Send>, // default: 0.0
    pub v_roughness: Arc<dyn Texture<Float> + Sync + Send>, // default: 0.0
    pub index: Arc<dyn Texture<Float> + Sync + Send>,
    pub bump_map: Option<Arc<dyn Texture<Float> + Send + Sync>>,
    pub remap_roughness: bool,
}

impl GlassMaterial {
    pub fn new(
        kr: Arc<dyn Texture<Spectrum> + Sync + Send>,
        kt: Arc<dyn Texture<Spectrum> + Sync + Send>,
        u_roughness: Arc<dyn Texture<Float> + Sync + Send>,
        v_roughness: Arc<dyn Texture<Float> + Sync + Send>,
        index: Arc<dyn Texture<Float> + Send + Sync>,
        bump_map: Option<Arc<dyn Texture<Float> + Sync + Send>>,
        remap_roughness: bool,
    ) -> Self {
        GlassMaterial {
            kr,
            kt,
            u_roughness,
            v_roughness,
            index,
            bump_map,
            remap_roughness,
        }
    }
    pub fn create(mp: &mut TextureParams) -> Arc<Material> {
        let kr = mp.get_spectrum_texture("Kr", Spectrum::new(1.0 as Float));
        let kt = mp.get_spectrum_texture("Kt", Spectrum::new(1.0 as Float));
        let roughu = mp.get_float_texture("uroughness", 0.0 as Float);
        let roughv = mp.get_float_texture("vroughness", 0.0 as Float);
        let bump_map = mp.get_float_texture_or_null("bumpmap");
        let remap_roughness: bool = mp.find_bool("remaproughness", true);
        let eta_option: Option<Arc<dyn Texture<Float> + Send + Sync>> =
            mp.get_float_texture_or_null("eta");
        if let Some(ref eta) = eta_option {
            Arc::new(Material::Glass(Box::new(GlassMaterial::new(
                kr,
                kt,
                roughu,
                roughv,
                eta.clone(),
                bump_map,
                remap_roughness,
            ))))
        } else {
            let eta: Arc<dyn Texture<Float> + Send + Sync> =
                mp.get_float_texture("index", 1.5 as Float);
            Arc::new(Material::Glass(Box::new(GlassMaterial::new(
                kr,
                kt,
                roughu,
                roughv,
                eta,
                bump_map,
                remap_roughness,
            ))))
        }
    }
    // Material
    pub fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        // arena: &mut Arena,
        mode: TransportMode,
        allow_multiple_lobes: bool,
        _material: Option<Arc<Material>>,
        scale_opt: Option<Spectrum>,
    ) {
        let mut use_scale: bool = false;
        let mut sc: Spectrum = Spectrum::default();
        if let Some(scale) = scale_opt {
            use_scale = true;
            sc = scale;
        }
        if let Some(ref bump) = self.bump_map {
            Material::bump(bump, si);
        }
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
        let eta: Float = self.index.evaluate(si);
        si.bsdf = Some(Bsdf::new(si, eta));
        if let Some(bsdf) = &mut si.bsdf {
            if is_specular && allow_multiple_lobes {
                if use_scale {
                    bsdf.add(Bxdf::FresnelSpec(FresnelSpecular::new(
                        r,
                        t,
                        1.0 as Float,
                        eta,
                        mode,
                        Some(sc),
                    )));
                } else {
                    bsdf.add(Bxdf::FresnelSpec(FresnelSpecular::new(
                        r,
                        t,
                        1.0 as Float,
                        eta,
                        mode,
                        None,
                    )));
                }
            } else {
                if self.remap_roughness {
                    urough = TrowbridgeReitzDistribution::roughness_to_alpha(urough);
                    vrough = TrowbridgeReitzDistribution::roughness_to_alpha(vrough);
                }
                if !r.is_black() {
                    let fresnel = Fresnel::Dielectric(FresnelDielectric {
                        eta_i: 1.0 as Float,
                        eta_t: eta,
                    });
                    if is_specular {
                        if use_scale {
                            bsdf.add(Bxdf::SpecRefl(SpecularReflection::new(
                                r,
                                fresnel,
                                Some(sc),
                            )));
                        } else {
                            bsdf.add(Bxdf::SpecRefl(SpecularReflection::new(r, fresnel, None)));
                        }
                    } else {
                        let distrib = MicrofacetDistribution::TrowbridgeReitz(
                            TrowbridgeReitzDistribution::new(urough, vrough, true),
                        );
                        if use_scale {
                            bsdf.add(Bxdf::MicrofacetRefl(MicrofacetReflection::new(
                                r,
                                distrib,
                                fresnel,
                                Some(sc),
                            )));
                        } else {
                            bsdf.add(Bxdf::MicrofacetRefl(MicrofacetReflection::new(
                                r, distrib, fresnel, None,
                            )));
                        }
                    }
                }
                if !t.is_black() {
                    if is_specular {
                        if use_scale {
                            bsdf.add(Bxdf::SpecTrans(SpecularTransmission::new(
                                t,
                                1.0,
                                eta,
                                mode,
                                Some(sc),
                            )));
                        } else {
                            bsdf.add(Bxdf::SpecTrans(SpecularTransmission::new(
                                t, 1.0, eta, mode, None,
                            )));
                        }
                    } else {
                        let distrib = MicrofacetDistribution::TrowbridgeReitz(
                            TrowbridgeReitzDistribution::new(urough, vrough, true),
                        );
                        if use_scale {
                            bsdf.add(Bxdf::MicrofacetTrans(MicrofacetTransmission::new(
                                t,
                                distrib,
                                1.0,
                                eta,
                                mode,
                                Some(sc),
                            )));
                        } else {
                            bsdf.add(Bxdf::MicrofacetTrans(MicrofacetTransmission::new(
                                t, distrib, 1.0, eta, mode, None,
                            )));
                        }
                    }
                }
            }
        }
    }
}
