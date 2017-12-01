//std
use std;
use std::sync::Arc;
// pbrt
use core::interaction::SurfaceInteraction;
use core::material::{Material, TransportMode};
use core::microfacet::TrowbridgeReitzDistribution;
use core::pbrt::{Float, Spectrum};
use core::reflection::{Bsdf, Bxdf, FresnelDielectric, FresnelSpecular, MicrofacetReflection,
                       SpecularReflection, SpecularTransmission};
use core::texture::Texture;

// see glass.h

/// Perfect or glossy specular reflection and transmission, weighted
/// by Fresnel terms for accurate angular-dependent variation.
pub struct GlassMaterial {
    pub kr: Arc<Texture<Spectrum> + Sync + Send>, // default: 1.0
    pub kt: Arc<Texture<Spectrum> + Sync + Send>, // default: 1.0
    pub u_roughness: Arc<Texture<Float> + Sync + Send>, // default: 0.0
    pub v_roughness: Arc<Texture<Float> + Sync + Send>, // default: 0.0
    pub index: Arc<Texture<Float> + Sync + Send>, // TODO: bump_map
    pub remap_roughness: bool,
}

impl GlassMaterial {
    pub fn bsdf(&self,
                si: &SurfaceInteraction,
                mode: TransportMode,
                allow_multiple_lobes: bool)
                -> Bsdf {
        let mut bxdfs: Vec<Box<Bxdf + Send + Sync>> = Vec::new();
        let eta: Float = self.index.evaluate(si);
        let mut urough: Float = self.u_roughness.evaluate(si);
        let mut vrough: Float = self.v_roughness.evaluate(si);
        let r: Spectrum = self.kr
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        let t: Spectrum = self.kt
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        let is_specular: bool = urough == 0.0 as Float && vrough == 0.0 as Float;
        if is_specular && allow_multiple_lobes {
            bxdfs.push(Box::new(FresnelSpecular::new(r, t, 1.0 as Float, eta, mode)));
        } else {
            if self.remap_roughness {
                urough = TrowbridgeReitzDistribution::roughness_to_alpha(urough);
                vrough = TrowbridgeReitzDistribution::roughness_to_alpha(vrough);
            }
            let distrib: Option<TrowbridgeReitzDistribution> = match is_specular {
                true => None,
                false => Some(TrowbridgeReitzDistribution::new(urough, vrough, true)),
            };
            if !r.is_black() {
                let fresnel = Arc::new(FresnelDielectric {
                                           eta_i: 1.0 as Float,
                                           eta_t: eta,
                                       });
                if is_specular {
                    bxdfs.push(Box::new(SpecularReflection::new(r, fresnel)));
                } else {
                    bxdfs.push(Box::new(MicrofacetReflection::new(r, distrib, fresnel)));
                }
            }
            if !t.is_black() {
                if is_specular {
                    bxdfs.push(Box::new(SpecularTransmission::new(t, 1.0, eta, mode)));
                } else {
                    // TODO: si->bsdf->Add(ARENA_ALLOC(arena, MicrofacetTransmission)(
                    // T, distrib, 1.f, eta, mode));
                }
            }
        }
        Bsdf::new(si, eta, bxdfs)
    }
}

impl Material for GlassMaterial {
    fn compute_scattering_functions(&self,
                                    si: &mut SurfaceInteraction,
                                    // arena: &mut Arena,
                                    mode: TransportMode,
                                    allow_multiple_lobes: bool) {
        si.bsdf = Some(Arc::new(self.bsdf(si, mode, allow_multiple_lobes)));
    }
}
