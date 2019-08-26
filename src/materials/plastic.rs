//std
use std;
use std::sync::Arc;
// pbrt
use crate::core::interaction::SurfaceInteraction;
use crate::core::material::{Material, TransportMode};
use crate::core::microfacet::TrowbridgeReitzDistribution;
use crate::core::paramset::TextureParams;
use crate::core::pbrt::{Float, Spectrum};
use crate::core::reflection::{Bsdf, Bxdf, FresnelDielectric, LambertianReflection, MicrofacetReflection};
use crate::core::texture::Texture;

// see plastic.h

/// Plastic can be modeled as a mixture of a diffuse and glossy
/// scattering function.
pub struct PlasticMaterial {
    pub kd: Arc<dyn Texture<Spectrum> + Sync + Send>, // default: 0.25
    pub ks: Arc<dyn Texture<Spectrum> + Sync + Send>, // default: 0.25
    pub roughness: Arc<dyn Texture<Float> + Sync + Send>, // default: 0.1
    pub bump_map: Option<Arc<dyn Texture<Float> + Send + Sync>>,
    pub remap_roughness: bool,
}

impl PlasticMaterial {
    pub fn new(
        kd: Arc<dyn Texture<Spectrum> + Send + Sync>,
        ks: Arc<dyn Texture<Spectrum> + Send + Sync>,
        roughness: Arc<dyn Texture<Float> + Sync + Send>,
        bump_map: Option<Arc<dyn Texture<Float> + Sync + Send>>,
        remap_roughness: bool,
    ) -> Self {
        PlasticMaterial {
            kd,
            ks,
            roughness,
            bump_map,
            remap_roughness,
        }
    }
    pub fn create(mp: &mut TextureParams) -> Arc<dyn Material + Send + Sync> {
        let kd = mp.get_spectrum_texture("Kd", Spectrum::new(0.25 as Float));
        let ks = mp.get_spectrum_texture("Ks", Spectrum::new(0.25 as Float));
        let roughness = mp.get_float_texture("roughness", 0.1 as Float);
        let bump_map = mp.get_float_texture_or_null("bumpmap");
        let remap_roughness: bool = mp.find_bool("remaproughness", true);
        Arc::new(PlasticMaterial::new(
            kd,
            ks,
            roughness,
            bump_map,
            remap_roughness,
        ))
    }
}

impl Material for PlasticMaterial {
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
        let mut bxdfs: Vec<Arc<dyn Bxdf + Send + Sync>> = Vec::new();
        // initialize diffuse component of plastic material
        let kd: Spectrum = self
            .kd
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        if !kd.is_black() {
            bxdfs.push(Arc::new(LambertianReflection::new(kd)));
        }
        // initialize specular component of plastic material
        let ks: Spectrum = self
            .ks
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        if !ks.is_black() {
            let fresnel = Arc::new(FresnelDielectric {
                eta_i: 1.5 as Float,
                eta_t: 1.0 as Float,
            });
            // create microfacet distribution _distrib_ for plastic material
            let mut rough: Float = self.roughness.evaluate(si);
            if self.remap_roughness {
                rough = TrowbridgeReitzDistribution::roughness_to_alpha(rough);
            }
            let distrib = Arc::new(TrowbridgeReitzDistribution::new(rough, rough, true));
            bxdfs.push(Arc::new(MicrofacetReflection::new(ks, distrib, fresnel)));
        }
        si.bsdf = Some(Arc::new(Bsdf::new(si, 1.0, bxdfs)));
    }
}
