//std
use std;
use std::sync::Arc;
// pbrt
use crate::core::interaction::SurfaceInteraction;
use crate::core::material::{Material, TransportMode};
use crate::core::paramset::TextureParams;
use crate::core::pbrt::{Float, Spectrum};
use crate::core::reflection::{Bsdf, Bxdf, Fresnel, FresnelNoOp, SpecularReflection};
use crate::core::texture::Texture;

// see mirror.h

/// A simple mirror, modeled with perfect specular reflection.
pub struct MirrorMaterial {
    pub kr: Arc<dyn Texture<Spectrum> + Sync + Send>, // default: 0.9
    pub bump_map: Option<Arc<dyn Texture<Float> + Send + Sync>>,
}

impl MirrorMaterial {
    pub fn new(
        kr: Arc<dyn Texture<Spectrum> + Send + Sync>,
        bump_map: Option<Arc<dyn Texture<Float> + Sync + Send>>,
    ) -> Self {
        MirrorMaterial { kr, bump_map }
    }
    pub fn create(mp: &mut TextureParams) -> Arc<dyn Material + Send + Sync> {
        let kr = mp.get_spectrum_texture("Kr", Spectrum::new(0.9 as Float));
        let bump_map = mp.get_float_texture_or_null("bumpmap");
        Arc::new(MirrorMaterial::new(kr, bump_map))
    }
}

impl Material for MirrorMaterial {
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        // arena: &mut Arena,
        _mode: TransportMode,
        _allow_multiple_lobes: bool,
        _material: Option<Arc<dyn Material + Send + Sync>>,
    ) -> Vec<Bxdf> {
        if let Some(ref bump) = self.bump_map {
            Self::bump(bump, si);
        }
        let mut bxdfs: Vec<Bxdf> = Vec::new();
        let r: Spectrum = self
            .kr
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        let fresnel = Fresnel::NoOp(FresnelNoOp {});
        bxdfs.push(Bxdf::SpecRefl(SpecularReflection::new(r, fresnel)));
        si.bsdf = Some(Arc::new(Bsdf::new(si, 1.0, Vec::new())));
        bxdfs
    }
}
