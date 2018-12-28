//std
use std;
use std::sync::Arc;
// pbrt
use core::interaction::SurfaceInteraction;
use core::material::{Material, TransportMode};
use core::paramset::TextureParams;
use core::pbrt::{Float, Spectrum};
use core::reflection::{Bsdf, Bxdf, FresnelNoOp, SpecularReflection};
use core::texture::Texture;

// see mirror.h

/// A simple mirror, modeled with perfect specular reflection.
pub struct MirrorMaterial {
    pub kr: Arc<Texture<Spectrum> + Sync + Send>, // default: 0.9
    pub bump_map: Option<Arc<Texture<Float> + Send + Sync>>,
}

impl MirrorMaterial {
    pub fn new(
        kr: Arc<Texture<Spectrum> + Send + Sync>,
        bump_map: Option<Arc<Texture<Float> + Sync + Send>>,
    ) -> Self {
        MirrorMaterial {
            kr: kr,
            bump_map: bump_map,
        }
    }
    pub fn create(mp: &mut TextureParams) -> Arc<Material + Send + Sync> {
        let kr = mp.get_spectrum_texture("Kr", Spectrum::new(0.9 as Float));
        let bump_map = mp.get_float_texture_or_null("bumpmap");
        Arc::new(MirrorMaterial::new(
            kr,
            bump_map,
        ))
    }
}

impl Material for MirrorMaterial {
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        // arena: &mut Arena,
        _mode: TransportMode,
        _allow_multiple_lobes: bool,
    ) {
        if let Some(ref bump) = self.bump_map {
            Self::bump(bump, si);
        }
        let mut bxdfs: Vec<Arc<Bxdf + Send + Sync>> = Vec::new();
        let r: Spectrum = self
            .kr
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        let fresnel = Arc::new(FresnelNoOp {});
        bxdfs.push(Arc::new(SpecularReflection::new(r, fresnel)));
        si.bsdf = Some(Arc::new(Bsdf::new(si, 1.0, bxdfs)));
    }
}
