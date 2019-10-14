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
        scale_opt: Option<Spectrum>,
    ) {
        let mut use_scale: bool = false;
        let mut sc: Spectrum = Spectrum::default();
        if let Some(scale) = scale_opt {
            use_scale = true;
            sc = scale;
        }
        if let Some(ref bump) = self.bump_map {
            Self::bump(bump, si);
        }
        let r: Spectrum = self
            .kr
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        si.bsdf = Some(Bsdf::new(si, 1.0));
        if let Some(bsdf) = &mut si.bsdf {
            let bxdf_idx: usize = 0;
            let fresnel = Fresnel::NoOp(FresnelNoOp {});
            if use_scale {
                bsdf.bxdfs[bxdf_idx] =
                    Bxdf::SpecRefl(SpecularReflection::new(r, fresnel, Some(sc)));
            } else {
                bsdf.bxdfs[bxdf_idx] = Bxdf::SpecRefl(SpecularReflection::new(r, fresnel, None));
            }
        }
    }
}
