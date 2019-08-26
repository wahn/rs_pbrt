//std
use std;
use std::sync::Arc;
// pbrt
use crate::core::interaction::SurfaceInteraction;
use crate::core::material::{Material, TransportMode};
use crate::core::paramset::TextureParams;
use crate::core::pbrt::clamp_t;
use crate::core::pbrt::{Float, Spectrum};
use crate::core::reflection::{Bsdf, Bxdf, LambertianReflection, OrenNayar};
use crate::core::texture::Texture;

// see matte.h

/// Describes a purely diffuse surface.
pub struct MatteMaterial {
    pub kd: Arc<dyn Texture<Spectrum> + Sync + Send>, // default: 0.5
    pub sigma: Arc<dyn Texture<Float> + Sync + Send>, // default: 0.0
    pub bump_map: Option<Arc<dyn Texture<Float> + Send + Sync>>,
}

impl MatteMaterial {
    pub fn new(
        kd: Arc<dyn Texture<Spectrum> + Send + Sync>,
        sigma: Arc<dyn Texture<Float> + Sync + Send>,
        bump_map: Option<Arc<dyn Texture<Float> + Sync + Send>>,
    ) -> Self {
        MatteMaterial {
            kd,
            sigma,
            bump_map,
        }
    }
    pub fn create(mp: &mut TextureParams) -> Arc<dyn Material + Send + Sync> {
        let kd: Arc<dyn Texture<Spectrum> + Sync + Send> =
            mp.get_spectrum_texture("Kd", Spectrum::new(0.5));
        let sigma: Arc<dyn Texture<Float> + Sync + Send> = mp.get_float_texture("sigma", 0.0);
        let bump_map = mp.get_float_texture_or_null("bumpmap");
        Arc::new(MatteMaterial::new(
            kd,
            sigma,
            bump_map,
        ))
    }
}

impl Material for MatteMaterial {
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
        let r: Spectrum = self
            .kd
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        let sig: Float = clamp_t(
            self.sigma.evaluate(si) as Float,
            0.0 as Float,
            90.0 as Float,
        );
        if !r.is_black() {
            if sig == 0.0 {
                bxdfs.push(Arc::new(LambertianReflection::new(r)));
            } else {
                bxdfs.push(Arc::new(OrenNayar::new(r, sig)));
            }
        }
        si.bsdf = Some(Arc::new(Bsdf::new(si, 1.0, bxdfs)));
    }
}
