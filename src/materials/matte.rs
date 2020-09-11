//std
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
    pub fn create(mp: &mut TextureParams) -> Arc<Material> {
        let kd: Arc<dyn Texture<Spectrum> + Sync + Send> =
            mp.get_spectrum_texture("Kd", Spectrum::new(0.5));
        let sigma: Arc<dyn Texture<Float> + Sync + Send> = mp.get_float_texture("sigma", 0.0);
        let bump_map = mp.get_float_texture_or_null("bumpmap");
        Arc::new(Material::Matte(Box::new(MatteMaterial::new(
            kd, sigma, bump_map,
        ))))
    }
    // Material
    pub fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        // arena: &mut Arena,
        _mode: TransportMode,
        _allow_multiple_lobes: bool,
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
        let r: Spectrum = self
            .kd
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        let sig: Float = clamp_t(
            self.sigma.evaluate(si) as Float,
            0.0 as Float,
            90.0 as Float,
        );
        si.bsdf = Some(Bsdf::new(si, 1.0));
        if let Some(bsdf) = &mut si.bsdf {
            if !r.is_black() {
                if sig == 0.0 {
                    if use_scale {
                        bsdf.add(Bxdf::LambertianRefl(LambertianReflection::new(r, Some(sc))));
                    } else {
                        bsdf.add(Bxdf::LambertianRefl(LambertianReflection::new(r, None)));
                    }
                } else if use_scale {
                    bsdf.add(Bxdf::OrenNayarRefl(OrenNayar::new(r, sig, Some(sc))));
                } else {
                    bsdf.add(Bxdf::OrenNayarRefl(OrenNayar::new(r, sig, None)));
                }
            }
        }
    }
}
