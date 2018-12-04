//std
use std;
use std::sync::Arc;
// pbrt
use core::interaction::SurfaceInteraction;
use core::material::{Material, TransportMode};
use core::paramset::TextureParams;
use core::pbrt::clamp_t;
use core::pbrt::{Float, Spectrum};
use core::reflection::{Bsdf, Bxdf, LambertianReflection, OrenNayar};
use core::texture::Texture;

// see matte.h

/// Describes a purely diffuse surface.
pub struct MatteMaterial {
    pub kd: Arc<Texture<Spectrum> + Sync + Send>, // default: 0.5
    pub sigma: Arc<Texture<Float> + Sync + Send>, // default: 0.0
                                                  // TODO: bump_map
}

impl MatteMaterial {
    pub fn new(
        kd: Arc<Texture<Spectrum> + Send + Sync>,
        sigma: Arc<Texture<Float> + Sync + Send>,
    ) -> Self {
        MatteMaterial {
            kd: kd,
            sigma: sigma,
        }
    }
    pub fn create(mp: &mut TextureParams) -> Arc<Material + Send + Sync> {
        let kd: Arc<Texture<Spectrum> + Sync + Send> =
            mp.get_spectrum_texture(String::from("Kd"), Spectrum::new(0.5));
        let sigma: Arc<Texture<Float> + Sync + Send> =
            mp.get_float_texture(String::from("sigma"), 0.0);
        Arc::new(MatteMaterial { kd, sigma })
    }
    pub fn bsdf(&self, si: &SurfaceInteraction) -> Bsdf {
        let mut bxdfs: Vec<Arc<Bxdf + Send + Sync>> = Vec::new();
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
        Bsdf::new(si, 1.0, bxdfs)
    }
}

impl Material for MatteMaterial {
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        // arena: &mut Arena,
        _mode: TransportMode,
        _allow_multiple_lobes: bool,
    ) {
        si.bsdf = Some(Arc::new(self.bsdf(si)));
    }
}
