//std
use std;
use std::sync::Arc;
// pbrt
use core::interaction::SurfaceInteraction;
use core::material::{Material, TransportMode};
use core::pbrt::{Float, Spectrum};
use core::reflection::{Bsdf, Bxdf, LambertianReflection, OrenNayar};
use core::texture::Texture;

// see matte.h

/// Describes a purely diffuse surface.
pub struct MatteMaterial {
    pub kd: Arc<Texture<Spectrum> + Sync + Send>, // default: 0.5
    pub sigma: Float, // default: 0.0
                      // TODO: bump_map
}

impl MatteMaterial {
    pub fn new(kd: Arc<Texture<Spectrum> + Send + Sync>, sigma: Float) -> Self {
        MatteMaterial {
            kd: kd,
            sigma: sigma,
        }
    }
    pub fn bsdf(&self, si: &SurfaceInteraction) -> Bsdf {
        let mut bxdfs: Vec<Arc<Bxdf + Send + Sync>> = Vec::new();
        let r: Spectrum = self.kd
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        if !r.is_black() {
            if self.sigma == 0.0 {
                bxdfs.push(Arc::new(LambertianReflection::new(r)));
            } else {
                bxdfs.push(Arc::new(OrenNayar::new(r, self.sigma)));
            }
        }
        Bsdf::new(si, 1.0, bxdfs)
    }
}

impl Material for MatteMaterial {
    fn compute_scattering_functions(&self,
                                    si: &mut SurfaceInteraction,
                                    // arena: &mut Arena,
                                    _mode: TransportMode,
                                    _allow_multiple_lobes: bool) {
        si.bsdf = Some(Arc::new(self.bsdf(si)));
    }
}
