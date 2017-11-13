//std
use std;
use std::sync::Arc;
// pbrt
use core::interaction::SurfaceInteraction;
use core::material::TransportMode;
use core::pbrt::{Float, Spectrum};
use core::reflection::{Bsdf, Bxdf, FresnelNoOp, SpecularReflection};
use materials::Material;
use textures::Texture;

// see mirror.h

/// A simple mirror, modeled with perfect specular reflection.
pub struct MirrorMaterial {
    pub kr: Arc<Texture<Spectrum> + Sync + Send>, // default: 0.9
                                                  // TODO: bump_map
}

impl MirrorMaterial {
    pub fn new(kr: Arc<Texture<Spectrum> + Send + Sync>) -> Self {
        MirrorMaterial { kr: kr }
    }
    pub fn bsdf(&self, si: &SurfaceInteraction) -> Bsdf {
        let mut bxdfs: Vec<Box<Bxdf + Send + Sync>> = Vec::new();
        let r: Spectrum = self.kr
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        let fresnel = Arc::new(FresnelNoOp {});
        bxdfs.push(Box::new(SpecularReflection::new(r, fresnel)));
        Bsdf::new(si, 1.5, bxdfs)
    }
}

impl Material for MirrorMaterial {
    fn compute_scattering_functions(&self,
                                    si: &mut SurfaceInteraction,
                                    // arena: &mut Arena,
                                    _mode: TransportMode,
                                    _allow_multiple_lobes: bool) {
        si.bsdf = Some(Arc::new(self.bsdf(si)));
    }
}
