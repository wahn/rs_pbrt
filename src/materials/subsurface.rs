//std
use std;
use std::sync::Arc;
// pbrt
use core::interaction::SurfaceInteraction;
use core::material::{Material, TransportMode};
use core::paramset::TextureParams;
use core::reflection::{Bsdf, Bxdf, LambertianReflection, OrenNayar};

// see subsurface.h

pub struct SubsurfaceMaterial {
}

impl SubsurfaceMaterial {
    pub fn new(
        ) -> Self {
        SubsurfaceMaterial {
        }
    }
    pub fn create(mp: &mut TextureParams) -> Arc<Material + Send + Sync> {
        Arc::new(SubsurfaceMaterial {})
    }
    pub fn bsdf(&self, si: &SurfaceInteraction) -> Bsdf {
        let mut bxdfs: Vec<Arc<Bxdf + Send + Sync>> = Vec::new();
        Bsdf::new(si, 1.0, bxdfs)
    }
}

impl Material for SubsurfaceMaterial {
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
