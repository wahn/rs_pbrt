//std
use std::collections::HashMap;
use std::sync::Arc;
// pbrt
use core::interaction::SurfaceInteraction;
use core::material::{Material, TransportMode};
use core::pbrt::Float;
use core::reflection::{Bsdf, Bxdf, FourierBSDFTable};
use core::texture::Texture;

// see fourier.h

pub struct FourierMaterial {
    pub bsdf_table: FourierBSDFTable,
    pub bump_map: Option<Arc<Texture<Float> + Sync + Send>>,
    pub loaded_bsdfs: HashMap<String, Arc<FourierBSDFTable>>,
}

impl FourierMaterial {
    pub fn bsdf(&self, si: &SurfaceInteraction) -> Bsdf {
        let // mut 
            bxdfs: Vec<Arc<Bxdf + Send + Sync>> = Vec::new();
        // WORK
        Bsdf::new(si, 1.0, bxdfs)
    }
}

impl Material for FourierMaterial {
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
