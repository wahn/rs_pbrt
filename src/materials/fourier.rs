//std
use std::sync::Arc;
// pbrt
use core::api::BsdfState;
use core::interaction::SurfaceInteraction;
use core::material::{Material, TransportMode};
use core::paramset::TextureParams;
use core::pbrt::Float;
use core::reflection::{Bsdf, Bxdf, FourierBSDFTable};
use core::texture::Texture;

// see fourier.h

pub struct FourierMaterial {
    pub bsdf_table: Arc<FourierBSDFTable>,
    pub bump_map: Option<Arc<Texture<Float> + Sync + Send>>,
}

impl FourierMaterial {
    pub fn new(
        bsdf_table: Arc<FourierBSDFTable>,
        bump_map: Option<Arc<Texture<Float> + Sync + Send>>
    ) -> Self {
        FourierMaterial {
            bump_map: bump_map,
            bsdf_table: bsdf_table,
        }
    }
    pub fn create(
        mp: &mut TextureParams,
        bsdf_state: &mut BsdfState,
    ) -> Arc<Material + Send + Sync> {
        let bump_map: Option<Arc<Texture<Float> + Send + Sync>> =
            mp.get_float_texture_or_null(String::from("bumpmap"));
        let bsdffile: String = mp.find_filename(String::from("bsdffile"), String::new());
        if let Some(bsdf_table) = bsdf_state.loaded_bsdfs.get(&bsdffile.clone()) {
            // use the BSDF table found
            Arc::new(FourierMaterial::new(bsdf_table.clone(), bump_map))
        } else {
            // read BSDF table from file
            let mut bsdf_table: FourierBSDFTable = FourierBSDFTable::default();
            println!(
                "reading {:?} returns {}",
                bsdffile,
                bsdf_table.read(&bsdffile)
            );
            Arc::new(FourierMaterial::new(Arc::new(bsdf_table), bump_map))
        }
    }
    pub fn bsdf(&self, si: &mut SurfaceInteraction) -> Bsdf {
        let mut bxdfs: Vec<Arc<Bxdf + Send + Sync>> = Vec::new();
        if let Some(ref bump_map) = self.bump_map {
            FourierMaterial::bump(bump_map, si);
        }
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
