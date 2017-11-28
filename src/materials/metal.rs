//std
use std::sync::Arc;
// pbrt
use core::interaction::SurfaceInteraction;
use core::material::TransportMode;
use core::microfacet::TrowbridgeReitzDistribution;
use core::pbrt::{Float, Spectrum};
use core::reflection::{Bsdf, Bxdf, FresnelConductor, MicrofacetReflection};
use materials::Material;
use textures::Texture;

pub struct MetalMaterial {
    pub eta: Arc<Texture<Spectrum> + Sync + Send>, // default: copper
    pub k: Arc<Texture<Spectrum> + Sync + Send>, // default: copper
    pub roughness: Arc<Texture<Float> + Sync + Send>, // default: 0.01
    pub u_roughness: Option<Arc<Texture<Float> + Sync + Send>>,
    pub v_roughness: Option<Arc<Texture<Float> + Sync + Send>>,
    // TODO: bump_map
    pub remap_roughness: bool,
}

impl MetalMaterial {
    pub fn new(eta: Arc<Texture<Spectrum> + Send + Sync>,
               k: Arc<Texture<Spectrum> + Send + Sync>,
               roughness: Arc<Texture<Float> + Sync + Send>,
               u_roughness: Option<Arc<Texture<Float> + Sync + Send>>,
               v_roughness: Option<Arc<Texture<Float> + Sync + Send>>,
               remap_roughness: bool)
               -> Self {
        MetalMaterial {
            eta: eta,
            k: k,
            roughness: roughness,
            u_roughness: u_roughness,
            v_roughness: v_roughness,
            remap_roughness: remap_roughness,
        }
    }
    pub fn bsdf(&self, si: &SurfaceInteraction) -> Bsdf {
        let mut bxdfs: Vec<Box<Bxdf + Send + Sync>> = Vec::new();
        let mut u_rough: Float;
        if let Some(ref u_roughness) = self.u_roughness {
            u_rough = u_roughness.evaluate(si);
        } else {
            u_rough = self.roughness.evaluate(si);
        }
        let mut v_rough: Float;
        if let Some(ref v_roughness) = self.v_roughness {
            v_rough = v_roughness.evaluate(si);
        } else {
            v_rough = self.roughness.evaluate(si);
        }
        if self.remap_roughness {
            u_rough = TrowbridgeReitzDistribution::roughness_to_alpha(u_rough);
            v_rough = TrowbridgeReitzDistribution::roughness_to_alpha(v_rough);
        }
        let fr_mf = Arc::new(FresnelConductor {
                                 eta_i: Spectrum::new(1.0 as Float),
                                 eta_t: self.eta.evaluate(si),
                                 k: self.k.evaluate(si),
                             });
        let distrib: Option<TrowbridgeReitzDistribution> =
            Some(TrowbridgeReitzDistribution::new(u_rough, v_rough, true));
        bxdfs
            .push(Box::new(MicrofacetReflection::new(Spectrum::new(1.0 as Float), distrib, fr_mf)));
        Bsdf::new(si, 1.0, bxdfs)
    }
}

impl Material for MetalMaterial {
    fn compute_scattering_functions(&self,
                                    si: &mut SurfaceInteraction,
                                    // arena: &mut Arena,
                                    _mode: TransportMode,
                                    _allow_multiple_lobes: bool) {
        si.bsdf = Some(Arc::new(self.bsdf(si)));
    }
}
