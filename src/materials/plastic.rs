//std
use std;
use std::sync::Arc;
// pbrt
use core::interaction::SurfaceInteraction;
use core::material::{Material, TransportMode};
use core::microfacet::TrowbridgeReitzDistribution;
use core::pbrt::{Float, Spectrum};
use core::reflection::{Bsdf, Bxdf, FresnelDielectric, LambertianReflection, MicrofacetReflection};
use textures::Texture;

// see plastic.h

/// Plastic can be modeled as a mixture of a diffuse and glossy
/// scattering function.
pub struct PlasticMaterial {
    pub kd: Arc<Texture<Spectrum> + Sync + Send>, // default: 0.25
    pub ks: Arc<Texture<Spectrum> + Sync + Send>, // default: 0.25
    pub roughness: Arc<Texture<Float> + Sync + Send>, // default: 0.1
    // TODO: bump_map
    pub remap_roughness: bool,
}

impl PlasticMaterial {
    pub fn new(kd: Arc<Texture<Spectrum> + Send + Sync>,
               ks: Arc<Texture<Spectrum> + Send + Sync>,
               roughness: Arc<Texture<Float> + Sync + Send>,
               remap_roughness: bool)
               -> Self {
        PlasticMaterial {
            kd: kd,
            ks: ks,
            roughness: roughness,
            remap_roughness: remap_roughness,
        }
    }
    pub fn bsdf(&self, si: &SurfaceInteraction) -> Bsdf {
        let mut bxdfs: Vec<Box<Bxdf + Send + Sync>> = Vec::new();
        // initialize diffuse component of plastic material
        let kd: Spectrum = self.kd
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        if !kd.is_black() {
            bxdfs.push(Box::new(LambertianReflection::new(kd)));
        }
        // initialize specular component of plastic material
        let ks: Spectrum = self.ks
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        if !ks.is_black() {
            let fresnel = Arc::new(FresnelDielectric {
                                       eta_i: 1.5 as Float,
                                       eta_t: 1.0 as Float,
                                   });
            // create microfacet distribution _distrib_ for plastic material
            let mut rough: Float = self.roughness.evaluate(si);
            if self.remap_roughness {
                rough = TrowbridgeReitzDistribution::roughness_to_alpha(rough);
            }
            let distrib: TrowbridgeReitzDistribution = TrowbridgeReitzDistribution {
                alpha_x: rough,
                alpha_y: rough,
                sample_visible_area: true,
            };
            bxdfs.push(Box::new(MicrofacetReflection::new(ks, Some(distrib), fresnel)));
        }
        Bsdf::new(si, 1.0, bxdfs)
    }
}

impl Material for PlasticMaterial {
    fn compute_scattering_functions(&self,
                                    si: &mut SurfaceInteraction,
                                    // arena: &mut Arena,
                                    _mode: TransportMode,
                                    _allow_multiple_lobes: bool) {
        si.bsdf = Some(Arc::new(self.bsdf(si)));
    }
}
