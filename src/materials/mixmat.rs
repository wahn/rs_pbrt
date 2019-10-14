//std
use std;
use std::sync::Arc;
// pbrt
use crate::core::interaction::SurfaceInteraction;
use crate::core::material::{Material, TransportMode};
use crate::core::pbrt::{Float, Spectrum};
use crate::core::reflection::Bxdf;
use crate::core::texture::Texture;

// see mixmat.h

/// The mix material takes two other materials and a texture and uses
/// the value returned by the texture to blend between the two
/// materials at the point being shaded.
pub struct MixMaterial {
    pub m1: Arc<dyn Material + Sync + Send>,
    pub m2: Arc<dyn Material + Sync + Send>,
    pub scale: Arc<dyn Texture<Spectrum> + Sync + Send>, // default: 0.5
}

impl MixMaterial {
    pub fn new(
        m1: Arc<dyn Material + Sync + Send>,
        m2: Arc<dyn Material + Sync + Send>,
        scale: Arc<dyn Texture<Spectrum> + Send + Sync>,
    ) -> Self {
        MixMaterial { m1, m2, scale }
    }
}

impl Material for MixMaterial {
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        // arena: &mut Arena,
        mode: TransportMode,
        allow_multiple_lobes: bool,
        _material: Option<Arc<dyn Material + Send + Sync>>,
        _scale: Option<Spectrum>,
    ) {
        let s1: Spectrum = self
            .scale
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        let s2: Spectrum =
            (Spectrum::new(1.0 as Float) - s1).clamp(0.0 as Float, std::f32::INFINITY as Float);
        let mut si2: SurfaceInteraction = SurfaceInteraction::new(
            &si.p,
            &si.p_error,
            &si.uv,
            &si.wo,
            &si.dpdu,
            &si.dpdv,
            &si.dndu,
            &si.dndv,
            si.time,
            si.shape,
        );
        self.m1.compute_scattering_functions(
            si,
            mode.clone(),
            allow_multiple_lobes,
            None,
            Some(s1),
        );
        self.m2.compute_scattering_functions(
            &mut si2,
            mode.clone(),
            allow_multiple_lobes,
            None,
            Some(s2),
        );
        let mut last_idx: usize = 0;
        if let Some(bsdf) = &mut si.bsdf {
            for bxdf_idx in 0..8 {
                bsdf.bxdfs[bxdf_idx] = match bsdf.bxdfs[bxdf_idx] {
                    Bxdf::Empty(_bxdf) => {
                        last_idx = bxdf_idx;
                        break;
                    },
                    _ => bsdf.bxdfs[bxdf_idx],
                };
            }
            if let Some(bsdf2) = si2.bsdf {
                for bxdf_idx in 0..8 {
                    bsdf.bxdfs[bxdf_idx + last_idx] = match bsdf2.bxdfs[bxdf_idx] {
                        Bxdf::Empty(_bxdf) => {
                            break;
                        },
                        _ => bsdf2.bxdfs[bxdf_idx],
                    };
                }
            }
        }
    }
}
