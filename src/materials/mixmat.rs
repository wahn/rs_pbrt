//std
use std;
use std::sync::Arc;
// pbrt
use crate::core::interaction::SurfaceInteraction;
use crate::core::material::{Material, TransportMode};
use crate::core::pbrt::{Float, Spectrum};
use crate::core::reflection::{Bsdf, Bxdf, BxdfType, ScaledBxDF};
use crate::core::texture::Texture;

// see mixmat.h

/// The mix material takes two other materials and a texture and uses
/// the value returned by the texture to blend between the two
/// materials at the point being shaded.
pub struct MixMaterial {
    pub m1: Arc<Material + Sync + Send>,
    pub m2: Arc<Material + Sync + Send>,
    pub scale: Arc<Texture<Spectrum> + Sync + Send>, // default: 0.5
}

impl MixMaterial {
    pub fn new(
        m1: Arc<Material + Sync + Send>,
        m2: Arc<Material + Sync + Send>,
        scale: Arc<Texture<Spectrum> + Send + Sync>,
    ) -> Self {
        MixMaterial {
            m1: m1,
            m2: m2,
            scale: scale,
        }
    }
}

impl Material for MixMaterial {
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        // arena: &mut Arena,
        mode: TransportMode,
        allow_multiple_lobes: bool,
        _material: Option<Arc<Material + Send + Sync>>,
    ) {
        let mut bxdfs: Vec<Arc<Bxdf + Send + Sync>> = Vec::new();
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
        self.m1
            .compute_scattering_functions(si, mode.clone(), allow_multiple_lobes, None);
        self.m2
            .compute_scattering_functions(&mut si2, mode.clone(), allow_multiple_lobes, None);
        let bsdf_flags: u8 = BxdfType::BsdfAll as u8;
        if let Some(ref bsdf1) = si.bsdf {
            let n1: u8 = bsdf1.num_components(bsdf_flags);
            if let Some(ref bsdf2) = si2.bsdf {
                let n2: u8 = bsdf2.num_components(bsdf_flags);
                for i in 0..n1 {
                    let bxdf = &bsdf1.bxdfs[i as usize];
                    bxdfs.push(Arc::new(ScaledBxDF::new(bxdf.clone(), s1)));
                }
                for i in 0..n2 {
                    let bxdf = &bsdf2.bxdfs[i as usize];
                    bxdfs.push(Arc::new(ScaledBxDF::new(bxdf.clone(), s2)));
                }
            }
        }
        si.bsdf = Some(Arc::new(Bsdf::new(si, 1.0, bxdfs)));
    }
}
