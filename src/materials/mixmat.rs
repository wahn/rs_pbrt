//std
use std;
use std::sync::Arc;
// pbrt
use crate::core::interaction::SurfaceInteraction;
use crate::core::material::{Material, TransportMode};
use crate::core::pbrt::{Float, Spectrum};
use crate::core::reflection::{
    Bxdf, Fresnel, FresnelConductor, FresnelDielectric, FresnelNoOp, NoBxdf, SpecularReflection,
};
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
        // find next empty slot
        if let Some(bsdf) = &si.bsdf {
            for bxdf_idx in 0..8 {
                match &bsdf.bxdfs[bxdf_idx] {
                    Bxdf::Empty(_bxdf) => {
                        last_idx = bxdf_idx;
                        break;
                    }
                    _ => {}
                }
            }
        }
        // get Bxdfs from si2 before it gets out of scope
        if si2.bsdf.is_some() {
            let bxdfs: [Bxdf; 8] = si2.bsdf.unwrap().bxdfs;
            if let Some(bsdf) = &mut si.bsdf {
                for bxdf_idx in 0..8 {
                    bsdf.bxdfs[bxdf_idx + last_idx] = match &bxdfs[bxdf_idx] {
                        Bxdf::Empty(_bxdf) => break,
                        Bxdf::SpecRefl(bxdf) =>  {
                            let fresnel = match &bxdf.fresnel {
                                Fresnel::Conductor(fresnel) => Fresnel::Conductor(FresnelConductor {
                                    eta_i: fresnel.eta_i,
                                    eta_t: fresnel.eta_t,
                                    k: fresnel.k,
                                }),
                                Fresnel::Dielectric(fresnel) => Fresnel::Dielectric(FresnelDielectric {
                                    eta_i: fresnel.eta_i,
                                    eta_t: fresnel.eta_t,
                                }),
                                _ => Fresnel::NoOp(FresnelNoOp {})
                            };
                            Bxdf::SpecRefl(SpecularReflection::new(
                            bxdf.r,
                            fresnel,
                            bxdf.sc_opt,
                        ))}
                        ,
                        // Bxdf::SpecTrans(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
                        // Bxdf::FresnelSpec(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
                        // Bxdf::LambertianRefl(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
                        // Bxdf::LambertianTrans(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
                        // Bxdf::OrenNayarRefl(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
                        // Bxdf::MicrofacetRefl(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
                        // Bxdf::MicrofacetTrans(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
                        // Bxdf::FresnelBlnd(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
                        // Bxdf::Fourier(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
                        // // Bxdf::Bssrdf(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
                        // Bxdf::DisDiff(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
                        // Bxdf::DisSS(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
                        // Bxdf::DisRetro(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
                        // Bxdf::DisSheen(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
                        // Bxdf::DisClearCoat(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
                        // Bxdf::Hair(bxdf) => bxdf.get_type() & t == bxdf.get_type(),
                        _ => Bxdf::Empty(NoBxdf::default()),
                    };
                }
            }
        }
    }
}
