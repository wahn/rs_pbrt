//std
use std;
use std::sync::Arc;
// others
use time::PreciseTime;
// pbrt
use crate::core::bssrdf::compute_beam_diffusion_bssrdf;
use crate::core::bssrdf::BssrdfTable;
use crate::core::bssrdf::TabulatedBssrdf;
use crate::core::interaction::SurfaceInteraction;
use crate::core::material::{Material, TransportMode};
use crate::core::medium::get_medium_scattering_properties;
use crate::core::microfacet::{MicrofacetDistribution, TrowbridgeReitzDistribution};
use crate::core::paramset::TextureParams;
use crate::core::pbrt::{Float, Spectrum};
use crate::core::reflection::{
    Bsdf, Bxdf, Fresnel, FresnelDielectric, FresnelSpecular, MicrofacetReflection,
    MicrofacetTransmission, SpecularReflection, SpecularTransmission,
};
use crate::core::texture::Texture;

// see subsurface.h

pub struct SubsurfaceMaterial {
    pub scale: Float,                                 // default: 1.0
    pub kr: Arc<dyn Texture<Spectrum> + Sync + Send>, // default: 1.0
    pub kt: Arc<dyn Texture<Spectrum> + Sync + Send>, // default: 1.0
    pub sigma_a: Arc<dyn Texture<Spectrum> + Sync + Send>,
    pub sigma_s: Arc<dyn Texture<Spectrum> + Sync + Send>,
    pub u_roughness: Arc<dyn Texture<Float> + Sync + Send>, // default: 0.0
    pub v_roughness: Arc<dyn Texture<Float> + Sync + Send>, // default: 0.0
    pub bump_map: Option<Arc<dyn Texture<Float> + Send + Sync>>,
    pub eta: Float,            // default: 1.33
    pub remap_roughness: bool, // default: true
    pub table: Arc<BssrdfTable>,
}

impl SubsurfaceMaterial {
    pub fn new(
        scale: Float,
        kr: Arc<dyn Texture<Spectrum> + Sync + Send>,
        kt: Arc<dyn Texture<Spectrum> + Sync + Send>,
        sigma_a: Arc<dyn Texture<Spectrum> + Sync + Send>,
        sigma_s: Arc<dyn Texture<Spectrum> + Sync + Send>,
        g: Float,
        eta: Float,
        u_roughness: Arc<dyn Texture<Float> + Sync + Send>,
        v_roughness: Arc<dyn Texture<Float> + Sync + Send>,
        bump_map: Option<Arc<dyn Texture<Float> + Sync + Send>>,
        remap_roughness: bool,
    ) -> Self {
        let mut table: BssrdfTable = BssrdfTable::new(100, 64);
        compute_beam_diffusion_bssrdf(g, eta, &mut table);
        SubsurfaceMaterial {
            scale,
            kr,
            kt,
            sigma_a,
            sigma_s,
            u_roughness,
            v_roughness,
            bump_map,
            eta,
            remap_roughness,
            table: Arc::new(table),
        }
    }
    pub fn create(mp: &mut TextureParams) -> Arc<Material> {
        let sig_a_rgb: [Float; 3] = [0.0011, 0.0024, 0.014];
        let sig_s_rgb: [Float; 3] = [2.55, 3.21, 3.77];
        let mut sig_a: Spectrum = Spectrum::from_rgb(&sig_a_rgb);
        let mut sig_s: Spectrum = Spectrum::from_rgb(&sig_s_rgb);
        let name: String = mp.find_string("name", String::from(""));
        let found: bool = get_medium_scattering_properties(&name, &mut sig_a, &mut sig_s);
        let mut g: Float = mp.find_float("g", 0.0 as Float);
        if name != String::from("") {
            if !found {
                println!(
                    "WARNING: Named material {:?} not found.  Using defaults.",
                    name
                );
            } else {
                // enforce g=0 (the database specifies reduced
                // scattering coefficients)
                g = 0.0;
            }
        }
        let scale: Float = mp.find_float("scale", 1.0 as Float);
        let eta: Float = mp.find_float("eta", 1.33 as Float);
        let sigma_a: Arc<dyn Texture<Spectrum> + Sync + Send> =
            mp.get_spectrum_texture("sigma_a", sig_a);
        let sigma_s: Arc<dyn Texture<Spectrum> + Sync + Send> =
            mp.get_spectrum_texture("sigma_s", sig_s);
        let kr: Arc<dyn Texture<Spectrum> + Sync + Send> =
            mp.get_spectrum_texture("Kr", Spectrum::new(1.0));
        let kt: Arc<dyn Texture<Spectrum> + Sync + Send> =
            mp.get_spectrum_texture("Kr", Spectrum::new(1.0));
        let roughu: Arc<dyn Texture<Float> + Sync + Send> =
            mp.get_float_texture("uroughness", 0.0 as Float);
        let roughv: Arc<dyn Texture<Float> + Sync + Send> =
            mp.get_float_texture("vroughness", 0.0 as Float);
        let bump_map = mp.get_float_texture_or_null("bumpmap");
        let remap_roughness: bool = mp.find_bool("remaproughness", true);
        let start = PreciseTime::now();
        let tmp = Arc::new(Material::Subsurface(SubsurfaceMaterial::new(
            scale,
            kr,
            kt,
            sigma_a,
            sigma_s,
            g,
            eta,
            roughu,
            roughv,
            bump_map,
            remap_roughness,
        )));
        let end = PreciseTime::now();
        println!(
            "{} seconds for SubsurfaceMaterial::new() ...",
            start.to(end)
        );
        tmp
    }
    // Material
    pub fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        // arena: &mut Arena,
        mode: TransportMode,
        allow_multiple_lobes: bool,
        material: Option<Arc<Material>>,
        scale_opt: Option<Spectrum>,
    ) {
        let mut use_scale: bool = false;
        let mut sc: Spectrum = Spectrum::default();
        if let Some(scale) = scale_opt {
            use_scale = true;
            sc = scale;
        }
        if let Some(ref bump) = self.bump_map {
            Material::bump(bump, si);
        }
        let mut bxdf_idx: usize = 0;
        // initialize BSDF for _SubsurfaceMaterial_
        let r: Spectrum = self
            .kr
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        let t: Spectrum = self
            .kt
            .evaluate(si)
            .clamp(0.0 as Float, std::f32::INFINITY as Float);
        let mut urough: Float = self.u_roughness.evaluate(si);
        let mut vrough: Float = self.v_roughness.evaluate(si);
        // initialize _bsdf_ for smooth or rough dielectric
        if r.is_black() && t.is_black() {
            return;
        }
        let is_specular: bool = urough == 0.0 as Float && vrough == 0.0 as Float;
        si.bsdf = Some(Bsdf::new(si, self.eta));
        if let Some(bsdf) = &mut si.bsdf {
            if is_specular && allow_multiple_lobes {
                if use_scale {
                    bsdf.bxdfs[bxdf_idx] = Bxdf::FresnelSpec(FresnelSpecular::new(
                        r,
                        t,
                        1.0 as Float,
                        self.eta,
                        mode,
                        Some(sc),
                    ));
                } else {
                    bsdf.bxdfs[bxdf_idx] = Bxdf::FresnelSpec(FresnelSpecular::new(
                        r,
                        t,
                        1.0 as Float,
                        self.eta,
                        mode,
                        None,
                    ));
                }
            } else {
                if self.remap_roughness {
                    urough = TrowbridgeReitzDistribution::roughness_to_alpha(urough);
                    vrough = TrowbridgeReitzDistribution::roughness_to_alpha(vrough);
                }
                if !r.is_black() {
                    let fresnel = Fresnel::Dielectric(FresnelDielectric {
                        eta_i: 1.0 as Float,
                        eta_t: self.eta,
                    });
                    if is_specular {
                        if use_scale {
                            bsdf.bxdfs[bxdf_idx] =
                                Bxdf::SpecRefl(SpecularReflection::new(r, fresnel, Some(sc)));
                            bxdf_idx += 1;
                        } else {
                            bsdf.bxdfs[bxdf_idx] =
                                Bxdf::SpecRefl(SpecularReflection::new(r, fresnel, None));
                            bxdf_idx += 1;
                        }
                    } else {
                        let distrib = MicrofacetDistribution::TrowbridgeReitz(
                            TrowbridgeReitzDistribution::new(urough, vrough, true),
                        );
                        if use_scale {
                            bsdf.bxdfs[bxdf_idx] = Bxdf::MicrofacetRefl(MicrofacetReflection::new(
                                r,
                                distrib,
                                fresnel,
                                Some(sc),
                            ));
                            bxdf_idx += 1;
                        } else {
                            bsdf.bxdfs[bxdf_idx] = Bxdf::MicrofacetRefl(MicrofacetReflection::new(
                                r, distrib, fresnel, None,
                            ));
                            bxdf_idx += 1;
                        }
                    }
                }
                if !t.is_black() {
                    if is_specular {
                        if use_scale {
                            bsdf.bxdfs[bxdf_idx] = Bxdf::SpecTrans(SpecularTransmission::new(
                                t,
                                1.0,
                                self.eta,
                                mode,
                                Some(sc),
                            ));
                        } else {
                            bsdf.bxdfs[bxdf_idx] = Bxdf::SpecTrans(SpecularTransmission::new(
                                t, 1.0, self.eta, mode, None,
                            ));
                        }
                    } else {
                        let distrib = MicrofacetDistribution::TrowbridgeReitz(
                            TrowbridgeReitzDistribution::new(urough, vrough, true),
                        );
                        if use_scale {
                            bsdf.bxdfs[bxdf_idx] =
                                Bxdf::MicrofacetTrans(MicrofacetTransmission::new(
                                    t,
                                    distrib,
                                    1.0,
                                    self.eta,
                                    mode,
                                    Some(sc),
                                ));
                        } else {
                            bsdf.bxdfs[bxdf_idx] = Bxdf::MicrofacetTrans(
                                MicrofacetTransmission::new(t, distrib, 1.0, self.eta, mode, None),
                            );
                        }
                    }
                }
            }
            let sig_a: Spectrum = self.scale
                * self
                    .sigma_a
                    .evaluate(si)
                    .clamp(0.0 as Float, std::f32::INFINITY as Float);
            let sig_s: Spectrum = self.scale
                * self
                    .sigma_s
                    .evaluate(si)
                    .clamp(0.0 as Float, std::f32::INFINITY as Float);
            si.bssrdf = Some(TabulatedBssrdf::new(
                si,
                material,
                mode,
                self.eta,
                &sig_a,
                &sig_s,
                self.table.clone(),
            ));
        }
    }
}
