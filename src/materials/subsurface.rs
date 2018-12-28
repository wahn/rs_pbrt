//std
use std;
use std::sync::Arc;
// pbrt
use core::bssrdf::compute_beam_diffusion_bssrdf;
use core::bssrdf::BSSRDFTable;
use core::interaction::SurfaceInteraction;
use core::material::{Material, TransportMode};
use core::medium::get_medium_scattering_properties;
use core::paramset::TextureParams;
use core::pbrt::{Float, Spectrum};
use core::reflection::{Bsdf, Bxdf, LambertianReflection, OrenNayar};
use core::texture::Texture;

// see subsurface.h

pub struct SubsurfaceMaterial {
    pub scale: Float,                             // default: 1.0
    pub kr: Arc<Texture<Spectrum> + Sync + Send>, // default: 1.0
    pub kt: Arc<Texture<Spectrum> + Sync + Send>, // default: 1.0
    pub sigma_a: Arc<Texture<Spectrum> + Sync + Send>,
    pub sigma_s: Arc<Texture<Spectrum> + Sync + Send>,
    pub u_roughness: Arc<Texture<Float> + Sync + Send>, // default: 0.0
    pub v_roughness: Arc<Texture<Float> + Sync + Send>, // default: 0.0
    pub bump_map: Option<Arc<Texture<Float> + Send + Sync>>,
    pub eta: Float,            // default: 1.33
    pub remap_roughness: bool, // default: true
    pub table: BSSRDFTable,
}

impl SubsurfaceMaterial {
    pub fn new(
        scale: Float,
        kr: Arc<Texture<Spectrum> + Sync + Send>,
        kt: Arc<Texture<Spectrum> + Sync + Send>,
        sigma_a: Arc<Texture<Spectrum> + Sync + Send>,
        sigma_s: Arc<Texture<Spectrum> + Sync + Send>,
        g: Float,
        eta: Float,
        u_roughness: Arc<Texture<Float> + Sync + Send>,
        v_roughness: Arc<Texture<Float> + Sync + Send>,
        bump_map: Option<Arc<Texture<Float> + Sync + Send>>,
        remap_roughness: bool,
    ) -> Self {
        let mut table: BSSRDFTable = BSSRDFTable::new(100, 64);
        compute_beam_diffusion_bssrdf(g, eta, &mut table);
        SubsurfaceMaterial {
            scale: scale,
            kr: kr,
            kt: kt,
            sigma_a: sigma_a,
            sigma_s: sigma_s,
            u_roughness: u_roughness,
            v_roughness: v_roughness,
            bump_map: bump_map,
            eta: eta,
            remap_roughness: remap_roughness,
            table: table,
        }
    }
    pub fn create(mp: &mut TextureParams) -> Arc<Material + Send + Sync> {
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
        let sigma_a: Arc<Texture<Spectrum> + Sync + Send> =
            mp.get_spectrum_texture("sigma_a", sig_a);
        let sigma_s: Arc<Texture<Spectrum> + Sync + Send> =
            mp.get_spectrum_texture("sigma_s", sig_s);
        let kr: Arc<Texture<Spectrum> + Sync + Send> =
            mp.get_spectrum_texture("Kr", Spectrum::new(1.0));
        let kt: Arc<Texture<Spectrum> + Sync + Send> =
            mp.get_spectrum_texture("Kr", Spectrum::new(1.0));
        let roughu: Arc<Texture<Float> + Sync + Send> =
            mp.get_float_texture("uroughness", 0.0 as Float);
        let roughv: Arc<Texture<Float> + Sync + Send> =
            mp.get_float_texture("vroughness", 0.0 as Float);
        let bump_map = mp.get_float_texture_or_null("bumpmap");
        let remap_roughness: bool = mp.find_bool("remaproughness", true);
        Arc::new(SubsurfaceMaterial::new(
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
        ))
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
        if let Some(ref bump_map) = self.bump_map {
            Self::bump(bump_map, si);
        }
        let mut bxdfs: Vec<Arc<Bxdf + Send + Sync>> = Vec::new();
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
        // si->bsdf = ARENA_ALLOC(arena, BSDF)(*si, eta);

        // if (R.IsBlack() && T.IsBlack()) return;

        // bool isSpecular = urough == 0 && vrough == 0;
        // if (isSpecular && allowMultipleLobes) {
        //     si->bsdf->Add(
        //         ARENA_ALLOC(arena, FresnelSpecular)(R, T, 1.f, eta, mode));
        // } else {
        //     if (remapRoughness) {
        //         urough = TrowbridgeReitzDistribution::RoughnessToAlpha(urough);
        //         vrough = TrowbridgeReitzDistribution::RoughnessToAlpha(vrough);
        //     }
        //     MicrofacetDistribution *distrib =
        //         isSpecular ? nullptr
        //                    : ARENA_ALLOC(arena, TrowbridgeReitzDistribution)(
        //                          urough, vrough);
        //     if (!R.IsBlack()) {
        //         Fresnel *fresnel = ARENA_ALLOC(arena, FresnelDielectric)(1.f, eta);
        //         if (isSpecular)
        //             si->bsdf->Add(
        //                 ARENA_ALLOC(arena, SpecularReflection)(R, fresnel));
        //         else
        //             si->bsdf->Add(ARENA_ALLOC(arena, MicrofacetReflection)(
        //                 R, distrib, fresnel));
        //     }
        //     if (!T.IsBlack()) {
        //         if (isSpecular)
        //             si->bsdf->Add(ARENA_ALLOC(arena, SpecularTransmission)(
        //                 T, 1.f, eta, mode));
        //         else
        //             si->bsdf->Add(ARENA_ALLOC(arena, MicrofacetTransmission)(
        //                 T, distrib, 1.f, eta, mode));
        //     }
        // }
        // Spectrum sig_a = scale * sigma_a->Evaluate(*si).Clamp();
        // Spectrum sig_s = scale * sigma_s->Evaluate(*si).Clamp();
        // si->bssrdf = ARENA_ALLOC(arena, TabulatedBSSRDF)(*si, this, mode, eta,
        //                                                  sig_a, sig_s, table);
        si.bsdf = Some(Arc::new(Bsdf::new(si, 1.0, bxdfs)));
    }
}
