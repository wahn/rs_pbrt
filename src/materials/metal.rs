//std
use std::sync::Arc;
// pbrt
use crate::core::interaction::SurfaceInteraction;
use crate::core::material::{Material, TransportMode};
use crate::core::microfacet::{MicrofacetDistribution, TrowbridgeReitzDistribution};
use crate::core::paramset::TextureParams;
use crate::core::pbrt::{Float, Spectrum};
use crate::core::reflection::{Bsdf, Bxdf, Fresnel, FresnelConductor, MicrofacetReflection};
use crate::core::texture::Texture;

pub const COPPER_SAMPLES: u8 = 56_u8;
pub const COPPER_WAVELENGTHS: [Float; COPPER_SAMPLES as usize] = [
    298.757_055_4,
    302.400_434_1,
    306.133_772_8,
    309.960_445,
    313.883_994_9,
    317.908_148_7,
    322.036_826,
    326.274_152_6,
    330.624_474_7,
    335.092_373,
    339.682_679_5,
    344.400_494_4,
    349.251_205_6,
    354.240_508_6,
    359.374_429,
    364.659_347_1,
    370.102_023_9,
    375.709_630_3,
    381.489_778_5,
    387.450_556_3,
    393.600_565_1,
    399.948_961_3,
    406.505_501_6,
    413.280_593_3,
    420.285_349_2,
    427.531_648_3,
    435.032_203_5,
    442.800_635_7,
    450.851_556_4,
    459.200_659_3,
    467.864_822_6,
    476.862_223_1,
    486.212_462_7,
    495.936_712,
    506.057_869_4,
    516.600_741_7,
    527.592_246_8,
    539.061_643_5,
    551.040_791_1,
    563.564_445_5,
    576.670_595_3,
    590.400_847_6,
    604.800_868_3,
    619.920_89,
    635.816_297_4,
    652.548_305_3,
    670.184_745_9,
    688.800_988_9,
    708.481_017_1,
    729.318_694_1,
    751.419_260_6,
    774.901_112_5,
    799.897_922_6,
    826.561_186_7,
    855.063_296_6,
    885.601_271_4,
];
pub const COPPER_N: [Float; COPPER_SAMPLES as usize] = [
    1.400_313, 1.38, 1.358_438, 1.34, 1.329_063, 1.325, 1.3325, 1.34, 1.334_375, 1.325, 1.317_812,
    1.31, 1.300_313, 1.29, 1.281_563, 1.27, 1.249_062, 1.225, 1.2, 1.18, 1.174_375, 1.175, 1.1775,
    1.18, 1.178_125, 1.175, 1.172_812, 1.17, 1.165_312, 1.16, 1.155_312, 1.15, 1.142_812, 1.135,
    1.131_562, 1.12, 1.092_437, 1.04, 0.950_375, 0.826, 0.645_875, 0.468, 0.35125, 0.272,
    0.230_813, 0.214, 0.20925, 0.213, 0.21625, 0.223, 0.2365, 0.25, 0.254_188, 0.26, 0.28, 0.3,
];
pub const COPPER_K: [Float; COPPER_SAMPLES as usize] = [
    1.662_125, 1.687, 1.703_313, 1.72, 1.744_563, 1.77, 1.791_625, 1.81, 1.822_125, 1.834, 1.85175,
    1.872, 1.89425, 1.916, 1.931_688, 1.95, 1.972_438, 2.015, 2.121_562, 2.21, 2.177_188, 2.13,
    2.160_063, 2.21, 2.249_938, 2.289, 2.326, 2.362, 2.397_625, 2.433, 2.469_187, 2.504, 2.535_875,
    2.564, 2.589_625, 2.605, 2.595_562, 2.583, 2.5765, 2.599, 2.678_062, 2.809, 3.01075, 3.24,
    3.458_187, 3.67, 3.863_125, 4.05, 4.239_563, 4.43, 4.619_563, 4.817, 5.034_125, 5.26,
    5.485_625, 5.717,
];

pub struct MetalMaterial {
    pub eta: Arc<dyn Texture<Spectrum> + Sync + Send>, // default: copper
    pub k: Arc<dyn Texture<Spectrum> + Sync + Send>,   // default: copper
    pub roughness: Arc<dyn Texture<Float> + Sync + Send>, // default: 0.01
    pub u_roughness: Option<Arc<dyn Texture<Float> + Sync + Send>>,
    pub v_roughness: Option<Arc<dyn Texture<Float> + Sync + Send>>,
    pub bump_map: Option<Arc<dyn Texture<Float> + Send + Sync>>,
    pub remap_roughness: bool,
}

impl MetalMaterial {
    pub fn new(
        eta: Arc<dyn Texture<Spectrum> + Send + Sync>,
        k: Arc<dyn Texture<Spectrum> + Send + Sync>,
        roughness: Arc<dyn Texture<Float> + Sync + Send>,
        u_roughness: Option<Arc<dyn Texture<Float> + Sync + Send>>,
        v_roughness: Option<Arc<dyn Texture<Float> + Sync + Send>>,
        bump_map: Option<Arc<dyn Texture<Float> + Sync + Send>>,
        remap_roughness: bool,
    ) -> Self {
        MetalMaterial {
            eta,
            k,
            roughness,
            u_roughness,
            v_roughness,
            bump_map,
            remap_roughness,
        }
    }
    pub fn create(mp: &mut TextureParams) -> Arc<Material> {
        let copper_n: Spectrum =
            Spectrum::from_sampled(&COPPER_WAVELENGTHS, &COPPER_N, COPPER_SAMPLES as i32);
        let eta: Arc<dyn Texture<Spectrum> + Send + Sync> =
            mp.get_spectrum_texture("eta", copper_n);
        let copper_k: Spectrum =
            Spectrum::from_sampled(&COPPER_WAVELENGTHS, &COPPER_K, COPPER_SAMPLES as i32);
        let k: Arc<dyn Texture<Spectrum> + Send + Sync> = mp.get_spectrum_texture("k", copper_k);
        let roughness: Arc<dyn Texture<Float> + Send + Sync> =
            mp.get_float_texture("roughness", 0.01 as Float);
        let u_roughness: Option<Arc<dyn Texture<Float> + Send + Sync>> =
            mp.get_float_texture_or_null("uroughness");
        let v_roughness: Option<Arc<dyn Texture<Float> + Send + Sync>> =
            mp.get_float_texture_or_null("vroughness");
        let bump_map = mp.get_float_texture_or_null("bumpmap");
        let remap_roughness: bool = mp.find_bool("remaproughness", true);
        Arc::new(Material::Metal(Box::new(MetalMaterial::new(
            eta,
            k,
            roughness,
            u_roughness,
            v_roughness,
            bump_map,
            remap_roughness,
        ))))
    }
    // Material
    pub fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        // arena: &mut Arena,
        _mode: TransportMode,
        _allow_multiple_lobes: bool,
        _material: Option<Arc<Material>>,
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
        let fr_mf = Fresnel::Conductor(FresnelConductor {
            eta_i: Spectrum::new(1.0 as Float),
            eta_t: self.eta.evaluate(si),
            k: self.k.evaluate(si),
        });
        let distrib = MicrofacetDistribution::TrowbridgeReitz(TrowbridgeReitzDistribution::new(
            u_rough, v_rough, true,
        ));
        si.bsdf = Some(Bsdf::new(si, 1.0));
        if let Some(bsdf) = &mut si.bsdf {
            let bxdf_idx: usize = 0;
            if use_scale {
                bsdf.bxdfs[bxdf_idx] = Bxdf::MicrofacetRefl(MicrofacetReflection::new(
                    Spectrum::new(1.0 as Float),
                    distrib,
                    fr_mf,
                    Some(sc),
                ));
            } else {
                bsdf.bxdfs[bxdf_idx] = Bxdf::MicrofacetRefl(MicrofacetReflection::new(
                    Spectrum::new(1.0 as Float),
                    distrib,
                    fr_mf,
                    None,
                ));
            }
        }
    }
}
