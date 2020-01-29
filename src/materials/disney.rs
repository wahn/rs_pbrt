use std::f32;
use std::sync::Arc;

use num::Zero;

use crate::core::geometry::{spherical_direction, vec3_abs_dot_vec3, vec3_dot_vec3};
use crate::core::geometry::{Point2f, Vector3f};
use crate::core::interaction::SurfaceInteraction;
use crate::core::material::{Material, TransportMode};
use crate::core::microfacet::{MicrofacetDistribution, TrowbridgeReitzDistribution};
use crate::core::paramset::TextureParams;
use crate::core::pbrt::{clamp_t, lerp};
use crate::core::pbrt::{Float, Spectrum};
use crate::core::reflection::reflect;
use crate::core::reflection::{abs_cos_theta, fr_schlick, vec3_same_hemisphere_vec3};
use crate::core::reflection::{
    Bsdf, Bxdf, BxdfType, DisneyFresnel, Fresnel, LambertianTransmission, MicrofacetReflection,
    MicrofacetTransmission, SpecularTransmission,
};
use crate::core::texture::Texture;

pub struct DisneyMaterial {
    color: Arc<dyn Texture<Spectrum> + Send + Sync>,
    // base_color: Arc<TextureFloat>,
    metallic: Arc<dyn Texture<Float> + Send + Sync>,
    eta: Arc<dyn Texture<Float> + Send + Sync>,
    roughness: Arc<dyn Texture<Float> + Send + Sync>,
    specular_tint: Arc<dyn Texture<Float> + Send + Sync>,
    anisotropic: Arc<dyn Texture<Float> + Send + Sync>,
    sheen: Arc<dyn Texture<Float> + Send + Sync>,
    sheen_tint: Arc<dyn Texture<Float> + Send + Sync>,
    clearcoat: Arc<dyn Texture<Float> + Send + Sync>,
    clearcoat_gloss: Arc<dyn Texture<Float> + Send + Sync>,
    spec_trans: Arc<dyn Texture<Float> + Send + Sync>,
    scatter_distance: Arc<dyn Texture<Spectrum> + Send + Sync>,
    flatness: Arc<dyn Texture<Float> + Send + Sync>,
    diff_trans: Arc<dyn Texture<Float> + Send + Sync>,
    bump_map: Option<Arc<dyn Texture<Float> + Send + Sync>>,
    thin: bool,
}

impl DisneyMaterial {
    pub fn create(mp: &mut TextureParams) -> Arc<Material> {
        let color = mp.get_spectrum_texture("color", Spectrum::from(0.5));
        let metallic = mp.get_float_texture("metallic", 0.0);
        let eta = mp.get_float_texture("eta", 1.5);
        let roughness = mp.get_float_texture("roughness", 0.5);
        let specular_tint = mp.get_float_texture("speculartint", 0.0);
        let anisotropic = mp.get_float_texture("anisotropic", 0.0);
        let sheen = mp.get_float_texture("sheen", 0.0);
        let sheen_tint = mp.get_float_texture("sheentint", 0.5);
        let clearcoat = mp.get_float_texture("clearcoat", 0.0);
        let clearcoat_gloss = mp.get_float_texture("clearcoatgloss", 1.0);
        let spec_trans = mp.get_float_texture("spectrans", 0.0);
        let scatter_distance = mp.get_spectrum_texture("scatterdistance", Spectrum::from(0.0));
        let thin = mp.find_bool("thin", false);
        let flatness = mp.get_float_texture("flatness", 0.0);
        let diff_trans = mp.get_float_texture("difftrans", 1.0);
        let bump_map = mp.get_float_texture_or_null("bumpmap");

        Arc::new(Material::Disney(DisneyMaterial {
            color,
            metallic,
            eta,
            roughness,
            specular_tint,
            anisotropic,
            sheen,
            sheen_tint,
            clearcoat,
            clearcoat_gloss,
            spec_trans,
            scatter_distance,
            flatness,
            diff_trans,
            bump_map,
            thin,
        }))
    }
    // Material
    pub fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        mode: TransportMode,
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
        // diffuse
        let c = self.color.evaluate(si).clamp(0.0, f32::INFINITY);
        let metallic_weight = self.metallic.evaluate(si);
        let e = self.eta.evaluate(si);
        let strans = self.spec_trans.evaluate(si);
        let diffuse_weight = (1.0 - metallic_weight) * (1.0 - strans);
        let dt = self.diff_trans.evaluate(si) / 2.0; // 0: all diffuse is reflected -> 1, transmitted
        let rough = self.roughness.evaluate(si);
        let lum = c.y();
        // normalize lum. to isolate hue+sat
        let c_tint = if lum > 0.0 {
            c / lum
        } else {
            Spectrum::new(1.0)
        };
        let sheen_weight = self.sheen.evaluate(si);
        let c_sheen = if sheen_weight > 0.0 {
            let stint = self.sheen_tint.evaluate(si);
            lerp(stint, Spectrum::new(1.0), c_tint)
        } else {
            Spectrum::zero()
        };
        let flat = self.flatness.evaluate(si);
        let sd = self.scatter_distance.evaluate(si);
        let aspect = Float::sqrt(1.0 - self.anisotropic.evaluate(si) * 0.9);
        let spec_tint = self.specular_tint.evaluate(si);
        let cc = self.clearcoat.evaluate(si);
        let gloss: Float = lerp(self.clearcoat_gloss.evaluate(si), 0.1, 0.001);
        si.bsdf = Some(Bsdf::new(si, 1.0));
        if let Some(bsdf) = &mut si.bsdf {
            let mut bxdf_idx: usize = 0;
            if diffuse_weight > 0.0 {
                if self.thin {
                    // Blend between DisneyDiffuse and fake subsurface based on flatness. Additionally,
                    // weight using diff_trans.
                    if use_scale {
                        bsdf.bxdfs[bxdf_idx] = Bxdf::DisDiff(DisneyDiffuse::new(
                            diffuse_weight * (1.0 - flat) * (1.0 - dt) * c,
                            Some(sc),
                        ));
                        bxdf_idx += 1;
                        bsdf.bxdfs[bxdf_idx] = Bxdf::DisSS(DisneyFakeSS::new(
                            diffuse_weight * flat * (1.0 - dt) * c,
                            rough,
                            Some(sc),
                        ));
                        bxdf_idx += 1;
                    } else {
                        bsdf.bxdfs[bxdf_idx] = Bxdf::DisDiff(DisneyDiffuse::new(
                            diffuse_weight * (1.0 - flat) * (1.0 - dt) * c,
                            None,
                        ));
                        bxdf_idx += 1;
                        bsdf.bxdfs[bxdf_idx] = Bxdf::DisSS(DisneyFakeSS::new(
                            diffuse_weight * flat * (1.0 - dt) * c,
                            rough,
                            None,
                        ));
                        bxdf_idx += 1;
                    }
                } else if sd.is_black() {
                    // No subsurface scattering; use regular (Fresnel modified) diffuse.
                    if use_scale {
                        bsdf.bxdfs[bxdf_idx] =
                            Bxdf::DisDiff(DisneyDiffuse::new(diffuse_weight * c, Some(sc)));
                        bxdf_idx += 1;
                    } else {
                        bsdf.bxdfs[bxdf_idx] =
                            Bxdf::DisDiff(DisneyDiffuse::new(diffuse_weight * c, None));
                        bxdf_idx += 1;
                    }
                } else {
                    // Use a BSSRDF instead.
                    if use_scale {
                        bsdf.bxdfs[bxdf_idx] = Bxdf::SpecTrans(SpecularTransmission::new(
                            Spectrum::from(1.0),
                            1.0,
                            e,
                            mode,
                            Some(sc),
                        ));
                        bxdf_idx += 1;
                    } else {
                        bsdf.bxdfs[bxdf_idx] = Bxdf::SpecTrans(SpecularTransmission::new(
                            Spectrum::from(1.0),
                            1.0,
                            e,
                            mode,
                            None,
                        ));
                        bxdf_idx += 1;
                    }
                    // TODO: BSSRDF
                }

                // Retro-reflection.
                if use_scale {
                    bsdf.bxdfs[bxdf_idx] =
                        Bxdf::DisRetro(DisneyRetro::new(diffuse_weight * c, rough, Some(sc)));
                    bxdf_idx += 1;
                } else {
                    bsdf.bxdfs[bxdf_idx] =
                        Bxdf::DisRetro(DisneyRetro::new(diffuse_weight * c, rough, None));
                    bxdf_idx += 1;
                }
                // Sheen (if enabled).
                if sheen_weight > 0.0 {
                    if use_scale {
                        bsdf.bxdfs[bxdf_idx] = Bxdf::DisSheen(DisneySheen::new(
                            diffuse_weight * sheen_weight * c_sheen,
                            Some(sc),
                        ));
                        bxdf_idx += 1;
                    } else {
                        bsdf.bxdfs[bxdf_idx] = Bxdf::DisSheen(DisneySheen::new(
                            diffuse_weight * sheen_weight * c_sheen,
                            None,
                        ));
                        bxdf_idx += 1;
                    }
                }
            }

            // Create the microfacet distribution for metallic and/or specular transmission.
            let ax = Float::max(0.001, sqr(rough) / aspect);
            let ay = Float::max(0.001, sqr(rough) * aspect);
            let distrib =
                MicrofacetDistribution::DisneyMicrofacet(DisneyMicrofacetDistribution::new(ax, ay));

            // Specular is Trowbridge-Reitz with a modified Fresnel function
            let cspec0 = lerp(
                metallic_weight,
                schlick_r0_from_eta(e) * lerp(spec_tint, Spectrum::new(1.0), c_tint),
                c,
            );
            let fresnel = Fresnel::Disney(DisneyFresnel::new(cspec0, metallic_weight, e));
            if use_scale {
                bsdf.bxdfs[bxdf_idx] =
                    Bxdf::MicrofacetRefl(MicrofacetReflection::new(c, distrib, fresnel, Some(sc)));
                bxdf_idx += 1;
            } else {
                bsdf.bxdfs[bxdf_idx] =
                    Bxdf::MicrofacetRefl(MicrofacetReflection::new(c, distrib, fresnel, None));
                bxdf_idx += 1;
            }
            // Clearcoat
            if cc > 0.0 {
                if use_scale {
                    bsdf.bxdfs[bxdf_idx] =
                        Bxdf::DisClearCoat(DisneyClearCoat::new(cc, gloss, Some(sc)));
                    bxdf_idx += 1;
                } else {
                    bsdf.bxdfs[bxdf_idx] =
                        Bxdf::DisClearCoat(DisneyClearCoat::new(cc, gloss, None));
                    bxdf_idx += 1;
                }
            }

            // BTDF
            if strans > 0.0 {
                // Walter et al.'s model, with the provided transmissive term scaled by sqrt(color), so
                // that after two refractions we're back to the provided color.
                let t = strans * c.sqrt();
                if self.thin {
                    // Scale roughness based on IOR (Burley 2015, Figure 15).
                    let rscaled = (0.65 * e - 0.35) * rough;
                    let ax = Float::max(0.001, sqr(rscaled) / aspect);
                    let ay = Float::max(0.001, sqr(rscaled) * aspect);
                    let scaled_distrib = MicrofacetDistribution::TrowbridgeReitz(
                        TrowbridgeReitzDistribution::new(ax, ay, true),
                    );
                    if use_scale {
                        bsdf.bxdfs[bxdf_idx] = Bxdf::MicrofacetTrans(MicrofacetTransmission::new(
                            t,
                            scaled_distrib,
                            1.0,
                            e,
                            mode,
                            Some(sc),
                        ));
                        bxdf_idx += 1;
                    } else {
                        bsdf.bxdfs[bxdf_idx] = Bxdf::MicrofacetTrans(MicrofacetTransmission::new(
                            t,
                            scaled_distrib,
                            1.0,
                            e,
                            mode,
                            None,
                        ));
                        bxdf_idx += 1;
                    }
                } else {
                    let distrib = MicrofacetDistribution::DisneyMicrofacet(
                        DisneyMicrofacetDistribution::new(ax, ay),
                    );
                    if use_scale {
                        bsdf.bxdfs[bxdf_idx] = Bxdf::MicrofacetTrans(MicrofacetTransmission::new(
                            t,
                            distrib,
                            1.0,
                            e,
                            mode,
                            Some(sc),
                        ));
                        bxdf_idx += 1;
                    } else {
                        bsdf.bxdfs[bxdf_idx] = Bxdf::MicrofacetTrans(MicrofacetTransmission::new(
                            t, distrib, 1.0, e, mode, None,
                        ));
                        bxdf_idx += 1;
                    }
                }
            }

            if self.thin {
                // Lambertian, weighted by (1.0 - diff_trans}
                if use_scale {
                    bsdf.bxdfs[bxdf_idx] =
                        Bxdf::LambertianTrans(LambertianTransmission::new(dt * c, Some(sc)));
                // bxdf_idx += 1;
                } else {
                    bsdf.bxdfs[bxdf_idx] =
                        Bxdf::LambertianTrans(LambertianTransmission::new(dt * c, None));
                    // bxdf_idx += 1;
                }
            }
        }
    }
}
// DisneyDiffuse

#[derive(Debug, Clone, Copy)]
pub struct DisneyDiffuse {
    pub r: Spectrum,
    pub sc_opt: Option<Spectrum>,
}

impl DisneyDiffuse {
    pub fn new(r: Spectrum, sc_opt: Option<Spectrum>) -> Self {
        DisneyDiffuse { r, sc_opt }
    }
    pub fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        let fo = schlick_weight(abs_cos_theta(wo));
        let fi = schlick_weight(abs_cos_theta(wi));

        // Diffuse fresnel - go from 1 at normal incidence to .5 at grazing.
        // Burley 2015, eq (4).
        if let Some(sc) = self.sc_opt {
            sc * self.r * f32::consts::FRAC_1_PI * (1.0 - fo / 2.0) * (1.0 - fi / 2.0)
        } else {
            self.r * f32::consts::FRAC_1_PI * (1.0 - fo / 2.0) * (1.0 - fi / 2.0)
        }
    }
    pub fn get_type(&self) -> u8 {
        BxdfType::BsdfReflection as u8 | BxdfType::BsdfDiffuse as u8
    }
}

// DisneyFakeSS

#[derive(Debug, Clone, Copy)]
pub struct DisneyFakeSS {
    pub r: Spectrum,
    pub roughness: Float,
    pub sc_opt: Option<Spectrum>,
}

impl DisneyFakeSS {
    pub fn new(r: Spectrum, roughness: Float, sc_opt: Option<Spectrum>) -> Self {
        DisneyFakeSS {
            r,
            roughness,
            sc_opt,
        }
    }
    pub fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        let mut wh = *wi + *wo;
        if wh.x == 0.0 && wh.y == 0.0 && wh.z == 0.0 {
            return Spectrum::from(0.0);
        }
        wh = wh.normalize();
        let cos_theta_d = vec3_dot_vec3(wi, &wh);

        // Fss90 used to "flatten" retroreflection based on roughness
        let fss90 = cos_theta_d * cos_theta_d * self.roughness;
        let fo = schlick_weight(abs_cos_theta(wo));
        let fi = schlick_weight(abs_cos_theta(wi));
        let fss = lerp(fo, 1.0, fss90) * lerp(fi, 1.0, fss90);
        // 1.25 scale is used to (roughly) preserve albedo
        let ss = 1.25 * (fss * (1.0 / (abs_cos_theta(wo) + abs_cos_theta(wi)) - 0.5) + 0.5);

        if let Some(sc) = self.sc_opt {
            sc * self.r * f32::consts::FRAC_1_PI * ss
        } else {
            self.r * f32::consts::FRAC_1_PI * ss
        }
    }
    pub fn get_type(&self) -> u8 {
        BxdfType::BsdfReflection as u8 | BxdfType::BsdfDiffuse as u8
    }
}

// DisneyRetro

#[derive(Debug, Clone, Copy)]
pub struct DisneyRetro {
    pub r: Spectrum,
    pub roughness: Float,
    pub sc_opt: Option<Spectrum>,
}

impl DisneyRetro {
    pub fn new(r: Spectrum, roughness: Float, sc_opt: Option<Spectrum>) -> Self {
        DisneyRetro {
            r,
            roughness,
            sc_opt,
        }
    }
    pub fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        let mut wh = *wi + *wo;
        if wh.x == 0.0 && wh.y == 0.0 && wh.z == 0.0 {
            return Spectrum::from(0.0);
        }
        wh = wh.normalize();
        let cos_theta_d = vec3_dot_vec3(wi, &wh);
        let fo = schlick_weight(abs_cos_theta(wo));
        let fi = schlick_weight(abs_cos_theta(wi));
        let rr = 2.0 * self.roughness * cos_theta_d * cos_theta_d;

        // Burley 2015, eq (4).
        if let Some(sc) = self.sc_opt {
            sc * self.r * f32::consts::FRAC_1_PI * rr * (fo + fi + fo * fi * (rr - 1.0))
        } else {
            self.r * f32::consts::FRAC_1_PI * rr * (fo + fi + fo * fi * (rr - 1.0))
        }
    }
    pub fn get_type(&self) -> u8 {
        BxdfType::BsdfReflection as u8 | BxdfType::BsdfDiffuse as u8
    }
}

// DisneySheen

#[derive(Debug, Clone, Copy)]
pub struct DisneySheen {
    pub r: Spectrum,
    pub sc_opt: Option<Spectrum>,
}

impl DisneySheen {
    pub fn new(r: Spectrum, sc_opt: Option<Spectrum>) -> Self {
        DisneySheen { r, sc_opt }
    }
    pub fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        let mut wh = *wi + *wo;
        if wh.x == 0.0 && wh.y == 0.0 && wh.z == 0.0 {
            return Spectrum::from(0.0);
        }
        wh = wh.normalize();
        let cos_theta_d = vec3_dot_vec3(wi, &wh);

        if let Some(sc) = self.sc_opt {
            sc * self.r * schlick_weight(cos_theta_d)
        } else {
            self.r * schlick_weight(cos_theta_d)
        }
    }
    pub fn get_type(&self) -> u8 {
        BxdfType::BsdfReflection as u8 | BxdfType::BsdfDiffuse as u8
    }
}

// DisneyClearCoat

#[derive(Debug, Clone, Copy)]
pub struct DisneyClearCoat {
    pub weight: Float,
    pub gloss: Float,
    pub sc_opt: Option<Spectrum>,
}

impl DisneyClearCoat {
    pub fn new(weight: Float, gloss: Float, sc_opt: Option<Spectrum>) -> Self {
        DisneyClearCoat {
            weight,
            gloss,
            sc_opt,
        }
    }
    pub fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        let mut wh = *wi + *wo;
        if wh.x == 0.0 && wh.y == 0.0 && wh.z == 0.0 {
            return Spectrum::from(0.0);
        }
        wh = wh.normalize();

        // Clearcoat has ior = 1.5 hardcoded -> F0 = 0.04. It then uses the
        // gtr1 distribution, which has even fatter tails than Trowbridge-Reitz
        // (which is GTR2).
        let dr = gtr1(abs_cos_theta(&wh), self.gloss);
        let fr = fr_schlick(0.04, vec3_dot_vec3(wo, &wh));
        // The geometric term always based on alpha = 0.25.
        let gr = smith_g_ggx(abs_cos_theta(wo), 0.25) * smith_g_ggx(abs_cos_theta(wi), 0.25);

        if let Some(sc) = self.sc_opt {
            sc * Spectrum::from(self.weight * gr * fr * dr / 4.0)
        } else {
            Spectrum::from(self.weight * gr * fr * dr / 4.0)
        }
    }
    pub fn sample_f(
        &self,
        wo: &Vector3f,
        wi: &mut Vector3f,
        u: &Point2f,
        pdf: &mut Float,
        _sampled_type: &mut u8,
    ) -> Spectrum {
        if wo.z == 0.0 {
            return Spectrum::zero();
        }

        let alpha2 = self.gloss * self.gloss;
        let cos_theta = Float::sqrt(Float::max(
            0.0,
            (1.0 - Float::powf(alpha2, 1.0 - u[0])) / (1.0 - alpha2),
        ));
        let sin_theta = Float::sqrt(Float::max(0.0, 1.0 - cos_theta * cos_theta));
        let phi = 2.0 * f32::consts::PI * u[1];
        let mut wh = spherical_direction(sin_theta, cos_theta, phi);
        if !vec3_same_hemisphere_vec3(wo, &wh) {
            wh = -wh;
        }
        *wi = reflect(wo, &wh);

        if !vec3_same_hemisphere_vec3(wo, wi) {
            return Spectrum::zero();
        }

        *pdf = self.pdf(wo, &wi);

        if let Some(sc) = self.sc_opt {
            sc * self.f(wo, wi)
        } else {
            self.f(wo, wi)
        }
    }
    pub fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        if !vec3_same_hemisphere_vec3(wo, wi) {
            return 0.0;
        }

        let mut wh = *wo + *wi;
        if wh.x == 0.0 && wh.y == 0.0 && wh.z == 0.0 {
            return 0.0;
        }
        wh = wh.normalize();

        // The sampling routine samples wh exactly from the gtr1 distribution.
        // Thus, the final value of the PDF is just the value of the
        // distribution for wh converted to a mesure with respect to the
        // surface normal.
        let dr = gtr1(abs_cos_theta(&wh), self.gloss);
        dr * abs_cos_theta(&wh) / (4.0 * vec3_dot_vec3(wo, &wh))
    }
    pub fn get_type(&self) -> u8 {
        BxdfType::BsdfReflection as u8 | BxdfType::BsdfGlossy as u8
    }
}

#[derive(Default, Clone, Copy)]
pub struct DisneyMicrofacetDistribution {
    pub inner: TrowbridgeReitzDistribution,
}

impl DisneyMicrofacetDistribution {
    pub fn new(alphax: Float, alphay: Float) -> DisneyMicrofacetDistribution {
        DisneyMicrofacetDistribution {
            inner: TrowbridgeReitzDistribution::new(alphax, alphay, true),
        }
    }
    pub fn d(&self, wh: &Vector3f) -> Float {
        self.inner.d(wh)
    }
    pub fn lambda(&self, wh: &Vector3f) -> Float {
        self.inner.lambda(wh)
    }
    pub fn g1(&self, w: &Vector3f) -> Float {
        1.0 as Float / (1.0 as Float + self.lambda(w))
    }
    pub fn g(&self, wi: &Vector3f, wo: &Vector3f) -> Float {
        // Disney uses the separable masking-shadowing model.
        self.g1(wi) * self.g1(wo)
    }
    pub fn pdf(&self, wo: &Vector3f, wh: &Vector3f) -> Float {
        if self.get_sample_visible_area() {
            self.d(wh) * self.g1(wo) * vec3_abs_dot_vec3(wo, wh) / abs_cos_theta(wo)
        } else {
            self.d(wh) * abs_cos_theta(wh)
        }
    }
    pub fn sample_wh(&self, wo: &Vector3f, u: &Point2f) -> Vector3f {
        self.inner.sample_wh(wo, u)
    }
    pub fn get_sample_visible_area(&self) -> bool {
        self.inner.get_sample_visible_area()
    }
}

/// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
///
/// The Schlick Fresnel approximation is:
///
/// R = R(0) + (1 - R(0)) (1 - cos theta)^5,
///
/// where R(0) is the reflectance at normal indicence.
#[inline]
fn schlick_weight(cos_theta: Float) -> Float {
    let m = clamp_t(1.0 - cos_theta, 0.0, 1.0);
    (m * m) * (m * m) * m
}

#[inline]
// For a dielectric, R(0) = (eta - 1)^2 / (eta + 1)^2, assuming we're
// coming from air.
fn schlick_r0_from_eta(eta: Float) -> Float {
    sqr(eta - 1.0) / sqr(eta + 1.0)
}

#[inline]
fn gtr1(cos_theta: Float, alpha: Float) -> Float {
    let alpha2 = alpha * alpha;

    (alpha2 - 1.0)
        / (f32::consts::PI * Float::log10(alpha2) * (1.0 + (alpha2 - 1.0) * cos_theta * cos_theta))
}

#[inline]
fn smith_g_ggx(cos_theta: Float, alpha: Float) -> Float {
    let alpha2 = alpha * alpha;
    let cos_theta2 = cos_theta * cos_theta;

    1.0 / (cos_theta + Float::sqrt(alpha2 + cos_theta2 - alpha2 * cos_theta2))
}

#[inline]
fn sqr(x: Float) -> Float {
    x * x
}
