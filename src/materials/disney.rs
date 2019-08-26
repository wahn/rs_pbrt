use std::f32;
use std::sync::Arc;

use num::Zero;

use crate::core::geometry::{spherical_direction, vec3_dot_vec3, Point2f, Vector3f};
use crate::core::interaction::SurfaceInteraction;
use crate::core::material::{Material, TransportMode};
use crate::core::microfacet::{MicrofacetDistribution, TrowbridgeReitzDistribution};
use crate::core::paramset::TextureParams;
use crate::core::pbrt::{clamp_t, lerp, Float, Spectrum};
use crate::core::reflection::{
    abs_cos_theta, vec3_same_hemisphere_vec3, MicrofacetReflection, MicrofacetTransmission,
};
use crate::core::reflection::{fr_dielectric, reflect, Fresnel};
use crate::core::reflection::{Bsdf, Bxdf, BxdfType, LambertianTransmission, SpecularTransmission};
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
    pub fn create(mp: &mut TextureParams) -> Arc<dyn Material + Send + Sync> {
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

        Arc::new(DisneyMaterial {
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
        })
    }
}

impl Material for DisneyMaterial {
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        mode: TransportMode,
        _allow_multiple_lobes: bool,
        _material: Option<Arc<dyn Material + Send + Sync>>,
    ) {
        if let Some(ref bump) = self.bump_map {
            Self::bump(bump, si);
        }

        let mut bxdfs: Vec<Arc<dyn Bxdf + Send + Sync>> = Vec::new();

        // Diffuse
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

        if diffuse_weight > 0.0 {
            if self.thin {
                let flat = self.flatness.evaluate(si);
                // Blend between DisneyDiffuse and fake subsurface based on flatness. Additionally,
                // weight using diff_trans.
                bxdfs.push(Arc::new(DisneyDiffuse::new(
                    diffuse_weight * (1.0 - flat) * (1.0 - dt) * c,
                )));
                bxdfs.push(Arc::new(DisneyFakeSS::new(
                    diffuse_weight * flat * (1.0 - dt) * c,
                    rough,
                )));
            } else {
                let sd = self.scatter_distance.evaluate(si);
                if sd.is_black() {
                    // No subsurface scattering; use regular (Fresnel modified) diffuse.
                    bxdfs.push(Arc::new(DisneyDiffuse::new(diffuse_weight * c)));
                } else {
                    // Use a BSSRDF instead.
                    bxdfs.push(Arc::new(SpecularTransmission::new(
                        Spectrum::from(1.0),
                        1.0,
                        e,
                        mode,
                    )));
                    // TODO: BSSRDF
                }
            }

            // Retro-reflection.
            bxdfs.push(Arc::new(DisneyRetro::new(diffuse_weight * c, rough)));

            // Sheen (if enabled).
            if sheen_weight > 0.0 {
                bxdfs.push(Arc::new(DisneySheen::new(
                    diffuse_weight * sheen_weight * c_sheen,
                )));
            }
        }

        // Create the microfacet distribution for metallic and/or specular transmission.
        let aspect = Float::sqrt(1.0 - self.anisotropic.evaluate(si) * 0.9);
        let ax = Float::max(0.001, sqr(rough) / aspect);
        let ay = Float::max(0.001, sqr(rough) * aspect);
        let distrib = Arc::new(DisneyMicrofacetDistribution::new(ax, ay));

        // Specular is Trowbridge-Reitz with a modified Fresnel function
        let spec_tint = self.specular_tint.evaluate(si);
        let cspec0 = lerp(
            metallic_weight,
            schlick_r0_from_eta(e) * lerp(spec_tint, Spectrum::new(1.0), c_tint),
            c,
        );
        let fresnel = Arc::new(DisneyFresnel::new(cspec0, metallic_weight, e));
        bxdfs.push(Arc::new(MicrofacetReflection::new(
            c,
            distrib.clone(),
            fresnel,
        )));

        // Clearcoat
        let cc = self.clearcoat.evaluate(si);
        if cc > 0.0 {
            bxdfs.push(Arc::new(DisneyClearCoat::new(
                cc,
                lerp(self.clearcoat_gloss.evaluate(si), 0.1, 0.001),
            )));
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
                let scaled_distrib = Arc::new(TrowbridgeReitzDistribution::new(ax, ay, true));
                bxdfs.push(Arc::new(MicrofacetTransmission::new(
                    t,
                    scaled_distrib,
                    1.0,
                    e,
                    mode,
                )));
            } else {
                bxdfs.push(Arc::new(MicrofacetTransmission::new(
                    t,
                    distrib.clone(),
                    1.0,
                    e,
                    mode,
                )));
            }
        }

        if self.thin {
            // Lambertian, weighted by (1.0 - diff_trans}
            bxdfs.push(Arc::new(LambertianTransmission::new(dt * c)));
        }

        si.bsdf = Some(Arc::new(Bsdf::new(si, 1.0, bxdfs)));
    }
}

// DisneyDiffuse
#[derive(Debug, Clone, Copy)]
struct DisneyDiffuse {
    r: Spectrum,
}

impl DisneyDiffuse {
    pub fn new(r: Spectrum) -> DisneyDiffuse {
        DisneyDiffuse { r }
    }
}

impl Bxdf for DisneyDiffuse {
    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        let fo = schlick_weight(abs_cos_theta(wo));
        let fi = schlick_weight(abs_cos_theta(wi));

        // Diffuse fresnel - go from 1 at normal incidence to .5 at grazing.
        // Burley 2015, eq (4).
        self.r * f32::consts::FRAC_1_PI * (1.0 - fo / 2.0) * (1.0 - fi / 2.0)
    }

    fn get_type(&self) -> u8 {
        BxdfType::BsdfReflection as u8 | BxdfType::BsdfDiffuse as u8
    }
}

// DisneyFakeSS
#[derive(Debug, Clone, Copy)]
struct DisneyFakeSS {
    r: Spectrum,
    roughness: Float,
}

impl DisneyFakeSS {
    pub fn new(r: Spectrum, roughness: Float) -> DisneyFakeSS {
        DisneyFakeSS { r, roughness }
    }
}

impl Bxdf for DisneyFakeSS {
    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
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

        self.r * f32::consts::FRAC_1_PI * ss
    }

    fn get_type(&self) -> u8 {
        BxdfType::BsdfReflection as u8 | BxdfType::BsdfDiffuse as u8
    }
}

// DisneyRetro
#[derive(Debug, Clone, Copy)]
struct DisneyRetro {
    r: Spectrum,
    roughness: Float,
}

impl DisneyRetro {
    pub fn new(r: Spectrum, roughness: Float) -> DisneyRetro {
        DisneyRetro { r, roughness }
    }
}

impl Bxdf for DisneyRetro {
    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
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
        self.r * f32::consts::FRAC_1_PI * rr * (fo + fi + fo * fi * (rr - 1.0))
    }

    fn get_type(&self) -> u8 {
        BxdfType::BsdfReflection as u8 | BxdfType::BsdfDiffuse as u8
    }
}

// DisneySheen
#[derive(Debug, Clone, Copy)]
struct DisneySheen {
    r: Spectrum,
}

impl DisneySheen {
    pub fn new(r: Spectrum) -> DisneySheen {
        DisneySheen { r }
    }
}

impl Bxdf for DisneySheen {
    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        let mut wh = *wi + *wo;
        if wh.x == 0.0 && wh.y == 0.0 && wh.z == 0.0 {
            return Spectrum::from(0.0);
        }
        wh = wh.normalize();
        let cos_theta_d = vec3_dot_vec3(wi, &wh);

        self.r * schlick_weight(cos_theta_d)
    }

    fn get_type(&self) -> u8 {
        BxdfType::BsdfReflection as u8 | BxdfType::BsdfDiffuse as u8
    }
}

// DisneyClearCoat
#[derive(Debug, Clone, Copy)]
struct DisneyClearCoat {
    weight: Float,
    gloss: Float,
}

impl DisneyClearCoat {
    pub fn new(weight: Float, gloss: Float) -> DisneyClearCoat {
        DisneyClearCoat { weight, gloss }
    }
}

impl Bxdf for DisneyClearCoat {
    fn f(&self, wo: &Vector3f, wi: &Vector3f) -> Spectrum {
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

        Spectrum::from(self.weight * gr * fr * dr / 4.0)
    }

    fn sample_f(
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

        self.f(wo, wi)
    }

    fn pdf(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
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

    fn get_type(&self) -> u8 {
        BxdfType::BsdfReflection as u8 | BxdfType::BsdfGlossy as u8
    }
}

// DisneyFresnel

/// Specialized Fresnel function used for the specular component, based on
/// a mixture between dielectric and the Schlick Fresnel approximation.
#[derive(Debug, Clone, Copy)]
struct DisneyFresnel {
    r0: Spectrum,
    metallic: Float,
    eta: Float,
}

impl DisneyFresnel {
    pub fn new(r0: Spectrum, metallic: Float, eta: Float) -> DisneyFresnel {
        DisneyFresnel { r0, metallic, eta }
    }
}

impl Fresnel for DisneyFresnel {
    fn evaluate(&self, cos_i: Float) -> Spectrum {
        lerp(
            self.metallic,
            Spectrum::from(fr_dielectric(cos_i, 1.0, self.eta)),
            fr_schlick_spectrum(self.r0, cos_i),
        )
    }
}

struct DisneyMicrofacetDistribution {
    inner: TrowbridgeReitzDistribution,
}

impl DisneyMicrofacetDistribution {
    fn new(alphax: Float, alphay: Float) -> DisneyMicrofacetDistribution {
        DisneyMicrofacetDistribution {
            inner: TrowbridgeReitzDistribution::new(alphax, alphay, true),
        }
    }
}

impl MicrofacetDistribution for DisneyMicrofacetDistribution {
    fn d(&self, wh: &Vector3f) -> Float {
        self.inner.d(wh)
    }

    fn lambda(&self, wh: &Vector3f) -> Float {
        self.inner.lambda(wh)
    }

    fn g(&self, wi: &Vector3f, wo: &Vector3f) -> Float {
        // Disney uses the separable masking-shadowing model.
        self.g1(wi) * self.g1(wo)
    }

    fn sample_wh(&self, wo: &Vector3f, u: &Point2f) -> Vector3f {
        self.inner.sample_wh(wo, u)
    }

    fn get_sample_visible_area(&self) -> bool {
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
fn fr_schlick(r0: Float, cos_theta: Float) -> Float {
    lerp(schlick_weight(cos_theta), r0, 1.0)
}

#[inline]
fn fr_schlick_spectrum(r0: Spectrum, cos_theta: Float) -> Spectrum {
    lerp(schlick_weight(cos_theta), r0, Spectrum::from(1.0))
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
