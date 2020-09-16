//! Reflection models based on microfacets that exhibit perfect
//! specular reflection and transmission have been effective at
//! modeling light scattering from a variety of glossy materials,
//! including metals, plastic, and frosted glass.

// std
use std::f32::consts::PI;
// pbrt
use crate::core::geometry::{spherical_direction, vec3_abs_dot_vec3f};
use crate::core::geometry::{Point2f, Vector3f, XYEnum};
use crate::core::pbrt::Float;
use crate::core::pbrt::{erf, erf_inv};
use crate::core::reflection::{
    abs_cos_theta, cos_2_phi, cos_2_theta, cos_phi, cos_theta, sin_2_phi, sin_phi, tan_2_theta,
    tan_theta, vec3_same_hemisphere_vec3,
};
use crate::materials::disney::DisneyMicrofacetDistribution;

// see microfacet.h

#[derive(Copy, Clone)]
pub enum MicrofacetDistribution {
    Beckmann(BeckmannDistribution),
    TrowbridgeReitz(TrowbridgeReitzDistribution),
    // see disney.rs
    DisneyMicrofacet(DisneyMicrofacetDistribution),
}

impl MicrofacetDistribution {
    pub fn d(&self, wh: &Vector3f) -> Float {
        match self {
            MicrofacetDistribution::Beckmann(distribution) => distribution.d(wh),
            MicrofacetDistribution::TrowbridgeReitz(distribution) => distribution.d(wh),
            MicrofacetDistribution::DisneyMicrofacet(distribution) => distribution.d(wh),
        }
    }
    pub fn lambda(&self, w: &Vector3f) -> Float {
        match self {
            MicrofacetDistribution::Beckmann(distribution) => distribution.lambda(w),
            MicrofacetDistribution::TrowbridgeReitz(distribution) => distribution.lambda(w),
            MicrofacetDistribution::DisneyMicrofacet(distribution) => distribution.lambda(w),
        }
    }
    pub fn g1(&self, w: &Vector3f) -> Float {
        match self {
            MicrofacetDistribution::Beckmann(distribution) => distribution.g1(w),
            MicrofacetDistribution::TrowbridgeReitz(distribution) => distribution.g1(w),
            MicrofacetDistribution::DisneyMicrofacet(distribution) => distribution.g1(w),
        }
    }
    pub fn g(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        match self {
            MicrofacetDistribution::Beckmann(distribution) => distribution.g(wo, wi),
            MicrofacetDistribution::TrowbridgeReitz(distribution) => distribution.g(wo, wi),
            MicrofacetDistribution::DisneyMicrofacet(distribution) => distribution.g(wo, wi),
        }
    }
    pub fn pdf(&self, wo: &Vector3f, wh: &Vector3f) -> Float {
        match self {
            MicrofacetDistribution::Beckmann(distribution) => distribution.pdf(wo, wh),
            MicrofacetDistribution::TrowbridgeReitz(distribution) => distribution.pdf(wo, wh),
            MicrofacetDistribution::DisneyMicrofacet(distribution) => distribution.pdf(wo, wh),
        }
    }
    pub fn sample_wh(&self, wo: &Vector3f, u: &Point2f) -> Vector3f {
        match self {
            MicrofacetDistribution::Beckmann(distribution) => distribution.sample_wh(wo, u),
            MicrofacetDistribution::TrowbridgeReitz(distribution) => distribution.sample_wh(wo, u),
            MicrofacetDistribution::DisneyMicrofacet(distribution) => distribution.sample_wh(wo, u),
        }
    }
    pub fn get_sample_visible_area(&self) -> bool {
        match self {
            MicrofacetDistribution::Beckmann(distribution) => {
                distribution.get_sample_visible_area()
            }
            MicrofacetDistribution::TrowbridgeReitz(distribution) => {
                distribution.get_sample_visible_area()
            }
            MicrofacetDistribution::DisneyMicrofacet(distribution) => {
                distribution.get_sample_visible_area()
            }
        }
    }
}

#[derive(Default, Copy, Clone)]
pub struct BeckmannDistribution {
    pub alpha_x: Float,
    pub alpha_y: Float,
    // inherited from class MicrofacetDistribution (see microfacet.h)
    pub sample_visible_area: bool,
}

impl BeckmannDistribution {
    pub fn new(alpha_x: Float, alpha_y: Float, sample_visible_area: bool) -> Self {
        BeckmannDistribution {
            alpha_x: alpha_x.max(0.001 as Float),
            alpha_y: alpha_y.max(0.001 as Float),
            sample_visible_area,
        }
    }
    pub fn roughness_to_alpha(roughness: Float) -> Float {
        let mut roughness = roughness;
        let limit: Float = 1e-3 as Float;
        if limit > roughness {
            roughness = limit;
        }
        let x: Float = roughness.ln(); // natural (base e) logarithm
        1.62142
            + 0.819_955 * x
            + 0.1734 * x * x
            + 0.017_120_1 * x * x * x
            + 0.000_640_711 * x * x * x * x
    }
    pub fn d(&self, wh: &Vector3f) -> Float {
        let tan_2_theta: Float = tan_2_theta(wh);
        if tan_2_theta.is_infinite() {
            return 0.0 as Float;
        }
        let cos_4_theta: Float = cos_2_theta(wh) * cos_2_theta(wh);
        (-tan_2_theta
            * (cos_2_phi(wh) / (self.alpha_x * self.alpha_x)
                + sin_2_phi(wh) / (self.alpha_y * self.alpha_y)))
            .exp()
            / (PI * self.alpha_x * self.alpha_y * cos_4_theta)
    }
    pub fn lambda(&self, w: &Vector3f) -> Float {
        let abs_tan_theta: Float = tan_theta(w).abs();
        if abs_tan_theta.is_infinite() {
            return 0.0;
        }
        // compute _alpha_ for direction _w_
        let alpha: Float = (cos_2_phi(w) * self.alpha_x * self.alpha_x
            + sin_2_phi(w) * self.alpha_y * self.alpha_y)
            .sqrt();
        let a: Float = 1.0 as Float / (alpha * abs_tan_theta);
        if a >= 1.6 as Float {
            return 0.0 as Float;
        }
        (1.0 as Float - 1.259 * a + 0.396 * a * a) / (3.535 * a + 2.181 * a * a)
    }
    pub fn g1(&self, w: &Vector3f) -> Float {
        1.0 as Float / (1.0 as Float + self.lambda(w))
    }
    pub fn g(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        1.0 as Float / (1.0 as Float + self.lambda(wo) + self.lambda(wi))
    }
    pub fn pdf(&self, wo: &Vector3f, wh: &Vector3f) -> Float {
        if self.get_sample_visible_area() {
            self.d(wh) * self.g1(wo) * vec3_abs_dot_vec3f(wo, wh) / abs_cos_theta(wo)
        } else {
            self.d(wh) * abs_cos_theta(wh)
        }
    }
    pub fn sample_wh(&self, wo: &Vector3f, u: &Point2f) -> Vector3f {
        if !self.sample_visible_area {
            // sample full distribution of normals for Beckmann
            // distribution

            // compute $\tan^2 \theta$ and $\phi$ for Beckmann
            // distribution sample
            let tan_2_theta: Float;
            let mut phi: Float;
            if self.alpha_x == self.alpha_y {
                let log_sample: Float = (1.0 as Float - u[XYEnum::X]).ln();
                assert!(!log_sample.is_infinite());
                tan_2_theta = -self.alpha_x * self.alpha_x * log_sample;
                phi = u[XYEnum::Y] * 2.0 as Float * PI;
            } else {
                // compute _tan_2_theta_ and _phi_ for anisotropic
                // Beckmann distribution
                let log_sample: Float = (1.0 as Float - u[XYEnum::X]).ln();
                assert!(!log_sample.is_infinite());
                phi = (self.alpha_y / self.alpha_x
                    * (2.0 as Float * PI * u[XYEnum::Y] + 0.5 * PI).tan())
                .atan();
                if u[XYEnum::Y] > 0.5 as Float {
                    phi += PI;
                }
                let sin_phi: Float = phi.sin();
                let cos_phi: Float = phi.cos();
                let alpha_x2: Float = self.alpha_x * self.alpha_x;
                let alpha_y2: Float = self.alpha_y * self.alpha_y;
                tan_2_theta =
                    -log_sample / (cos_phi * cos_phi / alpha_x2 + sin_phi * sin_phi / alpha_y2);
            }
            // map sampled Beckmann angles to normal direction _wh_
            let cos_theta: Float = 1.0 as Float / (1.0 as Float + tan_2_theta).sqrt();
            let sin_theta: Float =
                ((0.0 as Float).max(1.0 as Float - cos_theta * cos_theta)).sqrt();
            let mut wh: Vector3f = spherical_direction(sin_theta, cos_theta, phi);
            if !vec3_same_hemisphere_vec3(wo, &wh) {
                wh = -wh;
            }
            wh
        } else {
            // sample visible area of normals for Beckmann distribution
            let mut wh: Vector3f;
            let flip: bool = wo.z < 0.0 as Float;
            if flip {
                wh = beckmann_sample(
                    &-(*wo),
                    self.alpha_x,
                    self.alpha_y,
                    u[XYEnum::X],
                    u[XYEnum::Y],
                );
            } else {
                wh = beckmann_sample(wo, self.alpha_x, self.alpha_y, u[XYEnum::X], u[XYEnum::Y]);
            }
            if flip {
                wh = -wh;
            }
            wh
        }
    }
    pub fn get_sample_visible_area(&self) -> bool {
        self.sample_visible_area
    }
}

#[derive(Default, Copy, Clone)]
pub struct TrowbridgeReitzDistribution {
    pub alpha_x: Float,
    pub alpha_y: Float,
    // inherited from class MicrofacetDistribution (see microfacet.h)
    pub sample_visible_area: bool,
}

impl TrowbridgeReitzDistribution {
    pub fn new(alpha_x: Float, alpha_y: Float, sample_visible_area: bool) -> Self {
        TrowbridgeReitzDistribution {
            alpha_x: alpha_x.max(0.001 as Float),
            alpha_y: alpha_y.max(0.001 as Float),
            sample_visible_area,
        }
    }
    /// Microfacet distribution function: In comparison to the
    /// Beckmann-Spizzichino model, Trowbridge-Reitz has higher tails - it
    /// falls off to zero more slowly for directions far from the surface
    /// normal.
    pub fn roughness_to_alpha(roughness: Float) -> Float {
        let mut roughness = roughness;
        let limit: Float = 1e-3 as Float;
        if limit > roughness {
            roughness = limit;
        }
        let x: Float = roughness.ln(); // natural (base e) logarithm
        1.62142
            + 0.819_955 * x
            + 0.1734 * x * x
            + 0.017_120_1 * x * x * x
            + 0.000_640_711 * x * x * x * x
    }
    pub fn d(&self, wh: &Vector3f) -> Float {
        let tan_2_theta: Float = tan_2_theta(wh);
        if tan_2_theta.is_infinite() {
            return 0.0 as Float;
        }
        let cos_4_theta: Float = cos_2_theta(wh) * cos_2_theta(wh);
        let e: Float = (cos_2_phi(wh) / (self.alpha_x * self.alpha_x)
            + sin_2_phi(wh) / (self.alpha_y * self.alpha_y))
            * tan_2_theta;
        1.0 as Float
            / (PI
                * self.alpha_x
                * self.alpha_y
                * cos_4_theta
                * (1.0 as Float + e)
                * (1.0 as Float + e))
    }
    pub fn lambda(&self, w: &Vector3f) -> Float {
        let abs_tan_theta: Float = tan_theta(w).abs();
        if abs_tan_theta.is_infinite() {
            return 0.0;
        }
        // compute _alpha_ for direction _w_
        let alpha: Float = (cos_2_phi(w) * self.alpha_x * self.alpha_x
            + sin_2_phi(w) * self.alpha_y * self.alpha_y)
            .sqrt();
        let alpha_2_tan_2_theta: Float = (alpha * abs_tan_theta) * (alpha * abs_tan_theta);
        (-1.0 as Float + (1.0 as Float + alpha_2_tan_2_theta).sqrt()) / 2.0 as Float
    }
    pub fn g1(&self, w: &Vector3f) -> Float {
        1.0 as Float / (1.0 as Float + self.lambda(w))
    }
    pub fn g(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        1.0 as Float / (1.0 as Float + self.lambda(wo) + self.lambda(wi))
    }
    pub fn pdf(&self, wo: &Vector3f, wh: &Vector3f) -> Float {
        if self.get_sample_visible_area() {
            self.d(wh) * self.g1(wo) * vec3_abs_dot_vec3f(wo, wh) / abs_cos_theta(wo)
        } else {
            self.d(wh) * abs_cos_theta(wh)
        }
    }
    pub fn sample_wh(&self, wo: &Vector3f, u: &Point2f) -> Vector3f {
        let mut wh: Vector3f;
        if !self.sample_visible_area {
            let cos_theta;
            let mut phi: Float = (2.0 * PI) * u[XYEnum::Y];
            if self.alpha_x == self.alpha_y {
                let tan_theta2: Float =
                    self.alpha_x * self.alpha_x * u[XYEnum::X] / (1.0 - u[XYEnum::X]);
                cos_theta = 1.0 / (1.0 + tan_theta2).sqrt();
            } else {
                phi = (self.alpha_y / self.alpha_x * (2.0 * PI * u[XYEnum::Y] + 0.5 * PI).tan())
                    .atan();
                if u[XYEnum::Y] > 0.5 {
                    phi += PI;
                }
                let sin_phi: Float = phi.sin();
                let cos_phi: Float = phi.cos();
                let alphax2: Float = self.alpha_x * self.alpha_x;
                let alphay2: Float = self.alpha_y * self.alpha_y;
                let alpha2: Float =
                    1.0 / (cos_phi * cos_phi / alphax2 + sin_phi * sin_phi / alphay2);
                let tan_theta2: Float = alpha2 * u[XYEnum::X] / (1.0 - u[XYEnum::X]);
                cos_theta = 1.0 / (1.0 + tan_theta2).sqrt();
            }
            let sin_theta: Float = (0.0 as Float).max(1.0 - cos_theta * cos_theta).sqrt();
            wh = spherical_direction(sin_theta, cos_theta, phi);
            if !vec3_same_hemisphere_vec3(wo, &wh) {
                wh = -wh;
            }
        } else {
            let flip: bool = wo.z < 0.0;
            if flip {
                wh = trowbridge_reitz_sample(
                    &-(*wo),
                    self.alpha_x,
                    self.alpha_y,
                    u[XYEnum::X],
                    u[XYEnum::Y],
                );
                wh = -wh;
            } else {
                wh = trowbridge_reitz_sample(
                    wo,
                    self.alpha_x,
                    self.alpha_y,
                    u[XYEnum::X],
                    u[XYEnum::Y],
                );
            }
        }
        wh
    }
    pub fn get_sample_visible_area(&self) -> bool {
        self.sample_visible_area
    }
}

fn beckmann_sample_11(
    cos_theta_i: Float,
    u1: Float,
    u2: Float,
    slope_x: &mut Float,
    slope_y: &mut Float,
) {
    // special case (normal incidence)
    if cos_theta_i > 0.9999 {
        let r: Float = (-((1.0 - u1).ln())).sqrt();
        let phi: Float = 2.0 as Float * PI * u2;
        *slope_x = r * phi.cos();
        *slope_y = r * phi.sin();
        return;
    }

    // The original inversion routine from the paper contained
    // discontinuities, which causes issues for QMC integration and
    // techniques like Kelemen-style MLT. The following code performs
    // a numerical inversion with better behavior
    let sin_theta_i: Float = (0.0 as Float)
        .max(1.0 as Float - cos_theta_i * cos_theta_i)
        .sqrt();
    let tan_theta_i: Float = sin_theta_i / cos_theta_i;
    let cot_theta_i: Float = 1.0 as Float / tan_theta_i;

    // Search interval -- everything is parameterized in the Erf() domain

    let mut a: Float = -1.0;
    let mut c: Float = erf(cot_theta_i);
    let sample_x: Float = u1.max(1e-6 as Float);

    // We can do better (inverse of an approximation computed in
    // Mathematica)
    let theta_i: Float = cos_theta_i.acos();
    let fit: Float = 1.0 as Float + theta_i * (-0.876 + theta_i * (0.4265 - 0.0594 * theta_i));
    let mut b: Float = c - (1.0 as Float + c) * Float::powf(1.0 as Float - sample_x, fit);

    // normalization factor for the CDF
    let sqrt_pi_inv: Float = 1.0 as Float / PI.sqrt();
    let normalization: Float = 1.0 as Float
        / (1.0 as Float + c + sqrt_pi_inv * tan_theta_i * (-cot_theta_i * cot_theta_i).exp());

    for _it in 0..10 {
        // bisection criterion -- the oddly-looking Boolean expression
        // are intentional to check for NaNs at little additional cost
        if !(b >= a && b <= c) {
            b = 0.5 as Float * (a + c);
        }
        // evaluate the CDF and its derivative (i.e. the density
        // function)
        let inv_erf: Float = erf_inv(b);
        let value: Float = normalization
            * (1.0 as Float + b + sqrt_pi_inv * tan_theta_i * (-inv_erf * inv_erf).exp())
            - sample_x;
        let derivative: Float = normalization * (1.0 as Float - inv_erf * tan_theta_i);

        if value.abs() < 1e-5 as Float {
            break;
        }

        // update bisection intervals
        if value > 0.0 as Float {
            c = b;
        } else {
            a = b;
        }
        b -= value / derivative;
    }

    // now convert back into a slope value
    *slope_x = erf_inv(b);

    // simulate Y component
    *slope_y = erf_inv(2.0 as Float * u2.max(1e-6 as Float) - 1.0 as Float);

    assert!(!(*slope_x).is_infinite());
    assert!(!(*slope_x).is_nan());
    assert!(!(*slope_y).is_infinite());
    assert!(!(*slope_y).is_nan());
}

fn beckmann_sample(
    wi: &Vector3f,
    alpha_x: Float,
    alpha_y: Float,
    u1: Float,
    u2: Float,
) -> Vector3f {
    // 1. stretch wi
    let wi_stretched: Vector3f = Vector3f {
        x: alpha_x * wi.x,
        y: alpha_y * wi.y,
        z: wi.z,
    }
    .normalize();

    // 2. simulate P22_{wi}(x_slope, y_slope, 1, 1)
    let mut slope_x: Float = 0.0;
    let mut slope_y: Float = 0.0;
    beckmann_sample_11(cos_theta(&wi_stretched), u1, u2, &mut slope_x, &mut slope_y);

    // 3. rotate
    let tmp: Float = cos_phi(&wi_stretched) * slope_x - sin_phi(&wi_stretched) * slope_y;
    slope_y = sin_phi(&wi_stretched) * slope_x + cos_phi(&wi_stretched) * slope_y;
    slope_x = tmp;

    // 4. unstretch
    slope_x *= alpha_x;
    slope_y *= alpha_y;

    // 5. compute normal
    Vector3f {
        x: -slope_x,
        y: -slope_y,
        z: 1.0,
    }
    .normalize()
}

fn trowbridge_reitz_sample_11(
    cos_theta: Float,
    u1: Float,
    u2: Float,
    slope_x: &mut Float,
    slope_y: &mut Float,
) {
    // special case (normal incidence)
    if cos_theta > 0.9999 {
        let r: Float = (u1 / (1.0 - u1)).sqrt();
        let phi: Float = 6.283_185_307_18 * u2;
        *slope_x = r * phi.cos();
        *slope_y = r * phi.sin();
        return;
    }

    let sin_theta: Float = (0.0 as Float)
        .max(1.0 as Float - cos_theta * cos_theta)
        .sqrt();
    let tan_theta: Float = sin_theta / cos_theta;
    let a: Float = 1.0 / tan_theta;
    let g1: Float = 2.0 / (1.0 + (1.0 + 1.0 / (a * a)).sqrt());

    // sample slope_x
    let a: Float = 2.0 * u1 / g1 - 1.0;
    let mut tmp: Float = 1.0 / (a * a - 1.0);
    if tmp > 1e10 {
        tmp = 1e10;
    }
    let b: Float = tan_theta;
    let d: Float = (b * b * tmp * tmp - (a * a - b * b) * tmp)
        .max(0.0 as Float)
        .sqrt();
    let slope_x_1: Float = b * tmp - d;
    let slope_x_2: Float = b * tmp + d;
    if a < 0.0 || slope_x_2 > 1.0 / tan_theta {
        *slope_x = slope_x_1;
    } else {
        *slope_x = slope_x_2;
    }

    // sample slope_y
    let s: Float;
    let new_u2 = if u2 > 0.5 {
        s = 1.0;
        2.0 * (u2 - 0.5)
    } else {
        s = -1.0;
        2.0 * (0.5 - u2)
    };
    let z: Float = (new_u2 * (new_u2 * (new_u2 * 0.27385 - 0.73369) + 0.46341))
        / (new_u2 * (new_u2 * (new_u2 * 0.093_073 + 0.309_420) - 1.0) + 0.597_999);
    *slope_y = s * z * (1.0 + *slope_x * *slope_x).sqrt();

    assert!(!(*slope_y).is_infinite());
    assert!(!(*slope_y).is_nan());
}

fn trowbridge_reitz_sample(
    wi: &Vector3f,
    alpha_x: Float,
    alpha_y: Float,
    u1: Float,
    u2: Float,
) -> Vector3f {
    // 1. stretch wi
    let wi_stretched: Vector3f = Vector3f {
        x: alpha_x * wi.x,
        y: alpha_y * wi.y,
        z: wi.z,
    }
    .normalize();

    // 2. simulate P22_{wi}(x_slope, y_slope, 1, 1)
    let mut slope_x: Float = 0.0;
    let mut slope_y: Float = 0.0;
    trowbridge_reitz_sample_11(cos_theta(&wi_stretched), u1, u2, &mut slope_x, &mut slope_y);

    // 3. rotate
    let tmp: Float = cos_phi(&wi_stretched) * slope_x - sin_phi(&wi_stretched) * slope_y;
    slope_y = sin_phi(&wi_stretched) * slope_x + cos_phi(&wi_stretched) * slope_y;
    slope_x = tmp;

    // 4. unstretch
    slope_x *= alpha_x;
    slope_y *= alpha_y;

    // 5. compute normal
    Vector3f {
        x: -slope_x,
        y: -slope_y,
        z: 1.0,
    }
    .normalize()
}
