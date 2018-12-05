// std
use std::f32::consts::PI;
// pbrt
use core::geometry::{spherical_direction, vec3_abs_dot_vec3};
use core::geometry::{Point2f, Vector3f};
use core::pbrt::Float;
use core::reflection::{
    abs_cos_theta, cos_2_phi, cos_2_theta, cos_phi, cos_theta, sin_2_phi, sin_phi, tan_2_theta,
    tan_theta, vec3_same_hemisphere_vec3,
};

// see microfacet.h

pub trait MicrofacetDistribution {
    fn d(&self, wh: &Vector3f) -> Float;
    fn lambda(&self, w: &Vector3f) -> Float;
    fn g1(&self, w: &Vector3f) -> Float {
        1.0 as Float / (1.0 as Float + self.lambda(w))
    }
    fn g(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        1.0 as Float / (1.0 as Float + self.lambda(wo) + self.lambda(wi))
    }
    fn pdf(&self, wo: &Vector3f, wh: &Vector3f) -> Float {
        if self.get_sample_visible_area() {
            self.d(wh) * self.g1(wo) * vec3_abs_dot_vec3(wo, wh) / abs_cos_theta(wo)
        } else {
            self.d(wh) * abs_cos_theta(wh)
        }
    }
    fn sample_wh(&self, wo: &Vector3f, u: &Point2f) -> Vector3f;
    fn get_sample_visible_area(&self) -> bool;
}

pub struct TrowbridgeReitzDistribution {
    pub alpha_x: Float,
    pub alpha_y: Float,
    // inherited from class MicrofacetDistribution (see microfacet.h)
    pub sample_visible_area: bool,
}

impl TrowbridgeReitzDistribution {
    pub fn new(alpha_x: Float, alpha_y: Float, sample_visible_area: bool) -> Self {
        TrowbridgeReitzDistribution {
            alpha_x: alpha_x,
            alpha_y: alpha_y,
            sample_visible_area: sample_visible_area,
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
            + 0.819955 * x
            + 0.1734 * x * x
            + 0.0171201 * x * x * x
            + 0.000640711 * x * x * x * x
    }
}

impl MicrofacetDistribution for TrowbridgeReitzDistribution {
    fn d(&self, wh: &Vector3f) -> Float {
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

    fn lambda(&self, w: &Vector3f) -> Float {
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

    fn sample_wh(&self, wo: &Vector3f, u: &Point2f) -> Vector3f {
        let mut wh: Vector3f;
        if !self.sample_visible_area {
            let cos_theta;
            let mut phi: Float = (2.0 * PI) * u[1];
            if self.alpha_x == self.alpha_y {
                let tan_theta2: Float = self.alpha_x * self.alpha_x * u[0] / (1.0 - u[0]);
                cos_theta = 1.0 / (1.0 + tan_theta2).sqrt();
            } else {
                phi = (self.alpha_y / self.alpha_x * (2.0 * PI * u[1] + 0.5 * PI).tan()).atan();
                if u[1] > 0.5 {
                    phi += PI;
                }
                let sin_phi: Float = phi.sin();
                let cos_phi: Float = phi.cos();
                let alphax2: Float = self.alpha_x * self.alpha_x;
                let alphay2: Float = self.alpha_y * self.alpha_y;
                let alpha2: Float =
                    1.0 / (cos_phi * cos_phi / alphax2 + sin_phi * sin_phi / alphay2);
                let tan_theta2: Float = alpha2 * u[0] / (1.0 - u[0]);
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
                wh = trowbridge_reitz_sample(&-(*wo), self.alpha_x, self.alpha_y, u[0], u[1]);
                wh = -wh;
            } else {
                wh = trowbridge_reitz_sample(wo, self.alpha_x, self.alpha_y, u[0], u[1]);
            }
        }
        wh
    }

    fn get_sample_visible_area(&self) -> bool {
        self.sample_visible_area
    }
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
        let phi: Float = 6.28318530718 * u2;
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
    let new_u2: Float;
    if u2 > 0.5 {
        s = 1.0;
        new_u2 = 2.0 * (u2 - 0.5);
    } else {
        s = -1.0;
        new_u2 = 2.0 * (0.5 - u2);
    }
    let z: Float = (new_u2 * (new_u2 * (new_u2 * 0.27385 - 0.73369) + 0.46341))
        / (new_u2 * (new_u2 * (new_u2 * 0.093073 + 0.309420) - 1.0) + 0.597999);
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
    }.normalize();

    // 2. simulate P22_{wi}(x_slope, y_slope, 1, 1)
    let mut slope_x: Float = 0.0;
    let mut slope_y: Float = 0.0;
    trowbridge_reitz_sample_11(cos_theta(&wi_stretched), u1, u2, &mut slope_x, &mut slope_y);

    // 3. rotate
    let tmp: Float = cos_phi(&wi_stretched) * slope_x - sin_phi(&wi_stretched) * slope_y;
    slope_y = sin_phi(&wi_stretched) * slope_x + cos_phi(&wi_stretched) * slope_y;
    slope_x = tmp;

    // 4. unstretch
    slope_x = alpha_x * slope_x;
    slope_y = alpha_y * slope_y;

    // 5. compute normal
    Vector3f {
        x: -slope_x,
        y: -slope_y,
        z: 1.0,
    }.normalize()
}
