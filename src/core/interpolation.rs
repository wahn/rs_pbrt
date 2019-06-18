//! Spline-based interpolation to reconstruct BSDF values (instead of
//! using large lookup tables).

// std
use std::f32::consts::PI;
// pbrt
use crate::core::pbrt::find_interval;
use crate::core::pbrt::Float;
use crate::core::pbrt::INV_2_PI;

/// Calculates an offset and four weights for Catmull-Rom spline
/// interpolation.
pub fn catmull_rom_weights(
    nodes: &Vec<Float>,
    x: Float,
    offset: &mut i32,
    weights: &mut [Float; 4],
) -> bool {
    // return _false_ if _x_ is out of bounds
    if !(x >= *nodes.first().unwrap() && x <= *nodes.last().unwrap()) {
        return false;
    }
    // search for the interval _idx_ containing _x_
    let idx: i32 = find_interval(nodes.len() as i32, |index| nodes[index as usize] <= x);
    *offset = idx - 1;
    assert!(idx >= 0);
    let x0: Float = nodes[idx as usize];
    let x1: Float = nodes[(idx + 1) as usize];
    // compute the $t$ parameter and powers
    let t: Float = (x - x0) / (x1 - x0);
    let t2: Float = t * t;
    let t3: Float = t2 * t;
    // compute initial node weights $w_1$ and $w_2$
    weights[1] = 2.0 as Float * t3 - 3.0 as Float * t2 + 1.0 as Float;
    weights[2] = -2.0 as Float * t3 + 3.0 as Float * t2;
    // compute first node weight $w_0$
    if idx > 0_i32 {
        let w0: Float = (t3 - 2.0 as Float * t2 + t) * (x1 - x0) / (x1 - nodes[(idx - 1) as usize]);
        weights[0] = -w0;
        weights[2] += w0;
    } else {
        let w0: Float = t3 - 2.0 as Float * t2 + t;
        weights[0] = 0.0 as Float;
        weights[1] -= w0;
        weights[2] += w0;
    }
    // compute last node weight $w_3$
    if (idx + 2) < nodes.len() as i32 {
        let w3: Float = (t3 - t2) * (x1 - x0) / (nodes[(idx + 2) as usize] - x0);
        weights[1] -= w3;
        weights[3] = w3;
    } else {
        let w3: Float = t3 - t2;
        weights[1] -= w3;
        weights[2] += w3;
        weights[3] = 0.0 as Float;
    }
    true
}

/// Importance sampling of 2D functions via spline interpolants.
pub fn sample_catmull_rom_2d(
    nodes1: &Vec<Float>,
    nodes2: &Vec<Float>,
    values: &Vec<Float>,
    cdf: &Vec<Float>,
    alpha: Float,
    u: Float,
    fval: Option<&mut Float>,
    pdf: Option<&mut Float>,
) -> Float {
    let size2: i32 = nodes2.len() as i32;
    // local copy
    let mut u: Float = u;
    // determine offset and coefficients for the _alpha_ parameter
    let mut offset: i32 = 0;
    let mut weights: [Float; 4] = [0.0 as Float; 4];
    if !catmull_rom_weights(nodes1, alpha, &mut offset, &mut weights) {
        return 0.0 as Float;
    }
    // define a lambda function to interpolate table entries
    let interpolate = |array: &Vec<Float>, idx: i32| -> Float {
        let mut value: Float = 0.0;
        for i in 0..4 {
            if weights[i] != 0.0 as Float {
                let index: i32 = (offset + i as i32) * size2 + idx;
                assert!(index >= 0);
                value += array[index as usize] * weights[i];
            }
        }
        value
    };
    // map _u_ to a spline interval by inverting the interpolated _cdf_
    let maximum: Float = interpolate(cdf, size2 - 1_i32);
    u *= maximum;
    let idx: i32 = find_interval(size2, |index| interpolate(cdf, index) <= u);
    // look up node positions and interpolated function values
    let f0: Float = interpolate(values, idx);
    let f1: Float = interpolate(values, idx + 1);
    assert!(idx >= 0);
    let x0: Float = nodes2[idx as usize];
    let x1: Float = nodes2[(idx + 1) as usize];
    let width: Float = x1 - x0;
    // re-scale _u_ using the interpolated _cdf_
    u = (u - interpolate(cdf, idx)) / width;
    // approximate derivatives using finite differences of the interpolant
    let d0: Float;
    let d1: Float;
    if idx > 0_i32 {
        d0 = width * (f1 - interpolate(values, idx - 1)) / (x1 - nodes2[(idx - 1) as usize]);
    } else {
        d0 = f1 - f0;
    }
    if (idx + 2) < size2 {
        d1 = width * (interpolate(values, idx + 2) - f0) / (nodes2[(idx + 2) as usize] - x0);
    } else {
        d1 = f1 - f0;
    }

    // invert definite integral over spline segment and return solution

    // set initial guess for $t$ by importance sampling a linear interpolant
    let mut t: Float;
    if f0 != f1 {
        t = (f0
            - (0.0 as Float)
                .max(f0 * f0 + 2.0 as Float * u * (f1 - f0))
                .sqrt())
            / (f0 - f1);
    } else {
        t = u / f0;
    }
    let mut a: Float = 0.0;
    let mut b: Float = 1.0;
    let mut f_hat;
    let mut fhat;
    loop {
        // fall back to a bisection step when _t_ is out of bounds
        if !(t >= a && t <= b) {
            t = 0.5 as Float * (a + b);
        }
        // evaluate target function and its derivative in Horner form
        f_hat = t
            * (f0
                + t * (0.5 as Float * d0
                    + t * ((1.0 as Float / 3.0 as Float) * (-2.0 as Float * d0 - d1) + f1 - f0
                        + t * (0.25 as Float * (d0 + d1) + 0.5 as Float * (f0 - f1)))));
        fhat = f0
            + t * (d0
                + t * (-2.0 as Float * d0 - d1
                    + 3.0 as Float * (f1 - f0)
                    + t * (d0 + d1 + 2.0 as Float * (f0 - f1))));
        // stop the iteration if converged
        if (f_hat - u).abs() < 1e-6 as Float || b - a < 1e-6 as Float {
            break;
        }
        // update bisection bounds using updated _t_
        if (f_hat - u) < 0.0 as Float {
            a = t;
        } else {
            b = t;
        }
        // perform a Newton step
        t -= (f_hat - u) / fhat;
    }
    // return the sample position and function value
    if let Some(fval) = fval {
        *fval = fhat;
    }
    if let Some(pdf) = pdf {
        *pdf = fhat / maximum;
    }
    x0 + width * t
}

pub fn integrate_catmull_rom(
    n: i32,
    x: &Vec<Float>,
    offset: usize,
    values: &Vec<Float>,
    cdf: &mut Vec<Float>,
) -> Float {
    let mut sum: Float = 0.0;
    cdf[offset + 0] = 0.0 as Float;
    for i in 0..(n - 1) as usize {
        // look up $x_i$ and function values of spline segment _i_
        let x0: Float = x[i];
        let x1: Float = x[i + 1];
        let f0: Float = values[offset + i];
        let f1: Float = values[offset + i + 1];
        let width: Float = x1 - x0;
        // approximate derivatives using finite differences
        let d0: Float;
        let d1: Float;
        if i > 0 {
            d0 = width * (f1 - values[offset + i - 1]) / (x1 - x[i - 1]);
        } else {
            d0 = f1 - f0;
        }
        if i + 2 < n as usize {
            d1 = width * (values[offset + i + 2] - f0) / (x[i + 2] - x0);
        } else {
            d1 = f1 - f0;
        }
        // keep a running sum and build a cumulative distribution function
        sum += ((d0 - d1) * (1.0 as Float / 12.0 as Float) + (f0 + f1) * 0.5 as Float) * width;
        cdf[offset + i + 1] = sum;
    }
    sum
}

/// Evaluates the weighted sum of cosines.
pub fn fourier(a: &Vec<Float>, si: usize, m: i32, cos_phi: f64) -> Float {
    let mut value: f64 = 0.0;
    // initialize cosine iterates
    let mut cos_k_minus_one_phi: f64 = cos_phi;
    let mut cos_k_phi: f64 = 1.0;
    for k in 0..m as usize {
        // add the current summand and update the cosine iterates
        value += a[si + k] as f64 * cos_k_phi;
        let cos_k_plus_one_phi: f64 = 2.0 as f64 * cos_phi * cos_k_phi - cos_k_minus_one_phi;
        cos_k_minus_one_phi = cos_k_phi;
        cos_k_phi = cos_k_plus_one_phi;
    }
    value as Float
}

/// Returns the value of the Fourier expansion at the sampled
/// position.
pub fn sample_fourier(
    ak: &Vec<Float>,
    recip: &Vec<Float>,
    m: i32,
    u: Float,
    pdf: &mut Float,
    phi_ptr: &mut Float,
) -> Float {
    // local copy
    let mut u: Float = u;
    // pick a side and declare bisection variables
    let flip: bool;
    if u >= 0.5 as Float {
        flip = true;
        u = 1.0 as Float - 2.0 as Float * (u - 0.5 as Float);
    } else {
        flip = false;
        u *= 2.0 as Float;
    }
    let mut a: f64 = 0.0;
    let mut b: f64 = PI as f64;
    let mut phi: f64 = (0.5_f32 * PI) as f64;
    let mut cf: f64;
    let mut f: f64;
    loop {
        // evaluate $cf(\phi)$ and its derivative $f(\phi)$

        // initialize sine and cosine iterates
        let cos_phi: f64 = phi.cos();
        let sin_phi: f64 = ((0.0_f64).max(1.0_f64 - cos_phi * cos_phi)).sqrt();
        let mut cos_phi_prev: f64 = cos_phi;
        let mut cos_phi_cur: f64 = 1.0;
        let mut sin_phi_prev = -sin_phi;
        let mut sin_phi_cur: f64 = 0.0;
        // initialize _cf_ and _f_ with the first series term
        cf = ak[0] as f64 * phi;
        f = ak[0] as f64;
        for k in 1..m as usize {
            // compute next sine and cosine iterates
            let sin_phi_next: f64 = (2.0_f64 * cos_phi) * sin_phi_cur - sin_phi_prev;
            let cos_phi_next: f64 = (2.0_f64 * cos_phi) * cos_phi_cur - cos_phi_prev;
            sin_phi_prev = sin_phi_cur;
            sin_phi_cur = sin_phi_next;
            cos_phi_prev = cos_phi_cur;
            cos_phi_cur = cos_phi_next;
            // add the next series term to _cf_ and _f_
            cf += (ak[k] * recip[k]) as f64 * sin_phi_next;
            f += ak[k] as f64 * cos_phi_next;
        }
        cf -= (u * ak[0] * PI) as f64;
        // update bisection bounds using updated $\phi$
        if cf > 0.0 as f64 {
            b = phi;
        } else {
            a = phi;
        }
        // stop the Fourier bisection iteration if converged
        if cf.abs() < 1e-6 as f64 || b - a < 1e-6 as f64 {
            break;
        }
        // perform a Newton step given $f(\phi)$ and $cf(\phi)$
        phi -= cf / f;
        // fall back to a bisection step when $\phi$ is out of bounds
        if !(phi > a && phi < b) {
            phi = 0.5 as f64 * (a + b);
        }
    }
    // potentially flip $\phi$ and return the result
    if flip {
        phi = 2.0 as f64 * PI as f64 - phi;
    }
    *pdf = (INV_2_PI as f64 * f / ak[0] as f64) as Float;
    *phi_ptr = phi as Float;
    f as Float
}
