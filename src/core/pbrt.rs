//! Type definitions of Float and Spectrum, otherwise constants and
//! functions which can be used almost everywhere else in the code.

// std
use std;
use std::f32::consts::PI;
use std::ops::{Add, BitAnd, Div, Mul, Sub};
// others
use num;
// pbrt
use crate::core::spectrum::RGBSpectrum;

// see pbrt.h

pub type Spectrum = RGBSpectrum;

pub type Float = f32;

pub const MACHINE_EPSILON: Float = std::f32::EPSILON * 0.5;
pub const SHADOW_EPSILON: Float = 0.0001;
pub const INV_PI: Float = 0.318_309_886_183_790_671_54;
pub const INV_2_PI: Float = 0.159_154_943_091_895_335_77;
pub const INV_4_PI: Float = 0.079_577_471_545_947_667_88;
pub const PI_OVER_2: Float = 1.570_796_326_794_896_619_23;
pub const PI_OVER_4: Float = 0.785_398_163_397_448_309_61;
pub const SQRT_2: Float = 1.414_213_562_373_095_048_80;

/// Use **unsafe**
/// [std::mem::transmute_copy][transmute_copy]
/// to convert *f32* to *u32*.
///
/// [transmute_copy]: https://doc.rust-lang.org/std/mem/fn.transmute_copy.html
pub fn float_to_bits(f: f32) -> u32 {
    // uint64_t ui;
    // memcpy(&ui, &f, sizeof(double));
    // return ui;
    let rui: u32;
    unsafe {
        let ui: u32 = std::mem::transmute_copy(&f);
        rui = ui;
    }
    rui
}

/// Use **unsafe**
/// [std::mem::transmute_copy][transmute_copy]
/// to convert *u32* to *f32*.
///
/// [transmute_copy]: https://doc.rust-lang.org/std/mem/fn.transmute_copy.html
pub fn bits_to_float(ui: u32) -> f32 {
    // float f;
    // memcpy(&f, &ui, sizeof(uint32_t));
    // return f;
    let rf: f32;
    unsafe {
        let f: f32 = std::mem::transmute_copy(&ui);
        rf = f;
    }
    rf
}

/// Bump a floating-point value up to the next greater representable
/// floating-point value.
pub fn next_float_up(v: f32) -> f32 {
    if v.is_infinite() && v > 0.0 {
        v
    } else {
        let new_v = if v == -0.0 { 0.0 } else { v };
        let mut ui: u32 = float_to_bits(new_v);
        if new_v >= 0.0 {
            ui += 1;
        } else {
            ui -= 1;
        }
        bits_to_float(ui)
    }
}

/// Bump a floating-point value down to the next smaller representable
/// floating-point value.
pub fn next_float_down(v: f32) -> f32 {
    if v.is_infinite() && v < 0.0 {
        v
    } else {
        let new_v = if v == 0.0 { -0.0 } else { v };
        let mut ui: u32 = float_to_bits(new_v);
        if new_v > 0.0 {
            ui -= 1;
        } else {
            ui += 1;
        }
        bits_to_float(ui)
    }
}

/// Error propagation.
pub fn gamma(n: i32) -> Float {
    (n as Float * MACHINE_EPSILON) / (1.0 - n as Float * MACHINE_EPSILON)
}

/// Is used to write sRGB-compatible 8-bit image files.
pub fn gamma_correct(value: Float) -> Float {
    if value <= 0.003_130_8 {
        12.92 * value
    } else {
        1.055 as Float * value.powf((1.0 / 2.4) as Float) - 0.055
    }
}

/// Clamp the given value *val* to lie between the values *low* and *high*.
pub fn clamp_t<T>(val: T, low: T, high: T) -> T
where
    T: PartialOrd,
{
    let r: T;
    if val < low {
        r = low;
    } else if val > high {
        r = high;
    } else {
        r = val;
    }
    r
}

/// Computes the remainder of a/b. Provides the behavior that the
/// modulus of a negative number is always positive.
pub fn mod_t<T>(a: T, b: T) -> T
where
    T: num::Zero
        + Copy
        + PartialOrd
        + Add<T, Output = T>
        + Sub<T, Output = T>
        + Mul<T, Output = T>
        + Div<T, Output = T>,
{
    let result: T = a - (a / b) * b;
    if result < num::Zero::zero() {
        result + b
    } else {
        result
    }
}

/// Convert from angles expressed in degrees to radians.
pub fn radians(deg: Float) -> Float {
    (PI / 180.0) * deg
}

/// Convert from angles expressed in radians to degrees.
pub fn degrees(rad: Float) -> Float {
    (180.0 / PI) * rad
}

pub fn log_2(x: Float) -> Float {
    let inv_log2: Float = 1.442_695_040_888_963_387_004_650_940_071;
    x.ln() * inv_log2
}

/// Compute an integer base-2 logarithm function.
pub fn log_2_int_u32(v: u32) -> i32 {
    // C++: return 31 - __builtin_clz(v);
    31_i32 - v.leading_zeros() as i32
}

/// Compute an integer base-2 logarithm function.
pub fn log_2_int_i32(v: i32) -> i32 {
    log_2_int_u32(v as u32)
}

/// Compute an integer base-2 logarithm function.
pub fn log_2_int_u64(v: u64) -> i64 {
    // C++: return 63 - __builtin_clz(v);
    63_i64 - v.leading_zeros() as i64
}

/// Compute an integer base-2 logarithm function.
pub fn log_2_int_i64(v: i64) -> i64 {
    log_2_int_u64(v as u64)
}

/// Determine if a given integer is an exact power of 2.
pub fn is_power_of_2<T>(v: T) -> bool
where
    T: num::Zero + num::One + Copy + PartialOrd + BitAnd<T, Output = T> + Sub<T, Output = T>,
{
    // https://doc.rust-lang.org/std/primitive.u32.html#method.is_power_of_two
    (v > num::Zero::zero()) && !((v & (v - num::One::one())) > num::Zero::zero())
}

/// Round an integer up to the next higher (or equal) power of 2.
pub fn round_up_pow2_32(v: i32) -> i32 {
    let mut ret: i32 = v; // copy value
    ret -= 1_i32;
    ret |= ret >> 1;
    ret |= ret >> 2;
    ret |= ret >> 4;
    ret |= ret >> 8;
    ret |= ret >> 16;
    ret + 1
}

/// Round an integer up to the next higher (or equal) power of 2.
pub fn round_up_pow2_64(v: i64) -> i64 {
    let mut ret: i64 = v; // copy value
    ret -= 1_i64;
    ret |= ret >> 1;
    ret |= ret >> 2;
    ret |= ret >> 4;
    ret |= ret >> 8;
    ret |= ret >> 16;
    ret + 1
}

/// Helper function which emulates the behavior of std::upper_bound().
pub fn find_interval<P>(size: i32, pred: P) -> i32
where
    P: Fn(i32) -> bool,
{
    let mut first: i32 = 0;
    let mut len: i32 = size;
    while len > 0 {
        let half = len >> 1;
        let middle = first + half;
        // bisect range based on value of _pred_ at _middle_
        if pred(middle) {
            first = middle + 1;
            len -= half + 1;
        } else {
            len = half;
        }
    }
    clamp_t(first - 1, 0, size - 2)
}

/// Interpolate linearly between two provided values.
pub fn lerp<S, T>(t: S, a: T, b: T) -> T
where
    S: num::One,
    S: Sub<S, Output = S>,
    S: Copy,
    T: Add<T, Output = T>,
    T: Mul<S, Output = T>,
{
    let one: S = num::One::one();
    a * (one - t) + b * t
}

/// Find solution(s) of the quadratic equation at<sup>2</sup> + bt + c = 0.
pub fn quadratic(a: Float, b: Float, c: Float, t0: &mut Float, t1: &mut Float) -> bool {
    // find quadratic discriminant
    let discrim: f64 = (b as f64) * (b as f64) - 4.0 * (a as f64) * (c as f64);
    if discrim < 0.0 {
        false
    } else {
        let root_discrim: f64 = discrim.sqrt();
        // compute quadratic _t_ values
        let q = if b < 0.0 {
            -0.5 * (b as f64 - root_discrim)
        } else {
            -0.5 * (b as f64 + root_discrim)
        };
        *t0 = q as Float / a;
        *t1 = c / q as Float;
        if *t0 > *t1 {
            std::mem::swap(&mut (*t0), &mut (*t1))
        }
        true
    }
}

pub fn erf_inv(x: Float) -> Float {
    let clamped_x: Float = clamp_t(x, -0.99999, 0.99999);
    let mut w: Float = -((1.0 as Float - clamped_x) * (1.0 as Float + clamped_x)).ln();
    let mut p: Float;
    if w < 5.0 as Float {
        w -= 2.5 as Float;
        p = 2.810_226_36e-08;
        p = 3.432_739_39e-07 + p * w;
        p = -3.523_387_7e-06 + p * w;
        p = -4.391_506_54e-06 + p * w;
        p = 0.000_218_580_87 + p * w;
        p = -0.001_253_725_03 + p * w;
        p = -0.004_177_681_640 + p * w;
        p = 0.246_640_727 + p * w;
        p = 1.501_409_41 + p * w;
    } else {
        w = w.sqrt() - 3.0 as Float;
        p = -0.000_200_214_257;
        p = 0.000_100_950_558 + p * w;
        p = 0.001_349_343_22 + p * w;
        p = -0.003_673_428_44 + p * w;
        p = 0.005_739_507_73 + p * w;
        p = -0.007_622_461_3 + p * w;
        p = 0.009_438_870_47 + p * w;
        p = 1.001_674_06 + p * w;
        p = 2.832_976_82 + p * w;
    }
    p * clamped_x
}

pub fn erf(x: Float) -> Float {
    // constants
    let a1: Float = 0.254_829_592;
    let a2: Float = -0.284_496_736;
    let a3: Float = 1.421_413_741;
    let a4: Float = -1.453_152_027;
    let a5: Float = 1.061_405_429;
    let p: Float = 0.327_591_1;
    // save the sign of x
    let sign = if x < 0.0 as Float { -1.0 } else { 1.0 };
    let x: Float = x.abs();
    // A&S formula 7.1.26
    let t: Float = 1.0 as Float / (1.0 as Float + p * x);
    let y: Float =
        1.0 as Float - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * (-x * x).exp();
    sign * y
}
