// std
use std;
use std::f32::consts::PI;
use std::ops::{Add, BitAnd, Div, Mul, Sub};
// others
use num;
// pbrt
use core::spectrum::RGBSpectrum;

// see pbrt.h

pub type Spectrum = RGBSpectrum;

pub type Float = f32;

pub const MACHINE_EPSILON: Float = std::f32::EPSILON * 0.5;
pub const SHADOW_EPSILON: Float = 0.0001;
pub const INV_PI: Float = 0.31830988618379067154;
pub const INV_2_PI: Float = 0.15915494309189533577;
pub const PI_OVER_2: Float = 1.57079632679489661923;
pub const PI_OVER_4: Float = 0.78539816339744830961;

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
        let new_v: f32;
        if v == -0.0 {
            new_v = 0.0;
        } else {
            new_v = v;
        }
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
        let new_v: f32;
        if v == 0.0 {
            new_v = -0.0;
        } else {
            new_v = v;
        }
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
    if value <= 0.0031308 {
        12.92 * value
    } else {
        1.055 as Float * value.powf((1.0 / 2.4) as Float) - 0.055
    }
}

/// Clamp the given value *val* to lie between the values *low* and *high*.
pub fn clamp_t<T>(val: T, low: T, high: T) -> T
    where T: PartialOrd
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
    where T: num::Zero + Copy + PartialOrd +
    Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T> + Div<T, Output = T>
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

/// Determine if a given integer is an exact power of 2.
pub fn is_power_of_2<T>(v: T) -> bool
    where T: num::Zero + num::One + Copy + PartialOrd + BitAnd<T, Output = T> + Sub<T, Output = T>
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

/// Interpolate linearly between two provided values.
pub fn lerp(t: Float, v1: Float, v2: Float) -> Float {
    (1.0 as Float - t) * v1 + t * v2
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
        let q: f64;
        if b < 0.0 {
            q = -0.5 * (b as f64 - root_discrim);
        } else {
            q = -0.5 * (b as f64 + root_discrim);
        }
        *t0 = q as Float / a;
        *t1 = c / q as Float;
        if *t0 > *t1 {
            // std::swap(*t0, *t1);
            let swap = *t0;
            *t0 = *t1;
            *t1 = swap;
        }
        true
    }
}
