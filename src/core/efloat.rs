//! In addition to working out the error bounds algebraically, we can
//! also have the computer do this work for us as some computation is
//! being performed. This approach is known as *running error
//! analysis*.

// std
use std;
use std::ops::{Add, Div, Mul, Sub};
// pbrt
use core::pbrt::MACHINE_EPSILON;
use core::pbrt::{next_float_down, next_float_up};

// see efloat.h

/// Find solution(s) of the quadratic equation at<sup>2</sup> + bt + c = 0 using
/// *EFloat* instead of *Float* for error bounds.
pub fn quadratic_efloat(a: EFloat, b: EFloat, c: EFloat, t0: &mut EFloat, t1: &mut EFloat) -> bool {
    let discrim: f64 = b.v as f64 * b.v as f64 - 4.0f64 * a.v as f64 * c.v as f64;
    if discrim < 0.0 {
        false
    } else {
        let root_discrim: f64 = discrim.sqrt();
        let float_root_discrim: EFloat = EFloat::new(
            root_discrim as f32,
            MACHINE_EPSILON as f32 * root_discrim as f32,
        );
        // compute quadratic _t_ values
        let q: EFloat;
        if b.v < 0.0f32 {
            q = (b - float_root_discrim) * -0.5f32;
        } else {
            q = (b + float_root_discrim) * -0.5f32;
        }
        *t0 = q / a;
        *t1 = c / q;
        if (*t0).v > (*t1).v {
            let swap: EFloat = *t0;
            *t0 = *t1;
            *t1 = swap;
        }
        true
    }
}

/// **EFloat** keeps track of an interval that describes the
/// uncertainty of a value of interest.
#[derive(Debug, Default, Copy, Clone)]
pub struct EFloat {
    pub v: f32,
    pub low: f32,
    pub high: f32,
}

impl EFloat {
    pub fn new(v: f32, err: f32) -> Self {
        if err == 0.0 {
            EFloat {
                v: v,
                low: v,
                high: v,
            }
        } else {
            EFloat {
                v: v,
                low: next_float_down(v - err),
                high: next_float_up(v + err),
            }
        }
    }
    pub fn lower_bound(&self) -> f32 {
        self.low
    }
    pub fn upper_bound(&self) -> f32 {
        self.high
    }
}

impl PartialEq for EFloat {
    fn eq(&self, rhs: &EFloat) -> bool {
        self.v == rhs.v
    }
}

impl Add for EFloat {
    type Output = EFloat;
    fn add(self, rhs: EFloat) -> EFloat {
        // TODO: r.Check();
        EFloat {
            v: self.v + rhs.v,
            low: next_float_down(self.lower_bound() + rhs.lower_bound()),
            high: next_float_up(self.upper_bound() + rhs.upper_bound()),
        }
    }
}

impl Sub for EFloat {
    type Output = EFloat;
    fn sub(self, rhs: EFloat) -> EFloat {
        // TODO: r.Check();
        EFloat {
            v: self.v - rhs.v,
            low: next_float_down(self.lower_bound() - rhs.upper_bound()),
            high: next_float_up(self.upper_bound() - rhs.lower_bound()),
        }
    }
}

impl Mul for EFloat {
    type Output = EFloat;
    fn mul(self, rhs: EFloat) -> EFloat {
        let prod: [f32; 4] = [
            self.lower_bound() * rhs.lower_bound(),
            self.upper_bound() * rhs.lower_bound(),
            self.lower_bound() * rhs.upper_bound(),
            self.upper_bound() * rhs.upper_bound(),
        ];
        // TODO: r.Check();
        EFloat {
            v: self.v * rhs.v,
            low: next_float_down(prod[0].min(prod[1]).min(prod[2].min(prod[3]))),
            high: next_float_up(prod[0].max(prod[1]).max(prod[2].max(prod[3]))),
        }
    }
}

impl Mul<f32> for EFloat {
    type Output = EFloat;
    fn mul(self, rhs: f32) -> EFloat {
        EFloat::new(rhs, 0.0) * self
    }
}

impl Div for EFloat {
    type Output = EFloat;
    fn div(self, rhs: EFloat) -> EFloat {
        let div: [f32; 4] = [
            self.lower_bound() / rhs.lower_bound(),
            self.upper_bound() / rhs.lower_bound(),
            self.lower_bound() / rhs.upper_bound(),
            self.upper_bound() / rhs.upper_bound(),
        ];
        // TODO: r.Check();
        if rhs.low < 0.0 && rhs.high > 0.0 {
            // the interval we're dividing by straddles zero, so just
            // return an interval of everything
            EFloat {
                v: self.v / rhs.v,
                low: -std::f32::INFINITY,
                high: std::f32::INFINITY,
            }
        } else {
            EFloat {
                v: self.v / rhs.v,
                low: next_float_down(div[0].min(div[1]).min(div[2].min(div[3]))),
                high: next_float_up(div[0].max(div[1]).max(div[2].max(div[3]))),
            }
        }
    }
}
