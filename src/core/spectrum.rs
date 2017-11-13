// std
use std::ops::{Add, AddAssign, Div, Index, IndexMut, Mul, MulAssign, Sub};
// others
use num::Zero;
// pbrt
use core::pbrt::{Float, Spectrum};
use core::pbrt::clamp_t;

// see spectrum.h

// #[derive(Debug,Clone)]
// pub enum SpectrumType {
//     Reflectance,
//     Illuminant,
// }

#[derive(Debug,Default,Copy,Clone)]
pub struct RGBSpectrum {
    pub c: [Float; 3],
}

impl RGBSpectrum {
    pub fn new(v: Float) -> Self {
        // let n_spectrum_samples = 3; // RGB
        RGBSpectrum { c: [v, v, v] }
        // TODO: DCHECK(!HasNaNs());
    }
    pub fn rgb(r: Float, g: Float, b: Float) -> RGBSpectrum {
        RGBSpectrum { c: [r, g, b] }
    }
    pub fn from_rgb(rgb: &[Float; 3]) -> RGBSpectrum {
        let mut s: RGBSpectrum = RGBSpectrum::new(0.0 as Float);
        s.c[0] = rgb[0];
        s.c[1] = rgb[1];
        s.c[2] = rgb[2];
        // TODO: DCHECK(!s.HasNaNs());
        s
    }
    pub fn to_rgb(&self, rgb: &mut [Float; 3]) {
        rgb[0] = self.c[0];
        rgb[1] = self.c[1];
        rgb[2] = self.c[2];
    }
    pub fn to_xyz(&self, xyz: &mut [Float; 3]) {
        rgb_to_xyz(&self.c, xyz);
    }
    pub fn from_srgb(rgb: &[u8; 3]) -> RGBSpectrum {
        fn convert(v: u8) -> Float {
            let value = v as Float / 255.0;
            // see InverseGammaCorrect(Float value) in pbrt.h
            if value <= 0.04045 {
                value / 12.92
            } else {
                ((value + 0.055) * 1.0 / 1.055).powf(2.4)
            }
        }
        RGBSpectrum::rgb(convert(rgb[0]), convert(rgb[1]), convert(rgb[2]))
    }
    pub fn y(&self) -> Float {
        let y_weight: [Float; 3] = [0.212671, 0.715160, 0.072169];
        y_weight[0] * self.c[0] + y_weight[1] * self.c[1] + y_weight[2] * self.c[2]
    }
    // from CoefficientSpectrum
    pub fn is_black(&self) -> bool {
        for i in 0..3 {
            if self.c[i] != 0.0 as Float {
                return false;
            }
        }
        true
    }
    pub fn clamp(&self, low: Float, high: Float) -> RGBSpectrum {
        let mut ret: RGBSpectrum = RGBSpectrum::default();
        let n_spectrum_samples: usize = 3; // RGB
        for i in 0..n_spectrum_samples {
            ret.c[i] = clamp_t(self.c[i], low, high);
        }
        assert!(!ret.has_nans());
        ret
    }
    pub fn max_component_value(&self) -> Float {
        let mut m: Float = self.c[0];
        let n_spectrum_samples: usize = 3; // RGB
        for i in 1..n_spectrum_samples {
            m = m.max(self.c[i]);
        }
        m
    }
    pub fn has_nans(&self) -> bool {
        for i in 0..3 {
            if self.c[i].is_nan() {
                return true;
            }
        }
        false
    }
}

impl PartialEq for RGBSpectrum {
    fn eq(&self, rhs: &RGBSpectrum) -> bool {
        for i in 0..3 {
            if self.c[i] != rhs.c[i] {
                return false;
            }
        }
        true
    }
}

impl Add for RGBSpectrum {
    type Output = RGBSpectrum;
    fn add(self, rhs: RGBSpectrum) -> RGBSpectrum {
        RGBSpectrum {
            c: [self.c[0] + rhs.c[0],
                self.c[1] + rhs.c[1],
                self.c[2] + rhs.c[2]],
        }
    }
}

impl AddAssign for RGBSpectrum {
    fn add_assign(&mut self, rhs: RGBSpectrum) {
        // TODO: DCHECK(!s2.HasNaNs());
        self.c[0] += rhs.c[0];
        self.c[1] += rhs.c[1];
        self.c[2] += rhs.c[2];
    }
}

impl Mul for RGBSpectrum {
    type Output = RGBSpectrum;
    fn mul(self, rhs: RGBSpectrum) -> RGBSpectrum {
        RGBSpectrum {
            c: [self.c[0] * rhs.c[0],
                self.c[1] * rhs.c[1],
                self.c[2] * rhs.c[2]],
        }
    }
}

impl Mul<Float> for RGBSpectrum {
    type Output = RGBSpectrum;
    fn mul(self, rhs: Float) -> RGBSpectrum {
        RGBSpectrum { c: [self.c[0] * rhs, self.c[1] * rhs, self.c[2] * rhs] }
    }
}

impl MulAssign for RGBSpectrum {
    fn mul_assign(&mut self, rhs: RGBSpectrum) {
        // TODO: DCHECK(!HasNaNs());
        self.c[0] *= rhs.c[0];
        self.c[1] *= rhs.c[1];
        self.c[2] *= rhs.c[2];
    }
}

impl Sub for RGBSpectrum {
    type Output = RGBSpectrum;
    fn sub(self, rhs: RGBSpectrum) -> RGBSpectrum {
        RGBSpectrum {
            c: [self.c[0] - rhs.c[0],
                self.c[1] - rhs.c[1],
                self.c[2] - rhs.c[2]],
        }
    }
}

impl Div for RGBSpectrum {
    type Output = RGBSpectrum;
    fn div(self, rhs: RGBSpectrum) -> RGBSpectrum {
        RGBSpectrum {
            c: [self.c[0] / rhs.c[0],
                self.c[1] / rhs.c[1],
                self.c[2] / rhs.c[2]],
        }
    }
}

impl Div<Float> for RGBSpectrum {
    type Output = RGBSpectrum;
    fn div(self, rhs: Float) -> RGBSpectrum {
        assert_ne!(rhs, 0.0 as Float);
        assert!(!rhs.is_nan(), "rhs is NaN");
        let ret: RGBSpectrum =
            RGBSpectrum { c: [self.c[0] / rhs, self.c[1] / rhs, self.c[2] / rhs] };
        assert!(!ret.has_nans());
        ret
    }
}

impl Zero for RGBSpectrum {
    fn zero() -> RGBSpectrum {
        RGBSpectrum::new(0.0 as Float)
    }

    fn is_zero(&self) -> bool {
        self.is_black()
    }
}

impl Index<usize> for RGBSpectrum {
    type Output = Float;
    fn index(&self, index: usize) -> &Float {
        match index {
            0 => &self.c[0],
            1 => &self.c[1],
            2 => &self.c[2],
            _ => panic!("Check failed: i >= 0 && i <= 2"),
        }
    }
}

impl IndexMut<usize> for RGBSpectrum {
    fn index_mut(&mut self, index: usize) -> &mut Float {
        match index {
            0 => &mut self.c[0],
            1 => &mut self.c[1],
            2 => &mut self.c[2],
            _ => panic!("Check failed: i >= 0 && i <= 2"),
        }
    }
}

impl From<Float> for RGBSpectrum {
    fn from(f: Float) -> Self {
        RGBSpectrum::new(f)
    }
}

/// Calculate RGB coefficients from a XYZ representation.
pub fn xyz_to_rgb(xyz: &[Float; 3], rgb: &mut [Float; 3]) {
    rgb[0] = 3.240479 * xyz[0] - 1.537150 * xyz[1] - 0.498535 * xyz[2];
    rgb[1] = -0.969256 * xyz[0] + 1.875991 * xyz[1] + 0.041556 * xyz[2];
    rgb[2] = 0.055648 * xyz[0] - 0.204043 * xyz[1] + 1.057311 * xyz[2];
}

/// Calculate XYZ representation from RGB coefficients.
pub fn rgb_to_xyz(rgb: &[Float; 3], xyz: &mut [Float; 3]) {
    xyz[0] = 0.412453 * rgb[0] + 0.357580 * rgb[1] + 0.180423 * rgb[2];
    xyz[1] = 0.212671 * rgb[0] + 0.715160 * rgb[1] + 0.072169 * rgb[2];
    xyz[2] = 0.019334 * rgb[0] + 0.119193 * rgb[1] + 0.950227 * rgb[2];
}

/// Interpolate linearly between two provided values.
pub fn lerp_rgb(t: Float, s1: Spectrum, s2: Spectrum) -> Spectrum {
    s1 * (1.0 as Float - t) + s2 * t
}
