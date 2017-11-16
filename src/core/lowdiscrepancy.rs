// pbrt
use core::pbrt::Float;
use core::rng::ONE_MINUS_EPSILON;
use core::rng::Rng;
use core::sampling::shuffle;
use geometry::{Point2i, Point2f};

// see lowdiscrepancy.h

/// The bits of an integer quantity can be efficiently reversed with a
/// series of logical bit operations.
pub fn reverse_bits_32(n: u32) -> u32 {
    let mut n = (n << 16) | (n >> 16);
    n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8);
    n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4);
    n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2);
    n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1);
    n
}

/// The bits of a 64-bit value can be reversed by reversing the two
/// 32-bit components individually and then interchanging them.
pub fn reverse_bits_64(n: u64) -> u64 {
    let n0: u64 = reverse_bits_32(n as u32) as u64;
    let n1: u64 = reverse_bits_32((n >> 32) as u32) as u64;
    (n0 << 32) | n1
}


/// Compute the inverse of the radical inverse function.
pub fn inverse_radical_inverse(base: u8, inverse: u64, n_digits: u64) -> u64 {
    let mut inverse: u64 = inverse;
    let mut index: u64 = 0_64;
    for _i in 0..n_digits {
        let digit: u64 = inverse % base as u64;
        inverse /= base as u64;
        index = index * base as u64 + digit;
    }
    index
}

/// Takes a generator matrix *c*, a number of 1D samples to generate
/// *n*, and stores the corresponding samples in memory at the
/// location pointed to by *p*.
pub fn gray_code_sample_1d(c: [u32; 32], n: u32, scramble: u32, p: &mut [Float]) {
    let mut v: u32 = scramble;
    for i in 0..n as usize {
        // 1/2^32
        p[i] = (v as Float * 2.3283064365386963e-10 as Float).min(ONE_MINUS_EPSILON);
        v ^= c[(i + 1).trailing_zeros() as usize];
    }
}

/// Takes two generator matrices *c0* and *c1*, a number of 2D samples
/// to generate *n*, and stores the corresponding samples in memory at
/// the location pointed to by *p*.
pub fn gray_code_sample_2d(c0: &[u32], c1: &[u32], n: u32, scramble: &Point2i, p: &mut [Point2f]) {
    let mut v: [u32; 2] = [scramble.x as u32, scramble.y as u32];
    for i in 0..n as usize {
        p[i].x = (v[0] as Float * 2.3283064365386963e-10 as Float).min(ONE_MINUS_EPSILON);
        p[i].y = (v[1] as Float * 2.3283064365386963e-10 as Float).min(ONE_MINUS_EPSILON);
        v[0] ^= c0[(i + 1).trailing_zeros() as usize];
        v[1] ^= c1[(i + 1).trailing_zeros() as usize];
    }
}

/// Generates a number of scrambled 1D sample values using the Gray
/// code-based sampling machinery.
pub fn van_der_corput(n_samples_per_pixel_sample: i32,
                      n_pixel_samples: i32,
                      samples: &mut [Float],
                      rng: &mut Rng) {
    let scramble: u32 = rng.uniform_uint32();
    let c_van_der_corput: [u32; 32] = [0x80000000, 0x40000000, 0x20000000, 0x10000000, 0x8000000,
                                       0x4000000, 0x2000000, 0x1000000, 0x800000, 0x400000,
                                       0x200000, 0x100000, 0x80000, 0x40000, 0x20000, 0x10000,
                                       0x8000, 0x4000, 0x2000, 0x1000, 0x800, 0x400, 0x200, 0x100,
                                       0x80, 0x40, 0x20, 0x10, 0x8, 0x4, 0x2, 0x1];
    let total_samples: i32 = n_samples_per_pixel_sample * n_pixel_samples;
    gray_code_sample_1d(c_van_der_corput, total_samples as u32, scramble, samples);
    // randomly shuffle 1D sample points
    for i in 0..n_pixel_samples as usize {
        shuffle(&mut samples[(i * n_samples_per_pixel_sample as usize)..],
                n_samples_per_pixel_sample,
                1,
                rng);
    }
    shuffle(&mut samples[..],
            n_pixel_samples,
            n_samples_per_pixel_sample,
            rng);
}

/// Similar to *van_der_corput()*, but uses two generator matrices to
/// generate the first two dimensions of Sobol' points.
pub fn sobol_2d(n_samples_per_pixel_sample: i32,
                n_pixel_samples: i32,
                samples: &mut [Point2f],
                rng: &mut Rng) {
    let x: i32 = rng.uniform_uint32() as i32;
    let y: i32 = rng.uniform_uint32() as i32;
    let scramble: Point2i = Point2i { x: x, y: y };
    // define 2D Sobol$'$ generator matrices _c_sobol[2]_
    let c_sobol: [[u32; 32]; 2] =
        [[0x80000000_u32,
          0x40000000,
          0x20000000,
          0x10000000,
          0x8000000,
          0x4000000,
          0x2000000,
          0x1000000,
          0x800000,
          0x400000,
          0x200000,
          0x100000,
          0x80000,
          0x40000,
          0x20000,
          0x10000,
          0x8000,
          0x4000,
          0x2000,
          0x1000,
          0x800,
          0x400,
          0x200,
          0x100,
          0x80,
          0x40,
          0x20,
          0x10,
          0x8,
          0x4,
          0x2,
          0x1],
         [0x80000000, 0xc0000000, 0xa0000000, 0xf0000000, 0x88000000, 0xcc000000, 0xaa000000,
          0xff000000, 0x80800000, 0xc0c00000, 0xa0a00000, 0xf0f00000, 0x88880000, 0xcccc0000,
          0xaaaa0000, 0xffff0000, 0x80008000, 0xc000c000, 0xa000a000, 0xf000f000, 0x88008800,
          0xcc00cc00, 0xaa00aa00, 0xff00ff00, 0x80808080, 0xc0c0c0c0, 0xa0a0a0a0, 0xf0f0f0f0,
          0x88888888, 0xcccccccc, 0xaaaaaaaa, 0xffffffff]];
    gray_code_sample_2d(&c_sobol[0],
                        &c_sobol[1],
                        (n_samples_per_pixel_sample * n_pixel_samples) as u32,
                        &scramble,
                        &mut samples[..]);
    for i in 0..n_pixel_samples as usize {
        shuffle(&mut samples[(i * n_samples_per_pixel_sample as usize)..],
                n_samples_per_pixel_sample,
                1,
                rng);
    }
    shuffle(&mut samples[..],
            n_pixel_samples,
            n_samples_per_pixel_sample,
            rng);
}

// see lowdiscrepancy.cpp

/// Once we have an appropriate prime number, use it to compute the
/// radical inverse.
pub fn radical_inverse_specialized(base: u16, a: u64) -> Float {
    let inv_base: Float = 1.0 as Float / base as Float;
    let mut reversed_digits: u64 = 0_u64;
    let mut inv_base_n: Float = 1.0 as Float;
    let mut a: u64 = a; // shadowing input parameter
    while a != 0_u64 {
        let next: u64 = a / base as u64;
        let digit: u64 = a - next * base as u64;
        reversed_digits = reversed_digits * base as u64 + digit;
        inv_base_n *= inv_base;
        a = next;
    }
    assert!(reversed_digits as Float * inv_base_n < 1.00001 as Float);
    (reversed_digits as Float * inv_base_n).min(ONE_MINUS_EPSILON)
}

/// Map to an appropriate prime number and delegate to another
/// function to compute the radical inverse.
pub fn radical_inverse(base_index: u16, a: u64) -> Float {
    match base_index {
        0 => {
            // TODO: #ifndef PBRT_HAVE_HEX_FP_CONSTANTS
            // 0x1p-64 = (2.0 as Float).powi(-64 as i32)
            return reverse_bits_64(a) as Float * (2.0 as Float).powi(-64 as i32);
        }
        1 => {
            return radical_inverse_specialized(3_u16, a);
        }
        2 => {
            return radical_inverse_specialized(5_u16, a);
        }
        3 => {
            return radical_inverse_specialized(7_u16, a);
        }
        4 => {
            return radical_inverse_specialized(11_u16, a);
        }
        5 => {
            return radical_inverse_specialized(13_u16, a);
        }
        6 => {
            return radical_inverse_specialized(17_u16, a);
        }
        // WORK
        _ => {
            panic!("TODO: radical_inverse({:?}, {:?})", base_index, a);
        }
    };
}
