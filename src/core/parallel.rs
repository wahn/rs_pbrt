//! Using atomic operations on floating-point values. One example is
//! splatting pixel contributions.

// others
use atomic::{Atomic, Ordering};
// pbrt
use core::pbrt::Float;
use core::pbrt::{bits_to_float, float_to_bits};

// parallel.h

#[derive(Debug)]
pub struct AtomicFloat {
    pub bits: Atomic<u32>,
}

impl AtomicFloat {
    pub fn new(v: Float) -> AtomicFloat {
        AtomicFloat {
            bits: Atomic::new(float_to_bits(v)),
        }
    }
    pub fn add(&self, v: Float) {
        let mut old_bits: u32 = self.bits.load(Ordering::Relaxed);
        loop {
            let f: Float = bits_to_float(old_bits);
            let new_bits: u32 = float_to_bits(f + v);
            match self.bits.compare_exchange_weak(
                old_bits,
                new_bits,
                Ordering::SeqCst,
                Ordering::Relaxed,
            ) {
                Ok(_) => {
                    break;
                }
                Err(x) => {
                    old_bits = x;
                }
            }
        }
    }
}

impl Clone for AtomicFloat {
    fn clone(&self) -> Self {
        let bits: u32 = self.bits.load(Ordering::SeqCst);
        AtomicFloat {
            bits: Atomic::new(bits),
        }
    }
}

impl Default for AtomicFloat {
    fn default() -> AtomicFloat {
        AtomicFloat::new(0.0 as Float)
    }
}

impl<'a> From<&'a AtomicFloat> for Float {
    fn from(a: &'a AtomicFloat) -> Float {
        let bits: u32 = a.bits.load(Ordering::SeqCst);
        bits_to_float(bits) as Float
    }
}
