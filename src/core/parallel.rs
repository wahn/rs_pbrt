// std
use std::sync::Arc;
// others
use atomic::{Atomic, Ordering};
// pbrt
use core::pbrt::Float;
use core::pbrt::{bits_to_float, float_to_bits};

// parallel.h

#[derive(Debug, Clone)]
pub struct AtomicFloat {
    pub bits: Arc<Atomic<u32>>,
}

impl AtomicFloat {
    pub fn new(v: Float) -> AtomicFloat{
        AtomicFloat {
            bits: Arc::new(Atomic::new(float_to_bits(v))),
        }
    }
    pub fn add(&self, v: Float) {
        loop {
            let old_bits: u32 = self.bits.load(Ordering::SeqCst);
            let new_bits: u32 = float_to_bits(bits_to_float(old_bits) + v);
            let result = self.bits.compare_exchange_weak(
                old_bits,
                new_bits,
                Ordering::SeqCst,
                Ordering::Relaxed,
            );
            if result.is_ok() {
                break;
            }
        }
    }
}

impl Default for AtomicFloat {
    fn default() -> AtomicFloat {
        AtomicFloat::new(0.0 as Float)
    }
}

impl From<AtomicFloat> for Float {
    fn from(a: AtomicFloat) -> Float {
        let bits: u32 = a.bits.load(Ordering::SeqCst);
        bits_to_float(bits) as Float
    }
}
