// std
use std::sync::Arc;
// others
use atomic::{Atomic, Ordering};
// pbrt
use core::pbrt::Float;
use core::pbrt::{bits_to_float, float_to_bits};

// parallel.h

#[derive(Debug, Default, Clone)]
pub struct AtomicFloat {
    pub bits: Arc<Atomic<u32>>,
}

impl AtomicFloat {
    pub fn add(&self, v: Float) {
        let old_bits: u32 = self.bits.load(Ordering::SeqCst);
        self.bits
            .store(float_to_bits(bits_to_float(old_bits) + v), Ordering::SeqCst);
    }
}

impl From<AtomicFloat> for Float {
    fn from(a: AtomicFloat) -> Float {
        let bits: u32 = a.bits.load(Ordering::SeqCst);
        bits_to_float(bits) as Float
    }
}
