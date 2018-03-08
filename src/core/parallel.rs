// std
// use std;
// use std::sync::atomic::{AtomicU32, Ordering};
// pbrt
use core::pbrt::Float;
use core::pbrt::{bits_to_float, float_to_bits};

// parallel.h

#[derive(Debug,Default,Copy,Clone)]
pub struct AtomicFloat {
    bits: u32,
}

// impl AtomicFloat {
//     pub fn add(&self, v: Float) {
//         let old_bits: u32 = self.bits;
//         self.bits = float_to_bits(bits_to_float(old_bits) + v);
//         // loop {
//         // let new_bits: u32 = float_to_bits(bits_to_float(old_bits) + v);
//         //     let success = match self.bits.compare_exchange_weak(
//         //         old_bits,
//         //         new_bits,
//         //         Ordering::SeqCst,
//         //         Ordering::Relaxed,
//         //     ) {
//         //         Ok(_) => true,
//         //         Err(_) => false,
//         //     };
//         //     if success {
//         //         return;
//         //     }
//         // }
//     }
// }

impl From<AtomicFloat> for Float {
    fn from(a: AtomicFloat) -> Float {
        bits_to_float(a.bits) as Float
    }
}
