// pbrt
use core::pbrt::Float;
use core::pbrt::bits_to_float;

// parallel.h

#[derive(Debug,Default,Copy,Clone)]
pub struct AtomicFloat {
    bits: u32,
}

impl From<AtomicFloat> for Float {
    fn from(a: AtomicFloat) -> Float {
        bits_to_float(a.bits) as Float
    }
}
