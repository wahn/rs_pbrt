// pbrt
use core::pbrt::Float;

// see rng.h

pub const FLOAT_ONE_MINUS_EPSILON: Float = 0.99999994;
pub const PCG32_DEFAULT_STATE: u64 = 0x853c49e6748fea9b;
pub const PCG32_DEFAULT_STREAM: u64 = 0xda3e39cb94b95bdb;
pub const PCG32_MULT: u64 = 0x5851f42d4c957f2d;

/// Random number generator
#[derive(Debug,Default,Copy,Clone)]
pub struct Rng {
    state: u64,
    inc: u64,
}

impl Rng {
    pub fn new() -> Self {
        Rng {
            state: PCG32_DEFAULT_STATE,
            inc: PCG32_DEFAULT_STREAM,
        }
    }
    pub fn set_sequence(&mut self, initseq: u64) {
        self.state = 0_u64;
        let (shl, _overflow) = initseq.overflowing_shl(1);
        self.inc = shl | 1;
        self.uniform_uint32();
        let (add, _overflow) = self.state.overflowing_add(PCG32_DEFAULT_STATE);
        self.state = add;
        self.uniform_uint32();
    }
    pub fn uniform_uint32(&mut self) -> u32 {
        let oldstate: u64 = self.state;
        let (mul, _overflow) = oldstate.overflowing_mul(PCG32_MULT);
        let (add, _overflow) = mul.overflowing_add(self.inc);
        self.state = add;
        let (shr, _overflow) = oldstate.overflowing_shr(18);
        let combine = shr ^ oldstate;
        let (shr, _overflow) = combine.overflowing_shr(27);
        let xorshifted: u32 = shr as u32;
        let (shr, _overflow) = oldstate.overflowing_shr(59);
        let rot: u32 = shr as u32;
        // bitwise not in Rust is ! (not the ~ operator like in C)
        let (shr, _overflow) = xorshifted.overflowing_shr(rot);
        let (neg, _overflow) = rot.overflowing_neg();
        let (shl, _overflow) = xorshifted.overflowing_shl(neg & 31);
        shr | shl
    }
    pub fn uniform_uint32_bounded(&mut self, b: u32) -> u32 {
        // bitwise not in Rust is ! (not the ~ operator like in C)
        let threshold = (!b + 1) & b;
        loop {
            let r = self.uniform_uint32();
            if r >= threshold {
                return r % b;
            }
        }
    }
    pub fn uniform_float(&mut self) -> Float {
        (self.uniform_uint32() as Float * 2.3283064365386963e-10 as Float).min(FLOAT_ONE_MINUS_EPSILON)
    }
}
