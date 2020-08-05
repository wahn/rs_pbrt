//! Random Number Generator

use hexf::*;

// pbrt
use crate::core::pbrt::Float;

// see rng.h

//#ifndef PBRT_HAVE_HEX_FP_CONSTANTS
//pub const FLOAT_ONE_MINUS_EPSILON: Float = 0.99999994;
//#else
pub const FLOAT_ONE_MINUS_EPSILON: Float = hexf32!("0x1.fffffep-1");
//#endif
pub const PCG32_DEFAULT_STATE: u64 = 0x853c_49e6_748f_ea9b;
pub const PCG32_DEFAULT_STREAM: u64 = 0xda3e_39cb_94b9_5bdb;
pub const PCG32_MULT: u64 = 0x5851_f42d_4c95_7f2d;

/// Random number generator
#[derive(Debug, Default, Copy, Clone)]
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
        self.inc = initseq.wrapping_shl(1) | 1;
        self.uniform_uint32();
        self.state = self.state.wrapping_add(PCG32_DEFAULT_STATE);
        self.uniform_uint32();
    }
    pub fn uniform_uint32(&mut self) -> u32 {
        let oldstate: u64 = self.state;
        // C++: state = oldstate * PCG32_MULT + inc;
        self.state = oldstate.wrapping_mul(PCG32_MULT).wrapping_add(self.inc);
        // C++: uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
        let xorshifted: u32 = (oldstate.wrapping_shr(18) ^ oldstate).wrapping_shr(27) as u32;
        // C++: uint32_t rot = (uint32_t)(oldstate >> 59u);
        let rot: u32 = oldstate.wrapping_shr(59) as u32;
        // C++: return (xorshifted >> rot) | (xorshifted << ((~rot + 1u) & 31));
        xorshifted.wrapping_shr(rot)
            | xorshifted.wrapping_shl(rot.wrapping_neg().wrapping_add(1_u32) & 31)
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
        //#ifndef PBRT_HAVE_HEX_FP_CONSTANTS
        // (self.uniform_uint32() as Float * 2.3283064365386963e-10 as Float)
        //     .min(FLOAT_ONE_MINUS_EPSILON)
        //#else
        (self.uniform_uint32() as Float * hexf32!("0x1.0p-32") as Float)
            .min(FLOAT_ONE_MINUS_EPSILON)
        //#endif
    }
}
