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
    pub fn new(v: Float) -> AtomicFloat {
        AtomicFloat {
            bits: Arc::new(Atomic::new(float_to_bits(v))),
        }
    }
    pub fn add(&mut self, v: Float) {
        let arc: &mut Arc<Atomic<u32>> = &mut self.bits;
        let option: Option<&mut Atomic<u32>> = Arc::get_mut(arc);
        println!("Ready to unwrap? {:?}", option.is_some());
        let atom: &mut Atomic<u32> = option.unwrap();
        let mut old_bits: u32 = atom.load(Ordering::Relaxed);
        loop {
            let f: Float = bits_to_float(old_bits);
            print!("{:?} + {:?}: ", f, v);
            let new_bits: u32 = float_to_bits(f + v);
            match atom.compare_exchange_weak(
                old_bits,
                new_bits,
                Ordering::SeqCst,
                Ordering::Relaxed,
            ) {
                Ok(_) => {
                    println!("Ok");
                    break;
                },
                Err(x) => {
                    println!("Err({:?})", x);
                    old_bits = x;
                },
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
