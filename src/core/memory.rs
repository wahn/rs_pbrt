// std
use std;
// others
use num;
use std::ops::{Add, Index, IndexMut};

// see memory.h

const LOG_BLOCK_SIZE: usize = 3;
const BLOCK_SIZE: usize = 1 << LOG_BLOCK_SIZE;

fn round_up(x: usize) -> usize {
    (x + BLOCK_SIZE - 1) & !(BLOCK_SIZE - 1)
}

#[derive(Debug, Clone, Default)]
pub struct BlockedArray<T> {
    pub data: Vec<T>,
    pub u_res: usize,
    pub v_res: usize,
    pub u_blocks: usize,
    log_block_size: usize,
    block_size: usize,
}

impl<T> BlockedArray<T>
where
    T: num::Zero + Clone + Add<T, Output = T>,
{
    pub fn new(u_res: usize, v_res: usize) -> BlockedArray<T> {
        let data = vec![num::Zero::zero(); round_up(u_res) * round_up(v_res)];
        BlockedArray {
            u_res: u_res,
            v_res: v_res,
            u_blocks: round_up(u_res) >> LOG_BLOCK_SIZE,
            log_block_size: LOG_BLOCK_SIZE,
            block_size: BLOCK_SIZE,
            data: data,
        }
    }
    pub fn new_from(u_res: usize, v_res: usize, d: &[T]) -> BlockedArray<T> {
        let mut ba = Self::new(u_res, v_res);
        for u in 0..u_res {
            for v in 0..v_res {
                ba[(u, v)] = d[v * u_res + u].clone();
            }
        }
        ba
    }
    pub fn u_size(&self) -> usize {
        self.u_res
    }
    pub fn v_size(&self) -> usize {
        self.v_res
    }
    pub fn block_size(&self) -> usize {
        1 << self.log_block_size
    }
    pub fn block(&self, a: usize) -> usize {
        a >> self.log_block_size
    }
    pub fn offset(&self, a: usize) -> usize {
        a & (self.block_size() - 1)
    }
}

impl<T> Index<(usize, usize)> for BlockedArray<T>
where
    T: num::Zero + std::clone::Clone + Add<T, Output = T>,
{
    type Output = T;
    fn index(&self, i: (usize, usize)) -> &T {
        let (u, v) = i;
        let bu = self.block(u);
        let bv = self.block(v);
        let ou = self.offset(u);
        let ov = self.offset(v);
        let offset = self.block_size() * self.block_size() * (self.u_blocks * bv + bu)
            + self.block_size() * ov + ou;
        &self.data[offset]
    }
}

impl<T> IndexMut<(usize, usize)> for BlockedArray<T>
where
    T: num::Zero + std::clone::Clone + Add<T, Output = T>,
{
    fn index_mut(&mut self, i: (usize, usize)) -> &mut T {
        let (u, v) = i;
        let bu = self.block(u);
        let bv = self.block(v);
        let ou = self.offset(u);
        let ov = self.offset(v);
        let offset = self.block_size() * self.block_size() * (self.u_blocks * bv + bu)
            + self.block_size() * ov + ou;
        &mut self.data[offset]
    }
}
