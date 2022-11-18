//! Provides a queue of block indices that the sampler positions can
//! be initialized from for the worker threads. The queue itself is
//! not changed after creation we simply work through it with an
//! atomic counter to track the index of the next block to work on.

use std::sync::atomic::{AtomicUsize, Ordering};

// see github/tray_rust/src/sampler/block_queue.rs

/// The queue of blocks to be worked on shared immutably between worker threads.
pub struct BlockQueue {
    /// The block indices of blocks to work on for the image
    blocks: Vec<(u32, u32)>,
    /// Get the dimensions of an individual block
    dimensions: (u32, u32),
    /// Index of the next block to be worked on
    next: AtomicUsize,
}

impl BlockQueue {
    /// Create a block queue for the image with dimensions `img`.
    /// Panics if the image is not evenly broken into blocks of dimension `dim`
    pub fn new(img: (u32, u32), dim: (u32, u32), select_blocks: (usize, usize)) -> BlockQueue {
        if img.0 % dim.0 != 0 || img.1 % dim.1 != 0 {
            panic!(
                "Image with dimension {:?} not evenly divided by blocks of {:?}",
                img, dim
            );
        }
        let num_blocks = (img.0 / dim.0, img.1 / dim.1);
        // TODO: the .. operator precedence is very low so we need this paren here at the moment
        // once (hopefully) it's raised we can remove the parens
        let mut blocks: Vec<(u32, u32)> = (0..num_blocks.0 * num_blocks.1)
            .map(|i| (i % num_blocks.0, i / num_blocks.0))
            .collect();
        blocks.sort_by(|a, b| morton2(*a).cmp(&morton2(*b)));
        // If we're only rendering a subset of the blocks then filter our list down
        if select_blocks.1 > 0 {
            blocks = blocks
                .into_iter()
                .skip(select_blocks.0)
                .take(select_blocks.1)
                .collect();
        }
        if blocks.is_empty() {
            println!("Warning: This block queue is empty!");
        }
        BlockQueue {
            blocks,
            dimensions: dim,
            next: AtomicUsize::new(0),
        }
    }
    /// Get the dimensions of an individual block in the queue
    pub fn block_dim(&self) -> (u32, u32) {
        self.dimensions
    }
    /// Get an iterator to work through the queue
    pub fn iter(&self) -> BlockQueueIterator {
        BlockQueueIterator { queue: self }
    }
    /// Get the next block in the queue or None if the queue is finished
    pub fn next(&self) -> Option<(u32, u32)> {
        let i = self.next.fetch_add(1, Ordering::AcqRel);
        if i >= self.blocks.len() {
            None
        } else {
            Some(self.blocks[i])
        }
    }
    /// Get the length of the queue
    pub fn len(&self) -> usize {
        self.blocks.len()
    }
    /// Check if the queue is empty
    pub fn is_empty(&self) -> bool {
        self.next.load(Ordering::Acquire) >= self.blocks.len()
    }
}

/// Iterator to work through the queue safely
pub struct BlockQueueIterator<'a> {
    queue: &'a BlockQueue,
}

impl<'a> Iterator for BlockQueueIterator<'a> {
    type Item = (u32, u32);
    fn next(&mut self) -> Option<(u32, u32)> {
        self.queue.next()
    }
}

// see github/tray_rust/src/sampler/morton.rs

///! Provides utilities for 2D Morton code generation using Fabian
///! Giesen's Morton code decoding functions, see [his post on Morton
///! codes](https://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/)

/// Insert a 0 bit between each of the low 16 bits of x
fn part1_by1(mut x: u32) -> u32 {
    // x = ---- ---- ---- ---- fedc ba98 7654 3210
    x &= 0x0000_ffff;
    // x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x ^ (x << 8)) & 0x00ff_00ff;
    // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x ^ (x << 4)) & 0x0f0f_0f0f;
    // x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x ^ (x << 2)) & 0x3333_3333;
    // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    (x ^ (x << 1)) & 0x5555_5555
}
/// Compute the Morton code for the `(x, y)` position.
fn morton2(p: (u32, u32)) -> u32 {
    (part1_by1(p.1) << 1) + part1_by1(p.0)
}
