extern crate pbrt;
extern crate num_cpus;

// use pbrt::{Bounds2, Point2};

fn main() {
    // TODO: process command-line arguments
    // print welcome banner
    let date = "Jan  1 2017";
    let time = "13:58:59";
    let num_cores = num_cpus::get();
    println!("pbrt version 3 (build {} at {}) [Detected {} cores]",
             date,
             time,
             num_cores);
    // TODO: detect debug build
    println!("Copyright (c)1998-2016 Matt Pharr, Greg Humphreys, and Wenzel Jakob.");
    println!("The source code to pbrt (but *not* the book contents) is covered by the BSD \
              License.");
    println!("See the file LICENSE.txt for the conditions of the license.");
    // TODO: process scene description
}
