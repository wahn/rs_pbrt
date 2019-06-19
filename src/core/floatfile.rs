//! Read text files that store a series of floating-point values.

// std
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
// pbrt
use crate::core::pbrt::Float;

pub fn read_float_file(filename: &String, values: &mut Vec<Float>) -> bool {
    let path = Path::new(&filename);
    let result = File::open(path);
    if result.is_ok() {
        let f = result.unwrap();
        let reader = BufReader::new(f);
        for (line_number, line_result) in reader.lines().enumerate() {
            if line_result.is_ok() {
                let line = line_result.unwrap();
                if !line.is_empty() {
                    if line.chars().next() == Some('#') {
                        // ignore comments
                    } else {
                        for token in line.split_whitespace() {
                            match token.parse::<f32>() {
                                Ok(float) => values.push(float),
                                Err(_) => {
                                    println!(
                                    "WARNING: Unexpected text found at line {} of float file {:?}",
                                    line_number, filename
                                );
                                    continue;
                                }
                            }
                        }
                    }
                }
            } else {
                return false;
            }
        }
        true
    } else {
        println!("ERROR: Unable to open file {:?}", filename);
        false
    }
}
