extern crate pbrt;

use pbrt::{Bounds2i, BoxFilter, Film, Point2i, Vector2f};
use std::string::String;

fn main() {
    // see film.cpp CreateFilm()
    let film = Film {
        full_resolution: Point2i { x: 1280, y: 720 },
        diagonal: 35.0,
        filter: BoxFilter {
            radius: Vector2f { x: 0.5, y: 0.5 },
            inv_radius: Vector2f {
                x: 1.0 / 0.5,
                y: 1.0 / 0.5,
            },
        },
        filename: String::from("pbrt.exr"),
        cropped_pixel_bounds: Bounds2i {
            p_min: Point2i { x: 0, y: 0 },
            p_max: Point2i { x: 1280, y: 720 },
        },
    };

    println!("film = {:?}", film);
}
