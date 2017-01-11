extern crate pbrt;

use pbrt::{Bounds2i, BoxFilter, Film, Float, Point2i, Vector2f};

fn main() {
    // see api.cpp pbrtWorldEnd()

    // TODO: std::unique_ptr<Integrator> integrator(renderOptions->MakeIntegrator());
    // MakeCamera()
    // MakeFilter(FilterName, FilterParams)
    // see box.cpp CreateBoxFilter()
    let xw: Float = 0.5;
    let yw: Float = 0.5;
    let filter = BoxFilter {
        radius: Vector2f { x: xw, y: yw },
        inv_radius: Vector2f {
            x: 1.0 / xw,
            y: 1.0 / yw,
        },
    };
    // MakeFilm()
    // see film.cpp CreateFilm()
    let film = Film {
        full_resolution: Point2i { x: 1280, y: 720 },
        diagonal: 35.0,
        filter: filter,
        filename: String::from("pbrt.exr"),
        cropped_pixel_bounds: Bounds2i {
            p_min: Point2i { x: 0, y: 0 },
            p_max: Point2i { x: 1280, y: 720 },
        },
    };
    // TODO: std::unique_ptr<Scene> scene(renderOptions->MakeScene());
    // ...
}
