extern crate pbrt;

use pbrt::{Bounds2f, BoxFilter, Film, Float, Point2f, Point2i, Vector2f};

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
    let film: Film = Film::new(Point2i { x: 1280, y: 720 },
                               Bounds2f {
                                   p_min: Point2f { x: 0.0, y: 0.0 },
                                   p_max: Point2f { x: 1.0, y: 1.0 },
                               },
                               BoxFilter {
                                   radius: Vector2f { x: 0.5, y: 0.5 },
                                   inv_radius: Vector2f {
                                       x: 1.0 / 0.5,
                                       y: 1.0 / 0.5,
                                   },
                               },
                               35.0,
                               String::from("pbrt.exr"),
                               1.0,
                               std::f64::INFINITY);
    // TODO: std::unique_ptr<Scene> scene(renderOptions->MakeScene());
    // MakeAccelerator()
    // CreateBVHAccelerator()

}
