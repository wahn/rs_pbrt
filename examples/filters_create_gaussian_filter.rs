extern crate pbrt;

use pbrt::core::geometry::Vector2f;
use pbrt::core::pbrt::Float;
use pbrt::filters::gaussian::GaussianFilter;

fn main() {
    let xw: Float = 2.0;
    let yw: Float = 2.0;
    let alpha: Float = 2.0;
    let exp_x: Float = (-alpha * xw * xw).exp();
    let exp_y: Float = (-alpha * yw * yw).exp();
    let gaussian_filter = GaussianFilter {
        alpha: alpha,
        exp_x: exp_x,
        exp_y: exp_y,
        radius: Vector2f { x: xw, y: yw },
        inv_radius: Vector2f {
            x: 1.0 / xw,
            y: 1.0 / yw,
        },
    };

    println!("gaussian_filter = {:?}", gaussian_filter);
}
