extern crate pbrt;

use pbrt::core::geometry::Vector2f;
use pbrt::core::pbrt::Float;
use pbrt::filters::boxfilter::BoxFilter;

fn main() {
    let xw: Float = 0.5;
    let yw: Float = 0.5;
    let box_filter = BoxFilter {
        radius: Vector2f { x: xw, y: yw },
        inv_radius: Vector2f {
            x: 1.0 / xw,
            y: 1.0 / yw,
        },
    };

    println!("box_filter = {:?}", box_filter);
}
