extern crate pbrt;

use pbrt::core::geometry::Vector2f;
use pbrt::core::pbrt::Float;
use pbrt::filters::triangle::TriangleFilter;

fn main() {
    let xw: Float = 2.0;
    let yw: Float = 2.0;
    let triangle_filter = TriangleFilter {
        radius: Vector2f { x: xw, y: yw },
        inv_radius: Vector2f {
            x: 1.0 / xw,
            y: 1.0 / yw,
        },
    };

    println!("triangle_filter = {:?}", triangle_filter);
}
