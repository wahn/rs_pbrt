use pbrt::core::geometry::Vector2f;
use pbrt::core::pbrt::Float;
use pbrt::filters::sinc::LanczosSincFilter;

fn main() {
    let xw: Float = 4.0;
    let yw: Float = 4.0;
    let radius: Vector2f = Vector2f { x: xw, y: yw };
    let tau: Float = 3.0;
    let sinc_filter = LanczosSincFilter::new(&radius, tau);

    println!("sinc_filter = {:?}", sinc_filter);
}
