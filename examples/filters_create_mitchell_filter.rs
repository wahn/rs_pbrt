use pbrt::core::geometry::Vector2f;
use pbrt::core::pbrt::Float;
use pbrt::filters::mitchell::MitchellNetravali;

fn main() {
    let xwidth: Float = 2.0;
    let ywidth: Float = 2.0;
    let b: Float = 1.0 / 3.0;
    let c: Float = 1.0 / 3.0;
    let mitchell_filter = MitchellNetravali::new(xwidth, ywidth, b, c);

    println!("mitchell_filter = {:?}", mitchell_filter);
}
