// pbrt
use core::pbrt::find_interval;
use core::pbrt::Float;

/// Calculates an offset and four weights for Catmull-Rom spline
/// interpolation.
pub fn catmull_rom_weights(
    nodes: &Vec<Float>,
    x: Float,
    offset: &mut i32,
    weights: &mut [Float; 4],
) -> bool {
    // return _false_ if _x_ is out of bounds
    if !(x >= *nodes.first().unwrap() && x <= *nodes.last().unwrap()) {
        return false;
    }
    // search for the interval _idx_ containing _x_
    let idx: usize = find_interval(nodes.len(), |index| nodes[index as usize] <= x);
    *offset = idx as i32 - 1;
    let x0: Float = nodes[idx];
    let x1: Float = nodes[idx + 1];
    // compute the $t$ parameter and powers
    let t: Float = (x - x0) / (x1 - x0);
    let t2: Float = t * t;
    let t3: Float = t2 * t;
    // compute initial node weights $w_1$ and $w_2$
    weights[1] = 2.0 as Float * t3 - 3.0 as Float * t2 + 1.0 as Float;
    weights[2] = -2.0 as Float * t3 + 3.0 as Float * t2;
    // compute first node weight $w_0$
    if idx > 0_usize {
        let w0: Float = (t3 - 2.0 as Float * t2 + t) * (x1 - x0) / (x1 - nodes[idx - 1]);
        weights[0] = -w0;
        weights[2] += w0;
    } else {
        let w0: Float = t3 - 2.0 as Float * t2 + t;
        weights[0] = 0.0 as Float;
        weights[1] -= w0;
        weights[2] += w0;
    }
    // compute last node weight $w_3$
    if (idx + 2) < nodes.len() {
        let w3: Float = (t3 - t2) * (x1 - x0) / (nodes[idx + 2] - x0);
        weights[1] -= w3;
        weights[3] = w3;
    } else {
        let w3: Float = t3 - t2;
        weights[1] -= w3;
        weights[2] += w3;
        weights[3] = 0.0 as Float;
    }
    true
}
