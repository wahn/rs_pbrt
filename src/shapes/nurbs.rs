// pbrt
use core::pbrt::Float;

// see nurbs.cpp

#[derive(Debug, Default, Copy, Clone)]
pub struct Homogeneous3 {
    pub x: Float,
    pub y: Float,
    pub z: Float,
    pub w: Float,
}
