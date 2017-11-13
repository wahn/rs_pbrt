// std
use std::f32::consts::PI;
// pbrt
use core::interaction::SurfaceInteraction;
use core::pbrt::Float;
use geometry::{Point2f, Vector2f, Vector3f};
use geometry::vec3_dot_vec3;

pub mod checkerboard;
pub mod constant;
pub mod imagemap;

// see texture.h

pub trait TextureMapping2D {
    fn map(&self, si: &SurfaceInteraction, dstdx: &mut Vector2f, dstdy: &mut Vector2f) -> Point2f;
}

#[derive(Debug,Default,Copy,Clone)]
pub struct UVMapping2D {
    pub su: Float,
    pub sv: Float,
    pub du: Float,
    pub dv: Float,
}

impl TextureMapping2D for UVMapping2D {
    fn map(&self, si: &SurfaceInteraction, dstdx: &mut Vector2f, dstdy: &mut Vector2f) -> Point2f {
        // compute texture differentials for 2D identity mapping
        *dstdx = Vector2f {
            x: si.dudx * self.su,
            y: si.dvdx * self.sv,
        };
        *dstdy = Vector2f {
            x: si.dudy * self.su,
            y: si.dvdy * self.sv,
        };
        Point2f {
            x: si.uv[0] * self.su + self.du,
            y: si.uv[1] * self.sv + self.dv,
        }
    }
}

#[derive(Debug,Default,Copy,Clone)]
pub struct PlanarMapping2D {
    pub vs: Vector3f,
    pub vt: Vector3f,
    pub ds: Float,
    pub dt: Float,
}

impl TextureMapping2D for PlanarMapping2D {
    fn map(&self, si: &SurfaceInteraction, dstdx: &mut Vector2f, dstdy: &mut Vector2f) -> Point2f {
        let vec: Vector3f = Vector3f {
            x: si.p.x,
            y: si.p.y,
            z: si.p.z,
        };
        *dstdx = Vector2f {
            x: vec3_dot_vec3(si.dpdx, self.vs),
            y: vec3_dot_vec3(si.dpdx, self.vt),
        };
        *dstdy = Vector2f {
            x: vec3_dot_vec3(si.dpdy, self.vs),
            y: vec3_dot_vec3(si.dpdy, self.vt),
        };
        Point2f {
            x: self.ds + vec3_dot_vec3(vec, self.vs),
            y: self.dt + vec3_dot_vec3(vec, self.vt),
        }
    }
}

pub trait Texture<T> {
    fn evaluate(&self, si: &SurfaceInteraction) -> T;
}

pub fn lanczos(x: Float, tau: Float) -> Float {
    let mut x: Float = x;
    x = x.abs();
    if x < 1e-5 as Float {
        return 1.0 as Float;
    }
    if x > 1.0 as Float {
        return 0.0 as Float;
    }
    x *= PI;
    let s: Float = (x * tau).sin() / (x * tau);
    let lanczos: Float = x.sin() / x;
    s * lanczos
}
