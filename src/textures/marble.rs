// pbrt
use crate::core::geometry::{Point3f, Vector3f};
use crate::core::interaction::SurfaceInteraction;
use crate::core::pbrt::{Float, Spectrum};
use crate::core::texture::fbm;
use crate::core::texture::{Texture, TextureMapping3D};

// see marble.h

pub struct MarbleTexture {
    pub mapping: Box<dyn TextureMapping3D + Send + Sync>,
    pub octaves: i32,     // default: 8
    pub omega: Float,     // default: 0.5
    pub scale: Float,     // default: 1.0
    pub variation: Float, // default: 0.2
}

impl MarbleTexture {
    pub fn new(
        mapping: Box<dyn TextureMapping3D + Send + Sync>,
        octaves: i32,
        omega: Float,
        scale: Float,
        variation: Float,
    ) -> Self {
        MarbleTexture {
            mapping,
            omega,
            octaves,
            scale,
            variation,
        }
    }
}

impl Texture<Spectrum> for MarbleTexture {
    fn evaluate(&self, si: &SurfaceInteraction) -> Spectrum {
        let mut dpdx: Vector3f = Vector3f::default();
        let mut dpdy: Vector3f = Vector3f::default();
        let mut p: Point3f = self.mapping.map(si, &mut dpdx, &mut dpdy);
        p *= self.scale;
        let marble: Float = p.y
            + fbm(
                &p,
                &(dpdx * self.scale),
                &(dpdy * self.scale),
                self.omega,
                self.octaves,
            );
        let mut t: Float = 0.5 as Float + 0.5 as Float * marble.sin();
        let c: [[Float; 3]; 9] = [
            [0.58 as Float, 0.58 as Float, 0.6 as Float],
            [0.58 as Float, 0.58 as Float, 0.6 as Float],
            [0.58 as Float, 0.58 as Float, 0.6 as Float],
            [0.5 as Float, 0.5 as Float, 0.5 as Float],
            [0.6 as Float, 0.59 as Float, 0.58 as Float],
            [0.58 as Float, 0.58 as Float, 0.6 as Float],
            [0.58 as Float, 0.58 as Float, 0.6 as Float],
            [0.2 as Float, 0.2 as Float, 0.33 as Float],
            [0.58 as Float, 0.58 as Float, 0.6 as Float],
        ];
        let nseg: usize = 6;
        let first: usize = (t * nseg as Float).floor() as usize;
        t = t * nseg as Float - first as Float;
        let c0: Spectrum = Spectrum::from_rgb(&c[first]);
        let c1: Spectrum = Spectrum::from_rgb(&c[first + 1]);
        let c2: Spectrum = Spectrum::from_rgb(&c[first + 2]);
        let c3: Spectrum = Spectrum::from_rgb(&c[first + 3]);
        // Bezier spline evaluated with de Castilejau's algorithm
        let mut s0: Spectrum = c0 * (1.0 as Float - t) + c1 * t;
        let mut s1: Spectrum = c1 * (1.0 as Float - t) + c2 * t;
        let s2: Spectrum = c2 * (1.0 as Float - t) + c3 * t;
        s0 = s0 * (1.0 as Float - t) + s1 * t;
        s1 = s1 * (1.0 as Float - t) + s2 * t;
        // extra scale of 1.5 to increase variation among colors
        (s0 * (1.0 as Float - t) + s1 * t) * 1.5 as Float
    }
}
