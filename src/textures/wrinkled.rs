// pbrt
use crate::core::geometry::{Point3f, Vector3f};
use crate::core::interaction::SurfaceInteraction;
use crate::core::pbrt::Float;
use crate::core::texture::turbulence;
use crate::core::texture::{Texture, TextureMapping3D};

// see wrinkled.h

pub struct WrinkledTexture {
    pub mapping: Box<TextureMapping3D>,
    pub octaves: i32, // default: 8
    pub omega: Float, // default: 0.5
}

impl WrinkledTexture {
    pub fn new(
        mapping: Box<TextureMapping3D>,
        octaves: i32,
        omega: Float,
    ) -> Self {
        WrinkledTexture {
            mapping,
            omega,
            octaves,
        }
    }
}

impl<T> Texture<T> for WrinkledTexture
where
    T: From<Float>,
{
    fn evaluate(&self, si: &SurfaceInteraction) -> T {
        let mut dpdx: Vector3f = Vector3f::default();
        let mut dpdy: Vector3f = Vector3f::default();
        let p: Point3f = self.mapping.map(si, &mut dpdx, &mut dpdy);
        T::from(turbulence(&p, &dpdx, &dpdy, self.omega, self.octaves))
    }
}
