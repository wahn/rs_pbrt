// pbrt
use crate::core::geometry::{Point3f, Vector3f};
use crate::core::interaction::SurfaceInteraction;
use crate::core::pbrt::Float;
use crate::core::texture::fbm;
use crate::core::texture::{Texture, TextureMapping3D};

// see fbm.h

pub struct FBmTexture {
    pub mapping: Box<TextureMapping3D + Send + Sync>,
    pub omega: Float, // default: 0.5
    pub octaves: i32, // default: 8
}

impl FBmTexture {
    pub fn new(mapping: Box<TextureMapping3D + Send + Sync>, octaves: i32, omega: Float) -> Self {
        FBmTexture {
            mapping,
            omega,
            octaves,
        }
    }
}

impl<T> Texture<T> for FBmTexture
where
    T: From<Float>,
{
    fn evaluate(&self, si: &SurfaceInteraction) -> T {
        let mut dpdx: Vector3f = Vector3f::default();
        let mut dpdy: Vector3f = Vector3f::default();
        let p: Point3f = self.mapping.map(si, &mut dpdx, &mut dpdy);
        T::from(fbm(&p, &dpdx, &dpdy, self.omega, self.octaves))
    }
}
