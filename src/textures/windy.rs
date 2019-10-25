// pbrt
use crate::core::geometry::{Point3f, Vector3f};
use crate::core::interaction::SurfaceInteraction;
use crate::core::pbrt::Float;
use crate::core::texture::fbm;
use crate::core::texture::{Texture, TextureMapping3D};

// see windy.h

pub struct WindyTexture {
    pub mapping: Box<TextureMapping3D>,
}

impl WindyTexture {
    pub fn new(mapping: Box<TextureMapping3D>) -> Self {
        WindyTexture { mapping }
    }
}

impl<T> Texture<T> for WindyTexture
where
    T: From<Float>,
{
    fn evaluate(&self, si: &SurfaceInteraction) -> T {
        let mut dpdx: Vector3f = Vector3f::default();
        let mut dpdy: Vector3f = Vector3f::default();
        let p: Point3f = self.mapping.map(si, &mut dpdx, &mut dpdy);
        let wind_strength: Float = fbm(
            &(p * 0.1 as Float),
            &(dpdx * 0.1 as Float),
            &(dpdy * 0.1 as Float),
            0.5 as Float,
            3_i32,
        );
        let wave_height: Float = fbm(&p, &dpdx, &dpdy, 0.5 as Float, 6_i32);
        T::from(wind_strength.abs() * wave_height)
    }
}
