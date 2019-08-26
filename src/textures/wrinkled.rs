// pbrt
use crate::core::interaction::SurfaceInteraction;
use crate::core::pbrt::Float;
use crate::core::texture::Texture;

// see wrinkled.h

pub struct WrinkledTexture {}

impl WrinkledTexture {}

impl<T> Texture<T> for WrinkledTexture
where
    T: From<Float>,
{
    fn evaluate(&self, _si: &SurfaceInteraction) -> T {
        // WORK
        T::from(0.0 as Float)
    }
}
