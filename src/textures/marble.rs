// pbrt
use crate::core::interaction::SurfaceInteraction;
use crate::core::pbrt::Spectrum;
use crate::core::texture::Texture;

// see marble.h

pub struct MarbleTexture {
}

impl MarbleTexture {
}

impl Texture<Spectrum> for MarbleTexture {
    fn evaluate(&self, _si: &SurfaceInteraction) -> Spectrum {
        // WORK
        Spectrum::default()
    }
}
