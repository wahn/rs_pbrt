// std
use std::ops::Mul;
use std::sync::Arc;
// pbrt
use crate::core::interaction::SurfaceInteraction;
use crate::core::texture::Texture;

pub struct ScaleTexture<T> {
    pub tex1: Arc<Texture<T> + Send + Sync>,
    pub tex2: Arc<Texture<T> + Send + Sync>,
}

impl<T: Copy> ScaleTexture<T> {
    pub fn new(tex1: Arc<Texture<T> + Send + Sync>, tex2: Arc<Texture<T> + Send + Sync>) -> Self {
        ScaleTexture {
            tex1: tex1,
            tex2: tex2,
        }
    }
}

impl<T: Copy> Texture<T> for ScaleTexture<T>
where
    T: Mul<Output = T>,
{
    fn evaluate(&self, si: &SurfaceInteraction) -> T {
        self.tex1.evaluate(si) * self.tex2.evaluate(si)
    }
}
