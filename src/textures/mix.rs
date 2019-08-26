// std
use std::ops::{Add, Mul};
use std::sync::Arc;
// pbrt
use crate::core::interaction::SurfaceInteraction;
use crate::core::pbrt::Float;
use crate::core::texture::Texture;

pub struct MixTexture<T> {
    pub tex1: Arc<dyn Texture<T> + Send + Sync>,
    pub tex2: Arc<dyn Texture<T> + Send + Sync>,
    pub amount: Arc<dyn Texture<Float> + Send + Sync>,
}

impl<T: Copy> MixTexture<T> {
    pub fn new(
        tex1: Arc<dyn Texture<T> + Send + Sync>,
        tex2: Arc<dyn Texture<T> + Send + Sync>,
        amount: Arc<dyn Texture<Float> + Send + Sync>,
    ) -> Self {
        MixTexture { tex1, tex2, amount }
    }
}

impl<T: Copy> Texture<T> for MixTexture<T>
where
    T: Add<Output = T>,
    T: Mul<Output = T>,
    T: From<Float>,
{
    fn evaluate(&self, si: &SurfaceInteraction) -> T {
        let t1: T = self.tex1.evaluate(si);
        let t2: T = self.tex2.evaluate(si);
        let amt: Float = self.amount.evaluate(si);
        t1 * T::from(1.0 as Float - amt) + t2 * T::from(amt)
    }
}
