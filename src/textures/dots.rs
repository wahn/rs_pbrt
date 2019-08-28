// std
use std::sync::Arc;
// others
use num;
// pbrt
use crate::core::interaction::SurfaceInteraction;
use crate::core::texture::{Texture, TextureMapping3D};

// see dots.h

pub struct DotsTexture<T> {
    pub mapping: Box<dyn TextureMapping3D + Send + Sync>,
    pub outside_dot: Arc<dyn Texture<T> + Send + Sync>,
    pub inside_dot: Arc<dyn Texture<T> + Send + Sync>,
}

impl<T: Copy> DotsTexture<T> {
    pub fn new(
        mapping: Box<dyn TextureMapping3D + Send + Sync>,
        outside_dot: Arc<dyn Texture<T> + Send + Sync>,
        inside_dot: Arc<dyn Texture<T> + Send + Sync>,
    ) -> Self {
        DotsTexture {
            mapping,
            outside_dot,
            inside_dot,
        }
    }
}

impl<T: Copy> Texture<T> for DotsTexture<T>
where
    T: num::Zero,
{
    fn evaluate(&self, _si: &SurfaceInteraction) -> T {
        // WORK
        T::from(num::Zero::zero())
    }
}
