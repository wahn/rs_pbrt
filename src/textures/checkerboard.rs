// std
use std::sync::Arc;
// pbrt
use core::interaction::SurfaceInteraction;
use textures::{Texture, TextureMapping2D};
use geometry::{Point2f, Vector2f};

// checkerboard.h

pub struct Checkerboard2DTexture<T> {
    pub tex1: Arc<Texture<T> + Send + Sync>,
    pub tex2: Arc<Texture<T> + Send + Sync>,
    pub mapping: Box<TextureMapping2D + Send + Sync>,
    // TODO: const AAMethod aaMethod;
}

impl<T: Copy> Checkerboard2DTexture<T> {
    pub fn new(mapping: Box<TextureMapping2D + Send + Sync>,
               tex1: Arc<Texture<T> + Send + Sync>,
               tex2: Arc<Texture<T> + Send + Sync>// , TODO: aaMethod
    ) -> Self {
        Checkerboard2DTexture {
            tex1: tex1,
            tex2: tex2,
            mapping: mapping,
        }
    }
}

impl<T: Copy> Texture<T> for Checkerboard2DTexture<T> {
    fn evaluate(&self, si: &SurfaceInteraction) -> T {
        let mut dstdx: Vector2f = Vector2f::default();
        let mut dstdy: Vector2f = Vector2f::default();
        let st: Point2f = self.mapping.map(si, &mut dstdx, &mut dstdy);
        // TODO: if (aaMethod == AAMethod::None) {
        if (st.x.floor() as u32 + st.y.floor() as u32) % 2 == 0 {
            self.tex1.evaluate(si)
        } else {
            self.tex2.evaluate(si)
        }
    }
}
