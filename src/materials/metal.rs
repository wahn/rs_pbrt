//std
use std;
use std::sync::Arc;
// pbrt
use core::interaction::SurfaceInteraction;
use core::material::TransportMode;
use core::pbrt::{Float, Spectrum};
use core::reflection::{Bsdf, Bxdf, LambertianReflection, OrenNayar};
use materials::Material;
use textures::Texture;

pub struct MetalMaterial {
    pub eta: Arc<Texture<Spectrum> + Sync + Send>, // default: 0.5
}
