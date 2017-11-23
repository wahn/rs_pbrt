// std
use std::vec::Vec;
use std::collections::HashMap;
use std::string::String;
use std::sync::Arc;
// pbrt
use core::paramset::ParamSet;
use core::pbrt::Float;
use core::transform::Transform;
use shapes::Shape;
use textures::Texture;

pub fn create_ply_mesh(o2w: Transform,
                       w2o: Transform,
                       reverse_orientation: bool,
                       params: &ParamSet,
                       float_textures: HashMap<String, Arc<Texture<Float> + Send + Sync>>)
                       -> Vec<Arc<Shape + Send + Sync>> {
    // WORK
    Vec::new()
}
