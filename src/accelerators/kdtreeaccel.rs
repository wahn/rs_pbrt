// std
use std;
use std::sync::Arc;
// pbrt
use crate::core::geometry::{Bounds3f, Ray};
use crate::core::interaction::SurfaceInteraction;
use crate::core::light::AreaLight;
use crate::core::material::Material;
use crate::core::primitive::Primitive;

pub struct KdAccelNode {}

pub struct KdTreeAccel {
    pub primitives: Vec<Arc<dyn Primitive + Sync + Send>>,
    pub nodes: Vec<KdAccelNode>,
}

impl KdTreeAccel {}

impl Primitive for KdTreeAccel {
    fn world_bound(&self) -> Bounds3f {
        // WORK
        Bounds3f::default()
    }
    fn intersect(&self, ray: &mut Ray) -> Option<SurfaceInteraction> {
        // WORK
        None
    }
    fn intersect_p(&self, ray: &Ray) -> bool {
        // WORK
        false
    }
    fn get_material(&self) -> Option<Arc<dyn Material + Send + Sync>> {
        None
    }
    fn get_area_light(&self) -> Option<Arc<dyn AreaLight + Send + Sync>> {
        None
    }
}
