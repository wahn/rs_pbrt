// std
use std;
use std::sync::Arc;
// pbrt
use crate::core::geometry::{Bounds3f, Ray};
use crate::core::interaction::SurfaceInteraction;
use crate::core::light::AreaLight;
use crate::core::material::Material;
use crate::core::paramset::ParamSet;
use crate::core::pbrt::Float;
use crate::core::primitive::Primitive;

pub struct KdAccelNode {}

pub struct KdTreeAccel {
    pub primitives: Vec<Arc<dyn Primitive + Sync + Send>>,
    pub nodes: Vec<KdAccelNode>,
}

impl KdTreeAccel {
    pub fn new(p: Vec<Arc<dyn Primitive + Sync + Send>>) -> Self {
        KdTreeAccel {
            primitives: p,
            nodes: Vec::new(),
        }
    }
    pub fn create(prims: Vec<Arc<dyn Primitive + Send + Sync>>, ps: &ParamSet) -> Arc<KdTreeAccel> {
        let isect_cost: i32 = ps.find_one_int("intersectcost", 80);
        let trav_cost: i32 = ps.find_one_int("traversalcost", 1);
        let empty_bonus: Float = ps.find_one_float("emptybonus", 0.5 as Float);
        let max_prims: i32 = ps.find_one_int("maxprims", 1);
        let max_depth: i32 = ps.find_one_int("maxdepth", -1);
        Arc::new(KdTreeAccel::new(prims.clone()))
    }
}

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
