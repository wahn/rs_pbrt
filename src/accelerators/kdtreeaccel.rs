// std
use std;
use std::sync::Arc;
// pbrt
use crate::core::geometry::bnd3_union_bnd3;
use crate::core::geometry::{Bounds3f, Ray};
use crate::core::interaction::SurfaceInteraction;
use crate::core::light::AreaLight;
use crate::core::material::Material;
use crate::core::paramset::ParamSet;
use crate::core::pbrt::log_2_int_i32;
use crate::core::pbrt::Float;
use crate::core::primitive::Primitive;

pub struct KdAccelNode {}

pub struct KdTreeAccel {
    pub primitives: Vec<Arc<dyn Primitive + Sync + Send>>,
    pub nodes: Vec<KdAccelNode>,
    pub n_alloced_nodes: i32,
    pub next_free_node: i32,
    pub bounds: Bounds3f,
}

impl KdTreeAccel {
    pub fn new(
        p: Vec<Arc<dyn Primitive + Sync + Send>>,
        isect_cost: i32,
        trav_cost: i32,
        empty_bonus: Float,
        max_prims: i32,
        max_depth: i32,
    ) -> Self {
        let mut max_depth: i32 = max_depth;
        let mut bounds: Bounds3f = Bounds3f::default();
        // build kd-tree for accelerator
        let n_alloced_nodes: i32 = 0;
        let next_free_node: i32 = 0;
        if max_depth <= 0 {
            max_depth = (8.0 as Float + 1.3 as Float * log_2_int_i32(p.len() as i32) as Float)
                .round() as i32;
        }
        // compute bounds for kd-tree construction
        let mut prim_bounds: Vec<Bounds3f> = Vec::with_capacity(p.len());
        for i in 0..p.len() {
            let b: Bounds3f = p[i].world_bound();
            bounds = bnd3_union_bnd3(&bounds, &b);
            prim_bounds.push(b);
        }
        // allocate working memory for kd-tree construction
        // initialize _primNums_ for kd-tree construction
        // start recursive construction of kd-tree
        KdTreeAccel {
            primitives: p,
            nodes: Vec::new(),
            n_alloced_nodes,
            next_free_node,
            bounds,
        }
    }
    pub fn create(prims: Vec<Arc<dyn Primitive + Send + Sync>>, ps: &ParamSet) -> Arc<KdTreeAccel> {
        let isect_cost: i32 = ps.find_one_int("intersectcost", 80);
        let trav_cost: i32 = ps.find_one_int("traversalcost", 1);
        let empty_bonus: Float = ps.find_one_float("emptybonus", 0.5 as Float);
        let max_prims: i32 = ps.find_one_int("maxprims", 1);
        let max_depth: i32 = ps.find_one_int("maxdepth", -1);
        Arc::new(KdTreeAccel::new(
            prims.clone(),
            isect_cost,
            trav_cost,
            empty_bonus,
            max_prims,
            max_depth,
        ))
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
