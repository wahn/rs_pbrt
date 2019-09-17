// std
use std;
use std::sync::Arc;
// pbrt
use crate::core::geometry::bnd3_union_bnd3;
use crate::core::geometry::{Bounds3f, Ray, Vector3f};
use crate::core::interaction::SurfaceInteraction;
use crate::core::light::AreaLight;
use crate::core::material::Material;
use crate::core::paramset::ParamSet;
use crate::core::pbrt::log_2_int_i32;
use crate::core::pbrt::Float;
use crate::core::primitive::Primitive;

pub struct KdAccelNode {}

#[derive(Debug, Clone)]
pub enum EdgeType {
    Start,
    End,
}

#[derive(Debug, Clone)]
pub struct BoundEdge {
    pub t: Float,
    pub prim_num: usize,
    pub edge_type: EdgeType,
}

impl BoundEdge {
    pub fn new(t: Float, prim_num: usize, starting: bool) -> Self {
        let edge_type: EdgeType;
        if starting {
            edge_type = EdgeType::Start;
        } else {
            edge_type = EdgeType::End;
        }
        BoundEdge {
            t,
            prim_num,
            edge_type,
        }
    }
}

impl Default for BoundEdge {
    fn default() -> Self {
        BoundEdge {
            t: 0.0 as Float,
            prim_num: 0_usize,
            edge_type: EdgeType::Start,
        }
    }
}

pub struct KdTreeAccel {
    pub isect_cost: i32,
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
        let p_len: usize = p.len();
        let mut max_depth: i32 = max_depth;
        let mut bounds: Bounds3f = Bounds3f::default();
        // build kd-tree for accelerator
        let n_alloced_nodes: i32 = 0;
        let next_free_node: i32 = 0;
        if max_depth <= 0 {
            max_depth =
                (8.0 as Float + 1.3 as Float * log_2_int_i32(p_len as i32) as Float).round() as i32;
        }
        // compute bounds for kd-tree construction
        let mut prim_bounds: Vec<Bounds3f> = Vec::with_capacity(p_len);
        for i in 0..p_len {
            let b: Bounds3f = p[i].world_bound();
            bounds = bnd3_union_bnd3(&bounds, &b);
            prim_bounds.push(b);
        }
        // allocate working memory for kd-tree construction
        let mut edges: [Vec<BoundEdge>; 3] = [
            Vec::with_capacity(2 * p_len),
            Vec::with_capacity(2 * p_len),
            Vec::with_capacity(2 * p_len),
        ];
        // std::unique_ptr<int[]> prims0(new int[primitives.size()]);
        // std::unique_ptr<int[]> prims1(new int[(maxDepth + 1) * primitives.size()]);
        // initialize _prim_nums_ for kd-tree construction
        let mut prim_nums: Vec<usize> = Vec::with_capacity(p_len);
        for i in 0..p_len {
            prim_nums.push(i);
        }
        // start recursive construction of kd-tree
        let mut kd_tree: KdTreeAccel = KdTreeAccel {
            isect_cost,
            primitives: p,
            nodes: Vec::new(),
            n_alloced_nodes,
            next_free_node,
            bounds,
        };
        // build_tree(0, bounds, prim_bounds, prim_nums.get(), primitives.size(),
        //           maxDepth, edges, prims0.get(), prims1.get());
        KdTreeAccel::build_tree(
            &mut kd_tree,
            0 as i32,
            &bounds,
            &prim_bounds,
            &prim_nums[..],
            p_len,
            &mut edges,
        );
        kd_tree
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
    pub fn build_tree(
        &mut self,
        node_num: i32,
        node_bounds: &Bounds3f,
        all_prim_bounds: &Vec<Bounds3f>,
        prim_nums: &[usize],
        n_primitives: usize,
        edges: &mut [Vec<BoundEdge>; 3],
    ) {
        assert_eq!(node_num, self.next_free_node);
        if self.next_free_node == self.n_alloced_nodes {}
        self.next_free_node += 1;
        // ...
        // choose split axis position for interior node
        let mut best_axis: i32 = -1;
        let mut best_offset: i32 = -1;
        let mut best_cost: Float = std::f32::INFINITY;
        let old_cost: Float = self.isect_cost as Float * n_primitives as Float;
        let total_sa: Float = node_bounds.surface_area();
        let inv_total_sa: Float = 1.0 as Float / total_sa;
        let d: Vector3f = node_bounds.p_max - node_bounds.p_min;
        // choose which axis to split along
        let mut axis: u8 = node_bounds.maximum_extent();
        let mut retries: u8 = 0;
        // avoid 'goto retrySplit;'
        {
            // initialize edges for _axis_
            for i in 0..n_primitives {
                let pn: usize = prim_nums[i];
                let bounds: &Bounds3f = &all_prim_bounds[pn];
                edges[axis as usize][2 * i] = BoundEdge::new(bounds.p_min[axis], pn, true);
                edges[axis as usize][2 * i + 1] = BoundEdge::new(bounds.p_max[axis], pn, false);
            }
        }
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
