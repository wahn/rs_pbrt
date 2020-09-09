// std
use std::sync::Arc;
// pbrt
use crate::core::geometry::bnd3_union_bnd3f;
use crate::core::geometry::{Bounds3f, Ray, Vector3f, XYZEnum};
use crate::core::interaction::SurfaceInteraction;
use crate::core::light::Light;
use crate::core::material::Material;
use crate::core::paramset::ParamSet;
use crate::core::pbrt::log_2_int_i32;
use crate::core::pbrt::Float;
use crate::core::primitive::Primitive;

pub const MAX_TODO: usize = 64;

#[repr(C)]
union PrivateUnion {
    split: Float,
    one_primitive: i32,
    primitive_indices_offset: i32,
}

#[repr(C)]
pub union PublicUnion {
    pub flags: i32,
    pub n_prims: i32,
    pub above_child: i32,
}

#[derive(Copy, Clone, Default)]
pub struct KdToDo<'n> {
    pub node: Option<&'n KdAccelNode>,
    pub idx: usize,
    pub t_min: Float,
    pub t_max: Float,
}

pub struct KdAccelNode {
    priv_union: PrivateUnion,
    pub pub_union: PublicUnion,
}

impl KdAccelNode {
    pub fn init_leaf(&mut self, prim_nums: &[usize], np: usize, primitive_indices: &mut Vec<i32>) {
        self.pub_union.flags = 3;
        let n_prims: i32;
        unsafe {
            n_prims = self.pub_union.n_prims;
        }
        self.pub_union.n_prims = n_prims | ((np as i32) << 2);
        // store primitive ids for leaf node
        match np {
            0 => self.priv_union.one_primitive = 0_i32,
            1 => self.priv_union.one_primitive = prim_nums[0] as i32,
            _ => {
                self.priv_union.primitive_indices_offset = primitive_indices.len() as i32;
                for item in prim_nums.iter().take(np) {
                    primitive_indices.push(*item as i32);
                }
            }
        }
    }
    pub fn init_interior(&mut self, axis: i32, ac: i32, s: Float) {
        self.priv_union.split = s;
        self.pub_union.flags = axis;
        let above_child: i32;
        unsafe {
            above_child = self.pub_union.above_child;
        }
        self.pub_union.above_child = above_child | (ac << 2);
    }
    pub fn split_pos(&self) -> Float {
        let split: Float;
        unsafe {
            split = self.priv_union.split;
        }
        split
    }
    pub fn n_primitives(&self) -> i32 {
        let n_prims: i32;
        unsafe {
            n_prims = self.pub_union.n_prims;
        }
        n_prims >> 2
    }
    pub fn split_axis(&self) -> i32 {
        let flags: i32;
        unsafe {
            flags = self.pub_union.flags;
        }
        flags & 3
    }
    pub fn is_leaf(&self) -> bool {
        let flags: i32;
        unsafe {
            flags = self.pub_union.flags;
        }
        (flags & 3) == 3
    }
    pub fn above_child(&self) -> i32 {
        let above_child: i32;
        unsafe {
            above_child = self.pub_union.above_child;
        }
        above_child >> 2
    }
}

#[derive(Debug, PartialEq, PartialOrd)]
pub enum EdgeType {
    Start = 0,
    End = 1,
}

#[derive(Debug)]
pub struct BoundEdge {
    pub t: Float,
    pub prim_num: usize,
    pub edge_type: EdgeType,
}

impl BoundEdge {
    pub fn new(t: Float, prim_num: usize, starting: bool) -> Self {
        let edge_type = if starting {
            EdgeType::Start
        } else {
            EdgeType::End
        };
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
    pub traversal_cost: i32,
    pub max_prims: i32,
    pub empty_bonus: Float,
    pub primitives: Vec<Arc<Primitive>>,
    pub primitive_indices: Vec<i32>,
    pub nodes: Vec<KdAccelNode>,
    pub n_alloced_nodes: i32,
    pub next_free_node: i32,
    pub bounds: Bounds3f,
}

impl KdTreeAccel {
    pub fn new(
        p: Vec<Arc<Primitive>>,
        isect_cost: i32,
        traversal_cost: i32,
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
        for item in p.iter().take(p_len) {
            let b: Bounds3f = item.world_bound();
            bounds = bnd3_union_bnd3f(&bounds, &b);
            prim_bounds.push(b);
        }
        // allocate working memory for kd-tree construction
        let mut edges: [Vec<BoundEdge>; 3] = [
            Vec::with_capacity(2 * p_len),
            Vec::with_capacity(2 * p_len),
            Vec::with_capacity(2 * p_len),
        ];
        let mut prims0: Vec<usize> = Vec::with_capacity(p_len);
        let mut prims1: Vec<usize> = Vec::with_capacity((max_depth + 1) as usize * p_len);
        for _i in 0..((max_depth + 1) as usize * p_len) {
            prims1.push(0_usize);
        }
        // initialize _prim_nums_ for kd-tree construction
        let mut prim_nums: Vec<usize> = Vec::with_capacity(p_len);
        for i in 0..p_len {
            prims0.push(0_usize);
            prim_nums.push(i);
            // init all three edges Vecs
            edges[0].push(BoundEdge::default());
            edges[0].push(BoundEdge::default());
            edges[1].push(BoundEdge::default());
            edges[1].push(BoundEdge::default());
            edges[2].push(BoundEdge::default());
            edges[2].push(BoundEdge::default());
        }
        // start recursive construction of kd-tree
        let mut kd_tree: KdTreeAccel = KdTreeAccel {
            isect_cost,
            traversal_cost,
            max_prims,
            empty_bonus,
            primitives: p,
            primitive_indices: Vec::new(),
            nodes: Vec::new(),
            n_alloced_nodes,
            next_free_node,
            bounds,
        };
        KdTreeAccel::build_tree(
            &mut kd_tree,
            0 as i32,
            &bounds,
            &prim_bounds,
            &prim_nums[..],
            p_len,
            max_depth,
            &mut edges,
            &mut prims0[..],
            &mut prims1[..],
            0, // bad_refines
        );
        kd_tree
    }
    pub fn create(prims: Vec<Arc<Primitive>>, ps: &ParamSet) -> Primitive {
        let isect_cost: i32 = ps.find_one_int("intersectcost", 80);
        let trav_cost: i32 = ps.find_one_int("traversalcost", 1);
        let empty_bonus: Float = ps.find_one_float("emptybonus", 0.5 as Float);
        let max_prims: i32 = ps.find_one_int("maxprims", 1);
        let max_depth: i32 = ps.find_one_int("maxdepth", -1);
        Primitive::KdTree(Box::new(KdTreeAccel::new(
            prims,
            isect_cost,
            trav_cost,
            empty_bonus,
            max_prims,
            max_depth,
        )))
    }
    pub fn build_tree(
        &mut self,
        node_num: i32,
        node_bounds: &Bounds3f,
        all_prim_bounds: &[Bounds3f],
        prim_nums: &[usize],
        n_primitives: usize,
        depth: i32,
        edges: &mut [Vec<BoundEdge>; 3],
        prims0: &mut [usize],
        prims1: &mut [usize],
        bad_refines: i32,
    ) {
        let mut bad_refines: i32 = bad_refines;
        assert_eq!(node_num, self.next_free_node);
        if self.next_free_node == self.n_alloced_nodes {
            let n_new_alloc_nodes: i32 = std::cmp::max(2 * self.n_alloced_nodes, 512);
            if self.n_alloced_nodes > 0 {
                self.nodes
                    .resize_with(n_new_alloc_nodes as usize, || KdAccelNode {
                        priv_union: PrivateUnion {
                            one_primitive: 0_i32,
                        },
                        pub_union: PublicUnion { flags: 0_i32 },
                    });
            } else {
                let mut n: Vec<KdAccelNode> = Vec::with_capacity(n_new_alloc_nodes as usize);
                for _i in 0..n_new_alloc_nodes as usize {
                    n.push(KdAccelNode {
                        priv_union: PrivateUnion {
                            one_primitive: 0_i32,
                        },
                        pub_union: PublicUnion { flags: 0_i32 },
                    });
                }
                self.nodes = n;
            }
            self.n_alloced_nodes = n_new_alloc_nodes;
        }
        self.next_free_node += 1;
        // initialize leaf node if termination criteria met
        if n_primitives <= self.max_prims as usize || depth == 0 {
            self.nodes[node_num as usize].init_leaf(
                prim_nums,
                n_primitives,
                &mut self.primitive_indices,
            );
            return;
        }
        // choose split axis position for interior node
        let mut best_axis: i32 = -1;
        let mut best_axis_i: XYZEnum = match best_axis {
            0 => XYZEnum::X,
            1 => XYZEnum::Y,
            _ => XYZEnum::Z,
        };
        let mut best_offset: i32 = -1;
        let mut best_cost: Float = std::f32::INFINITY;
        let old_cost: Float = self.isect_cost as Float * n_primitives as Float;
        let total_sa: Float = node_bounds.surface_area();
        let inv_total_sa: Float = 1.0 as Float / total_sa;
        let d: Vector3f = node_bounds.p_max - node_bounds.p_min;
        // choose which axis to split along
        let mut axis: u8 = node_bounds.maximum_extent();
        let mut axis_i: XYZEnum = match axis {
            0 => XYZEnum::X,
            1 => XYZEnum::Y,
            _ => XYZEnum::Z,
        };
        let mut retries: u8 = 0;
        // avoid 'goto retrySplit;'
        loop {
            // trim edges to 2 * n_primitives
            edges[axis as usize].resize_with(2 * n_primitives, BoundEdge::default);
            // initialize edges for _axis_
            for (i, item) in prim_nums.iter().enumerate().take(n_primitives) {
                let pn: usize = *item;
                let bounds: &Bounds3f = &all_prim_bounds[pn];
                edges[axis as usize][2 * i] = BoundEdge::new(bounds.p_min[axis_i], pn, true);
                edges[axis as usize][2 * i + 1] = BoundEdge::new(bounds.p_max[axis_i], pn, false);
            }
            // sort _edges_ for _axis_
            edges[axis as usize].sort_unstable_by(|e0, e1| {
                if e0.t == e1.t {
                    e0.edge_type.partial_cmp(&e1.edge_type).unwrap()
                } else {
                    e0.t.partial_cmp(&e1.t).unwrap()
                }
            });
            // for i in 0..n_primitives {
            //     println!("{:?}", edges[axis as usize][2 * i]);
            //     println!("{:?}", edges[axis as usize][2 * i + 1]);
            // }

            // compute cost of all splits for _axis_ to find best
            let mut n_below: usize = 0;
            let mut n_above: usize = n_primitives;
            for i in 0..(2 * n_primitives) {
                if edges[axis as usize][i].edge_type == EdgeType::End {
                    n_above -= 1;
                }
                let edge_t: Float = edges[axis as usize][i].t;
                if edge_t > node_bounds.p_min[axis_i] && edge_t < node_bounds.p_max[axis_i] {
                    // compute cost for split at _i_th edge

                    // compute child surface areas for split at _edge_t_
                    let other_axis_0: u8 = (axis + 1) % 3;
                    let other_axis_1: u8 = (axis + 2) % 3;
                    let other_axis_0_i: XYZEnum = match other_axis_0 {
                        0 => XYZEnum::X,
                        1 => XYZEnum::Y,
                        _ => XYZEnum::Z,
                    };
                    let other_axis_1_i: XYZEnum = match other_axis_1 {
                        0 => XYZEnum::X,
                        1 => XYZEnum::Y,
                        _ => XYZEnum::Z,
                    };
                    let below_sa: Float = 2.0 as Float
                        * (d[other_axis_0_i] * d[other_axis_1_i]
                            + (edge_t - node_bounds.p_min[axis_i])
                                * (d[other_axis_0_i] + d[other_axis_1_i]));
                    let above_sa: Float = 2.0 as Float
                        * (d[other_axis_0_i] * d[other_axis_1_i]
                            + (node_bounds.p_max[axis_i] - edge_t)
                                * (d[other_axis_0_i] + d[other_axis_1_i]));
                    let p_below: Float = below_sa * inv_total_sa;
                    let p_above: Float = above_sa * inv_total_sa;
                    let eb = if n_above == 0 || n_below == 0 {
                        self.empty_bonus
                    } else {
                        0.0 as Float
                    };
                    let cost: Float = self.traversal_cost as Float
                        + self.isect_cost as Float
                            * (1.0 as Float - eb)
                            * (p_below * n_below as Float + p_above * n_above as Float);
                    // update best split if this is lowest cost so far
                    if cost < best_cost {
                        best_cost = cost;
                        best_axis = axis as i32;
                        best_axis_i = match best_axis {
                            0 => XYZEnum::X,
                            1 => XYZEnum::Y,
                            _ => XYZEnum::Z,
                        };
                        best_offset = i as i32;
                    }
                }
                if edges[axis as usize][i].edge_type == EdgeType::Start {
                    n_below += 1;
                }
            }
            assert!(
                n_below == n_primitives && n_above == 0,
                "{} == {}? && {} == 0?",
                n_below,
                n_primitives,
                n_above
            );
            // create leaf if no good splits were found
            if best_axis == -1 && retries < 2 {
                retries += 1;
                axis = (axis + 1) % 3;
                axis_i = match axis {
                    0 => XYZEnum::X,
                    1 => XYZEnum::Y,
                    _ => XYZEnum::Z,
                };
            // goto retrySplit;
            } else {
                break;
            }
        }
        if best_cost > old_cost {
            bad_refines += 1;
        }
        if (best_cost > 4.0 as Float * old_cost && n_primitives < 16)
            || best_axis == -1
            || bad_refines == 3
        {
            self.nodes[node_num as usize].init_leaf(
                prim_nums,
                n_primitives,
                &mut self.primitive_indices,
            );
            return;
        }
        // classify primitives with respect to split
        let mut n0: usize = 0;
        let mut n1: usize = 0;
        for i in 0..best_offset as usize {
            if edges[best_axis as usize][i].edge_type == EdgeType::Start {
                prims0[n0] = edges[best_axis as usize][i].prim_num;
                n0 += 1;
            }
        }
        for i in ((best_offset + 1) as usize)..(2 * n_primitives) {
            if edges[best_axis as usize][i].edge_type == EdgeType::End {
                prims1[n1] = edges[best_axis as usize][i].prim_num;
                n1 += 1;
            }
        }
        // recursively initialize children nodes
        let t_split: Float = edges[best_axis as usize][best_offset as usize].t;
        let mut bounds0: Bounds3f = *node_bounds;
        let mut bounds1: Bounds3f = *node_bounds;
        bounds0.p_max[best_axis_i] = t_split;
        bounds1.p_min[best_axis_i] = t_split;
        // copy prims0
        let mut prim_nums: Vec<usize> = Vec::with_capacity(prims0.len());
        for i in 0..prims0.len() {
            prim_nums.push(prims0[i]);
        }
        self.build_tree(
            node_num + 1,
            &bounds0,
            all_prim_bounds,
            &prim_nums[..],
            n0,
            depth - 1,
            edges,
            prims0,
            &mut prims1[n_primitives..],
            bad_refines,
        );
        let above_child: i32 = self.next_free_node;
        self.nodes[node_num as usize].init_interior(best_axis, above_child, t_split);
        // copy prims1
        let mut prim_nums: Vec<usize> = Vec::with_capacity(prims1.len());
        for i in 0..prims1.len() {
            prim_nums.push(prims1[i]);
        }
        self.build_tree(
            above_child,
            &bounds1,
            all_prim_bounds,
            &prim_nums[..],
            n1,
            depth - 1,
            edges,
            prims0,
            &mut prims1[n_primitives..],
            bad_refines,
        );
    }
    // Primitive
    pub fn world_bound(&self) -> Bounds3f {
        self.bounds
    }
    pub fn intersect(&self, ray: &mut Ray, isect: &mut SurfaceInteraction) -> bool {
        // TODO: ProfilePhase p(Prof::AccelIntersect);
        if self.nodes.is_empty() {
            return false;
        }
        // compute initial parametric range of ray inside kd-tree extent
        let mut t_min: Float = 0.0;
        let mut t_max: Float = 0.0;
        if !self.bounds.intersect_b(&ray, &mut t_min, &mut t_max) {
            return false;
        }
        // prepare to traverse kd-tree for ray
        let inv_dir: Vector3f = Vector3f {
            x: 1.0 / ray.d.x,
            y: 1.0 / ray.d.y,
            z: 1.0 / ray.d.z,
        };
        let mut todo: [KdToDo; MAX_TODO] = [KdToDo::default(); MAX_TODO];
        let mut todo_pos: usize = 0;
        // traverse kd-tree nodes in order for ray
        let mut hit: bool = false;
        let mut node_idx: usize = 0;
        let mut node_opt: Option<&KdAccelNode> = self.nodes.get(node_idx);
        while let Some(node) = node_opt {
            // bail out if we found a hit closer than the current node
            if ray.t_max < t_min {
                break;
            }
            if !node.is_leaf() {
                // process kd-tree interior node

                // compute parametric distance along ray to split plane
                let axis: u8 = node.split_axis() as u8;
                let axis_i: XYZEnum = match axis {
                    0 => XYZEnum::X,
                    1 => XYZEnum::Y,
                    _ => XYZEnum::Z,
                };
                let t_plane: Float = (node.split_pos() - ray.o[axis_i]) * inv_dir[axis_i];
                // get node children pointers for ray
                let below_first: bool = (ray.o[axis_i] < node.split_pos())
                    || (ray.o[axis_i] == node.split_pos() && ray.d[axis_i] <= 0.0 as Float);
                let first_child: Option<&KdAccelNode>;
                let second_child: Option<&KdAccelNode>;
                let first_idx: usize;
                let second_idx: usize;
                if below_first {
                    first_idx = node_idx + 1;
                    first_child = self.nodes.get(first_idx);
                    second_idx = node.above_child() as usize;
                    second_child = self.nodes.get(second_idx);
                } else {
                    first_idx = node.above_child() as usize;
                    first_child = self.nodes.get(first_idx);
                    second_idx = node_idx + 1;
                    second_child = self.nodes.get(second_idx);
                }
                // advance to next child node, possibly enqueue other child
                if t_plane > t_max || t_plane <= 0.0 as Float {
                    node_opt = first_child;
                    node_idx = first_idx;
                } else if t_plane < t_min {
                    node_opt = second_child;
                    node_idx = second_idx;
                } else {
                    // enqueue _second_child_ in todo list
                    todo[todo_pos].node = second_child;
                    todo[todo_pos].idx = second_idx;
                    todo[todo_pos].t_min = t_plane;
                    todo[todo_pos].t_max = t_max;
                    todo_pos += 1;
                    node_opt = first_child;
                    node_idx = first_idx;
                    t_max = t_plane;
                }
            } else {
                // check for intersections inside leaf node
                let n_primitives: i32 = node.n_primitives();
                if n_primitives == 1 {
                    let one_primitive: i32;
                    unsafe {
                        one_primitive = node.priv_union.one_primitive;
                    }
                    let p: &Arc<Primitive> = &self.primitives[one_primitive as usize];
                    // check one primitive inside leaf node
                    if p.intersect(ray, isect) {
                        hit = true;
                    }
                } else {
                    for i in 0..n_primitives {
                        let primitive_indices_offset: i32;
                        unsafe {
                            primitive_indices_offset = node.priv_union.primitive_indices_offset;
                        }
                        let index: usize = self.primitive_indices
                            [(primitive_indices_offset + i) as usize]
                            as usize;
                        let p: &Arc<Primitive> = &self.primitives[index];
                        // check one primitive inside leaf node
                        if p.intersect(ray, isect) {
                            hit = true;
                        }
                    }
                }
                // grab next node to process from todo list
                if todo_pos > 0 {
                    todo_pos -= 1;
                    node_opt = todo[todo_pos].node;
                    node_idx = todo[todo_pos].idx;
                    t_min = todo[todo_pos].t_min;
                    t_max = todo[todo_pos].t_max;
                } else {
                    break;
                }
            }
        }
        hit
    }
    pub fn intersect_p(&self, ray: &Ray) -> bool {
        // TODO: ProfilePhase p(Prof::AccelIntersectP);
        if self.nodes.is_empty() {
            return false;
        }
        // compute initial parametric range of ray inside kd-tree extent
        let mut t_min: Float = 0.0;
        let mut t_max: Float = 0.0;
        if !self.bounds.intersect_b(&ray, &mut t_min, &mut t_max) {
            return false;
        }
        // prepare to traverse kd-tree for ray
        let inv_dir: Vector3f = Vector3f {
            x: 1.0 / ray.d.x,
            y: 1.0 / ray.d.y,
            z: 1.0 / ray.d.z,
        };
        let mut todo: [KdToDo; MAX_TODO] = [KdToDo::default(); MAX_TODO];
        let mut todo_pos: usize = 0;
        let mut node_idx: usize = 0;
        let mut node_opt: Option<&KdAccelNode> = self.nodes.get(node_idx);
        while let Some(node) = node_opt {
            if node.is_leaf() {
                // check for shadow ray intersections inside leaf node
                let n_primitives: i32 = node.n_primitives();
                if n_primitives == 1 {
                    let one_primitive: i32;
                    unsafe {
                        one_primitive = node.priv_union.one_primitive;
                    }
                    let p: &Arc<Primitive> = &self.primitives[one_primitive as usize];
                    if p.intersect_p(ray) {
                        return true;
                    }
                } else {
                    for i in 0..n_primitives {
                        let primitive_indices_offset: i32;
                        unsafe {
                            primitive_indices_offset = node.priv_union.primitive_indices_offset;
                        }
                        let primitive_index: usize = self.primitive_indices
                            [(primitive_indices_offset + i) as usize]
                            as usize;
                        let prim: &Arc<Primitive> = &self.primitives[primitive_index];
                        if prim.intersect_p(ray) {
                            return true;
                        }
                    }
                }
                // grab next node to process from todo list
                if todo_pos > 0 {
                    todo_pos -= 1;
                    node_opt = todo[todo_pos].node;
                    node_idx = todo[todo_pos].idx;
                    t_min = todo[todo_pos].t_min;
                    t_max = todo[todo_pos].t_max;
                } else {
                    break;
                }
            } else {
                // process kd-tree interior node

                // compute parametric distance along ray to split plane
                let axis: u8 = node.split_axis() as u8;
                let axis_i: XYZEnum = match axis {
                    0 => XYZEnum::X,
                    1 => XYZEnum::Y,
                    _ => XYZEnum::Z,
                };
                let t_plane: Float = (node.split_pos() - ray.o[axis_i]) * inv_dir[axis_i];
                // get node children pointers for ray
                let below_first: bool = (ray.o[axis_i] < node.split_pos())
                    || (ray.o[axis_i] == node.split_pos() && ray.d[axis_i] <= 0.0 as Float);
                let first_child: Option<&KdAccelNode>;
                let second_child: Option<&KdAccelNode>;
                let first_idx: usize;
                let second_idx: usize;
                if below_first {
                    first_idx = node_idx + 1;
                    first_child = self.nodes.get(first_idx);
                    second_idx = node.above_child() as usize;
                    second_child = self.nodes.get(second_idx);
                } else {
                    first_idx = node.above_child() as usize;
                    first_child = self.nodes.get(first_idx);
                    second_idx = node_idx + 1;
                    second_child = self.nodes.get(second_idx);
                }
                // advance to next child node, possibly enqueue other child
                if t_plane > t_max || t_plane <= 0.0 as Float {
                    node_opt = first_child;
                    node_idx = first_idx;
                } else if t_plane < t_min {
                    node_opt = second_child;
                    node_idx = second_idx;
                } else {
                    // enqueue _second_child_ in todo list
                    todo[todo_pos].node = second_child;
                    todo[todo_pos].idx = second_idx;
                    todo[todo_pos].t_min = t_plane;
                    todo[todo_pos].t_max = t_max;
                    todo_pos += 1;
                    node_opt = first_child;
                    node_idx = first_idx;
                    t_max = t_plane;
                }
            }
        }
        false
    }
    pub fn get_material(&self) -> Option<Arc<Material>> {
        None
    }
    pub fn get_area_light(&self) -> Option<Arc<Light>> {
        None
    }
}
