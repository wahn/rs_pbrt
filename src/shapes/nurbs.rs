// std
use std;
// pbrt
use core::geometry::{Point3f, Vector3f};
use core::pbrt::Float;

// see nurbs.cpp

pub fn knot_offset(knot: &Vec<Float>, order: i32, np: i32, t: Float) -> usize {
    let first_knot: usize = (order - 1_i32) as usize;
    let mut knot_offset: usize = first_knot;
    while t > knot[knot_offset + 1_usize] {
        knot_offset += 1;
    }
    // np == last_knot
    assert!(knot_offset < np as usize);
    assert!(t >= knot[knot_offset] && t <= knot[knot_offset + 1]);
    knot_offset
}

#[derive(Debug, Default, Copy, Clone)]
pub struct Homogeneous3 {
    pub x: Float,
    pub y: Float,
    pub z: Float,
    pub w: Float,
}

pub fn nurbs_evaluate(
    order: i32,
    knot: &Vec<Float>,
    cp: &Vec<Homogeneous3>,
    cp_start: i32,
    np: i32,
    cp_stride: i32,
    t: Float,
    deriv_opt: Option<&mut Vector3f>,
) -> Homogeneous3 {
    let mut alpha: Float;
    let knot_offset: usize = knot_offset(knot, order, np, t);
    let cp_offset: usize = knot_offset + 1 - order as usize;
    assert!(cp_offset < np as usize);
    let mut cp_work: Vec<Homogeneous3> = Vec::with_capacity(order as usize);
    for i in 0_usize..order as usize {
        cp_work.push(cp[(cp_start + (cp_offset + i) as i32 * cp_stride) as usize]);
    }
    for i in 0..(order - 2) {
        for j in 0..(order - 1 - i) {
            alpha = (knot[knot_offset + 1 + j as usize] - t)
                / (knot[knot_offset + 1 + j as usize]
                    - knot[(knot_offset as i32 + (j + 2 + i - order)) as usize]);
            assert!(alpha >= 0.0 as Float && alpha <= 1.0 as Float);
            let one_minus_alpha: Float = 1.0 as Float - alpha;
            cp_work[j as usize].x =
                cp_work[j as usize].x * alpha + cp_work[(j + 1) as usize].x * one_minus_alpha;
            cp_work[j as usize].y =
                cp_work[j as usize].y * alpha + cp_work[(j + 1) as usize].y * one_minus_alpha;
            cp_work[j as usize].z =
                cp_work[j as usize].z * alpha + cp_work[(j + 1) as usize].z * one_minus_alpha;
            cp_work[j as usize].w =
                cp_work[j as usize].w * alpha + cp_work[(j + 1) as usize].w * one_minus_alpha;
        }
    }
    alpha = (knot[knot_offset + 1] - t) / (knot[knot_offset + 1] - knot[knot_offset + 0]);
    assert!(alpha >= 0.0 as Float && alpha <= 1.0 as Float);
    let one_minus_alpha: Float = 1.0 as Float - alpha;
    let val: Homogeneous3 = Homogeneous3 {
        x: cp_work[0].x * alpha + cp_work[1].x * one_minus_alpha,
        y: cp_work[0].y * alpha + cp_work[1].y * one_minus_alpha,
        z: cp_work[0].z * alpha + cp_work[1].z * one_minus_alpha,
        w: cp_work[0].w * alpha + cp_work[1].w * one_minus_alpha,
    };
    if let Some(deriv) = deriv_opt {
        let factor: Float = (order - 1) as Float / (knot[knot_offset + 1] - knot[knot_offset + 0]);
        let delta: Homogeneous3 = Homogeneous3 {
            x: (cp_work[1].x - cp_work[0].x) * factor,
            y: (cp_work[1].y - cp_work[0].y) * factor,
            z: (cp_work[1].z - cp_work[0].z) * factor,
            w: (cp_work[1].w - cp_work[0].w) * factor,
        };
        deriv.x = delta.x / val.w - (val.x * delta.w / (val.w * val.w));
        deriv.y = delta.y / val.w - (val.y * delta.w / (val.w * val.w));
        deriv.z = delta.z / val.w - (val.z * delta.w / (val.w * val.w));
    }
    val
}

pub fn nurbs_evaluate_surface(
    u_order: i32,
    u_knot: &Vec<Float>,
    ucp: i32,
    u: Float,
    v_order: i32,
    v_knot: &Vec<Float>,
    vcp: i32,
    v: Float,
    cp: &Vec<Homogeneous3>,
    dpdu_opt: Option<&mut Vector3f>,
    dpdv_opt: Option<&mut Vector3f>,
) -> Point3f {
    let mut iso: Vec<Homogeneous3> = Vec::with_capacity(std::cmp::max(u_order, v_order) as usize);
    let u_offset: usize = knot_offset(u_knot, u_order, ucp, u);
    let u_first_cp: usize = u_offset + 1 - u_order as usize;
    assert!(u_first_cp + u_order as usize - 1 < ucp as usize);
    for i in 0..u_order {
        iso.push(nurbs_evaluate(
            v_order,
            v_knot,
            &cp,
            u_first_cp as i32 + i,
            vcp,
            ucp,
            v,
            None,
        ));
    }
    let v_offset: usize = knot_offset(v_knot, v_order, vcp, v);
    let v_first_cp: usize = v_offset + 1 - v_order as usize;
    assert!(v_first_cp + v_order as usize - 1 < vcp as usize);
    let p: Homogeneous3 = nurbs_evaluate(
        u_order, u_knot, &iso, -(u_first_cp as i32),
        ucp, 1, u, dpdu_opt,
    );
    if let Some(dpdv) = dpdv_opt {
        for i in 0..v_order {
            iso[i as usize] = nurbs_evaluate(
                u_order,
                u_knot,
                &cp,
                (v_first_cp as i32 + i) * ucp,
                ucp,
                1,
                u,
                None,
            );
        }
        nurbs_evaluate(
            v_order, v_knot, &iso, -(v_first_cp as i32),
            vcp, 1, v, Some(dpdv),
        );
    }
    Point3f {
        x: p.x / p.w,
        y: p.y / p.w,
        z: p.z / p.w,
    }
}
