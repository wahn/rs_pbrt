//! In general, transformations make it possible to work in the most
//! convenient coordinate space.
//!
//! ## 4 x 4 Matrices
//!
//! The **Matrix4x4** structure provides a low-level representation of
//! 4 x 4 matrices. It is an integral part of the **Transform** class.
//!
//! ## Transformations
//!
//! In general a transformation is a mapping from points to points and
//! from vectors to vectors. When a new **Transform** is created, it
//! defaults to the *identity transformation* - the transformation
//! that maps each point and each vector to itself.
//!
//! ### Translations
//!
//! One of the simplest transformations is the translation
//! transformation. Translations only affect points, leaving vectors
//! unchanged.
//!
//! ### Scaling
//!
//! Another basic transformations is the scale transformation. We can
//! differentiate between **uniform** scaling, where all three scale
//! factors have the same value, and **nonuniform** scaling, where
//! they may have different values.
//!
//! ### X, Y, And Z Axis Rotations
//!
//! Another useful type of transformation is the rotation
//! transformation.
//!
//! ### Rotation Around an Arbitrary Axis
//!
//! We also provide a routine to compute the transformation that
//! represents rotation around an arbitrary axis.
//!
//! ### The Look-At Transformation
//!
//! The *look-at* transformation is particularly useful for placing a
//! camera in the scene. The caller specifies the desired position of
//! the camera, a point the camera is looking at, and an "up" vector
//! that orients the camera along the viewing direction implied by the
//! first two parameters. All of these values are given in world space
//! coordinates. The look-at construction then gives a transformation
//! between camera space and world space.
//!
//! ## Quaternions
//!
//! Quaternions were originally invented by Sir William Hamilton in
//! 1843 as a generalization of complex numbers (2 dimensions) to four
//! dimensions.

// std
use std::cell::Cell;
use std::f32::consts::PI;
use std::ops::{Add, Mul};
// pbrt
use crate::core::geometry::{
    bnd3_union_bnd3f, bnd3_union_pnt3f, nrm_faceforward_nrm, vec3_cross_vec3, vec3_dot_vec3f,
};
use crate::core::geometry::{Bounds3f, Normal3f, Point3f, Ray, RayDifferential, Vector3f};
use crate::core::interaction::SurfaceInteraction;
use crate::core::pbrt::Float;
use crate::core::pbrt::{clamp_t, gamma, lerp, radians};
use crate::core::quaternion::Quaternion;
use crate::core::quaternion::{quat_dot_quat, quat_normalize, quat_slerp};

// see transform.h

#[derive(Debug, Copy, Clone)]
pub struct Matrix4x4 {
    pub m: [[Float; 4]; 4],
}

impl Default for Matrix4x4 {
    fn default() -> Self {
        Matrix4x4 {
            m: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        }
    }
}

impl Matrix4x4 {
    pub fn new(
        t00: Float,
        t01: Float,
        t02: Float,
        t03: Float,
        t10: Float,
        t11: Float,
        t12: Float,
        t13: Float,
        t20: Float,
        t21: Float,
        t22: Float,
        t23: Float,
        t30: Float,
        t31: Float,
        t32: Float,
        t33: Float,
    ) -> Self {
        Matrix4x4 {
            m: [
                [t00, t01, t02, t03],
                [t10, t11, t12, t13],
                [t20, t21, t22, t23],
                [t30, t31, t32, t33],
            ],
        }
    }
    pub fn transpose(m: &Matrix4x4) -> Matrix4x4 {
        Matrix4x4 {
            m: [
                [m.m[0][0], m.m[1][0], m.m[2][0], m.m[3][0]],
                [m.m[0][1], m.m[1][1], m.m[2][1], m.m[3][1]],
                [m.m[0][2], m.m[1][2], m.m[2][2], m.m[3][2]],
                [m.m[0][3], m.m[1][3], m.m[2][3], m.m[3][3]],
            ],
        }
    }
    pub fn inverse(m: &Matrix4x4) -> Matrix4x4 {
        let mut indxc = vec![0; 4];
        let mut indxr = vec![0; 4];
        let mut ipiv = vec![0; 4];
        let mut minv: Matrix4x4 = Matrix4x4::new(
            m.m[0][0], m.m[0][1], m.m[0][2], m.m[0][3], m.m[1][0], m.m[1][1], m.m[1][2], m.m[1][3],
            m.m[2][0], m.m[2][1], m.m[2][2], m.m[2][3], m.m[3][0], m.m[3][1], m.m[3][2], m.m[3][3],
        );
        for i in 0..4 {
            let mut irow = 0;
            let mut icol = 0;
            let mut big: Float = 0.0;
            // choose pivot
            for j in 0..4 {
                if ipiv[j] != 1 {
                    for (k, item) in ipiv.iter().enumerate().take(4) {
                        if *item == 0 {
                            let abs: Float = (minv.m[j][k]).abs();
                            if abs >= big {
                                big = abs;
                                irow = j;
                                icol = k;
                            }
                        } else if *item > 1 {
                            println!("Singular matrix in MatrixInvert");
                        }
                    }
                }
            }
            ipiv[icol] += 1;
            // swap rows _irow_ and _icol_ for pivot
            if irow != icol {
                for k in 0..4 {
                    // C++: std::swap(minv[irow][k], minv[icol][k]);
                    let swap = minv.m[irow][k];
                    minv.m[irow][k] = minv.m[icol][k];
                    minv.m[icol][k] = swap;
                }
            }
            indxr[i] = irow;
            indxc[i] = icol;
            if minv.m[icol][icol] == 0.0 {
                println!("Singular matrix in MatrixInvert");
            }
            // set $m[icol][icol]$ to one by scaling row _icol_ appropriately
            let pivinv: Float = 1.0 / minv.m[icol][icol];
            minv.m[icol][icol] = 1.0;
            for j in 0..4 {
                minv.m[icol][j] *= pivinv;
            }
            // subtract this row from others to zero out their columns
            for j in 0..4 {
                if j != icol {
                    let save: Float = minv.m[j][icol];
                    minv.m[j][icol] = 0.0;
                    for k in 0..4 {
                        minv.m[j][k] -= minv.m[icol][k] * save;
                    }
                }
            }
        }
        // swap columns to reflect permutation
        for i in 0..4 {
            let j = 3 - i;
            if indxr[j] != indxc[j] {
                for k in 0..4 {
                    // C++: std::swap(minv[k][indxr[j]], minv[k][indxc[j]]);
                    minv.m[k].swap(indxr[j], indxc[j])
                }
            }
        }
        minv
    }
}

impl PartialEq for Matrix4x4 {
    fn eq(&self, rhs: &Matrix4x4) -> bool {
        for i in 0..4 {
            for j in 0..4 {
                if self.m[i][j] != rhs.m[i][j] {
                    return false;
                }
            }
        }
        true
    }
}

// see transform.cpp

/// Finds the closed-form solution of a 2x2 linear system.
pub fn solve_linear_system_2x2(
    a: [[Float; 2]; 2],
    b: [Float; 2],
    x0: &mut Float,
    x1: &mut Float,
) -> bool {
    let det: Float = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    if det.abs() < 1e-10 as Float {
        return false;
    }
    *x0 = (a[1][1] * b[0] - a[0][1] * b[1]) / det;
    *x1 = (a[0][0] * b[1] - a[1][0] * b[0]) / det;
    if (*x0).is_nan() || (*x1).is_nan() {
        return false;
    }
    true
}

/// The product of two matrices.
pub fn mtx_mul(m1: &Matrix4x4, m2: &Matrix4x4) -> Matrix4x4 {
    let mut r: Matrix4x4 = Matrix4x4::default();
    for i in 0..4 {
        for j in 0..4 {
            r.m[i][j] = m1.m[i][0] * m2.m[0][j]
                + m1.m[i][1] * m2.m[1][j]
                + m1.m[i][2] * m2.m[2][j]
                + m1.m[i][3] * m2.m[3][j];
        }
    }
    r
}

#[derive(Debug, Copy, Clone)]
pub struct Transform {
    pub m: Matrix4x4,
    pub m_inv: Matrix4x4,
}

impl Default for Transform {
    fn default() -> Self {
        Transform {
            m: Matrix4x4::default(),
            m_inv: Matrix4x4::default(),
        }
    }
}

impl Transform {
    pub fn new(
        t00: Float,
        t01: Float,
        t02: Float,
        t03: Float,
        t10: Float,
        t11: Float,
        t12: Float,
        t13: Float,
        t20: Float,
        t21: Float,
        t22: Float,
        t23: Float,
        t30: Float,
        t31: Float,
        t32: Float,
        t33: Float,
    ) -> Self {
        Transform {
            m: Matrix4x4::new(
                t00, t01, t02, t03, t10, t11, t12, t13, t20, t21, t22, t23, t30, t31, t32, t33,
            ),
            m_inv: Matrix4x4::inverse(&Matrix4x4::new(
                t00, t01, t02, t03, t10, t11, t12, t13, t20, t21, t22, t23, t30, t31, t32, t33,
            )),
        }
    }
    pub fn inverse(t: &Transform) -> Transform {
        Transform {
            m: t.m_inv,
            m_inv: t.m,
        }
    }
    pub fn is_identity(&self) -> bool {
        self.m.m[0][0] == 1.0 as Float
            && self.m.m[0][1] == 0.0 as Float
            && self.m.m[0][2] == 0.0 as Float
            && self.m.m[0][3] == 0.0 as Float
            && self.m.m[1][0] == 0.0 as Float
            && self.m.m[1][1] == 1.0 as Float
            && self.m.m[1][2] == 0.0 as Float
            && self.m.m[1][3] == 0.0 as Float
            && self.m.m[2][0] == 0.0 as Float
            && self.m.m[2][1] == 0.0 as Float
            && self.m.m[2][2] == 1.0 as Float
            && self.m.m[2][3] == 0.0 as Float
            && self.m.m[3][0] == 0.0 as Float
            && self.m.m[3][1] == 0.0 as Float
            && self.m.m[3][2] == 0.0 as Float
            && self.m.m[3][3] == 1.0 as Float
    }
    pub fn swaps_handedness(&self) -> bool {
        let det: Float = self.m.m[0][0]
            * (self.m.m[1][1] * self.m.m[2][2] - self.m.m[1][2] * self.m.m[2][1])
            - self.m.m[0][1] * (self.m.m[1][0] * self.m.m[2][2] - self.m.m[1][2] * self.m.m[2][0])
            + self.m.m[0][2] * (self.m.m[1][0] * self.m.m[2][1] - self.m.m[1][1] * self.m.m[2][0]);
        det < 0.0 as Float
    }
    pub fn translate(delta: &Vector3f) -> Transform {
        Transform {
            m: Matrix4x4::new(
                1.0, 0.0, 0.0, delta.x, 0.0, 1.0, 0.0, delta.y, 0.0, 0.0, 1.0, delta.z, 0.0, 0.0,
                0.0, 1.0,
            ),
            m_inv: Matrix4x4::new(
                1.0, 0.0, 0.0, -delta.x, 0.0, 1.0, 0.0, -delta.y, 0.0, 0.0, 1.0, -delta.z, 0.0,
                0.0, 0.0, 1.0,
            ),
        }
    }
    pub fn scale(x: Float, y: Float, z: Float) -> Transform {
        Transform {
            m: Matrix4x4::new(
                x, 0.0, 0.0, 0.0, 0.0, y, 0.0, 0.0, 0.0, 0.0, z, 0.0, 0.0, 0.0, 0.0, 1.0,
            ),
            m_inv: Matrix4x4::new(
                1.0 / x,
                0.0,
                0.0,
                0.0,
                0.0,
                1.0 / y,
                0.0,
                0.0,
                0.0,
                0.0,
                1.0 / z,
                0.0,
                0.0,
                0.0,
                0.0,
                1.0,
            ),
        }
    }
    pub fn rotate_x(theta: Float) -> Transform {
        let sin_theta: Float = radians(theta).sin();
        let cos_theta: Float = radians(theta).cos();
        let m = Matrix4x4::new(
            1.0, 0.0, 0.0, 0.0, 0.0, cos_theta, -sin_theta, 0.0, 0.0, sin_theta, cos_theta, 0.0,
            0.0, 0.0, 0.0, 1.0,
        );
        Transform {
            m,
            m_inv: Matrix4x4::transpose(&m),
        }
    }
    pub fn rotate_y(theta: Float) -> Transform {
        let sin_theta: Float = radians(theta).sin();
        let cos_theta: Float = radians(theta).cos();
        let m = Matrix4x4::new(
            cos_theta, 0.0, sin_theta, 0.0, 0.0, 1.0, 0.0, 0.0, -sin_theta, 0.0, cos_theta, 0.0,
            0.0, 0.0, 0.0, 1.0,
        );
        Transform {
            m,
            m_inv: Matrix4x4::transpose(&m),
        }
    }
    pub fn rotate_z(theta: Float) -> Transform {
        let sin_theta: Float = radians(theta).sin();
        let cos_theta: Float = radians(theta).cos();
        let m = Matrix4x4::new(
            cos_theta, -sin_theta, 0.0, 0.0, sin_theta, cos_theta, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 1.0,
        );
        Transform {
            m,
            m_inv: Matrix4x4::transpose(&m),
        }
    }
    pub fn rotate(theta: Float, axis: &Vector3f) -> Transform {
        let a: Vector3f = axis.normalize();
        let sin_theta: Float = radians(theta).sin();
        let cos_theta: Float = radians(theta).cos();
        let mut m = Matrix4x4::default();
        // compute rotation of first basis vector
        m.m[0][0] = a.x * a.x + (1.0 - a.x * a.x) * cos_theta;
        m.m[0][1] = a.x * a.y * (1.0 - cos_theta) - a.z * sin_theta;
        m.m[0][2] = a.x * a.z * (1.0 - cos_theta) + a.y * sin_theta;
        m.m[0][3] = 0.0;
        // compute rotations of second basis vectors
        m.m[1][0] = a.x * a.y * (1.0 - cos_theta) + a.z * sin_theta;
        m.m[1][1] = a.y * a.y + (1.0 - a.y * a.y) * cos_theta;
        m.m[1][2] = a.y * a.z * (1.0 - cos_theta) - a.x * sin_theta;
        m.m[1][3] = 0.0;
        // compute rotations of third basis vectors
        m.m[2][0] = a.x * a.z * (1.0 - cos_theta) - a.y * sin_theta;
        m.m[2][1] = a.y * a.z * (1.0 - cos_theta) + a.x * sin_theta;
        m.m[2][2] = a.z * a.z + (1.0 - a.z * a.z) * cos_theta;
        m.m[2][3] = 0.0;
        Transform {
            m,
            m_inv: Matrix4x4::transpose(&m),
        }
    }
    pub fn look_at(pos: &Point3f, look: &Point3f, up: &Vector3f) -> Transform {
        let mut camera_to_world = Matrix4x4::default();
        // initialize fourth column of viewing matrix
        camera_to_world.m[0][3] = pos.x;
        camera_to_world.m[1][3] = pos.y;
        camera_to_world.m[2][3] = pos.z;
        camera_to_world.m[3][3] = 1.0;
        // initialize first three columns of viewing matrix
        let dir: Vector3f = (*look - *pos).normalize();
        if vec3_cross_vec3(&up.normalize(), &dir).length() == 0.0 {
            println!(
                "\"up\" vector ({}, {}, {}) and viewing direction ({}, {}, {}) passed to \
                 LookAt are pointing in the same direction.  Using the identity \
                 transformation.",
                up.x, up.y, up.z, dir.x, dir.y, dir.z
            );
            Transform::default()
        } else {
            let left: Vector3f = vec3_cross_vec3(&up.normalize(), &dir).normalize();
            let new_up: Vector3f = vec3_cross_vec3(&dir, &left);
            camera_to_world.m[0][0] = left.x;
            camera_to_world.m[1][0] = left.y;
            camera_to_world.m[2][0] = left.z;
            camera_to_world.m[3][0] = 0.0;
            camera_to_world.m[0][1] = new_up.x;
            camera_to_world.m[1][1] = new_up.y;
            camera_to_world.m[2][1] = new_up.z;
            camera_to_world.m[3][1] = 0.0;
            camera_to_world.m[0][2] = dir.x;
            camera_to_world.m[1][2] = dir.y;
            camera_to_world.m[2][2] = dir.z;
            camera_to_world.m[3][2] = 0.0;
            Transform {
                m: Matrix4x4::inverse(&camera_to_world),
                m_inv: camera_to_world,
            }
        }
    }
    pub fn orthographic(z_near: Float, z_far: Float) -> Transform {
        let translate: Transform = Transform::translate(&Vector3f {
            x: 0.0,
            y: 0.0,
            z: -z_near,
        });
        let scale: Transform = Transform::scale(1.0, 1.0, 1.0 / (z_far - z_near));
        scale * translate
    }
    pub fn perspective(fov: Float, n: Float, f: Float) -> Transform {
        // perform projective divide for perspective projection
        let persp = Matrix4x4::new(
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            f / (f - n),
            -f * n / (f - n),
            0.0,
            0.0,
            1.0,
            0.0,
        );
        // scale canonical perspective view to specified field of view
        let inv_tan_ang: Float = 1.0 / (radians(fov) / 2.0).tan();
        let scale: Transform = Transform::scale(inv_tan_ang, inv_tan_ang, 1.0);
        let persp_trans: Transform = Transform {
            m: persp,
            m_inv: Matrix4x4::inverse(&persp),
        };
        scale * persp_trans
    }
    pub fn transform_point(&self, p: &Point3f) -> Point3f {
        let x: Float = p.x;
        let y: Float = p.y;
        let z: Float = p.z;
        let xp: Float =
            self.m.m[0][0] * x + self.m.m[0][1] * y + self.m.m[0][2] * z + self.m.m[0][3];
        let yp: Float =
            self.m.m[1][0] * x + self.m.m[1][1] * y + self.m.m[1][2] * z + self.m.m[1][3];
        let zp: Float =
            self.m.m[2][0] * x + self.m.m[2][1] * y + self.m.m[2][2] * z + self.m.m[2][3];
        let wp: Float =
            self.m.m[3][0] * x + self.m.m[3][1] * y + self.m.m[3][2] * z + self.m.m[3][3];
        assert!(wp != 0.0, "wp = {:?} != 0.0", wp);
        if wp == 1.0 as Float {
            Point3f {
                x: xp,
                y: yp,
                z: zp,
            }
        } else {
            let inv: Float = 1.0 as Float / wp;
            Point3f {
                x: inv * xp,
                y: inv * yp,
                z: inv * zp,
            }
        }
    }
    pub fn transform_vector(&self, v: &Vector3f) -> Vector3f {
        let x: Float = v.x;
        let y: Float = v.y;
        let z: Float = v.z;
        Vector3f {
            x: self.m.m[0][0] * x + self.m.m[0][1] * y + self.m.m[0][2] * z,
            y: self.m.m[1][0] * x + self.m.m[1][1] * y + self.m.m[1][2] * z,
            z: self.m.m[2][0] * x + self.m.m[2][1] * y + self.m.m[2][2] * z,
        }
    }
    pub fn transform_normal(&self, n: &Normal3f) -> Normal3f {
        let x: Float = n.x;
        let y: Float = n.y;
        let z: Float = n.z;
        Normal3f {
            x: self.m_inv.m[0][0] * x + self.m_inv.m[1][0] * y + self.m_inv.m[2][0] * z,
            y: self.m_inv.m[0][1] * x + self.m_inv.m[1][1] * y + self.m_inv.m[2][1] * z,
            z: self.m_inv.m[0][2] * x + self.m_inv.m[1][2] * y + self.m_inv.m[2][2] * z,
        }
    }
    pub fn transform_ray(&self, r: &Ray) -> Ray {
        let mut o_error: Vector3f = Vector3f::default();
        let mut o: Point3f = self.transform_point_with_error(&r.o, &mut o_error);
        let d: Vector3f = self.transform_vector(&r.d);
        // offset ray origin to edge of error bounds and compute _tMax_
        let length_squared: Float = d.length_squared();
        let mut t_max: Float = r.t_max.get();
        if length_squared > 0.0 as Float {
            let dt: Float = vec3_dot_vec3f(&d.abs(), &o_error) / length_squared;
            o += d * dt;
            t_max -= dt;
        }
        if let Some(rd) = r.differential {
            let diff: RayDifferential = RayDifferential {
                rx_origin: self.transform_point(&rd.rx_origin),
                ry_origin: self.transform_point(&rd.ry_origin),
                rx_direction: self.transform_vector(&rd.rx_direction),
                ry_direction: self.transform_vector(&rd.ry_direction),
            };
            if let Some(ref medium_arc) = r.medium {
                Ray {
                    o,
                    d,
                    t_max: Cell::new(t_max),
                    time: r.time,
                    differential: Some(diff),
                    medium: Some(medium_arc.clone()),
                }
            } else {
                Ray {
                    o,
                    d,
                    t_max: Cell::new(t_max),
                    time: r.time,
                    differential: Some(diff),
                    medium: None,
                }
            }
        } else if let Some(ref medium_arc) = r.medium {
            Ray {
                o,
                d,
                t_max: Cell::new(t_max),
                time: r.time,
                differential: None,
                medium: Some(medium_arc.clone()),
            }
        } else {
            Ray {
                o,
                d,
                t_max: Cell::new(t_max),
                time: r.time,
                differential: None,
                medium: None,
            }
        }
    }
    pub fn transform_bounds(&self, b: &Bounds3f) -> Bounds3f {
        let m: Transform = *self;
        let p: Point3f = self.transform_point(&Point3f {
            x: b.p_min.x,
            y: b.p_min.y,
            z: b.p_min.z,
        });
        let mut ret: Bounds3f = Bounds3f { p_min: p, p_max: p };
        ret = bnd3_union_pnt3f(
            &ret,
            &m.transform_point(&Point3f {
                x: b.p_max.x,
                y: b.p_min.y,
                z: b.p_min.z,
            }),
        );
        ret = bnd3_union_pnt3f(
            &ret,
            &m.transform_point(&Point3f {
                x: b.p_min.x,
                y: b.p_max.y,
                z: b.p_min.z,
            }),
        );
        ret = bnd3_union_pnt3f(
            &ret,
            &m.transform_point(&Point3f {
                x: b.p_min.x,
                y: b.p_min.y,
                z: b.p_max.z,
            }),
        );
        ret = bnd3_union_pnt3f(
            &ret,
            &m.transform_point(&Point3f {
                x: b.p_min.x,
                y: b.p_max.y,
                z: b.p_max.z,
            }),
        );
        ret = bnd3_union_pnt3f(
            &ret,
            &m.transform_point(&Point3f {
                x: b.p_max.x,
                y: b.p_max.y,
                z: b.p_min.z,
            }),
        );
        ret = bnd3_union_pnt3f(
            &ret,
            &m.transform_point(&Point3f {
                x: b.p_max.x,
                y: b.p_min.y,
                z: b.p_max.z,
            }),
        );
        ret = bnd3_union_pnt3f(
            &ret,
            &m.transform_point(&Point3f {
                x: b.p_max.x,
                y: b.p_max.y,
                z: b.p_max.z,
            }),
        );
        ret
    }
    pub fn transform_point_with_error(&self, p: &Point3f, p_error: &mut Vector3f) -> Point3f {
        let x: Float = p.x;
        let y: Float = p.y;
        let z: Float = p.z;
        // compute transformed coordinates from point _pt_
        let xp: Float =
            self.m.m[0][0] * x + self.m.m[0][1] * y + self.m.m[0][2] * z + self.m.m[0][3];
        let yp: Float =
            self.m.m[1][0] * x + self.m.m[1][1] * y + self.m.m[1][2] * z + self.m.m[1][3];
        let zp: Float =
            self.m.m[2][0] * x + self.m.m[2][1] * y + self.m.m[2][2] * z + self.m.m[2][3];
        let wp: Float =
            self.m.m[3][0] * x + self.m.m[3][1] * y + self.m.m[3][2] * z + self.m.m[3][3];
        // compute absolute error for transformed point
        let x_abs_sum: Float = (self.m.m[0][0] * x).abs()
            + (self.m.m[0][1] * y).abs()
            + (self.m.m[0][2] * z).abs()
            + self.m.m[0][3].abs();
        let y_abs_sum: Float = (self.m.m[1][0] * x).abs()
            + (self.m.m[1][1] * y).abs()
            + (self.m.m[1][2] * z).abs()
            + self.m.m[1][3].abs();
        let z_abs_sum: Float = (self.m.m[2][0] * x).abs()
            + (self.m.m[2][1] * y).abs()
            + (self.m.m[2][2] * z).abs()
            + self.m.m[2][3].abs();
        *p_error = Vector3f {
            x: x_abs_sum,
            y: y_abs_sum,
            z: z_abs_sum,
        } * gamma(3i32);
        assert!(wp != 0.0, "wp = {:?} != 0.0", wp);
        if wp == 1. {
            Point3f {
                x: xp,
                y: yp,
                z: zp,
            }
        } else {
            let inv: Float = 1.0 as Float / wp;
            Point3f {
                x: inv * xp,
                y: inv * yp,
                z: inv * zp,
            }
        }
    }
    pub fn transform_point_with_abs_error(
        &self,
        pt: &Point3f,
        pt_error: &Vector3f,
        abs_error: &mut Vector3f,
    ) -> Point3f {
        let x: Float = pt.x;
        let y: Float = pt.y;
        let z: Float = pt.z;
        // compute transformed coordinates from point _pt_
        let xp: Float =
            self.m.m[0][0] * x + self.m.m[0][1] * y + self.m.m[0][2] * z + self.m.m[0][3];
        let yp: Float =
            self.m.m[1][0] * x + self.m.m[1][1] * y + self.m.m[1][2] * z + self.m.m[1][3];
        let zp: Float =
            self.m.m[2][0] * x + self.m.m[2][1] * y + self.m.m[2][2] * z + self.m.m[2][3];
        let wp: Float =
            self.m.m[3][0] * x + self.m.m[3][1] * y + self.m.m[3][2] * z + self.m.m[3][3];
        abs_error.x = (gamma(3i32) + 1.0 as Float)
            * (self.m.m[0][0].abs() * pt_error.x
                + self.m.m[0][1].abs() * pt_error.y
                + self.m.m[0][2].abs() * pt_error.z)
            + gamma(3i32)
                * ((self.m.m[0][0] * x).abs()
                    + (self.m.m[0][1] * y).abs()
                    + (self.m.m[0][2] * z).abs()
                    + self.m.m[0][3].abs());
        abs_error.y = (gamma(3i32) + 1.0 as Float)
            * (self.m.m[1][0].abs() * pt_error.x
                + self.m.m[1][1].abs() * pt_error.y
                + self.m.m[1][2].abs() * pt_error.z)
            + gamma(3i32)
                * ((self.m.m[1][0] * x).abs()
                    + (self.m.m[1][1] * y).abs()
                    + (self.m.m[1][2] * z).abs()
                    + self.m.m[1][3].abs());
        abs_error.z = (gamma(3i32) + 1.0 as Float)
            * (self.m.m[2][0].abs() * pt_error.x
                + self.m.m[2][1].abs() * pt_error.y
                + self.m.m[2][2].abs() * pt_error.z)
            + gamma(3i32)
                * ((self.m.m[2][0] * x).abs()
                    + (self.m.m[2][1] * y).abs()
                    + (self.m.m[2][2] * z).abs()
                    + self.m.m[2][3].abs());
        assert!(wp != 0.0, "wp = {:?} != 0.0", wp);
        if wp == 1. {
            Point3f {
                x: xp,
                y: yp,
                z: zp,
            }
        } else {
            let inv: Float = 1.0 as Float / wp;
            Point3f {
                x: inv * xp,
                y: inv * yp,
                z: inv * zp,
            }
        }
    }
    pub fn transform_vector_with_error(&self, v: &Vector3f, abs_error: &mut Vector3f) -> Vector3f {
        let x: Float = v.x;
        let y: Float = v.y;
        let z: Float = v.z;
        let gamma: Float = gamma(3i32);
        abs_error.x = gamma
            * ((self.m.m[0][0] * v.x).abs()
                + (self.m.m[0][1] * v.y).abs()
                + (self.m.m[0][2] * v.z).abs());
        abs_error.y = gamma
            * ((self.m.m[1][0] * v.x).abs()
                + (self.m.m[1][1] * v.y).abs()
                + (self.m.m[1][2] * v.z).abs());
        abs_error.z = gamma
            * ((self.m.m[2][0] * v.x).abs()
                + (self.m.m[2][1] * v.y).abs()
                + (self.m.m[2][2] * v.z).abs());
        Vector3f {
            x: self.m.m[0][0] * x + self.m.m[0][1] * y + self.m.m[0][2] * z,
            y: self.m.m[1][0] * x + self.m.m[1][1] * y + self.m.m[1][2] * z,
            z: self.m.m[2][0] * x + self.m.m[2][1] * y + self.m.m[2][2] * z,
        }
    }
    pub fn transform_ray_with_error(
        &self,
        r: &Ray,
        o_error: &mut Vector3f,
        d_error: &mut Vector3f,
    ) -> Ray {
        let mut o: Point3f = self.transform_point_with_error(&r.o, o_error);
        let d: Vector3f = self.transform_vector_with_error(&r.d, d_error);
        let length_squared: Float = d.length_squared();
        if length_squared > 0.0 {
            let dt: Float = vec3_dot_vec3f(&d.abs(), &*o_error) / length_squared;
            o += d * dt;
        }
        Ray {
            o,
            d,
            t_max: Cell::new(r.t_max.get()),
            time: r.time,
            differential: None,
            medium: r.medium.clone(),
        }
    }
    pub fn transform_surface_interaction(&self, si: &mut SurfaceInteraction) {
        let mut ret: SurfaceInteraction = SurfaceInteraction::default();
        {
            // transform _p_ and _pError_ in _SurfaceInteraction_
            ret.common.p = self.transform_point_with_abs_error(
                &si.common.p,
                &si.common.p_error,
                &mut ret.common.p_error,
            );
            // transform remaining members of _SurfaceInteraction_
            ret.common.n = self.transform_normal(&si.common.n).normalize();
            ret.common.wo = self.transform_vector(&si.common.wo).normalize();
            ret.common.time = si.common.time;
        }
        ret.uv = si.uv;
        ret.shape = None; // TODO? si.shape;
        ret.dpdu = self.transform_vector(&si.dpdu);
        ret.dpdv = self.transform_vector(&si.dpdv);
        ret.dndu = self.transform_normal(&si.dndu);
        ret.dndv = self.transform_normal(&si.dndv);
        ret.shading.n = self.transform_normal(&si.shading.n).normalize();
        ret.shading.dpdu = self.transform_vector(&si.shading.dpdu);
        ret.shading.dpdv = self.transform_vector(&si.shading.dpdv);
        ret.shading.dndu = self.transform_normal(&si.shading.dndu);
        ret.shading.dndv = self.transform_normal(&si.shading.dndv);
        ret.dudx = Cell::new(si.dudx.get());
        ret.dvdx = Cell::new(si.dvdx.get());
        ret.dudy = Cell::new(si.dudy.get());
        ret.dvdy = Cell::new(si.dvdy.get());
        ret.dpdx = Cell::new(si.dpdx.get());
        ret.dpdy = Cell::new(si.dpdy.get());
        // if let Some(bsdf) = &si.bsdf {
        //     if let Some(mut bsdf2) = ret.bsdf {
        //         for bxdf_idx in 0..8 {
        //             bsdf2.bxdfs[bxdf_idx] = match bsdf.bxdfs[bxdf_idx] {
        //                 _ => bsdf.bxdfs[bxdf_idx],
        //             };
        //         }
        //     }
        // }
        // ret.bssrdf = si.bssrdf.clone();
        ret.primitive = None; // TODO? si.primitive;
        ret.shading.n = nrm_faceforward_nrm(&ret.shading.n, &ret.common.n);
        // TODO: ret.faceIndex = si.faceIndex;
        *si = ret;
    }
}

impl PartialEq for Transform {
    fn eq(&self, rhs: &Transform) -> bool {
        rhs.m == self.m && rhs.m_inv == self.m_inv
    }
}

impl Mul for Transform {
    type Output = Transform;
    fn mul(self, rhs: Transform) -> Transform {
        Transform {
            m: mtx_mul(&self.m, &rhs.m),
            m_inv: mtx_mul(&rhs.m_inv, &self.m_inv),
        }
    }
}

#[derive(Debug, Default, Copy, Clone)]
pub struct DerivativeTerm {
    kc: Float,
    kx: Float,
    ky: Float,
    kz: Float,
}

impl DerivativeTerm {
    pub fn eval(&self, p: &Point3f) -> Float {
        self.kc + self.kx * p.x + self.ky * p.y + self.kz * p.z
    }
}

#[derive(Debug, Default, Copy, Clone)]
pub struct AnimatedTransform {
    start_transform: Transform,
    end_transform: Transform,
    start_time: Float,
    end_time: Float,
    actually_animated: bool,
    t: [Vector3f; 2],
    r: [Quaternion; 2],
    s: [Matrix4x4; 2],
    has_rotation: bool,
    c1: [DerivativeTerm; 3],
    c2: [DerivativeTerm; 3],
    c3: [DerivativeTerm; 3],
    c4: [DerivativeTerm; 3],
    c5: [DerivativeTerm; 3],
}

impl AnimatedTransform {
    pub fn new(
        start_transform: &Transform,
        start_time: Float,
        end_transform: &Transform,
        end_time: Float,
    ) -> Self {
        let mut at: AnimatedTransform = AnimatedTransform::default();
        at.start_transform = *start_transform;
        at.end_transform = *end_transform;
        at.start_time = start_time;
        at.end_time = end_time;
        at.actually_animated = *start_transform != *end_transform;
        AnimatedTransform::decompose(&start_transform.m, &mut at.t[0], &mut at.r[0], &mut at.s[0]);
        AnimatedTransform::decompose(&end_transform.m, &mut at.t[1], &mut at.r[1], &mut at.s[1]);
        // flip _r[1]_ if needed to select shortest path
        if quat_dot_quat(&at.r[0], &at.r[1]) < 0.0 {
            at.r[1] = -at.r[1];
        }
        at.has_rotation = quat_dot_quat(&at.r[0], &at.r[1]) < 0.9995;
        // compute terms of motion derivative function
        if at.has_rotation {
            let cos_theta: Float = quat_dot_quat(&at.r[0], &at.r[1]);
            let theta: Float = (clamp_t(cos_theta, -1.0, 1.0)).acos();
            let qperp: Quaternion = quat_normalize(&(at.r[1] - at.r[0] * cos_theta));
            let t0x: Float = at.t[0].x;
            let t0y: Float = at.t[0].y;
            let t0z: Float = at.t[0].z;
            let t1x: Float = at.t[1].x;
            let t1y: Float = at.t[1].y;
            let t1z: Float = at.t[1].z;
            let q1x: Float = at.r[0].v.x;
            let q1y: Float = at.r[0].v.y;
            let q1z: Float = at.r[0].v.z;
            let q1w: Float = at.r[0].w;
            let qperpx: Float = qperp.v.x;
            let qperpy: Float = qperp.v.y;
            let qperpz: Float = qperp.v.z;
            let qperpw: Float = qperp.w;
            let s000: Float = at.s[0].m[0][0];
            let s001: Float = at.s[0].m[0][1];
            let s002: Float = at.s[0].m[0][2];
            let s010: Float = at.s[0].m[1][0];
            let s011: Float = at.s[0].m[1][1];
            let s012: Float = at.s[0].m[1][2];
            let s020: Float = at.s[0].m[2][0];
            let s021: Float = at.s[0].m[2][1];
            let s022: Float = at.s[0].m[2][2];
            let s100: Float = at.s[1].m[0][0];
            let s101: Float = at.s[1].m[0][1];
            let s102: Float = at.s[1].m[0][2];
            let s110: Float = at.s[1].m[1][0];
            let s111: Float = at.s[1].m[1][1];
            let s112: Float = at.s[1].m[1][2];
            let s120: Float = at.s[1].m[2][0];
            let s121: Float = at.s[1].m[2][1];
            let s122: Float = at.s[1].m[2][2];
            at.c1[0] = DerivativeTerm {
                kc: -t0x + t1x,
                kx: (-1.0 + q1y * q1y + q1z * q1z + qperpy * qperpy + qperpz * qperpz) * s000
                    + q1w * q1z * s010
                    - qperpx * qperpy * s010
                    + qperpw * qperpz * s010
                    - q1w * q1y * s020
                    - qperpw * qperpy * s020
                    - qperpx * qperpz * s020
                    + s100
                    - q1y * q1y * s100
                    - q1z * q1z * s100
                    - qperpy * qperpy * s100
                    - qperpz * qperpz * s100
                    - q1w * q1z * s110
                    + qperpx * qperpy * s110
                    - qperpw * qperpz * s110
                    + q1w * q1y * s120
                    + qperpw * qperpy * s120
                    + qperpx * qperpz * s120
                    + q1x * (-(q1y * s010) - q1z * s020 + q1y * s110 + q1z * s120),
                ky: (-1.0 + q1y * q1y + q1z * q1z + qperpy * qperpy + qperpz * qperpz) * s001
                    + q1w * q1z * s011
                    - qperpx * qperpy * s011
                    + qperpw * qperpz * s011
                    - q1w * q1y * s021
                    - qperpw * qperpy * s021
                    - qperpx * qperpz * s021
                    + s101
                    - q1y * q1y * s101
                    - q1z * q1z * s101
                    - qperpy * qperpy * s101
                    - qperpz * qperpz * s101
                    - q1w * q1z * s111
                    + qperpx * qperpy * s111
                    - qperpw * qperpz * s111
                    + q1w * q1y * s121
                    + qperpw * qperpy * s121
                    + qperpx * qperpz * s121
                    + q1x * (-(q1y * s011) - q1z * s021 + q1y * s111 + q1z * s121),
                kz: (-1.0 + q1y * q1y + q1z * q1z + qperpy * qperpy + qperpz * qperpz) * s002
                    + q1w * q1z * s012
                    - qperpx * qperpy * s012
                    + qperpw * qperpz * s012
                    - q1w * q1y * s022
                    - qperpw * qperpy * s022
                    - qperpx * qperpz * s022
                    + s102
                    - q1y * q1y * s102
                    - q1z * q1z * s102
                    - qperpy * qperpy * s102
                    - qperpz * qperpz * s102
                    - q1w * q1z * s112
                    + qperpx * qperpy * s112
                    - qperpw * qperpz * s112
                    + q1w * q1y * s122
                    + qperpw * qperpy * s122
                    + qperpx * qperpz * s122
                    + q1x * (-(q1y * s012) - q1z * s022 + q1y * s112 + q1z * s122),
            };
            at.c2[0] = DerivativeTerm {
                kc: 0.0,
                kx: -(qperpy * qperpy * s000) - qperpz * qperpz * s000 + qperpx * qperpy * s010
                    - qperpw * qperpz * s010
                    + qperpw * qperpy * s020
                    + qperpx * qperpz * s020
                    + q1y * q1y * (s000 - s100)
                    + q1z * q1z * (s000 - s100)
                    + qperpy * qperpy * s100
                    + qperpz * qperpz * s100
                    - qperpx * qperpy * s110
                    + qperpw * qperpz * s110
                    - qperpw * qperpy * s120
                    - qperpx * qperpz * s120
                    + 2.0 * q1x * qperpy * s010 * theta
                    - 2.0 * q1w * qperpz * s010 * theta
                    + 2.0 * q1w * qperpy * s020 * theta
                    + 2.0 * q1x * qperpz * s020 * theta
                    + q1y
                        * (q1x * (-s010 + s110)
                            + q1w * (-s020 + s120)
                            + 2.0 * (-2.0 * qperpy * s000 + qperpx * s010 + qperpw * s020) * theta)
                    + q1z
                        * (q1w * (s010 - s110) + q1x * (-s020 + s120)
                            - 2.0 * (2.0 * qperpz * s000 + qperpw * s010 - qperpx * s020) * theta),
                ky: -(qperpy * qperpy * s001) - qperpz * qperpz * s001 + qperpx * qperpy * s011
                    - qperpw * qperpz * s011
                    + qperpw * qperpy * s021
                    + qperpx * qperpz * s021
                    + q1y * q1y * (s001 - s101)
                    + q1z * q1z * (s001 - s101)
                    + qperpy * qperpy * s101
                    + qperpz * qperpz * s101
                    - qperpx * qperpy * s111
                    + qperpw * qperpz * s111
                    - qperpw * qperpy * s121
                    - qperpx * qperpz * s121
                    + 2.0 * q1x * qperpy * s011 * theta
                    - 2.0 * q1w * qperpz * s011 * theta
                    + 2.0 * q1w * qperpy * s021 * theta
                    + 2.0 * q1x * qperpz * s021 * theta
                    + q1y
                        * (q1x * (-s011 + s111)
                            + q1w * (-s021 + s121)
                            + 2.0 * (-2.0 * qperpy * s001 + qperpx * s011 + qperpw * s021) * theta)
                    + q1z
                        * (q1w * (s011 - s111) + q1x * (-s021 + s121)
                            - 2.0 * (2.0 * qperpz * s001 + qperpw * s011 - qperpx * s021) * theta),
                kz: -(qperpy * qperpy * s002) - qperpz * qperpz * s002 + qperpx * qperpy * s012
                    - qperpw * qperpz * s012
                    + qperpw * qperpy * s022
                    + qperpx * qperpz * s022
                    + q1y * q1y * (s002 - s102)
                    + q1z * q1z * (s002 - s102)
                    + qperpy * qperpy * s102
                    + qperpz * qperpz * s102
                    - qperpx * qperpy * s112
                    + qperpw * qperpz * s112
                    - qperpw * qperpy * s122
                    - qperpx * qperpz * s122
                    + 2.0 * q1x * qperpy * s012 * theta
                    - 2.0 * q1w * qperpz * s012 * theta
                    + 2.0 * q1w * qperpy * s022 * theta
                    + 2.0 * q1x * qperpz * s022 * theta
                    + q1y
                        * (q1x * (-s012 + s112)
                            + q1w * (-s022 + s122)
                            + 2.0 * (-2.0 * qperpy * s002 + qperpx * s012 + qperpw * s022) * theta)
                    + q1z
                        * (q1w * (s012 - s112) + q1x * (-s022 + s122)
                            - 2.0 * (2.0 * qperpz * s002 + qperpw * s012 - qperpx * s022) * theta),
            };
            at.c3[0] = DerivativeTerm {
                kc: 0.0,
                kx: -2.0
                    * (q1x * qperpy * s010 - q1w * qperpz * s010
                        + q1w * qperpy * s020
                        + q1x * qperpz * s020
                        - q1x * qperpy * s110
                        + q1w * qperpz * s110
                        - q1w * qperpy * s120
                        - q1x * qperpz * s120
                        + q1y
                            * (-2.0 * qperpy * s000
                                + qperpx * s010
                                + qperpw * s020
                                + 2.0 * qperpy * s100
                                - qperpx * s110
                                - qperpw * s120)
                        + q1z
                            * (-2.0 * qperpz * s000 - qperpw * s010
                                + qperpx * s020
                                + 2.0 * qperpz * s100
                                + qperpw * s110
                                - qperpx * s120))
                    * theta,
                ky: -2.0
                    * (q1x * qperpy * s011 - q1w * qperpz * s011
                        + q1w * qperpy * s021
                        + q1x * qperpz * s021
                        - q1x * qperpy * s111
                        + q1w * qperpz * s111
                        - q1w * qperpy * s121
                        - q1x * qperpz * s121
                        + q1y
                            * (-2.0 * qperpy * s001
                                + qperpx * s011
                                + qperpw * s021
                                + 2.0 * qperpy * s101
                                - qperpx * s111
                                - qperpw * s121)
                        + q1z
                            * (-2.0 * qperpz * s001 - qperpw * s011
                                + qperpx * s021
                                + 2.0 * qperpz * s101
                                + qperpw * s111
                                - qperpx * s121))
                    * theta,
                kz: -2.0
                    * (q1x * qperpy * s012 - q1w * qperpz * s012
                        + q1w * qperpy * s022
                        + q1x * qperpz * s022
                        - q1x * qperpy * s112
                        + q1w * qperpz * s112
                        - q1w * qperpy * s122
                        - q1x * qperpz * s122
                        + q1y
                            * (-2.0 * qperpy * s002
                                + qperpx * s012
                                + qperpw * s022
                                + 2.0 * qperpy * s102
                                - qperpx * s112
                                - qperpw * s122)
                        + q1z
                            * (-2.0 * qperpz * s002 - qperpw * s012
                                + qperpx * s022
                                + 2.0 * qperpz * s102
                                + qperpw * s112
                                - qperpx * s122))
                    * theta,
            };
            at.c4[0] = DerivativeTerm {
                kc: 0.0,
                kx: -(q1x * qperpy * s010) + q1w * qperpz * s010
                    - q1w * qperpy * s020
                    - q1x * qperpz * s020
                    + q1x * qperpy * s110
                    - q1w * qperpz * s110
                    + q1w * qperpy * s120
                    + q1x * qperpz * s120
                    + 2.0 * q1y * q1y * s000 * theta
                    + 2.0 * q1z * q1z * s000 * theta
                    - 2.0 * qperpy * qperpy * s000 * theta
                    - 2.0 * qperpz * qperpz * s000 * theta
                    + 2.0 * qperpx * qperpy * s010 * theta
                    - 2.0 * qperpw * qperpz * s010 * theta
                    + 2.0 * qperpw * qperpy * s020 * theta
                    + 2.0 * qperpx * qperpz * s020 * theta
                    + q1y
                        * (-(qperpx * s010) - qperpw * s020
                            + 2.0 * qperpy * (s000 - s100)
                            + qperpx * s110
                            + qperpw * s120
                            - 2.0 * q1x * s010 * theta
                            - 2.0 * q1w * s020 * theta)
                    + q1z
                        * (2.0 * qperpz * s000 + qperpw * s010
                            - qperpx * s020
                            - 2.0 * qperpz * s100
                            - qperpw * s110
                            + qperpx * s120
                            + 2.0 * q1w * s010 * theta
                            - 2.0 * q1x * s020 * theta),
                ky: -(q1x * qperpy * s011) + q1w * qperpz * s011
                    - q1w * qperpy * s021
                    - q1x * qperpz * s021
                    + q1x * qperpy * s111
                    - q1w * qperpz * s111
                    + q1w * qperpy * s121
                    + q1x * qperpz * s121
                    + 2.0 * q1y * q1y * s001 * theta
                    + 2.0 * q1z * q1z * s001 * theta
                    - 2.0 * qperpy * qperpy * s001 * theta
                    - 2.0 * qperpz * qperpz * s001 * theta
                    + 2.0 * qperpx * qperpy * s011 * theta
                    - 2.0 * qperpw * qperpz * s011 * theta
                    + 2.0 * qperpw * qperpy * s021 * theta
                    + 2.0 * qperpx * qperpz * s021 * theta
                    + q1y
                        * (-(qperpx * s011) - qperpw * s021
                            + 2.0 * qperpy * (s001 - s101)
                            + qperpx * s111
                            + qperpw * s121
                            - 2.0 * q1x * s011 * theta
                            - 2.0 * q1w * s021 * theta)
                    + q1z
                        * (2.0 * qperpz * s001 + qperpw * s011
                            - qperpx * s021
                            - 2.0 * qperpz * s101
                            - qperpw * s111
                            + qperpx * s121
                            + 2.0 * q1w * s011 * theta
                            - 2.0 * q1x * s021 * theta),
                kz: -(q1x * qperpy * s012) + q1w * qperpz * s012
                    - q1w * qperpy * s022
                    - q1x * qperpz * s022
                    + q1x * qperpy * s112
                    - q1w * qperpz * s112
                    + q1w * qperpy * s122
                    + q1x * qperpz * s122
                    + 2.0 * q1y * q1y * s002 * theta
                    + 2.0 * q1z * q1z * s002 * theta
                    - 2.0 * qperpy * qperpy * s002 * theta
                    - 2.0 * qperpz * qperpz * s002 * theta
                    + 2.0 * qperpx * qperpy * s012 * theta
                    - 2.0 * qperpw * qperpz * s012 * theta
                    + 2.0 * qperpw * qperpy * s022 * theta
                    + 2.0 * qperpx * qperpz * s022 * theta
                    + q1y
                        * (-(qperpx * s012) - qperpw * s022
                            + 2.0 * qperpy * (s002 - s102)
                            + qperpx * s112
                            + qperpw * s122
                            - 2.0 * q1x * s012 * theta
                            - 2.0 * q1w * s022 * theta)
                    + q1z
                        * (2.0 * qperpz * s002 + qperpw * s012
                            - qperpx * s022
                            - 2.0 * qperpz * s102
                            - qperpw * s112
                            + qperpx * s122
                            + 2.0 * q1w * s012 * theta
                            - 2.0 * q1x * s022 * theta),
            };
            at.c5[0] = DerivativeTerm {
                kc: 0.0,
                kx: 2.0
                    * (qperpy * qperpy * s000 + qperpz * qperpz * s000 - qperpx * qperpy * s010
                        + qperpw * qperpz * s010
                        - qperpw * qperpy * s020
                        - qperpx * qperpz * s020
                        - qperpy * qperpy * s100
                        - qperpz * qperpz * s100
                        + q1y * q1y * (-s000 + s100)
                        + q1z * q1z * (-s000 + s100)
                        + qperpx * qperpy * s110
                        - qperpw * qperpz * s110
                        + q1y * (q1x * (s010 - s110) + q1w * (s020 - s120))
                        + qperpw * qperpy * s120
                        + qperpx * qperpz * s120
                        + q1z * (-(q1w * s010) + q1x * s020 + q1w * s110 - q1x * s120))
                    * theta,
                ky: 2.0
                    * (qperpy * qperpy * s001 + qperpz * qperpz * s001 - qperpx * qperpy * s011
                        + qperpw * qperpz * s011
                        - qperpw * qperpy * s021
                        - qperpx * qperpz * s021
                        - qperpy * qperpy * s101
                        - qperpz * qperpz * s101
                        + q1y * q1y * (-s001 + s101)
                        + q1z * q1z * (-s001 + s101)
                        + qperpx * qperpy * s111
                        - qperpw * qperpz * s111
                        + q1y * (q1x * (s011 - s111) + q1w * (s021 - s121))
                        + qperpw * qperpy * s121
                        + qperpx * qperpz * s121
                        + q1z * (-(q1w * s011) + q1x * s021 + q1w * s111 - q1x * s121))
                    * theta,
                kz: 2.0
                    * (qperpy * qperpy * s002 + qperpz * qperpz * s002 - qperpx * qperpy * s012
                        + qperpw * qperpz * s012
                        - qperpw * qperpy * s022
                        - qperpx * qperpz * s022
                        - qperpy * qperpy * s102
                        - qperpz * qperpz * s102
                        + q1y * q1y * (-s002 + s102)
                        + q1z * q1z * (-s002 + s102)
                        + qperpx * qperpy * s112
                        - qperpw * qperpz * s112
                        + q1y * (q1x * (s012 - s112) + q1w * (s022 - s122))
                        + qperpw * qperpy * s122
                        + qperpx * qperpz * s122
                        + q1z * (-(q1w * s012) + q1x * s022 + q1w * s112 - q1x * s122))
                    * theta,
            };
            at.c1[1] = DerivativeTerm {
                kc: -t0y + t1y,
                kx: -(qperpx * qperpy * s000) - qperpw * qperpz * s000 - s010
                    + q1z * q1z * s010
                    + qperpx * qperpx * s010
                    + qperpz * qperpz * s010
                    - q1y * q1z * s020
                    + qperpw * qperpx * s020
                    - qperpy * qperpz * s020
                    + qperpx * qperpy * s100
                    + qperpw * qperpz * s100
                    + q1w * q1z * (-s000 + s100)
                    + q1x * q1x * (s010 - s110)
                    + s110
                    - q1z * q1z * s110
                    - qperpx * qperpx * s110
                    - qperpz * qperpz * s110
                    + q1x * (q1y * (-s000 + s100) + q1w * (s020 - s120))
                    + q1y * q1z * s120
                    - qperpw * qperpx * s120
                    + qperpy * qperpz * s120,
                ky: -(qperpx * qperpy * s001) - qperpw * qperpz * s001 - s011
                    + q1z * q1z * s011
                    + qperpx * qperpx * s011
                    + qperpz * qperpz * s011
                    - q1y * q1z * s021
                    + qperpw * qperpx * s021
                    - qperpy * qperpz * s021
                    + qperpx * qperpy * s101
                    + qperpw * qperpz * s101
                    + q1w * q1z * (-s001 + s101)
                    + q1x * q1x * (s011 - s111)
                    + s111
                    - q1z * q1z * s111
                    - qperpx * qperpx * s111
                    - qperpz * qperpz * s111
                    + q1x * (q1y * (-s001 + s101) + q1w * (s021 - s121))
                    + q1y * q1z * s121
                    - qperpw * qperpx * s121
                    + qperpy * qperpz * s121,
                kz: -(qperpx * qperpy * s002) - qperpw * qperpz * s002 - s012
                    + q1z * q1z * s012
                    + qperpx * qperpx * s012
                    + qperpz * qperpz * s012
                    - q1y * q1z * s022
                    + qperpw * qperpx * s022
                    - qperpy * qperpz * s022
                    + qperpx * qperpy * s102
                    + qperpw * qperpz * s102
                    + q1w * q1z * (-s002 + s102)
                    + q1x * q1x * (s012 - s112)
                    + s112
                    - q1z * q1z * s112
                    - qperpx * qperpx * s112
                    - qperpz * qperpz * s112
                    + q1x * (q1y * (-s002 + s102) + q1w * (s022 - s122))
                    + q1y * q1z * s122
                    - qperpw * qperpx * s122
                    + qperpy * qperpz * s122,
            };
            at.c2[1] = DerivativeTerm {
                kc: 0.0,
                kx: qperpx * qperpy * s000 + qperpw * qperpz * s000 + q1z * q1z * s010
                    - qperpx * qperpx * s010
                    - qperpz * qperpz * s010
                    - q1y * q1z * s020
                    - qperpw * qperpx * s020
                    + qperpy * qperpz * s020
                    - qperpx * qperpy * s100
                    - qperpw * qperpz * s100
                    + q1x * q1x * (s010 - s110)
                    - q1z * q1z * s110
                    + qperpx * qperpx * s110
                    + qperpz * qperpz * s110
                    + q1y * q1z * s120
                    + qperpw * qperpx * s120
                    - qperpy * qperpz * s120
                    + 2.0 * q1z * qperpw * s000 * theta
                    + 2.0 * q1y * qperpx * s000 * theta
                    - 4.0 * q1z * qperpz * s010 * theta
                    + 2.0 * q1z * qperpy * s020 * theta
                    + 2.0 * q1y * qperpz * s020 * theta
                    + q1x
                        * (q1w * s020 + q1y * (-s000 + s100) - q1w * s120
                            + 2.0 * qperpy * s000 * theta
                            - 4.0 * qperpx * s010 * theta
                            - 2.0 * qperpw * s020 * theta)
                    + q1w
                        * (-(q1z * s000) + q1z * s100 + 2.0 * qperpz * s000 * theta
                            - 2.0 * qperpx * s020 * theta),
                ky: qperpx * qperpy * s001 + qperpw * qperpz * s001 + q1z * q1z * s011
                    - qperpx * qperpx * s011
                    - qperpz * qperpz * s011
                    - q1y * q1z * s021
                    - qperpw * qperpx * s021
                    + qperpy * qperpz * s021
                    - qperpx * qperpy * s101
                    - qperpw * qperpz * s101
                    + q1x * q1x * (s011 - s111)
                    - q1z * q1z * s111
                    + qperpx * qperpx * s111
                    + qperpz * qperpz * s111
                    + q1y * q1z * s121
                    + qperpw * qperpx * s121
                    - qperpy * qperpz * s121
                    + 2.0 * q1z * qperpw * s001 * theta
                    + 2.0 * q1y * qperpx * s001 * theta
                    - 4.0 * q1z * qperpz * s011 * theta
                    + 2.0 * q1z * qperpy * s021 * theta
                    + 2.0 * q1y * qperpz * s021 * theta
                    + q1x
                        * (q1w * s021 + q1y * (-s001 + s101) - q1w * s121
                            + 2.0 * qperpy * s001 * theta
                            - 4.0 * qperpx * s011 * theta
                            - 2.0 * qperpw * s021 * theta)
                    + q1w
                        * (-(q1z * s001) + q1z * s101 + 2.0 * qperpz * s001 * theta
                            - 2.0 * qperpx * s021 * theta),
                kz: qperpx * qperpy * s002 + qperpw * qperpz * s002 + q1z * q1z * s012
                    - qperpx * qperpx * s012
                    - qperpz * qperpz * s012
                    - q1y * q1z * s022
                    - qperpw * qperpx * s022
                    + qperpy * qperpz * s022
                    - qperpx * qperpy * s102
                    - qperpw * qperpz * s102
                    + q1x * q1x * (s012 - s112)
                    - q1z * q1z * s112
                    + qperpx * qperpx * s112
                    + qperpz * qperpz * s112
                    + q1y * q1z * s122
                    + qperpw * qperpx * s122
                    - qperpy * qperpz * s122
                    + 2.0 * q1z * qperpw * s002 * theta
                    + 2.0 * q1y * qperpx * s002 * theta
                    - 4.0 * q1z * qperpz * s012 * theta
                    + 2.0 * q1z * qperpy * s022 * theta
                    + 2.0 * q1y * qperpz * s022 * theta
                    + q1x
                        * (q1w * s022 + q1y * (-s002 + s102) - q1w * s122
                            + 2.0 * qperpy * s002 * theta
                            - 4.0 * qperpx * s012 * theta
                            - 2.0 * qperpw * s022 * theta)
                    + q1w
                        * (-(q1z * s002) + q1z * s102 + 2.0 * qperpz * s002 * theta
                            - 2.0 * qperpx * s022 * theta),
            };
            at.c3[1] = DerivativeTerm {
                kc: 0.0,
                kx: 2.0
                    * (-(q1x * qperpy * s000) - q1w * qperpz * s000
                        + 2.0 * q1x * qperpx * s010
                        + q1x * qperpw * s020
                        + q1w * qperpx * s020
                        + q1x * qperpy * s100
                        + q1w * qperpz * s100
                        - 2.0 * q1x * qperpx * s110
                        - q1x * qperpw * s120
                        - q1w * qperpx * s120
                        + q1z
                            * (2.0 * qperpz * s010 - qperpy * s020 + qperpw * (-s000 + s100)
                                - 2.0 * qperpz * s110
                                + qperpy * s120)
                        + q1y * (-(qperpx * s000) - qperpz * s020 + qperpx * s100 + qperpz * s120))
                    * theta,
                ky: 2.0
                    * (-(q1x * qperpy * s001) - q1w * qperpz * s001
                        + 2.0 * q1x * qperpx * s011
                        + q1x * qperpw * s021
                        + q1w * qperpx * s021
                        + q1x * qperpy * s101
                        + q1w * qperpz * s101
                        - 2.0 * q1x * qperpx * s111
                        - q1x * qperpw * s121
                        - q1w * qperpx * s121
                        + q1z
                            * (2.0 * qperpz * s011 - qperpy * s021 + qperpw * (-s001 + s101)
                                - 2.0 * qperpz * s111
                                + qperpy * s121)
                        + q1y * (-(qperpx * s001) - qperpz * s021 + qperpx * s101 + qperpz * s121))
                    * theta,
                kz: 2.0
                    * (-(q1x * qperpy * s002) - q1w * qperpz * s002
                        + 2.0 * q1x * qperpx * s012
                        + q1x * qperpw * s022
                        + q1w * qperpx * s022
                        + q1x * qperpy * s102
                        + q1w * qperpz * s102
                        - 2.0 * q1x * qperpx * s112
                        - q1x * qperpw * s122
                        - q1w * qperpx * s122
                        + q1z
                            * (2.0 * qperpz * s012 - qperpy * s022 + qperpw * (-s002 + s102)
                                - 2.0 * qperpz * s112
                                + qperpy * s122)
                        + q1y * (-(qperpx * s002) - qperpz * s022 + qperpx * s102 + qperpz * s122))
                    * theta,
            };
            at.c4[1] = DerivativeTerm {
                kc: 0.0,
                kx: -(q1x * qperpy * s000) - q1w * qperpz * s000
                    + 2.0 * q1x * qperpx * s010
                    + q1x * qperpw * s020
                    + q1w * qperpx * s020
                    + q1x * qperpy * s100
                    + q1w * qperpz * s100
                    - 2.0 * q1x * qperpx * s110
                    - q1x * qperpw * s120
                    - q1w * qperpx * s120
                    + 2.0 * qperpx * qperpy * s000 * theta
                    + 2.0 * qperpw * qperpz * s000 * theta
                    + 2.0 * q1x * q1x * s010 * theta
                    + 2.0 * q1z * q1z * s010 * theta
                    - 2.0 * qperpx * qperpx * s010 * theta
                    - 2.0 * qperpz * qperpz * s010 * theta
                    + 2.0 * q1w * q1x * s020 * theta
                    - 2.0 * qperpw * qperpx * s020 * theta
                    + 2.0 * qperpy * qperpz * s020 * theta
                    + q1y
                        * (-(qperpx * s000) - qperpz * s020 + qperpx * s100 + qperpz * s120
                            - 2.0 * q1x * s000 * theta)
                    + q1z
                        * (2.0 * qperpz * s010 - qperpy * s020 + qperpw * (-s000 + s100)
                            - 2.0 * qperpz * s110
                            + qperpy * s120
                            - 2.0 * q1w * s000 * theta
                            - 2.0 * q1y * s020 * theta),
                ky: -(q1x * qperpy * s001) - q1w * qperpz * s001
                    + 2.0 * q1x * qperpx * s011
                    + q1x * qperpw * s021
                    + q1w * qperpx * s021
                    + q1x * qperpy * s101
                    + q1w * qperpz * s101
                    - 2.0 * q1x * qperpx * s111
                    - q1x * qperpw * s121
                    - q1w * qperpx * s121
                    + 2.0 * qperpx * qperpy * s001 * theta
                    + 2.0 * qperpw * qperpz * s001 * theta
                    + 2.0 * q1x * q1x * s011 * theta
                    + 2.0 * q1z * q1z * s011 * theta
                    - 2.0 * qperpx * qperpx * s011 * theta
                    - 2.0 * qperpz * qperpz * s011 * theta
                    + 2.0 * q1w * q1x * s021 * theta
                    - 2.0 * qperpw * qperpx * s021 * theta
                    + 2.0 * qperpy * qperpz * s021 * theta
                    + q1y
                        * (-(qperpx * s001) - qperpz * s021 + qperpx * s101 + qperpz * s121
                            - 2.0 * q1x * s001 * theta)
                    + q1z
                        * (2.0 * qperpz * s011 - qperpy * s021 + qperpw * (-s001 + s101)
                            - 2.0 * qperpz * s111
                            + qperpy * s121
                            - 2.0 * q1w * s001 * theta
                            - 2.0 * q1y * s021 * theta),
                kz: -(q1x * qperpy * s002) - q1w * qperpz * s002
                    + 2.0 * q1x * qperpx * s012
                    + q1x * qperpw * s022
                    + q1w * qperpx * s022
                    + q1x * qperpy * s102
                    + q1w * qperpz * s102
                    - 2.0 * q1x * qperpx * s112
                    - q1x * qperpw * s122
                    - q1w * qperpx * s122
                    + 2.0 * qperpx * qperpy * s002 * theta
                    + 2.0 * qperpw * qperpz * s002 * theta
                    + 2.0 * q1x * q1x * s012 * theta
                    + 2.0 * q1z * q1z * s012 * theta
                    - 2.0 * qperpx * qperpx * s012 * theta
                    - 2.0 * qperpz * qperpz * s012 * theta
                    + 2.0 * q1w * q1x * s022 * theta
                    - 2.0 * qperpw * qperpx * s022 * theta
                    + 2.0 * qperpy * qperpz * s022 * theta
                    + q1y
                        * (-(qperpx * s002) - qperpz * s022 + qperpx * s102 + qperpz * s122
                            - 2.0 * q1x * s002 * theta)
                    + q1z
                        * (2.0 * qperpz * s012 - qperpy * s022 + qperpw * (-s002 + s102)
                            - 2.0 * qperpz * s112
                            + qperpy * s122
                            - 2.0 * q1w * s002 * theta
                            - 2.0 * q1y * s022 * theta),
            };
            at.c5[1] = DerivativeTerm {
                kc: 0.,
                kx: -2.0
                    * (qperpx * qperpy * s000 + qperpw * qperpz * s000 + q1z * q1z * s010
                        - qperpx * qperpx * s010
                        - qperpz * qperpz * s010
                        - q1y * q1z * s020
                        - qperpw * qperpx * s020
                        + qperpy * qperpz * s020
                        - qperpx * qperpy * s100
                        - qperpw * qperpz * s100
                        + q1w * q1z * (-s000 + s100)
                        + q1x * q1x * (s010 - s110)
                        - q1z * q1z * s110
                        + qperpx * qperpx * s110
                        + qperpz * qperpz * s110
                        + q1x * (q1y * (-s000 + s100) + q1w * (s020 - s120))
                        + q1y * q1z * s120
                        + qperpw * qperpx * s120
                        - qperpy * qperpz * s120)
                    * theta,
                ky: -2.0
                    * (qperpx * qperpy * s001 + qperpw * qperpz * s001 + q1z * q1z * s011
                        - qperpx * qperpx * s011
                        - qperpz * qperpz * s011
                        - q1y * q1z * s021
                        - qperpw * qperpx * s021
                        + qperpy * qperpz * s021
                        - qperpx * qperpy * s101
                        - qperpw * qperpz * s101
                        + q1w * q1z * (-s001 + s101)
                        + q1x * q1x * (s011 - s111)
                        - q1z * q1z * s111
                        + qperpx * qperpx * s111
                        + qperpz * qperpz * s111
                        + q1x * (q1y * (-s001 + s101) + q1w * (s021 - s121))
                        + q1y * q1z * s121
                        + qperpw * qperpx * s121
                        - qperpy * qperpz * s121)
                    * theta,
                kz: -2.0
                    * (qperpx * qperpy * s002 + qperpw * qperpz * s002 + q1z * q1z * s012
                        - qperpx * qperpx * s012
                        - qperpz * qperpz * s012
                        - q1y * q1z * s022
                        - qperpw * qperpx * s022
                        + qperpy * qperpz * s022
                        - qperpx * qperpy * s102
                        - qperpw * qperpz * s102
                        + q1w * q1z * (-s002 + s102)
                        + q1x * q1x * (s012 - s112)
                        - q1z * q1z * s112
                        + qperpx * qperpx * s112
                        + qperpz * qperpz * s112
                        + q1x * (q1y * (-s002 + s102) + q1w * (s022 - s122))
                        + q1y * q1z * s122
                        + qperpw * qperpx * s122
                        - qperpy * qperpz * s122)
                    * theta,
            };
            at.c1[2] = DerivativeTerm {
                kc: -t0z + t1z,
                kx: (qperpw * qperpy * s000
                    - qperpx * qperpz * s000
                    - q1y * q1z * s010
                    - qperpw * qperpx * s010
                    - qperpy * qperpz * s010
                    - s020
                    + q1y * q1y * s020
                    + qperpx * qperpx * s020
                    + qperpy * qperpy * s020
                    - qperpw * qperpy * s100
                    + qperpx * qperpz * s100
                    + q1x * q1z * (-s000 + s100)
                    + q1y * q1z * s110
                    + qperpw * qperpx * s110
                    + qperpy * qperpz * s110
                    + q1w * (q1y * (s000 - s100) + q1x * (-s010 + s110))
                    + q1x * q1x * (s020 - s120)
                    + s120
                    - q1y * q1y * s120
                    - qperpx * qperpx * s120
                    - qperpy * qperpy * s120),
                ky: (qperpw * qperpy * s001
                    - qperpx * qperpz * s001
                    - q1y * q1z * s011
                    - qperpw * qperpx * s011
                    - qperpy * qperpz * s011
                    - s021
                    + q1y * q1y * s021
                    + qperpx * qperpx * s021
                    + qperpy * qperpy * s021
                    - qperpw * qperpy * s101
                    + qperpx * qperpz * s101
                    + q1x * q1z * (-s001 + s101)
                    + q1y * q1z * s111
                    + qperpw * qperpx * s111
                    + qperpy * qperpz * s111
                    + q1w * (q1y * (s001 - s101) + q1x * (-s011 + s111))
                    + q1x * q1x * (s021 - s121)
                    + s121
                    - q1y * q1y * s121
                    - qperpx * qperpx * s121
                    - qperpy * qperpy * s121),
                kz: (qperpw * qperpy * s002
                    - qperpx * qperpz * s002
                    - q1y * q1z * s012
                    - qperpw * qperpx * s012
                    - qperpy * qperpz * s012
                    - s022
                    + q1y * q1y * s022
                    + qperpx * qperpx * s022
                    + qperpy * qperpy * s022
                    - qperpw * qperpy * s102
                    + qperpx * qperpz * s102
                    + q1x * q1z * (-s002 + s102)
                    + q1y * q1z * s112
                    + qperpw * qperpx * s112
                    + qperpy * qperpz * s112
                    + q1w * (q1y * (s002 - s102) + q1x * (-s012 + s112))
                    + q1x * q1x * (s022 - s122)
                    + s122
                    - q1y * q1y * s122
                    - qperpx * qperpx * s122
                    - qperpy * qperpy * s122),
            };
            at.c2[2] = DerivativeTerm {
                kc: 0.0,
                kx: (q1w * q1y * s000 - q1x * q1z * s000 - qperpw * qperpy * s000
                    + qperpx * qperpz * s000
                    - q1w * q1x * s010
                    - q1y * q1z * s010
                    + qperpw * qperpx * s010
                    + qperpy * qperpz * s010
                    + q1x * q1x * s020
                    + q1y * q1y * s020
                    - qperpx * qperpx * s020
                    - qperpy * qperpy * s020
                    - q1w * q1y * s100
                    + q1x * q1z * s100
                    + qperpw * qperpy * s100
                    - qperpx * qperpz * s100
                    + q1w * q1x * s110
                    + q1y * q1z * s110
                    - qperpw * qperpx * s110
                    - qperpy * qperpz * s110
                    - q1x * q1x * s120
                    - q1y * q1y * s120
                    + qperpx * qperpx * s120
                    + qperpy * qperpy * s120
                    - 2.0 * q1y * qperpw * s000 * theta
                    + 2.0 * q1z * qperpx * s000 * theta
                    - 2.0 * q1w * qperpy * s000 * theta
                    + 2.0 * q1x * qperpz * s000 * theta
                    + 2.0 * q1x * qperpw * s010 * theta
                    + 2.0 * q1w * qperpx * s010 * theta
                    + 2.0 * q1z * qperpy * s010 * theta
                    + 2.0 * q1y * qperpz * s010 * theta
                    - 4.0 * q1x * qperpx * s020 * theta
                    - 4.0 * q1y * qperpy * s020 * theta),
                ky: (q1w * q1y * s001 - q1x * q1z * s001 - qperpw * qperpy * s001
                    + qperpx * qperpz * s001
                    - q1w * q1x * s011
                    - q1y * q1z * s011
                    + qperpw * qperpx * s011
                    + qperpy * qperpz * s011
                    + q1x * q1x * s021
                    + q1y * q1y * s021
                    - qperpx * qperpx * s021
                    - qperpy * qperpy * s021
                    - q1w * q1y * s101
                    + q1x * q1z * s101
                    + qperpw * qperpy * s101
                    - qperpx * qperpz * s101
                    + q1w * q1x * s111
                    + q1y * q1z * s111
                    - qperpw * qperpx * s111
                    - qperpy * qperpz * s111
                    - q1x * q1x * s121
                    - q1y * q1y * s121
                    + qperpx * qperpx * s121
                    + qperpy * qperpy * s121
                    - 2.0 * q1y * qperpw * s001 * theta
                    + 2.0 * q1z * qperpx * s001 * theta
                    - 2.0 * q1w * qperpy * s001 * theta
                    + 2.0 * q1x * qperpz * s001 * theta
                    + 2.0 * q1x * qperpw * s011 * theta
                    + 2.0 * q1w * qperpx * s011 * theta
                    + 2.0 * q1z * qperpy * s011 * theta
                    + 2.0 * q1y * qperpz * s011 * theta
                    - 4.0 * q1x * qperpx * s021 * theta
                    - 4.0 * q1y * qperpy * s021 * theta),
                kz: (q1w * q1y * s002 - q1x * q1z * s002 - qperpw * qperpy * s002
                    + qperpx * qperpz * s002
                    - q1w * q1x * s012
                    - q1y * q1z * s012
                    + qperpw * qperpx * s012
                    + qperpy * qperpz * s012
                    + q1x * q1x * s022
                    + q1y * q1y * s022
                    - qperpx * qperpx * s022
                    - qperpy * qperpy * s022
                    - q1w * q1y * s102
                    + q1x * q1z * s102
                    + qperpw * qperpy * s102
                    - qperpx * qperpz * s102
                    + q1w * q1x * s112
                    + q1y * q1z * s112
                    - qperpw * qperpx * s112
                    - qperpy * qperpz * s112
                    - q1x * q1x * s122
                    - q1y * q1y * s122
                    + qperpx * qperpx * s122
                    + qperpy * qperpy * s122
                    - 2.0 * q1y * qperpw * s002 * theta
                    + 2.0 * q1z * qperpx * s002 * theta
                    - 2.0 * q1w * qperpy * s002 * theta
                    + 2.0 * q1x * qperpz * s002 * theta
                    + 2.0 * q1x * qperpw * s012 * theta
                    + 2.0 * q1w * qperpx * s012 * theta
                    + 2.0 * q1z * qperpy * s012 * theta
                    + 2.0 * q1y * qperpz * s012 * theta
                    - 4.0 * q1x * qperpx * s022 * theta
                    - 4.0 * q1y * qperpy * s022 * theta),
            };
            at.c3[2] = DerivativeTerm {
                kc: 0.0,
                kx: -2.0
                    * (-(q1w * qperpy * s000)
                        + q1x * qperpz * s000
                        + q1x * qperpw * s010
                        + q1w * qperpx * s010
                        - 2.0 * q1x * qperpx * s020
                        + q1w * qperpy * s100
                        - q1x * qperpz * s100
                        - q1x * qperpw * s110
                        - q1w * qperpx * s110
                        + q1z * (qperpx * s000 + qperpy * s010 - qperpx * s100 - qperpy * s110)
                        + 2.0 * q1x * qperpx * s120
                        + q1y
                            * (qperpz * s010 - 2.0 * qperpy * s020 + qperpw * (-s000 + s100)
                                - qperpz * s110
                                + 2.0 * qperpy * s120))
                    * theta,
                ky: -2.0
                    * (-(q1w * qperpy * s001)
                        + q1x * qperpz * s001
                        + q1x * qperpw * s011
                        + q1w * qperpx * s011
                        - 2.0 * q1x * qperpx * s021
                        + q1w * qperpy * s101
                        - q1x * qperpz * s101
                        - q1x * qperpw * s111
                        - q1w * qperpx * s111
                        + q1z * (qperpx * s001 + qperpy * s011 - qperpx * s101 - qperpy * s111)
                        + 2.0 * q1x * qperpx * s121
                        + q1y
                            * (qperpz * s011 - 2.0 * qperpy * s021 + qperpw * (-s001 + s101)
                                - qperpz * s111
                                + 2.0 * qperpy * s121))
                    * theta,
                kz: -2.0
                    * (-(q1w * qperpy * s002)
                        + q1x * qperpz * s002
                        + q1x * qperpw * s012
                        + q1w * qperpx * s012
                        - 2.0 * q1x * qperpx * s022
                        + q1w * qperpy * s102
                        - q1x * qperpz * s102
                        - q1x * qperpw * s112
                        - q1w * qperpx * s112
                        + q1z * (qperpx * s002 + qperpy * s012 - qperpx * s102 - qperpy * s112)
                        + 2.0 * q1x * qperpx * s122
                        + q1y
                            * (qperpz * s012 - 2.0 * qperpy * s022 + qperpw * (-s002 + s102)
                                - qperpz * s112
                                + 2.0 * qperpy * s122))
                    * theta,
            };
            at.c4[2] = DerivativeTerm {
                kc: 0.0,
                kx: q1w * qperpy * s000
                    - q1x * qperpz * s000
                    - q1x * qperpw * s010
                    - q1w * qperpx * s010
                    + 2.0 * q1x * qperpx * s020
                    - q1w * qperpy * s100
                    + q1x * qperpz * s100
                    + q1x * qperpw * s110
                    + q1w * qperpx * s110
                    - 2.0 * q1x * qperpx * s120
                    - 2.0 * qperpw * qperpy * s000 * theta
                    + 2.0 * qperpx * qperpz * s000 * theta
                    - 2.0 * q1w * q1x * s010 * theta
                    + 2.0 * qperpw * qperpx * s010 * theta
                    + 2.0 * qperpy * qperpz * s010 * theta
                    + 2.0 * q1x * q1x * s020 * theta
                    + 2.0 * q1y * q1y * s020 * theta
                    - 2.0 * qperpx * qperpx * s020 * theta
                    - 2.0 * qperpy * qperpy * s020 * theta
                    + q1z
                        * (-(qperpx * s000) - qperpy * s010 + qperpx * s100 + qperpy * s110
                            - 2.0 * q1x * s000 * theta)
                    + q1y
                        * (-(qperpz * s010)
                            + 2.0 * qperpy * s020
                            + qperpw * (s000 - s100)
                            + qperpz * s110
                            - 2.0 * qperpy * s120
                            + 2.0 * q1w * s000 * theta
                            - 2.0 * q1z * s010 * theta),
                ky: q1w * qperpy * s001
                    - q1x * qperpz * s001
                    - q1x * qperpw * s011
                    - q1w * qperpx * s011
                    + 2.0 * q1x * qperpx * s021
                    - q1w * qperpy * s101
                    + q1x * qperpz * s101
                    + q1x * qperpw * s111
                    + q1w * qperpx * s111
                    - 2.0 * q1x * qperpx * s121
                    - 2.0 * qperpw * qperpy * s001 * theta
                    + 2.0 * qperpx * qperpz * s001 * theta
                    - 2.0 * q1w * q1x * s011 * theta
                    + 2.0 * qperpw * qperpx * s011 * theta
                    + 2.0 * qperpy * qperpz * s011 * theta
                    + 2.0 * q1x * q1x * s021 * theta
                    + 2.0 * q1y * q1y * s021 * theta
                    - 2.0 * qperpx * qperpx * s021 * theta
                    - 2.0 * qperpy * qperpy * s021 * theta
                    + q1z
                        * (-(qperpx * s001) - qperpy * s011 + qperpx * s101 + qperpy * s111
                            - 2.0 * q1x * s001 * theta)
                    + q1y
                        * (-(qperpz * s011)
                            + 2.0 * qperpy * s021
                            + qperpw * (s001 - s101)
                            + qperpz * s111
                            - 2.0 * qperpy * s121
                            + 2.0 * q1w * s001 * theta
                            - 2.0 * q1z * s011 * theta),
                kz: q1w * qperpy * s002
                    - q1x * qperpz * s002
                    - q1x * qperpw * s012
                    - q1w * qperpx * s012
                    + 2.0 * q1x * qperpx * s022
                    - q1w * qperpy * s102
                    + q1x * qperpz * s102
                    + q1x * qperpw * s112
                    + q1w * qperpx * s112
                    - 2.0 * q1x * qperpx * s122
                    - 2.0 * qperpw * qperpy * s002 * theta
                    + 2.0 * qperpx * qperpz * s002 * theta
                    - 2.0 * q1w * q1x * s012 * theta
                    + 2.0 * qperpw * qperpx * s012 * theta
                    + 2.0 * qperpy * qperpz * s012 * theta
                    + 2.0 * q1x * q1x * s022 * theta
                    + 2.0 * q1y * q1y * s022 * theta
                    - 2.0 * qperpx * qperpx * s022 * theta
                    - 2.0 * qperpy * qperpy * s022 * theta
                    + q1z
                        * (-(qperpx * s002) - qperpy * s012 + qperpx * s102 + qperpy * s112
                            - 2.0 * q1x * s002 * theta)
                    + q1y
                        * (-(qperpz * s012)
                            + 2.0 * qperpy * s022
                            + qperpw * (s002 - s102)
                            + qperpz * s112
                            - 2.0 * qperpy * s122
                            + 2.0 * q1w * s002 * theta
                            - 2.0 * q1z * s012 * theta),
            };
            at.c5[2] = DerivativeTerm {
                kc: 0.,
                kx: 2.0
                    * (qperpw * qperpy * s000 - qperpx * qperpz * s000 + q1y * q1z * s010
                        - qperpw * qperpx * s010
                        - qperpy * qperpz * s010
                        - q1y * q1y * s020
                        + qperpx * qperpx * s020
                        + qperpy * qperpy * s020
                        + q1x * q1z * (s000 - s100)
                        - qperpw * qperpy * s100
                        + qperpx * qperpz * s100
                        + q1w * (q1y * (-s000 + s100) + q1x * (s010 - s110))
                        - q1y * q1z * s110
                        + qperpw * qperpx * s110
                        + qperpy * qperpz * s110
                        + q1y * q1y * s120
                        - qperpx * qperpx * s120
                        - qperpy * qperpy * s120
                        + q1x * q1x * (-s020 + s120))
                    * theta,
                ky: 2.0
                    * (qperpw * qperpy * s001 - qperpx * qperpz * s001 + q1y * q1z * s011
                        - qperpw * qperpx * s011
                        - qperpy * qperpz * s011
                        - q1y * q1y * s021
                        + qperpx * qperpx * s021
                        + qperpy * qperpy * s021
                        + q1x * q1z * (s001 - s101)
                        - qperpw * qperpy * s101
                        + qperpx * qperpz * s101
                        + q1w * (q1y * (-s001 + s101) + q1x * (s011 - s111))
                        - q1y * q1z * s111
                        + qperpw * qperpx * s111
                        + qperpy * qperpz * s111
                        + q1y * q1y * s121
                        - qperpx * qperpx * s121
                        - qperpy * qperpy * s121
                        + q1x * q1x * (-s021 + s121))
                    * theta,
                kz: 2.0
                    * (qperpw * qperpy * s002 - qperpx * qperpz * s002 + q1y * q1z * s012
                        - qperpw * qperpx * s012
                        - qperpy * qperpz * s012
                        - q1y * q1y * s022
                        + qperpx * qperpx * s022
                        + qperpy * qperpy * s022
                        + q1x * q1z * (s002 - s102)
                        - qperpw * qperpy * s102
                        + qperpx * qperpz * s102
                        + q1w * (q1y * (-s002 + s102) + q1x * (s012 - s112))
                        - q1y * q1z * s112
                        + qperpw * qperpx * s112
                        + qperpy * qperpz * s112
                        + q1y * q1y * s122
                        - qperpx * qperpx * s122
                        - qperpy * qperpy * s122
                        + q1x * q1x * (-s022 + s122))
                    * theta,
            };
        }
        at
    }
    pub fn decompose(m: &Matrix4x4, t: &mut Vector3f, rquat: &mut Quaternion, s: &mut Matrix4x4) {
        // extract translation from transformation matrix
        t.x = m.m[0][3];
        t.y = m.m[1][3];
        t.z = m.m[2][3];
        // compute new transformation matrix _m_ without translation
        let mut matrix: Matrix4x4 = *m;
        for i in 0..3 {
            matrix.m[i][3] = 0.0;
            matrix.m[3][i] = 0.0;
        }
        matrix.m[3][3] = 1.0;
        // extract rotation _r_ from transformation matrix
        let mut norm: Float;
        let mut count: u8 = 0;
        let mut r: Matrix4x4 = matrix;
        loop {
            // compute next matrix _rnext_ in series
            let mut rnext: Matrix4x4 = Matrix4x4::default();
            let rit: Matrix4x4 = Matrix4x4::inverse(&Matrix4x4::transpose(&r));
            for i in 0..4 {
                for j in 0..4 {
                    rnext.m[i][j] = 0.5 * (r.m[i][j] + rit.m[i][j]);
                }
            }
            // compute norm of difference between _r_ and _rnext_
            norm = 0.0;
            for i in 0..3 {
                let n: Float = (r.m[i][0] - rnext.m[i][0]).abs()
                    + (r.m[i][1] - rnext.m[i][1]).abs()
                    + (r.m[i][2] - rnext.m[i][2]).abs();
                norm = norm.max(n);
            }
            r = rnext;
            count += 1;
            if count >= 100 || norm <= 0.0001 {
                break;
            }
        }
        // XXX TODO FIXME deal with flip...
        let transform: Transform = Transform {
            m: r,
            m_inv: Matrix4x4::inverse(&r.clone()),
        };
        *rquat = Quaternion::new(transform);

        // compute scale _S_ using rotation and original matrix
        *s = mtx_mul(&Matrix4x4::inverse(&r), &*m);
    }
    pub fn interpolate(&self, time: Float, t: &mut Transform) {
        // handle boundary conditions for matrix interpolation
        if !self.actually_animated || time <= self.start_time {
            *t = self.start_transform;
            return;
        }
        if time >= self.end_time {
            *t = self.end_transform;
            return;
        }
        let dt: Float = (time - self.start_time) / (self.end_time - self.start_time);
        // interpolate translation at _dt_
        let trans: Vector3f = self.t[0] * (1.0 as Float - dt) + self.t[1] * dt;

        // interpolate rotation at _dt_
        let rotate: Quaternion = quat_slerp(dt, &self.r[0], &self.r[1]);

        // interpolate scale at _dt_
        let mut scale: Matrix4x4 = Matrix4x4::default();
        for i in 0..3 {
            for j in 0..3 {
                scale.m[i][j] = lerp(dt, self.s[0].m[i][j], self.s[1].m[i][j]);
            }
        }

        // compute interpolated matrix as product of interpolated components
        *t = Transform::translate(&trans)
            * rotate.to_transform()
            * Transform {
                m: scale,
                m_inv: Matrix4x4::inverse(&scale),
            };
    }
    pub fn transform_ray(&self, r: &Ray) -> Ray {
        if !self.actually_animated || r.time <= self.start_time {
            self.start_transform.transform_ray(r)
        } else if r.time >= self.end_time {
            self.end_transform.transform_ray(r)
        } else {
            let mut t: Transform = Transform::default();
            self.interpolate(r.time, &mut t);
            t.transform_ray(r)
        }
    }
    pub fn transform_point(&self, time: Float, p: &Point3f) -> Point3f {
        if !self.actually_animated || time <= self.start_time {
            self.start_transform.transform_point(p)
        } else if time >= self.end_time {
            self.end_transform.transform_point(p)
        } else {
            let mut t: Transform = Transform::default();
            self.interpolate(time, &mut t);
            t.transform_point(p)
        }
    }
    pub fn transform_vector(&self, time: Float, v: &Vector3f) -> Vector3f {
        if !self.actually_animated || time <= self.start_time {
            self.start_transform.transform_vector(v)
        } else if time >= self.end_time {
            self.end_transform.transform_vector(v)
        } else {
            let mut t: Transform = Transform::default();
            self.interpolate(time, &mut t);
            t.transform_vector(v)
        }
    }
    pub fn motion_bounds(&self, b: &Bounds3f) -> Bounds3f {
        if !self.actually_animated {
            return self.start_transform.transform_bounds(b);
        }
        if !self.has_rotation {
            return bnd3_union_bnd3f(
                &self.start_transform.transform_bounds(b),
                &self.end_transform.transform_bounds(b),
            );
        }
        // return motion bounds accounting for animated rotation
        let mut bounds: Bounds3f = Bounds3f::default();
        for corner in 0..8 {
            bounds = bnd3_union_bnd3f(&bounds, &self.bound_point_motion(&b.corner(corner)));
        }
        bounds
    }
    pub fn bound_point_motion(&self, p: &Point3f) -> Bounds3f {
        if !self.actually_animated {
            // set both to start
            let bounds: Bounds3f = Bounds3f::new(
                self.start_transform.transform_point(p),
                self.start_transform.transform_point(p),
            );
            return bounds;
        }
        let mut bounds: Bounds3f = Bounds3f::new(
            self.start_transform.transform_point(p),
            self.end_transform.transform_point(p),
        );
        let cos_theta: Float = quat_dot_quat(&self.r[0], &self.r[1]);
        let theta: Float = clamp_t(cos_theta, -1.0 as Float, 1.0 as Float).acos();
        for c in 0..3 {
            // find any motion derivative zeros for the component _c_
            let mut zeros: [Float; 8] = [0.0 as Float; 8];
            let mut n_zeros: u8 = 0;
            interval_find_zeros(
                self.c1[c].eval(p),
                self.c2[c].eval(p),
                self.c3[c].eval(p),
                self.c4[c].eval(p),
                self.c5[c].eval(p),
                theta,
                Interval::new(0.0 as Float, 1.0 as Float),
                &mut zeros,
                &mut n_zeros,
                8_usize,
            );
            // expand bounding box for any motion derivative zeros found
            for item in &zeros {
                let pz: Point3f =
                    self.transform_point(lerp(*item, self.start_time, self.end_time), p);
                bounds = bnd3_union_pnt3f(&bounds, &pz);
            }
        }
        bounds
    }
}

#[derive(Debug, Default, Copy, Clone)]
pub struct Interval {
    pub low: Float,
    pub high: Float,
}

impl Interval {
    pub fn new(v0: Float, v1: Float) -> Self {
        Interval {
            low: v0.min(v1),
            high: v0.max(v1),
        }
    }
}

impl Add for Interval {
    type Output = Interval;
    fn add(self, rhs: Interval) -> Interval {
        Interval {
            low: self.low + rhs.low,
            high: self.high + rhs.high,
        }
    }
}

impl Mul for Interval {
    type Output = Interval;
    fn mul(self, rhs: Interval) -> Interval {
        let min_rhs_low: Float = (self.low * rhs.low).min(self.high * rhs.low);
        let min_rhs_high: Float = (self.low * rhs.high).min(self.high * rhs.high);
        let max_rhs_low: Float = (self.low * rhs.low).max(self.high * rhs.low);
        let max_rhs_high: Float = (self.low * rhs.high).max(self.high * rhs.high);
        let low: Float = min_rhs_low.min(min_rhs_high);
        let high: Float = max_rhs_low.max(max_rhs_high);
        Interval { low, high }
    }
}

pub fn interval_sin(i: Interval) -> Interval {
    assert!(i.low >= 0.0 as Float);
    assert!(i.high <= 2.0001 as Float * PI);
    let mut sin_low: Float = i.low.sin();
    let mut sin_high: Float = i.high.sin();
    if sin_low > sin_high {
        std::mem::swap(&mut sin_low, &mut sin_high);
    }
    if i.low < PI / 2.0 as Float && i.high > PI / 2.0 as Float {
        sin_high = 1.0 as Float;
    }
    if i.low < (3.0 as Float / 2.0 as Float) * PI && i.high > (3.0 as Float / 2.0 as Float) * PI {
        sin_low = -1.0 as Float;
    }
    Interval {
        low: sin_low,
        high: sin_high,
    }
}

pub fn interval_cos(i: Interval) -> Interval {
    assert!(i.low >= 0.0 as Float);
    assert!(i.high <= 2.0001 as Float * PI);
    let mut cos_low: Float = i.low.cos();
    let mut cos_high: Float = i.high.cos();
    if cos_low > cos_high {
        std::mem::swap(&mut cos_low, &mut cos_high);
    }
    if i.low < PI && i.high > PI {
        cos_low = -1.0 as Float;
    }
    Interval {
        low: cos_low,
        high: cos_high,
    }
}

pub fn interval_find_zeros(
    c1: Float,
    c2: Float,
    c3: Float,
    c4: Float,
    c5: Float,
    theta: Float,
    t_interval: Interval,
    zeros: &mut [Float; 8],
    zero_count: &mut u8,
    depth: usize,
) {
    // evaluate motion derivative in interval form, return if no zeros
    let two_theta: Float = 2.0 as Float * theta;
    let range: Interval = Interval::new(c1, c1)
        + (Interval::new(c2, c2) + Interval::new(c3, c3) * t_interval)
            * interval_cos(Interval::new(two_theta, two_theta) * t_interval)
        + (Interval::new(c4, c4) + Interval::new(c5, c5) * t_interval)
            * interval_sin(Interval::new(two_theta, two_theta) * t_interval);
    if range.low > 0.0 as Float || range.high < 0.0 as Float || range.low == range.high {
        return;
    }
    if depth > 0_usize {
        // split _t_interval_ and check both resulting intervals
        let mid: Float = (t_interval.low + t_interval.high) * 0.5 as Float;
        interval_find_zeros(
            c1,
            c2,
            c3,
            c4,
            c5,
            theta,
            Interval::new(t_interval.low, mid),
            zeros,
            zero_count,
            depth - 1_usize,
        );
        interval_find_zeros(
            c1,
            c2,
            c3,
            c4,
            c5,
            theta,
            Interval::new(mid, t_interval.high),
            zeros,
            zero_count,
            depth - 1_usize,
        );
    } else {
        // use Newton's method to refine zero
        let mut t_newton: Float = (t_interval.low + t_interval.high) * 0.5 as Float;
        for _i in 0..4 {
            let f_newton: Float = c1
                + (c2 + c3 * t_newton) * (2.0 as Float * theta * t_newton).cos()
                + (c4 + c5 * t_newton) * (2.0 as Float * theta * t_newton).sin();
            let f_prime_newton: Float = (c3 + 2.0 as Float * (c4 + c5 * t_newton) * theta)
                * (2.0 as Float * t_newton * theta).cos()
                + (c5 - 2.0 as Float * (c2 + c3 * t_newton) * theta)
                    * (2.0 as Float * t_newton * theta).sin();
            if f_newton == 0.0 as Float || f_prime_newton == 0.0 as Float {
                break;
            }
            t_newton -= f_newton / f_prime_newton;
        }
        if t_newton >= t_interval.low - 1e-3 as Float && t_newton < t_interval.high + 1e-3 as Float
        {
            zeros[*zero_count as usize] = t_newton;
            *zero_count += 1;
        }
    }
}
