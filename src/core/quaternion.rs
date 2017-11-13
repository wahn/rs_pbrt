// std
use std::ops::{Add, Div, Mul, Neg, Sub};
// pbrt
use core::pbrt::Float;
use core::pbrt::clamp_t;
use core::transform::{Matrix4x4, Transform};
use geometry::Vector3f;
use geometry::vec3_dot_vec3;

// see quaternion.h

#[derive(Debug,Copy,Clone)]
pub struct Quaternion {
    pub v: Vector3f,
    pub w: Float,
}

impl Default for Quaternion {
    fn default() -> Self {
        Quaternion {
            v: Vector3f::default(),
            w: 1.0,
        }
    }
}

impl Quaternion {
    pub fn new(t: Transform) -> Self {
        let m: Matrix4x4 = t.m.clone();
        let trace: Float = m.m[0][0] + m.m[1][1] + m.m[2][2];
        // if (trace > 0.f) {
        if trace > 0.0 {
            // compute w from matrix trace, then xyz
            // 4w^2 = m[0][0] + m[1][1] + m[2][2] + m[3][3] (but m[3][3] == 1)
            let mut s: Float = (trace + 1.0).sqrt();
            let w: Float = s / 2.0;
            s = 0.5 / s;
            Quaternion {
                v: Vector3f {
                    x: (m.m[2][1] - m.m[1][2]) * s,
                    y: (m.m[0][2] - m.m[2][0]) * s,
                    z: (m.m[1][0] - m.m[0][1]) * s,
                },
                w: w,
            }
        } else {
            // compute largest of $x$, $y$, or $z$, then remaining components
            let nxt: [usize; 3] = [1, 2, 0];
            let mut q: [Float; 3] = [0.0; 3];
            let mut i: usize = 0;
            if m.m[1][1] > m.m[0][0] {
                i = 1;
            }
            if m.m[2][2] > m.m[i][i] {
                i = 2;
            }
            let j = nxt[i];
            let k = nxt[j];
            let mut s: Float = ((m.m[i][i] - (m.m[j][j] + m.m[k][k])) + 1.0).sqrt();
            q[i] = s * 0.5;
            if s != 0.0 {
                s = 0.5 / s;
            }
            let w: Float = (m.m[k][j] - m.m[j][k]) * s;
            q[j] = (m.m[j][i] + m.m[i][j]) * s;
            q[k] = (m.m[k][i] + m.m[i][k]) * s;
            Quaternion {
                v: Vector3f {
                    x: q[0],
                    y: q[1],
                    z: q[2],
                },
                w: w,
            }
        }
    }
    pub fn to_transform(&self) -> Transform {
        let xx: Float = self.v.x * self.v.x;
        let yy: Float = self.v.y * self.v.y;
        let zz: Float = self.v.z * self.v.z;
        let xy: Float = self.v.x * self.v.y;
        let xz: Float = self.v.x * self.v.z;
        let yz: Float = self.v.y * self.v.z;
        let wx: Float = self.v.x * self.w;
        let wy: Float = self.v.y * self.w;
        let wz: Float = self.v.z * self.w;

        let mut m: Matrix4x4 = Matrix4x4::default();
        m.m[0][0] = 1.0 as Float - 2.0 as Float * (yy + zz);
        m.m[0][1] = 2.0 as Float * (xy + wz);
        m.m[0][2] = 2.0 as Float * (xz - wy);
        m.m[1][0] = 2.0 as Float * (xy - wz);
        m.m[1][1] = 1.0 as Float - 2.0 as Float * (xx + zz);
        m.m[1][2] = 2.0 as Float * (yz + wx);
        m.m[2][0] = 2.0 as Float * (xz + wy);
        m.m[2][1] = 2.0 as Float * (yz - wx);
        m.m[2][2] = 1.0 as Float - 2.0 as Float * (xx + yy);

        // transpose since we are left-handed
        Transform {
            m: Matrix4x4::transpose(m),
            m_inv: m,
        }
    }
}

impl Add for Quaternion {
    type Output = Quaternion;
    fn add(self, rhs: Quaternion) -> Quaternion {
        Quaternion {
            v: self.v + rhs.v,
            w: self.w + rhs.w,
        }
    }
}

impl Sub for Quaternion {
    type Output = Quaternion;
    fn sub(self, rhs: Quaternion) -> Quaternion {
        Quaternion {
            v: self.v - rhs.v,
            w: self.w - rhs.w,
        }
    }
}

impl Mul<Float> for Quaternion {
    type Output = Quaternion;
    fn mul(self, rhs: Float) -> Quaternion {
        Quaternion {
            v: self.v * rhs,
            w: self.w * rhs,
        }
    }
}

impl Div<Float> for Quaternion {
    type Output = Quaternion;
    fn div(self, rhs: Float) -> Quaternion {
        Quaternion {
            v: self.v / rhs,
            w: self.w / rhs,
        }
    }
}

impl Neg for Quaternion {
    type Output = Quaternion;
    fn neg(self) -> Quaternion {
        Quaternion {
            v: -self.v,
            w: -self.w,
        }
    }
}

pub fn quat_slerp(t: Float, q1: Quaternion, q2: Quaternion) -> Quaternion {
    let cos_theta: Float = quat_dot_quat(q1, q2);
    if cos_theta > 0.9995 as Float {
        quat_normalize(q1 * (1.0 as Float - t) + q2 * t)
    } else {
        let theta: Float = clamp_t(cos_theta, -1.0 as Float, 1.0 as Float).acos();
        let thetap: Float = theta * t;
        let qperp: Quaternion = quat_normalize(q2 - q1 * cos_theta);
        q1 * thetap.cos() + qperp * thetap.sin()
    }
}

/// The inner product of two quaterions.
pub fn quat_dot_quat(q1: Quaternion, q2: Quaternion) -> Float {
    vec3_dot_vec3(q1.v, q2.v) + q1.w * q2.w
}

/// A quaternion can be normalized by dividing by its length.
pub fn quat_normalize(q: Quaternion) -> Quaternion {
    q / quat_dot_quat(q, q).sqrt()
}
