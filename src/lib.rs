//! # pbrt
//!
//! ## Vectors
//!
//! **pbrt** provides both 2D and 3D **vector** classes. Both are
//! parameterized by the type of the underlying vector element, thus
//! making it easy to instantiate vectors of both integer and
//! floating-point types.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::Vector3;
//!
//! fn main() {
//!     let int_null = Vector3 { x: 0, y: 0, z: 0 };
//!     let float_null = Vector3 {
//!         x: 0.0,
//!         y: 0.0,
//!         z: 0.0,
//!     };
//!
//!     println!("int   {:?}", int_null);
//!     println!("float {:?}", float_null);
//! }
//! ```
//!
//! ## Points
//!
//! A **point** is a zero-dimensional location in 2D or 3D space. The
//! **Point2** and **Point3** classes in **pbrt** represent points in
//! the obvious way: using x, y, z (in 3D) coordinates with respect to
//! a coordinate system. Although the same representation is used for
//! vectors, the fact that a point represents a position whereas a
//! vector represents a direction leads to a number of important
//! differences in how they are treated.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::Point3;
//!
//! fn main() {
//!     let int_origin = Point3 { x: 0, y: 0, z: 0 };
//!     let float_origin = Point3 {
//!         x: 0.0,
//!         y: 0.0,
//!         z: 0.0,
//!     };
//!
//!     println!("int   {:?}", int_origin);
//!     println!("float {:?}", float_origin);
//! }
//! ```
//!
//! ## Normals
//!
//! A surface **normal** (or just normal) is a vector that is
//! perpendicular to a surface at a particular position. It can be
//! defined as the cross product of any two nonparallel vectors that
//! are tangent to the surface at a point. Although normals are
//! superficially similar to vectors, it is important to distinguish
//! between the two of them: because normals are defined in terms of
//! their relationship to a particular surface, they behave
//! differently than vectors in some situations, particularly when
//! applying transformations.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::Normal3;
//!
//! fn main() {
//!     let int_null = Normal3 { x: 0, y: 0, z: 0 };
//!     let float_null = Normal3 {
//!         x: 0.0,
//!         y: 0.0,
//!         z: 0.0,
//!     };
//!
//!     println!("int   {:?}", int_null);
//!     println!("float {:?}", float_null);
//! }
//! ```
//!
//! ## Rays
//!
//! A **ray** is a semi-infinite line specified by its origin and
//! direction. **pbrt** represents a **Ray** with a **Point3f** for
//! the origin and a **Vector3f** for the direction. We only need rays
//! with floating-point origins and directions, so **Ray** isn't a
//! template class parameterized by an arbitrary type, as points,
//! vectors, and normals were.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::{Ray, Point3f, Vector3f};
//!
//! fn main() {
//!     let origin = Point3f {
//!         x: -5.5,
//!         y: 2.75,
//!         z: 0.0,
//!     };
//!     let direction = Vector3f {
//!         x: 1.0,
//!         y: -8.75,
//!         z: 2.25,
//!     };
//!     let ray = Ray {
//!         o: origin,
//!         d: direction,
//!     };
//!
//!     println!("{:?}", ray);
//! }
//! ```
//!
//! ## Bounding Boxes
//!
//! Many parts of the system operate on axis-aligned regions of
//! space. For example, multi-threading in **pbrt** is implemented by
//! subdividing the image into rectangular tiles that can be processed
//! independently, and the bounding volume hierarchy uses 3D boxes to
//! bound geometric primitives in the scene. The **Bounds2** and
//! **Bounds3** template classes are used to represent the extent of
//! these sort of regions. Both are parameterized by a type T that is
//! used to represent the coordinates of its extents.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::{Bounds3, Point3};
//!
//! fn main() {
//!     let int_origin = Point3 { x: 0, y: 0, z: 0 };
//!     let int_xyz111 = Point3 { x: 1, y: 1, z: 1 };
//!     let float_origin = Point3 {
//!         x: 0.0,
//!         y: 0.0,
//!         z: 0.0,
//!     };
//!     let float_xyz111 = Point3 {
//!         x: 1.0,
//!         y: 1.0,
//!         z: 1.0,
//!     };
//!     let int_unit_cube = Bounds3 {
//!         p_min: int_origin,
//!         p_max: int_xyz111,
//!     };
//!     let float_unit_cube = Bounds3 {
//!         p_min: float_origin,
//!         p_max: float_xyz111,
//!     };
//!
//!     println!("int   {:?}", int_unit_cube);
//!     println!("float {:?}", float_unit_cube);
//! }
//! ```
//!
//! ## 4 x 4 Matrices
//!
//! The **Matrix4x4** structure provides a low-level representation of
//! 4 x 4 matrices. It is an integral part of the **Transform** class.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::Matrix4x4;
//!
//! fn main() {
//!     let identity: Matrix4x4 = Matrix4x4::default();
//!     let m: Matrix4x4 = Matrix4x4::new(0.0,
//!                                       0.1,
//!                                       0.2,
//!                                       0.3,
//!                                       1.0,
//!                                       1.1,
//!                                       1.2,
//!                                       1.3,
//!                                       2.0,
//!                                       2.1,
//!                                       2.2,
//!                                       2.3,
//!                                       3.0,
//!                                       3.1,
//!                                       3.2,
//!                                       3.3);
//!
//!     println!("identity matrix = {:?}", identity);
//!     println!("m = {:?}", m);
//! }
//! ```
//!
//! ## Transformations
//!
//! In general a transformation is a mapping from points to points and
//! from vectors to vectors. When a new **Transform** is created, it
//! defaults to the *identity transformation* - the transformation
//! that maps each point and each vector to itself.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::Transform;
//!
//! fn main() {
//!     let identity: Transform = Transform::default();
//!     let t: Transform = Transform::new(2.0,
//!                                       0.0,
//!                                       0.0,
//!                                       -1.25,
//!                                       0.0,
//!                                       -8.0,
//!                                       0.0,
//!                                       3.5,
//!                                       0.0,
//!                                       0.0,
//!                                       1.0,
//!                                       7.875,
//!                                       0.0,
//!                                       0.0,
//!                                       0.0,
//!                                       1.0);
//!
//!     println!("identity matrix = {:?}", identity);
//!     println!("t = {:?}", t);
//! }
//! ```
//!
//! ## Translations
//!
//! One of the simplest transformations is the translation
//! transformation. Translations only affect points, leaving vectors
//! unchanged.
//!
//! ```rust
//! extern crate pbrt;
//! 
//! use pbrt::{Transform, Vector3f};
//! 
//! fn main() {
//!     let t: Transform = Transform::translate(Vector3f {
//!         x: -1.25,
//!         y: 3.5,
//!         z: 7.875,
//!     });
//! 
//!     println!("t = {:?}", t);
//! }
//! ```

extern crate num;

use std::ops::{Add, Mul};
use std::default::Default;

pub type Float = f64;

// see geometry.h

pub type Point3f = Point3<Float>;
pub type Vector3f = Vector3<Float>;

#[derive(Debug)]
pub struct Vector2<T> {
    pub x: T,
    pub y: T,
}

impl<T> Vector2<T> {
    pub fn length_squared(&self) -> T
        where T: Copy + Add<T, Output = T> + Mul<T, Output = T>
    {
        self.x * self.x + self.y * self.y
    }
    pub fn length(&self) -> T
        where T: num::Float
    {
        self.length_squared().sqrt()
    }
}

pub fn vec2_dot<T>(v1: Vector2<T>, v2: Vector2<T>) -> T
    where T: Copy + Add<T, Output = T> + Mul<T, Output = T>
{
    v1.x * v2.x + v1.y * v2.y
}

#[derive(Debug)]
pub struct Vector3<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T> Vector3<T> {
    pub fn length_squared(&self) -> T
        where T: Copy + Add<T, Output = T> + Mul<T, Output = T>
    {
        self.x * self.x + self.y * self.y + self.z * self.z
    }
    pub fn length(&self) -> T
        where T: num::Float
    {
        self.length_squared().sqrt()
    }
}

pub fn vec3_dot<T>(v1: Vector3<T>, v2: Vector3<T>) -> T
    where T: Copy + Add<T, Output = T> + Mul<T, Output = T>
{
    v1.x * v2.x + v1.y * v2.y + v1.z * v2.z
}

#[derive(Debug)]
pub struct Point2<T> {
    pub x: T,
    pub y: T,
}

#[derive(Debug)]
pub struct Point3<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

#[derive(Debug)]
pub struct Normal3<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T> Normal3<T> {
    pub fn length_squared(&self) -> T
        where T: Copy + Add<T, Output = T> + Mul<T, Output = T>
    {
        self.x * self.x + self.y * self.y + self.z * self.z
    }
    pub fn length(&self) -> T
        where T: num::Float
    {
        self.length_squared().sqrt()
    }
}

#[derive(Debug)]
pub struct Bounds2<T> {
    pub p_min: Point2<T>,
    pub p_max: Point2<T>,
}

#[derive(Debug)]
pub struct Bounds3<T> {
    pub p_min: Point3<T>,
    pub p_max: Point3<T>,
}

#[derive(Debug)]
pub struct Ray {
    /// origin
    pub o: Point3f,
    /// direction
    pub d: Vector3f,
}

// see transform.h

#[derive(Debug)]
pub struct Matrix4x4 {
    pub m: [[Float; 4]; 4],
}

impl Default for Matrix4x4 {
    fn default() -> Matrix4x4 {
        Matrix4x4 {
            m: [[1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0]],
        }
    }
}

impl Matrix4x4 {
    pub fn new(t00: Float,
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
               t33: Float)
               -> Matrix4x4 {
        Matrix4x4 {
            m: [[t00, t01, t02, t03],
                [t10, t11, t12, t13],
                [t20, t21, t22, t23],
                [t30, t31, t32, t33]],
        }
    }
    pub fn inverse(m: Matrix4x4) -> Matrix4x4 {
        let mut indxc = vec![0; 4];
        let mut indxr = vec![0; 4];
        let mut ipiv = vec![0; 4];
        let mut minv: Matrix4x4 = Matrix4x4::new(m.m[0][0],
                                                 m.m[0][1],
                                                 m.m[0][2],
                                                 m.m[0][3],
                                                 m.m[1][0],
                                                 m.m[1][1],
                                                 m.m[1][2],
                                                 m.m[1][3],
                                                 m.m[2][0],
                                                 m.m[2][1],
                                                 m.m[2][2],
                                                 m.m[2][3],
                                                 m.m[3][0],
                                                 m.m[3][1],
                                                 m.m[3][2],
                                                 m.m[3][3]);
        for i in 0..4 {
            let mut irow = 0;
            let mut icol = 0;
            let mut big: Float = 0.0;
            // choose pivot
            for j in 0..4 {
                if ipiv[j] != 1 {
                    for k in 0..4 {
                        if ipiv[k] == 0 {
                            let abs: Float = (minv.m[j][k]).abs();
                            if abs >= big {
                                big = abs;
                                irow = j;
                                icol = k;
                            }
                        } else {
                            if ipiv[k] > 1 {
                                println!("Singular matrix in MatrixInvert");
                            }
                        }
                    }
                }
            }
            ipiv[icol] = ipiv[icol] + 1;
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
                minv.m[icol][j] = minv.m[icol][j] * pivinv;
            }
            // subtract this row from others to zero out their columns
            for j in 0..4 {
                if j != icol {
                    let save: Float = minv.m[j][icol];
                    minv.m[j][icol] = 0.0;
                    for k in 0..4 {
                        minv.m[j][k] = minv.m[j][k] - (minv.m[icol][k] * save);
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
                    let swap = minv.m[k][indxr[j]];
                    minv.m[k][indxr[j]] = minv.m[k][indxc[j]];
                    minv.m[k][indxc[j]] = swap;
                }
            }
        }
        minv
    }
}

#[derive(Debug)]
pub struct Transform {
    pub m: Matrix4x4,
    pub m_inv: Matrix4x4,
}

impl Default for Transform {
    fn default() -> Transform {
        Transform {
            m: Matrix4x4::default(),
            m_inv: Matrix4x4::default(),
        }
    }
}

impl Transform {
    pub fn new(t00: Float,
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
               t33: Float)
               -> Transform {
        Transform {
            m: Matrix4x4::new(t00,
                              t01,
                              t02,
                              t03,
                              t10,
                              t11,
                              t12,
                              t13,
                              t20,
                              t21,
                              t22,
                              t23,
                              t30,
                              t31,
                              t32,
                              t33),
            m_inv: Matrix4x4::inverse(Matrix4x4::new(t00,
                                                     t01,
                                                     t02,
                                                     t03,
                                                     t10,
                                                     t11,
                                                     t12,
                                                     t13,
                                                     t20,
                                                     t21,
                                                     t22,
                                                     t23,
                                                     t30,
                                                     t31,
                                                     t32,
                                                     t33)),
        }
    }
    pub fn translate(delta: Vector3f) -> Transform {
        Transform {
            m: Matrix4x4::new(1.0,
                              0.0,
                              0.0,
                              delta.x,
                              0.0,
                              1.0,
                              0.0,
                              delta.y,
                              0.0,
                              0.0,
                              1.0,
                              delta.z,
                              0.0,
                              0.0,
                              0.0,
                              1.0),
            m_inv: Matrix4x4::new(1.0,
                                  0.0,
                                  0.0,
                                  -delta.x,
                                  0.0,
                                  1.0,
                                  0.0,
                                  -delta.y,
                                  0.0,
                                  0.0,
                                  1.0,
                                  -delta.z,
                                  0.0,
                                  0.0,
                                  0.0,
                                  1.0),
        }
    }
}
