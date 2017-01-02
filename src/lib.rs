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
//! ### Translations
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
//!
//! ### Scaling
//!
//! Another basic transformations is the scale transformation. We can
//! differentiate between **uniform** scaling, where all three scale
//! factors have the same value, and **nonuniform** scaling, where
//! they may have different values.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::Transform;
//!
//! fn main() {
//!     let t: Transform = Transform::scale(2.0, -8.0, 1.0);
//!
//!     println!("t = {:?}", t);
//! }
//! ```
//!
//! ### X, Y, And Z Axis Rotations
//!
//! Another useful type of transformation is the rotation
//! transformation.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::{Float, Transform};
//!
//! fn main() {
//!     let theta_x: Float = 30.0;
//!     let theta_y: Float = 45.0;
//!     let theta_z: Float = 60.0;
//!     let t_x: Transform = Transform::rotate_x(theta_x);
//!     let t_y: Transform = Transform::rotate_y(theta_y);
//!     let t_z: Transform = Transform::rotate_z(theta_z);
//!
//!     println!("Transform::rotate_x({}) = {:?}", theta_x, t_x);
//!     println!("Transform::rotate_y({}) = {:?}", theta_y, t_y);
//!     println!("Transform::rotate_z({}) = {:?}", theta_z, t_z);
//! }
//! ```
//!
//! ### Rotation Around an Arbitrary Axis
//!
//! We also provide a routine to compute the transformation that
//! represents rotation around an arbitrary axis.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::{Float, Transform, Vector3f};
//!
//! fn main() {
//!     let theta: Float = 30.0;
//!     let axis = Vector3f {
//!         x: 1.0,
//!         y: 2.0,
//!         z: 3.0,
//!     };
//!     let t: Transform = Transform::rotate(theta, axis);
//!
//!     println!("Transform::rotate({}, {:?}) = {:?}", theta, axis, t);
//! }
//! ```
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
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::{Point3f, Transform, Vector3f};
//!
//! fn main() {
//!     // LookAt 2 2 5  0 -.4 0  0 1 0 (see spheres-differentials-texfilt.pbrt)
//!     let pos = Point3f {
//!         x: 2.0,
//!         y: 2.0,
//!         z: 5.0,
//!     };
//!     let look = Point3f {
//!         x: 0.0,
//!         y: -0.4,
//!         z: 0.0,
//!     };
//!     let up = Vector3f {
//!         x: 0.0,
//!         y: 1.0,
//!         z: 0.0,
//!     };
//!     let t: Transform = Transform::look_at(pos, look, up);
//!
//!     println!("Transform::look_at({:?}, {:?}, {:?}) = {:?}",
//!              pos,
//!              look,
//!              up,
//!              t);
//! }
//! ```
//!
//! ## Shapes
//!
//! Careful abstraction of geometric shapes in a ray tracer is a key
//! component of a clean system design, and shapes are the ideal
//! candidate for an object-oriented approach. All geometric
//! primitives implement a common interface, and the rest of the
//! renderer can use this interface without needing any details about
//! the underlying shape. This makes it possible to separate the
//! geometric and the shading subsystem of pbrt.
//!
//! ### Spheres
//!
//! Spheres are a special case of a general type of surfaces called
//! quadrics. They are the simplest type of curved surfaces that is
//! useful to a ray tracer and are a good starting point for general
//! ray intersection routines.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::Sphere;
//!
//! fn main() {
//!     let default_sphere: Sphere = Sphere::default();
//!     let sphere: Sphere = Sphere::new(2.0, -0.5, 0.75, 270.0);
//!
//!     println!("default sphere = {:?}", default_sphere);
//!     println!("sphere = {:?}", sphere);
//! }
//! ```
//!
//! ### Cones
//!
//! TODO
//!
//! ### Disks
//!
//! TODO
//!
//! ### Cylinders
//!
//! TODO
//!
//! ### Hyperboloids
//!
//! TODO
//!
//! ### Paraboloids
//!
//! TODO
//!

extern crate num;

use std::ops::{Add, Sub, Mul, Div};
use std::default::Default;
use std::f64::consts::PI;

pub type Float = f64;

// see pbrt.h

pub fn clamp<T>(val: T, low: T, high: T) -> T
    where T: PartialOrd
{
    let r: T;
    if val < low {
        r = low;
    } else if val > high {
        r = high;
    } else {
        r = val;
    }
    r
}

pub fn radians(deg: Float) -> Float {
    (PI / 180.0) * deg
}

pub fn degrees(rad: Float) -> Float {
    (180.0 / PI) * rad
}

// see geometry.h

pub type Point3f = Point3<Float>;
pub type Vector3f = Vector3<Float>;

#[derive(Debug,Copy,Clone)]
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

#[derive(Debug,Copy,Clone)]
pub struct Vector3<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T> Div<T> for Vector3<T>
    where T: Copy + Div<T, Output = T>
{
    type Output = Vector3<T>;
    fn div(self, rhs: T) -> Vector3<T>
        where T: Copy + Div<T, Output = T>
    {
        Vector3::<T> {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
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

pub fn vec3_cross<T>(v1: Vector3<T>, v2: Vector3<T>) -> Vector3<T>
    where T: Copy + Sub<T, Output = T> + Mul<T, Output = T>
{
    let v1x: T = v1.x;
    let v1y: T = v1.y;
    let v1z: T = v1.z;
    let v2x: T = v2.x;
    let v2y: T = v2.y;
    let v2z: T = v2.z;
    Vector3::<T> {
        x: (v1y * v2z) - (v1z * v2y),
        y: (v1z * v2x) - (v1x * v2z),
        z: (v1x * v2y) - (v1y * v2x),
    }
}

pub fn normalize<T>(v: Vector3<T>) -> Vector3<T>
    where T: num::Float + Copy + Div<T, Output = T>
{
    v / v.length()
}

#[derive(Debug,Copy,Clone)]
pub struct Point2<T> {
    pub x: T,
    pub y: T,
}

#[derive(Debug,Copy,Clone)]
pub struct Point3<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T> Sub<Point3<T>> for Point3<T>
    where T: Sub<T, Output = T>
{
    type Output = Vector3<T>;
    fn sub(self, rhs: Point3<T>) -> Vector3<T>
        where T: Sub<T, Output = T>
    {
        Vector3::<T> {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

#[derive(Debug,Copy,Clone)]
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

#[derive(Debug,Copy,Clone)]
pub struct Bounds2<T> {
    pub p_min: Point2<T>,
    pub p_max: Point2<T>,
}

#[derive(Debug,Copy,Clone)]
pub struct Bounds3<T> {
    pub p_min: Point3<T>,
    pub p_max: Point3<T>,
}

#[derive(Debug,Copy,Clone)]
pub struct Ray {
    /// origin
    pub o: Point3f,
    /// direction
    pub d: Vector3f,
}

// see transform.h

#[derive(Debug,Copy,Clone)]
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
    pub fn transpose(m: Matrix4x4) -> Matrix4x4 {
        Matrix4x4 {
            m: [[m.m[0][0], m.m[1][0], m.m[2][0], m.m[3][0]],
                [m.m[0][1], m.m[1][1], m.m[2][1], m.m[3][1]],
                [m.m[0][2], m.m[1][2], m.m[2][2], m.m[3][2]],
                [m.m[0][3], m.m[1][3], m.m[2][3], m.m[3][3]]],
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

#[derive(Debug,Copy,Clone)]
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
    pub fn scale(x: Float, y: Float, z: Float) -> Transform {
        Transform {
            m: Matrix4x4::new(x,
                              0.0,
                              0.0,
                              0.0,
                              0.0,
                              y,
                              0.0,
                              0.0,
                              0.0,
                              0.0,
                              z,
                              0.0,
                              0.0,
                              0.0,
                              0.0,
                              1.0),
            m_inv: Matrix4x4::new(1.0 / x,
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
                                  1.0),
        }
    }
    pub fn rotate_x(theta: Float) -> Transform {
        let sin_theta: Float = radians(theta).sin();
        let cos_theta: Float = radians(theta).cos();
        let m = Matrix4x4::new(1.0,
                               0.0,
                               0.0,
                               0.0,
                               0.0,
                               cos_theta,
                               -sin_theta,
                               0.0,
                               0.0,
                               sin_theta,
                               cos_theta,
                               0.0,
                               0.0,
                               0.0,
                               0.0,
                               1.0);
        Transform {
            m: m,
            m_inv: Matrix4x4::transpose(m),
        }
    }
    pub fn rotate_y(theta: Float) -> Transform {
        let sin_theta: Float = radians(theta).sin();
        let cos_theta: Float = radians(theta).cos();
        let m = Matrix4x4::new(cos_theta,
                               0.0,
                               sin_theta,
                               0.0,
                               0.0,
                               1.0,
                               0.0,
                               0.0,
                               -sin_theta,
                               0.0,
                               cos_theta,
                               0.0,
                               0.0,
                               0.0,
                               0.0,
                               1.0);
        Transform {
            m: m,
            m_inv: Matrix4x4::transpose(m),
        }
    }
    pub fn rotate_z(theta: Float) -> Transform {
        let sin_theta: Float = radians(theta).sin();
        let cos_theta: Float = radians(theta).cos();
        let m = Matrix4x4::new(cos_theta,
                               -sin_theta,
                               0.0,
                               0.0,
                               sin_theta,
                               cos_theta,
                               0.0,
                               0.0,
                               0.0,
                               0.0,
                               1.0,
                               0.0,
                               0.0,
                               0.0,
                               0.0,
                               1.0);
        Transform {
            m: m,
            m_inv: Matrix4x4::transpose(m),
        }
    }
    pub fn rotate(theta: Float, axis: Vector3f) -> Transform {
        let a: Vector3f = normalize(axis);
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
            m: m,
            m_inv: Matrix4x4::transpose(m),
        }
    }
    pub fn look_at(pos: Point3f, look: Point3f, up: Vector3f) -> Transform {
        let mut camera_to_world = Matrix4x4::default();
        // initialize fourth column of viewing matrix
        camera_to_world.m[0][3] = pos.x;
        camera_to_world.m[1][3] = pos.y;
        camera_to_world.m[2][3] = pos.z;
        camera_to_world.m[3][3] = 1.0;
        // initialize first three columns of viewing matrix
        let dir: Vector3f = normalize(look - pos);
        if vec3_cross(normalize(up), dir).length() == 0.0 {
            println!("\"up\" vector ({}, {}, {}) and viewing direction ({}, {}, {}) passed to \
                      LookAt are pointing in the same direction.  Using the identity \
                      transformation.",
                     up.x,
                     up.y,
                     up.z,
                     dir.x,
                     dir.y,
                     dir.z);
            Transform::default()
        } else {
            let left: Vector3f = normalize(vec3_cross(normalize(up), dir));
            let new_up: Vector3f = vec3_cross(dir, left);
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
                m: Matrix4x4::inverse(camera_to_world),
                m_inv: camera_to_world,
            }
        }
    }
}

// see sphere.h

#[derive(Debug,Copy,Clone)]
pub struct Sphere {
    radius: Float,
    z_min: Float,
    z_max: Float,
    theta_min: Float,
    theta_max: Float,
    phi_max: Float,
}

impl Default for Sphere {
    fn default() -> Sphere {
        Sphere {
            radius: 1.0,
            z_min: -1.0,
            z_max: 1.0,
            theta_min: (-1.0 as Float).acos(),
            theta_max: (1.0 as Float).acos(),
            phi_max: 360.0,
        }
    }
}

impl Sphere {
    pub fn new(radius: Float, z_min: Float, z_max: Float, phi_max: Float) -> Sphere {
        Sphere {
            radius: radius,
            z_min: clamp(z_min.min(z_max), -radius, radius),
            z_max: clamp(z_min.max(z_max), -radius, radius),
            theta_min: clamp(z_min.min(z_max) / radius, -1.0, 1.0).acos(),
            theta_max: clamp(z_min.max(z_max) / radius, -1.0, 1.0).acos(),
            phi_max: phi_max,
        }
    }
}
