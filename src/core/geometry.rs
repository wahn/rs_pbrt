//! Almost all nontrivial graphics programs are built on a foundation
//! of geometric classes. These classes represent mathematical
//! constructs like points, vectors, and rays.
//!
//! # Points
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
//! use rs_pbrt::core::geometry::{Point3i, Point3f};
//!
//!     let int_origin = Point3i { x: 0, y: 0, z: 0 };
//!     let float_origin = Point3f {
//!         x: 0.0,
//!         y: 0.0,
//!         z: 0.0,
//!     };
//!
//!     println!("int   {:?}", int_origin);
//!     println!("float {:?}", float_origin);
//! ```
//!
//! # Vectors
//!
//! **pbrt** provides both 2D and 3D **vector** classes. Both are
//! parameterized by the type of the underlying vector element, thus
//! making it easy to instantiate vectors of both integer and
//! floating-point types.
//!
//! ```rust
//! use rs_pbrt::core::geometry::{Vector3i, Vector3f};
//!
//!     let int_null = Vector3i { x: 0, y: 0, z: 0 };
//!     let float_null = Vector3f {
//!         x: 0.0,
//!         y: 0.0,
//!         z: 0.0,
//!     };
//!
//!     println!("int   {:?}", int_null);
//!     println!("float {:?}", float_null);
//! ```
//!
//! # Normals
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
//! use rs_pbrt::core::geometry::Normal3f;
//!
//!     let float_null = Normal3f {
//!         x: 0.0,
//!         y: 0.0,
//!         z: 0.0,
//!     };
//!
//!     println!("float {:?}", float_null);
//! ```
//!
//! # Rays
//!
//! A **ray** is a semi-infinite line specified by its origin and
//! direction. **pbrt** represents a **Ray** with a **Point3f** for
//! the origin and a **Vector3f** for the direction. We only need rays
//! with floating-point origins and directions, so **Ray** isn't a
//! template class parameterized by an arbitrary type, as points,
//! vectors, and normals were.
//!
//! ```rust
//! use rs_pbrt::core::geometry::{Ray, Point3f, Vector3f};
//! use std::cell::Cell;
//!
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
//!         t_max: Cell::new(std::f32::INFINITY),
//!         time: 0.0,
//!         medium: None,
//!         differential: None,
//!     };
//! ```
//!
//! ## RayDifferentials
//!
//! **RayDifferential** is a subclass of **Ray** that contains
//! additional information about two auxiliary rays. These extra rays
//! represent camera rays offset by one sample in the *x* and *y*
//! direction from the main ray on the film plane. By determining the
//! area that these three rays project on an object being shaded, the
//! **Texture** can estimate an area to average over for proper
//! antialiasing.
//!
//! In Rust we don't have inheritance, therefore we use an
//! **`Option<RayDifferential>`** in the **Ray** struct, which means the
//! additional information can be present (or not).
//!
//! # Bounding Boxes
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
//!
//! ```rust
//! use rs_pbrt::core::geometry::{Bounds3i, Bounds3f, Point3i, Point3f};
//!
//!     let int_origin = Point3i { x: 0, y: 0, z: 0 };
//!     let int_xyz111 = Point3i { x: 1, y: 1, z: 1 };
//!     let float_origin = Point3f {
//!         x: 0.0,
//!         y: 0.0,
//!         z: 0.0,
//!     };
//!     let float_xyz111 = Point3f {
//!         x: 1.0,
//!         y: 1.0,
//!         z: 1.0,
//!     };
//!     let int_unit_cube = Bounds3i {
//!         p_min: int_origin,
//!         p_max: int_xyz111,
//!     };
//!     let float_unit_cube = Bounds3f {
//!         p_min: float_origin,
//!         p_max: float_xyz111,
//!     };
//!
//!     println!("int   {:?}", int_unit_cube);
//!     println!("float {:?}", float_unit_cube);
//! ```

// std
use std::cell::Cell;
use std::f32::consts::PI;
use std::ops;
use std::ops::{Index, IndexMut};
use std::sync::Arc;
// others
use strum::IntoEnumIterator;
use strum_macros::EnumIter;
// pbrt
use crate::core::medium::Medium;
use crate::core::pbrt::Float;
use crate::core::pbrt::{clamp_t, gamma, lerp, next_float_down, next_float_up};

// see geometry.h

// pub type Point2f = Point2<Float>;
// pub type Point2i = Point2<i32>;
// pub type Point3f = Point3<Float>;
// pub type Point3i = Point3<i32>;
// pub type Vector2f = Vector2<Float>;
// pub type Vector2i = Vector2<i32>;
// pub type Vector3f = Vector3<Float>;
// pub type Vector3i = Vector3<i32>;
// pub type Normal3f = Normal3<Float>;

#[derive(EnumIter, Debug, Copy, Clone)]
#[repr(u8)]
pub enum XYEnum {
    X = 0,
    Y = 1,
}

#[derive(EnumIter, Debug, Copy, Clone)]
#[repr(u8)]
pub enum MinMaxEnum {
    Min = 0,
    Max = 1,
}

#[derive(EnumIter, Debug, Copy, Clone)]
#[repr(u8)]
pub enum XYZEnum {
    X = 0,
    Y = 1,
    Z = 2,
}

#[derive(Debug, Default, Copy, Clone)]
pub struct Vector2f {
    pub x: Float,
    pub y: Float,
}

impl Vector2f {
    pub fn has_nans(&self) -> bool {
        self.x.is_nan() || self.y.is_nan()
    }
    pub fn length_squared(&self) -> Float {
        self.x * self.x + self.y * self.y
    }
    pub fn length(&self) -> Float {
        self.length_squared().sqrt()
    }
}

#[derive(Debug, Default, Copy, Clone)]
pub struct Vector2i {
    pub x: i32,
    pub y: i32,
}

impl Vector2i {
    pub fn length_squared(&self) -> i32 {
        self.x * self.x + self.y * self.y
    }
}

impl Index<XYEnum> for Vector2f {
    type Output = Float;
    fn index(&self, index: XYEnum) -> &Float {
        match index {
            XYEnum::X => &self.x,
            _ => &self.y,
        }
    }
}

impl Index<XYEnum> for Vector2i {
    type Output = i32;
    fn index(&self, index: XYEnum) -> &i32 {
        match index {
            XYEnum::X => &self.x,
            _ => &self.y,
        }
    }
}

impl IndexMut<XYEnum> for Vector2f {
    fn index_mut(&mut self, index: XYEnum) -> &mut Float {
        match index {
            XYEnum::X => &mut self.x,
            _ => &mut self.y,
        }
    }
}

impl IndexMut<XYEnum> for Vector2i {
    fn index_mut(&mut self, index: XYEnum) -> &mut i32 {
        match index {
            XYEnum::X => &mut self.x,
            _ => &mut self.y,
        }
    }
}

// impl<T> AddAssign<Vector2<T>> for Vector2<T>
// where
//     T: AddAssign,
// {
//     fn add_assign(&mut self, rhs: Vector2<T>) {
//         self.x += rhs.x;
//         self.y += rhs.y;
//     }
// }

// impl<T> Add for Vector2<T>
// where
//     T: Copy + Add<T, Output = T>,
// {
//     type Output = Vector2<T>;
//     fn add(self, rhs: Vector2<T>) -> Vector2<T> {
//         Vector2::<T> {
//             x: self.x + rhs.x,
//             y: self.y + rhs.y,
//         }
//     }
// }

// impl<T> SubAssign<Vector2<T>> for Vector2<T>
// where
//     T: SubAssign,
// {
//     fn sub_assign(&mut self, rhs: Vector2<T>) {
//         self.x -= rhs.x;
//         self.y -= rhs.y;
//     }
// }

// impl<T> Sub for Vector2<T>
// where
//     T: Copy + Sub<T, Output = T>,
// {
//     type Output = Vector2<T>;
//     fn sub(self, rhs: Vector2<T>) -> Vector2<T> {
//         Vector2::<T> {
//             x: self.x - rhs.x,
//             y: self.y - rhs.y,
//         }
//     }
// }

// impl<T> MulAssign<T> for Vector2<T>
// where
//     T: Copy + MulAssign,
// {
//     fn mul_assign(&mut self, rhs: T) {
//         self.x *= rhs;
//         self.y *= rhs;
//     }
// }

// work around bug
// https://github.com/rust-lang/rust/issues/40395
// impl Div<Float> for Vector2<f32> {
//     type Output = Vector2<f32>;
//     fn div(self, rhs: Float) -> Vector2<f32> {
//         assert_ne!(rhs, 0.0 as Float);
//         let inv: Float = 1.0 as Float / rhs;
//         Vector2::<f32> {
//             x: self.x * inv,
//             y: self.y * inv,
//         }
//     }
// }

// work around bug
// https://github.com/rust-lang/rust/issues/40395
// impl DivAssign<Float> for Vector2<f32> {
//     fn div_assign(&mut self, rhs: Float) {
//         assert_ne!(rhs, 0.0 as Float);
//         let inv: Float = 1.0 as Float / rhs;
//         self.x *= inv;
//         self.y *= inv;
//     }
// }

impl From<Vector2f> for Point2f {
    fn from(v: Vector2f) -> Self {
        Point2f { x: v.x, y: v.y }
    }
}

impl From<Vector2i> for Point2i {
    fn from(v: Vector2i) -> Self {
        Point2i { x: v.x, y: v.y }
    }
}

/// Product of the Euclidean magnitudes of the two vectors and the
/// cosine of the angle between them. A return value of zero means
/// both vectors are orthogonal, a value if one means they are
/// codirectional.
pub fn vec2_dotf(v1: &Vector2f, v2: &Vector2f) -> Float {
    v1.x * v2.x + v1.y * v2.y
}

/// Product of the Euclidean magnitudes of the two vectors and the
/// cosine of the angle between them. A return value of zero means
/// both vectors are orthogonal, a value if one means they are
/// codirectional.
pub fn vec2_doti(v1: &Vector2i, v2: &Vector2i) -> i32 {
    v1.x * v2.x + v1.y * v2.y
}

#[derive(Debug, Default, Copy, Clone, PartialEq)]
pub struct Vector3f {
    pub x: Float,
    pub y: Float,
    pub z: Float,
}

impl Vector3f {
    pub fn has_nans(&self) -> bool {
        self.x.is_nan() || self.y.is_nan() || self.z.is_nan()
    }
    pub fn abs(&self) -> Vector3f {
        Vector3f {
            x: self.x.abs(),
            y: self.y.abs(),
            z: self.z.abs(),
        }
    }
    pub fn length_squared(&self) -> Float {
        self.x * self.x + self.y * self.y + self.z * self.z
    }
    pub fn length(&self) -> Float {
        self.length_squared().sqrt()
    }
    /// Compute a new vector pointing in the same direction but with unit
    /// length.
    pub fn normalize(&self) -> Vector3f {
        *self / self.length()
    }
}

#[derive(Debug, Default, Copy, Clone, PartialEq)]
pub struct Vector3i {
    pub x: i32,
    pub y: i32,
    pub z: i32,
}

impl Vector3i {
    pub fn abs(&self) -> Vector3i {
        Vector3i {
            x: self.x.abs(),
            y: self.y.abs(),
            z: self.z.abs(),
        }
    }
    pub fn length_squared(&self) -> i32 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }
}

// impl<T> AddAssign<Vector3<T>> for Vector3<T>
// where
//     T: AddAssign,
// {
//     fn add_assign(&mut self, rhs: Vector3<T>) {
//         self.x += rhs.x;
//         self.y += rhs.y;
//         self.z += rhs.z;
//     }
// }

// impl<T> Add for Vector3<T>
// where
//     T: Copy + Add<T, Output = T>,
// {
//     type Output = Vector3<T>;
//     fn add(self, rhs: Vector3<T>) -> Vector3<T> {
//         Vector3::<T> {
//             x: self.x + rhs.x,
//             y: self.y + rhs.y,
//             z: self.z + rhs.z,
//         }
//     }
// }

// impl<T> SubAssign<Vector3<T>> for Vector3<T>
// where
//     T: SubAssign,
// {
//     fn sub_assign(&mut self, rhs: Vector3<T>) {
//         self.x -= rhs.x;
//         self.y -= rhs.y;
//         self.z -= rhs.z;
//     }
// }

// impl<T> Sub for Vector3<T>
// where
//     T: Copy + Sub<T, Output = T>,
// {
//     type Output = Vector3<T>;
//     fn sub(self, rhs: Vector3<T>) -> Vector3<T> {
//         Vector3::<T> {
//             x: self.x - rhs.x,
//             y: self.y - rhs.y,
//             z: self.z - rhs.z,
//         }
//     }
// }

// impl<T> Mul<T> for Vector3<T>
// where
//     T: Copy + Mul<T, Output = T>,
// {
//     type Output = Vector3<T>;
//     fn mul(self, rhs: T) -> Vector3<T>
//     where
//         T: Copy + Mul<T, Output = T>,
//     {
//         Vector3::<T> {
//             x: self.x * rhs,
//             y: self.y * rhs,
//             z: self.z * rhs,
//         }
//     }
// }

// impl<T> MulAssign<T> for Vector3<T>
// where
//     T: Copy + MulAssign,
// {
//     fn mul_assign(&mut self, rhs: T) {
//         self.x *= rhs;
//         self.y *= rhs;
//         self.z *= rhs;
//     }
// }

// work around bug
// https://github.com/rust-lang/rust/issues/40395
// impl Div<Float> for Vector3<f32> {
//     type Output = Vector3<f32>;
//     fn div(self, rhs: Float) -> Vector3<f32> {
//         assert_ne!(rhs, 0.0 as Float);
//         let inv: Float = 1.0 as Float / rhs;
//         Vector3::<f32> {
//             x: self.x * inv,
//             y: self.y * inv,
//             z: self.z * inv,
//         }
//     }
// }

// work around bug
// https://github.com/rust-lang/rust/issues/40395
// impl DivAssign<Float> for Vector3<f32> {
//     fn div_assign(&mut self, rhs: Float) {
//         assert_ne!(rhs, 0.0 as Float);
//         let inv: Float = 1.0 as Float / rhs;
//         self.x *= inv;
//         self.y *= inv;
//         self.z *= inv;
//     }
// }

impl_op!(-|a: Vector2f| -> Vector2f { Vector2f { x: -a.x, y: -a.y } });

impl_op!(-|a: Vector3f| -> Vector3f {
    Vector3f {
        x: -a.x,
        y: -a.y,
        z: -a.z,
    }
});

impl_op!(-|a: Normal3f| -> Normal3f {
    Normal3f {
        x: -a.x,
        y: -a.y,
        z: -a.z,
    }
});

// impl<T> Neg for Vector3<T>
// where
//     T: Copy + Neg<Output = T>,
// {
//     type Output = Vector3<T>;
//     fn neg(self) -> Vector3<T> {
//         Vector3::<T> {
//             x: -self.x,
//             y: -self.y,
//             z: -self.z,
//         }
//     }
// }

impl Index<XYZEnum> for Vector3f {
    type Output = Float;
    fn index(&self, index: XYZEnum) -> &Float {
        match index {
            XYZEnum::X => &self.x,
            XYZEnum::Y => &self.y,
            _ => &self.z,
        }
    }
}

impl Index<XYZEnum> for Vector3i {
    type Output = i32;
    fn index(&self, index: XYZEnum) -> &i32 {
        match index {
            XYZEnum::X => &self.x,
            XYZEnum::Y => &self.y,
            _ => &self.z,
        }
    }
}

// impl<T> IndexMut<XYZEnum> for Vector3<T> {
//     fn index_mut(&mut self, index: XYZEnum) -> &mut T {
//         match index {
//             XYZEnum::X => &mut self.x,
//             XYZEnum::Y => &mut self.y,
//             _ => &mut self.z,
//         }
//     }
// }

impl From<Point3f> for Vector3f {
    fn from(p: Point3f) -> Self {
        Vector3f {
            x: p.x,
            y: p.y,
            z: p.z,
        }
    }
}

impl From<Normal3f> for Vector3f {
    fn from(n: Normal3f) -> Self {
        Vector3f {
            x: n.x,
            y: n.y,
            z: n.z,
        }
    }
}

/// Product of the Euclidean magnitudes of the two vectors and the
/// cosine of the angle between them. A return value of zero means
/// both vectors are orthogonal, a value if one means they are
/// codirectional.
pub fn vec3_dot_vec3f(v1: &Vector3f, v2: &Vector3f) -> Float {
    v1.x * v2.x + v1.y * v2.y + v1.z * v2.z
}

/// Product of the Euclidean magnitudes of the two vectors and the
/// cosine of the angle between them. A return value of zero means
/// both vectors are orthogonal, a value if one means they are
/// codirectional.
pub fn vec3_dot_vec3i(v1: &Vector3i, v2: &Vector3i) -> i32 {
    v1.x * v2.x + v1.y * v2.y + v1.z * v2.z
}

/// Product of the Euclidean magnitudes of a vector (and a normal) and
/// the cosine of the angle between them. A return value of zero means
/// both are orthogonal, a value if one means they are codirectional.
pub fn vec3_dot_nrmf(v1: &Vector3f, n2: &Normal3f) -> Float {
    // DCHECK(!v1.HasNaNs() && !n2.HasNaNs());
    v1.x * n2.x + v1.y * n2.y + v1.z * n2.z
}

/// Product of the Euclidean magnitudes of a vector (and a normal) and
/// the cosine of the angle between them. A return value of zero means
/// both are orthogonal, a value if one means they are codirectional.
pub fn vec3_dot_nrmi(v1: &Vector3i, n2: &Normal3i) -> i32 {
    // DCHECK(!v1.HasNaNs() && !n2.HasNaNs());
    v1.x * n2.x + v1.y * n2.y + v1.z * n2.z
}

/// Computes the absolute value of the dot product.
pub fn vec3_abs_dot_vec3f(v1: &Vector3f, v2: &Vector3f) -> Float {
    vec3_dot_vec3f(v1, v2).abs()
}

/// Computes the absolute value of the dot product.
pub fn vec3_abs_dot_vec3i(v1: &Vector3i, v2: &Vector3i) -> i32 {
    vec3_dot_vec3i(v1, v2).abs()
}

/// Computes the absolute value of the dot product.
pub fn vec3_abs_dot_nrmf(v1: &Vector3f, n2: &Normal3f) -> Float {
    vec3_dot_nrmf(v1, n2).abs()
}

/// Computes the absolute value of the dot product.
pub fn vec3_abs_dot_nrmi(v1: &Vector3i, n2: &Normal3i) -> i32 {
    vec3_dot_nrmi(v1, n2).abs()
}

/// Given two vectors in 3D, the cross product is a vector that is
/// perpendicular to both of them.
pub fn vec3_cross_vec3(v1: &Vector3f, v2: &Vector3f) -> Vector3f {
    let v1x: f64 = v1.x as f64;
    let v1y: f64 = v1.y as f64;
    let v1z: f64 = v1.z as f64;
    let v2x: f64 = v2.x as f64;
    let v2y: f64 = v2.y as f64;
    let v2z: f64 = v2.z as f64;
    Vector3f {
        x: ((v1y * v2z) - (v1z * v2y)) as Float,
        y: ((v1z * v2x) - (v1x * v2z)) as Float,
        z: ((v1x * v2y) - (v1y * v2x)) as Float,
    }
}

/// Given a vectors and a normal in 3D, the cross product is a vector
/// that is perpendicular to both of them.
pub fn vec3_cross_nrm(v1: &Vector3f, v2: &Normal3f) -> Vector3f {
    let v1x: f64 = v1.x as f64;
    let v1y: f64 = v1.y as f64;
    let v1z: f64 = v1.z as f64;
    let v2x: f64 = v2.x as f64;
    let v2y: f64 = v2.y as f64;
    let v2z: f64 = v2.z as f64;
    Vector3f {
        x: ((v1y * v2z) - (v1z * v2y)) as Float,
        y: ((v1z * v2x) - (v1x * v2z)) as Float,
        z: ((v1x * v2y) - (v1y * v2x)) as Float,
    }
}

/// Return the largest coordinate value.
pub fn vec3_max_componentf(v: &Vector3f) -> Float {
    v.x.max(v.y.max(v.z))
}

/// Return the largest coordinate value.
pub fn vec3_max_componenti(v: &Vector3i) -> i32 {
    v.x.max(v.y.max(v.z))
}

/// Return the index of the component with the largest value.
pub fn vec3_max_dimensionf(v: &Vector3f) -> usize {
    if v.x > v.y {
        if v.x > v.z {
            0_usize
        } else {
            2_usize
        }
    } else if v.y > v.z {
        1_usize
    } else {
        2_usize
    }
}

/// Return the index of the component with the largest value.
pub fn vec3_max_dimensioni(v: &Vector3i) -> usize {
    if v.x > v.y {
        if v.x > v.z {
            0_usize
        } else {
            2_usize
        }
    } else if v.y > v.z {
        1_usize
    } else {
        2_usize
    }
}

/// Permute the coordinate values according to the povided
/// permutation.
pub fn vec3_permutef(v: &Vector3f, x: usize, y: usize, z: usize) -> Vector3f {
    let v3: [Float; 3] = [v.x, v.y, v.z];
    let xp: Float = v3[x];
    let yp: Float = v3[y];
    let zp: Float = v3[z];
    Vector3f {
        x: xp,
        y: yp,
        z: zp,
    }
}

/// Permute the coordinate values according to the povided
/// permutation.
pub fn vec3_permutei(v: &Vector3i, x: usize, y: usize, z: usize) -> Vector3i {
    let v3: [i32; 3] = [v.x, v.y, v.z];
    let xp: i32 = v3[x];
    let yp: i32 = v3[y];
    let zp: i32 = v3[z];
    Vector3i {
        x: xp,
        y: yp,
        z: zp,
    }
}

/// Construct a local coordinate system given only a single 3D vector.
pub fn vec3_coordinate_system(v1: &Vector3f, v2: &mut Vector3f, v3: &mut Vector3f) {
    if v1.x.abs() > v1.y.abs() {
        *v2 = Vector3f {
            x: -v1.z,
            y: 0.0 as Float,
            z: v1.x,
        } / (v1.x * v1.x + v1.z * v1.z).sqrt();
    } else {
        *v2 = Vector3f {
            x: 0.0 as Float,
            y: v1.z,
            z: -v1.y,
        } / (v1.y * v1.y + v1.z * v1.z).sqrt();
    }
    *v3 = vec3_cross_vec3(v1, &*v2);
}

#[derive(Debug, Default, Copy, Clone)]
pub struct Point2f {
    pub x: Float,
    pub y: Float,
}

impl Point2f {
    pub fn has_nans(&self) -> bool {
        self.x.is_nan() || self.y.is_nan()
    }
}

#[derive(Debug, Default, Copy, Clone)]
pub struct Point2i {
    pub x: i32,
    pub y: i32,
}

// impl<T> PartialEq for Point2<T>
// where
//     T: std::cmp::PartialOrd,
// {
//     fn eq(&self, rhs: &Point2<T>) -> bool {
//         self.x == rhs.x && self.y == rhs.y
//     }
// }

// impl<T> Add<Point2<T>> for Point2<T>
// where
//     T: Add<T, Output = T>,
// {
//     type Output = Point2<T>;
//     fn add(self, rhs: Point2<T>) -> Point2<T> {
//         Point2::<T> {
//             x: self.x + rhs.x,
//             y: self.y + rhs.y,
//         }
//     }
// }

// impl<T> Add<Vector2<T>> for Point2<T>
// where
//     T: Add<T, Output = T>,
// {
//     type Output = Point2<T>;
//     fn add(self, rhs: Vector2<T>) -> Point2<T> {
//         Point2::<T> {
//             x: self.x + rhs.x,
//             y: self.y + rhs.y,
//         }
//     }
// }

// impl<T> Sub<Point2<T>> for Point2<T>
// where
//     T: Sub<T, Output = T>,
// {
//     type Output = Vector2<T>;
//     fn sub(self, rhs: Point2<T>) -> Vector2<T> {
//         Vector2::<T> {
//             x: self.x - rhs.x,
//             y: self.y - rhs.y,
//         }
//     }
// }

// impl<T> Sub<Vector2<T>> for Point2<T>
// where
//     T: Sub<T, Output = T>,
// {
//     type Output = Point2<T>;
//     fn sub(self, rhs: Vector2<T>) -> Point2<T> {
//         Point2::<T> {
//             x: self.x - rhs.x,
//             y: self.y - rhs.y,
//         }
//     }
// }

// impl<T> Mul<T> for Point2<T>
// where
//     T: Copy + Mul<T, Output = T>,
// {
//     type Output = Point2<T>;
//     fn mul(self, rhs: T) -> Point2<T>
//     where
//         T: Copy + Mul<T, Output = T>,
//     {
//         Point2::<T> {
//             x: self.x * rhs,
//             y: self.y * rhs,
//         }
//     }
// }

// impl<T> Neg for Vector2<T>
// where
//     T: Copy + Neg<Output = T>,
// {
//     type Output = Vector2<T>;
//     fn neg(self) -> Vector2<T> {
//         Vector2::<T> {
//             x: -self.x,
//             y: -self.y,
//         }
//     }
// }

impl Index<XYEnum> for Point2f {
    type Output = Float;
    fn index(&self, index: XYEnum) -> &Float {
        match index {
            XYEnum::X => &self.x,
            _ => &self.y,
        }
    }
}

impl Index<XYEnum> for Point2i {
    type Output = i32;
    fn index(&self, index: XYEnum) -> &i32 {
        match index {
            XYEnum::X => &self.x,
            _ => &self.y,
        }
    }
}

impl IndexMut<XYEnum> for Point2f {
    fn index_mut(&mut self, index: XYEnum) -> &mut Float {
        match index {
            XYEnum::X => &mut self.x,
            _ => &mut self.y,
        }
    }
}

impl IndexMut<XYEnum> for Point2i {
    fn index_mut(&mut self, index: XYEnum) -> &mut i32 {
        match index {
            XYEnum::X => &mut self.x,
            _ => &mut self.y,
        }
    }
}

/// Apply floor operation component-wise.
pub fn pnt2_floor(p: Point2f) -> Point2f {
    Point2f {
        x: p.x.floor(),
        y: p.y.floor(),
    }
}

/// Apply ceil operation component-wise.
pub fn pnt2_ceil(p: Point2f) -> Point2f {
    Point2f {
        x: p.x.ceil(),
        y: p.y.ceil(),
    }
}

/// Apply std::cmp::min operation component-wise.
pub fn pnt2_min_pnt2i(pa: Point2i, pb: Point2i) -> Point2i {
    Point2i {
        x: std::cmp::min(pa.x, pb.x),
        y: std::cmp::min(pa.y, pb.y),
    }
}

/// Apply std::cmp::max operation component-wise.
pub fn pnt2_max_pnt2i(pa: Point2i, pb: Point2i) -> Point2i {
    Point2i {
        x: std::cmp::max(pa.x, pb.x),
        y: std::cmp::max(pa.y, pb.y),
    }
}

/// Given a bounding box and a point, the **bnd2_union_pnt2()**
/// function returns a new bounding box that encompasses that point as
/// well as the original box.
pub fn bnd2_union_pnt2(b: &Bounds2f, p: Point2f) -> Bounds2f {
    let p_min: Point2f = Point2f {
        x: b.p_min.x.min(p.x),
        y: b.p_min.y.min(p.y),
    };
    let p_max: Point2f = Point2f {
        x: b.p_max.x.max(p.x),
        y: b.p_max.y.max(p.y),
    };
    Bounds2f { p_min, p_max }
}

/// Determine if a given point is inside the bounding box.
pub fn pnt2_inside_bnd2f(pt: Point2f, b: &Bounds2f) -> bool {
    pt.x >= b.p_min.x && pt.x <= b.p_max.x && pt.y >= b.p_min.y && pt.y <= b.p_max.y
}

/// Determine if a given point is inside the bounding box.
pub fn pnt2_inside_bnd2i(pt: Point2i, b: &Bounds2i) -> bool {
    pt.x >= b.p_min.x && pt.x <= b.p_max.x && pt.y >= b.p_min.y && pt.y <= b.p_max.y
}

/// Is a 2D point inside a 2D bound?
pub fn pnt2_inside_exclusivef(pt: Point2f, b: &Bounds2f) -> bool {
    pt.x >= b.p_min.x && pt.x < b.p_max.x && pt.y >= b.p_min.y && pt.y < b.p_max.y
}

/// Is a 2D point inside a 2D bound?
pub fn pnt2_inside_exclusivei(pt: Point2i, b: &Bounds2i) -> bool {
    pt.x >= b.p_min.x && pt.x < b.p_max.x && pt.y >= b.p_min.y && pt.y < b.p_max.y
}

/// Pads the bounding box by a constant factor in both dimensions.
pub fn bnd2_expand(b: &Bounds2f, delta: Float) -> Bounds2f {
    Bounds2f {
        p_min: b.p_min - Vector2f { x: delta, y: delta },
        p_max: b.p_max + Vector2f { x: delta, y: delta },
    }
}

#[derive(Debug, Default, Copy, Clone)]
pub struct Point3f {
    pub x: Float,
    pub y: Float,
    pub z: Float,
}

#[derive(Debug, Default, Copy, Clone)]
pub struct Point3i {
    pub x: i32,
    pub y: i32,
    pub z: i32,
}

impl Point3f {
    pub fn has_nans(&self) -> bool {
        self.x.is_nan() || self.y.is_nan() || self.z.is_nan()
    }
}

// impl<T> AddAssign<Point3<T>> for Point3<T>
// where
//     T: AddAssign,
// {
//     fn add_assign(&mut self, rhs: Point3<T>) {
//         self.x += rhs.x;
//         self.y += rhs.y;
//         self.z += rhs.z;
//     }
// }

impl_op_ex!(+|a: &Point3f, b: &Point3f| -> Point3f {
    Point3f {
        x: a.x + b.x,
        y: a.y + b.y,
        z: a.z + b.z,
    }
});

impl_op_ex!(+|a: &Point2i, b: &Point2i| -> Point2i {
    Point2i {
        x: a.x + b.x,
        y: a.y + b.y,
    }
});

impl_op_ex!(+|a: &Point2f, b: &Point2f| -> Point2f {
    Point2f {
        x: a.x + b.x,
        y: a.y + b.y,
    }
});

impl_op_ex!(+|a: &Vector3f, b: &Vector3f| -> Vector3f {
    Vector3f {
        x: a.x + b.x,
        y: a.y + b.y,
        z: a.z + b.z,
    }
});

impl_op_ex!(+|a: &Normal3f, b: &Normal3f| -> Normal3f {
    Normal3f {
        x: a.x + b.x,
        y: a.y + b.y,
        z: a.z + b.z,
    }
});

impl_op_ex!(-|a: &Normal3f, b: &Normal3f| -> Normal3f {
    Normal3f {
        x: a.x - b.x,
        y: a.y - b.y,
        z: a.z - b.z,
    }
});

impl_op_ex!(-|a: &Vector3f, b: &Vector3f| -> Vector3f {
    Vector3f {
        x: a.x - b.x,
        y: a.y - b.y,
        z: a.z - b.z,
    }
});

impl_op_ex!(+|a: &Point3f, b: &Vector3f| -> Point3f {
    Point3f {
        x: a.x + b.x,
        y: a.y + b.y,
        z: a.z + b.z,
    }
});

impl_op_ex!(+|a: &Point3i, b: &Vector3i| -> Point3i {
    Point3i {
        x: a.x + b.x,
        y: a.y + b.y,
        z: a.z + b.z,
    }
});

impl_op_ex!(+|a: &Point2f, b: &Vector2f| -> Point2f {
    Point2f {
        x: a.x + b.x,
        y: a.y + b.y,
    }
});

impl_op_ex!(-|a: &Point3f, b: &Point3f| -> Vector3f {
    Vector3f {
        x: a.x - b.x,
        y: a.y - b.y,
        z: a.z - b.z,
    }
});

impl_op_ex!(-|a: &Point3i, b: &Point3i| -> Vector3i {
    Vector3i {
        x: a.x - b.x,
        y: a.y - b.y,
        z: a.z - b.z,
    }
});

impl_op_ex!(-|a: &Point2f, b: &Point2f| -> Vector2f {
    Vector2f {
        x: a.x - b.x,
        y: a.y - b.y,
    }
});

impl_op_ex!(-|a: &Point2i, b: &Point2i| -> Vector2i {
    Vector2i {
        x: a.x - b.x,
        y: a.y - b.y,
    }
});

impl_op_ex!(-|a: &Point2f, b: &Vector2f| -> Point2f {
    Point2f {
        x: a.x - b.x,
        y: a.y - b.y,
    }
});

impl_op_ex!(-|a: &Point2i, b: &Vector2i| -> Point2i {
    Point2i {
        x: a.x - b.x,
        y: a.y - b.y,
    }
});

impl_op_ex!(-|a: &Point3f, b: &Vector3f| -> Point3f {
    Point3f {
        x: a.x - b.x,
        y: a.y - b.y,
        z: a.z - b.z,
    }
});

// impl<T> Add<Point3<T>> for Point3<T>
// where
//     T: Add<T, Output = T>,
// {
//     type Output = Point3<T>;
//     fn add(self, rhs: Point3<T>) -> Point3<T> {
//         Point3::<T> {
//             x: self.x + rhs.x,
//             y: self.y + rhs.y,
//             z: self.z + rhs.z,
//         }
//     }
// }

// impl<T> Add<Vector3<T>> for Point3<T>
// where
//     T: Add<T, Output = T>,
// {
//     type Output = Point3<T>;
//     fn add(self, rhs: Vector3<T>) -> Point3<T> {
//         Point3::<T> {
//             x: self.x + rhs.x,
//             y: self.y + rhs.y,
//             z: self.z + rhs.z,
//         }
//     }
// }

// impl<T> AddAssign<Vector3<T>> for Point3<T>
// where
//     T: AddAssign,
// {
//     fn add_assign(&mut self, rhs: Vector3<T>) {
//         self.x += rhs.x;
//         self.y += rhs.y;
//         self.z += rhs.z;
//     }
// }

// impl<T> Sub<Vector3<T>> for Point3<T>
// where
//     T: Sub<T, Output = T>,
// {
//     type Output = Point3<T>;
//     fn sub(self, rhs: Vector3<T>) -> Point3<T> {
//         Point3::<T> {
//             x: self.x - rhs.x,
//             y: self.y - rhs.y,
//             z: self.z - rhs.z,
//         }
//     }
// }

impl_op_ex!(*|a: &Point2f, b: Float| -> Point2f {
    Point2f {
        x: a.x * b,
        y: a.y * b,
    }
});

impl_op_ex!(*|a: &Point3f, b: Float| -> Point3f {
    Point3f {
        x: a.x * b,
        y: a.y * b,
        z: a.z * b,
    }
});

impl_op_ex!(*|a: &Normal3f, b: Float| -> Normal3f {
    Normal3f {
        x: a.x * b,
        y: a.y * b,
        z: a.z * b,
    }
});

impl_op_ex!(*|a: &Vector3f, b: Float| -> Vector3f {
    Vector3f {
        x: a.x * b,
        y: a.y * b,
        z: a.z * b,
    }
});

impl_op_ex!(/|a: &Point3f, b: Float| -> Point3f {
    assert_ne!(b, 0.0 as Float);
    let inv: Float = 1.0 as Float / b;
    Point3f {
        x: a.x * inv,
        y: a.y * inv,
        z: a.z * inv,
    }
});

impl_op_ex!(/|a: &Vector3f, b: Float| -> Vector3f {
    assert_ne!(b, 0.0 as Float);
    let inv: Float = 1.0 as Float / b;
    Vector3f {
        x: a.x * inv,
        y: a.y * inv,
        z: a.z * inv,
    }
});

impl_op_ex!(/|a: &Vector2f, b: Float| -> Vector2f {
    assert_ne!(b, 0.0 as Float);
    let inv: Float = 1.0 as Float / b;
    Vector2f {
        x: a.x * inv,
        y: a.y * inv,
    }
});

impl_op_ex!(/|a: &Normal3f, b: Float| -> Normal3f {
    assert_ne!(b, 0.0 as Float);
    let inv: Float = 1.0 as Float / b;
    Normal3f {
        x: a.x * inv,
        y: a.y * inv,
        z: a.z * inv,
    }
});

// impl<T> Mul<T> for Point3<T>
// where
//     T: Copy + Mul<T, Output = T>,
// {
//     type Output = Point3<T>;
//     fn mul(self, rhs: T) -> Point3<T>
//     where
//         T: Copy + Mul<T, Output = T>,
//     {
//         Point3::<T> {
//             x: self.x * rhs,
//             y: self.y * rhs,
//             z: self.z * rhs,
//         }
//     }
// }

// impl<T> MulAssign<T> for Point3<T>
// where
//     T: Copy + MulAssign,
// {
//     fn mul_assign(&mut self, rhs: T) {
//         self.x *= rhs;
//         self.y *= rhs;
//         self.z *= rhs;
//     }
// }

// work around bug
// https://github.com/rust-lang/rust/issues/40395
// impl Div<Float> for Point3<f32> {
//     type Output = Point3<f32>;
//     fn div(self, rhs: Float) -> Point3<f32> {
//         assert_ne!(rhs, 0.0 as Float);
//         let inv: Float = 1.0 as Float / rhs;
//         Point3::<f32> {
//             x: self.x * inv,
//             y: self.y * inv,
//             z: self.z * inv,
//         }
//     }
// }

impl_op!(*= |a: &mut Point2f, b: Float| {
    a.x *= b;
    a.y *= b;
});

impl_op!(*= |a: &mut Point3f, b: Float| {
    a.x *= b;
    a.y *= b;
    a.z *= b;
});

impl_op!(+= |a: &mut Point3f, b: Point3f| {
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
});

impl_op!(+= |a: &mut Point3f, b: Vector3f| {
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
});

impl_op!(+= |a: &mut Vector3f, b: Vector3f| {
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
});

impl_op!(*= |a: &mut Vector2f, b: Float| {
    a.x *= b;
    a.y *= b;
});

impl_op!(*= |a: &mut Vector3f, b: Float| {
    a.x *= b;
    a.y *= b;
    a.z *= b;
});

impl_op!(*= |a: &mut Normal3f, b: Float| {
    a.x *= b;
    a.y *= b;
    a.z *= b;
});

impl_op!(/= |a: &mut Point3f, b: Float| {
    assert_ne!(b, 0.0 as Float);
    let inv: Float = 1.0 as Float / b;
    a.x *= inv;
    a.y *= inv;
    a.z *= inv;
});

impl_op!(/= |a: &mut Vector3f, b: Float| {
    assert_ne!(b, 0.0 as Float);
    let inv: Float = 1.0 as Float / b;
    a.x *= inv;
    a.y *= inv;
    a.z *= inv;
});

// work around bug
// https://github.com/rust-lang/rust/issues/40395
// impl DivAssign<Float> for Point3<f32> {
//     fn div_assign(&mut self, rhs: Float) {
//         assert_ne!(rhs, 0.0 as Float);
//         let inv: Float = 1.0 as Float / rhs;
//         self.x *= inv;
//         self.y *= inv;
//         self.z *= inv;
//     }
// }

impl Index<XYZEnum> for Point3f {
    type Output = Float;
    fn index(&self, index: XYZEnum) -> &Float {
        match index {
            XYZEnum::X => &self.x,
            XYZEnum::Y => &self.y,
            _ => &self.z,
        }
    }
}

impl Index<XYZEnum> for Point3i {
    type Output = i32;
    fn index(&self, index: XYZEnum) -> &i32 {
        match index {
            XYZEnum::X => &self.x,
            XYZEnum::Y => &self.y,
            _ => &self.z,
        }
    }
}

impl IndexMut<XYZEnum> for Point3f {
    fn index_mut(&mut self, index: XYZEnum) -> &mut Float {
        match index {
            XYZEnum::X => &mut self.x,
            XYZEnum::Y => &mut self.y,
            _ => &mut self.z,
        }
    }
}

impl IndexMut<XYZEnum> for Point3i {
    fn index_mut(&mut self, index: XYZEnum) -> &mut i32 {
        match index {
            XYZEnum::X => &mut self.x,
            XYZEnum::Y => &mut self.y,
            _ => &mut self.z,
        }
    }
}

/// Permute the coordinate values according to the povided
/// permutation.
pub fn pnt3_permutef(v: &Point3f, x: usize, y: usize, z: usize) -> Point3f {
    let v3: [Float; 3] = [v.x, v.y, v.z];
    let xp: Float = v3[x];
    let yp: Float = v3[y];
    let zp: Float = v3[z];
    Point3f {
        x: xp,
        y: yp,
        z: zp,
    }
}

/// Permute the coordinate values according to the povided
/// permutation.
pub fn pnt3_permutei(v: &Point3i, x: usize, y: usize, z: usize) -> Point3i {
    let v3: [i32; 3] = [v.x, v.y, v.z];
    let xp: i32 = v3[x];
    let yp: i32 = v3[y];
    let zp: i32 = v3[z];
    Point3i {
        x: xp,
        y: yp,
        z: zp,
    }
}

/// Interpolate linearly between two provided points.
pub fn pnt3_lerp(t: Float, p0: &Point3f, p1: &Point3f) -> Point3f {
    *p0 * (1.0 as Float - t) as Float + *p1 * t
}

/// Apply floor operation component-wise.
pub fn pnt3_floor(p: &Point3f) -> Point3f {
    Point3f {
        x: p.x.floor(),
        y: p.y.floor(),
        z: p.z.floor(),
    }
}

/// Apply ceil operation component-wise.
pub fn pnt3_ceil(p: &Point3f) -> Point3f {
    Point3f {
        x: p.x.ceil(),
        y: p.y.ceil(),
        z: p.z.ceil(),
    }
}

/// Apply abs operation component-wise.
pub fn pnt3_abs(p: &Point3f) -> Point3f {
    Point3f {
        x: p.x.abs(),
        y: p.y.abs(),
        z: p.z.abs(),
    }
}

/// The distance between two points is the length of the vector
/// between them.
pub fn pnt3_distancef(p1: &Point3f, p2: &Point3f) -> Float {
    (p1 - p2).length()
}

/// The distance squared between two points is the length of the
/// vector between them squared.
pub fn pnt3_distance_squaredf(p1: &Point3f, p2: &Point3f) -> Float {
    (p1 - p2).length_squared()
}

/// When tracing spawned rays leaving the intersection point p, we
/// offset their origins enough to ensure that they are past the
/// boundary of the error box and thus won't incorrectly re-intersect
/// the surface.
pub fn pnt3_offset_ray_origin(
    p: &Point3f,
    p_error: &Vector3f,
    n: &Normal3f,
    w: &Vector3f,
) -> Point3f {
    let d: Float = nrm_dot_vec3f(&nrm_absf(n), p_error);
    let mut offset: Vector3f = Vector3f::from(*n) * d;
    if vec3_dot_nrmf(w, n) < 0.0 as Float {
        offset = -offset;
    }
    let mut po: Point3f = *p + offset;
    // round offset point _po_ away from _p_
    for i in XYZEnum::iter() {
        if offset[i] > 0.0 as Float {
            po[i] = next_float_up(po[i]);
        } else if offset[i] < 0.0 as Float {
            po[i] = next_float_down(po[i]);
        }
    }
    po
}

/// Calculate appropriate direction vector from two angles.
pub fn spherical_direction(sin_theta: Float, cos_theta: Float, phi: Float) -> Vector3f {
    Vector3f {
        x: sin_theta * phi.cos(),
        y: sin_theta * phi.sin(),
        z: cos_theta,
    }
}

/// Take three basis vectors representing the x, y, and z axes and
/// return the appropriate direction vector with respect to the
/// coordinate frame defined by them.
pub fn spherical_direction_vec3(
    sin_theta: Float,
    cos_theta: Float,
    phi: Float,
    x: &Vector3f,
    y: &Vector3f,
    z: &Vector3f,
) -> Vector3f {
    *x * (sin_theta * phi.cos()) + *y * (sin_theta * phi.sin()) + *z * cos_theta
}

/// Conversion of a direction to spherical angles. Note that
/// **spherical_theta()** assumes that the vector **v** has been
/// normalized before being passed in.
pub fn spherical_theta(v: &Vector3f) -> Float {
    clamp_t(v.z, -1.0 as Float, 1.0 as Float).acos()
}

/// Conversion of a direction to spherical angles.
pub fn spherical_phi(v: &Vector3f) -> Float {
    let p: Float = v.y.atan2(v.x);
    if p < 0.0 as Float {
        p + 2.0 as Float * PI
    } else {
        p
    }
}

#[derive(Debug, Default, Copy, Clone, PartialEq)]
pub struct Normal3f {
    pub x: Float,
    pub y: Float,
    pub z: Float,
}

impl Normal3f {
    pub fn length_squared(&self) -> Float {
        self.x * self.x + self.y * self.y + self.z * self.z
    }
    pub fn length(&self) -> Float {
        self.length_squared().sqrt()
    }
    /// Compute a new normal pointing in the same direction but with unit
    /// length.
    pub fn normalize(&self) -> Normal3f {
        *self / self.length()
    }
}

#[derive(Debug, Default, Copy, Clone)]
pub struct Normal3i {
    pub x: i32,
    pub y: i32,
    pub z: i32,
}

// impl<T> Add for Normal3<T>
// where
//     T: Copy + Add<T, Output = T>,
// {
//     type Output = Normal3<T>;
//     fn add(self, rhs: Normal3<T>) -> Normal3<T> {
//         Normal3::<T> {
//             x: self.x + rhs.x,
//             y: self.y + rhs.y,
//             z: self.z + rhs.z,
//         }
//     }
// }

// impl<T> Sub for Normal3<T>
// where
//     T: Copy + Sub<T, Output = T>,
// {
//     type Output = Normal3<T>;
//     fn sub(self, rhs: Normal3<T>) -> Normal3<T> {
//         Normal3::<T> {
//             x: self.x - rhs.x,
//             y: self.y - rhs.y,
//             z: self.z - rhs.z,
//         }
//     }
// }

// impl<T> Mul<T> for Normal3<T>
// where
//     T: Copy + Mul<T, Output = T>,
// {
//     type Output = Normal3<T>;
//     fn mul(self, rhs: T) -> Normal3<T>
//     where
//         T: Copy + Mul<T, Output = T>,
//     {
//         Normal3::<T> {
//             x: self.x * rhs,
//             y: self.y * rhs,
//             z: self.z * rhs,
//         }
//     }
// }

// impl<T> MulAssign<T> for Normal3<T>
// where
//     T: Copy + MulAssign,
// {
//     fn mul_assign(&mut self, rhs: T) {
//         self.x *= rhs;
//         self.y *= rhs;
//         self.z *= rhs;
//     }
// }

// impl<T> Neg for Normal3<T>
// where
//     T: Copy + Neg<Output = T>,
// {
//     type Output = Normal3<T>;
//     fn neg(self) -> Normal3<T> {
//         Normal3::<T> {
//             x: -self.x,
//             y: -self.y,
//             z: -self.z,
//         }
//     }
// }

impl Index<XYZEnum> for Normal3f {
    type Output = Float;
    fn index(&self, index: XYZEnum) -> &Float {
        match index {
            XYZEnum::X => &self.x,
            XYZEnum::Y => &self.y,
            _ => &self.z,
        }
    }
}

// impl<T> IndexMut<XYZEnum> for Normal3<T> {
//     fn index_mut(&mut self, index: XYZEnum) -> &mut T {
//         match index {
//             XYZEnum::X => &mut self.x,
//             XYZEnum::Y => &mut self.y,
//             _ => &mut self.z,
//         }
//     }
// }

// impl<T> Normal3<T> {
//     pub fn length_squared(&self) -> T
//     where
//         T: Copy + Add<T, Output = T> + Mul<T, Output = T>,
//     {
//         self.x * self.x + self.y * self.y + self.z * self.z
//     }
//     pub fn length(&self) -> T
//     where
//         T: num::Float,
//     {
//         self.length_squared().sqrt()
//     }
// }

// impl<T> PartialEq for Normal3<T>
// where
//     T: std::cmp::PartialOrd,
// {
//     fn eq(&self, rhs: &Normal3<T>) -> bool {
//         self.x == rhs.x && self.y == rhs.y && self.z == rhs.z
//     }
// }

// work around bug
// https://github.com/rust-lang/rust/issues/40395
// impl Div<Float> for Normal3<f32> {
//     type Output = Normal3<f32>;
//     fn div(self, rhs: Float) -> Normal3<f32> {
//         assert_ne!(rhs, 0.0 as Float);
//         let inv: Float = 1.0 as Float / rhs;
//         Normal3::<f32> {
//             x: self.x * inv,
//             y: self.y * inv,
//             z: self.z * inv,
//         }
//     }
// }

impl From<Vector3f> for Normal3f {
    fn from(v: Vector3f) -> Self {
        Normal3f {
            x: v.x,
            y: v.y,
            z: v.z,
        }
    }
}

/// Given a normal and a vector in 3D, the cross product is a vector
/// that is perpendicular to both of them.
pub fn nrm_cross_vec3(n1: &Normal3f, v2: &Vector3f) -> Vector3f {
    let n1x: f64 = n1.x as f64;
    let n1y: f64 = n1.y as f64;
    let n1z: f64 = n1.z as f64;
    let v2x: f64 = v2.x as f64;
    let v2y: f64 = v2.y as f64;
    let v2z: f64 = v2.z as f64;
    Vector3f {
        x: ((n1y * v2z) - (n1z * v2y)) as Float,
        y: ((n1z * v2x) - (n1x * v2z)) as Float,
        z: ((n1x * v2y) - (n1y * v2x)) as Float,
    }
}

/// Product of the Euclidean magnitudes of a normal (and another
/// normal) and the cosine of the angle between them. A return value
/// of zero means both are orthogonal, a value if one means they are
/// codirectional.
pub fn nrm_dot_nrmf(n1: &Normal3f, n2: &Normal3f) -> Float {
    n1.x * n2.x + n1.y * n2.y + n1.z * n2.z
}

/// Product of the Euclidean magnitudes of a normal (and another
/// normal) and the cosine of the angle between them. A return value
/// of zero means both are orthogonal, a value if one means they are
/// codirectional.
pub fn nrm_dot_nrmi(n1: &Normal3i, n2: &Normal3i) -> i32 {
    n1.x * n2.x + n1.y * n2.y + n1.z * n2.z
}

/// Product of the Euclidean magnitudes of a normal (and a vector) and
/// the cosine of the angle between them. A return value of zero means
/// both are orthogonal, a value if one means they are codirectional.
pub fn nrm_dot_vec3f(n1: &Normal3f, v2: &Vector3f) -> Float {
    n1.x * v2.x + n1.y * v2.y + n1.z * v2.z
}

/// Product of the Euclidean magnitudes of a normal (and a vector) and
/// the cosine of the angle between them. A return value of zero means
/// both are orthogonal, a value if one means they are codirectional.
pub fn nrm_dot_vec3i(n1: &Normal3i, v2: &Vector3i) -> i32 {
    n1.x * v2.x + n1.y * v2.y + n1.z * v2.z
}

/// Computes the absolute value of the dot product.
pub fn nrm_abs_dot_vec3f(n1: &Normal3f, v2: &Vector3f) -> Float {
    nrm_dot_vec3f(n1, v2).abs()
}

/// Computes the absolute value of the dot product.
pub fn nrm_abs_dot_vec3i(n1: &Normal3i, v2: &Vector3i) -> i32 {
    nrm_dot_vec3i(n1, v2).abs()
}

/// Return normal with the absolute value of each coordinate.
pub fn nrm_absf(n: &Normal3f) -> Normal3f {
    Normal3f {
        x: n.x.abs(),
        y: n.y.abs(),
        z: n.z.abs(),
    }
}

/// Return normal with the absolute value of each coordinate.
pub fn nrm_absi(n: &Normal3i) -> Normal3i {
    Normal3i {
        x: n.x.abs(),
        y: n.y.abs(),
        z: n.z.abs(),
    }
}

/// Flip a surface normal so that it lies in the same hemisphere as a
/// given vector.
pub fn nrm_faceforward_vec3(n: &Normal3f, v: &Vector3f) -> Normal3f {
    if nrm_dot_vec3f(n, v) < 0.0 as Float {
        -(*n)
    } else {
        *n
    }
}

/// Flip a surface normal so that it lies in the same hemisphere as a
/// given normal.
pub fn nrm_faceforward_nrm(n: &Normal3f, n2: &Normal3f) -> Normal3f {
    if nrm_dot_nrmf(n, n2) < 0.0 as Float {
        -(*n)
    } else {
        *n
    }
}

// pub type Bounds2f = Bounds2<Float>;
// pub type Bounds2i = Bounds2<i32>;
// pub type Bounds3f = Bounds3<Float>;
// pub type Bounds3i = Bounds3<i32>;

#[derive(Debug, Default, Copy, Clone)]
pub struct Bounds2f {
    pub p_min: Point2f,
    pub p_max: Point2f,
}

#[derive(Debug, Default, Copy, Clone)]
pub struct Bounds2i {
    pub p_min: Point2i,
    pub p_max: Point2i,
}

impl Bounds2f {
    pub fn diagonal(&self) -> Vector2f {
        self.p_max - self.p_min
    }
    pub fn area(&self) -> Float {
        let d: Vector2f = self.p_max - self.p_min;
        d.x * d.y
    }
    pub fn lerp(&self, t: Point2f) -> Point2f {
        Point2f {
            x: lerp(t.x, self.p_min.x, self.p_max.x),
            y: lerp(t.y, self.p_min.y, self.p_max.y),
        }
    }
    pub fn offset(&self, p: Point2f) -> Vector2f {
        let mut o: Vector2f = p - self.p_min;
        if self.p_max.x > self.p_min.x {
            o.x /= self.p_max.x - self.p_min.x;
        }
        if self.p_max.y > self.p_min.y {
            o.y /= self.p_max.y - self.p_min.y;
        }
        o
    }
}

impl Bounds2i {
    pub fn new(p1: Point2i, p2: Point2i) -> Self {
        let p_min: Point2i = Point2i {
            x: std::cmp::min(p1.x, p2.x),
            y: std::cmp::min(p1.y, p2.y),
        };
        let p_max: Point2i = Point2i {
            x: std::cmp::max(p1.x, p2.x),
            y: std::cmp::max(p1.y, p2.y),
        };
        Bounds2i { p_min, p_max }
    }
    pub fn diagonal(&self) -> Vector2i {
        self.p_max - self.p_min
    }
    pub fn area(&self) -> i32 {
        let d: Vector2i = self.p_max - self.p_min;
        d.x * d.y
    }
}

pub struct Bounds2Iterator<'a> {
    p: Point2i,
    bounds: &'a Bounds2i,
}

impl<'a> Iterator for Bounds2Iterator<'a> {
    type Item = Point2i;

    fn next(&mut self) -> Option<Point2i> {
        self.p.x += 1;
        if self.p.x == self.bounds.p_max.x {
            self.p.x = self.bounds.p_min.x;
            self.p.y += 1;
        }
        if self.p.y == self.bounds.p_max.y {
            None
        } else {
            Some(self.p)
        }
    }
}

impl<'a> IntoIterator for &'a Bounds2i {
    type Item = Point2i;
    type IntoIter = Bounds2Iterator<'a>;

    fn into_iter(self) -> Self::IntoIter {
        Bounds2Iterator {
            // need to start 1 before p_min.x as next() will be called
            // to get the first element
            p: Point2i {
                x: self.p_min.x - 1,
                y: self.p_min.y,
            },
            bounds: self,
        }
    }
}

/// The intersection of two bounding boxes can be found by computing
/// the maximum of their two respective minimum coordinates and the
/// minimum of their maximum coordinates.
pub fn bnd2_intersect_bnd2i(b1: &Bounds2i, b2: &Bounds2i) -> Bounds2i {
    Bounds2i {
        p_min: Point2i {
            x: std::cmp::max(b1.p_min.x, b2.p_min.x),
            y: std::cmp::max(b1.p_min.y, b2.p_min.y),
        },
        p_max: Point2i {
            x: std::cmp::min(b1.p_max.x, b2.p_max.x),
            y: std::cmp::min(b1.p_max.y, b2.p_max.y),
        },
    }
}

#[derive(Debug, Copy, Clone)]
pub struct Bounds3f {
    pub p_min: Point3f,
    pub p_max: Point3f,
}

#[derive(Debug, Copy, Clone)]
pub struct Bounds3i {
    pub p_min: Point3i,
    pub p_max: Point3i,
}

// work around bug
// https://github.com/rust-lang/rust/issues/40395
impl Default for Bounds3f {
    fn default() -> Bounds3f {
        let min_num: Float = std::f32::MIN;
        let max_num: Float = std::f32::MAX;
        // Bounds3f
        Bounds3f {
            p_min: Point3f {
                x: max_num,
                y: max_num,
                z: max_num,
            },
            p_max: Point3f {
                x: min_num,
                y: min_num,
                z: min_num,
            },
        }
    }
}

impl Bounds3f {
    pub fn new(p1: Point3f, p2: Point3f) -> Self {
        let p_min: Point3f = Point3f {
            x: p1.x.min(p2.x),
            y: p1.y.min(p2.y),
            z: p1.z.min(p2.z),
        };
        let p_max: Point3f = Point3f {
            x: p1.x.max(p2.x),
            y: p1.y.max(p2.y),
            z: p1.z.max(p2.z),
        };
        Bounds3f { p_min, p_max }
    }
    pub fn corner(&self, corner: u8) -> Point3f {
        // assert!(corner >= 0_u8);
        assert!(corner < 8_u8);
        let x: Float = if corner & 1 == 0 {
            self.p_min.x
        } else {
            self.p_max.x
        };
        let y: Float = if corner & 2 == 0 {
            self.p_min.y
        } else {
            self.p_max.y
        };
        let z: Float = if corner & 4 == 0 {
            self.p_min.z
        } else {
            self.p_max.z
        };
        Point3f { x, y, z }
    }
    pub fn diagonal(&self) -> Vector3f {
        self.p_max - self.p_min
    }
    pub fn surface_area(&self) -> Float {
        let d: Vector3f = self.diagonal();
        // 2 * (d.x * d.y + d.x * d.z + d.y * d.z)
        let r: Float = d.x * d.y + d.x * d.z + d.y * d.z;
        r + r // avoid '2 *'
    }
    pub fn maximum_extent(&self) -> u8 {
        let d: Vector3f = self.diagonal();
        if d.x > d.y && d.x > d.z {
            0_u8
        } else if d.y > d.z {
            1_u8
        } else {
            2_u8
        }
    }
    pub fn offset(&self, p: &Point3f) -> Vector3f {
        let mut o: Vector3f = p - self.p_min;
        if self.p_max.x > self.p_min.x {
            o.x /= self.p_max.x - self.p_min.x;
        }
        if self.p_max.y > self.p_min.y {
            o.y /= self.p_max.y - self.p_min.y;
        }
        if self.p_max.z > self.p_min.z {
            o.z /= self.p_max.z - self.p_min.z;
        }
        o
    }
    pub fn bounding_sphere(b: &Bounds3f, center: &mut Point3f, radius: &mut Float) {
        let p_min: Point3f = b.p_min as Point3f;
        let p_max: Point3f = b.p_max as Point3f;
        let sum: Point3f = p_min + p_max;
        *center = sum / 2.0;
        let center_copy: Point3f = *center as Point3f;
        let is_inside: bool = pnt3_inside_bnd3(&center_copy, b);
        if is_inside {
            *radius = pnt3_distancef(&center_copy, &p_max);
        } else {
            *radius = 0.0;
        }
    }
}

impl Bounds3i {
    pub fn new(p1: Point3i, p2: Point3i) -> Self {
        let p_min: Point3i = Point3i {
            x: p1.x.min(p2.x),
            y: p1.y.min(p2.y),
            z: p1.z.min(p2.z),
        };
        let p_max: Point3i = Point3i {
            x: p1.x.max(p2.x),
            y: p1.y.max(p2.y),
            z: p1.z.max(p2.z),
        };
        Bounds3i { p_min, p_max }
    }
    pub fn corner(&self, corner: u8) -> Point3i {
        // assert!(corner >= 0_u8);
        assert!(corner < 8_u8);
        let x: i32 = if corner & 1 == 0 {
            self.p_min.x
        } else {
            self.p_max.x
        };
        let y: i32 = if corner & 2 == 0 {
            self.p_min.y
        } else {
            self.p_max.y
        };
        let z: i32 = if corner & 4 == 0 {
            self.p_min.z
        } else {
            self.p_max.z
        };
        Point3i { x, y, z }
    }
    pub fn diagonal(&self) -> Vector3i {
        self.p_max - self.p_min
    }
    pub fn surface_area(&self) -> i32 {
        let d: Vector3i = self.diagonal();
        // 2 * (d.x * d.y + d.x * d.z + d.y * d.z)
        let r: i32 = d.x * d.y + d.x * d.z + d.y * d.z;
        r + r // avoid '2 *'
    }
    pub fn maximum_extent(&self) -> u8 {
        let d: Vector3i = self.diagonal();
        if d.x > d.y && d.x > d.z {
            0_u8
        } else if d.y > d.z {
            1_u8
        } else {
            2_u8
        }
    }
    pub fn offset(&self, p: &Point3i) -> Vector3i {
        let mut o: Vector3i = *p - self.p_min;
        if self.p_max.x > self.p_min.x {
            o.x /= self.p_max.x - self.p_min.x;
        }
        if self.p_max.y > self.p_min.y {
            o.y /= self.p_max.y - self.p_min.y;
        }
        if self.p_max.z > self.p_min.z {
            o.z /= self.p_max.z - self.p_min.z;
        }
        o
    }
    pub fn bounding_sphere(b: &Bounds3f, center: &mut Point3f, radius: &mut Float) {
        let p_min: Point3f = b.p_min as Point3f;
        let p_max: Point3f = b.p_max as Point3f;
        let sum: Point3f = p_min + p_max;
        *center = sum / 2.0;
        let center_copy: Point3f = *center as Point3f;
        let is_inside: bool = pnt3_inside_bnd3(&center_copy, b);
        if is_inside {
            *radius = pnt3_distancef(&center_copy, &p_max);
        } else {
            *radius = 0.0;
        }
    }
}

impl Bounds3f {
    pub fn lerp(&self, t: &Point3f) -> Point3f {
        Point3f {
            x: lerp(t.x, self.p_min.x as Float, self.p_max.x as Float),
            y: lerp(t.y, self.p_min.y as Float, self.p_max.y as Float),
            z: lerp(t.z, self.p_min.z as Float, self.p_max.z as Float),
        }
    }
    pub fn intersect_b(&self, ray: &Ray, hitt0: &mut Float, hitt1: &mut Float) -> bool {
        let mut t0: Float = 0.0;
        let mut t1: Float = ray.t_max.get();
        for i in XYZEnum::iter() {
            // update interval for _i_th bounding box slab
            let inv_ray_dir: Float = 1.0 as Float / ray.d[i];
            let mut t_near: Float = (self.p_min[i] - ray.o[i]) * inv_ray_dir;
            let mut t_far: Float = (self.p_max[i] - ray.o[i]) * inv_ray_dir;
            // update parametric interval from slab intersection $t$ values
            if t_near > t_far {
                std::mem::swap(&mut t_near, &mut t_far);
            }
            // update _t_far_ to ensure robust ray--bounds intersection
            t_far *= 1.0 as Float + 2.0 as Float * gamma(3_i32);
            if t_near > t0 {
                t0 = t_near;
            }
            if t_far < t1 {
                t1 = t_far;
            }
            if t0 > t1 {
                return false;
            }
        }
        *hitt0 = t0;
        *hitt1 = t1;
        true
    }
    pub fn intersect_p(&self, ray: &Ray, inv_dir: &Vector3f, dir_is_neg: &[u8; 3]) -> bool {
        let dir_is_neg_0: MinMaxEnum = match dir_is_neg[0] {
            0 => MinMaxEnum::Min,
            _ => MinMaxEnum::Max,
        };
        let dir_is_not_neg_0: MinMaxEnum = match dir_is_neg[0] {
            0 => MinMaxEnum::Max,
            _ => MinMaxEnum::Min,
        };
        let dir_is_neg_1: MinMaxEnum = match dir_is_neg[1] {
            0 => MinMaxEnum::Min,
            _ => MinMaxEnum::Max,
        };
        let dir_is_not_neg_1: MinMaxEnum = match dir_is_neg[1] {
            0 => MinMaxEnum::Max,
            _ => MinMaxEnum::Min,
        };
        let dir_is_neg_2: MinMaxEnum = match dir_is_neg[2] {
            0 => MinMaxEnum::Min,
            _ => MinMaxEnum::Max,
        };
        let dir_is_not_neg_2: MinMaxEnum = match dir_is_neg[2] {
            0 => MinMaxEnum::Max,
            _ => MinMaxEnum::Min,
        };
        // check for ray intersection against $x$ and $y$ slabs
        let mut t_min: Float = (self[dir_is_neg_0].x - ray.o.x) * inv_dir.x;
        let mut t_max: Float = (self[dir_is_not_neg_0].x - ray.o.x) * inv_dir.x;
        let ty_min: Float = (self[dir_is_neg_1].y - ray.o.y) * inv_dir.y;
        let mut ty_max: Float = (self[dir_is_not_neg_1].y - ray.o.y) * inv_dir.y;
        // update _t_max_ and _ty_max_ to ensure robust bounds intersection
        t_max *= 1.0 + 2.0 * gamma(3_i32);
        ty_max *= 1.0 + 2.0 * gamma(3_i32);
        if t_min > ty_max || ty_min > t_max {
            return false;
        }
        if ty_min > t_min {
            t_min = ty_min;
        }
        if ty_max < t_max {
            t_max = ty_max;
        }
        // check for ray intersection against $z$ slab
        let tz_min: Float = (self[dir_is_neg_2].z - ray.o.z) * inv_dir.z;
        let mut tz_max: Float = (self[dir_is_not_neg_2].z - ray.o.z) * inv_dir.z;
        // update _tz_max_ to ensure robust bounds intersection
        tz_max *= 1.0 + 2.0 * gamma(3_i32);
        if t_min > tz_max || tz_min > t_max {
            return false;
        }
        if tz_min > t_min {
            t_min = tz_min;
        }
        if tz_max < t_max {
            t_max = tz_max;
        }
        (t_min < ray.t_max.get()) && (t_max > 0.0)
    }
}

impl Index<MinMaxEnum> for Bounds3f {
    type Output = Point3f;
    fn index(&self, i: MinMaxEnum) -> &Point3f {
        match i {
            MinMaxEnum::Min => &self.p_min,
            _ => &self.p_max,
        }
    }
}

// /// Minimum squared distance from point to box; returns zero if point
// /// is inside.
// pub fn pnt3_distance_squared_bnd3(p: Point3f, b: Bounds3f) -> Float {
//     let dx: Float = (b.p_min.x - p.x).max(num::Zero::zero()).max(p.x - b.p_max.x);
//     let dy: Float = (b.p_min.y - p.y).max(num::Zero::zero()).max(p.y - b.p_max.y);
//     let dz: Float = (b.p_min.z - p.z).max(num::Zero::zero()).max(p.z - b.p_max.z);
//     dx * dx + dy * dy + dz * dz
// }

// /// Minimum distance from point to box; returns zero if point is
// /// inside.
// pub fn pnt3_distance_bnd3(p: Point3f, b: Bounds3f) -> Float {
//     pnt3_distance_squared_bnd3(p, b).sqrt()
// }

/// Given a bounding box and a point, the **bnd3_union_pnt3()**
/// function returns a new bounding box that encompasses that point as
/// well as the original box.
pub fn bnd3_union_pnt3f(b: &Bounds3f, p: &Point3f) -> Bounds3f {
    let p_min: Point3f = Point3f {
        x: b.p_min.x.min(p.x),
        y: b.p_min.y.min(p.y),
        z: b.p_min.z.min(p.z),
    };
    let p_max: Point3f = Point3f {
        x: b.p_max.x.max(p.x),
        y: b.p_max.y.max(p.y),
        z: b.p_max.z.max(p.z),
    };
    Bounds3f { p_min, p_max }
}

/// Construct a new box that bounds the space encompassed by two other
/// bounding boxes.
pub fn bnd3_union_bnd3f(b1: &Bounds3f, b2: &Bounds3f) -> Bounds3f {
    let p_min: Point3f = Point3f {
        x: b1.p_min.x.min(b2.p_min.x),
        y: b1.p_min.y.min(b2.p_min.y),
        z: b1.p_min.z.min(b2.p_min.z),
    };
    let p_max: Point3f = Point3f {
        x: b1.p_max.x.max(b2.p_max.x),
        y: b1.p_max.y.max(b2.p_max.y),
        z: b1.p_max.z.max(b2.p_max.z),
    };
    Bounds3f { p_min, p_max }
}

/// Determine if a given point is inside the bounding box.
pub fn pnt3_inside_bnd3(p: &Point3f, b: &Bounds3f) -> bool {
    p.x >= b.p_min.x
        && p.x <= b.p_max.x
        && p.y >= b.p_min.y
        && p.y <= b.p_max.y
        && p.z >= b.p_min.z
        && p.z <= b.p_max.z
}

/// Is a 3D point inside a 3D bound?
pub fn pnt3i_inside_exclusive(p: &Point3i, b: &Bounds3i) -> bool {
    p.x >= b.p_min.x
        && p.x < b.p_max.x
        && p.y >= b.p_min.y
        && p.y < b.p_max.y
        && p.z >= b.p_min.z
        && p.z < b.p_max.z
}

/// Is a 3D point inside a 3D bound?
pub fn pnt3f_inside_exclusive(p: &Point3f, b: &Bounds3f) -> bool {
    p.x >= b.p_min.x
        && p.x < b.p_max.x
        && p.y >= b.p_min.y
        && p.y < b.p_max.y
        && p.z >= b.p_min.z
        && p.z < b.p_max.z
}

/// Pads the bounding box by a constant factor in all dimensions.
pub fn bnd3_expand(b: &Bounds3f, delta: Float) -> Bounds3f {
    Bounds3f::new(
        b.p_min
            - Vector3f {
                x: delta,
                y: delta,
                z: delta,
            },
        b.p_max
            + Vector3f {
                x: delta,
                y: delta,
                z: delta,
            },
    )
}

#[derive(Default, Clone)]
pub struct Ray {
    /// origin
    pub o: Point3f,
    /// direction
    pub d: Vector3f,
    /// limits the ray to a segment along its infinite extent
    pub t_max: Cell<Float>,
    /// used for animations
    pub time: Float,
    pub medium: Option<Arc<Medium>>,
    /// in C++: 'class RayDifferential : public Ray'
    pub differential: Option<RayDifferential>,
}

impl Ray {
    // Point3f operator()(Float t) const { return o + d * t; }
    pub fn position(&self, t: Float) -> Point3f {
        self.o + self.d * t
    }
    // from class RayDifferential
    pub fn scale_differentials(&mut self, s: Float) {
        if let Some(d) = self.differential.iter_mut().next() {
            d.rx_origin = self.o + (d.rx_origin - self.o) * s;
            d.ry_origin = self.o + (d.ry_origin - self.o) * s;
            d.rx_direction = self.d + (d.rx_direction - self.d) * s;
            d.ry_direction = self.d + (d.ry_direction - self.d) * s;
        }
    }
}

#[derive(Debug, Default, Copy, Clone)]
pub struct RayDifferential {
    pub rx_origin: Point3f,
    pub ry_origin: Point3f,
    pub rx_direction: Vector3f,
    pub ry_direction: Vector3f,
}
