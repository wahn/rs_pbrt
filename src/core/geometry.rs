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
//! extern crate pbrt;
//!
//! use pbrt::core::geometry::Point3;
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
//! # Vectors
//!
//! **pbrt** provides both 2D and 3D **vector** classes. Both are
//! parameterized by the type of the underlying vector element, thus
//! making it easy to instantiate vectors of both integer and
//! floating-point types.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::core::geometry::Vector3;
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
//! extern crate pbrt;
//!
//! use pbrt::core::geometry::Normal3;
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
//! extern crate pbrt;
//!
//! use pbrt::core::geometry::{Ray, Point3f, Vector3f};
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
//!         t_max: std::f32::INFINITY,
//!         time: 0.0,
//!         medium: None,
//!         differential: None,
//!     };
//! }
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
//! **Option<RayDifferential>** in the **Ray** struct, which means the
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
//! extern crate pbrt;
//!
//! use pbrt::core::geometry::{Bounds3, Point3};
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

// std
use std;
use std::f32::consts::PI;
use std::ops::{Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub};
use std::sync::Arc;
// others
use num;
// pbrt
use core::medium::Medium;
use core::pbrt::Float;
use core::pbrt::{clamp_t, gamma, lerp, next_float_down, next_float_up};

// see geometry.h

pub type Point2f = Point2<Float>;
pub type Point2i = Point2<i32>;
pub type Point3f = Point3<Float>;
pub type Point3i = Point3<i32>;
pub type Vector2f = Vector2<Float>;
pub type Vector2i = Vector2<i32>;
pub type Vector3f = Vector3<Float>;
pub type Vector3i = Vector3<i32>;
pub type Normal3f = Normal3<Float>;

#[derive(Debug, Default, Copy, Clone)]
pub struct Vector2<T> {
    pub x: T,
    pub y: T,
}

impl<T> Vector2<T> {
    pub fn length_squared(&self) -> T
    where
        T: Copy + Add<T, Output = T> + Mul<T, Output = T>,
    {
        self.x * self.x + self.y * self.y
    }
    pub fn length(&self) -> T
    where
        T: num::Float,
    {
        self.length_squared().sqrt()
    }
}

impl<T> Index<u8> for Vector2<T> {
    type Output = T;
    fn index(&self, index: u8) -> &T {
        match index {
            0 => &self.x,
            1 => &self.y,
            _ => panic!("Check failed: i >= 0 && i <= 1"),
        }
    }
}

impl<T> IndexMut<u8> for Vector2<T> {
    fn index_mut(&mut self, index: u8) -> &mut T {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            _ => panic!("Check failed: i >= 0 && i <= 1"),
        }
    }
}

impl<T> MulAssign<T> for Vector2<T>
where
    T: Copy + MulAssign,
{
    fn mul_assign(&mut self, rhs: T) {
        self.x *= rhs;
        self.y *= rhs;
    }
}

/// Product of the Euclidean magnitudes of the two vectors and the
/// cosine of the angle between them. A return value of zero means
/// both vectors are orthogonal, a value if one means they are
/// codirectional.
pub fn vec2_dot<T>(v1: &Vector2<T>, v2: &Vector2<T>) -> T
where
    T: Copy + Add<T, Output = T> + Mul<T, Output = T>,
{
    v1.x * v2.x + v1.y * v2.y
}

#[derive(Debug, Default, Copy, Clone, PartialEq)]
pub struct Vector3<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T> Vector3<T> {
    pub fn abs(&self) -> Vector3<T>
    where
        T: num::Float,
    {
        Vector3::<T> {
            x: self.x.abs(),
            y: self.y.abs(),
            z: self.z.abs(),
        }
    }
    pub fn length_squared(&self) -> T
    where
        T: Copy + Add<T, Output = T> + Mul<T, Output = T>,
    {
        self.x * self.x + self.y * self.y + self.z * self.z
    }
    pub fn length(&self) -> T
    where
        T: num::Float,
    {
        self.length_squared().sqrt()
    }
}

impl<T> AddAssign<Vector3<T>> for Vector3<T>
where
    T: AddAssign,
{
    fn add_assign(&mut self, rhs: Vector3<T>) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl<T> Add for Vector3<T>
where
    T: Copy + Add<T, Output = T>,
{
    type Output = Vector3<T>;
    fn add(self, rhs: Vector3<T>) -> Vector3<T> {
        Vector3::<T> {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl<T> Sub for Vector3<T>
where
    T: Copy + Sub<T, Output = T>,
{
    type Output = Vector3<T>;
    fn sub(self, rhs: Vector3<T>) -> Vector3<T> {
        Vector3::<T> {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl<T> Mul<T> for Vector3<T>
where
    T: Copy + Mul<T, Output = T>,
{
    type Output = Vector3<T>;
    fn mul(self, rhs: T) -> Vector3<T>
    where
        T: Copy + Mul<T, Output = T>,
    {
        Vector3::<T> {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl<T> MulAssign<T> for Vector3<T>
where
    T: Copy + MulAssign,
{
    fn mul_assign(&mut self, rhs: T) {
        self.x *= rhs;
        self.y *= rhs;
        self.z *= rhs;
    }
}

// work around bug
// https://github.com/rust-lang/rust/issues/40395
impl Div<Float> for Vector3<f32> {
    type Output = Vector3<f32>;
    fn div(self, rhs: Float) -> Vector3<f32> {
        assert_ne!(rhs, 0.0 as Float);
        let inv: Float = 1.0 as Float / rhs;
        Vector3::<f32> {
            x: self.x * inv,
            y: self.y * inv,
            z: self.z * inv,
        }
    }
}

// work around bug
// https://github.com/rust-lang/rust/issues/40395
impl DivAssign<Float> for Vector3<f32> {
    fn div_assign(&mut self, rhs: Float) {
        assert_ne!(rhs, 0.0 as Float);
        let inv: Float = 1.0 as Float / rhs;
        self.x *= inv;
        self.y *= inv;
        self.z *= inv;
    }
}

impl<T> Neg for Vector3<T>
where
    T: Copy + Neg<Output = T>,
{
    type Output = Vector3<T>;
    fn neg(self) -> Vector3<T> {
        Vector3::<T> {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl<T> Index<u8> for Vector3<T> {
    type Output = T;
    fn index(&self, index: u8) -> &T {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Check failed: i >= 0 && i <= 2"),
        }
    }
}

impl<T> IndexMut<u8> for Vector3<T> {
    fn index_mut(&mut self, index: u8) -> &mut T {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => panic!("Check failed: i >= 0 && i <= 2"),
        }
    }
}

impl<T> From<Point3<T>> for Vector3<T> {
    fn from(p: Point3<T>) -> Self {
        Vector3::<T> {
            x: p.x,
            y: p.y,
            z: p.z,
        }
    }
}

impl<T> From<Normal3<T>> for Vector3<T> {
    fn from(n: Normal3<T>) -> Self {
        Vector3::<T> {
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
pub fn vec3_dot_vec3<T>(v1: &Vector3<T>, v2: &Vector3<T>) -> T
where
    T: Copy + Add<T, Output = T> + Mul<T, Output = T>,
{
    v1.x * v2.x + v1.y * v2.y + v1.z * v2.z
}

/// Product of the Euclidean magnitudes of a vector (and a normal) and
/// the cosine of the angle between them. A return value of zero means
/// both are orthogonal, a value if one means they are codirectional.
pub fn vec3_dot_nrm<T>(v1: &Vector3<T>, n2: &Normal3<T>) -> T
where
    T: Copy + Add<T, Output = T> + Mul<T, Output = T>,
{
    // DCHECK(!v1.HasNaNs() && !n2.HasNaNs());
    return v1.x * n2.x + v1.y * n2.y + v1.z * n2.z;
}

/// Computes the absolute value of the dot product.
pub fn vec3_abs_dot_vec3<T>(v1: &Vector3<T>, v2: &Vector3<T>) -> T
where
    T: num::Float,
{
    vec3_dot_vec3(v1, v2).abs()
}

/// Computes the absolute value of the dot product.
pub fn vec3_abs_dot_nrm<T>(v1: &Vector3<T>, n2: &Normal3<T>) -> T
where
    T: num::Float,
{
    vec3_dot_nrm(v1, n2).abs()
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

/// Compute a new vector pointing in the same direction but with unit
/// length.
pub fn vec3_normalize(v: &Vector3f) -> Vector3f {
    *v / v.length()
}

/// Return the largest coordinate value.
pub fn vec3_max_component<T>(v: &Vector3<T>) -> T
where
    T: num::Float,
{
    v.x.max(v.y.max(v.z))
}

/// Return the index of the component with the largest value.
pub fn vec3_max_dimension<T>(v: &Vector3<T>) -> usize
where
    T: std::cmp::PartialOrd,
{
    if v.x > v.y {
        if v.x > v.z {
            0_usize
        } else {
            2_usize
        }
    } else {
        if v.y > v.z {
            1_usize
        } else {
            2_usize
        }
    }
}

/// Permute the coordinate values according to the povided
/// permutation.
pub fn vec3_permute<T>(v: &Vector3<T>, x: usize, y: usize, z: usize) -> Vector3<T>
where
    T: Copy,
{
    let v3: Vec<T> = vec![v.x, v.y, v.z];
    let xp: T = v3[x];
    let yp: T = v3[y];
    let zp: T = v3[z];
    Vector3::<T> {
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
pub struct Point2<T> {
    pub x: T,
    pub y: T,
}

impl<T> PartialEq for Point2<T>
where
    T: std::cmp::PartialOrd,
{
    fn eq(&self, rhs: &Point2<T>) -> bool {
        if self.x == rhs.x && self.y == rhs.y {
            true
        } else {
            false
        }
    }
    fn ne(&self, rhs: &Point2<T>) -> bool {
        !self.eq(rhs)
    }
}

impl<T> Add<Point2<T>> for Point2<T>
where
    T: Add<T, Output = T>,
{
    type Output = Point2<T>;
    fn add(self, rhs: Point2<T>) -> Point2<T> {
        // TODO: DCHECK(!v.HasNaNs());
        Point2::<T> {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl<T> Add<Vector2<T>> for Point2<T>
where
    T: Add<T, Output = T>,
{
    type Output = Point2<T>;
    fn add(self, rhs: Vector2<T>) -> Point2<T> {
        // TODO: DCHECK(!v.HasNaNs());
        Point2::<T> {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl<T> Sub<Point2<T>> for Point2<T>
where
    T: Sub<T, Output = T>,
{
    type Output = Vector2<T>;
    fn sub(self, rhs: Point2<T>) -> Vector2<T> {
        Vector2::<T> {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

impl<T> Sub<Vector2<T>> for Point2<T>
where
    T: Sub<T, Output = T>,
{
    type Output = Point2<T>;
    fn sub(self, rhs: Vector2<T>) -> Point2<T> {
        Point2::<T> {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

impl<T> Mul<T> for Point2<T>
where
    T: Copy + Mul<T, Output = T>,
{
    type Output = Point2<T>;
    fn mul(self, rhs: T) -> Point2<T>
    where
        T: Copy + Mul<T, Output = T>,
    {
        Point2::<T> {
            x: self.x * rhs,
            y: self.y * rhs,
        }
    }
}

impl<T> Neg for Vector2<T>
where
    T: Copy + Neg<Output = T>,
{
    type Output = Vector2<T>;
    fn neg(self) -> Vector2<T> {
        Vector2::<T> {
            x: -self.x,
            y: -self.y,
        }
    }
}

impl<T> Index<u8> for Point2<T> {
    type Output = T;
    fn index(&self, index: u8) -> &T {
        match index {
            0 => &self.x,
            1 => &self.y,
            _ => panic!("Check failed: i >= 0 && i <= 1"),
        }
    }
}

impl<T> IndexMut<u8> for Point2<T> {
    fn index_mut(&mut self, index: u8) -> &mut T {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            _ => panic!("Check failed: i >= 0 && i <= 1"),
        }
    }
}

/// Apply floor operation component-wise.
pub fn pnt2_floor<T>(p: &Point2<T>) -> Point2<T>
where
    T: num::Float,
{
    Point2 {
        x: p.x.floor(),
        y: p.y.floor(),
    }
}

/// Apply ceil operation component-wise.
pub fn pnt2_ceil<T>(p: &Point2<T>) -> Point2<T>
where
    T: num::Float,
{
    Point2 {
        x: p.x.ceil(),
        y: p.y.ceil(),
    }
}

/// Apply std::cmp::min operation component-wise.
pub fn pnt2_min_pnt2<T>(pa: Point2<T>, pb: Point2<T>) -> Point2<T>
where
    T: Ord,
{
    Point2 {
        x: std::cmp::min(pa.x, pb.x),
        y: std::cmp::min(pa.y, pb.y),
    }
}

/// Apply std::cmp::max operation component-wise.
pub fn pnt2_max_pnt2<T>(pa: Point2<T>, pb: Point2<T>) -> Point2<T>
where
    T: Ord,
{
    Point2 {
        x: std::cmp::max(pa.x, pb.x),
        y: std::cmp::max(pa.y, pb.y),
    }
}

/// Is a 2D point inside a 2D bound?
pub fn pnt2_inside_exclusive<T>(pt: &Point2<T>, b: &Bounds2<T>) -> bool
where
    T: PartialOrd,
{
    pt.x >= b.p_min.x && pt.x < b.p_max.x && pt.y >= b.p_min.y && pt.y < b.p_max.y
}

#[derive(Debug, Default, Copy, Clone)]
pub struct Point3<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T> AddAssign<Point3<T>> for Point3<T>
where
    T: AddAssign,
{
    fn add_assign(&mut self, rhs: Point3<T>) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl<T> Add<Point3<T>> for Point3<T>
where
    T: Add<T, Output = T>,
{
    type Output = Point3<T>;
    fn add(self, rhs: Point3<T>) -> Point3<T> {
        // TODO: DCHECK(!v.HasNaNs());
        Point3::<T> {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl<T> Add<Vector3<T>> for Point3<T>
where
    T: Add<T, Output = T>,
{
    type Output = Point3<T>;
    fn add(self, rhs: Vector3<T>) -> Point3<T> {
        // TODO: DCHECK(!v.HasNaNs());
        Point3::<T> {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl<T> AddAssign<Vector3<T>> for Point3<T>
where
    T: AddAssign,
{
    fn add_assign(&mut self, rhs: Vector3<T>) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl<T> Sub<Point3<T>> for Point3<T>
where
    T: Sub<T, Output = T>,
{
    type Output = Vector3<T>;
    fn sub(self, rhs: Point3<T>) -> Vector3<T> {
        Vector3::<T> {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl<T> Sub<Vector3<T>> for Point3<T>
where
    T: Sub<T, Output = T>,
{
    type Output = Point3<T>;
    fn sub(self, rhs: Vector3<T>) -> Point3<T> {
        Point3::<T> {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl<T> Mul<T> for Point3<T>
where
    T: Copy + Mul<T, Output = T>,
{
    type Output = Point3<T>;
    fn mul(self, rhs: T) -> Point3<T>
    where
        T: Copy + Mul<T, Output = T>,
    {
        Point3::<T> {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl<T> MulAssign<T> for Point3<T>
where
    T: Copy + MulAssign,
{
    fn mul_assign(&mut self, rhs: T) {
        self.x *= rhs;
        self.y *= rhs;
        self.z *= rhs;
    }
}

// work around bug
// https://github.com/rust-lang/rust/issues/40395
impl Div<Float> for Point3<f32> {
    type Output = Point3<f32>;
    fn div(self, rhs: Float) -> Point3<f32> {
        assert_ne!(rhs, 0.0 as Float);
        let inv: Float = 1.0 as Float / rhs;
        Point3::<f32> {
            x: self.x * inv,
            y: self.y * inv,
            z: self.z * inv,
        }
    }
}

// work around bug
// https://github.com/rust-lang/rust/issues/40395
impl DivAssign<Float> for Point3<f32> {
    fn div_assign(&mut self, rhs: Float) {
        assert_ne!(rhs, 0.0 as Float);
        let inv: Float = 1.0 as Float / rhs;
        self.x *= inv;
        self.y *= inv;
        self.z *= inv;
    }
}

impl<T> Index<u8> for Point3<T> {
    type Output = T;
    fn index(&self, index: u8) -> &T {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("Check failed: i >= 0 && i <= 2"),
        }
    }
}

impl<T> IndexMut<u8> for Point3<T> {
    fn index_mut(&mut self, index: u8) -> &mut T {
        match index {
            0 => &mut self.x,
            1 => &mut self.y,
            2 => &mut self.z,
            _ => panic!("Check failed: i >= 0 && i <= 2"),
        }
    }
}

/// Permute the coordinate values according to the povided
/// permutation.
pub fn pnt3_permute<T>(v: &Point3<T>, x: usize, y: usize, z: usize) -> Point3<T>
where
    T: Copy,
{
    let v3: Vec<T> = vec![v.x, v.y, v.z];
    let xp: T = v3[x];
    let yp: T = v3[y];
    let zp: T = v3[z];
    Point3::<T> {
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
pub fn pnt3_floor<T>(p: &Point3<T>) -> Point3<T>
where
    T: num::Float,
{
    Point3 {
        x: p.x.floor(),
        y: p.y.floor(),
        z: p.z.floor(),
    }
}

/// Apply ceil operation component-wise.
pub fn pnt3_ceil<T>(p: &Point3<T>) -> Point3<T>
where
    T: num::Float,
{
    Point3 {
        x: p.x.ceil(),
        y: p.y.ceil(),
        z: p.z.ceil(),
    }
}

/// Apply abs operation component-wise.
pub fn pnt3_abs<T>(p: &Point3<T>) -> Point3<T>
where
    T: num::Float,
{
    Point3 {
        x: p.x.abs(),
        y: p.y.abs(),
        z: p.z.abs(),
    }
}

/// The distance between two points is the length of the vector
/// between them.
pub fn pnt3_distance<T>(p1: &Point3<T>, p2: &Point3<T>) -> T
where
    T: num::Float + Sub<T, Output = T>,
{
    (*p1 - *p2).length()
}

/// The distance squared between two points is the length of the
/// vector between them squared.
pub fn pnt3_distance_squared<T>(p1: &Point3<T>, p2: &Point3<T>) -> T
where
    T: num::Float + Sub<T, Output = T>,
{
    (*p1 - *p2).length_squared()
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
    //     Float d = Dot(Abs(n), pError);
    let d: Float = nrm_dot_vec3(&nrm_abs(n), p_error);
    // #ifdef PBRT_FLOAT_AS_DOUBLE
    //     // We have tons of precision; for now bump up the offset a bunch just
    //     // to be extra sure that we start on the right side of the surface
    //     // (In case of any bugs in the epsilons code...)
    //     d *= 1024.;
    // #endif
    let mut offset: Vector3f = Vector3f::from(*n) * d;
    if vec3_dot_nrm(w, n) < 0.0 as Float {
        offset = -offset;
    }
    let mut po: Point3f = *p + offset;
    // round offset point _po_ away from _p_
    for i in 0..3 {
        if offset[i] > 0.0 as Float {
            po[i] = next_float_up(po[i]);
        } else {
            if offset[i] < 0.0 as Float {
                po[i] = next_float_down(po[i]);
            }
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

#[derive(Debug, Default, Copy, Clone)]
pub struct Normal3<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T> Add for Normal3<T>
where
    T: Copy + Add<T, Output = T>,
{
    type Output = Normal3<T>;
    fn add(self, rhs: Normal3<T>) -> Normal3<T> {
        Normal3::<T> {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl<T> Sub for Normal3<T>
where
    T: Copy + Sub<T, Output = T>,
{
    type Output = Normal3<T>;
    fn sub(self, rhs: Normal3<T>) -> Normal3<T> {
        Normal3::<T> {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl<T> Mul<T> for Normal3<T>
where
    T: Copy + Mul<T, Output = T>,
{
    type Output = Normal3<T>;
    fn mul(self, rhs: T) -> Normal3<T>
    where
        T: Copy + Mul<T, Output = T>,
    {
        Normal3::<T> {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl<T> MulAssign<T> for Normal3<T>
where
    T: Copy + MulAssign,
{
    fn mul_assign(&mut self, rhs: T) {
        self.x *= rhs;
        self.y *= rhs;
        self.z *= rhs;
    }
}

impl<T> Neg for Normal3<T>
where
    T: Copy + Neg<Output = T>,
{
    type Output = Normal3<T>;
    fn neg(self) -> Normal3<T> {
        Normal3::<T> {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl<T> Normal3<T> {
    pub fn length_squared(&self) -> T
    where
        T: Copy + Add<T, Output = T> + Mul<T, Output = T>,
    {
        self.x * self.x + self.y * self.y + self.z * self.z
    }
    pub fn length(&self) -> T
    where
        T: num::Float,
    {
        self.length_squared().sqrt()
    }
}

impl<T> PartialEq for Normal3<T>
where
    T: std::cmp::PartialOrd,
{
    fn eq(&self, rhs: &Normal3<T>) -> bool {
        if self.x == rhs.x && self.y == rhs.y && self.z == rhs.z {
            true
        } else {
            false
        }
    }
    fn ne(&self, rhs: &Normal3<T>) -> bool {
        if self.x != rhs.x || self.y != rhs.y || self.z != rhs.z {
            true
        } else {
            false
        }
    }
}

// work around bug
// https://github.com/rust-lang/rust/issues/40395
impl Div<Float> for Normal3<f32> {
    type Output = Normal3<f32>;
    fn div(self, rhs: Float) -> Normal3<f32> {
        assert_ne!(rhs, 0.0 as Float);
        let inv: Float = 1.0 as Float / rhs;
        Normal3::<f32> {
            x: self.x * inv,
            y: self.y * inv,
            z: self.z * inv,
        }
    }
}

impl<T> From<Vector3<T>> for Normal3<T> {
    fn from(v: Vector3<T>) -> Self {
        // TODO: DCHECK(!v.HasNaNs());
        Normal3::<T> {
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

/// Compute a new normal pointing in the same direction but with unit
/// length.
pub fn nrm_normalize(n: &Normal3f) -> Normal3f {
    *n / n.length()
}

/// Product of the Euclidean magnitudes of a normal (and another
/// normal) and the cosine of the angle between them. A return value
/// of zero means both are orthogonal, a value if one means they are
/// codirectional.
pub fn nrm_dot_nrm<T>(n1: &Normal3<T>, n2: &Normal3<T>) -> T
where
    T: Copy + Add<T, Output = T> + Mul<T, Output = T>,
{
    // TODO: DCHECK(!n1.HasNaNs() && !n2.HasNaNs());
    n1.x * n2.x + n1.y * n2.y + n1.z * n2.z
}

/// Product of the Euclidean magnitudes of a normal (and a vector) and
/// the cosine of the angle between them. A return value of zero means
/// both are orthogonal, a value if one means they are codirectional.
pub fn nrm_dot_vec3<T>(n1: &Normal3<T>, v2: &Vector3<T>) -> T
where
    T: Copy + Add<T, Output = T> + Mul<T, Output = T>,
{
    // TODO: DCHECK(!n1.HasNaNs() && !v2.HasNaNs());
    n1.x * v2.x + n1.y * v2.y + n1.z * v2.z
}

/// Computes the absolute value of the dot product.
pub fn nrm_abs_dot_vec3<T>(n1: &Normal3<T>, v2: &Vector3<T>) -> T
where
    T: num::Float,
{
    nrm_dot_vec3(n1, v2).abs()
}

/// Return normal with the absolute value of each coordinate.
pub fn nrm_abs<T>(n: &Normal3<T>) -> Normal3<T>
where
    T: num::Float,
{
    Normal3::<T> {
        x: n.x.abs(),
        y: n.y.abs(),
        z: n.z.abs(),
    }
}

/// Flip a surface normal so that it lies in the same hemisphere as a
/// given vector.
pub fn nrm_faceforward_vec3(n: &Normal3f, v: &Vector3f) -> Normal3f {
    if nrm_dot_vec3(n, v) < 0.0 as Float {
        -(*n)
    } else {
        *n
    }
}

/// Flip a surface normal so that it lies in the same hemisphere as a
/// given normal.
pub fn nrm_faceforward_nrm(n: &Normal3f, n2: &Normal3f) -> Normal3f {
    if nrm_dot_nrm(n, n2) < 0.0 as Float {
        -(*n)
    } else {
        *n
    }
}

pub type Bounds2f = Bounds2<Float>;
pub type Bounds2i = Bounds2<i32>;
pub type Bounds3f = Bounds3<Float>;
pub type Bounds3i = Bounds3<i32>;

#[derive(Debug, Default, Copy, Clone)]
pub struct Bounds2<T> {
    pub p_min: Point2<T>,
    pub p_max: Point2<T>,
}

impl<T> Bounds2<T> {
    pub fn new(p1: Point2<T>, p2: Point2<T>) -> Self
    where
        T: Copy + Ord,
    {
        let p_min: Point2<T> = Point2::<T> {
            x: std::cmp::min(p1.x, p2.x),
            y: std::cmp::min(p1.y, p2.y),
        };
        let p_max: Point2<T> = Point2::<T> {
            x: std::cmp::max(p1.x, p2.x),
            y: std::cmp::max(p1.y, p2.y),
        };
        Bounds2::<T> {
            p_min: p_min,
            p_max: p_max,
        }
    }
    pub fn diagonal(&self) -> Vector2<T>
    where
        T: Copy + Sub<T, Output = T>,
    {
        self.p_max - self.p_min
    }
    pub fn area(&self) -> T
    where
        T: Copy + Sub<T, Output = T> + Mul<T, Output = T>,
    {
        let d: Vector2<T> = self.p_max - self.p_min;
        d.x * d.y
    }
}

impl Bounds2<Float> {
    pub fn lerp(&self, t: &Point2f) -> Point2f {
        Point2f {
            x: lerp(t.x, self.p_min.x, self.p_max.x),
            y: lerp(t.y, self.p_min.y, self.p_max.y),
        }
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
pub fn bnd2_intersect_bnd2<T>(b1: &Bounds2<T>, b2: &Bounds2<T>) -> Bounds2<T>
where
    T: Copy + Ord,
{
    Bounds2::<T> {
        p_min: Point2::<T> {
            x: std::cmp::max(b1.p_min.x, b2.p_min.x),
            y: std::cmp::max(b1.p_min.y, b2.p_min.y),
        },
        p_max: Point2::<T> {
            x: std::cmp::min(b1.p_max.x, b2.p_max.x),
            y: std::cmp::min(b1.p_max.y, b2.p_max.y),
        },
    }
}

#[derive(Debug, Copy, Clone)]
pub struct Bounds3<T> {
    pub p_min: Point3<T>,
    pub p_max: Point3<T>,
}

// work around bug
// https://github.com/rust-lang/rust/issues/40395
impl Default for Bounds3<f32> {
    fn default() -> Bounds3<f32> {
        let min_num: Float = std::f32::MIN;
        let max_num: Float = std::f32::MAX;
        // Bounds3f
        Bounds3::<f32> {
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

impl<T> Bounds3<T> {
    pub fn new(p1: Point3<T>, p2: Point3<T>) -> Self
    where
        T: num::Float,
    {
        let p_min: Point3<T> = Point3::<T> {
            x: p1.x.min(p2.x),
            y: p1.y.min(p2.y),
            z: p1.z.min(p2.z),
        };
        let p_max: Point3<T> = Point3::<T> {
            x: p1.x.max(p2.x),
            y: p1.y.max(p2.y),
            z: p1.z.max(p2.z),
        };
        Bounds3::<T> {
            p_min: p_min,
            p_max: p_max,
        }
    }
    pub fn corner(&self, corner: u8) -> Point3<T>
    where
        T: Copy,
    {
        // assert!(corner >= 0_u8);
        assert!(corner < 8_u8);
        let x: T;
        if corner & 1 == 0 {
            x = self.p_min.x;
        } else {
            x = self.p_max.x;
        }
        let y: T;
        if corner & 2 == 0 {
            y = self.p_min.y;
        } else {
            y = self.p_max.y;
        }
        let z: T;
        if corner & 4 == 0 {
            z = self.p_min.z;
        } else {
            z = self.p_max.z;
        }
        Point3::<T> { x: x, y: y, z: z }
    }
    pub fn diagonal(&self) -> Vector3<T>
    where
        T: Copy + Sub<T, Output = T>,
    {
        self.p_max - self.p_min
    }
    pub fn surface_area(&self) -> T
    where
        T: Copy + Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T>,
    {
        let d: Vector3<T> = self.diagonal();
        // 2 * (d.x * d.y + d.x * d.z + d.y * d.z)
        let r: T = d.x * d.y + d.x * d.z + d.y * d.z;
        r + r // avoid '2 *'
    }
    pub fn maximum_extent(&self) -> u8
    where
        T: Copy + std::cmp::PartialOrd + Sub<T, Output = T>,
    {
        let d: Vector3<T> = self.diagonal();
        if d.x > d.y && d.x > d.z {
            0_u8
        } else if d.y > d.z {
            1_u8
        } else {
            2_u8
        }
    }
    pub fn offset(&self, p: &Point3<T>) -> Vector3<T>
    where
        T: Copy + std::cmp::PartialOrd + Sub<T, Output = T> + DivAssign<T>,
    {
        let mut o: Vector3<T> = *p - self.p_min;
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
            *radius = pnt3_distance(&center_copy, &p_max);
        } else {
            *radius = 0.0;
        }
    }
}

impl Bounds3<Float> {
    pub fn lerp(&self, t: &Point3f) -> Point3f {
        Point3f {
            x: lerp(t.x, self.p_min.x as Float, self.p_max.x as Float),
            y: lerp(t.y, self.p_min.y as Float, self.p_max.y as Float),
            z: lerp(t.z, self.p_min.z as Float, self.p_max.z as Float),
        }
    }
    pub fn intersect_p(&self, ray: &Ray, inv_dir: &Vector3f, dir_is_neg: [u8; 3]) -> bool {
        // check for ray intersection against $x$ and $y$ slabs
        let mut t_min: Float = (self[dir_is_neg[0]].x - ray.o.x) * inv_dir.x;
        let mut t_max: Float = (self[1_u8 - dir_is_neg[0]].x - ray.o.x) * inv_dir.x;
        let ty_min: Float = (self[dir_is_neg[1]].y - ray.o.y) * inv_dir.y;
        let mut ty_max: Float = (self[1_u8 - dir_is_neg[1]].y - ray.o.y) * inv_dir.y;
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
        let tz_min: Float = (self[dir_is_neg[2]].z - ray.o.z) * inv_dir.z;
        let mut tz_max: Float = (self[1_u8 - dir_is_neg[2]].z - ray.o.z) * inv_dir.z;
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
        (t_min < ray.t_max) && (t_max > 0.0)
    }
}

impl<T> Index<u8> for Bounds3<T> {
    type Output = Point3<T>;

    fn index(&self, i: u8) -> &Point3<T> {
        match i {
            0 => &self.p_min,
            1 => &self.p_max,
            _ => panic!("Invalid index!"),
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
pub fn bnd3_union_pnt3<T>(b: &Bounds3<T>, p: &Point3<T>) -> Bounds3<T>
where
    T: num::Float,
{
    let p_min: Point3<T> = Point3::<T> {
        x: b.p_min.x.min(p.x),
        y: b.p_min.y.min(p.y),
        z: b.p_min.z.min(p.z),
    };
    let p_max: Point3<T> = Point3::<T> {
        x: b.p_max.x.max(p.x),
        y: b.p_max.y.max(p.y),
        z: b.p_max.z.max(p.z),
    };
    Bounds3::new(p_min, p_max)
}

/// Construct a new box that bounds the space encompassed by two other
/// bounding boxes.
pub fn bnd3_union_bnd3<T>(b1: &Bounds3<T>, b2: &Bounds3<T>) -> Bounds3<T>
where
    T: num::Float,
{
    let p_min: Point3<T> = Point3::<T> {
        x: b1.p_min.x.min(b2.p_min.x),
        y: b1.p_min.y.min(b2.p_min.y),
        z: b1.p_min.z.min(b2.p_min.z),
    };
    let p_max: Point3<T> = Point3::<T> {
        x: b1.p_max.x.max(b2.p_max.x),
        y: b1.p_max.y.max(b2.p_max.y),
        z: b1.p_max.z.max(b2.p_max.z),
    };
    Bounds3::new(p_min, p_max)
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

/// Pads the bounding box by a constant factor in both dimensions.
pub fn bnd3_expand(b: &Bounds3f, delta: Float) -> Bounds3f {
    Bounds3f::new(
        b.p_min - Vector3f {
            x: delta,
            y: delta,
            z: delta,
        },
        b.p_max + Vector3f {
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
    pub t_max: Float,
    /// used for animations
    pub time: Float,
    pub medium: Option<Arc<Medium + Send + Sync>>,
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
