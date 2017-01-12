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
//! ## Quaternions
//!
//! Quaternions were originally invented by Sir William Hamilton in
//! 1843 as a generalization of complex numbers (2 dimensions) to four
//! dimensions.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::Quaternion;
//!
//! fn main() {
//!     let default_quaternion: Quaternion = Quaternion::default();
//!
//!     println!("default quaternion = {:?}", default_quaternion);
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
//! ## Filters
//!
//! ### Box Filter
//!
//! One of the most commonly used filters in graphics is the box
//! filter. The box filter equally weights all samples within a square
//! region of the image. Although computational efficient, it's just
//! about the worst filter possible.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::{BoxFilter, Float, Vector2f};
//!
//! fn main() {
//!     // see box.cpp CreateBoxFilter()
//!     let xw: Float = 0.5;
//!     let yw: Float = 0.5;
//!     let box_filter = BoxFilter {
//!         radius: Vector2f { x: xw, y: yw },
//!         inv_radius: Vector2f {
//!             x: 1.0 / xw,
//!             y: 1.0 / yw,
//!         },
//!     };
//!
//!     println!("box_filter = {:?}", box_filter);
//! }
//! ```
//! ## Film
//!
//! The type of film or sensor in a camera has a dramatic effect on
//! the way that incident light is transformed into colors in an
//! image. In **pbrt**, the **Film** class models the sensing device
//! in the simulated camera. After the radiance is found for each
//! camera ray, the **Film** implementation determines the sample's
//! contribution to the pixel around the point on the film plane where
//! the camera ray began and updates its representation of the
//! image. When the main rendering loop exits, the **Film** writes the
//! final image to file.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::{Bounds2i, BoxFilter, Film, Point2i, Vector2f};
//! use std::string::String;
//!
//! fn main() {
//!     // see film.cpp CreateFilm()
//!     let film = Film {
//!         full_resolution: Point2i { x: 1280, y: 720 },
//!         diagonal: 35.0,
//!         filter: BoxFilter {
//!             radius: Vector2f { x: 0.5, y: 0.5 },
//!             inv_radius: Vector2f {
//!                 x: 1.0 / 0.5,
//!                 y: 1.0 / 0.5,
//!             },
//!         },
//!         filename: String::from("pbrt.exr"),
//!         cropped_pixel_bounds: Bounds2i {
//!             p_min: Point2i { x: 0, y: 0 },
//!             p_max: Point2i { x: 1280, y: 720 },
//!         },
//!     };
//!
//!     println!("film = {:?}", film);
//! }
//! ```
//!
//! ## Direct Lighting
//!
//! The **DirectLightingIntegrator** accounts only for direct lighting
//! - light that has traveled directly from a light source to the
//! point being shaded - and ignores indirect illumination from
//! objects that are not themselfes emissive, except for basic
//! specular reflection and transmission effects.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::DirectLightingIntegrator;
//!
//! fn main() {
//!     // see directlighting.cpp CreateDirectLightingIntegrator()
//!     let integrator = DirectLightingIntegrator::default();
//!
//!     println!("integrator = {:?}", integrator);
//! }
//! ```

extern crate num;

use std::ops::{Add, Sub, Mul, Div, Neg};
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

pub fn quadratic(a: Float, b: Float, c: Float, t0: &mut Float, t1: &mut Float) -> bool {
    // find quadratic discriminant
    let discrim: f64 = (b as f64) * (b as f64) - 4.0 * (a as f64) * (c as f64);
    if discrim < 0.0 {
        false
    } else {
        let rootDiscrim: f64 = discrim.sqrt();
        // compute quadratic _t_ values
        let mut q: f64 = 0.0;
        if b < 0.0 {
            q = -0.5 * (b - rootDiscrim);
        } else {
            q = -0.5 * (b + rootDiscrim);
        }
        *t0 = q / a;
        *t1 = c / q;
        if *t0 > *t1 {
            // std::swap(*t0, *t1);
            let swap = *t0;
            *t0 = *t1;
            *t1 = swap;
        }
        true
    }
}

// see geometry.h

pub type Point2f = Point2<Float>;
pub type Point2i = Point2<i32>;
pub type Point3f = Point3<Float>;
pub type Point3i = Point3<i32>;
pub type Vector2f = Vector2<Float>;
pub type Vector3f = Vector3<Float>;
pub type Normal3f = Normal3<Float>;

#[derive(Debug,Default,Copy,Clone)]
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

#[derive(Debug,Default,Copy,Clone)]
pub struct Vector3<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T> Sub for Vector3<T>
    where T: Copy + Sub<T, Output = T>
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
    where T: Copy + Mul<T, Output = T>
{
    type Output = Vector3<T>;
    fn mul(self, rhs: T) -> Vector3<T>
        where T: Copy + Mul<T, Output = T>
    {
        Vector3::<T> {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
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

impl<T> Neg for Vector3<T>
    where T: Copy + Neg<Output = T>
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

pub fn vec3_normalize<T>(v: Vector3<T>) -> Vector3<T>
    where T: num::Float + Copy + Div<T, Output = T>
{
    v / v.length()
}

#[derive(Debug,Default,Copy,Clone)]
pub struct Point2<T> {
    pub x: T,
    pub y: T,
}

#[derive(Debug,Default,Copy,Clone)]
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

#[derive(Debug,Default,Copy,Clone)]
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

pub type Bounds2f = Bounds2<Float>;
pub type Bounds2i = Bounds2<i32>;
pub type Bounds3f = Bounds3<Float>;
pub type Bounds3i = Bounds3<i32>;

#[derive(Debug,Default,Copy,Clone)]
pub struct Bounds2<T> {
    pub p_min: Point2<T>,
    pub p_max: Point2<T>,
}

#[derive(Debug,Default,Copy,Clone)]
pub struct Bounds3<T> {
    pub p_min: Point3<T>,
    pub p_max: Point3<T>,
}

#[derive(Debug,Default,Copy,Clone)]
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

pub fn mtx_mul(m1: Matrix4x4, m2: Matrix4x4) -> Matrix4x4 {
    let mut r: Matrix4x4 = Matrix4x4::default();
    for i in 0..4 {
        for j in 0..4 {
            r.m[i][j] = m1.m[i][0] * m2.m[0][j] + m1.m[i][1] * m2.m[1][j] +
                        m1.m[i][2] * m2.m[2][j] + m1.m[i][3] * m2.m[3][j];
        }
    }
    r
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

impl Mul for Transform {
    type Output = Transform;
    fn mul(self, rhs: Transform) -> Transform {
        Transform {
            m: mtx_mul(self.m, rhs.m),
            m_inv: mtx_mul(rhs.m_inv, self.m_inv),
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
    pub fn inverse(t: Transform) -> Transform {
        Transform {
            m: t.m_inv,
            m_inv: t.m,
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
        let a: Vector3f = vec3_normalize(axis);
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
        let dir: Vector3f = vec3_normalize(look - pos);
        if vec3_cross(vec3_normalize(up), dir).length() == 0.0 {
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
            let left: Vector3f = vec3_normalize(vec3_cross(vec3_normalize(up), dir));
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
    pub fn perspective(fov: Float, n: Float, f: Float) -> Transform {
        // perform projective divide for perspective projection
        let persp = Matrix4x4::new(1.0,
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
                                   0.0);
        // scale canonical perspective view to specified field of view
        let inv_tan_ang: Float = 1.0 / (radians(fov) / 2.0).tan();
        let scale: Transform = Transform::scale(inv_tan_ang, inv_tan_ang, 1.0);
        let persp_trans: Transform = Transform {
            m: persp,
            m_inv: Matrix4x4::inverse(persp),
        };
        scale * persp_trans
    }
}

#[derive(Debug,Default,Copy,Clone)]
pub struct DerivativeTerm {
    kc: Float,
    kx: Float,
    ky: Float,
    kz: Float,
}

#[derive(Debug,Default,Copy,Clone)]
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
    pub fn new(start_transform: &Transform,
               start_time: Float,
               end_transform: &Transform,
               end_time: Float)
               -> AnimatedTransform {
        let mut at: AnimatedTransform = AnimatedTransform::default();
        at.start_transform = start_transform.clone();
        at.end_transform = end_transform.clone();
        at.start_time = start_time;
        at.end_time = end_time;
        AnimatedTransform::decompose(&start_transform.m, &mut at.t[0], &mut at.r[0], &mut at.s[0]);
        AnimatedTransform::decompose(&end_transform.m, &mut at.t[1], &mut at.r[1], &mut at.s[1]);
        // flip _r[1]_ if needed to select shortest path
        if quat_dot(at.r[0], at.r[1]) < 0.0 {
            at.r[1] = -at.r[1];
        }
        at.has_rotation = quat_dot(at.r[0], at.r[1]) < 0.9995;
        // compute terms of motion derivative function
        if at.has_rotation {
            let cos_theta: Float = quat_dot(at.r[0], at.r[1]);
            let theta: Float = (clamp(cos_theta, -1.0, 1.0)).acos();
            let qperp: Quaternion = quat_normalize(at.r[1] - at.r[0] * cos_theta);
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
                kx: (-1.0 + q1y * q1y + q1z * q1z + qperpy * qperpy + qperpz * qperpz) * s000 +
                    q1w * q1z * s010 - qperpx * qperpy * s010 +
                    qperpw * qperpz * s010 - q1w * q1y * s020 -
                    qperpw * qperpy * s020 -
                    qperpx * qperpz * s020 + s100 -
                    q1y * q1y * s100 - q1z * q1z * s100 -
                    qperpy * qperpy * s100 -
                    qperpz * qperpz * s100 - q1w * q1z * s110 +
                    qperpx * qperpy * s110 -
                    qperpw * qperpz * s110 + q1w * q1y * s120 +
                    qperpw * qperpy * s120 +
                    qperpx * qperpz * s120 +
                    q1x * (-(q1y * s010) - q1z * s020 + q1y * s110 + q1z * s120),
                ky: (-1.0 + q1y * q1y + q1z * q1z + qperpy * qperpy + qperpz * qperpz) * s001 +
                    q1w * q1z * s011 - qperpx * qperpy * s011 +
                    qperpw * qperpz * s011 - q1w * q1y * s021 -
                    qperpw * qperpy * s021 -
                    qperpx * qperpz * s021 + s101 -
                    q1y * q1y * s101 - q1z * q1z * s101 -
                    qperpy * qperpy * s101 -
                    qperpz * qperpz * s101 - q1w * q1z * s111 +
                    qperpx * qperpy * s111 -
                    qperpw * qperpz * s111 + q1w * q1y * s121 +
                    qperpw * qperpy * s121 +
                    qperpx * qperpz * s121 +
                    q1x * (-(q1y * s011) - q1z * s021 + q1y * s111 + q1z * s121),
                kz: (-1.0 + q1y * q1y + q1z * q1z + qperpy * qperpy + qperpz * qperpz) * s002 +
                    q1w * q1z * s012 - qperpx * qperpy * s012 +
                    qperpw * qperpz * s012 - q1w * q1y * s022 -
                    qperpw * qperpy * s022 -
                    qperpx * qperpz * s022 + s102 -
                    q1y * q1y * s102 - q1z * q1z * s102 -
                    qperpy * qperpy * s102 -
                    qperpz * qperpz * s102 - q1w * q1z * s112 +
                    qperpx * qperpy * s112 -
                    qperpw * qperpz * s112 + q1w * q1y * s122 +
                    qperpw * qperpy * s122 +
                    qperpx * qperpz * s122 +
                    q1x * (-(q1y * s012) - q1z * s022 + q1y * s112 + q1z * s122),
            };
            at.c2[0] = DerivativeTerm {
                kc: 0.0,
                kx: -(qperpy * qperpy * s000) - qperpz * qperpz * s000 + qperpx * qperpy * s010 -
                    qperpw * qperpz * s010 +
                    qperpw * qperpy * s020 + qperpx * qperpz * s020 +
                    q1y * q1y * (s000 - s100) + q1z * q1z * (s000 - s100) +
                    qperpy * qperpy * s100 +
                    qperpz * qperpz * s100 - qperpx * qperpy * s110 +
                    qperpw * qperpz * s110 -
                    qperpw * qperpy * s120 - qperpx * qperpz * s120 +
                    2.0 * q1x * qperpy * s010 * theta -
                    2.0 * q1w * qperpz * s010 * theta +
                    2.0 * q1w * qperpy * s020 * theta +
                    2.0 * q1x * qperpz * s020 * theta +
                    q1y *
                    (q1x * (-s010 + s110) + q1w * (-s020 + s120) +
                     2.0 * (-2.0 * qperpy * s000 + qperpx * s010 + qperpw * s020) * theta) +
                    q1z *
                    (q1w * (s010 - s110) + q1x * (-s020 + s120) -
                     2.0 * (2.0 * qperpz * s000 + qperpw * s010 - qperpx * s020) * theta),
                ky: -(qperpy * qperpy * s001) - qperpz * qperpz * s001 + qperpx * qperpy * s011 -
                    qperpw * qperpz * s011 +
                    qperpw * qperpy * s021 + qperpx * qperpz * s021 +
                    q1y * q1y * (s001 - s101) + q1z * q1z * (s001 - s101) +
                    qperpy * qperpy * s101 +
                    qperpz * qperpz * s101 - qperpx * qperpy * s111 +
                    qperpw * qperpz * s111 -
                    qperpw * qperpy * s121 - qperpx * qperpz * s121 +
                    2.0 * q1x * qperpy * s011 * theta -
                    2.0 * q1w * qperpz * s011 * theta +
                    2.0 * q1w * qperpy * s021 * theta +
                    2.0 * q1x * qperpz * s021 * theta +
                    q1y *
                    (q1x * (-s011 + s111) + q1w * (-s021 + s121) +
                     2.0 * (-2.0 * qperpy * s001 + qperpx * s011 + qperpw * s021) * theta) +
                    q1z *
                    (q1w * (s011 - s111) + q1x * (-s021 + s121) -
                     2.0 * (2.0 * qperpz * s001 + qperpw * s011 - qperpx * s021) * theta),
                kz: -(qperpy * qperpy * s002) - qperpz * qperpz * s002 + qperpx * qperpy * s012 -
                    qperpw * qperpz * s012 +
                    qperpw * qperpy * s022 + qperpx * qperpz * s022 +
                    q1y * q1y * (s002 - s102) + q1z * q1z * (s002 - s102) +
                    qperpy * qperpy * s102 +
                    qperpz * qperpz * s102 - qperpx * qperpy * s112 +
                    qperpw * qperpz * s112 -
                    qperpw * qperpy * s122 - qperpx * qperpz * s122 +
                    2.0 * q1x * qperpy * s012 * theta -
                    2.0 * q1w * qperpz * s012 * theta +
                    2.0 * q1w * qperpy * s022 * theta +
                    2.0 * q1x * qperpz * s022 * theta +
                    q1y *
                    (q1x * (-s012 + s112) + q1w * (-s022 + s122) +
                     2.0 * (-2.0 * qperpy * s002 + qperpx * s012 + qperpw * s022) * theta) +
                    q1z *
                    (q1w * (s012 - s112) + q1x * (-s022 + s122) -
                     2.0 * (2.0 * qperpz * s002 + qperpw * s012 - qperpx * s022) * theta),
            };
            at.c3[0] = DerivativeTerm {
                kc: 0.0,
                kx: -2.0 *
                    (q1x * qperpy * s010 - q1w * qperpz * s010 + q1w * qperpy * s020 +
                     q1x * qperpz * s020 - q1x * qperpy * s110 +
                     q1w * qperpz * s110 - q1w * qperpy * s120 -
                     q1x * qperpz * s120 +
                     q1y *
                     (-2.0 * qperpy * s000 + qperpx * s010 + qperpw * s020 +
                      2.0 * qperpy * s100 - qperpx * s110 - qperpw * s120) +
                     q1z *
                     (-2.0 * qperpz * s000 - qperpw * s010 + qperpx * s020 +
                      2.0 * qperpz * s100 + qperpw * s110 - qperpx * s120)) *
                    theta,
                ky: -2.0 *
                    (q1x * qperpy * s011 - q1w * qperpz * s011 + q1w * qperpy * s021 +
                     q1x * qperpz * s021 - q1x * qperpy * s111 +
                     q1w * qperpz * s111 - q1w * qperpy * s121 -
                     q1x * qperpz * s121 +
                     q1y *
                     (-2.0 * qperpy * s001 + qperpx * s011 + qperpw * s021 +
                      2.0 * qperpy * s101 - qperpx * s111 - qperpw * s121) +
                     q1z *
                     (-2.0 * qperpz * s001 - qperpw * s011 + qperpx * s021 +
                      2.0 * qperpz * s101 + qperpw * s111 - qperpx * s121)) *
                    theta,
                kz: -2.0 *
                    (q1x * qperpy * s012 - q1w * qperpz * s012 + q1w * qperpy * s022 +
                     q1x * qperpz * s022 - q1x * qperpy * s112 +
                     q1w * qperpz * s112 - q1w * qperpy * s122 -
                     q1x * qperpz * s122 +
                     q1y *
                     (-2.0 * qperpy * s002 + qperpx * s012 + qperpw * s022 +
                      2.0 * qperpy * s102 - qperpx * s112 - qperpw * s122) +
                     q1z *
                     (-2.0 * qperpz * s002 - qperpw * s012 + qperpx * s022 +
                      2.0 * qperpz * s102 + qperpw * s112 - qperpx * s122)) *
                    theta,
            };
            at.c4[0] = DerivativeTerm {
                kc: 0.0,
                kx: -(q1x * qperpy * s010) + q1w * qperpz * s010 - q1w * qperpy * s020 -
                    q1x * qperpz * s020 + q1x * qperpy * s110 -
                    q1w * qperpz * s110 + q1w * qperpy * s120 +
                    q1x * qperpz * s120 + 2.0 * q1y * q1y * s000 * theta +
                    2.0 * q1z * q1z * s000 * theta -
                    2.0 * qperpy * qperpy * s000 * theta -
                    2.0 * qperpz * qperpz * s000 * theta +
                    2.0 * qperpx * qperpy * s010 * theta -
                    2.0 * qperpw * qperpz * s010 * theta +
                    2.0 * qperpw * qperpy * s020 * theta +
                    2.0 * qperpx * qperpz * s020 * theta +
                    q1y *
                    (-(qperpx * s010) - qperpw * s020 + 2.0 * qperpy * (s000 - s100) +
                     qperpx * s110 + qperpw * s120 -
                     2.0 * q1x * s010 * theta - 2.0 * q1w * s020 * theta) +
                    q1z *
                    (2.0 * qperpz * s000 + qperpw * s010 - qperpx * s020 - 2.0 * qperpz * s100 -
                     qperpw * s110 + qperpx * s120 +
                     2.0 * q1w * s010 * theta - 2.0 * q1x * s020 * theta),
                ky: -(q1x * qperpy * s011) + q1w * qperpz * s011 - q1w * qperpy * s021 -
                    q1x * qperpz * s021 + q1x * qperpy * s111 -
                    q1w * qperpz * s111 + q1w * qperpy * s121 +
                    q1x * qperpz * s121 + 2.0 * q1y * q1y * s001 * theta +
                    2.0 * q1z * q1z * s001 * theta -
                    2.0 * qperpy * qperpy * s001 * theta -
                    2.0 * qperpz * qperpz * s001 * theta +
                    2.0 * qperpx * qperpy * s011 * theta -
                    2.0 * qperpw * qperpz * s011 * theta +
                    2.0 * qperpw * qperpy * s021 * theta +
                    2.0 * qperpx * qperpz * s021 * theta +
                    q1y *
                    (-(qperpx * s011) - qperpw * s021 + 2.0 * qperpy * (s001 - s101) +
                     qperpx * s111 + qperpw * s121 -
                     2.0 * q1x * s011 * theta - 2.0 * q1w * s021 * theta) +
                    q1z *
                    (2.0 * qperpz * s001 + qperpw * s011 - qperpx * s021 - 2.0 * qperpz * s101 -
                     qperpw * s111 + qperpx * s121 +
                     2.0 * q1w * s011 * theta - 2.0 * q1x * s021 * theta),
                kz: -(q1x * qperpy * s012) + q1w * qperpz * s012 - q1w * qperpy * s022 -
                    q1x * qperpz * s022 + q1x * qperpy * s112 -
                    q1w * qperpz * s112 + q1w * qperpy * s122 +
                    q1x * qperpz * s122 + 2.0 * q1y * q1y * s002 * theta +
                    2.0 * q1z * q1z * s002 * theta -
                    2.0 * qperpy * qperpy * s002 * theta -
                    2.0 * qperpz * qperpz * s002 * theta +
                    2.0 * qperpx * qperpy * s012 * theta -
                    2.0 * qperpw * qperpz * s012 * theta +
                    2.0 * qperpw * qperpy * s022 * theta +
                    2.0 * qperpx * qperpz * s022 * theta +
                    q1y *
                    (-(qperpx * s012) - qperpw * s022 + 2.0 * qperpy * (s002 - s102) +
                     qperpx * s112 + qperpw * s122 -
                     2.0 * q1x * s012 * theta - 2.0 * q1w * s022 * theta) +
                    q1z *
                    (2.0 * qperpz * s002 + qperpw * s012 - qperpx * s022 - 2.0 * qperpz * s102 -
                     qperpw * s112 + qperpx * s122 +
                     2.0 * q1w * s012 * theta - 2.0 * q1x * s022 * theta),
            };
            at.c5[0] = DerivativeTerm {
                kc: 0.0,
                kx: 2.0 *
                    (qperpy * qperpy * s000 + qperpz * qperpz * s000 - qperpx * qperpy * s010 +
                     qperpw * qperpz * s010 -
                     qperpw * qperpy * s020 - qperpx * qperpz * s020 -
                     qperpy * qperpy * s100 -
                     qperpz * qperpz * s100 + q1y * q1y * (-s000 + s100) +
                     q1z * q1z * (-s000 + s100) + qperpx * qperpy * s110 -
                     qperpw * qperpz * s110 +
                     q1y * (q1x * (s010 - s110) + q1w * (s020 - s120)) +
                     qperpw * qperpy * s120 + qperpx * qperpz * s120 +
                     q1z * (-(q1w * s010) + q1x * s020 + q1w * s110 - q1x * s120)) *
                    theta,
                ky: 2.0 *
                    (qperpy * qperpy * s001 + qperpz * qperpz * s001 - qperpx * qperpy * s011 +
                     qperpw * qperpz * s011 -
                     qperpw * qperpy * s021 - qperpx * qperpz * s021 -
                     qperpy * qperpy * s101 -
                     qperpz * qperpz * s101 + q1y * q1y * (-s001 + s101) +
                     q1z * q1z * (-s001 + s101) + qperpx * qperpy * s111 -
                     qperpw * qperpz * s111 +
                     q1y * (q1x * (s011 - s111) + q1w * (s021 - s121)) +
                     qperpw * qperpy * s121 + qperpx * qperpz * s121 +
                     q1z * (-(q1w * s011) + q1x * s021 + q1w * s111 - q1x * s121)) *
                    theta,
                kz: 2.0 *
                    (qperpy * qperpy * s002 + qperpz * qperpz * s002 - qperpx * qperpy * s012 +
                     qperpw * qperpz * s012 -
                     qperpw * qperpy * s022 - qperpx * qperpz * s022 -
                     qperpy * qperpy * s102 -
                     qperpz * qperpz * s102 + q1y * q1y * (-s002 + s102) +
                     q1z * q1z * (-s002 + s102) + qperpx * qperpy * s112 -
                     qperpw * qperpz * s112 +
                     q1y * (q1x * (s012 - s112) + q1w * (s022 - s122)) +
                     qperpw * qperpy * s122 + qperpx * qperpz * s122 +
                     q1z * (-(q1w * s012) + q1x * s022 + q1w * s112 - q1x * s122)) *
                    theta,
            };
            at.c1[1] = DerivativeTerm {
                kc: -t0y + t1y,
                kx: -(qperpx * qperpy * s000) - qperpw * qperpz * s000 - s010 + q1z * q1z * s010 +
                    qperpx * qperpx * s010 +
                    qperpz * qperpz * s010 - q1y * q1z * s020 +
                    qperpw * qperpx * s020 - qperpy * qperpz * s020 +
                    qperpx * qperpy * s100 + qperpw * qperpz * s100 +
                    q1w * q1z * (-s000 + s100) +
                    q1x * q1x * (s010 - s110) + s110 - q1z * q1z * s110 -
                    qperpx * qperpx * s110 - qperpz * qperpz * s110 +
                    q1x * (q1y * (-s000 + s100) + q1w * (s020 - s120)) +
                    q1y * q1z * s120 - qperpw * qperpx * s120 +
                    qperpy * qperpz * s120,
                ky: -(qperpx * qperpy * s001) - qperpw * qperpz * s001 - s011 + q1z * q1z * s011 +
                    qperpx * qperpx * s011 +
                    qperpz * qperpz * s011 - q1y * q1z * s021 +
                    qperpw * qperpx * s021 - qperpy * qperpz * s021 +
                    qperpx * qperpy * s101 + qperpw * qperpz * s101 +
                    q1w * q1z * (-s001 + s101) +
                    q1x * q1x * (s011 - s111) + s111 - q1z * q1z * s111 -
                    qperpx * qperpx * s111 - qperpz * qperpz * s111 +
                    q1x * (q1y * (-s001 + s101) + q1w * (s021 - s121)) +
                    q1y * q1z * s121 - qperpw * qperpx * s121 +
                    qperpy * qperpz * s121,
                kz: -(qperpx * qperpy * s002) - qperpw * qperpz * s002 - s012 + q1z * q1z * s012 +
                    qperpx * qperpx * s012 +
                    qperpz * qperpz * s012 - q1y * q1z * s022 +
                    qperpw * qperpx * s022 - qperpy * qperpz * s022 +
                    qperpx * qperpy * s102 + qperpw * qperpz * s102 +
                    q1w * q1z * (-s002 + s102) +
                    q1x * q1x * (s012 - s112) + s112 - q1z * q1z * s112 -
                    qperpx * qperpx * s112 - qperpz * qperpz * s112 +
                    q1x * (q1y * (-s002 + s102) + q1w * (s022 - s122)) +
                    q1y * q1z * s122 - qperpw * qperpx * s122 +
                    qperpy * qperpz * s122,
            };
            at.c2[1] = DerivativeTerm {
                kc: 0.0,
                kx: qperpx * qperpy * s000 + qperpw * qperpz * s000 + q1z * q1z * s010 -
                    qperpx * qperpx * s010 - qperpz * qperpz * s010 -
                    q1y * q1z * s020 - qperpw * qperpx * s020 +
                    qperpy * qperpz * s020 -
                    qperpx * qperpy * s100 - qperpw * qperpz * s100 +
                    q1x * q1x * (s010 - s110) - q1z * q1z * s110 +
                    qperpx * qperpx * s110 +
                    qperpz * qperpz * s110 +
                    q1y * q1z * s120 + qperpw * qperpx * s120 -
                    qperpy * qperpz * s120 + 2.0 * q1z * qperpw * s000 * theta +
                    2.0 * q1y * qperpx * s000 * theta -
                    4.0 * q1z * qperpz * s010 * theta +
                    2.0 * q1z * qperpy * s020 * theta +
                    2.0 * q1y * qperpz * s020 * theta +
                    q1x *
                    (q1w * s020 + q1y * (-s000 + s100) - q1w * s120 + 2.0 * qperpy * s000 * theta -
                     4.0 * qperpx * s010 * theta -
                     2.0 * qperpw * s020 * theta) +
                    q1w *
                    (-(q1z * s000) + q1z * s100 + 2.0 * qperpz * s000 * theta -
                     2.0 * qperpx * s020 * theta),
                ky: qperpx * qperpy * s001 + qperpw * qperpz * s001 + q1z * q1z * s011 -
                    qperpx * qperpx * s011 - qperpz * qperpz * s011 -
                    q1y * q1z * s021 - qperpw * qperpx * s021 +
                    qperpy * qperpz * s021 -
                    qperpx * qperpy * s101 - qperpw * qperpz * s101 +
                    q1x * q1x * (s011 - s111) - q1z * q1z * s111 +
                    qperpx * qperpx * s111 +
                    qperpz * qperpz * s111 +
                    q1y * q1z * s121 + qperpw * qperpx * s121 -
                    qperpy * qperpz * s121 + 2.0 * q1z * qperpw * s001 * theta +
                    2.0 * q1y * qperpx * s001 * theta -
                    4.0 * q1z * qperpz * s011 * theta +
                    2.0 * q1z * qperpy * s021 * theta +
                    2.0 * q1y * qperpz * s021 * theta +
                    q1x *
                    (q1w * s021 + q1y * (-s001 + s101) - q1w * s121 + 2.0 * qperpy * s001 * theta -
                     4.0 * qperpx * s011 * theta -
                     2.0 * qperpw * s021 * theta) +
                    q1w *
                    (-(q1z * s001) + q1z * s101 + 2.0 * qperpz * s001 * theta -
                     2.0 * qperpx * s021 * theta),
                kz: qperpx * qperpy * s002 + qperpw * qperpz * s002 + q1z * q1z * s012 -
                    qperpx * qperpx * s012 - qperpz * qperpz * s012 -
                    q1y * q1z * s022 - qperpw * qperpx * s022 +
                    qperpy * qperpz * s022 -
                    qperpx * qperpy * s102 - qperpw * qperpz * s102 +
                    q1x * q1x * (s012 - s112) - q1z * q1z * s112 +
                    qperpx * qperpx * s112 +
                    qperpz * qperpz * s112 +
                    q1y * q1z * s122 + qperpw * qperpx * s122 -
                    qperpy * qperpz * s122 + 2.0 * q1z * qperpw * s002 * theta +
                    2.0 * q1y * qperpx * s002 * theta -
                    4.0 * q1z * qperpz * s012 * theta +
                    2.0 * q1z * qperpy * s022 * theta +
                    2.0 * q1y * qperpz * s022 * theta +
                    q1x *
                    (q1w * s022 + q1y * (-s002 + s102) - q1w * s122 + 2.0 * qperpy * s002 * theta -
                     4.0 * qperpx * s012 * theta -
                     2.0 * qperpw * s022 * theta) +
                    q1w *
                    (-(q1z * s002) + q1z * s102 + 2.0 * qperpz * s002 * theta -
                     2.0 * qperpx * s022 * theta),
            };
            at.c3[1] =
                DerivativeTerm {
                    kc: 0.0,
                    kx: 2.0 *
                        (-(q1x * qperpy * s000) - q1w * qperpz * s000 + 2.0 * q1x * qperpx * s010 +
                         q1x * qperpw * s020 + q1w * qperpx * s020 +
                         q1x * qperpy * s100 + q1w * qperpz * s100 -
                         2.0 * q1x * qperpx * s110 -
                         q1x * qperpw * s120 - q1w * qperpx * s120 +
                         q1z *
                         (2.0 * qperpz * s010 - qperpy * s020 + qperpw * (-s000 + s100) -
                          2.0 * qperpz * s110 + qperpy * s120) +
                         q1y * (-(qperpx * s000) - qperpz * s020 + qperpx * s100 + qperpz * s120)) *
                        theta,
                    ky: 2.0 *
                        (-(q1x * qperpy * s001) - q1w * qperpz * s001 + 2.0 * q1x * qperpx * s011 +
                         q1x * qperpw * s021 + q1w * qperpx * s021 +
                         q1x * qperpy * s101 + q1w * qperpz * s101 -
                         2.0 * q1x * qperpx * s111 -
                         q1x * qperpw * s121 - q1w * qperpx * s121 +
                         q1z *
                         (2.0 * qperpz * s011 - qperpy * s021 + qperpw * (-s001 + s101) -
                          2.0 * qperpz * s111 + qperpy * s121) +
                         q1y * (-(qperpx * s001) - qperpz * s021 + qperpx * s101 + qperpz * s121)) *
                        theta,
                    kz: 2.0 *
                        (-(q1x * qperpy * s002) - q1w * qperpz * s002 + 2.0 * q1x * qperpx * s012 +
                         q1x * qperpw * s022 + q1w * qperpx * s022 +
                         q1x * qperpy * s102 + q1w * qperpz * s102 -
                         2.0 * q1x * qperpx * s112 -
                         q1x * qperpw * s122 - q1w * qperpx * s122 +
                         q1z *
                         (2.0 * qperpz * s012 - qperpy * s022 + qperpw * (-s002 + s102) -
                          2.0 * qperpz * s112 + qperpy * s122) +
                         q1y * (-(qperpx * s002) - qperpz * s022 + qperpx * s102 + qperpz * s122)) *
                        theta,
                };
            at.c4[1] = DerivativeTerm {
                kc: 0.0,
                kx: -(q1x * qperpy * s000) - q1w * qperpz * s000 + 2.0 * q1x * qperpx * s010 +
                    q1x * qperpw * s020 + q1w * qperpx * s020 +
                    q1x * qperpy * s100 +
                    q1w * qperpz * s100 - 2.0 * q1x * qperpx * s110 -
                    q1x * qperpw * s120 -
                    q1w * qperpx * s120 + 2.0 * qperpx * qperpy * s000 * theta +
                    2.0 * qperpw * qperpz * s000 * theta +
                    2.0 * q1x * q1x * s010 * theta +
                    2.0 * q1z * q1z * s010 * theta -
                    2.0 * qperpx * qperpx * s010 * theta -
                    2.0 * qperpz * qperpz * s010 * theta +
                    2.0 * q1w * q1x * s020 * theta -
                    2.0 * qperpw * qperpx * s020 * theta +
                    2.0 * qperpy * qperpz * s020 * theta +
                    q1y *
                    (-(qperpx * s000) - qperpz * s020 + qperpx * s100 + qperpz * s120 -
                     2.0 * q1x * s000 * theta) +
                    q1z *
                    (2.0 * qperpz * s010 - qperpy * s020 + qperpw * (-s000 + s100) -
                     2.0 * qperpz * s110 + qperpy * s120 -
                     2.0 * q1w * s000 * theta - 2.0 * q1y * s020 * theta),
                ky: -(q1x * qperpy * s001) - q1w * qperpz * s001 + 2.0 * q1x * qperpx * s011 +
                    q1x * qperpw * s021 + q1w * qperpx * s021 +
                    q1x * qperpy * s101 +
                    q1w * qperpz * s101 - 2.0 * q1x * qperpx * s111 -
                    q1x * qperpw * s121 -
                    q1w * qperpx * s121 + 2.0 * qperpx * qperpy * s001 * theta +
                    2.0 * qperpw * qperpz * s001 * theta +
                    2.0 * q1x * q1x * s011 * theta +
                    2.0 * q1z * q1z * s011 * theta -
                    2.0 * qperpx * qperpx * s011 * theta -
                    2.0 * qperpz * qperpz * s011 * theta +
                    2.0 * q1w * q1x * s021 * theta -
                    2.0 * qperpw * qperpx * s021 * theta +
                    2.0 * qperpy * qperpz * s021 * theta +
                    q1y *
                    (-(qperpx * s001) - qperpz * s021 + qperpx * s101 + qperpz * s121 -
                     2.0 * q1x * s001 * theta) +
                    q1z *
                    (2.0 * qperpz * s011 - qperpy * s021 + qperpw * (-s001 + s101) -
                     2.0 * qperpz * s111 + qperpy * s121 -
                     2.0 * q1w * s001 * theta - 2.0 * q1y * s021 * theta),
                kz: -(q1x * qperpy * s002) - q1w * qperpz * s002 + 2.0 * q1x * qperpx * s012 +
                    q1x * qperpw * s022 + q1w * qperpx * s022 +
                    q1x * qperpy * s102 +
                    q1w * qperpz * s102 - 2.0 * q1x * qperpx * s112 -
                    q1x * qperpw * s122 -
                    q1w * qperpx * s122 + 2.0 * qperpx * qperpy * s002 * theta +
                    2.0 * qperpw * qperpz * s002 * theta +
                    2.0 * q1x * q1x * s012 * theta +
                    2.0 * q1z * q1z * s012 * theta -
                    2.0 * qperpx * qperpx * s012 * theta -
                    2.0 * qperpz * qperpz * s012 * theta +
                    2.0 * q1w * q1x * s022 * theta -
                    2.0 * qperpw * qperpx * s022 * theta +
                    2.0 * qperpy * qperpz * s022 * theta +
                    q1y *
                    (-(qperpx * s002) - qperpz * s022 + qperpx * s102 + qperpz * s122 -
                     2.0 * q1x * s002 * theta) +
                    q1z *
                    (2.0 * qperpz * s012 - qperpy * s022 + qperpw * (-s002 + s102) -
                     2.0 * qperpz * s112 + qperpy * s122 -
                     2.0 * q1w * s002 * theta - 2.0 * q1y * s022 * theta),
            };
            at.c5[1] = DerivativeTerm {
                kc: 0.,
                kx: -2.0 *
                    (qperpx * qperpy * s000 + qperpw * qperpz * s000 + q1z * q1z * s010 -
                     qperpx * qperpx * s010 - qperpz * qperpz * s010 -
                     q1y * q1z * s020 - qperpw * qperpx * s020 +
                     qperpy * qperpz * s020 -
                     qperpx * qperpy * s100 - qperpw * qperpz * s100 +
                     q1w * q1z * (-s000 + s100) +
                     q1x * q1x * (s010 - s110) - q1z * q1z * s110 +
                     qperpx * qperpx * s110 + qperpz * qperpz * s110 +
                     q1x * (q1y * (-s000 + s100) + q1w * (s020 - s120)) +
                     q1y * q1z * s120 + qperpw * qperpx * s120 -
                     qperpy * qperpz * s120) * theta,
                ky: -2.0 *
                    (qperpx * qperpy * s001 + qperpw * qperpz * s001 + q1z * q1z * s011 -
                     qperpx * qperpx * s011 - qperpz * qperpz * s011 -
                     q1y * q1z * s021 - qperpw * qperpx * s021 +
                     qperpy * qperpz * s021 -
                     qperpx * qperpy * s101 - qperpw * qperpz * s101 +
                     q1w * q1z * (-s001 + s101) +
                     q1x * q1x * (s011 - s111) - q1z * q1z * s111 +
                     qperpx * qperpx * s111 + qperpz * qperpz * s111 +
                     q1x * (q1y * (-s001 + s101) + q1w * (s021 - s121)) +
                     q1y * q1z * s121 + qperpw * qperpx * s121 -
                     qperpy * qperpz * s121) * theta,
                kz: -2.0 *
                    (qperpx * qperpy * s002 + qperpw * qperpz * s002 + q1z * q1z * s012 -
                     qperpx * qperpx * s012 - qperpz * qperpz * s012 -
                     q1y * q1z * s022 - qperpw * qperpx * s022 +
                     qperpy * qperpz * s022 -
                     qperpx * qperpy * s102 - qperpw * qperpz * s102 +
                     q1w * q1z * (-s002 + s102) +
                     q1x * q1x * (s012 - s112) - q1z * q1z * s112 +
                     qperpx * qperpx * s112 + qperpz * qperpz * s112 +
                     q1x * (q1y * (-s002 + s102) + q1w * (s022 - s122)) +
                     q1y * q1z * s122 + qperpw * qperpx * s122 -
                     qperpy * qperpz * s122) * theta,
            };
            at.c1[2] = DerivativeTerm {
                kc: -t0z + t1z,
                kx: (qperpw * qperpy * s000 - qperpx * qperpz * s000 - q1y * q1z * s010 -
                     qperpw * qperpx * s010 - qperpy * qperpz * s010 - s020 +
                     q1y * q1y * s020 +
                     qperpx * qperpx * s020 + qperpy * qperpy * s020 -
                     qperpw * qperpy * s100 +
                     qperpx * qperpz * s100 + q1x * q1z * (-s000 + s100) +
                     q1y * q1z * s110 +
                     qperpw * qperpx * s110 + qperpy * qperpz * s110 +
                     q1w * (q1y * (s000 - s100) + q1x * (-s010 + s110)) +
                     q1x * q1x * (s020 - s120) + s120 - q1y * q1y * s120 -
                     qperpx * qperpx * s120 - qperpy * qperpy * s120),
                ky: (qperpw * qperpy * s001 - qperpx * qperpz * s001 - q1y * q1z * s011 -
                     qperpw * qperpx * s011 - qperpy * qperpz * s011 - s021 +
                     q1y * q1y * s021 +
                     qperpx * qperpx * s021 + qperpy * qperpy * s021 -
                     qperpw * qperpy * s101 +
                     qperpx * qperpz * s101 + q1x * q1z * (-s001 + s101) +
                     q1y * q1z * s111 +
                     qperpw * qperpx * s111 + qperpy * qperpz * s111 +
                     q1w * (q1y * (s001 - s101) + q1x * (-s011 + s111)) +
                     q1x * q1x * (s021 - s121) + s121 - q1y * q1y * s121 -
                     qperpx * qperpx * s121 - qperpy * qperpy * s121),
                kz: (qperpw * qperpy * s002 - qperpx * qperpz * s002 - q1y * q1z * s012 -
                     qperpw * qperpx * s012 - qperpy * qperpz * s012 - s022 +
                     q1y * q1y * s022 +
                     qperpx * qperpx * s022 + qperpy * qperpy * s022 -
                     qperpw * qperpy * s102 +
                     qperpx * qperpz * s102 + q1x * q1z * (-s002 + s102) +
                     q1y * q1z * s112 +
                     qperpw * qperpx * s112 + qperpy * qperpz * s112 +
                     q1w * (q1y * (s002 - s102) + q1x * (-s012 + s112)) +
                     q1x * q1x * (s022 - s122) + s122 - q1y * q1y * s122 -
                     qperpx * qperpx * s122 - qperpy * qperpy * s122),
            };
            at.c2[2] = DerivativeTerm {
                kc: 0.0,
                kx: (q1w * q1y * s000 - q1x * q1z * s000 - qperpw * qperpy * s000 +
                     qperpx * qperpz * s000 - q1w * q1x * s010 -
                     q1y * q1z * s010 + qperpw * qperpx * s010 +
                     qperpy * qperpz * s010 + q1x * q1x * s020 +
                     q1y * q1y * s020 - qperpx * qperpx * s020 -
                     qperpy * qperpy * s020 - q1w * q1y * s100 +
                     q1x * q1z * s100 + qperpw * qperpy * s100 -
                     qperpx * qperpz * s100 + q1w * q1x * s110 +
                     q1y * q1z * s110 - qperpw * qperpx * s110 -
                     qperpy * qperpz * s110 - q1x * q1x * s120 -
                     q1y * q1y * s120 + qperpx * qperpx * s120 +
                     qperpy * qperpy * s120 -
                     2.0 * q1y * qperpw * s000 * theta +
                     2.0 * q1z * qperpx * s000 * theta -
                     2.0 * q1w * qperpy * s000 * theta +
                     2.0 * q1x * qperpz * s000 * theta +
                     2.0 * q1x * qperpw * s010 * theta +
                     2.0 * q1w * qperpx * s010 * theta +
                     2.0 * q1z * qperpy * s010 * theta +
                     2.0 * q1y * qperpz * s010 * theta -
                     4.0 * q1x * qperpx * s020 * theta -
                     4.0 * q1y * qperpy * s020 * theta),
                ky: (q1w * q1y * s001 - q1x * q1z * s001 - qperpw * qperpy * s001 +
                     qperpx * qperpz * s001 - q1w * q1x * s011 -
                     q1y * q1z * s011 + qperpw * qperpx * s011 +
                     qperpy * qperpz * s011 + q1x * q1x * s021 +
                     q1y * q1y * s021 - qperpx * qperpx * s021 -
                     qperpy * qperpy * s021 - q1w * q1y * s101 +
                     q1x * q1z * s101 + qperpw * qperpy * s101 -
                     qperpx * qperpz * s101 + q1w * q1x * s111 +
                     q1y * q1z * s111 - qperpw * qperpx * s111 -
                     qperpy * qperpz * s111 - q1x * q1x * s121 -
                     q1y * q1y * s121 + qperpx * qperpx * s121 +
                     qperpy * qperpy * s121 -
                     2.0 * q1y * qperpw * s001 * theta +
                     2.0 * q1z * qperpx * s001 * theta -
                     2.0 * q1w * qperpy * s001 * theta +
                     2.0 * q1x * qperpz * s001 * theta +
                     2.0 * q1x * qperpw * s011 * theta +
                     2.0 * q1w * qperpx * s011 * theta +
                     2.0 * q1z * qperpy * s011 * theta +
                     2.0 * q1y * qperpz * s011 * theta -
                     4.0 * q1x * qperpx * s021 * theta -
                     4.0 * q1y * qperpy * s021 * theta),
                kz: (q1w * q1y * s002 - q1x * q1z * s002 - qperpw * qperpy * s002 +
                     qperpx * qperpz * s002 - q1w * q1x * s012 -
                     q1y * q1z * s012 + qperpw * qperpx * s012 +
                     qperpy * qperpz * s012 + q1x * q1x * s022 +
                     q1y * q1y * s022 - qperpx * qperpx * s022 -
                     qperpy * qperpy * s022 - q1w * q1y * s102 +
                     q1x * q1z * s102 + qperpw * qperpy * s102 -
                     qperpx * qperpz * s102 + q1w * q1x * s112 +
                     q1y * q1z * s112 - qperpw * qperpx * s112 -
                     qperpy * qperpz * s112 - q1x * q1x * s122 -
                     q1y * q1y * s122 + qperpx * qperpx * s122 +
                     qperpy * qperpy * s122 -
                     2.0 * q1y * qperpw * s002 * theta +
                     2.0 * q1z * qperpx * s002 * theta -
                     2.0 * q1w * qperpy * s002 * theta +
                     2.0 * q1x * qperpz * s002 * theta +
                     2.0 * q1x * qperpw * s012 * theta +
                     2.0 * q1w * qperpx * s012 * theta +
                     2.0 * q1z * qperpy * s012 * theta +
                     2.0 * q1y * qperpz * s012 * theta -
                     4.0 * q1x * qperpx * s022 * theta -
                     4.0 * q1y * qperpy * s022 * theta),
            };
            at.c3[2] =
                DerivativeTerm {
                    kc: 0.0,
                    kx: -2.0 *
                        (-(q1w * qperpy * s000) + q1x * qperpz * s000 + q1x * qperpw * s010 +
                         q1w * qperpx * s010 - 2.0 * q1x * qperpx * s020 +
                         q1w * qperpy * s100 - q1x * qperpz * s100 -
                         q1x * qperpw * s110 - q1w * qperpx * s110 +
                         q1z * (qperpx * s000 + qperpy * s010 - qperpx * s100 - qperpy * s110) +
                         2.0 * q1x * qperpx * s120 +
                         q1y *
                         (qperpz * s010 - 2.0 * qperpy * s020 + qperpw * (-s000 + s100) -
                          qperpz * s110 + 2.0 * qperpy * s120)) * theta,
                    ky: -2.0 *
                        (-(q1w * qperpy * s001) + q1x * qperpz * s001 + q1x * qperpw * s011 +
                         q1w * qperpx * s011 - 2.0 * q1x * qperpx * s021 +
                         q1w * qperpy * s101 - q1x * qperpz * s101 -
                         q1x * qperpw * s111 - q1w * qperpx * s111 +
                         q1z * (qperpx * s001 + qperpy * s011 - qperpx * s101 - qperpy * s111) +
                         2.0 * q1x * qperpx * s121 +
                         q1y *
                         (qperpz * s011 - 2.0 * qperpy * s021 + qperpw * (-s001 + s101) -
                          qperpz * s111 + 2.0 * qperpy * s121)) * theta,
                    kz: -2.0 *
                        (-(q1w * qperpy * s002) + q1x * qperpz * s002 + q1x * qperpw * s012 +
                         q1w * qperpx * s012 - 2.0 * q1x * qperpx * s022 +
                         q1w * qperpy * s102 - q1x * qperpz * s102 -
                         q1x * qperpw * s112 - q1w * qperpx * s112 +
                         q1z * (qperpx * s002 + qperpy * s012 - qperpx * s102 - qperpy * s112) +
                         2.0 * q1x * qperpx * s122 +
                         q1y *
                         (qperpz * s012 - 2.0 * qperpy * s022 + qperpw * (-s002 + s102) -
                          qperpz * s112 + 2.0 * qperpy * s122)) * theta,
                };
            at.c4[2] = DerivativeTerm {
                kc: 0.0,
                kx: q1w * qperpy * s000 - q1x * qperpz * s000 - q1x * qperpw * s010 -
                    q1w * qperpx * s010 + 2.0 * q1x * qperpx * s020 -
                    q1w * qperpy * s100 + q1x * qperpz * s100 +
                    q1x * qperpw * s110 + q1w * qperpx * s110 -
                    2.0 * q1x * qperpx * s120 -
                    2.0 * qperpw * qperpy * s000 * theta +
                    2.0 * qperpx * qperpz * s000 * theta -
                    2.0 * q1w * q1x * s010 * theta +
                    2.0 * qperpw * qperpx * s010 * theta +
                    2.0 * qperpy * qperpz * s010 * theta +
                    2.0 * q1x * q1x * s020 * theta +
                    2.0 * q1y * q1y * s020 * theta -
                    2.0 * qperpx * qperpx * s020 * theta -
                    2.0 * qperpy * qperpy * s020 * theta +
                    q1z *
                    (-(qperpx * s000) - qperpy * s010 + qperpx * s100 + qperpy * s110 -
                     2.0 * q1x * s000 * theta) +
                    q1y *
                    (-(qperpz * s010) + 2.0 * qperpy * s020 + qperpw * (s000 - s100) +
                     qperpz * s110 - 2.0 * qperpy * s120 +
                     2.0 * q1w * s000 * theta - 2.0 * q1z * s010 * theta),
                ky: q1w * qperpy * s001 - q1x * qperpz * s001 - q1x * qperpw * s011 -
                    q1w * qperpx * s011 + 2.0 * q1x * qperpx * s021 -
                    q1w * qperpy * s101 + q1x * qperpz * s101 +
                    q1x * qperpw * s111 + q1w * qperpx * s111 -
                    2.0 * q1x * qperpx * s121 -
                    2.0 * qperpw * qperpy * s001 * theta +
                    2.0 * qperpx * qperpz * s001 * theta -
                    2.0 * q1w * q1x * s011 * theta +
                    2.0 * qperpw * qperpx * s011 * theta +
                    2.0 * qperpy * qperpz * s011 * theta +
                    2.0 * q1x * q1x * s021 * theta +
                    2.0 * q1y * q1y * s021 * theta -
                    2.0 * qperpx * qperpx * s021 * theta -
                    2.0 * qperpy * qperpy * s021 * theta +
                    q1z *
                    (-(qperpx * s001) - qperpy * s011 + qperpx * s101 + qperpy * s111 -
                     2.0 * q1x * s001 * theta) +
                    q1y *
                    (-(qperpz * s011) + 2.0 * qperpy * s021 + qperpw * (s001 - s101) +
                     qperpz * s111 - 2.0 * qperpy * s121 +
                     2.0 * q1w * s001 * theta - 2.0 * q1z * s011 * theta),
                kz: q1w * qperpy * s002 - q1x * qperpz * s002 - q1x * qperpw * s012 -
                    q1w * qperpx * s012 + 2.0 * q1x * qperpx * s022 -
                    q1w * qperpy * s102 + q1x * qperpz * s102 +
                    q1x * qperpw * s112 + q1w * qperpx * s112 -
                    2.0 * q1x * qperpx * s122 -
                    2.0 * qperpw * qperpy * s002 * theta +
                    2.0 * qperpx * qperpz * s002 * theta -
                    2.0 * q1w * q1x * s012 * theta +
                    2.0 * qperpw * qperpx * s012 * theta +
                    2.0 * qperpy * qperpz * s012 * theta +
                    2.0 * q1x * q1x * s022 * theta +
                    2.0 * q1y * q1y * s022 * theta -
                    2.0 * qperpx * qperpx * s022 * theta -
                    2.0 * qperpy * qperpy * s022 * theta +
                    q1z *
                    (-(qperpx * s002) - qperpy * s012 + qperpx * s102 + qperpy * s112 -
                     2.0 * q1x * s002 * theta) +
                    q1y *
                    (-(qperpz * s012) + 2.0 * qperpy * s022 + qperpw * (s002 - s102) +
                     qperpz * s112 - 2.0 * qperpy * s122 +
                     2.0 * q1w * s002 * theta - 2.0 * q1z * s012 * theta),
            };
            at.c5[2] = DerivativeTerm {
                kc: 0.,
                kx: 2.0 *
                    (qperpw * qperpy * s000 - qperpx * qperpz * s000 + q1y * q1z * s010 -
                     qperpw * qperpx * s010 - qperpy * qperpz * s010 -
                     q1y * q1y * s020 + qperpx * qperpx * s020 +
                     qperpy * qperpy * s020 +
                     q1x * q1z * (s000 - s100) - qperpw * qperpy * s100 +
                     qperpx * qperpz * s100 +
                     q1w * (q1y * (-s000 + s100) + q1x * (s010 - s110)) -
                     q1y * q1z * s110 + qperpw * qperpx * s110 +
                     qperpy * qperpz * s110 + q1y * q1y * s120 -
                     qperpx * qperpx * s120 - qperpy * qperpy * s120 +
                     q1x * q1x * (-s020 + s120)) * theta,
                ky: 2.0 *
                    (qperpw * qperpy * s001 - qperpx * qperpz * s001 + q1y * q1z * s011 -
                     qperpw * qperpx * s011 - qperpy * qperpz * s011 -
                     q1y * q1y * s021 + qperpx * qperpx * s021 +
                     qperpy * qperpy * s021 +
                     q1x * q1z * (s001 - s101) - qperpw * qperpy * s101 +
                     qperpx * qperpz * s101 +
                     q1w * (q1y * (-s001 + s101) + q1x * (s011 - s111)) -
                     q1y * q1z * s111 + qperpw * qperpx * s111 +
                     qperpy * qperpz * s111 + q1y * q1y * s121 -
                     qperpx * qperpx * s121 - qperpy * qperpy * s121 +
                     q1x * q1x * (-s021 + s121)) * theta,
                kz: 2.0 *
                    (qperpw * qperpy * s002 - qperpx * qperpz * s002 + q1y * q1z * s012 -
                     qperpw * qperpx * s012 - qperpy * qperpz * s012 -
                     q1y * q1y * s022 + qperpx * qperpx * s022 +
                     qperpy * qperpy * s022 +
                     q1x * q1z * (s002 - s102) - qperpw * qperpy * s102 +
                     qperpx * qperpz * s102 +
                     q1w * (q1y * (-s002 + s102) + q1x * (s012 - s112)) -
                     q1y * q1z * s112 + qperpw * qperpx * s112 +
                     qperpy * qperpz * s112 + q1y * q1y * s122 -
                     qperpx * qperpx * s122 - qperpy * qperpy * s122 +
                     q1x * q1x * (-s022 + s122)) * theta,
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
        let mut matrix: Matrix4x4 = m.clone();
        for i in 0..3 {
            matrix.m[i][3] = 0.0;
            matrix.m[3][i] = 0.0;
        }
        matrix.m[3][3] = 1.0;
        // extract rotation _r_ from transformation matrix
        let mut norm: Float;
        let mut count: u8 = 0;
        let mut r: Matrix4x4 = matrix.clone();
        loop {
            // compute next matrix _rnext_ in series
            let mut rnext: Matrix4x4 = Matrix4x4::default();
            let rit: Matrix4x4 = Matrix4x4::inverse(Matrix4x4::transpose(r));
            for i in 0..4 {
                for j in 0..4 {
                    rnext.m[i][j] = 0.5 * (r.m[i][j] + rit.m[i][j]);
                }
            }
            // compute norm of difference between _r_ and _rnext_
            norm = 0.0;
            for i in 0..3 {
                let n: Float = (r.m[i][0] - rnext.m[i][0]).abs() +
                               (r.m[i][1] - rnext.m[i][1]).abs() +
                               (r.m[i][2] - rnext.m[i][2]).abs();
                norm = norm.max(n);
            }
            r = rnext.clone();
            count += 1;
            if count >= 100 || norm <= 0.0001 {
                break;
            }
        }
        // XXX TODO FIXME deal with flip...
        let transform: Transform = Transform {
            m: r.clone(),
            m_inv: Matrix4x4::inverse(r.clone()),
        };
        *rquat = Quaternion::new(transform);

        // compute scale _S_ using rotation and original matrix
        *s = mtx_mul(Matrix4x4::inverse(r), *m);
    }
}

// see quaternion.h

#[derive(Debug,Copy,Clone)]
pub struct Quaternion {
    pub v: Vector3f,
    pub w: Float,
}

impl Default for Quaternion {
    fn default() -> Quaternion {
        Quaternion {
            v: Vector3f::default(),
            w: 1.0,
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

impl Quaternion {
    pub fn new(t: Transform) -> Quaternion {
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
}

pub fn quat_dot(q1: Quaternion, q2: Quaternion) -> Float {
    vec3_dot(q1.v, q2.v) + q1.w * q2.w
}

pub fn quat_normalize(q: Quaternion) -> Quaternion {
    q / quat_dot(q, q)
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

// see zerotwosequence.h

#[derive(Debug,Default,Copy,Clone)]
pub struct ZeroTwoSequenceSampler {
    pub samples_per_pixel: i64,
    pub n_sampled_dimensions: i64,
}

// see box.h

#[derive(Debug,Default,Copy,Clone)]
pub struct BoxFilter {
    pub radius: Vector2f,
    pub inv_radius: Vector2f,
}

// see film.h

#[derive(Debug,Default,Clone)]
pub struct Film {
    pub full_resolution: Point2i,
    pub diagonal: Float,
    pub filter: BoxFilter, // TODO: Filter
    pub filename: String,
    pub cropped_pixel_bounds: Bounds2i,
}

// see perspective.h

#[derive(Debug,Default,Clone)]
pub struct PerspectiveCamera {
    // inherited from Camera (see camera.h)
    pub camera_to_world: AnimatedTransform,
    pub shutter_open: Float,
    pub shutter_close: Float,
    pub film: Film,
    // TODO: const Medium *medium;
    // inherited from ProjectiveCamera (see camera.h)
    camera_to_screen: Transform,
    raster_to_camera: Transform,
    screen_to_raster: Transform,
    raster_to_screen: Transform,
    lens_radius: Float,
    focal_distance: Float,
    // private data (see perspective.h)
    // dx_camera: Vector3f,
    // dy_camera: Vector3f,
    // a: Float,
}

impl PerspectiveCamera {
    pub fn new(camera_to_world: AnimatedTransform,
               camera_to_screen: Transform,
               screen_window: Bounds2f,
               shutter_open: Float,
               shutter_close: Float,
               lens_radius: Float,
               focal_distance: Float,
               fov: Float,
               film: Film /* const Medium *medium */)
               -> PerspectiveCamera {
        // see perspective.cpp
        // compute differential changes in origin for perspective camera rays
        // let dxCamera = (RasterToCamera(Point3f(1, 0, 0)) - RasterToCamera(Point3f(0, 0, 0)));
        // let dyCamera = (RasterToCamera(Point3f(0, 1, 0)) - RasterToCamera(Point3f(0, 0, 0)));
        // compute image plane bounds at $z=1$ for _PerspectiveCamera_
        // Point2i res = film->fullResolution;
        // Point3f pMin = RasterToCamera(Point3f(0, 0, 0));
        // Point3f pMax = RasterToCamera(Point3f(res.x, res.y, 0));
        // pMin /= pMin.z;
        // pMax /= pMax.z;
        // A = std::abs((pMax.x - pMin.x) * (pMax.y - pMin.y));
        // see camera.h
        // compute projective camera screen transformations
        let scale1 = Transform::scale(film.full_resolution.x as Float,
                                      film.full_resolution.y as Float,
                                      1.0);
        let scale2 = Transform::scale(1.0 / (screen_window.p_max.x - screen_window.p_min.x),
                                      1.0 / (screen_window.p_min.y - screen_window.p_max.y),
                                      1.0);
        let translate = Transform::translate(Vector3f {
            x: -screen_window.p_min.x,
            y: -screen_window.p_max.y,
            z: 0.0,
        });
        let screen_to_raster = scale1 * scale2 * translate;
        let raster_to_screen = Transform::inverse(screen_to_raster);
        let raster_to_camera = Transform::inverse(camera_to_screen) * raster_to_screen;

        PerspectiveCamera {
            camera_to_world: camera_to_world,
            shutter_open: shutter_open,
            shutter_close: shutter_close,
            film: film,
            camera_to_screen: camera_to_screen,
            raster_to_camera: raster_to_camera,
            screen_to_raster: screen_to_raster,
            raster_to_screen: raster_to_screen,
            lens_radius: lens_radius,
            focal_distance: focal_distance,
            // dx_camera: dx_camera,
            // dy_camera: dy_camera,
            // a: a,
        }
    }
}

// see directlighting.h

#[derive(Debug,Clone)]
pub enum LightStrategy {
    UniformSampleAll,
    UniformSampleOne,
}

#[derive(Debug,Clone)]
pub struct DirectLightingIntegrator {
    // inherited from SamplerIntegrator (see integrator.h)
    camera: PerspectiveCamera, // std::shared_ptr<const Camera> camera;
    // TODO: std::shared_ptr<Sampler> sampler;
    pixel_bounds: Bounds2i,
    // see directlighting.h
    // TODO: const LightStrategy strategy;
    strategy: LightStrategy,
    max_depth: i64,
    n_light_samples: Vec<i64>,
}

impl Default for DirectLightingIntegrator {
    fn default() -> DirectLightingIntegrator {
        DirectLightingIntegrator {
            camera: PerspectiveCamera::default(),
            pixel_bounds: Bounds2i {
                p_min: Point2i { x: 0, y: 0 },
                p_max: Point2i { x: 1280, y: 720 },
            },
            strategy: LightStrategy::UniformSampleAll,
            max_depth: 5,
            n_light_samples: Vec::new(),
        }
    }
}

// see paramset.h

pub struct ParamSetItem<T> {
    pub name: String,
    pub values: Vec<T>,
    pub n_values: usize,
    pub looked_up: bool, // false
}

pub struct ParamSet {
    bools: Vec<ParamSetItem<bool>>,
    ints: Vec<ParamSetItem<i64>>,
    floats: Vec<ParamSetItem<Float>>,
    point2fs: Vec<ParamSetItem<Point2f>>,
    vector2fs: Vec<ParamSetItem<Vector2f>>,
    point3fs: Vec<ParamSetItem<Point3f>>,
    vector3fs: Vec<ParamSetItem<Vector3f>>,
    normals: Vec<ParamSetItem<Normal3f>>,
    // TODO: std::vector<std::shared_ptr<ParamSetItem<Spectrum>>> spectra;
    strings: Vec<ParamSetItem<String>>,
    textures: Vec<ParamSetItem<String>>,
    // TODO: static std::map<std::string, Spectrum> cachedSpectra;
}

// see api.cpp

pub struct RenderOptions {
    transform_start_time: Float,
    transform_end_time: Float,
    filter_name: String,
    filter_params: ParamSet,
    film_name: String, // "box"
    film_params: ParamSet,
    sampler_name: String, // "halton";
    sampler_params: ParamSet,
    accelerator_name: String, // "bvh";
    accelerator_params: ParamSet,
    integrator_name: String, // "path";
    integrator_params: ParamSet,
    camera_name: String, // "perspective";
    camera_params: ParamSet,
    // TODO: TransformSet camera_to_world;
    // TODO: std::map<std::string, std::shared_ptr<Medium>> namedMedia;
    // TODO: std::vector<std::shared_ptr<Light>> lights;
    // TODO: std::vector<std::shared_ptr<Primitive>> primitives;
    // TODO: std::map<std::string, std::vector<std::shared_ptr<Primitive>>> instances;
    // TODO: std::vector<std::shared_ptr<Primitive>> *currentInstance = nullptr;
    have_scattering_media: bool, // false
}
