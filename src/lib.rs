//! # pbrt
//!
//! [Rust][rust] crate to implement at least parts of the [PBRT
//! book][book]'s C++ code. You can find a copy of the current code
//! [here][repo].
//!
//! [rust]: https://www.rust-lang.org/en-US
//! [book]: http://www.pbrt.org
//! [repo]: https://github.com/wahn/rs_pbrt
//!
//! ## Scene
//!
//! As the scene file is parsed, objects are created that represent
//! the lights and geometric primitives in the scene. These are all
//! stored in the **Scene** object.
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
//!         t_max: std::f64::INFINITY,
//!         time: 0.0,
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
//! ## Surface Interaction
//!
//! The geometry of a particular point on a surface is represented by
//! a **SurfaceInteraction**. Having this abstraction lets most of the
//! system work with points on surfaces without needing to consider
//! the particular type of geometric shape the points lie on; the
//! **SurfaceInteraction** abstraction supplies enough information
//! about the surface point to allow the shading and geometric
//! operations in the rest of **pbrt** to be implemented generically.
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
//! use pbrt::{Sphere, Transform, Vector3f};
//!
//! fn main() {
//!     let default_sphere: Sphere = Sphere::default();
//!     let translate: Transform = Transform::translate(Vector3f {
//!         x: -1.3,
//!         y: 0.0,
//!         z: 0.0,
//!     });
//!     let inverse: Transform = Transform::inverse(translate);
//!     let sphere: Sphere = Sphere::new(translate, inverse, false, false, 2.0, -0.5, 0.75, 270.0);
//!
//!     println!("default sphere = {:?}", default_sphere);
//!     println!("sphere = {:?}", sphere);
//! }
//! ```
//!
//! ### Triangle Meshes
//!
//! While a natural representation would be to have a **Triangle**
//! shape implementation where each triangle stored the positions of
//! its three vertices, a more memory-efficient representation is to
//! separately store entire triangle meshes with an array of vertex
//! positions where each individual triangle just stores three offsets
//! into this array for its three vertices.
//!
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::{Point2f, Point3f, Transform, TriangleMesh, Vector3f};
//!
//! fn main() {
//!     let translate: Transform = Transform::translate(Vector3f {
//!         x: 0.25,
//!         y: 0.0,
//!         z: 0.0,
//!     });
//!     let inverse: Transform = Transform::inverse(translate);
//!     let n_triangles: usize = 2;
//!     let vertex_indices: Vec<usize> = vec![0_usize, 2, 1, 0, 3, 2];
//!     let n_vertices: usize = 4;
//!     let p: Vec<Point3f> = vec![Point3f {
//!                                    x: -100.0,
//!                                    y: -1.0,
//!                                    z: -100.0,
//!                                },
//!                                Point3f {
//!                                    x: 400.0,
//!                                    y: -1.0,
//!                                    z: -100.0,
//!                                },
//!                                Point3f {
//!                                    x: 400.0,
//!                                    y: -1.0,
//!                                    z: 400.0,
//!                                },
//!                                Point3f {
//!                                    x: -100.0,
//!                                    y: -1.0,
//!                                    z: 400.0,
//!                                }];
//!     let s: Vec<Vector3f> = Vec::new();
//!     let n: Vec<Vector3f> = Vec::new();
//!     let uv: Vec<Point2f> = vec![Point2f { x: 0.0, y: 0.0 },
//!                                 Point2f { x: 1.0, y: 0.0 },
//!                                 Point2f { x: 0.0, y: 1.0 },
//!                                 Point2f { x: 1.0, y: 1.0 }];
//!     let triangle_mesh: TriangleMesh = TriangleMesh::new(translate,
//!                                                         inverse,
//!                                                         false,
//!                                                         false,
//!                                                         n_triangles,
//!                                                         vertex_indices,
//!                                                         n_vertices,
//!                                                         p,
//!                                                         s,
//!                                                         n,
//!                                                         uv);
//!
//!     println!("translate = {:?}", translate);
//!     println!("inverse = {:?}", inverse);
//!     println!("triangle_mesh = {:?}", triangle_mesh);
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
//! use pbrt::{Bounds2f, BoxFilter, Film, Point2f, Point2i, Vector2f};
//! use std::string::String;
//!
//! fn main() {
//!     // see film.cpp CreateFilm()
//!     let film: Film = Film::new(Point2i { x: 1280, y: 720 },
//!                                Bounds2f {
//!                                    p_min: Point2f { x: 0.0, y: 0.0 },
//!                                    p_max: Point2f { x: 1.0, y: 1.0 },
//!                                },
//!                                BoxFilter {
//!                                    radius: Vector2f { x: 0.5, y: 0.5 },
//!                                    inv_radius: Vector2f {
//!                                        x: 1.0 / 0.5,
//!                                        y: 1.0 / 0.5,
//!                                    },
//!                                },
//!                                35.0,
//!                                String::from("pbrt.exr"),
//!                                1.0,
//!                                std::f64::INFINITY);
//!
//!     println!("film = {:?}", film);
//! }
//! ```
//!
//! ## Light Sources
//!
//! ### Point Lights
//!
//! TODO
//!
//! ### Distant Lights
//!
//! A distant light, also known as directional light, describes an
//! emitter that deposits illumination from the same direction at
//! every point in space.
//!
//! ```rust
//! extern crate pbrt;
//!
//! use pbrt::{DistantLight, Point3f, Spectrum, Transform, Vector3f};
//!
//! fn main() {
//!     let l: Spectrum = Spectrum::new(3.141593);
//!     let sc: Spectrum = Spectrum::new(1.0);
//!     let from: Point3f = Point3f {
//!         x: 0.0,
//!         y: 10.0,
//!         z: 0.0,
//!     };
//!     let to: Point3f = Point3f {
//!         x: 0.0,
//!         y: 0.0,
//!         z: 0.0,
//!     };
//!     let dir: Vector3f = from - to;
//!     let light_to_world: Transform = Transform::default();
//!     let lsc: Spectrum = l * sc;
//!     let distant_light: DistantLight = DistantLight::new(&light_to_world, &lsc, &dir);
//!     println!("distant_light = {:?}", distant_light);
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

use std::cmp::PartialEq;
use std::ops::{Add, AddAssign, Sub, Mul, MulAssign, Div, DivAssign, Neg, Index};
use std::default::Default;
use std::f64::consts::PI;
use std::mem;
use std::sync::Arc;

pub type Float = f64;

// see scene.h

#[derive(Clone)]
pub struct Scene<'a> {
    pub lights: Vec<DistantLight>, // TODO: Light
    pub infinite_lights: Vec<DistantLight>, // TODO: Light
    aggregate: &'a BVHAccel<'a>, // TODO: Primitive,
    world_bound: Bounds3f,
}

impl<'a> Scene<'a> {
    pub fn new(aggregate: &'a BVHAccel<'a>, lights: Vec<DistantLight>) -> Self {
        let world_bound: Bounds3f = aggregate.world_bound();
        let scene: Scene = Scene {
            lights: Vec::new(),
            infinite_lights: Vec::new(),
            aggregate: aggregate,
            world_bound: world_bound,
        };
        let mut changed_lights = Vec::new();
        let mut infinite_lights = Vec::new();
        for mut light in lights {
            light.preprocess(&scene);
            changed_lights.push(light);
            let check: i32 = light.flags & LightFlags::Infinite as i32;
            if check == LightFlags::Infinite as i32 {
                infinite_lights.push(light);
            }
        }
        Scene {
            lights: changed_lights,
            infinite_lights: infinite_lights,
            aggregate: aggregate,
            world_bound: world_bound,
        }
    }
    pub fn world_bound(&self) -> &Bounds3f {
        &self.world_bound
    }
}

// see pbrt.h

pub type Spectrum = RGBSpectrum;

const MACHINE_EPSILON: Float = std::f64::EPSILON * 0.5;

/// Use **unsafe**
/// [std::mem::transmute_copy][transmute_copy]
/// to convert *f32* to *u32*.
///
/// [transmute_copy]: https://doc.rust-lang.org/std/mem/fn.transmute_copy.html
pub fn float_to_bits(f: f32) -> u32 {
    // uint64_t ui;
    // memcpy(&ui, &f, sizeof(double));
    // return ui;
    let rui: u32;
    unsafe {
        let ui: u32 = std::mem::transmute_copy(&f);
        rui = ui;
    }
    rui
}

/// Use **unsafe**
/// [std::mem::transmute_copy][transmute_copy]
/// to convert *u32* to *f32*.
///
/// [transmute_copy]: https://doc.rust-lang.org/std/mem/fn.transmute_copy.html
pub fn bits_to_float(ui: u32) -> f32 {
    // float f;
    // memcpy(&f, &ui, sizeof(uint32_t));
    // return f;
    let rf: f32;
    unsafe {
        let f: f32 = std::mem::transmute_copy(&ui);
        rf = f;
    }
    rf
}

/// Bump a floating-point value up to the next greater representable
/// floating-point value.
pub fn next_float_up(v: f32) -> f32 {
    if v.is_infinite() && v > 0.0 {
        v
    } else {
        let new_v: f32;
        if v == -0.0 {
            new_v = 0.0;
        } else {
            new_v = v;
        }
        let mut ui: u32 = float_to_bits(v);
        if new_v > 0.0 {
            ui += 1;
        } else {
            ui -= 1;
        }
        bits_to_float(ui)
    }
}

/// Bump a floating-point value down to the next smaller representable
/// floating-point value.
pub fn next_float_down(v: f32) -> f32 {
    if v.is_infinite() && v < 0.0 {
        v
    } else {
        let new_v: f32;
        if v == 0.0 {
            new_v = -0.0;
        } else {
            new_v = v;
        }
        let mut ui: u32 = float_to_bits(v);
        if new_v > 0.0 {
            ui -= 1;
        } else {
            ui += 1;
        }
        bits_to_float(ui)
    }
}

/// Error propagation.
pub fn gamma(n: i32) -> Float {
    (n as Float * MACHINE_EPSILON) / (1.0 - n as Float * MACHINE_EPSILON)
}


/// Clamp the given value *val* to lie between the values *low* and *high*.
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

/// Convert from angles expressed in degrees to radians.
pub fn radians(deg: Float) -> Float {
    (PI / 180.0) * deg
}

/// Convert from angles expressed in radians to degrees.
pub fn degrees(rad: Float) -> Float {
    (180.0 / PI) * rad
}

/// Round an integer up to the next higher (or equal) power of 2.
pub fn round_up_pow2_32(v: &mut i32) -> i32 {
    *v -= 1_i32;
    *v |= *v >> 1;
    *v |= *v >> 2;
    *v |= *v >> 4;
    *v |= *v >> 8;
    *v |= *v >> 16;
    *v + 1
}

/// Round an integer up to the next higher (or equal) power of 2.
pub fn round_up_pow2_64(v: &mut i64) -> i64 {
    *v -= 1_i64;
    *v |= *v >> 1;
    *v |= *v >> 2;
    *v |= *v >> 4;
    *v |= *v >> 8;
    *v |= *v >> 16;
    *v + 1
}

/// Find solution(s) of the quadratic equation at^2 + bt + c = 0.
pub fn quadratic(a: Float, b: Float, c: Float, t0: &mut Float, t1: &mut Float) -> bool {
    // find quadratic discriminant
    let discrim: f64 = (b as f64) * (b as f64) - 4.0 * (a as f64) * (c as f64);
    if discrim < 0.0 {
        false
    } else {
        let root_discrim: f64 = discrim.sqrt();
        // compute quadratic _t_ values
        let q: f64;
        if b < 0.0 {
            q = -0.5 * (b - root_discrim);
        } else {
            q = -0.5 * (b + root_discrim);
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

// see efloat.h

/// Find solution(s) of the quadratic equation at^2 + bt + c = 0 using
/// *EFloat* instead of *Float* for error bounds.
pub fn quadratic_efloat(a: EFloat, b: EFloat, c: EFloat, t0: &mut EFloat, t1: &mut EFloat) -> bool {
    let discrim: f64 = b.v as f64 * b.v as f64 - 4.0f64 * a.v as f64 * c.v as f64;
    if discrim < 0.0 {
        false
    } else {
        let root_discrim: f64 = discrim.sqrt();
        let float_root_discrim: EFloat = EFloat::new(root_discrim as f32,
                                                     MACHINE_EPSILON as f32 * root_discrim as f32);
        // compute quadratic _t_ values
        let q: EFloat;
        if b.v < 0.0f32 {
            q = (b - float_root_discrim) * -0.5f32;
        } else {
            q = (b + float_root_discrim) * -0.5f32;
        }
        *t0 = q / a;
        *t1 = c / q;
        if (*t0).v > (*t1).v {
            let swap: EFloat = *t0;
            *t0 = *t1;
            *t1 = swap;
        }
        true
    }
}

#[derive(Debug,Default,Copy,Clone)]
pub struct EFloat {
    v: f32,
    low: f32,
    high: f32,
}

impl EFloat {
    pub fn new(v: f32, err: f32) -> Self {
        if err == 0.0 {
            EFloat {
                v: v,
                low: v,
                high: v,
            }
        } else {
            EFloat {
                v: v,
                low: next_float_down(v - err),
                high: next_float_up(v + err),
            }
        }
    }
    pub fn lower_bound(&self) -> f32 {
        self.low
    }
    pub fn upper_bound(&self) -> f32 {
        self.high
    }
}

impl PartialEq for EFloat {
    fn eq(&self, rhs: &EFloat) -> bool {
        self.v == rhs.v
    }
}

impl Add for EFloat {
    type Output = EFloat;
    fn add(self, rhs: EFloat) -> EFloat {
        // TODO: r.Check();
        EFloat {
            v: self.v + rhs.v,
            low: next_float_down(self.lower_bound() + rhs.lower_bound()),
            high: next_float_up(self.upper_bound() + rhs.upper_bound()),
        }
    }
}

impl Sub for EFloat {
    type Output = EFloat;
    fn sub(self, rhs: EFloat) -> EFloat {
        // TODO: r.Check();
        EFloat {
            v: self.v - rhs.v,
            low: next_float_down(self.lower_bound() - rhs.upper_bound()),
            high: next_float_up(self.upper_bound() - rhs.lower_bound()),
        }
    }
}

impl Mul for EFloat {
    type Output = EFloat;
    fn mul(self, rhs: EFloat) -> EFloat {
        let prod: [f32; 4] = [self.lower_bound() * rhs.lower_bound(),
                              self.upper_bound() * rhs.lower_bound(),
                              self.lower_bound() * rhs.upper_bound(),
                              self.upper_bound() * rhs.upper_bound()];
        // TODO: r.Check();
        EFloat {
            v: self.v * rhs.v,
            low: next_float_down(prod[0].min(prod[1]).min(prod[2].min(prod[3]))),
            high: next_float_up(prod[0].max(prod[1]).max(prod[2].max(prod[3]))),
        }
    }
}

impl Mul<f32> for EFloat {
    type Output = EFloat;
    fn mul(self, rhs: f32) -> EFloat {
        EFloat::new(rhs, 0.0) * self
    }
}

impl Div for EFloat {
    type Output = EFloat;
    fn div(self, rhs: EFloat) -> EFloat {
        let div: [f32; 4] = [self.lower_bound() / rhs.lower_bound(),
                             self.upper_bound() / rhs.lower_bound(),
                             self.lower_bound() / rhs.upper_bound(),
                             self.upper_bound() / rhs.upper_bound()];
        // TODO: r.Check();
        if rhs.low < 0.0 && rhs.high > 0.0 {
            // the interval we're dividing by straddles zero, so just
            // return an interval of everything
            EFloat {
                v: self.v / rhs.v,
                low: -std::f32::INFINITY,
                high: std::f32::INFINITY,
            }
        } else {
            EFloat {
                v: self.v / rhs.v,
                low: next_float_down(div[0].min(div[1]).min(div[2].min(div[3]))),
                high: next_float_up(div[0].max(div[1]).max(div[2].max(div[3]))),
            }
        }
    }
}

// parallel.h

#[derive(Debug,Default,Copy,Clone)]
struct AtomicFloat {
    bits: u32,
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

/// Product of the Euclidean magnitudes of the two vectors and the
/// cosine of the angle between them. A return value of zero means
/// both vectors are orthogonal, a value if one means they are
/// codirectional.
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

impl<T> Add for Vector3<T>
    where T: Copy + Add<T, Output = T>
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
    pub fn abs(&self) -> Vector3<T>
        where T: num::Float
    {
        Vector3::<T> {
            x: self.x.abs(),
            y: self.y.abs(),
            z: self.z.abs(),
        }
    }
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

impl<T> Index<u8> for Vector3<T>
{
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

/// Product of the Euclidean magnitudes of the two vectors and the
/// cosine of the angle between them. A return value of zero means
/// both vectors are orthogonal, a value if one means they are
/// codirectional.
pub fn vec3_dot<T>(v1: Vector3<T>, v2: Vector3<T>) -> T
    where T: Copy + Add<T, Output = T> + Mul<T, Output = T>
{
    v1.x * v2.x + v1.y * v2.y + v1.z * v2.z
}

/// Given two vectors in 3D, the cross product is a vector that is
/// perpendicular to both of them.
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

/// Compute a new vector pointing in the same direction but with unit
/// length.
pub fn vec3_normalize<T>(v: Vector3<T>) -> Vector3<T>
    where T: num::Float + Copy + Div<T, Output = T>
{
    v / v.length()
}

/// Return the largest coordinate value.
pub fn vec3_max_component<T>(v: Vector3<T>) -> T
    where T: num::Float
{
    v.x.max(v.y.max(v.z))
}

/// Return the index of the component with the largest value.
pub fn vec3_max_dimension<T>(v: Vector3<T>) -> usize
    where T: std::cmp::PartialOrd
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
pub fn vec3_permute<T>(v: Vector3<T>, x: usize, y: usize, z: usize) -> Vector3<T>
    where T: Copy
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
pub fn vec3_coordinate_system<T>(v1: &Vector3<T>, v2: &mut Vector3<T>, v3: &mut Vector3<T>)
    where T: num::Float
{
    if v1.x.abs() > v1.y.abs() {
        *v2 = Vector3::<T> {
            x: -v1.z,
            y: num::Zero::zero(),
            z: v1.x,
        } / (v1.x * v1.x + v1.z * v1.z).sqrt();
    } else {
        *v2 = Vector3::<T> {
            x: num::Zero::zero(),
            y: v1.z,
            z: -v1.y,
        } / (v1.y * v1.y + v1.z * v1.z).sqrt();
    }
    *v3 = vec3_cross(*v1, *v2);
}

#[derive(Debug,Default,Copy,Clone)]
pub struct Point2<T> {
    pub x: T,
    pub y: T,
}

impl<T> Add<Point2<T>> for Point2<T>
    where T: Add<T, Output = T>
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
    where T: Add<T, Output = T>
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
    where T: Sub<T, Output = T>
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
    where T: Sub<T, Output = T>
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
    where T: Copy + Mul<T, Output = T>
{
    type Output = Point2<T>;
    fn mul(self, rhs: T) -> Point2<T>
        where T: Copy + Mul<T, Output = T>
    {
        Point2::<T> {
            x: self.x * rhs,
            y: self.y * rhs,
        }
    }
}

/// Apply floor operation component-wise.
pub fn pnt2_floor<T>(p: Point2<T>) -> Point2<T>
    where T: num::Float
{
    Point2 {
        x: p.x.floor(),
        y: p.y.floor(),
    }
}

/// Apply ceil operation component-wise.
pub fn pnt2_ceil<T>(p: Point2<T>) -> Point2<T>
    where T: num::Float
{
    Point2 {
        x: p.x.ceil(),
        y: p.y.ceil(),
    }
}

#[derive(Debug,Default,Copy,Clone)]
pub struct Point3<T> {
    pub x: T,
    pub y: T,
    pub z: T,
}

impl<T> Add<Point3<T>> for Point3<T>
    where T: Add<T, Output = T>
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
    where T: Add<T, Output = T>
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
    where T: AddAssign
{
    fn add_assign(&mut self, rhs: Vector3<T>) {
        self.x += rhs.x;
        self.y += rhs.y;
        self.z += rhs.z;
    }
}

impl<T> Sub<Point3<T>> for Point3<T>
    where T: Sub<T, Output = T>
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
    where T: Sub<T, Output = T>
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
    where T: Copy + Mul<T, Output = T>
{
    type Output = Point3<T>;
    fn mul(self, rhs: T) -> Point3<T>
        where T: Copy + Mul<T, Output = T>
    {
        Point3::<T> {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl<T> MulAssign<T> for Point3<T>
    where T: Copy + MulAssign
{
    fn mul_assign(&mut self, rhs: T) {
        self.x *= rhs;
        self.y *= rhs;
        self.z *= rhs;
    }
}

impl<T> Div<T> for Point3<T>
    where T: Copy + Div<T, Output = T>
{
    type Output = Point3<T>;
    fn div(self, rhs: T) -> Point3<T>
        where T: Copy + Div<T, Output = T>
    {
        Point3::<T> {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

impl<T> DivAssign<T> for Point3<T>
    where T: Copy + DivAssign
{
    fn div_assign(&mut self, rhs: T) {
        self.x /= rhs;
        self.y /= rhs;
        self.z /= rhs;
    }
}

impl<T> Index<u8> for Point3<T>
{
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

/// Permute the coordinate values according to the povided
/// permutation.
pub fn pnt3_permute<T>(v: Point3<T>, x: usize, y: usize, z: usize) -> Point3<T>
    where T: Copy
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

/// Apply floor operation component-wise.
pub fn pnt3_floor<T>(p: Point3<T>) -> Point3<T>
    where T: num::Float
{
    Point3 {
        x: p.x.floor(),
        y: p.y.floor(),
        z: p.z.floor(),
    }
}

/// Apply ceil operation component-wise.
pub fn pnt3_ceil<T>(p: Point3<T>) -> Point3<T>
    where T: num::Float
{
    Point3 {
        x: p.x.ceil(),
        y: p.y.ceil(),
        z: p.z.ceil(),
    }
}

/// The distance between two points is the length of the vector
/// between them.
pub fn pnt3_distance<T>(p1: Point3<T>, p2: Point3<T>) -> T
    where T: num::Float + Sub<T, Output = T>
{
    (p1 - p2).length()
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

impl<T> Bounds2<T> {
    pub fn new(p1: Point2<T>, p2: Point2<T>) -> Self
        where T: num::Float
    {
        let p_min: Point2<T> = Point2::<T> {
            x: p1.x.min(p2.x),
            y: p1.y.min(p2.y),
        };
        let p_max: Point2<T> = Point2::<T> {
            x: p1.x.max(p2.x),
            y: p1.y.max(p2.y),
        };
        Bounds2::<T> {
            p_min: p_min,
            p_max: p_max,
        }
    }
}

#[derive(Debug,Copy,Clone)]
pub struct Bounds3<T> {
    pub p_min: Point3<T>,
    pub p_max: Point3<T>,
}

// work around bug
// https://github.com/rust-lang/rust/issues/40395
impl Default for /* Bounds3f */ Bounds3<f64> {
    fn default() -> /* Bounds3f */ Bounds3<f64> {
        let min_num: Float = std::f64::MIN;
        let max_num: Float = std::f64::MAX;
        /* Bounds3f */ Bounds3::<f64> {
            p_min: Point3f {
                x: max_num,
                y: max_num,
                z: max_num,
            }, p_max: Point3f {
                x: min_num,
                y: min_num,
                z: min_num,
            },
        }
    }
}

impl<T> Bounds3<T> {
    pub fn new(p1: Point3<T>, p2: Point3<T>) -> Self
        where T: num::Float
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
    pub fn diagonal(&self) -> Vector3<T>
        where T: Copy + Sub<T, Output = T>
    {
        self.p_max - self.p_min
    }
    pub fn surface_area(&self) -> T
        where T: Copy + Add<T, Output = T> + Sub<T, Output = T> + Mul<T, Output = T>
    {
        let d: Vector3<T> = self.diagonal();
        // 2 * (d.x * d.y + d.x * d.z + d.y * d.z)
        let r: T = d.x * d.y + d.x * d.z + d.y * d.z;
        r + r // avoid '2 *'
    }
    pub fn maximum_extent(&self) -> u8
        where T: Copy + std::cmp::PartialOrd + Sub<T, Output = T>
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
    pub fn offset(&self, p: Point3<T>) -> Vector3<T>
        where T: Copy + std::cmp::PartialOrd + Sub<T, Output = T> + DivAssign<T>
    {
        let mut o: Vector3<T> = p - self.p_min;
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
        let is_inside: bool = pnt3_inside_bnd3(center_copy, *b);
        if is_inside {
            *radius = pnt3_distance(center_copy, p_max);
        } else {
            *radius = 0.0;
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
pub fn bnd3_union_pnt3<T>(b: Bounds3<T>, p: Point3<T>) -> Bounds3<T>
    where T: num::Float
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
pub fn bnd3_union_bnd3<T>(b1: Bounds3<T>, b2: Bounds3<T>) -> Bounds3<T>
    where T: num::Float
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
pub fn pnt3_inside_bnd3(p: Point3f, b: Bounds3f) -> bool {
    p.x >= b.p_min.x && p.x <= b.p_max.x && p.y >= b.p_min.y && p.y <= b.p_max.y &&
    p.z >= b.p_min.z && p.z <= b.p_max.z
}

#[derive(Debug,Default,Copy,Clone)]
pub struct Ray {
    /// origin
    pub o: Point3f,
    /// direction
    pub d: Vector3f,
    /// limits the ray to a segment along its infinite extent
    pub t_max: Float,
    /// used for animations
    pub time: Float,
}

impl Ray {
    // Point3f operator()(Float t) const { return o + d * t; }
    fn position(&self, t: Float) -> Point3f {
        self.o + self.d * t
    }
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
               -> Self {
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

/// The product of two matrices.
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
               -> Self {
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
    pub fn transform_point(&self, p: Point3<Float>) -> Point3<Float> {
        let x: Float = p.x;
        let y: Float = p.y;
        let z: Float = p.z;
        let xp: Float = self.m.m[0][0] * x + self.m.m[0][1] * y + self.m.m[0][2] * z +
                        self.m.m[0][3];
        let yp: Float = self.m.m[1][0] * x + self.m.m[1][1] * y + self.m.m[1][2] * z +
                        self.m.m[1][3];
        let zp: Float = self.m.m[2][0] * x + self.m.m[2][1] * y + self.m.m[2][2] * z +
                        self.m.m[2][3];
        let wp: Float = self.m.m[3][0] * x + self.m.m[3][1] * y + self.m.m[3][2] * z +
                        self.m.m[3][3];
        assert!(wp != 0.0, "wp = {:?} != 0.0", wp);
        if wp == 1. {
            Point3::<Float> {
                x: xp,
                y: yp,
                z: zp,
            }
        } else {
            Point3::<Float> {
                x: xp / wp,
                y: yp / wp,
                z: zp / wp,
            }
        }
    }
    pub fn transform_vector(&self, v: Vector3<Float>) -> Vector3<Float> {
        let x: Float = v.x;
        let y: Float = v.y;
        let z: Float = v.z;
        Vector3::<Float> {
            x: self.m.m[0][0] * x + self.m.m[0][1] * y + self.m.m[0][2] * z,
            y: self.m.m[1][0] * x + self.m.m[1][1] * y + self.m.m[1][2] * z,
            z: self.m.m[2][0] * x + self.m.m[2][1] * y + self.m.m[2][2] * z,
        }
    }
    pub fn transform_bounds(&self, b: Bounds3f) -> Bounds3f {
        let m: Transform = *self;
        let p: Point3f = self.transform_point(Point3f {
            x: b.p_min.x,
            y: b.p_min.y,
            z: b.p_min.z,
        });
        let mut ret: Bounds3f = Bounds3f {
            p_min: p,
            p_max: p,
        };
        ret = bnd3_union_pnt3(ret,
                              m.transform_point(Point3f {
                                  x: b.p_max.x,
                                  y: b.p_min.y,
                                  z: b.p_min.z,
                              }));
        ret = bnd3_union_pnt3(ret,
                              m.transform_point(Point3f {
                                  x: b.p_min.x,
                                  y: b.p_max.y,
                                  z: b.p_min.z,
                              }));
        ret = bnd3_union_pnt3(ret,
                              m.transform_point(Point3f {
                                  x: b.p_min.x,
                                  y: b.p_min.y,
                                  z: b.p_max.z,
                              }));
        ret = bnd3_union_pnt3(ret,
                              m.transform_point(Point3f {
                                  x: b.p_min.x,
                                  y: b.p_max.y,
                                  z: b.p_max.z,
                              }));
        ret = bnd3_union_pnt3(ret,
                              m.transform_point(Point3f {
                                  x: b.p_max.x,
                                  y: b.p_max.y,
                                  z: b.p_min.z,
                              }));
        ret = bnd3_union_pnt3(ret,
                              m.transform_point(Point3f {
                                  x: b.p_max.x,
                                  y: b.p_min.y,
                                  z: b.p_max.z,
                              }));
        ret = bnd3_union_pnt3(ret,
                              m.transform_point(Point3f {
                                  x: b.p_max.x,
                                  y: b.p_max.y,
                                  z: b.p_max.z,
                              }));
        ret
    }
    pub fn transform_point_with_error(&self,
                                      p: Point3<Float>,
                                      p_error: &mut Vector3<Float>)
                                      -> Point3<Float> {
        let x: Float = p.x;
        let y: Float = p.y;
        let z: Float = p.z;
        // compute transformed coordinates from point _pt_
        let xp: Float = self.m.m[0][0] * x + self.m.m[0][1] * y + self.m.m[0][2] * z +
                        self.m.m[0][3];
        let yp: Float = self.m.m[1][0] * x + self.m.m[1][1] * y + self.m.m[1][2] * z +
                        self.m.m[1][3];
        let zp: Float = self.m.m[2][0] * x + self.m.m[2][1] * y + self.m.m[2][2] * z +
                        self.m.m[2][3];
        let wp: Float = self.m.m[3][0] * x + self.m.m[3][1] * y + self.m.m[3][2] * z +
                        self.m.m[3][3];
        // compute absolute error for transformed point
        let x_abs_sum: Float = self.m.m[0][0].abs() * x + self.m.m[0][1].abs() * y +
                               self.m.m[0][2].abs() * z +
                               self.m.m[0][3].abs();
        let y_abs_sum: Float = self.m.m[1][0].abs() * x + self.m.m[1][1].abs() * y +
                               self.m.m[1][2].abs() * z +
                               self.m.m[1][3].abs();
        let z_abs_sum: Float = self.m.m[2][0].abs() * x + self.m.m[2][1].abs() * y +
                               self.m.m[2][2].abs() * z +
                               self.m.m[2][3].abs();
        *p_error = Vector3::<Float> {
            x: x_abs_sum,
            y: y_abs_sum,
            z: z_abs_sum,
        } * gamma(3i32);
        assert!(wp != 0.0, "wp = {:?} != 0.0", wp);
        if wp == 1. {
            Point3::<Float> {
                x: xp,
                y: yp,
                z: zp,
            }
        } else {
            Point3::<Float> {
                x: xp / wp,
                y: yp / wp,
                z: zp / wp,
            }
        }
    }
    pub fn transform_vector_with_error(&self,
                                       v: Vector3<Float>,
                                       abs_error: &mut Vector3<Float>)
                                       -> Vector3<Float> {
        let x: Float = v.x;
        let y: Float = v.y;
        let z: Float = v.z;
        let gamma: Float = gamma(3i32);
        abs_error.x = gamma *
                      ((self.m.m[0][0] * v.x).abs() + (self.m.m[0][1] * v.y).abs() +
                       (self.m.m[0][2] * v.z).abs());
        abs_error.y = gamma *
                      ((self.m.m[1][0] * v.x).abs() + (self.m.m[1][1] * v.y).abs() +
                       (self.m.m[1][2] * v.z).abs());
        abs_error.z = gamma *
                      ((self.m.m[2][0] * v.x).abs() + (self.m.m[2][1] * v.y).abs() +
                       (self.m.m[2][2] * v.z).abs());
        Vector3::<Float> {
            x: self.m.m[0][0] * x + self.m.m[0][1] * y + self.m.m[0][2] * z,
            y: self.m.m[1][0] * x + self.m.m[1][1] * y + self.m.m[1][2] * z,
            z: self.m.m[2][0] * x + self.m.m[2][1] * y + self.m.m[2][2] * z,
        }
    }
    pub fn transform_ray_with_error(&self,
                                    r: &Ray,
                                    o_error: &mut Vector3<Float>,
                                    d_error: &mut Vector3<Float>)
                                    -> Ray {
        let mut o: Point3f = self.transform_point_with_error(r.o, o_error);
        let d: Vector3f = self.transform_vector_with_error(r.d, d_error);
        let length_squared: Float = d.length_squared();
        if length_squared > 0.0 {
            let dt: Float = vec3_dot(d.abs(), *o_error) / length_squared;
            o += d * dt;
        }
        Ray {
            o: o,
            d: d,
            t_max: r.t_max,
            time: 0.0,
        }
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
               -> Self {
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
}

/// The inner product of two quaterions.
pub fn quat_dot(q1: Quaternion, q2: Quaternion) -> Float {
    vec3_dot(q1.v, q2.v) + q1.w * q2.w
}

/// A quaternion can be normalized by dividing by its length.
pub fn quat_normalize(q: Quaternion) -> Quaternion {
    q / quat_dot(q, q)
}

// see interaction.h

#[derive(Debug,Default,Copy,Clone)]
pub struct SurfaceInteraction {
    // Interaction Public Data
    p: Point3f,
    time: Float,
    p_error: Vector3f,
    wo: Vector3f,
    n: Normal3f,
    // TODO: MediumInterface mediumInterface;
    // SurfaceInteraction Public Data
    uv: Point2f,
    dpdu: Vector3f,
    dpdv: Vector3f,
    dndu: Normal3f,
    dndv: Normal3f, /* const Shape *shape = nullptr;
                     * struct {
                     *     Normal3f n;
                     *     Vector3f dpdu, dpdv;
                     *     Normal3f dndu, dndv;
                     * } shading;
                     * const Primitive *primitive = nullptr;
                     * BSDF *bsdf = nullptr;
                     * BSSRDF *bssrdf = nullptr;
                     * mutable Vector3f dpdx, dpdy;
                     * mutable Float dudx = 0, dvdx = 0, dudy = 0, dvdy = 0; */
}

impl SurfaceInteraction {
    pub fn new(p: Point3f,
               p_error: Vector3f,
               uv: Point2f,
               wo: Vector3f,
               dpdu: Vector3f,
               dpdv: Vector3f,
               dndu: Normal3f,
               dndv: Normal3f,
               time: Float /* ,
                            * sh: &Shape */)
               -> Self {
        let nv: Vector3f = vec3_normalize(vec3_cross(dpdu, dpdv));
        let n: Normal3f = Normal3f {
            x: nv.x,
            y: nv.y,
            z: nv.z,
        };
        SurfaceInteraction {
            p: p,
            time: time,
            p_error: p_error,
            wo: wo,
            n: n,
            uv: uv,
            dpdu: dpdu,
            dpdv: dpdv,
            dndu: dndu,
            dndv: dndv,
        }
    }
}

// see primitive.h

pub trait Primitive {
    fn world_bound(&self) -> Bounds3f;
    fn intersect(&self,
                 r: &Ray,
                 t_hit: &mut Float,
                 isect: &mut SurfaceInteraction /* , bool testAlphaTexture */)
                 -> bool;
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
    // inherited from class Shape (see shape.h)
    object_to_world: Transform,
    world_to_object: Transform,
    reverse_orientation: bool,
    transform_swaps_handedness: bool,
}

impl Default for Sphere {
    fn default() -> Sphere {
        Sphere {
            // Shape
            object_to_world: Transform::default(),
            world_to_object: Transform::default(),
            reverse_orientation: false,
            transform_swaps_handedness: false,
            // Sphere
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
    pub fn new(object_to_world: Transform,
               world_to_object: Transform,
               reverse_orientation: bool,
               transform_swaps_handedness: bool,
               radius: Float,
               z_min: Float,
               z_max: Float,
               phi_max: Float)
               -> Self {
        Sphere {
            // Shape
            object_to_world: object_to_world,
            world_to_object: world_to_object,
            reverse_orientation: reverse_orientation,
            transform_swaps_handedness: transform_swaps_handedness,
            // Sphere
            radius: radius,
            z_min: clamp(z_min.min(z_max), -radius, radius),
            z_max: clamp(z_min.max(z_max), -radius, radius),
            theta_min: clamp(z_min.min(z_max) / radius, -1.0, 1.0).acos(),
            theta_max: clamp(z_min.max(z_max) / radius, -1.0, 1.0).acos(),
            phi_max: phi_max,
        }
    }
    pub fn object_bound(&self) -> Bounds3f {
        Bounds3f {
            p_min: Point3f {
                x: -self.radius,
                y: -self.radius,
                z: self.z_min,
            },
            p_max: Point3f {
                x: self.radius,
                y: self.radius,
                z: self.z_max,
            },
        }
    }
}

impl Primitive for Sphere {
    fn world_bound(&self) -> Bounds3f {
        // in C++: Bounds3f Shape::WorldBound() const { return (*ObjectToWorld)(ObjectBound()); }
        self.object_to_world.transform_bounds(self.object_bound())
    }
    fn intersect(&self,
                 r: &Ray,
                 t_hit: &mut Float,
                 isect: &mut SurfaceInteraction /* , bool testAlphaTexture */)
                 -> bool {
        // transform _Ray_ to object space
        let mut o_err: Vector3f = Vector3f::default();
        let mut d_err: Vector3f = Vector3f::default();
        let ray: Ray = self.world_to_object.transform_ray_with_error(r, &mut o_err, &mut d_err);

        // compute quadratic sphere coefficients

        // initialize _EFloat_ ray coordinate values
        let ox = EFloat::new(ray.o.x as f32, o_err.x as f32);
        let oy = EFloat::new(ray.o.y as f32, o_err.y as f32);
        let oz = EFloat::new(ray.o.z as f32, o_err.z as f32);
        let dx = EFloat::new(ray.d.x as f32, d_err.x as f32);
        let dy = EFloat::new(ray.d.y as f32, d_err.y as f32);
        let dz = EFloat::new(ray.d.z as f32, d_err.z as f32);
        let a: EFloat = dx * dx + dy * dy + dz * dz;
        let b: EFloat = (dx * ox + dy * oy + dz * oz) * 2.0f32;
        let c: EFloat = ox * ox + oy * oy + oz * oz -
                        EFloat::new(self.radius as f32, 0.0) * EFloat::new(self.radius as f32, 0.0);

        // solve quadratic equation for _t_ values
        let mut t0: EFloat = EFloat::default();
        let mut t1: EFloat = EFloat::default();
        if !quadratic_efloat(a, b, c, &mut t0, &mut t1) {
            return false;
        }
        // check quadric shape _t0_ and _t1_ for nearest intersection
        if t0.upper_bound() > ray.t_max as f32 || t1.lower_bound() <= 0.0f32 {
            return false;
        }
        let mut t_shape_hit: EFloat = t0;
        if t_shape_hit.lower_bound() <= 0.0f32 {
            t_shape_hit = t1;
            if t_shape_hit.upper_bound() > ray.t_max as f32 {
                return false;
            }
        }
        // compute sphere hit position and $\phi$
        let mut p_hit: Point3f = ray.position(t_shape_hit.v as f64);
        // refine sphere intersection point
        p_hit *= self.radius / pnt3_distance(p_hit, Point3f::default());
        if p_hit.x == 0.0 && p_hit.y == 0.0 {
            p_hit.x = 1e-5_f64 * self.radius;
        }
        let mut phi: Float = p_hit.y.atan2(p_hit.x);
        if phi < 0.0 {
            phi += 2.0_f64 * PI;
        }
        // test sphere intersection against clipping parameters
        if (self.z_min > -self.radius && p_hit.z < self.z_min) ||
           (self.z_max < self.radius && p_hit.z > self.z_max) || phi > self.phi_max {
            if t_shape_hit == t1 {
                return false;
            }
            if t1.upper_bound() > ray.t_max as f32 {
                return false;
            }
            t_shape_hit = t1;
            // compute sphere hit position and $\phi$
            p_hit = ray.position(t_shape_hit.v as f64);

            // refine sphere intersection point
            p_hit *= self.radius / pnt3_distance(p_hit, Point3f::default());
            if p_hit.x == 0.0 && p_hit.y == 0.0 {
                p_hit.x = 1e-5_f64 * self.radius;
            }
            phi = p_hit.y.atan2(p_hit.x);
            if phi < 0.0 {
                phi += 2.0_f64 * PI;
            }
            if (self.z_min > -self.radius && p_hit.z < self.z_min) ||
               (self.z_max < self.radius && p_hit.z > self.z_max) ||
               phi > self.phi_max {
                return false;
            }
        }
        // find parametric representation of sphere hit
        let u: Float = phi / self.phi_max;
        let theta: Float = clamp(p_hit.z / self.radius, -1.0, 1.0).acos();
        let v: Float = (theta - self.theta_min) / (self.theta_max - self.theta_min);
        // compute sphere $\dpdu$ and $\dpdv$
        let z_radius: Float = (p_hit.x * p_hit.x + p_hit.y * p_hit.y).sqrt();
        let inv_z_radius: Float = 1.0 / z_radius;
        let cos_phi: Float = p_hit.x * inv_z_radius;
        let sin_phi: Float = p_hit.y * inv_z_radius;
        let dpdu: Vector3f = Vector3f {
            x: -self.phi_max * p_hit.y,
            y: self.phi_max * p_hit.x,
            z: 0.0,
        };
        let dpdv: Vector3f = Vector3f {
            x: p_hit.z * cos_phi,
            y: p_hit.z * sin_phi,
            z: -self.radius * theta.sin(),
        } * (self.theta_max - self.theta_min);
        // compute sphere $\dndu$ and $\dndv$
        let d2_p_duu: Vector3f = Vector3f {
            x: p_hit.x,
            y: p_hit.y,
            z: 0.0,
        } * -self.phi_max * self.phi_max;
        let d2_p_duv: Vector3f = Vector3f {
            x: -sin_phi,
            y: cos_phi,
            z: 0.0,
        } * (self.theta_max - self.theta_min) * p_hit.z *
                                 self.phi_max;
        let d2_p_dvv: Vector3f = Vector3f {
            x: p_hit.x,
            y: p_hit.y,
            z: p_hit.z,
        } * -(self.theta_max - self.theta_min) *
                                 (self.theta_max - self.theta_min);
        // compute coefficients for fundamental forms
        let ec: Float = vec3_dot(dpdu, dpdu);
        let fc: Float = vec3_dot(dpdu, dpdv);
        let gc: Float = vec3_dot(dpdv, dpdv);
        let nc: Vector3f = vec3_normalize(vec3_cross(dpdu, dpdv));
        let el: Float = vec3_dot(nc, d2_p_duu);
        let fl: Float = vec3_dot(nc, d2_p_duv);
        let gl: Float = vec3_dot(nc, d2_p_dvv);
        // compute $\dndu$ and $\dndv$ from fundamental form coefficients
        let inv_egf2: Float = 1.0 / (ec * gc - fc * fc);
        let dndu = dpdu * (fl * fc - el * gc) * inv_egf2 + dpdv * (el * fc - fl * ec) * inv_egf2;
        let dndu = Normal3f {
            x: dndu.x,
            y: dndu.y,
            z: dndu.z,
        };
        let dndv = dpdu * (gl * fc - fl * gc) * inv_egf2 + dpdv * (fl * fc - gl * ec) * inv_egf2;
        let dndv = Normal3f {
            x: dndv.x,
            y: dndv.y,
            z: dndv.z,
        };
        // compute error bounds for sphere intersection
        let p_error: Vector3f = Vector3f {
                x: p_hit.x,
                y: p_hit.y,
                z: p_hit.z,
            }
            .abs() * gamma(5_i32);
        // initialize _SurfaceInteraction_ from parametric information
        // *isect = (*ObjectToWorld)(SurfaceInteraction(pHit, p_error, Point2f(u, v),
        //                                              -ray.d, dpdu, dpdv, dndu, dndv,
        //                                              ray.time, this));
        let uv_hit: Point2f = Point2f { x: u, y: v };
        let wo: Vector3f = -ray.d;
        let si: SurfaceInteraction =
            SurfaceInteraction::new(p_hit, p_error, uv_hit, wo, dpdu, dpdv, dndu, dndv, ray.time);
        println!("isect = {:?}", isect); // TMP
        isect.p = si.p;
        isect.time = si.time;
        isect.p_error = si.p_error;
        isect.wo = si.wo;
        isect.n = si.n;
        isect.uv = si.uv;
        println!("isect = {:?}", isect); // TMP
        *t_hit = t_shape_hit.v as f64;
        true
    }
}

// see triangle.h

#[derive(Debug,Clone)]
pub struct TriangleMesh {
    /// the total number of triangles in the mesh
    pub n_triangles: usize,
    /// vector of vertex indices
    pub vertex_indices: Vec<usize>,
    /// the total number of vertices in the mesh
    pub n_vertices: usize,
    /// vector of *n_vertices* vertex positions
    pub p: Vec<Point3f>,
    /// an optional vector of normal vectors (can be empty)
    pub n: Vec<Vector3f>,
    /// an optional vector of tangent vectors (can be empty)
    pub s: Vec<Vector3f>,
    /// an optional vector of paramtric (u, v) values (texture coordinates)
    pub uv: Vec<Point2f>,
    // TODO: std::shared_ptr<Texture<Float>> alphaMask, shadowAlphaMask;
    // inherited from class Shape (see shape.h)
    pub object_to_world: Transform, // TODO: not pub?
    pub world_to_object: Transform, // TODO: not pub?
    reverse_orientation: bool,
    pub transform_swaps_handedness: bool, // TODO: not pub?
}

impl TriangleMesh {
    pub fn new(object_to_world: Transform,
               world_to_object: Transform,
               reverse_orientation: bool,
               transform_swaps_handedness: bool,
               n_triangles: usize,
               vertex_indices: Vec<usize>,
               n_vertices: usize,
               p: Vec<Point3f>,
               s: Vec<Vector3f>,
               n: Vec<Vector3f>,
               uv: Vec<Point2f>)
               -> Self {
        TriangleMesh {
            // Shape
            object_to_world: object_to_world,
            world_to_object: world_to_object,
            reverse_orientation: reverse_orientation,
            transform_swaps_handedness: transform_swaps_handedness,
            // TriangleMesh
            n_triangles: n_triangles,
            vertex_indices: vertex_indices,
            n_vertices: n_vertices,
            p: p,
            n: n,
            s: s,
            uv: uv,
        }
    }
}

#[derive(Debug,Copy,Clone)]
pub struct Triangle<'a> {
    mesh: &'a TriangleMesh,
    id: usize,
    // inherited from class Shape (see shape.h)
    object_to_world: Transform,
    world_to_object: Transform,
    reverse_orientation: bool,
    transform_swaps_handedness: bool,
}

impl<'a> Triangle<'a> {
    pub fn new(object_to_world: Transform,
               world_to_object: Transform,
               reverse_orientation: bool,
               mesh: &'a TriangleMesh,
               tri_number: usize)
               -> Self {
        Triangle {
            mesh: mesh,
            id: tri_number,
            object_to_world: object_to_world,
            world_to_object: world_to_object,
            reverse_orientation: reverse_orientation,
            transform_swaps_handedness: false,
        }
    }
    pub fn object_bound(&self) -> Bounds3f {
        let p0: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 0]];
        let p1: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 1]];
        let p2: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 2]];
        bnd3_union_pnt3(Bounds3f::new(self.world_to_object.transform_point(p0),
                                      self.world_to_object.transform_point(p1)),
                        self.world_to_object.transform_point(p2))
    }
    pub fn get_uvs(&self) -> [Point2f; 3] {
        if self.mesh.uv.is_empty() {
            [Point2f { x: 0.0, y: 0.0 }, Point2f { x: 1.0, y: 0.0 }, Point2f { x: 1.0, y: 1.0 }]
        } else {
            [self.mesh.uv[self.mesh.vertex_indices[self.id * 3 + 0]],
             self.mesh.uv[self.mesh.vertex_indices[self.id * 3 + 1]],
             self.mesh.uv[self.mesh.vertex_indices[self.id * 3 + 2]]]
        }
    }
}

impl<'a> Primitive for Triangle<'a> {
    fn world_bound(&self) -> Bounds3f {
        let p0: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 0]];
        let p1: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 1]];
        let p2: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 2]];
        bnd3_union_pnt3(Bounds3f::new(p0, p1), p2)
    }
    fn intersect(&self,
                 ray: &Ray,
                 t_hit: &mut Float,
                 isect: &mut SurfaceInteraction /* , bool testAlphaTexture */)
                 -> bool {
        // get triangle vertices in _p0_, _p1_, and _p2_
        let p0: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 0]];
        let p1: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 1]];
        let p2: Point3f = self.mesh.p[self.mesh.vertex_indices[self.id * 3 + 2]];
        // translate vertices based on ray origin
        let mut p0t: Point3f = p0 -
                               Vector3f {
            x: ray.o.x,
            y: ray.o.y,
            z: ray.o.z,
        };
        let mut p1t: Point3f = p1 -
                               Vector3f {
            x: ray.o.x,
            y: ray.o.y,
            z: ray.o.z,
        };
        let mut p2t: Point3f = p2 -
                               Vector3f {
            x: ray.o.x,
            y: ray.o.y,
            z: ray.o.z,
        };
        // permute components of triangle vertices and ray direction
        let kz: usize = vec3_max_dimension(ray.d.abs());
        let mut kx: usize = kz + 1;
        if kx == 3 {
            kx = 0;
        }
        let mut ky: usize = kx + 1;
        if ky == 3 {
            ky = 0;
        }
        let d: Vector3f = vec3_permute(ray.d, kx, ky, kz);
        p0t = pnt3_permute(p0t, kx, ky, kz);
        p1t = pnt3_permute(p1t, kx, ky, kz);
        p2t = pnt3_permute(p2t, kx, ky, kz);
        // apply shear transformation to translated vertex positions
        let sx: Float = -d.x / d.z;
        let sy: Float = -d.y / d.z;
        let sz: Float = 1.0 / d.z;
        p0t.x += sx * p0t.z;
        p0t.y += sy * p0t.z;
        p1t.x += sx * p1t.z;
        p1t.y += sy * p1t.z;
        p2t.x += sx * p2t.z;
        p2t.y += sy * p2t.z;
        // compute edge function coefficients _e0_, _e1_, and _e2_
        let e0: Float = p1t.x * p2t.y - p1t.y * p2t.x;
        let e1: Float = p2t.x * p0t.y - p2t.y * p0t.x;
        let e2: Float = p0t.x * p1t.y - p0t.y * p1t.x;
        // TODO: fall back to double precision test at triangle edges
        if mem::size_of::<Float>() == mem::size_of::<f32>() {
            println!("[Triangle::intersect()]: TODO fall back to double precision test at \
                      triangle edges");
        }
        // if (sizeof(Float) == sizeof(float) &&
        //     (e0 == 0.0f || e1 == 0.0f || e2 == 0.0f)) {
        //     double p2txp1ty = (double)p2t.x * (double)p1t.y;
        //     double p2typ1tx = (double)p2t.y * (double)p1t.x;
        //     e0 = (float)(p2typ1tx - p2txp1ty);
        //     double p0txp2ty = (double)p0t.x * (double)p2t.y;
        //     double p0typ2tx = (double)p0t.y * (double)p2t.x;
        //     e1 = (float)(p0typ2tx - p0txp2ty);
        //     double p1txp0ty = (double)p1t.x * (double)p0t.y;
        //     double p1typ0tx = (double)p1t.y * (double)p0t.x;
        //     e2 = (float)(p1typ0tx - p1txp0ty);
        // }
        // perform triangle edge and determinant tests
        if (e0 < 0.0 || e1 < 0.0 || e2 < 0.0) && (e0 > 0.0 || e1 > 0.0 || e2 > 0.0) {
            return false;
        }
        let det: Float = e0 + e1 + e2;
        if det == 0.0 {
            return false;
        }
        // compute scaled hit distance to triangle and test against ray $t$ range
        p0t.z *= sz;
        p1t.z *= sz;
        p2t.z *= sz;
        let t_scaled: Float = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
        if det < 0.0 && (t_scaled >= 0.0 || t_scaled < ray.t_max * det) {
            return false;
        } else if det > 0.0 && (t_scaled <= 0.0 || t_scaled > ray.t_max * det) {
            return false;
        }
        // compute barycentric coordinates and $t$ value for triangle intersection
        let inv_det: Float = 1.0 / det;
        let b0: Float = e0 * inv_det;
        let b1: Float = e1 * inv_det;
        let b2: Float = e2 * inv_det;
        let t: Float = t_scaled * inv_det;

        // ensure that computed triangle $t$ is conservatively greater than zero

        // compute $\delta_z$ term for triangle $t$ error bounds
        let max_zt: Float = vec3_max_component(Vector3f {
                x: p0t.z,
                y: p1t.z,
                z: p2t.z,
            }
            .abs());
        let delta_z: Float = gamma(3_i32) * max_zt;
        // compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
        let max_xt: Float = vec3_max_component(Vector3f {
                x: p0t.x,
                y: p1t.x,
                z: p2t.x,
            }
            .abs());
        let max_yt: Float = vec3_max_component(Vector3f {
                x: p0t.y,
                y: p1t.y,
                z: p2t.y,
            }
            .abs());
        let delta_x: Float = gamma(5) * (max_xt + max_zt);
        let delta_y: Float = gamma(5) * (max_yt + max_zt);
        // compute $\delta_e$ term for triangle $t$ error bounds
        let delta_e: Float = 2.0 *
                             (gamma(2) * max_xt * max_yt + delta_y * max_xt + delta_x * max_yt);
        // compute $\delta_t$ term for triangle $t$ error bounds and check _t_
        let max_e: Float = vec3_max_component(Vector3f {
                x: e0,
                y: e1,
                z: e2,
            }
            .abs());
        let delta_t: Float =
            3.0 * (gamma(3) * max_e * max_zt + delta_e * max_zt + delta_z * max_e) * inv_det.abs();
        if t <= delta_t {
            return false;
        }
        // compute triangle partial derivatives
        let uv: [Point2f; 3] = self.get_uvs();
        // compute deltas for triangle partial derivatives
        let duv02: Vector2f = uv[0] - uv[2];
        let duv12: Vector2f = uv[1] - uv[2];
        let dp02: Vector3f = p0 - p2;
        let dp12: Vector3f = p1 - p2;
        let determinant: Float = duv02.x * duv12.y - duv02.y * duv12.x;
        let degenerate_uv: bool = determinant.abs() < 1e-8 as Float;
        // Vector3f dpdu, dpdv;
        let mut dpdu: Vector3f = Vector3f::default();
        let mut dpdv: Vector3f = Vector3f::default();
        if !degenerate_uv {
            let invdet: Float = 1.0 / determinant;
            dpdu = (dp02 * duv12.y - dp12 * duv02.y) * invdet;
            dpdv = (dp02 * -duv12.x + dp12 * duv02.x) * invdet;
        }
        if degenerate_uv || vec3_cross(dpdu, dpdv).length_squared() == 0.0 {
            // handle zero determinant for triangle partial derivative matrix
            vec3_coordinate_system(&vec3_normalize(vec3_cross(p2 - p0, p1 - p0)),
                                   &mut dpdu,
                                   &mut dpdv);
        }
        // compute error bounds for triangle intersection
        let x_abs_sum: Float = (b0 * p0.x).abs() + (b1 * p1.x).abs() + (b2 * p2.x).abs();
        let y_abs_sum: Float = (b0 * p0.y).abs() + (b1 * p1.y).abs() + (b2 * p2.y).abs();
        let z_abs_sum: Float = (b0 * p0.z).abs() + (b1 * p1.z).abs() + (b2 * p2.z).abs();
        let p_error: Vector3f = Vector3f {
            x: x_abs_sum,
            y: y_abs_sum,
            z: z_abs_sum,
        } * gamma(7);
        // interpolate $(u,v)$ parametric coordinates and hit point
        let p_hit: Point3f = p0 * b0 + p1 * b1 + p2 * b2;
        let uv_hit: Point2f = uv[0] * b0 + uv[1] * b1 + uv[2] * b2;
        // TODO: test intersection against alpha texture, if present
        // if (testAlphaTexture && mesh->alphaMask) {
        //     SurfaceInteraction isectLocal(p_hit, Vector3f(0, 0, 0), uv_hit, -ray.d,
        //                                   dpdu, dpdv, Normal3f(0, 0, 0),
        //                                   Normal3f(0, 0, 0), ray.time, this);
        //     if (mesh->alphaMask->Evaluate(isectLocal) == 0) return false;
        // }
        // fill in _SurfaceInteraction_ from triangle hit
        let dndu: Normal3f = Normal3f::default();
        let dndv: Normal3f = Normal3f::default();
        let wo: Vector3f = -ray.d;
        let si: SurfaceInteraction =
            SurfaceInteraction::new(p_hit, p_error, uv_hit, wo, dpdu, dpdv, dndu, dndv, ray.time);
        println!("isect = {:?}", isect); // TMP
        isect.p = si.p;
        isect.time = si.time;
        isect.p_error = si.p_error;
        isect.wo = si.wo;
        isect.n = si.n;
        isect.uv = si.uv;
        println!("isect = {:?}", isect); // TMP
        // WORK
        *t_hit = t;
        // TODO: ++nHits;
        true
    }
}

// see spectrum.h

#[derive(Debug,Default,Copy,Clone)]
pub struct RGBSpectrum {
    c: [Float; 3],
}

impl Mul for RGBSpectrum {
    type Output = RGBSpectrum;
    fn mul(self, rhs: RGBSpectrum) -> RGBSpectrum {
        RGBSpectrum { c: [self.c[0] * rhs.c[0], self.c[1] * rhs.c[1], self.c[2] * rhs.c[2]] }
    }
}

impl RGBSpectrum {
    pub fn new(v: Float) -> Self {
        let n_spectrum_samples = 3; // RGB
        RGBSpectrum { c: [v, v, v] }
        // TODO: DCHECK(!HasNaNs());
    }
}

// see bvh.h

#[derive(Debug,Clone)]
pub enum SplitMethod {
    SAH,
    HLBVH,
    Middle,
    EqualCounts,
}

#[derive(Debug,Default,Copy,Clone)]
pub struct BVHPrimitiveInfo {
    primitive_number: usize,
    bounds: Bounds3f,
    centroid: Point3f,
}

impl BVHPrimitiveInfo {
    pub fn new(primitive_number: usize, bounds: Bounds3f) -> Self {
        BVHPrimitiveInfo {
            primitive_number: primitive_number,
            bounds: bounds,
            centroid: bounds.p_min * 0.5 + bounds.p_max * 0.5,
        }
    }
}

#[derive(Debug,Clone)]
enum BVHLink {
    Empty,
    More(Box<BVHBuildNode>),
}

#[derive(Debug,Clone)]
pub struct BVHBuildNode {
    bounds: Bounds3f,
    child1: BVHLink,
    child2: BVHLink,
    split_axis: u8,
    first_prim_offset: usize,
    n_primitives: usize,
}

impl Default for BVHBuildNode {
    fn default() -> BVHBuildNode {
        BVHBuildNode {
            bounds: Bounds3f::default(),
            child1: BVHLink::Empty,
            child2: BVHLink::Empty,
            split_axis: 0_u8,
            first_prim_offset: 0_usize,
            n_primitives: 0_usize,
        }
    }
}

impl BVHBuildNode {
    pub fn init_leaf(&mut self, first: usize, n: usize, b: &Bounds3f) {
        self.first_prim_offset = first;
        self.n_primitives = n;
        self.bounds = *b;
        self.child1 = BVHLink::Empty;
        self.child2 = BVHLink::Empty;
    }
    pub fn init_interior(&mut self, axis: u8, c0: Box<BVHBuildNode>, c1: Box<BVHBuildNode>) {
        self.n_primitives = 0;
        self.bounds = bnd3_union_bnd3(c0.bounds, c1.bounds);
        self.child1 = BVHLink::More(c0);
        self.child2 = BVHLink::More(c1);
        self.split_axis = axis;
    }
}

#[derive(Debug,Copy,Clone)]
struct BucketInfo {
    count: usize,
    bounds: Bounds3f,
}

impl Default for BucketInfo {
    fn default() -> BucketInfo {
        BucketInfo {
            count: 0_usize,
            bounds: Bounds3f::default(),
        }
    }
}

#[derive(Debug,Default,Copy,Clone)]
pub struct LinearBVHNode {
    bounds: Bounds3f,
    // in C++ a union { int primitivesOffset;     // leaf
    //                  int secondChildOffset; }; // interior
    offset: usize,
    n_primitives: usize,
    axis: u8,
    // TODO? pad
}

// BVHAccel -> Aggregate -> Primitive
pub struct BVHAccel<'a> {
    max_prims_in_node: usize,
    split_method: SplitMethod,
    pub primitives: Vec<&'a Primitive>,
    pub nodes: Vec<LinearBVHNode>,
}

impl<'a> BVHAccel<'a> {
    pub fn new(p: Vec<&'a Primitive>, max_prims_in_node: usize, split_method: SplitMethod) -> Self {
        let bvh = Arc::new(BVHAccel {
            max_prims_in_node: std::cmp::min(max_prims_in_node, 255),
            split_method: split_method.clone(),
            primitives: p,
            nodes: Vec::new(),
        });
        let num_prims = bvh.primitives.len();
        let mut primitive_info = vec![BVHPrimitiveInfo::default(); num_prims];
        for i in 0..num_prims {
            let world_bound = bvh.primitives[i].world_bound();
            primitive_info[i] = BVHPrimitiveInfo::new(i, world_bound);
        }
        // TODO: if (splitMethod == SplitMethod::HLBVH)
        let mut total_nodes: usize = 0;
        let mut ordered_prims: Vec<&'a Primitive> = Vec::with_capacity(num_prims);
        let root = BVHAccel::recursive_build(bvh.clone(), // instead of self
                                             // arena,
                                             &mut primitive_info,
                                             0,
                                             num_prims,
                                             &mut total_nodes,
                                             &mut ordered_prims);
        // flatten first
        let mut nodes = vec![LinearBVHNode::default(); total_nodes];
        let mut offset: usize = 0;
        BVHAccel::flatten_bvh_tree(&root, &mut nodes, &mut offset);
        assert!(nodes.len() == total_nodes);
        // primitives.swap(orderedPrims);
        let bvh_ordered_prims = Arc::new(BVHAccel {
            max_prims_in_node: std::cmp::min(max_prims_in_node, 255),
            split_method: split_method.clone(),
            primitives: ordered_prims,
            nodes: nodes,
        });
        let unwrapped = Arc::try_unwrap(bvh_ordered_prims);
        unwrapped.ok().unwrap()
    }
    pub fn world_bound(&self) -> Bounds3f {
        if self.nodes.len() > 0 {
            self.nodes[0].bounds
        } else {
            Bounds3f::default()
        }
    }
    pub fn recursive_build(bvh: Arc<BVHAccel<'a>>,
                           // arena,
                           primitive_info: &mut Vec<BVHPrimitiveInfo>,
                           start: usize,
                           end: usize,
                           total_nodes: &mut usize,
                           ordered_prims: &mut Vec<&'a Primitive>)
                           -> Box<BVHBuildNode> {
        assert_ne!(start, end);
        let mut node: Box<BVHBuildNode> = Box::new(BVHBuildNode::default());
        *total_nodes += 1_usize;
        // compute bounds of all primitives in BVH node
        let mut bounds: Bounds3f = Bounds3f::default();
        for i in start..end {
            bounds = bnd3_union_bnd3(bounds, primitive_info[i].bounds);
        }
        let n_primitives: usize = end - start;
        if n_primitives == 1 {
            // create leaf _BVHBuildNode_
            let first_prim_offset: usize = ordered_prims.len();
            for i in start..end {
                let prim_num: usize = primitive_info[i].primitive_number;
                ordered_prims.push(bvh.primitives[prim_num]);
            }
            node.init_leaf(first_prim_offset, n_primitives, &bounds);
            return node;
        } else {
            // compute bound of primitive centroids, choose split dimension _dim_
            let mut centroid_bounds: Bounds3f = Bounds3f::default();
            for i in start..end {
                centroid_bounds = bnd3_union_pnt3(centroid_bounds, primitive_info[i].centroid);
            }
            let dim: u8 = centroid_bounds.maximum_extent();
            // partition primitives into two sets and build children
            let mut mid: usize = (start + end) / 2_usize;
            if centroid_bounds.p_max[dim] == centroid_bounds.p_min[dim] {
                // create leaf _BVHBuildNode_
                let first_prim_offset: usize = ordered_prims.len();
                for i in start..end {
                    let prim_num: usize = primitive_info[i].primitive_number;
                    ordered_prims.push(bvh.primitives[prim_num]);
                }
                node.init_leaf(first_prim_offset, n_primitives, &bounds);
                return node;
            } else {
                // partition primitives based on _splitMethod_
                match bvh.split_method {
                    SplitMethod::Middle => {
                        // TODO
                    }
                    SplitMethod::EqualCounts => {
                        // TODO
                    }
                    SplitMethod::SAH | SplitMethod::HLBVH => {
                        if n_primitives <= 2 {
                            mid = (start + end) / 2;
                            if start != end - 1 {
                                if primitive_info[end - 1].centroid[dim] <
                                    primitive_info[start].centroid[dim] {
                                        primitive_info.swap(start, end - 1);
                                }
                            }
                        } else {
                            // allocate _BucketInfo_ for SAH partition buckets
                            let n_buckets: usize = 12;
                            let mut buckets: [BucketInfo; 12] = [BucketInfo::default(); 12];
                            // initialize _BucketInfo_ for SAH partition buckets
                            for i in start..end {
                                let mut b: usize = n_buckets *
                                    centroid_bounds.offset(primitive_info[i].centroid)[dim]
                                    as usize;
                                if b == n_buckets {
                                    b = n_buckets - 1;
                                }
                                assert!(b >= 0_usize, "b >= 0");
                                assert!(b < n_buckets, "b < {}", n_buckets);
                                buckets[b].count += 1;
                                buckets[b].bounds = bnd3_union_bnd3(buckets[b].bounds,
                                                                    primitive_info[i].bounds);
                            }
                            // compute costs for splitting after each bucket
                            let mut cost: [Float; 11] = [0.0; 11];
                            for i in 0..(n_buckets - 1) {
                                let mut b0: Bounds3f = Bounds3f::default();
                                let mut b1: Bounds3f = Bounds3f::default();
                                let mut count0: usize = 0;
                                let mut count1: usize = 0;
                                for j in 0..(i + 1) {
                                    b0 = bnd3_union_bnd3(b0, buckets[j].bounds);
                                    count0 += buckets[j].count;
                                }
                                for j in (i + 1)..n_buckets {
                                    b1 = bnd3_union_bnd3(b1, buckets[j].bounds);
                                    count1 += buckets[j].count;
                                }
                                cost[i] = 1.0 +
                                          (count0 as Float * b0.surface_area() +
                                           count1 as Float * b1.surface_area()) /
                                          bounds.surface_area();
                            }
                            // find bucket to split at that minimizes SAH metric
                            let mut min_cost: Float = cost[0];
                            let mut min_cost_split_bucket: usize = 0;
                            for i in 0..(n_buckets - 1) {
                                if cost[i] < min_cost {
                                    min_cost = cost[i];
                                    min_cost_split_bucket = i;
                                }
                            }
                            // either create leaf or split primitives
                            // at selected SAH bucket
                            let leaf_cost: Float = n_primitives as Float;
                            if n_primitives > bvh.max_prims_in_node || min_cost < leaf_cost {
                                let (left, right): (Vec<BVHPrimitiveInfo>, Vec<BVHPrimitiveInfo>) =
                                    primitive_info[start..end].into_iter()
                                    .partition(|&pi| {
                                        let mut b: usize = n_buckets *
                                            centroid_bounds.offset(pi.centroid)[dim] as usize;
                                        if b == n_buckets {b = n_buckets - 1;}
                                        assert!(b >= 0_usize, "b >= 0");
                                        assert!(b < n_buckets, "b < {}", n_buckets);
                                        b <= min_cost_split_bucket
                                    });
                                mid = start + left.len();
                                let combined = [left, right].concat();
                                if combined.len() == primitive_info.len() {
                                    primitive_info.copy_from_slice(combined.as_slice());
                                } else {
                                    println!("TODO: combined.len() != primitive_info.len(); {} != {}",
                                    combined.len(), primitive_info.len());
                                }
                            } else {
                                // create leaf _BVHBuildNode_
                                let first_prim_offset: usize = ordered_prims.len();
                                for i in start..end {
                                    let prim_num: usize = primitive_info[i].primitive_number;
                                    ordered_prims.push(bvh.primitives[prim_num]);
                                }
                                node.init_leaf(first_prim_offset, n_primitives, &bounds);
                                return node;
                            }
                        }
                    }
                }
                // make sure we get result for c1 before c0
                let c1 = BVHAccel::recursive_build(bvh.clone(),
                                                   primitive_info,
                                                   mid,
                                                   end,
                                                   total_nodes,
                                                   ordered_prims);
                let c0 = BVHAccel::recursive_build(bvh.clone(),
                                                   primitive_info,
                                                   start,
                                                   mid,
                                                   total_nodes,
                                                   ordered_prims);
                node.init_interior(dim,
                                   c0,
                                   c1);
            }
        }
        return node;
    }
    fn flatten_bvh_tree(node: &Box<BVHBuildNode>, nodes: &mut Vec<LinearBVHNode>, offset: &mut usize) -> usize {
        let my_offset: usize = *offset;
        *offset += 1;
        if node.n_primitives > 0 {
            // leaf
            let linear_node = LinearBVHNode {
                bounds: node.bounds,
                offset: node.first_prim_offset,
                n_primitives: node.n_primitives,
                axis: 0_u8,
            };
            nodes[my_offset] = linear_node;
        } else {
            // interior
            let child1 = match node.clone().child1 {
                BVHLink::More(c) => c,
                _ => Box::new(BVHBuildNode::default()),
            };
            let child2 = match node.clone().child2 {
                BVHLink::More(c) => c,
                _ => Box::new(BVHBuildNode::default()),
            };
            BVHAccel::flatten_bvh_tree(&child1, nodes, offset);
            let linear_node = LinearBVHNode {
                bounds: node.bounds,
                offset: BVHAccel::flatten_bvh_tree(&child2, nodes, offset),
                n_primitives: 0_usize,
                axis: node.split_axis,
            };
            nodes[my_offset] = linear_node;
        }
        my_offset
    }
}

impl<'a> Primitive for BVHAccel<'a> {
    fn world_bound(&self) -> Bounds3f {
        // WORK
        Bounds3f::default()
    }
    fn intersect(&self,
                 r: &Ray,
                 t_hit: &mut Float,
                 isect: &mut SurfaceInteraction /* , bool testAlphaTexture */)
                 -> bool {
        // WORK
        false
    }
}

// see zerotwosequence.h

#[derive(Debug,Default,Clone)]
pub struct ZeroTwoSequenceSampler {
    pub samples_per_pixel: i64,
    pub n_sampled_dimensions: i64,
    // inherited from class Sampler (see sampler.h)
    pub samples_2d_array_sizes: Vec<i32>, // TODO: not pub?
    pub samples_2d_array: Vec<Point2f>, // TODO: not pub?
}

impl ZeroTwoSequenceSampler {
    pub fn round_count(&self, count: i32) -> i32 {
        let mut mut_count: i32 = count;
        round_up_pow2_32(&mut mut_count)
    }
    // sampler Interface
    pub fn request_2d_array(&mut self, n: i32) {
        assert_eq!(self.round_count(n), n);
        self.samples_2d_array_sizes.push(n);
        let size: usize = (n * self.samples_per_pixel as i32) as usize;
        let mut additional_points: Vec<Point2f> = vec![Point2f::default(); size];
        self.samples_2d_array.append(&mut additional_points);
    }
}

// see box.h

#[derive(Debug,Default,Copy,Clone)]
pub struct BoxFilter {
    pub radius: Vector2f,
    pub inv_radius: Vector2f,
}

// see film.h

#[derive(Debug,Default,Copy,Clone)]
pub struct Pixel {
    xyz: [Float; 3],
    filter_weight_sum: Float,
    splat_xyz: [AtomicFloat; 3],
    pad: Float,
}

#[derive(Debug,Default,Clone)]
pub struct Film {
    // Film Public Data
    /// The overall resolution of the image in pixels
    pub full_resolution: Point2i,
    /// The length of the diagonal of the film's physical area (specified in mm, stored in meters)
    pub diagonal: Float,
    /// A filter function
    pub filter: BoxFilter, // TODO: Filter
    /// The filename of the output image
    pub filename: String,
    /// A crop window that may specify a subset of the image to render
    pub cropped_pixel_bounds: Bounds2i,

    // Film Private Data
    scale: Float,
    max_sample_luminance: Float,
}

impl Film {
    pub fn new(resolution: Point2i,
               crop_window: Bounds2f,
               filt: BoxFilter,
               diagonal: Float,
               filename: String,
               scale: Float,
               max_sample_luminance: Float)
               -> Self {
        let cropped_pixel_bounds: Bounds2i = Bounds2i {
            p_min: Point2i {
                x: (resolution.x as Float * crop_window.p_min.x).ceil() as i32,
                y: (resolution.y as Float * crop_window.p_min.y).ceil() as i32,
            },
            p_max: Point2i {
                x: (resolution.x as Float * crop_window.p_max.x).ceil() as i32,
                y: (resolution.y as Float * crop_window.p_max.y).ceil() as i32,
            },
        };
        // TODO: allocate film image storage
        // TODO: precompute filter weight table
        Film {
            full_resolution: resolution,
            diagonal: diagonal * 0.001,
            filter: filt,
            filename: filename,
            cropped_pixel_bounds: cropped_pixel_bounds,
            scale: scale,
            max_sample_luminance: max_sample_luminance,
        }
    }
    pub fn get_sample_bounds(&self) -> Bounds2i {
        let f: Point2f = pnt2_floor(Point2f {
            x: self.cropped_pixel_bounds.p_min.x as Float,
            y: self.cropped_pixel_bounds.p_min.y as Float,
        } + Vector2f { x: 0.5, y: 0.5 } - self.filter.radius);
        let c: Point2f = pnt2_ceil(Point2f {
            x: self.cropped_pixel_bounds.p_max.x as Float,
            y: self.cropped_pixel_bounds.p_max.y as Float,
        } + Vector2f { x: 0.5, y: 0.5 } + self.filter.radius);
        let float_bounds: Bounds2f = Bounds2f::new(f, c);
        Bounds2i {
            p_min: Point2i {
                x: float_bounds.p_min.x as i32,
                y: float_bounds.p_min.y as i32,
            },
            p_max: Point2i {
                x: float_bounds.p_max.x as i32,
                y: float_bounds.p_max.y as i32,
            },
        }
    }
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
    dx_camera: Vector3f,
    dy_camera: Vector3f,
    a: Float,
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
               -> Self {
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
        // see perspective.cpp
        // compute differential changes in origin for perspective camera rays
        let dx_camera: Vector3f = raster_to_camera.transform_point(Point3f {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        }) -
                                  raster_to_camera.transform_point(Point3f {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        });
        let dy_camera: Vector3f = raster_to_camera.transform_point(Point3f {
            x: 0.0,
            y: 1.0,
            z: 0.0,
        }) -
                                  raster_to_camera.transform_point(Point3f {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        });
        // compute image plane bounds at $z=1$ for _PerspectiveCamera_
        let res: Point2i = film.full_resolution;
        let mut p_min: Point3f = raster_to_camera.transform_point(Point3f {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        });
        // Point3f p_max = RasterToCamera(Point3f(res.x, res.y, 0));
        let mut p_max: Point3f = raster_to_camera.transform_point(Point3f {
            x: res.x as Float,
            y: res.y as Float,
            z: 0.0,
        });
        p_min /= p_min.z;
        p_max /= p_max.z;
        let a: Float = ((p_max.x - p_min.x) * (p_max.y - p_min.y)).abs();

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
            dx_camera: dx_camera,
            dy_camera: dy_camera,
            a: a,
        }
    }
}

// see light.h

#[repr(i32)]
pub enum LightFlags {
    DeltaPosition = 1,
    DeltaDirection = 2,
    Area = 4,
    Infinite = 8
}

// see distant.h

#[derive(Debug,Copy,Clone)]
pub struct DistantLight {
    // private data (see distant.h)
    l: Spectrum,
    w_light: Vector3f,
    world_center: Point3f,
    world_radius: Float,
    // inherited from class Light (see light.h)
    flags: i32,
    n_samples: i32, // const?
    // TODO: const MediumInterface mediumInterface;
    // TODO: const Transform LightToWorld, WorldToLight;
}

impl DistantLight {
    pub fn new(light_to_world: &Transform, l: &Spectrum, w_light: &Vector3f) -> Self {
        DistantLight {
            l: *l,
            w_light: vec3_normalize(light_to_world.transform_vector(*w_light)),
            world_center: Point3f::default(),
            world_radius: 0.0,
            flags: LightFlags::DeltaDirection as i32,
            n_samples: 1_i32,
        }
    }
    /// Some of the **DistanceLight** methods need to know the bounds
    /// of the scene. Because lights are created before the scene
    /// geometry, these bounds aren't available when the
    /// **DistanceLight** constructor runs. Therefore,
    /// **DistanceLight** implements the optional *preprocess()*
    /// method to get the bound. This method is called at the end of
    /// the **Scene** constructor.
    pub fn preprocess(&mut self, scene: &Scene) {
        Bounds3f::bounding_sphere(scene.world_bound(),
                                  &mut self.world_center,
                                  &mut self.world_radius);
    }
}

// see directlighting.h

#[derive(Debug,Clone,PartialEq)]
pub enum LightStrategy {
    UniformSampleAll,
    UniformSampleOne,
}

#[derive(Debug,Clone)]
pub struct DirectLightingIntegrator {
    // inherited from SamplerIntegrator (see integrator.h)
    camera: PerspectiveCamera, // std::shared_ptr<const Camera> camera;
    sampler: ZeroTwoSequenceSampler, // std::shared_ptr<Sampler> sampler;
    pixel_bounds: Bounds2i,
    // see directlighting.h
    // TODO: const LightStrategy strategy;
    strategy: LightStrategy,
    max_depth: i64,
    n_light_samples: Vec<i32>,
}

impl Default for DirectLightingIntegrator {
    fn default() -> DirectLightingIntegrator {
        DirectLightingIntegrator {
            camera: PerspectiveCamera::default(),
            pixel_bounds: Bounds2i {
                p_min: Point2i { x: 0, y: 0 },
                p_max: Point2i { x: 1280, y: 720 },
            },
            sampler: ZeroTwoSequenceSampler::default(),
            strategy: LightStrategy::UniformSampleAll,
            max_depth: 5,
            n_light_samples: Vec::new(),
        }
    }
}

impl DirectLightingIntegrator {
    pub fn new(strategy: LightStrategy,
               max_depth: i64,
               camera: PerspectiveCamera,
               sampler: ZeroTwoSequenceSampler,
               pixel_bounds: Bounds2i)
               -> Self {
        DirectLightingIntegrator {
            camera: camera,
            sampler: sampler,
            pixel_bounds: pixel_bounds,
            strategy: strategy,
            max_depth: max_depth,
            n_light_samples: Vec::new(),
        }
    }
    pub fn preprocess(&mut self, scene: &Scene) { // , sampler: &mut ZeroTwoSequenceSampler
        if self.strategy == LightStrategy::UniformSampleAll {
            // compute number of samples to use for each light
            for li in 0..scene.lights.len() {
                let ref light = scene.lights[li];
                self.n_light_samples.push(self.sampler.round_count(light.n_samples));
            }
            // request samples for sampling all lights
            for i in 0..self.max_depth {
                for j in 0..scene.lights.len() {
                    self.sampler.request_2d_array(self.n_light_samples[j]);
                    self.sampler.request_2d_array(self.n_light_samples[j]);
                }
            }
        }
    }
    pub fn render(&mut self, scene: &Scene) {
        // SamplerIntegrator::Render (integrator.cpp)
        self.preprocess(scene); // , &mut self.sampler
        let sample_bounds: Bounds2i = self.camera.film.get_sample_bounds();
        println!("sample_bounds = {:?}", sample_bounds);
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
    textures: Vec<ParamSetItem<String>>, // TODO: static std::map<std::string, Spectrum> cachedSpectra;
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
