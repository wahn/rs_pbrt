//! All filter implementations are derived from an abstract **Filter**
//! class, which provides the interface for the functions used in
//! filtering.
//!
//! - BoxFilter
//! - GaussianFilter
//! - MitchellFilter
//! - LanczosSincFilter
//! - TriangleFilter
//!
//! ## Box Filter
//!
//! One of the most commonly used filters in graphics is the box
//! filter. The box filter equally weights all samples within a square
//! region of the image. Although computational efficient, it's just
//! about the worst filter possible.
//!
//! ```rust
//! use rs_pbrt::core::pbrt::Float;
//! use rs_pbrt::filters::boxfilter::BoxFilter;
//! use rs_pbrt::core::geometry::Vector2f;
//!
//! fn main() {
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
//!
//! ## Gaussian Filter
//!
//! Unlike the box and triangle filters, the Gaussian filter gives a
//! reasonably good result in practice. The Gaussian filter does tend
//! to cause slight blurring of the final image compared to some of
//! the other filters, but this blurring can actually help mask any
//! remaining aliasing in the image.
//!
//! ```rust
//! use rs_pbrt::core::pbrt::Float;
//! use rs_pbrt::filters::gaussian::GaussianFilter;
//! use rs_pbrt::core::geometry::Vector2f;
//!
//! fn main() {
//!     let xw: Float = 2.0;
//!     let yw: Float = 2.0;
//!     let alpha: Float = 2.0;
//!     let exp_x: Float = (-alpha * xw * xw).exp();
//!     let exp_y: Float = (-alpha * yw * yw).exp();
//!     let gaussian_filter = GaussianFilter {
//!         alpha: alpha,
//!         exp_x: exp_x,
//!         exp_y: exp_y,
//!         radius: Vector2f { x: xw, y: yw },
//!         inv_radius: Vector2f {
//!             x: 1.0 / xw,
//!             y: 1.0 / yw,
//!         },
//!     };
//!
//!     println!("gaussian_filter = {:?}", gaussian_filter);
//! }
//! ```
//!
//! ## MitchellFilter
//!
//! ```rust
//! use rs_pbrt::core::geometry::Vector2f;
//! use rs_pbrt::core::pbrt::Float;
//! use rs_pbrt::filters::mitchell::MitchellNetravali;
//!
//! fn main() {
//!     let xwidth: Float = 2.0;
//!     let ywidth: Float = 2.0;
//!     let b: Float = 1.0 / 3.0;
//!     let c: Float = 1.0 / 3.0;
//!     let mitchell_filter = MitchellNetravali::new(xwidth, ywidth, b, c);
//!
//!     println!("mitchell_filter = {:?}", mitchell_filter);
//! }
//! ```
//!
//! ## LanczosSincFilter
//!
//! ```rust
//! use rs_pbrt::core::geometry::Vector2f;
//! use rs_pbrt::core::pbrt::Float;
//! use rs_pbrt::filters::sinc::LanczosSincFilter;
//!
//! fn main() {
//!     let xw: Float = 4.0;
//!     let yw: Float = 4.0;
//!     let radius: Vector2f = Vector2f { x: xw, y: yw };
//!     let tau: Float = 3.0;
//!     let sinc_filter = LanczosSincFilter::new(radius, tau);
//!
//!     println!("sinc_filter = {:?}", sinc_filter);
//! }
//! ```
//!
//! ## TriangleFilter
//!
//! ```rust
//! use rs_pbrt::core::geometry::Vector2f;
//! use rs_pbrt::core::pbrt::Float;
//! use rs_pbrt::filters::triangle::TriangleFilter;
//!
//! fn main() {
//!     let xw: Float = 2.0;
//!     let yw: Float = 2.0;
//!     let triangle_filter = TriangleFilter {
//!         radius: Vector2f { x: xw, y: yw },
//!         inv_radius: Vector2f {
//!             x: 1.0 / xw,
//!             y: 1.0 / yw,
//!         },
//!     };
//!
//!     println!("triangle_filter = {:?}", triangle_filter);
//! }
//! ```
//!

pub mod boxfilter;
pub mod gaussian;
pub mod mitchell;
pub mod sinc;
pub mod triangle;
