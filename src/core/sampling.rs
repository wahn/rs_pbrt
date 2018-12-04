//! Draw random samples from a chosen probability distribution.

// std
use std::f32::consts::PI;
use std::sync::Arc;
// pbrt
use core::geometry::{Point2f, Vector2f, Vector3f};
use core::pbrt::clamp_t;
use core::pbrt::Float;
use core::pbrt::{INV_2_PI, INV_4_PI, INV_PI, PI_OVER_2, PI_OVER_4};
use core::rng::Rng;

// see sampling.h

#[derive(Debug, Default, Clone)]
pub struct Distribution1D {
    pub func: Vec<Float>,
    pub cdf: Vec<Float>,
    pub func_int: Float,
}

impl Distribution1D {
    pub fn new(f: Vec<Float>) -> Self {
        let n: usize = f.len();
        // compute integral of step function at $x_i$
        let mut cdf: Vec<Float> = Vec::new();
        cdf.push(0.0 as Float);
        for i in 1..(n + 1) {
            let previous: Float = cdf[i - 1];
            cdf.push(previous + f[i - 1] / n as Float);
        }
        // transform step function integral into CDF
        let func_int: Float = cdf[n];
        if func_int == 0.0 as Float {
            for i in 1..(n + 1) {
                cdf[i] = i as Float / n as Float;
            }
        } else {
            for i in 1..(n + 1) {
                cdf[i] /= func_int;
            }
        }
        Distribution1D {
            func: f,
            cdf: cdf,
            func_int: func_int,
        }
    }
    pub fn count(&self) -> usize {
        self.func.len()
    }
    pub fn sample_continuous(
        &self,
        u: Float,
        pdf: Option<&mut Float>,
        off: Option<&mut usize>,
    ) -> Float {
        // find surrounding CDF segments and _offset_
        // int offset = find_interval((int)cdf.size(),
        //                           [&](int index) { return cdf[index] <= u; });

        // see pbrt.h (int FindInterval(int size, const Predicate &pred) {...})
        let mut first: usize = 0;
        let mut len: usize = self.cdf.len();
        while len > 0 as usize {
            let half: usize = len >> 1;
            let middle: usize = first + half;
            // bisect range based on value of _pred_ at _middle_
            if self.cdf[middle] <= u {
                first = middle + 1;
                len -= half + 1;
            } else {
                len = half;
            }
        }
        let offset: usize = clamp_t(
            first as isize - 1_isize,
            0 as isize,
            self.cdf.len() as isize - 2_isize,
        ) as usize;
        if let Some(off_ref) = off {
            *off_ref = offset;
        }
        // compute offset along CDF segment
        let mut du: Float = u - self.cdf[offset];
        if (self.cdf[offset + 1] - self.cdf[offset]) > 0.0 as Float {
            assert!(self.cdf[offset + 1] > self.cdf[offset]);
            du /= self.cdf[offset + 1] - self.cdf[offset];
        }
        assert!(!du.is_nan());
        // compute PDF for sampled offset
        if pdf.is_some() {
            if self.func_int > 0.0 as Float {
                *pdf.unwrap() = self.func[offset] / self.func_int;
            } else {
                *pdf.unwrap() = 0.0;
            }
        }
        // return $x\in{}[0,1)$ corresponding to sample
        (offset as Float + du) / self.count() as Float
    }
    pub fn sample_discrete(
        &self,
        u: Float,
        pdf: Option<&mut Float>, /* TODO: Float *uRemapped = nullptr */
    ) -> usize {
        // find surrounding CDF segments and _offset_
        // let offset: usize = find_interval(cdf.size(),
        //                           [&](int index) { return cdf[index] <= u; });

        // see pbrt.h (int FindInterval(int size, const Predicate &pred) {...})
        let mut first: usize = 0;
        let mut len: usize = self.cdf.len();
        while len > 0 as usize {
            let half: usize = len >> 1;
            let middle: usize = first + half;
            // bisect range based on value of _pred_ at _middle_
            if self.cdf[middle] <= u {
                first = middle + 1;
                len -= half + 1;
            } else {
                len = half;
            }
        }
        let offset: usize = clamp_t(
            first as isize - 1_isize,
            0 as isize,
            self.cdf.len() as isize - 2_isize,
        ) as usize;
        if pdf.is_some() {
            if self.func_int > 0.0 as Float {
                *pdf.unwrap() = self.func[offset] / (self.func_int * self.func.len() as Float);
            } else {
                *pdf.unwrap() = 0.0;
            }
        }
        // TODO: if (uRemapped)
        //     *uRemapped = (u - cdf[offset]) / (cdf[offset + 1] - cdf[offset]);
        // if (uRemapped) CHECK(*uRemapped >= 0.f && *uRemapped <= 1.f);
        offset
    }
    pub fn discrete_pdf(&self, index: usize) -> Float {
        assert!(index < self.func.len());
        self.func[index] / (self.func_int * self.func.len() as Float)
    }
}

#[derive(Debug, Default, Clone)]
pub struct Distribution2D {
    pub p_conditional_v: Vec<Arc<Distribution1D>>,
    pub p_marginal: Arc<Distribution1D>,
}

impl Distribution2D {
    pub fn new(func: Vec<Float>, nu: i32, nv: i32) -> Self {
        let mut p_conditional_v: Vec<Arc<Distribution1D>> = Vec::new();
        for v in 0..nv {
            // compute conditional sampling distribution for $\tilde{v}$
            let f: Vec<Float> = func[(v * nu) as usize..((v + 1) * nu) as usize].to_vec();
            p_conditional_v.push(Arc::new(Distribution1D::new(f)));
        }
        // compute marginal sampling distribution $p[\tilde{v}]$
        let mut marginal_func: Vec<Float> = Vec::with_capacity(nv as usize);
        for v in 0..nv {
            marginal_func.push(p_conditional_v[v as usize].func_int);
        }
        let p_marginal: Arc<Distribution1D> = Arc::new(Distribution1D::new(marginal_func));
        Distribution2D {
            p_conditional_v: p_conditional_v,
            p_marginal: p_marginal,
        }
    }
    pub fn sample_continuous(&self, u: &Point2f, pdf: &mut Float) -> Point2f {
        let mut pdfs: [Float; 2] = [0.0 as Float; 2];
        let mut v: usize = 0_usize;
        let d1: Float = self
            .p_marginal
            .sample_continuous(u[1], Some(&mut (pdfs[1])), Some(&mut v));
        let d0: Float = self.p_conditional_v[v].sample_continuous(u[0], Some(&mut (pdfs[0])), None);
        *pdf = pdfs[0] * pdfs[1];
        Point2f { x: d0, y: d1 }
    }
    pub fn pdf(&self, p: &Point2f) -> Float {
        let iu: usize = clamp_t(
            (p[0] * self.p_conditional_v[0].count() as Float) as usize,
            0_usize,
            self.p_conditional_v[0].count() - 1_usize,
        );
        let iv: usize = clamp_t(
            (p[1] * self.p_marginal.count() as Float) as usize,
            0_usize,
            self.p_marginal.count() - 1_usize,
        );
        self.p_conditional_v[iv].func[iu] / self.p_marginal.func_int
    }
}

/// Randomly permute an array of *count* sample values, each of which
/// has *n_dimensions* dimensions.
pub fn shuffle<T>(samp: &mut [T], count: i32, n_dimensions: i32, rng: &mut Rng) {
    for i in 0..count {
        let other: i32 = i + rng.uniform_uint32_bounded((count - i) as u32) as i32;
        for j in 0..n_dimensions {
            samp.swap(
                (n_dimensions * i + j) as usize,
                (n_dimensions * other + j) as usize,
            );
        }
    }
}

/// Cosine-weighted hemisphere sampling using Malley's method.
pub fn cosine_sample_hemisphere(u: &Point2f) -> Vector3f {
    let d: Point2f = concentric_sample_disk(u);
    let z: Float = (0.0 as Float)
        .max(1.0 as Float - d.x * d.x - d.y * d.y)
        .sqrt();
    Vector3f {
        x: d.x,
        y: d.y,
        z: z,
    }
}

/// Returns a weight of cos_theta / PI.
pub fn cosine_hemisphere_pdf(cos_theta: Float) -> Float {
    cos_theta * INV_PI
}

/// Reducing the variance according to Veach's heuristic.
pub fn power_heuristic(nf: u8, f_pdf: Float, ng: u8, g_pdf: Float) -> Float {
    let f: Float = nf as Float * f_pdf;
    let g: Float = ng as Float * g_pdf;
    (f * f) / (f * f + g * g)
}

// see sampling.cpp

/// Uniformly sample rays in a hemisphere. Choose a direction.
pub fn uniform_sample_hemisphere(u: &Point2f) -> Vector3f {
    let z: Float = u[0_u8];
    let r: Float = (0.0 as Float).max(1.0 as Float - z * z).sqrt();
    let phi: Float = 2.0 as Float * PI * u[1_u8];
    Vector3f {
        x: r * phi.cos(),
        y: r * phi.sin(),
        z: z,
    }
}

/// Uniformly sample rays in a hemisphere. Probability density
/// function (PDF).
pub fn uniform_hemisphere_pdf() -> Float {
    INV_2_PI
}

/// Uniformly sample rays in a full sphere. Choose a direction.
pub fn uniform_sample_sphere(u: &Point2f) -> Vector3f {
    let z: Float = 1.0 as Float - 2.0 as Float * u[0];
    let r: Float = (0.0 as Float).max(1.0 as Float - z * z).sqrt();
    let phi: Float = 2.0 as Float * PI * u[1];
    Vector3f {
        x: r * phi.cos(),
        y: r * phi.sin(),
        z: z,
    }
}

/// Probability density function (PDF) of a sphere.
pub fn uniform_sphere_pdf() -> Float {
    INV_4_PI
}

/// Uniformly distribute samples over a unit disk.
pub fn concentric_sample_disk(u: &Point2f) -> Point2f {
    // map uniform random numbers to $[-1,1]^2$
    let u_offset: Point2f = *u * 2.0 as Float - Vector2f { x: 1.0, y: 1.0 };
    // handle degeneracy at the origin
    if u_offset.x == 0.0 as Float && u_offset.y == 0.0 as Float {
        return Point2f::default();
    }
    // apply concentric mapping to point
    let theta: Float;
    let r: Float;
    if u_offset.x.abs() > u_offset.y.abs() {
        r = u_offset.x;
        theta = PI_OVER_4 * (u_offset.y / u_offset.x);
    } else {
        r = u_offset.y;
        theta = PI_OVER_2 - PI_OVER_4 * (u_offset.x / u_offset.y);
    }
    Point2f {
        x: theta.cos(),
        y: theta.sin(),
    } * r
}

/// Uniformly sample rays in a cone of directions. Probability density
/// function (PDF).
pub fn uniform_cone_pdf(cos_theta_max: Float) -> Float {
    1.0 as Float / (2.0 as Float * PI * (1.0 as Float - cos_theta_max))
}

/// Samples in a cone of directions about the (0, 0, 1) axis.
pub fn uniform_sample_cone(u: &Point2f, cos_theta_max: Float) -> Vector3f {
    let cos_theta: Float = (1.0 as Float - u[0]) + u[0] * cos_theta_max;
    let sin_theta: Float = (1.0 as Float - cos_theta * cos_theta).sqrt();
    let phi: Float = u[1] * 2.0 as Float * PI;
    Vector3f {
        x: phi.cos() * sin_theta,
        y: phi.sin() * sin_theta,
        z: cos_theta,
    }
}

/// Uniformly distributing samples over isosceles right triangles
/// actually works for any triangle.
pub fn uniform_sample_triangle(u: &Point2f) -> Point2f {
    let su0: Float = u[0].sqrt();
    Point2f {
        x: 1.0 as Float - su0,
        y: u[1] * su0,
    }
}
