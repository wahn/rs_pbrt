//! To limit the potential number of texels that need to be accessed,
//! filtering methods use an image pyramid of increasingly lower
//! resolution prefiltered versions of the original image to
//! accelerate their operation.

// std
use std;
use std::ops::{Add, AddAssign, Div, Mul};
// others
use num;
// pbrt
use crate::core::geometry::{Point2f, Point2i, Vector2f};
use crate::core::memory::BlockedArray;
use crate::core::pbrt::{clamp_t, is_power_of_2, lerp, mod_t, round_up_pow2_32};
use crate::core::pbrt::{Float, Spectrum};
use crate::core::texture::lanczos;

// see mipmap.h

const WEIGHT_LUT_SIZE: usize = 128;

#[derive(Debug, Clone)]
pub enum ImageWrap {
    Repeat,
    Black,
    Clamp,
}

#[derive(Debug, Default, Copy, Clone)]
pub struct ResampleWeight {
    pub first_texel: i32,
    pub weight: [Float; 4],
}

pub struct MipMap<T> {
    // MIPMap Private Data
    pub do_trilinear: bool,
    pub max_anisotropy: Float,
    pub wrap_mode: ImageWrap,
    pub resolution: Point2i,
    pub pyramid: Vec<BlockedArray<T>>,
    // TODO: static Float weightLut[WeightLUTSize];
    pub weight_lut: [Float; WEIGHT_LUT_SIZE],
}

impl<T> MipMap<T>
where
    T: std::default::Default
        + num::Zero
        + std::clone::Clone
        + Add<T, Output = T>
        + AddAssign
        + Clampable
        + Copy
        + Div<Float, Output = T>
        + Mul<T, Output = T>
        + Mul<Float, Output = T>,
{
    pub fn new(
        res: Point2i,
        img: &[T],
        do_trilinear: bool,
        max_anisotropy: Float,
        wrap_mode: ImageWrap,
    ) -> Self {
        let mut resolution = res;
        let mut resampled_image: Vec<T> = Vec::new();
        if !is_power_of_2(resolution.x) || !is_power_of_2(resolution.y) {
            // resample image to power-of-two resolution
            let res_pow_2: Point2i = Point2i {
                x: round_up_pow2_32(resolution.x),
                y: round_up_pow2_32(resolution.y),
            };
            // println!(
            //     "Resampling MIPMap from {:?} to {:?}. Ratio= {:?}",
            //     resolution,
            //     res_pow_2,
            //     (res_pow_2.x * res_pow_2.y) as Float / (resolution.x * resolution.y) as Float
            // );
            // resample image in $s$ direction
            let s_weights: Vec<ResampleWeight> =
                MipMap::<T>::resample_weights(resolution.x, res_pow_2.x);
            // TODO: resampled_image.reset(new T[resPow2[0] * resPow2[1]]);
            resampled_image = vec![T::default(); (res_pow_2.x * res_pow_2.y) as usize];
            // apply _s_weights_ to zoom in $s$ direction
            // TODO: ParallelFor([&](int t) {
            for t in 0..resolution.y {
                // chunk size 16
                for s in 0..res_pow_2.x {
                    // compute texel $(s,t)$ in $s$-zoomed image
                    resampled_image[(t * res_pow_2.x + s) as usize] = T::default();
                    for j in 0..4 {
                        let mut orig_s: i32 = s_weights[s as usize].first_texel + j as i32;
                        orig_s = match wrap_mode {
                            ImageWrap::Repeat => mod_t(orig_s, resolution.x),
                            ImageWrap::Clamp => clamp_t(orig_s, 0_i32, resolution.x - 1_i32),
                            _ => orig_s,
                        };
                        if orig_s >= 0_i32 && orig_s < resolution.x {
                            resampled_image[(t * res_pow_2.x + s) as usize] += img
                                [(t * resolution.x + orig_s) as usize]
                                * s_weights[s as usize].weight[j];
                        }
                    }
                }
            }
            // TODO: }, resolution[1], 16);
            // resample image in $t$ direction
            let t_weights: Vec<ResampleWeight> =
                MipMap::<T>::resample_weights(resolution.y, res_pow_2.y);
            // std::vector<T *> resample_bufs;
            // int nThreads = MaxThreadIndex();
            // for (int i = 0; i < nThreads; ++i)
            //     resample_bufs.push_back(new T[resPow2[1]]);
            // let resampled_bufs: Vec<T> = vec![T::default(); res_pow_2.y as usize]; // single-threaded
            let mut work_data: Vec<T> = vec![T::default(); res_pow_2.y as usize]; // single-threaded
                                                                                  // TODO: ParallelFor([&](int s) {
                                                                                  // T *work_data = resample_bufs[ThreadIndex];
            for s in 0..res_pow_2.x {
                // chunk size 32
                for t in 0..res_pow_2.y {
                    work_data[t as usize] = T::default();
                    for j in 0..4 {
                        let mut offset: i32 = t_weights[t as usize].first_texel + j as i32;
                        offset = match wrap_mode {
                            ImageWrap::Repeat => mod_t(offset, resolution.y),
                            ImageWrap::Clamp => clamp_t(offset, 0_i32, resolution.y - 1_i32),
                            _ => offset,
                        };
                        if offset >= 0_i32 && offset < resolution.y {
                            work_data[t as usize] += resampled_image
                                [(offset * res_pow_2.x + s) as usize]
                                * t_weights[t as usize].weight[j];
                        }
                    }
                }
                for t in 0..res_pow_2.y {
                    resampled_image[(t * res_pow_2.x + s) as usize] = Clampable::clamp(
                        work_data[t as usize],
                        0.0 as Float,
                        std::f32::INFINITY as Float,
                    );
                }
            }
            // TODO: }, resPow2[0], 32);
            // for (auto ptr : resample_bufs) delete[] ptr;
            resolution = res_pow_2;
        }
        let mut mipmap = MipMap::<T> {
            do_trilinear,
            max_anisotropy,
            wrap_mode,
            resolution,
            pyramid: Vec::new(),
            weight_lut: [0.0 as Float; WEIGHT_LUT_SIZE],
        };
        // initialize levels of MipMap for image
        let n_levels = 1 + (std::cmp::max(resolution.x, resolution.y) as Float).log2() as usize;
        // initialize most detailed level of MipMap
        let img_data: &[T] = if resampled_image.is_empty() {
            img
        } else {
            &resampled_image[..]
        };
        mipmap.pyramid.push(BlockedArray::new_from(
            resolution.x as usize,
            resolution.y as usize,
            img_data,
        ));
        for i in 1..n_levels {
            // initialize $i$th MipMap level from $i-1$st level
            let s_res = std::cmp::max(1, mipmap.pyramid[i - 1].u_size() / 2);
            let t_res = std::cmp::max(1, mipmap.pyramid[i - 1].v_size() / 2);
            let mut ba = BlockedArray::<T>::new(s_res, t_res);
            // filter 4 texels from finer level of pyramid
            for t in 0..t_res {
                for s in 0..s_res {
                    let (si, ti) = (s as isize, t as isize);
                    ba[(s, t)] = (*mipmap.texel(i - 1, 2 * si, 2 * ti)
                        + *mipmap.texel(i - 1, 2 * si + 1, 2 * ti)
                        + *mipmap.texel(i - 1, 2 * si, 2 * ti + 1)
                        + *mipmap.texel(i - 1, 2 * si + 1, 2 * ti + 1))
                        as T
                        * 0.25 as Float;
                }
            }
            mipmap.pyramid.push(ba);
        }
        // initialize EWA filter weights if needed
        if mipmap.weight_lut[0] == 0.0 as Float {
            for i in 0..WEIGHT_LUT_SIZE {
                let alpha: Float = 2.0 as Float;
                let r2: Float = i as Float / (WEIGHT_LUT_SIZE - 1) as Float;
                mipmap.weight_lut[i] = (-alpha * r2).exp() - (-alpha).exp();
            }
        }
        // TODO: mipMapMemory += (4 * resolution[0] * resolution[1] * sizeof(T)) / 3;
        mipmap
    }
    pub fn width(&self) -> i32 {
        self.resolution.x
    }
    pub fn height(&self) -> i32 {
        self.resolution.y
    }
    pub fn levels(&self) -> usize {
        self.pyramid.len()
    }
    pub fn texel(&self, level: usize, s: isize, t: isize) -> &T {
        let l = &self.pyramid[level];
        let (u_size, v_size) = (l.u_size() as isize, l.v_size() as isize);
        let (ss, tt): (usize, usize) = match self.wrap_mode {
            ImageWrap::Repeat => (
                mod_t(s as usize, u_size as usize),
                mod_t(t as usize, v_size as usize),
            ),
            ImageWrap::Clamp => (
                clamp_t(s, 0, u_size - 1) as usize,
                clamp_t(t, 0, v_size - 1) as usize,
            ),
            ImageWrap::Black => {
                // TODO: let black: T = num::Zero::zero();
                if s < 0 || s >= u_size || t < 0 || t >= v_size {
                    // TODO: return &black;
                    (
                        clamp_t(s, 0, u_size - 1) as usize,
                        clamp_t(t, 0, v_size - 1) as usize,
                    ) // TMP
                } else {
                    (s as usize, t as usize)
                }
            }
        };
        &l[(ss, tt)]
    }
    pub fn lookup_pnt_flt(&self, st: Point2f, width: Float) -> T {
        // TODO: ++nTrilerpLookups;
        // TODO: ProfilePhase p(Prof::TexFiltTrilerp);
        // compute MIPMap level for trilinear filtering
        let level: Float = self.levels() as Float - 1.0 as Float + width.max(1e-8 as Float).log2();
        // perform trilinear interpolation at appropriate MIPMap level
        if level < 0.0 as Float {
            self.triangle(0_usize, st)
        } else if level >= self.levels() as Float - 1.0 as Float {
            *self.texel(self.levels() - 1, 0_isize, 0_isize)
        } else {
            let i_level: usize = level.floor() as usize;
            let delta: Float = level - i_level as Float;
            lerp(
                delta,
                self.triangle(i_level, st),
                self.triangle(i_level + 1_usize, st),
            )
        }
    }
    pub fn lookup_pnt_vec_vec(&self, st: Point2f, dst0: &mut Vector2f, dst1: &mut Vector2f) -> T {
        if self.do_trilinear {
            let width: Float = dst0
                .x
                .abs()
                .max(dst0.y.abs())
                .max(dst1.x.abs().max(dst1.y.abs()));
            return self.lookup_pnt_flt(st, width);
        }
        // TODO: ++nEWALookups;
        // TODO: ProfilePhase p(Prof::TexFiltEWA);
        // compute ellipse minor and major axes
        if dst0.length_squared() < dst1.length_squared() {
            // std::swap(dst0, dst1);
            let swap: Vector2f = Vector2f {
                x: dst0.x,
                y: dst0.y,
            };
            // dst0 = dst1
            dst0.x = dst1.x;
            dst0.y = dst1.y;
            // dst1 = dst0
            dst1.x = swap.x;
            dst1.y = swap.y;
        }
        let major_length: Float = dst0.length();
        let mut minor_length: Float = dst1.length();
        // clamp ellipse eccentricity if too large
        if minor_length * self.max_anisotropy < major_length && minor_length > 0.0 as Float {
            let scale: Float = major_length / (minor_length * self.max_anisotropy);
            *dst1 *= scale;
            minor_length *= scale;
        }
        if minor_length == 0.0 as Float {
            return self.triangle(0, st);
        }
        // choose level of detail for EWA lookup and perform EWA filtering
        let lod: Float = (0.0 as Float)
            .max(self.levels() as Float - 1.0 as Float + minor_length.log2() as Float);
        let ilod: usize = lod.floor() as usize;
        let col2: T = self.ewa(ilod + 1, st.clone(), dst0.clone(), dst1.clone());
        let col1: T = self.ewa(ilod, st.clone(), dst0.clone(), dst1.clone());
        let ret: T = lerp(lod - ilod as Float, col1, col2);
        ret
    }
    fn resample_weights(old_res: i32, new_res: i32) -> Vec<ResampleWeight> {
        assert!(new_res >= old_res);
        let mut wt: Vec<ResampleWeight> = Vec::with_capacity(new_res as usize);
        let filterwidth: Float = 2.0 as Float;
        for i in 0..new_res {
            // compute image resampling weights for _i_th texel
            let center: Float = (i as Float + 0.5 as Float) * old_res as Float / new_res as Float;
            let mut rw: ResampleWeight = ResampleWeight::default();
            rw.first_texel = ((center - filterwidth) + 0.5 as Float).floor() as i32;
            for j in 0..4 {
                let pos: Float = rw.first_texel as Float + j as Float + 0.5 as Float;
                rw.weight[j] = lanczos((pos - center) / filterwidth, 2.0 as Float);
            }
            // normalize filter weights for texel resampling
            let inv_sum_wts: Float =
                1.0 as Float / (rw.weight[0] + rw.weight[1] + rw.weight[2] + rw.weight[3]);
            for j in 0..4 {
                rw.weight[j] *= inv_sum_wts;
            }
            wt.push(rw); // add to vector
        }
        wt
    }
    fn triangle(&self, level: usize, st: Point2f) -> T {
        let level: usize = clamp_t(level, 0_usize, self.levels() - 1_usize);
        let s: Float = st.x * self.pyramid[level].u_size() as Float - 0.5;
        let t: Float = st.y * self.pyramid[level].v_size() as Float - 0.5;
        let s0: isize = s.floor() as isize;
        let t0: isize = t.floor() as isize;
        let ds: Float = s - s0 as Float;
        let dt: Float = t - t0 as Float;
        let tmp1: T = *self.texel(level, s0 + 1, t0 + 1) * (ds * dt);
        let tmp2: T = *self.texel(level, s0 + 1, t0) * (ds * (1.0 - dt));
        let tmp3: T = *self.texel(level, s0, t0 + 1) * ((1.0 - ds) * dt);
        let tmp4: T = *self.texel(level, s0, t0) * ((1.0 - ds) * (1.0 - dt));
        tmp4 + tmp3 + tmp2 + tmp1
    }
    fn ewa(&self, level: usize, st: Point2f, dst0: Vector2f, dst1: Vector2f) -> T {
        if level >= self.levels() {
            return *self.texel(self.levels() - 1, 0, 0);
        }
        // convert EWA coordinates to appropriate scale for level
        let mut new_st: Vector2f = Vector2f { x: st.x, y: st.y };
        new_st.x = new_st.x * self.pyramid[level].u_size() as Float - 0.5 as Float;
        new_st.y = new_st.y * self.pyramid[level].v_size() as Float - 0.5 as Float;
        let mut new_dst0: Vector2f = Vector2f {
            x: dst0.x,
            y: dst0.y,
        };
        let mut new_dst1: Vector2f = Vector2f {
            x: dst1.x,
            y: dst1.y,
        };
        new_dst0.x *= self.pyramid[level].u_size() as Float;
        new_dst0.y *= self.pyramid[level].v_size() as Float;
        new_dst1.x *= self.pyramid[level].u_size() as Float;
        new_dst1.y *= self.pyramid[level].v_size() as Float;
        // compute ellipse coefficients to bound EWA filter region
        let mut a: Float = new_dst0.y * new_dst0.y + new_dst1.y * new_dst1.y + 1.0 as Float;
        let mut b: Float = -2.0 as Float * (new_dst0.x * new_dst0.y + new_dst1.x * new_dst1.y);
        let mut c: Float = new_dst0.x * new_dst0.x + new_dst1.x * new_dst1.x + 1.0 as Float;
        let inv_f: Float = 1.0 as Float / (a * c - b * b * 0.25 as Float);
        a *= inv_f;
        b *= inv_f;
        c *= inv_f;
        // compute the ellipse's $(s,t)$ bounding box in texture space
        let det: Float = -b * b + 4.0 as Float * a * c;
        let inv_det: Float = 1.0 as Float / det;
        let u_sqrt: Float = (det * c).sqrt();
        let v_sqrt: Float = (a * det).sqrt();
        let s0: isize = (new_st.x - 2.0 as Float * inv_det * u_sqrt).ceil() as isize;
        let s1: isize = (new_st.x + 2.0 as Float * inv_det * u_sqrt).floor() as isize;
        let t0: isize = (new_st.y - 2.0 as Float * inv_det * v_sqrt).ceil() as isize;
        let t1: isize = (new_st.y + 2.0 as Float * inv_det * v_sqrt).floor() as isize;
        // scan over ellipse bound and compute quadratic equation
        let mut sum: T = T::default();
        let mut sum_wts: Float = 0.0;
        for it in t0..=t1 {
            let tt: Float = it as Float - new_st.y;
            for is in s0..=s1 {
                let ss: Float = is as Float - new_st.x;
                // compute squared radius and filter texel if inside ellipse
                let r2: Float = a * ss * ss + b * ss * tt + c * tt * tt;
                if r2 < 1.0 as Float {
                    let index: usize = std::cmp::min(
                        (r2 * WEIGHT_LUT_SIZE as Float) as usize,
                        WEIGHT_LUT_SIZE - 1,
                    );
                    let weight: Float = self.weight_lut[index];
                    sum += *self.texel(level, is as isize, it as isize) * weight;
                    sum_wts += weight;
                }
            }
        }
        sum / sum_wts
    }
}

pub trait Clampable {
    fn clamp(self, min: Float, max: Float) -> Self;
}

impl Clampable for Float {
    fn clamp(self, min: Float, max: Float) -> Float {
        clamp_t(self, min, max)
    }
}

impl Clampable for Spectrum {
    fn clamp(self, min: Float, max: Float) -> Spectrum {
        Spectrum::rgb(
            clamp_t(self.c[0], min, max),
            clamp_t(self.c[1], min, max),
            clamp_t(self.c[2], min, max),
        )
    }
}
