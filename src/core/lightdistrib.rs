// std
use std;
use std::sync::{Arc, RwLock};
use std::sync::atomic::{AtomicU64, Ordering};
// pbrt
use core::geometry::{Bounds3f, Normal3f, Point2f, Point3f, Point3i, Vector3f};
use core::integrator::compute_light_power_distribution;
use core::interaction::InteractionCommon;
use core::light::VisibilityTester;
use core::lowdiscrepancy::radical_inverse;
use core::pbrt::{Float, Spectrum};
use core::pbrt::clamp_t;
use core::scene::Scene;
use core::sampling::Distribution1D;

// see lightdistrib.h

/// LightDistribution defines a general interface for classes that
/// provide probability distributions for sampling light sources at a
/// given point in space.
pub trait LightDistribution {
    /// Given a point |p| in space, this method returns a (hopefully
    /// effective) sampling distribution for light sources at that
    /// point.
    fn lookup(&self, p: &Point3f) -> Arc<Distribution1D>;
}

#[derive(Debug, Default)]
struct HashEntry {
    packed_pos: AtomicU64,
    distribution: RwLock<Option<Arc<Distribution1D>>>,
}

pub struct UniformLightDistribution {
    pub distrib: Arc<Distribution1D>,
}

/// The simplest possible implementation of LightDistribution: this
/// returns a uniform distribution over all light sources, ignoring
/// the provided point. This approach works well for very simple
/// scenes, but is quite ineffective for scenes with more than a
/// handful of light sources. (This was the sampling method originally
/// used for the PathIntegrator and the VolPathIntegrator in the
/// printed book, though without the UniformLightDistribution class.)
impl UniformLightDistribution {
    pub fn new(scene: &Scene) -> Self {
        let prob: Vec<Float> = vec![1.0 as Float; scene.lights.len()];
        UniformLightDistribution {
            distrib: Arc::new(Distribution1D::new(prob)),
        }
    }
}

impl LightDistribution for UniformLightDistribution {
    fn lookup(&self, _p: &Point3f) -> Arc<Distribution1D> {
        self.distrib.clone()
    }
}

/// PowerLightDistribution returns a distribution with sampling
/// probability proportional to the total emitted power for each
/// light. (It also ignores the provided point |p|.)  This approach
/// works well for scenes where there the most powerful lights are
/// also the most important contributors to lighting in the scene, but
/// doesn't do well if there are many lights and if different lights
/// are relatively important in some areas of the scene and
/// unimportant in others. (This was the default sampling method used
/// for the BDPT integrator and MLT integrator in the printed book,
/// though also without the PowerLightDistribution class.)
pub struct PowerLightDistribution {
    pub distrib: Option<Arc<Distribution1D>>,
}

impl PowerLightDistribution {
    pub fn new(scene: &Scene) -> Self {
        PowerLightDistribution {
            distrib: compute_light_power_distribution(scene),
        }
    }
}

impl LightDistribution for PowerLightDistribution {
    fn lookup(&self, _p: &Point3f) -> Arc<Distribution1D> {
        if let Some(ref distrib) = self.distrib {
            distrib.clone()
        } else {
            // WARNING: this should only happen if scene.lights.is_empty()
            let prob: Vec<Float> = Vec::new();
            Arc::new(Distribution1D::new(prob))
        }
    }
}

/// A spatially-varying light distribution that adjusts the
/// probability of sampling a light source based on an estimate of its
/// contribution to a region of space.  A fixed voxel grid is imposed
/// over the scene bounds and a sampling distribution is computed as
/// needed for each voxel.
pub struct SpatialLightDistribution {
    pub scene: Scene,
    pub n_voxels: [i32; 3],
    hash_table: Arc<Vec<HashEntry>>,
    pub hash_table_size: usize,
}

impl SpatialLightDistribution {
    pub fn new(scene: &Scene, max_voxels: u32) -> Self {
        // compute the number of voxels so that the widest scene
        // bounding box dimension has maxVoxels voxels and the other
        // dimensions have a number of voxels so that voxels are
        // roughly cube shaped.
        let b: Bounds3f = scene.world_bound();
        let diag: Vector3f = b.diagonal();
        let bmax: Float = diag[b.maximum_extent()];
        let mut n_voxels: [i32; 3] = [0_i32; 3];
        for i in 0..3 {
            n_voxels[i] = std::cmp::max(
                1 as i32,
                (diag[i as u8] / bmax * max_voxels as Float).round() as i32,
            );
            // in the Lookup() method, we require that 20 or fewer
            // bits be sufficient to represent each coordinate
            // value. It's fairly hard to imagine that this would ever
            // be a problem.
            assert!(n_voxels[i] < (1 << 20));
        }
        let hash_table_size: usize = (4 as i32 * n_voxels[0] * n_voxels[1] * n_voxels[2]) as usize;
        let mut hash_table: Vec<HashEntry> = Vec::new();
        // let null: *mut Distribution1D = std::ptr::null_mut();
        for _i in 0..hash_table_size {
            let hash_entry: HashEntry = HashEntry {
                packed_pos: AtomicU64::new(INVALID_PACKED_POS),
                distribution: RwLock::new(None),
            };
            hash_table.push(hash_entry);
        }
        println!(
            "SpatialLightDistribution: scene bounds {:?}, voxel res ({:?}, {:?}, {:?})",
            b, n_voxels[0], n_voxels[1], n_voxels[2]
        );
        SpatialLightDistribution {
            scene: scene.clone(),
            n_voxels: n_voxels,
            hash_table: Arc::new(hash_table),
            hash_table_size: hash_table_size,
        }
    }
    /// Compute the sampling distribution for the voxel with integer
    /// coordiantes given by "pi".
    pub fn compute_distribution(&self, pi: &Point3i) -> Distribution1D {
        // Compute the world-space bounding box of the voxel
        // corresponding to |pi|.
        let p0: Point3f = Point3f {
            x: pi[0] as Float / self.n_voxels[0] as Float,
            y: pi[1] as Float / self.n_voxels[1] as Float,
            z: pi[2] as Float / self.n_voxels[2] as Float,
        };
        let p1: Point3f = Point3f {
            x: (pi[0] + 1) as Float / self.n_voxels[0] as Float,
            y: (pi[1] + 1) as Float / self.n_voxels[1] as Float,
            z: (pi[2] + 1) as Float / self.n_voxels[2] as Float,
        };
        let voxel_bounds: Bounds3f = Bounds3f {
            p_min: self.scene.world_bound().lerp(&p0),
            p_max: self.scene.world_bound().lerp(&p1),
        };
        // Compute the sampling distribution. Sample a number of
        // points inside voxelBounds using a 3D Halton sequence; at
        // each one, sample each light source and compute a weight
        // based on Li/pdf for the light's sample (ignoring visibility
        // between the point in the voxel and the point on the light
        // source) as an approximation to how much the light is likely
        // to contribute to illumination in the voxel.
        let n_samples: usize = 128;
        let mut light_contrib: Vec<Float> = vec![0.0 as Float; self.scene.lights.len()];
        for i in 0..n_samples {
            let po: Point3f = voxel_bounds.lerp(&Point3f {
                x: radical_inverse(0, i as u64),
                y: radical_inverse(1, i as u64),
                z: radical_inverse(2, i as u64),
            });
            let time: Float = 0.0;
            let intr: InteractionCommon = InteractionCommon {
                p: po,
                time: time,
                p_error: Vector3f::default(),
                wo: Vector3f {
                    x: 1.0,
                    y: 0.0,
                    z: 0.0,
                },
                n: Normal3f::default(),
            };
            // Use the next two Halton dimensions to sample a point on the
            // light source.
            let u: Point2f = Point2f {
                x: radical_inverse(3, i as u64),
                y: radical_inverse(4, i as u64),
            };
            for j in 0..self.scene.lights.len() {
                let mut pdf: Float = 0.0 as Float;
                let mut wi: Vector3f = Vector3f::default();
                let mut vis: VisibilityTester = VisibilityTester::default();
                let li: Spectrum =
                    self.scene.lights[j].sample_li(&intr, &u, &mut wi, &mut pdf, &mut vis);
                if pdf > 0.0 as Float {
                    // TODO: look at tracing shadow rays / computing
                    // beam transmittance. Probably shouldn't give
                    // those full weight but instead e.g. have an
                    // occluded shadow ray scale down the contribution
                    // by 10 or something.
                    light_contrib[j] += li.y() / pdf;
                }
            }
        }
        // We don't want to leave any lights with a zero probability;
        // it's possible that a light contributes to points in the
        // voxel even though we didn't find such a point when sampling
        // above. Therefore, compute a minimum (small) weight and
        // ensure that all lights are given at least the corresponding
        // probability.
        let sum_contrib: Float = light_contrib.iter().sum();
        let avg_contrib: Float = sum_contrib / (n_samples * light_contrib.len()) as Float;
        let min_contrib: Float;
        if avg_contrib > 0.0 as Float {
            min_contrib = 0.001 * avg_contrib;
        } else {
            min_contrib = 1.0 as Float;
        }
        for i in 0..light_contrib.len() {
            // println!("Voxel pi = {:?}, light {:?} contrib = {:?}",
            //          pi, i, light_contrib[i]);
            light_contrib[i] = light_contrib[i].max(min_contrib);
        }
        // println!("Initialized light distribution in voxel pi= {:?}, avg_contrib = {:?}",
        //          pi, avg_contrib);
        // Compute a sampling distribution from the accumulated contributions.
        Distribution1D::new(light_contrib)
    }
}

impl LightDistribution for SpatialLightDistribution {
    fn lookup(&self, p: &Point3f) -> Arc<Distribution1D> {
        // TODO: ProfilePhase _(Prof::LightDistribLookup);
        // TODO: ++nLookups;

        // first, compute integer voxel coordinates for the given
        // point |p| with respect to the overall voxel grid.
        let offset: Vector3f = self.scene.world_bound().offset(&p); // offset in [0,1].
        let mut pi: Point3i = Point3i::default();
        for i in 0..3 {
            // the clamp should almost never be necessary, but is
            // there to be robust to computed intersection points
            // being slightly outside the scene bounds due to
            // floating-point roundoff error.
            pi[i] = clamp_t(
                (offset[i] * self.n_voxels[i as usize] as Float) as i32,
                0_i32,
                self.n_voxels[i as usize] - 1_i32,
            );
        }
        // pack the 3D integer voxel coordinates into a single 64-bit value.
        let packed_pos: u64 = ((pi[0] as u64) << 40) | ((pi[1] as u64) << 20) | (pi[2] as u64);
        assert_ne!(packed_pos, INVALID_PACKED_POS);
        // Compute a hash value from the packed voxel coordinates. We
        // could just take packed_Pos mod the hash table size, but
        // since packed_Pos isn't necessarily well distributed on its
        // own, it's worthwhile to do a little work to make sure that
        // its bits values are individually fairly random. For details
        // of and motivation for the following, see:
        // http://zimbry.blogspot.ch/2011/09/better-bit-mixing-improving-on.html
        let mut hash: u64 = packed_pos;
        hash ^= hash >> 31;
        // hash *= 0x7fb5d329728ea185;
        let (mul, _overflow) = hash.overflowing_mul(0x7fb5d329728ea185);
        hash = mul;
        hash ^= hash >> 27;
        // hash *= 0x81dadef4bc2dd44d;
        let (mul, _overflow) = hash.overflowing_mul(0x81dadef4bc2dd44d);
        hash = mul;
        hash ^= hash >> 33;
        hash %= self.hash_table_size as u64;
        // // hash ^= hash >> 31;
        // let (shr, _overflow) = hash.overflowing_shr(31);
        // hash ^= shr;
        // // hash ^= hash >> 27;
        // let (shr, _overflow) = hash.overflowing_shr(27);
        // hash ^= shr;
        // // hash ^= hash >> 33;
        // let (shr, _overflow) = hash.overflowing_shr(33);
        // hash ^= shr;
        // // hash %= self.hash_table_size as u64;
        // let (rem, _overflow) = hash.overflowing_rem(self.hash_table_size as u64);
        // hash = rem;
        // BELOW: comparison is useless due to type limits
        // assert!(hash >= 0_u64, "hash needs to be greater or equal zero");
        // Now, see if the hash table already has an entry for the
        // voxel. We'll use quadratic probing when the hash table
        // entry is already used for another value; step stores the
        // square root of the probe step.
        let mut step: u64 = 1;
        // TODO: int nProbes = 0;
        loop {
            // TODO: ++nProbes;
            let entry: &HashEntry = &self.hash_table[hash as usize];
            // does the hash table entry at offset |hash| match the current point?
            let entry_packed_pos: u64 = entry.packed_pos.load(Ordering::Acquire);
            if entry_packed_pos == packed_pos {
                // Yes! Most of the time, there should already by a light
                // sampling distribution available.
                let option: &Option<Arc<Distribution1D>> = &*entry.distribution.read().unwrap();
                if let Some(ref dist) = *option {
                    // We have a valid sampling distribution.
                    return Arc::clone(dist);
                }
            } else if entry_packed_pos != INVALID_PACKED_POS {
                // The hash table entry we're checking has already
                // been allocated for another voxel. Advance to the
                // next entry with quadratic probing.
                hash += step * step;
                if hash >= self.hash_table_size as u64 {
                    hash %= self.hash_table_size as u64;
                }
                step += 1_u64;
            } else {
                // We have found an invalid entry. (Though this may
                // have changed since the load into entryPackedPos
                // above.)  Use an atomic compare/exchange to try to
                // claim this entry for the current position.
                let invalid: u64 = INVALID_PACKED_POS;
                let success = match entry.packed_pos.compare_exchange_weak(
                    invalid,
                    packed_pos,
                    Ordering::SeqCst,
                    Ordering::Relaxed,
                ) {
                    Ok(_) => true,
                    Err(_) => false,
                };
                if success {
                    // Success; we've claimed this position for this
                    // voxel's distribution. Now compute the sampling
                    // distribution and add it to the hash table.
                    let dist: Distribution1D = self.compute_distribution(&pi);
                    let arc_dist: Arc<Distribution1D> = Arc::new(dist.clone());
                    let mut distribution = entry.distribution.write().unwrap();
                    *distribution = Some(arc_dist.clone());
                    return arc_dist;
                }
            }
        }
    }
}

// see lightdistrib.cpp

const INVALID_PACKED_POS: u64 = 0xffffffffffffffff;

/// Decides based on the name and the number of scene lights which
/// light distribution to return.
pub fn create_light_sample_distribution(
    name: String,
    scene: &Scene,
) -> Option<Arc<LightDistribution + Send + Sync>> {
    if name == String::from("uniform") || scene.lights.len() == 1 {
        return Some(Arc::new(UniformLightDistribution::new(scene)));
    } else if name == String::from("power") {
        return Some(Arc::new(PowerLightDistribution::new(scene)));
    } else if name == String::from("spatial") {
        return Some(Arc::new(SpatialLightDistribution::new(scene, 64)));
    } else {
        println!(
            "Light sample distribution type \"{:?}\" unknown. Using \"spatial\".",
            name
        );
        return Some(Arc::new(SpatialLightDistribution::new(scene, 64)));
    }
}
