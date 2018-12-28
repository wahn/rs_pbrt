//std
use std;
use std::f32::consts::PI;
// pbrt
use core::interpolation::integrate_catmull_rom;
use core::pbrt::{Float, Spectrum};

pub struct BSSRDFTable {
    pub n_rho_samples: i32,
    pub n_radius_samples: i32,
    pub rho_samples: Vec<Float>,
    pub radius_samples: Vec<Float>,
    pub profile: Vec<Float>,
    pub rho_eff: Vec<Float>,
    pub profile_cdf: Vec<Float>,
}

impl BSSRDFTable {
    pub fn new(n_rho_samples: i32, n_radius_samples: i32) -> Self {
        BSSRDFTable {
            n_rho_samples: n_rho_samples,
            n_radius_samples: n_radius_samples,
            rho_samples: Vec::with_capacity(n_rho_samples as usize),
            radius_samples: Vec::with_capacity(n_radius_samples as usize),
            profile: Vec::with_capacity((n_radius_samples * n_rho_samples) as usize),
            rho_eff: Vec::with_capacity(n_rho_samples as usize),
            profile_cdf: Vec::with_capacity((n_radius_samples * n_rho_samples) as usize),
        }
    }
    pub fn eval_profile(&self, rho_index: i32, radius_index: i32) -> Float {
        self.profile[(rho_index * self.n_radius_samples + radius_index) as usize]
    }
}

pub fn beam_diffusion_ms(sigma_s: Float, sigma_a: Float, g: Float, eta: Float, r: Float) -> Float {
    // WORK
    0.0 as Float
}

pub fn beam_diffusion_ss(sigma_s: Float, sigma_a: Float, g: Float, eta: Float, r: Float) -> Float {
    // WORK
    0.0 as Float
}

pub fn compute_beam_diffusion_bssrdf(g: Float, eta: Float, t: &mut BSSRDFTable) {
    // choose radius values of the diffusion profile discretization
    t.radius_samples[0] = 0.0 as Float;
    t.radius_samples[1] = 2.5e-3 as Float;
    for i in 2..t.n_radius_samples as usize {
        t.radius_samples[i] = t.radius_samples[i - 1] * 1.2 as Float;
    }
    // choose albedo values of the diffusion profile discretization
    for i in 0..t.n_rho_samples as usize {
        t.rho_samples[i] = (1.0 as Float
            - (-8.0 as Float * i as Float / (t.n_rho_samples as Float - 1.0 as Float)).exp())
            / (1.0 as Float - (-8.0 as Float).exp());
    }
    // ParallelFor([&](int i) {
    for i in 0..t.n_rho_samples as usize {
        // compute the diffusion profile for the _i_th albedo sample

        // compute scattering profile for chosen albedo $\rho$
        for j in 0..t.n_radius_samples as usize {
            //         Float rho = t.rho_samples[i], r = t.radius_samples[j];
            let rho: Float = t.rho_samples[i];
            let r: Float = t.radius_samples[j];
            t.profile[i * t.n_radius_samples as usize + j] = 2.0 as Float
                * PI
                * r
                * (beam_diffusion_ss(rho, 1.0 as Float - rho, g, eta, r)
                    + beam_diffusion_ms(rho, 1.0 as Float - rho, g, eta, r));
        }
        // compute effective albedo $\rho_{\roman{eff}}$ and CDF for
        // importance sampling
        t.rho_eff[i] = integrate_catmull_rom(
            t.n_radius_samples,
            &t.radius_samples,
            &t.profile,
            i * t.n_radius_samples as usize,
            &mut t.profile_cdf,
            i * t.n_radius_samples as usize,
        );
    }
    // }, t.n_rho_samples);
}
