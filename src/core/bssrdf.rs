//std
use std;
use std::f32::consts::PI;
use std::sync::Arc;
// pbrt
use core::geometry::{Normal3f, Point2f, Point3f, Ray, Vector3f};
use core::interaction::SurfaceInteraction;
use core::interpolation::integrate_catmull_rom;
use core::material::{Material, TransportMode};
use core::medium::phase_hg;
use core::pbrt::INV_4_PI;
use core::pbrt::{Float, Spectrum};
use core::reflection::fr_dielectric;
use core::scene::Scene;

pub trait Bssrdf {
    fn s(&self, pi: &SurfaceInteraction, wi: &Vector3f) -> Spectrum;
    fn sample_s(
        &self,
        scene: &Scene,
        u1: Float,
        u2: &Point2f,
        si: &mut SurfaceInteraction,
        pdf: &mut Float,
    ) -> Spectrum;
}

pub struct TabulatedBssrdf {
    // BSSRDF Protected Data
    // pub po: SurfaceInteraction,
    pub eta: Float,
    // SeparableBSSRDF Private Data
    pub ns: Normal3f,
    pub ss: Vector3f,
    pub ts: Vector3f,
    pub material: Arc<Material>,
    pub mode: TransportMode,
    // TabulatedBSSRDF Private Data
    pub table: Arc<BssrdfTable>,
    pub sigma_t: Spectrum,
    pub rho: Spectrum,
}

impl TabulatedBssrdf {
    pub fn new(
        // TODO: po: SurfaceInteraction,
        material: Arc<Material>,
        mode: TransportMode,
        eta: Float,
        sigma_a: &Spectrum,
        sigma_s: &Spectrum,
        table: Arc<BssrdfTable>,
    ) -> Self {
        let sigma_t: Spectrum = *sigma_a + *sigma_s;
        let mut rho: Spectrum = Spectrum::new(0.0 as Float);
        for c in 0..3 {
            if sigma_t[c] != 0.0 as Float {
                rho.c[c] = sigma_s[c] / sigma_t[c];
            } else {
                rho.c[c] = 0.0 as Float;
            }
        }
        TabulatedBssrdf {
            // TODO: po:
            eta: eta,
            ns: Normal3f::default(), // TODO: po.shading.n
            ss: Vector3f::default(), // TODO: normalize(po.shading.dpdu)
            ts: Vector3f::default(), // TODO: cross(ns, ss)
            material: material,
            mode: mode,
            table: table.clone(),
            sigma_t: sigma_t,
            rho: rho,
        }
    }
}

pub struct BssrdfTable {
    pub n_rho_samples: i32,
    pub n_radius_samples: i32,
    pub rho_samples: Vec<Float>,
    pub radius_samples: Vec<Float>,
    pub profile: Vec<Float>,
    pub rho_eff: Vec<Float>,
    pub profile_cdf: Vec<Float>,
}

impl BssrdfTable {
    pub fn new(n_rho_samples: i32, n_radius_samples: i32) -> Self {
        BssrdfTable {
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

pub fn fresnel_moment1(eta: Float) -> Float {
    let eta2: Float = eta * eta;
    let eta3: Float = eta2 * eta;
    let eta4: Float = eta3 * eta;
    let eta5: Float = eta4 * eta;
    if eta < 1.0 as Float {
        0.45966 as Float - 1.73965 as Float * eta + 3.37668 as Float * eta2 - 3.904945 * eta3
            + 2.49277 as Float * eta4
            - 0.68441 as Float * eta5
    } else {
        -4.61686 as Float + 11.1136 as Float * eta - 10.4646 as Float * eta2
            + 5.11455 as Float * eta3
            - 1.27198 as Float * eta4
            + 0.12746 as Float * eta5
    }
}

pub fn fresnel_moment2(eta: Float) -> Float {
    let eta2: Float = eta * eta;
    let eta3: Float = eta2 * eta;
    let eta4: Float = eta3 * eta;
    let eta5: Float = eta4 * eta;
    if eta < 1.0 as Float {
        0.27614 as Float - 0.87350 as Float * eta + 1.12077 as Float * eta2
            - 0.65095 as Float * eta3
            + 0.07883 as Float * eta4
            + 0.04860 as Float * eta5
    } else {
        let r_eta = 1.0 as Float / eta;
        let r_eta2 = r_eta * r_eta;
        let r_eta3 = r_eta2 * r_eta;
        -547.033 as Float + 45.3087 as Float * r_eta3 - 218.725 as Float * r_eta2
            + 458.843 as Float * r_eta
            + 404.557 as Float * eta
            - 189.519 as Float * eta2
            + 54.9327 as Float * eta3
            - 9.00603 as Float * eta4
            + 0.63942 as Float * eta5
    }
}

pub fn beam_diffusion_ms(sigma_s: Float, sigma_a: Float, g: Float, eta: Float, r: Float) -> Float {
    let n_samples: i32 = 100;
    let mut ed: Float = 0.0;

    // precompute information for dipole integrand

    // compute reduced scattering coefficients $\sigmaps, \sigmapt$
    // and albedo $\rhop$
    let sigmap_s: Float = sigma_s * (1.0 as Float - g);
    let sigmap_t: Float = sigma_a + sigmap_s;
    let rhop: Float = sigmap_s / sigmap_t;
    // compute non-classical diffusion coefficient $D_\roman{G}$ using
    // Equation (15.24)
    let d_g: Float = (2.0 as Float * sigma_a + sigmap_s) / (3.0 as Float * sigmap_t * sigmap_t);
    // compute effective transport coefficient $\sigmatr$ based on $D_\roman{G}$
    let sigma_tr: Float = (sigma_a / d_g).sqrt();
    // determine linear extrapolation distance $\depthextrapolation$
    // using Equation (15.28)
    let fm1: Float = fresnel_moment1(eta);
    let fm2: Float = fresnel_moment2(eta);
    let ze: Float = -2.0 as Float * d_g * (1.0 as Float + 3.0 as Float * fm2)
        / (1.0 as Float - 2.0 as Float * fm1);
    // determine exitance scale factors using Equations (15.31) and (15.32)
    let c_phi: Float = 0.25 as Float * (1.0 as Float - 2.0 as Float * fm1);
    let c_e = 0.5 as Float * (1.0 as Float - 3.0 as Float * fm2);
    // for (int i = 0; i < n_samples; ++i) {
    for i in 0..n_samples {
        // sample real point source depth $\depthreal$
        let zr: Float =
            -(1.0 as Float - (i as Float + 0.5 as Float).ln() / n_samples as Float) / sigmap_t;
        // evaluate dipole integrand $E_{\roman{d}}$ at $\depthreal$ and add to _ed_
        let zv: Float = -zr + 2.0 as Float * ze;
        let dr: Float = (r * r + zr * zr).sqrt();
        let dv: Float = (r * r + zv * zv).sqrt();
        // compute dipole fluence rate $\dipole(r)$ using Equation (15.27)
        let phi_d: Float =
            INV_4_PI / d_g * ((-sigma_tr * dr).exp() / dr - (-sigma_tr * dv).exp() / dv);
        // compute dipole vector irradiance $-\N{}\cdot\dipoleE(r)$
        // using Equation (15.27)
        let ed_n: Float = INV_4_PI
            * (zr * (1.0 as Float + sigma_tr * dr) * (-sigma_tr * dr).exp() / (dr * dr * dr)
                - zv * (1.0 as Float + sigma_tr * dv) * (-sigma_tr * dv).exp() / (dv * dv * dv));
        // add contribution from dipole for depth $\depthreal$ to _ed_
        let e: Float = phi_d * c_phi + ed_n * c_e;
        let kappa: Float = 1.0 as Float - (-2.0 as Float * sigmap_t * (dr + zr)).exp();
        ed += kappa * rhop * rhop * e;
    }
    ed / n_samples as Float
}

pub fn beam_diffusion_ss(sigma_s: Float, sigma_a: Float, g: Float, eta: Float, r: Float) -> Float {
    // compute material parameters and minimum $t$ below the critical angle
    let sigma_t: Float = sigma_a + sigma_s;
    let rho: Float = sigma_s / sigma_t;
    let t_crit: Float = r * (eta * eta - 1.0 as Float).sqrt();
    let mut ess: Float = 0.0 as Float;
    let n_samples: i32 = 100;
    for i in 0..n_samples {
        // evaluate single scattering integrand and add to _ess_
        let ti: Float = t_crit
            - (1.0 as Float - (i as Float + 0.5 as Float) / n_samples as Float).ln() / sigma_t;
        // determine length $d$ of connecting segment and $\cos\theta_\roman{o}$
        let d: Float = (r * r + ti * ti).sqrt();
        let cos_theta_o: Float = ti / d;
        // add contribution of single scattering at depth $t$
        ess += rho * (-sigma_t * (d + t_crit)).exp() / (d * d)
            * phase_hg(cos_theta_o, g)
            * (1.0 as Float - fr_dielectric(-cos_theta_o, 1.0 as Float, eta))
            * (cos_theta_o).abs();
    }
    ess / n_samples as Float
}

pub fn compute_beam_diffusion_bssrdf(g: Float, eta: Float, t: &mut BssrdfTable) {
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
