//std
use std;
use std::f32::consts::PI;
use std::sync::Arc;
// pbrt
use core::geometry::{
    nrm_cross_vec3, nrm_dot_nrm, nrm_dot_vec3, pnt3_distance, vec3_dot_nrm, vec3_dot_vec3,
};
use core::geometry::{Normal3f, Point2f, Point3f, Ray, Vector3f};
use core::interaction::{InteractionCommon, SurfaceInteraction};
use core::interpolation::{catmull_rom_weights, integrate_catmull_rom, sample_catmull_rom_2d};
use core::material::{Material, TransportMode};
use core::medium::phase_hg;
use core::pbrt::clamp_t;
use core::pbrt::INV_4_PI;
use core::pbrt::{Float, Spectrum};
use core::primitive::Primitive;
use core::reflection::{cos_theta, fr_dielectric};
use core::reflection::{Bsdf, Bxdf, BxdfType};
use core::scene::Scene;

pub trait Bssrdf {
    fn s(&self, pi: &SurfaceInteraction, wi: &Vector3f) -> Spectrum;
    fn sample_s(
        &self,
        // the next three (extra) parameters are used for SeparableBssrdfAdapter
        sc: Arc<SeparableBssrdf + Sync + Send>,
        mode: TransportMode,
        eta: Float,
        // done
        scene: &Scene,
        u1: Float,
        u2: &Point2f,
        pdf: &mut Float,
    ) -> (Spectrum, Option<SurfaceInteraction>);
}

pub trait SeparableBssrdf {
    fn sw(&self, w: &Vector3f) -> Spectrum;
    fn sp(&self, pi: &SurfaceInteraction) -> Spectrum;
    fn pdf_sp(&self, pi: &SurfaceInteraction) -> Float;
    fn sample_sp(
        &self,
        scene: &Scene,
        u1: Float,
        u2: &Point2f,
        pi: &mut SurfaceInteraction,
        pdf: &mut Float,
    ) -> Spectrum;
    fn sr(&self, r: Float) -> Spectrum;
    fn pdf_sr(&self, ch: usize, r: Float) -> Float;
    fn sample_sr(&self, ch: usize, u: Float) -> Float;
}

pub struct TabulatedBssrdf {
    // BSSRDF Protected Data
    pub po_p: Point3f,   // pub po: &SurfaceInteraction,
    pub po_time: Float,  // TMP
    pub po_wo: Vector3f, // TMP
    pub eta: Float,
    // SeparableBSSRDF Private Data
    pub ns: Normal3f,
    pub ss: Vector3f,
    pub ts: Vector3f,
    pub material: Arc<Material + Send + Sync>,
    pub mode: TransportMode,
    // TabulatedBSSRDF Private Data
    pub table: Arc<BssrdfTable>,
    pub sigma_t: Spectrum,
    pub rho: Spectrum,
}

impl TabulatedBssrdf {
    pub fn new(
        po: &SurfaceInteraction,
        material_opt: Option<Arc<Material + Send + Sync>>,
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
        let ns: Normal3f = po.shading.n;
        let ss: Vector3f = po.shading.dpdu.normalize();
        if let Some(material) = material_opt {
            TabulatedBssrdf {
                po_p: po.p,
                po_time: po.time,
                po_wo: po.wo,
                eta: eta,
                ns: ns,
                ss: ss,
                ts: nrm_cross_vec3(&ns, &ss),
                material: material,
                mode: mode,
                table: table.clone(),
                sigma_t: sigma_t,
                rho: rho,
            }
        } else {
            panic!("TabulatedBssrdf needs Material pointer")
        }
    }
}

impl Bssrdf for TabulatedBssrdf {
    fn s(&self, pi: &SurfaceInteraction, wi: &Vector3f) -> Spectrum {
        // ProfilePhase pp(Prof::BSSRDFEvaluation);
        let ft: Float = fr_dielectric(cos_theta(&self.po_wo), 1.0 as Float, self.eta);
        self.sp(pi) * self.sw(wi) * (1.0 as Float - ft)
    }
    fn sample_s(
        &self,
        // the next three (extra) parameters are used for SeparableBssrdfAdapter
        sc: Arc<SeparableBssrdf + Sync + Send>,
        mode: TransportMode,
        eta: Float,
        // done
        scene: &Scene,
        u1: Float,
        u2: &Point2f,
        pdf: &mut Float,
    ) -> (Spectrum, Option<SurfaceInteraction>) {
        // ProfilePhase pp(Prof::BSSRDFSampling);
        let mut si: SurfaceInteraction = SurfaceInteraction::default();
        let sp: Spectrum = self.sample_sp(scene, u1, u2, &mut si, pdf);
        if !sp.is_black() {
            // initialize material model at sampled surface interaction
            let mut bxdfs: Vec<Arc<Bxdf + Send + Sync>> = Vec::new();
            bxdfs.push(Arc::new(SeparableBssrdfAdapter::new(sc, mode, eta)));
            si.bsdf = Some(Arc::new(Bsdf::new(&si, 1.0, bxdfs)));
            si.wo = Vector3f::from(si.shading.n);
            (sp, Some(si))
        } else {
            (sp, None)
        }
    }
}

impl SeparableBssrdf for TabulatedBssrdf {
    fn sw(&self, w: &Vector3f) -> Spectrum {
        let c: Float = 1.0 as Float - 2.0 as Float * fresnel_moment1(1.0 as Float / self.eta);
        Spectrum::new(
            (1.0 as Float - fr_dielectric(cos_theta(w), 1.0 as Float, self.eta)) / (c * PI),
        )
    }
    fn sp(&self, pi: &SurfaceInteraction) -> Spectrum {
        self.sr(pnt3_distance(&self.po_p, &pi.p))
    }
    fn pdf_sp(&self, pi: &SurfaceInteraction) -> Float {
        // express $\pti-\pto$ and $\bold{n}_i$ with respect to local coordinates at $\pto$
        let d: Vector3f = self.po_p - pi.p;
        let d_local: Vector3f = Vector3f {
            x: vec3_dot_vec3(&self.ss, &d),
            y: vec3_dot_vec3(&self.ts, &d),
            z: nrm_dot_vec3(&self.ns, &d),
        };
        let n_local: Normal3f = Normal3f {
            x: vec3_dot_nrm(&self.ss, &pi.n),
            y: vec3_dot_nrm(&self.ts, &pi.n),
            z: nrm_dot_nrm(&self.ns, &pi.n),
        };
        // compute BSSRDF profile radius under projection along each axis
        let r_proj: [Float; 3] = [
            (d_local.y * d_local.y + d_local.z * d_local.z).sqrt(),
            (d_local.z * d_local.z + d_local.x * d_local.x).sqrt(),
            (d_local.x * d_local.x + d_local.y * d_local.y).sqrt(),
        ];
        // return combined probability from all BSSRDF sampling strategies
        let mut pdf: Float = 0.0;
        let axis_prob: [Float; 3] = [0.25 as Float, 0.25 as Float, 0.5 as Float];
        let ch_prob: Float = 1.0 as Float / 3.0 as Float;
        for axis in 0..3_usize {
            for ch in 0..3_usize {
                pdf += self.pdf_sr(ch, r_proj[axis])
                    * n_local[axis as u8].abs()
                    * ch_prob
                    * axis_prob[axis];
            }
        }
        pdf
    }
    fn sample_sp(
        &self,
        scene: &Scene,
        u1: Float,
        u2: &Point2f,
        pi: &mut SurfaceInteraction,
        pdf: &mut Float,
    ) -> Spectrum {
        // ProfilePhase pp(Prof::BSSRDFEvaluation);
        let mut u1: Float = u1; // shadowing input parameter

        // choose projection axis for BSSRDF sampling
        let vx: Vector3f;
        let vy: Vector3f;
        let vz: Vector3f;
        if u1 < 0.5 as Float {
            vx = self.ss;
            vy = self.ts;
            vz = Vector3f::from(self.ns);
            u1 *= 2.0 as Float;
        } else if u1 < 0.75 as Float {
            // prepare for sampling rays with respect to _self.ss_
            vx = self.ts;
            vy = Vector3f::from(self.ns);
            vz = self.ss;
            u1 = (u1 - 0.5 as Float) * 4.0 as Float;
        } else {
            // prepare for sampling rays with respect to _self.ts_
            vx = Vector3f::from(self.ns);
            vy = self.ss;
            vz = self.ts;
            u1 = (u1 - 0.75 as Float) * 4.0 as Float;
        }
        // choose spectral channel for BSSRDF sampling
        let ch: usize = clamp_t((u1 * 3.0 as Float) as usize, 0_usize, 2_usize);
        u1 = u1 * 3.0 as Float - ch as Float;
        // sample BSSRDF profile in polar coordinates
        let r: Float = self.sample_sr(ch, u2.x);
        if r < 0.0 as Float {
            return Spectrum::default();
        }
        let phi: Float = 2.0 as Float * PI * u2.y;
        // compute BSSRDF profile bounds and intersection height
        let r_max: Float = self.sample_sr(ch, 0.999 as Float);
        if r >= r_max {
            return Spectrum::default();
        }
        let l: Float = 2.0 as Float * (r_max * r_max - r * r).sqrt();
        // compute BSSRDF sampling ray segment
        let mut base: InteractionCommon = InteractionCommon::default();
        base.p = self.po_p + (vx * phi.cos() + vy * phi.sin()) * r - vz * (l * 0.5 as Float);
        base.time = self.po_time;
        let p_target: Point3f = base.p + vz * l;

        // intersect BSSRDF sampling ray against the scene geometry

        // declare _IntersectionChain_ and linked list
        // struct IntersectionChain {
        //     SurfaceInteraction si;
        //     IntersectionChain *next = nullptr;
        // };
        // IntersectionChain *chain = ARENA_ALLOC(arena, IntersectionChain)();

        // accumulate chain of intersections along ray
        // IntersectionChain *ptr = chain;
        let mut chain: Vec<SurfaceInteraction> = Vec::new();
        let mut n_found: usize = 0;
        loop {
            let mut r: Ray = base.spawn_ray_to_pnt(&p_target);
            if r.d == Vector3f::default() {
                break;
            }
            if let Some(si) = scene.intersect(&mut r) {
                // base = ptr->si;
                base.p = si.p;
                base.time = si.time;
                base.p_error = si.p_error;
                base.wo = si.wo;
                base.n = si.n;
                // TODO: si.medium_interface;
                base.medium_interface = None;
                // append admissible intersection to _IntersectionChain_
                if let Some(geo_prim) = si.primitive {
                    if let Some(material) = geo_prim.get_material() {
                        //     if (ptr->si.primitive->GetMaterial() == this->material) {
                        if Arc::ptr_eq(&material, &self.material) {
                            //         IntersectionChain *next = ARENA_ALLOC(arena, IntersectionChain)();
                            //         ptr->next = next;
                            //         ptr = next;
                            chain.push(si.clone());
                            n_found += 1;
                        }
                    }
                }
            } else {
                break;
            }
        }

        // randomly choose one of several intersections during BSSRDF sampling
        if n_found == 0_usize {
            return Spectrum::default();
        }
        let selected: usize = clamp_t(
            (u1 * n_found as Float) as usize,
            0_usize,
            (n_found - 1) as usize,
        );
        // while (selected-- > 0) chain = chain->next;
        // *pi = chain->si;
        let selected_si: &SurfaceInteraction = &chain[selected];
        pi.p = selected_si.p;
        pi.time = selected_si.time;
        pi.p_error = selected_si.p_error;
        pi.wo = selected_si.wo;
        pi.n = selected_si.n;
        if let Some(ref medium_interface) = selected_si.medium_interface {
            pi.medium_interface = Some(medium_interface.clone());
        } else {
            pi.medium_interface = None;
        }
        pi.uv = selected_si.uv;
        pi.dpdu = selected_si.dpdu;
        pi.dpdv = selected_si.dpdv;
        pi.dndu = selected_si.dndu;
        pi.dndv = selected_si.dndv;
        pi.dpdx = selected_si.dpdx;
        pi.dpdy = selected_si.dpdy;
        pi.dudx = selected_si.dudx;
        pi.dvdx = selected_si.dvdx;
        pi.dudy = selected_si.dudy;
        pi.dvdy = selected_si.dvdy;
        pi.shading = selected_si.shading;
        // no primitive!
        if let Some(ref bsdf) = selected_si.bsdf {
            pi.bsdf = Some(bsdf.clone());
        } else {
            pi.bsdf = None;
        }
        if let Some(ref bssrdf) = selected_si.bssrdf {
            pi.bssrdf = Some(bssrdf.clone());
        } else {
            pi.bssrdf = None;
        }
        // no shape!
        // compute sample PDF and return the spatial BSSRDF term $\sp$
        *pdf = self.pdf_sp(selected_si) / n_found as Float;
        self.sp(selected_si)
    }
    fn sr(&self, r: Float) -> Spectrum {
        let mut sr: Spectrum = Spectrum::default();
        for ch in 0..3_usize {
            // convert $r$ into unitless optical radius $r_{\roman{optical}}$
            let r_optical: Float = r * self.sigma_t.c[ch];
            // compute spline weights to interpolate BSSRDF on channel _ch_
            let mut rho_offset: i32 = 0;
            let mut radius_offset: i32 = 0;
            let mut rho_weights: [Float; 4] = [0.0 as Float; 4];
            let mut radius_weights: [Float; 4] = [0.0 as Float; 4];
            if !catmull_rom_weights(
                &self.table.rho_samples,
                self.rho.c[ch],
                &mut rho_offset,
                &mut rho_weights,
            ) || !catmull_rom_weights(
                &self.table.radius_samples,
                r_optical,
                &mut radius_offset,
                &mut radius_weights,
            ) {
                continue;
            }
            // set BSSRDF value _Sr[ch]_ using tensor spline interpolation
            let mut srf: Float = 0.0;
            for i in 0..4_usize {
                for j in 0..4_usize {
                    let weight: Float = rho_weights[i] * radius_weights[j];
                    if weight != 0.0 as Float {
                        srf += weight
                            * self
                                .table
                                .eval_profile(rho_offset + i as i32, radius_offset + j as i32);
                    }
                }
            }
            // cancel marginal PDF factor from tabulated BSSRDF profile
            if r_optical != 0.0 as Float {
                srf /= 2.0 as Float * PI * r_optical;
            }
            sr.c[ch] = srf;
        }
        // transform BSSRDF value into world space units
        sr *= self.sigma_t * self.sigma_t;
        sr.clamp(0.0 as Float, std::f32::INFINITY as Float)
    }
    fn pdf_sr(&self, ch: usize, r: Float) -> Float {
        // convert $r$ into unitless optical radius $r_{\roman{optical}}$
        let r_optical: Float = r * self.sigma_t[ch];
        // compute spline weights to interpolate BSSRDF density on channel _ch_
        let mut rho_offset: i32 = 0;
        let mut radius_offset: i32 = 0;
        let mut rho_weights: [Float; 4] = [0.0 as Float; 4];
        let mut radius_weights: [Float; 4] = [0.0 as Float; 4];
        if !catmull_rom_weights(
            &self.table.rho_samples,
            self.rho.c[ch],
            &mut rho_offset,
            &mut rho_weights,
        ) || !catmull_rom_weights(
            &self.table.radius_samples,
            r_optical,
            &mut radius_offset,
            &mut radius_weights,
        ) {
            return 0.0 as Float;
        }
        // return BSSRDF profile density for channel _ch_
        let mut sr: Float = 0.0;
        let mut rho_eff: Float = 0.0;
        for i in 0..4_usize {
            if rho_weights[i] == 0.0 as Float {
                continue;
            }
            rho_eff += self.table.rho_eff[rho_offset as usize + i] * rho_weights[i];
            for j in 0..4_usize {
                if radius_weights[j] == 0.0 as Float {
                    continue;
                }
                sr += self
                    .table
                    .eval_profile(rho_offset + i as i32, radius_offset + j as i32)
                    * rho_weights[i]
                    * radius_weights[j];
            }
        }
        // cancel marginal PDF factor from tabulated BSSRDF profile
        if r_optical != 0.0 as Float {
            sr /= 2.0 as Float * PI * r_optical;
        }
        (0.0 as Float).max(sr * self.sigma_t[ch] * self.sigma_t[ch] / rho_eff)
    }
    fn sample_sr(&self, ch: usize, u: Float) -> Float {
        if self.sigma_t[ch] == 0.0 as Float {
            return -1.0 as Float;
        }
        sample_catmull_rom_2d(
            &self.table.rho_samples,
            &self.table.radius_samples,
            &self.table.profile,
            &self.table.profile_cdf,
            self.rho[ch],
            u,
            None,
            None,
        ) / self.sigma_t[ch]
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
        // initialize all Vec<Float> vectors to zero
        let rho_samples: Vec<Float> = vec![0.0 as Float; n_rho_samples as usize];
        let radius_samples: Vec<Float> = vec![0.0 as Float; n_radius_samples as usize];
        let profile: Vec<Float> = vec![0.0 as Float; (n_radius_samples * n_rho_samples) as usize];
        let rho_eff: Vec<Float> = vec![0.0 as Float; n_rho_samples as usize];
        let profile_cdf: Vec<Float> =
            vec![0.0 as Float; (n_radius_samples * n_rho_samples) as usize];
        BssrdfTable {
            n_rho_samples: n_rho_samples,
            n_radius_samples: n_radius_samples,
            rho_samples: rho_samples,
            radius_samples: radius_samples,
            profile: profile,
            rho_eff: rho_eff,
            profile_cdf: profile_cdf,
        }
    }
    pub fn eval_profile(&self, rho_index: i32, radius_index: i32) -> Float {
        self.profile[(rho_index * self.n_radius_samples + radius_index) as usize]
    }
}

pub struct SeparableBssrdfAdapter {
    // pub bssrdf: &'b (SeparableBssrdf + Send + Sync),
    pub bssrdf: Arc<SeparableBssrdf + Sync + Send>,
    mode: TransportMode,
    eta2: Float,
}

impl SeparableBssrdfAdapter {
    pub fn new(
        bssrdf: Arc<SeparableBssrdf + Sync + Send>,
        mode: TransportMode,
        eta: Float,
    ) -> Self {
        SeparableBssrdfAdapter {
            bssrdf: bssrdf,
            mode: mode,
            eta2: eta * eta,
        }
    }
}

impl Bxdf for SeparableBssrdfAdapter {
    fn f(&self, _wo: &Vector3f, wi: &Vector3f) -> Spectrum {
        let mut f: Spectrum = self.bssrdf.sw(wi);
        // update BSSRDF transmission term to account for adjoint light transport
        if self.mode == TransportMode::Radiance {
            f *= Spectrum::new(self.eta2);
        }
        f
    }
    fn get_type(&self) -> u8 {
        BxdfType::BsdfDiffuse as u8 | BxdfType::BsdfReflection as u8
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
            -(1.0 as Float - (i as Float + 0.5 as Float) / n_samples as Float).ln() / sigmap_t;
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
        let prev_radius_sample: Float = t.radius_samples[i - 1];
        t.radius_samples[i] = prev_radius_sample * 1.2 as Float;
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
            i * t.n_radius_samples as usize,
            &t.profile,
            &mut t.profile_cdf,
        );
    }
    // }, t.n_rho_samples);
}
