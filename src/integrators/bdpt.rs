// std
use std;
use std::sync::Arc;
// pbrt
use core::camera::{Camera, CameraSample};
use core::geometry::{Bounds2i, Normal3f, Point2f, Point3f, Ray, Vector3f};
use core::geometry::{nrm_abs_dot_vec3, pnt3_offset_ray_origin, vec3_abs_dot_nrm};
use core::light::{Light, LightFlags};
use core::material::TransportMode;
use core::interaction::{Interaction, SurfaceInteraction};
use core::pbrt::{Float, Spectrum};
use core::reflection::BxdfType;
use core::sampler::Sampler;
use core::sampling::Distribution1D;
use core::scene::Scene;

// see bdpt.h

#[derive(Default)]
pub struct EndpointInteraction<'a> {
    // Interaction Public Data
    pub p: Point3f,
    pub time: Float,
    pub p_error: Vector3f,
    pub wo: Vector3f,
    pub n: Normal3f,
    // EndpointInteraction Public Data
    pub camera: Option<&'a Box<Camera + Send + Sync>>,
    pub light: Option<Arc<Light + Send + Sync>>,
}

impl<'a> EndpointInteraction<'a> {
    pub fn new(p: &Point3f, time: Float) -> Self {
        EndpointInteraction {
            p: *p,
            time,
            ..Default::default()
        }
    }
    pub fn new_camera(camera: &'a Box<Camera + Send + Sync>, ray: &Ray) -> Self {
        let mut ei: EndpointInteraction = EndpointInteraction::new(&ray.o, ray.time);
        ei.camera = Some(camera);
        ei
    }
    pub fn new_light(light: Arc<Light + Send + Sync>, ray: &Ray, nl: &Normal3f) -> Self {
        let mut ei: EndpointInteraction = EndpointInteraction::new(&ray.o, ray.time);
        ei.light = Some(light);
        ei.n = *nl;
        ei
    }
    pub fn new_ray(ray: &Ray) -> Self {
        let mut ei: EndpointInteraction = EndpointInteraction::new(&ray.o, ray.time);
        ei.n = Normal3f::from(-ray.d);
        ei
    }
}

impl<'a> Interaction for EndpointInteraction<'a> {
    fn is_surface_interaction(&self) -> bool {
        self.n != Normal3f::default()
    }
    fn is_medium_interaction(&self) -> bool {
        !self.is_surface_interaction()
    }
    fn spawn_ray(&self, d: &Vector3f) -> Ray {
        let o: Point3f = pnt3_offset_ray_origin(&self.p, &self.p_error, &self.n, d);
        Ray {
            o: o,
            d: *d,
            t_max: std::f32::INFINITY,
            time: self.time,
            differential: None,
        }
    }
    fn get_p(&self) -> Point3f {
        self.p.clone()
    }
    fn get_time(&self) -> Float {
        self.time
    }
    fn get_p_error(&self) -> Vector3f {
        self.p_error.clone()
    }
    fn get_wo(&self) -> Vector3f {
        self.wo.clone()
    }
    fn get_n(&self) -> Normal3f {
        self.n.clone()
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum VertexType {
    Camera,
    Light,
    Surface,
    Medium,
}

pub struct Vertex<'a, 'p, 's> {
    vertex_type: VertexType,
    beta: Spectrum,
    ei: Option<EndpointInteraction<'a>>,
    si: Option<SurfaceInteraction<'p, 's>>,
    delta: bool,
    pdf_fwd: Float,
    pdf_rev: Float,
}

impl<'a, 'p, 's> Vertex<'a, 'p, 's> {
    pub fn new(vertex_type: VertexType, ei: EndpointInteraction<'a>, beta: &Spectrum) -> Self {
        Vertex {
            vertex_type: vertex_type,
            beta: *beta,
            ei: Some(ei),
            si: None,
            delta: false,
            pdf_fwd: 0.0 as Float,
            pdf_rev: 0.0 as Float,
        }
    }
    pub fn create_camera(
        camera: &'a Box<Camera + Send + Sync>,
        ray: &Ray,
        beta: &Spectrum,
    ) -> Vertex<'a, 'p, 's> {
        Vertex::new(
            VertexType::Camera,
            EndpointInteraction::new_camera(camera, ray),
            beta,
        )
    }
    pub fn create_surface_interaction(
        si: SurfaceInteraction<'p, 's>,
        beta: &Spectrum,
        pdf: Float,
        prev: &Vertex,
    ) -> Vertex<'a, 'p, 's> {
        let mut v: Vertex = Vertex {
            vertex_type: VertexType::Surface,
            beta: *beta,
            ei: None,
            si: Some(si),
            delta: false,
            pdf_fwd: 0.0 as Float,
            pdf_rev: 0.0 as Float,
        };
        v.pdf_fwd = prev.convert_density(pdf, &v);
        v
    }
    pub fn create_light_interaction(
        ei: EndpointInteraction<'a>,
        beta: &Spectrum,
        pdf: Float,
    ) -> Vertex<'a, 'p, 's> {
        let mut v: Vertex = Vertex::new(VertexType::Light, ei, beta);
        v.pdf_fwd = pdf;
        v
    }
    pub fn create_light(
        light: Arc<Light + Send + Sync>,
        ray: &Ray,
        nl: &Normal3f,
        le: &Spectrum,
        pdf: Float,
    ) -> Vertex<'a, 'p, 's> {
        Vertex::new(
            VertexType::Light,
            EndpointInteraction::new_light(light, ray, nl),
            le,
        )
    }
    pub fn p(&self) -> Point3f {
        match self.vertex_type {
            VertexType::Medium => Point3f::default(),
            VertexType::Surface => {
                if let Some(ref si) = self.si {
                    si.p
                } else {
                    Point3f::default()
                }
            }
            _ => {
                if let Some(ref ei) = self.ei {
                    ei.p
                } else {
                    Point3f::default()
                }
            }
        }
    }
    pub fn time(&self) -> Float {
        match self.vertex_type {
            VertexType::Medium => Float::default(),
            VertexType::Surface => {
                if let Some(ref si) = self.si {
                    si.time
                } else {
                    Float::default()
                }
            }
            _ => {
                if let Some(ref ei) = self.ei {
                    ei.time
                } else {
                    Float::default()
                }
            }
        }
    }
    pub fn ng(&self) -> Normal3f {
        match self.vertex_type {
            VertexType::Medium => Normal3f::default(),
            VertexType::Surface => {
                if let Some(ref si) = self.si {
                    si.n
                } else {
                    Normal3f::default()
                }
            }
            _ => {
                if let Some(ref ei) = self.ei {
                    ei.n
                } else {
                    Normal3f::default()
                }
            }
        }
    }
    pub fn is_on_surface(&self) -> bool {
        self.ng() != Normal3f::default()
    }
    pub fn is_infinite_light(&self) -> bool {
        if self.vertex_type != VertexType::Light {
            return false;
        } else if let Some(ref ei) = self.ei {
            if let Some(ref light) = ei.light {
                let check: u8 = light.get_flags() & LightFlags::Infinite as u8;
                if check == LightFlags::Infinite as u8 {
                    return true;
                }
                let check: u8 = light.get_flags() & LightFlags::DeltaDirection as u8;
                if check == LightFlags::DeltaDirection as u8 {
                    return true;
                }
            }
        }
        false
    }
    pub fn convert_density(&self, pdf: Float, next: &Vertex) -> Float {
        // return solid angle density if _next_ is an infinite area light
        if next.is_infinite_light() {
            return pdf;
        }
        let w: Vector3f = next.p() - self.p();
        if w.length_squared() == 0.0 as Float {
            return 0.0 as Float;
        }
        let inv_dist_2: Float = 1.0 as Float / w.length_squared();
        let mut pdf: Float = pdf; // shadow
        if next.is_on_surface() {
            pdf *= nrm_abs_dot_vec3(&next.ng(), &(w * inv_dist_2.sqrt()));
        }
        pdf * inv_dist_2
    }
}

/// Bidirectional Path Tracing (Global Illumination)
pub struct BDPTIntegrator {
    pub max_depth: u32,
    visualize_strategies: bool,
    visualize_weights: bool,
    pub pixel_bounds: Bounds2i,
    light_sample_strategy: String, // "power"
}

impl BDPTIntegrator {
    pub fn new(
        // TODO: sampler
        // TODO: camera
        max_depth: u32,
        visualize_strategies: bool,
        visualize_weights: bool,
        pixel_bounds: Bounds2i,
        light_sample_strategy: String,
    ) -> Self {
        BDPTIntegrator {
            max_depth: max_depth,
            visualize_strategies: visualize_strategies,
            visualize_weights: visualize_weights,
            pixel_bounds: pixel_bounds,
            light_sample_strategy: light_sample_strategy,
        }
    }
    pub fn get_light_sample_strategy(&self) -> String {
        self.light_sample_strategy.clone()
    }
}

// BDPT Utility Functions

pub fn correct_shading_normal(
    isect: &SurfaceInteraction,
    wo: &Vector3f,
    &wi: &Vector3f,
    mode: TransportMode,
) -> Float {
    if mode == TransportMode::Importance {
        let num: Float = vec3_abs_dot_nrm(&wo, &isect.shading.n) * vec3_abs_dot_nrm(&wi, &isect.n);
        let denom: Float =
            vec3_abs_dot_nrm(&wo, &isect.n) * vec3_abs_dot_nrm(&wi, &isect.shading.n);
        // wi is occasionally perpendicular to isect.shading.n; this
        // is fine, but we don't want to return an infinite or NaN
        // value in that case.
        if denom == 0.0 as Float {
            return 0.0 as Float;
        }
        num / denom
    } else {
        1.0 as Float
    }
}

pub fn generate_camera_subpath<'a>(
    scene: &'a Scene,
    sampler: &mut Box<Sampler + Send + Sync>,
    max_depth: u32,
    camera: &'a Box<Camera + Send + Sync>,
    p_film: &Point2f,
    path: &'a mut Vec<Vertex<'a, 'a, 'a>>,
) -> (usize, Point3f, Float) {
    if max_depth == 0 {
        return (0_usize, Point3f::default(), Float::default());
    }
    // TODO: ProfilePhase _(Prof::BDPTGenerateSubpath);
    // sample initial ray for camera subpath
    let mut camera_sample: CameraSample = CameraSample::default();
    camera_sample.p_film = *p_film;
    camera_sample.time = sampler.get_1d();
    camera_sample.p_lens = sampler.get_2d();
    let mut ray: Ray = Ray::default();
    let mut beta: Spectrum =
        Spectrum::new(camera.generate_ray_differential(&camera_sample, &mut ray));
    ray.scale_differentials(1.0 as Float / (sampler.get_samples_per_pixel() as Float).sqrt());
    // generate first vertex on camera subpath and start random walk
    let vertex: Vertex = Vertex::create_camera(camera, &ray, &beta);
    // get extra info
    let p: Point3f = vertex.p();
    let time: Float = vertex.time();
    // store vertex
    path.push(vertex);
    let (_pdf_pos, pdf_dir) = camera.pdf_we(&ray);
    (
        random_walk(
            scene,
            &mut ray,
            sampler,
            &mut beta,
            pdf_dir,
            max_depth - 1_u32,
            TransportMode::Radiance,
            path,
            None,
        ) + 1_usize,
        p,
        time,
    )
}

pub fn generate_light_subpath<'a>(
    scene: &'a Scene,
    sampler: &mut Box<Sampler + Send + Sync>,
    max_depth: u32,
    time: Float,
    light_distr: Arc<Distribution1D>,
    // TODO: light_to_index
    path: &'a mut Vec<Vertex<'a, 'a, 'a>>,
) -> usize {
    let mut n_vertices: usize = 0_usize;
    if max_depth == 0_u32 {
        return 0_usize;
    }
    // TODO: ProfilePhase _(Prof::BDPTGenerateSubpath);
    // sample initial ray for light subpath
    let mut light_pdf: Option<Float> = Some(0.0 as Float);
    let light_num: usize = light_distr.sample_discrete(sampler.get_1d(), light_pdf.as_mut());
    let ref light = scene.lights[light_num];
    let mut ray: Ray = Ray::default();
    let mut n_light: Normal3f = Normal3f::default();
    let mut pdf_pos: Float = 0.0 as Float;
    let mut pdf_dir: Float = 0.0 as Float;
    let le: Spectrum = light.sample_le(
        &sampler.get_2d(),
        &sampler.get_2d(),
        time,
        &mut ray,
        &mut n_light,
        &mut pdf_pos,
        &mut pdf_dir,
    );
    if pdf_pos == 0.0 as Float || pdf_dir == 0.0 as Float || le.is_black() {
        return 0_usize;
    }
    if let Some(light_pdf) = light_pdf {
        // generate first vertex on light subpath and start random walk
        let vertex: Vertex =
            Vertex::create_light(light.clone(), &ray, &n_light, &le, pdf_pos * light_pdf);
        let is_infinite_light: bool = vertex.is_infinite_light();
        path.push(vertex);
        let mut beta: Spectrum =
            le * nrm_abs_dot_vec3(&n_light, &ray.d) / (light_pdf * pdf_pos * pdf_dir);
        println!(
            "Starting light subpath. Ray: {:?}, Le {:?}, beta {:?}, pdf_pos {:?}, pdf_dir {:?}",
            ray, le, beta, pdf_pos, pdf_dir
        );
        let n_vertices: usize;
        if is_infinite_light {
            n_vertices = random_walk(
                scene,
                &mut ray,
                sampler,
                &mut beta,
                pdf_dir,
                max_depth - 1,
                TransportMode::Importance,
                path,
                Some(pdf_pos),
            );
        } else {
            n_vertices = random_walk(
                scene,
                &mut ray,
                sampler,
                &mut beta,
                pdf_dir,
                max_depth - 1,
                TransportMode::Importance,
                path,
                None,
            );
        }
    }
    n_vertices + 1
}

pub fn random_walk<'a>(
    scene: &'a Scene,
    ray: &mut Ray,
    sampler: &mut Box<Sampler + Send + Sync>,
    beta: &mut Spectrum,
    pdf: Float,
    max_depth: u32,
    mode: TransportMode,
    path: &'a mut Vec<Vertex<'a, 'a, 'a>>,
    density_info: Option<Float>,
) -> usize {
    let mut bounces: usize = 0_usize;
    if max_depth == 0_u32 {
        return bounces;
    }
    // declare variables for forward and reverse probability densities
    let mut pdf_fwd: Float = pdf;
    let mut pdf_rev: Float = 0.0 as Float;
    loop {
        // attempt to create the next subpath vertex in _path_
        println!(
            "Random walk. Bounces {:?}, beta {:?}, pdf_fwd {:?}, pdf_rev {:?}",
            bounces, beta, pdf_fwd, pdf_rev
        );
        // TODO: Handle MediumInteraction
        // trace a ray and sample the medium, if any
        if let Some(mut isect) = scene.intersect(ray) {
            // compute scattering functions for _mode_ and skip over medium
            // boundaries
            isect.compute_scattering_functions(ray /*, arena, */, true, mode.clone());

            // if (!isect.bsdf) {
            //     ray = isect.SpawnRay(ray.d);
            //     continue;
            // }

            // initialize _vertex_ with surface intersection information
            let mut vertex: Vertex = Vertex::create_surface_interaction(
                isect.clone(),
                &beta,
                pdf_fwd,
                &path[bounces as usize],
            );
            bounces += 1;
            if bounces as u32 >= max_depth {
                break;
            }
            if let Some(ref bsdf) = isect.clone().bsdf {
                // sample BSDF at current vertex and compute reverse probability
                let mut wi: Vector3f = Vector3f::default();
                let bsdf_flags: u8 = BxdfType::BsdfAll as u8;
                let mut sampled_type: u8 = u8::max_value(); // != 0
                let f: Spectrum = bsdf.sample_f(
                    &isect.wo,
                    &mut wi,
                    &sampler.get_2d(),
                    &mut pdf_fwd,
                    bsdf_flags,
                    &mut sampled_type,
                );
                println!(
                    "Random walk sampled dir {:?} f: {:?}, pdf_fwd: {:?}",
                    wi, f, pdf_fwd
                );
                if f.is_black() || pdf_fwd == 0.0 as Float {
                    break;
                }
                *beta *= f * vec3_abs_dot_nrm(&wi, &isect.shading.n) / pdf_fwd;
                println!("Random walk beta now {:?}", beta);
                pdf_rev = bsdf.pdf(&wi, &isect.wo, bsdf_flags);
                if (sampled_type & BxdfType::BsdfSpecular as u8) != 0_u8 {
                    vertex.delta = true;
                    pdf_rev = 0.0 as Float;
                    pdf_fwd = 0.0 as Float;
                }
                *beta *=
                    Spectrum::new(correct_shading_normal(&isect, &isect.wo, &wi, mode.clone()));
                println!(
                    "Random walk beta after shading normal correction {:?}",
                    beta
                );
                let new_ray = isect.spawn_ray(&wi);
                *ray = new_ray;
                // compute reverse area density at preceding vertex
                let mut new_pdf_rev;
                {
                    let prev: &Vertex = &path[(bounces - 1) as usize];
                    new_pdf_rev = vertex.convert_density(pdf_rev, prev);
                }
                // update previous vertex
                path[(bounces - 1) as usize].pdf_rev = new_pdf_rev;
                // store new vertex
                path.push(vertex);
            } else {
                let new_ray = isect.spawn_ray(&ray.d);
                *ray = new_ray;
                continue;
            }
        } else {
            // capture escaped rays when tracing from the camera
            if mode.clone() == TransportMode::Radiance {
                let vertex: Vertex = Vertex::create_light_interaction(
                    EndpointInteraction::new_ray(ray),
                    &beta,
                    pdf_fwd,
                );
                // store new vertex
                path.push(vertex);
                bounces += 1;
            }
            break;
        }
    }
    // correct subpath sampling densities for infinite area lights
    if let Some(pdf_pos) = density_info {
        // set spatial density of _path[1]_ for infinite area light
        if bounces > 0 {
            path[1].pdf_fwd = pdf_pos;
            if path[1].is_on_surface() {
                path[1].pdf_fwd *= vec3_abs_dot_nrm(&ray.d, &path[1].ng());
            }
        }
        // set spatial density of _path[0]_ for infinite area light
        // path[0].pdf_fwd = infinite_light_density(scene, light_distr, light_to_index, ray.d);
    }
    bounces
}
