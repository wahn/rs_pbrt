// std
use std;
// pbrt
use core::camera::{Camera, CameraSample};
use core::geometry::{Bounds2i, Normal3f, Point2f, Point3f, Ray, Vector3f};
use core::geometry::pnt3_offset_ray_origin;
use core::light::Light;
use core::material::TransportMode;
use core::interaction::{Interaction, SurfaceInteraction};
use core::pbrt::{Float, Spectrum};
use core::sampler::Sampler;
use core::scene::Scene;

// see bdpt.h

#[derive(Default)]
pub struct EndpointInteraction<'c, 'l> {
    // Interaction Public Data
    pub p: Point3f,
    pub time: Float,
    pub p_error: Vector3f,
    pub wo: Vector3f,
    pub n: Normal3f,
    // EndpointInteraction Public Data
    pub camera: Option<&'c Box<Camera + Send + Sync>>,
    pub light: Option<&'l Box<Light + Send + Sync>>,
}

impl<'c, 'l> EndpointInteraction<'c, 'l> {
    pub fn new(p: &Point3f, time: Float) -> Self {
        EndpointInteraction {
            p: *p,
            time,
            ..Default::default()
        }
    }
    pub fn new_camera(camera: &'c Box<Camera + Send + Sync>, ray: &Ray) -> Self {
        let mut ei: EndpointInteraction = EndpointInteraction::new(&ray.o, ray.time);
        ei.camera = Some(camera);
        ei
    }
    pub fn new_ray(ray: &Ray) -> Self {
        let mut ei: EndpointInteraction = EndpointInteraction::new(&ray.o, ray.time);
        ei.n = Normal3f::from(-ray.d);
        ei
    }
}

impl<'c, 'l> Interaction for EndpointInteraction<'c, 'l> {
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

#[derive(Debug, Clone)]
pub enum VertexType {
    Camera,
    Light,
    Surface,
    Medium,
}

pub struct Vertex<'c, 'l, 'p, 's> {
    vertex_type: VertexType,
    beta: Spectrum,
    ei: Option<EndpointInteraction<'c, 'l>>,
    si: Option<SurfaceInteraction<'p, 's>>,
    delta: bool,
    pdf_fwd: Float,
    pdf_rev: Float,
}

impl<'c, 'l, 'p, 's> Vertex<'c, 'l, 'p, 's> {
    pub fn new(vertex_type: VertexType, ei: EndpointInteraction<'c, 'l>, beta: &Spectrum) -> Self {
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
        camera: &'c Box<Camera + Send + Sync>,
        ray: &Ray,
        beta: &Spectrum,
    ) -> Vertex<'c, 'l, 'p, 's> {
        Vertex::new(
            VertexType::Camera,
            EndpointInteraction::new_camera(camera, ray),
            beta,
        )
    }
    pub fn create_surface(
        si: SurfaceInteraction<'p, 's>,
        beta: &Spectrum,
        pdf: Float,
        prev: &Vertex,
    ) -> Vertex<'c, 'l, 'p, 's> {
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
    pub fn create_light(
        ei: EndpointInteraction<'c, 'l>,
        beta: &Spectrum,
        pdf: Float,
    ) -> Vertex<'c, 'l, 'p, 's> {
        let mut v: Vertex = Vertex::new(VertexType::Light, ei, beta);
        v.pdf_fwd = pdf;
        v
    }
    pub fn convert_density(&self, pdf: Float, next: &Vertex) -> Float {
        // WORK
        0.0 as Float
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

pub fn generate_camera_subpath(
    scene: &Scene,
    sampler: &mut Box<Sampler + Send + Sync>,
    max_depth: u32,
    camera: &Box<Camera + Send + Sync>,
    p_film: &Point2f,
) -> u32 {
    if max_depth == 0 {
        return 0_u32;
    }
    // TODO: ProfilePhase _(Prof::BDPTGenerateSubpath);
    // sample initial ray for camera subpath
    let mut camera_sample: CameraSample = CameraSample::default();
    camera_sample.p_film = *p_film;
    camera_sample.time = sampler.get_1d();
    camera_sample.p_lens = sampler.get_2d();
    let mut ray: Ray = Ray::default();
    let beta: Spectrum = Spectrum::new(camera.generate_ray_differential(&camera_sample, &mut ray));
    ray.scale_differentials(1.0 as Float / (sampler.get_samples_per_pixel() as Float).sqrt());
    // generate first vertex on camera subpath and start random walk
    let vertex: Vertex = Vertex::create_camera(camera, &ray, &beta);
    let mut path: Vec<Vertex> = Vec::with_capacity(max_depth as usize);
    path.push(vertex);
    let (pdf_pos, pdf_dir) = camera.pdf_we(&ray);
    random_walk(
        scene,
        &mut ray,
        sampler,
        beta,
        pdf_dir,
        max_depth - 1_u32,
        TransportMode::Radiance,
        &mut path,
    );
    // WORK
    0_u32
}

pub fn random_walk(
    scene: &Scene,
    ray: &mut Ray,
    sampler: &mut Box<Sampler + Send + Sync>,
    beta: Spectrum,
    pdf: Float,
    max_depth: u32,
    mode: TransportMode,
    path: &mut Vec<Vertex>,
) -> u32 {
    let mut bounces: u32 = 0_u32;
    if max_depth == 0_u32 {
        return bounces;
    }
    // declare variables for forward and reverse probability densities
    let pdf_fwd: Float = pdf;
    let pdf_rev: Float = 0.0 as Float;
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
            if let Some(ref bsdf) = isect.bsdf {
                // initialize _vertex_ with surface intersection information
                // vertex = Vertex::CreateSurface(isect, beta, pdf_fwd, prev);

                // TODO: error[E0623]: lifetime mismatch
                // path[bounces as usize] =
                //     Vertex::create_surface(isect, &beta, pdf_fwd, &path[(bounces - 1) as usize]);

                // if (++bounces >= maxDepth) break;
                // sample BSDF at current vertex and compute reverse probability
                // Vector3f wi, wo = isect.wo;
                // BxDFType type;
                // Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf_fwd,
                //                                   BSDF_ALL, &type);
                // VLOG(2) << "Random walk sampled dir " << wi << " f: " << f <<
                //     ", pdf_fwd: " << pdf_fwd;
                // if (f.IsBlack() || pdf_fwd == 0.f) break;
                // beta *= f * AbsDot(wi, isect.shading.n) / pdf_fwd;
                // VLOG(2) << "Random walk beta now " << beta;
                // pdfRev = isect.bsdf->Pdf(wi, wo, BSDF_ALL);
                // if (type & BSDF_SPECULAR) {
                //     vertex.delta = true;
                //     pdfRev = pdf_fwd = 0;
                // }
                // beta *= CorrectShadingNormal(isect, wo, wi, mode.clone());
                // VLOG(2) << "Random walk beta after shading normal correction " << beta;
                // ray = isect.SpawnRay(wi);
                // compute reverse area density at preceding vertex
                // prev.pdfRev = vertex.ConvertDensity(pdfRev, prev);
                // WORK
            } else {
                let new_ray = isect.spawn_ray(&ray.d);
                *ray = new_ray;
                continue;
            }
        } else {
            // capture escaped rays when tracing from the camera
            if mode.clone() == TransportMode::Radiance {
                let vertex: Vertex =
                    Vertex::create_light(EndpointInteraction::new_ray(ray), &beta, pdf_fwd);
                bounces += 1;
            }
            break;
        }
    }
    bounces
}
