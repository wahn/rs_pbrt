// std
use std;
// pbrt
use core::camera::{Camera, CameraSample};
use core::geometry::{Bounds2i, Normal3f, Point2f, Point3f, Ray, Vector3f};
use core::geometry::pnt3_offset_ray_origin;
use core::light::Light;
use core::interaction::Interaction;
use core::pbrt::{Float, Spectrum};
use core::sampler::Sampler;

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
    pub light: Option<Box<Light + Send + Sync>>,
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
        // WORK: How do we store a pointer to a camera here?
        ei.camera = Some(camera);
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

#[derive(Debug, Clone)]
pub enum VertexType {
    Camera,
    Light,
    Surface,
    Medium,
}

pub struct Vertex<'a> {
    vertex_type: VertexType,
    beta: Spectrum,
    ei: Option<EndpointInteraction<'a>>,
    delta: bool,
    pdf_fwd: Float,
    pdf_rev: Float,
}

impl<'a> Vertex<'a> {
    pub fn new(vertex_type: VertexType, ei: EndpointInteraction<'a>, beta: &Spectrum) -> Self {
        Vertex {
            vertex_type: vertex_type,
            beta: *beta,
            ei: Some(ei),
            delta: false,
            pdf_fwd: 0.0 as Float,
            pdf_rev: 0.0 as Float,
        }
    }
    pub fn create_camera(camera: &'a Box<Camera + Send + Sync>, ray: &Ray, beta: &Spectrum) -> Vertex<'a> {
        Vertex::new(VertexType::Camera, EndpointInteraction::new_camera(camera, ray), beta)
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
    // TODO: return random_walk(...) bdpt.cpp:142
    // WORK
    0_u32
}
