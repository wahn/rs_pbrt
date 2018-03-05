// std
use std;
use std::sync::Arc;
// pbrt
use core::camera::{Camera, CameraSample};
use core::geometry::{Bounds2i, Normal3f, Point2f, Point3f, Ray, Vector3f};
use core::geometry::{nrm_abs_dot_vec3, pnt3_offset_ray_origin, vec3_abs_dot_nrm, vec3_normalize};
use core::light::{Light, LightFlags, VisibilityTester};
use core::material::TransportMode;
use core::primitive::Primitive;
use core::interaction::{Interaction, InteractionCommon, SurfaceInteraction};
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
    pub light: Option<&'a Arc<Light + Send + Sync>>,
}

impl<'a> EndpointInteraction<'a> {
    pub fn new(p: &Point3f, time: Float) -> Self {
        EndpointInteraction {
            p: *p,
            time,
            ..Default::default()
        }
    }
    pub fn new_interaction_from_camera(
        it: &InteractionCommon,
        camera: &'a Box<Camera + Send + Sync>,
    ) -> Self {
        let mut ei: EndpointInteraction = EndpointInteraction::new(&it.p, it.time);
        ei.camera = Some(camera);
        ei
    }
    pub fn new_camera(camera: &'a Box<Camera + Send + Sync>, ray: &Ray) -> Self {
        let mut ei: EndpointInteraction = EndpointInteraction::new(&ray.o, ray.time);
        ei.camera = Some(camera);
        ei
    }
    pub fn new_light(light: &'a Arc<Light + Send + Sync>, ray: &Ray, nl: &Normal3f) -> Self {
        let mut ei: EndpointInteraction = EndpointInteraction::new(&ray.o, ray.time);
        ei.light = Some(light);
        ei.n = *nl;
        ei
    }
    pub fn new_interaction_from_light(
        it: &InteractionCommon,
        light: &'a Arc<Light + Send + Sync>,
    ) -> Self {
        let mut ei: EndpointInteraction = EndpointInteraction::default();
        ei.p = it.p;
        ei.time = it.time;
        ei.p_error = it.p_error;
        ei.wo = it.wo;
        ei.n = it.n;
        ei.light = Some(light);
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
    pub fn create_camera_from_ray(
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
    pub fn create_camera_from_interaction(
        camera: &'a Box<Camera + Send + Sync>,
        it: &InteractionCommon,
        beta: &Spectrum,
    ) -> Vertex<'a, 'p, 's> {
        Vertex::new(
            VertexType::Camera,
            EndpointInteraction::new_interaction_from_camera(it, camera),
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
        light: &'a Arc<Light + Send + Sync>,
        ray: &Ray,
        nl: &Normal3f,
        le: &Spectrum,
        _pdf: Float,
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
    pub fn ns(&self) -> Normal3f {
        match self.vertex_type {
            VertexType::Medium => {}
            VertexType::Surface => {
                if let Some(ref si) = self.si {
                    return si.shading.n;
                } else {
                    return Normal3f::default();
                }
            }
            VertexType::Light => {
                if let Some(ref ei) = self.ei {
                    return ei.n;
                } else {
                    return Normal3f::default();
                }
            }
            _ => {
                return Normal3f::default();
            }
        }
        Normal3f::default()
    }
    pub fn is_on_surface(&self) -> bool {
        self.ng() != Normal3f::default()
    }
    pub fn f(&self, next: &Vertex, mode: TransportMode) -> Spectrum {
        let mut wi: Vector3f = next.p() - self.p();
        if wi.length_squared() == 0.0 as Float {
            return Spectrum::default();
        }
        wi = vec3_normalize(&wi);
        match self.vertex_type {
            VertexType::Surface => {
                if let Some(ref si) = self.si {
                    if let Some(ref bsdf) = si.bsdf {
                        let bsdf_flags: u8 = BxdfType::BsdfAll as u8;
                        return bsdf.f(&si.wo, &wi, bsdf_flags)
                            * correct_shading_normal(si, &si.wo, &wi, mode);
                    } else {
                        return Spectrum::default();
                    }
                } else {
                    return Spectrum::default();
                }
            }
            VertexType::Medium => {
                // TODO: return mi.phase->p(mi.wo, wi);
                return Spectrum::default();
            }
            _ => {
                panic!("Vertex::f(): Unimplemented");
            }
        }
    }
    pub fn is_connectible(&self) -> bool {
        match self.vertex_type {
            VertexType::Medium => true,
            VertexType::Light => {
                if let Some(ref ei) = self.ei {
                    if let Some(ref light) = ei.light {
                        let check: u8 = light.get_flags() & LightFlags::DeltaDirection as u8;
                        check == 0 as u8
                    } else {
                        false
                    }
                } else {
                    false
                }
            }
            VertexType::Camera => true,
            VertexType::Surface => {
                if let Some(ref si) = self.si {
                    if let Some(ref bsdf) = si.bsdf {
                        let bsdf_flags: u8 = BxdfType::BsdfDiffuse as u8
                            | BxdfType::BsdfGlossy as u8
                            | BxdfType::BsdfReflection as u8
                            | BxdfType::BsdfTransmission as u8;
                        bsdf.num_components(bsdf_flags) > 0
                    } else {
                        false
                    }
                } else {
                    false
                }
            }
        }
    }
    pub fn is_light(&self) -> bool {
        if self.vertex_type == VertexType::Light {
            return true;
        } else {
            if self.vertex_type == VertexType::Surface {
                if let Some(ref si) = self.si {
                    if let Some(primitive) = si.primitive {
                        if let Some(_area_light) = primitive.get_area_light() {
                            return true;
                        }
                    }
                }
            }
        }
        false
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
    pub fn le(&self, scene: &Scene, v: &Vertex) -> Spectrum {
        if !self.is_light() {
            return Spectrum::default();
        }
        let mut w: Vector3f = v.p() - self.p();
        if w.length_squared() == 0.0 as Float {
            return Spectrum::default();
        }
        w = vec3_normalize(&w);
        if self.is_infinite_light() {
            // return emitted radiance for infinite light sources
            let mut le: Spectrum = Spectrum::default();
            for light in &scene.infinite_lights {
                let mut ray: Ray = Ray {
                    o: self.p(),
                    d: -w,
                    t_max: Float::default(),
                    time: Float::default(),
                    differential: None,
                };
                le += light.le(&mut ray);
            }
            return le;
        } else {
            if let Some(ref si) = self.si {
                if let Some(primitive) = si.primitive {
                    if let Some(light) = primitive.get_area_light() {
                        let mut iref: InteractionCommon = InteractionCommon::default();
                        iref.p = si.p;
                        iref.time = si.time;
                        iref.p_error = si.p_error;
                        iref.wo = si.wo;
                        iref.n = si.n;
                        return light.l(&iref, &w);
                    }
                }
            }
        }
        Spectrum::default()
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
    pub fn pdf_light_origin(
        &self,
        scene: &Scene,
        v: &Vertex,
        light_distr: &Arc<Distribution1D>,
    ) -> Float {
        let mut w: Vector3f = v.p() - self.p();
        if w.length_squared() == 0.0 as Float {
            return 0.0 as Float;
        }
        w = vec3_normalize(&w);
        if self.is_infinite_light() {
            // return solid angle density for infinite light sources
            return infinite_light_density(scene, &light_distr, &w);
        } else {
            // return solid angle density for non-infinite light sources
            //         Float pdf_pos, pdf_dir, pdf_choice = 0;
            let mut pdf_pos: Float = 0.0;
            let mut pdf_dir: Float = 0.0;
            let mut pdf_choice: Float = 0.0;
            // get pointer _light_ to the light source at the vertex
            assert!(self.is_light());
            if self.vertex_type == VertexType::Light {
                // a real light source (not geometry emitting light)
                if let Some(ref ei) = self.ei {
                    if let Some(ref light_ref) = ei.light {
                        // find light in light vector
                        for i in 0..scene.lights.len() {
                            let light = &scene.lights[i];
                            // use ** (alloc::arc::Arc<Light> **)
                            let pr = &**light_ref as *const _ as *const usize;
                            let pl = &*light as *const _ as *const usize;
                            if pr == pl {
                                // compute the discrete probability of sampling _light_, _pdf_choice_
                                pdf_choice = light_distr.discrete_pdf(i);
                                light.pdf_le(
                                    &Ray {
                                        o: self.p(),
                                        d: w,
                                        t_max: std::f32::INFINITY,
                                        time: self.time(),
                                        differential: None,
                                    },
                                    &self.ng(),
                                    &mut pdf_pos,
                                    &mut pdf_dir,
                                );
                                break;
                            }
                        }
                        assert!(pdf_choice != 0.0 as Float);
                        return pdf_pos * pdf_choice;
                    }
                }
            } else {
                // area light from primitive
                if let Some(ref si) = self.si {
                    if let Some(primitive) = si.primitive {
                        if let Some(area_light) = primitive.get_area_light() {
                            // find area light in light vector
                            for i in 0..scene.lights.len() {
                                let light = &scene.lights[i];
                                let pa = &*area_light as *const _ as *const usize;
                                let pl = &*light as *const _ as *const usize;
                                if pa == pl {
                                    // compute the discrete probability of sampling _light_, _pdf_choice_
                                    pdf_choice = light_distr.discrete_pdf(i);
                                    light.pdf_le(
                                        &Ray {
                                            o: self.p(),
                                            d: w,
                                            t_max: std::f32::INFINITY,
                                            time: self.time(),
                                            differential: None,
                                        },
                                        &self.ng(),
                                        &mut pdf_pos,
                                        &mut pdf_dir,
                                    );
                                    break;
                                }
                            }
                            assert!(pdf_choice != 0.0 as Float);
                            return pdf_pos * pdf_choice;
                        }
                    }
                }
            }
        }
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
    path: &mut Vec<Vertex<'a, 'a, 'a>>,
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
    let vertex: Vertex = Vertex::create_camera_from_ray(camera, &ray, &beta);
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
    light_distr: &Arc<Distribution1D>,
    // TODO: light_to_index
    path: &mut Vec<Vertex<'a, 'a, 'a>>,
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
        let vertex: Vertex = Vertex::create_light(light, &ray, &n_light, &le, pdf_pos * light_pdf);
        let is_infinite_light: bool = vertex.is_infinite_light();
        path.push(vertex);
        let mut beta: Spectrum =
            le * nrm_abs_dot_vec3(&n_light, &ray.d) / (light_pdf * pdf_pos * pdf_dir);
        // println!(
        //     "Starting light subpath. Ray: {:?}, Le {:?}, beta {:?}, pdf_pos {:?}, pdf_dir {:?}",
        //     ray, le, beta, pdf_pos, pdf_dir
        // );
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
    path: &mut Vec<Vertex<'a, 'a, 'a>>,
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
        // println!(
        //     "Random walk. Bounces {:?}, beta {:?}, pdf_fwd {:?}, pdf_rev {:?}",
        //     bounces, beta, pdf_fwd, pdf_rev
        // );
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
                // store new vertex
                path.push(vertex);
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
                // println!(
                //     "Random walk sampled dir {:?} f: {:?}, pdf_fwd: {:?}",
                //     wi, f, pdf_fwd
                // );
                if f.is_black() || pdf_fwd == 0.0 as Float {
                    break;
                }
                *beta *= f * vec3_abs_dot_nrm(&wi, &isect.shading.n) / pdf_fwd;
                // println!("Random walk beta now {:?}", beta);
                pdf_rev = bsdf.pdf(&wi, &isect.wo, bsdf_flags);
                if (sampled_type & BxdfType::BsdfSpecular as u8) != 0_u8 {
                    vertex.delta = true;
                    pdf_rev = 0.0 as Float;
                    pdf_fwd = 0.0 as Float;
                }
                *beta *=
                    Spectrum::new(correct_shading_normal(&isect, &isect.wo, &wi, mode.clone()));
                // println!(
                //     "Random walk beta after shading normal correction {:?}",
                //     beta
                // );
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

pub fn connect_bdpt<'a>(
    scene: &'a Scene,
    light_vertices: &'a Vec<Vertex<'a, 'a, 'a>>,
    camera_vertices: &'a Vec<Vertex<'a, 'a, 'a>>,
    s: usize,
    t: usize,
    light_distr: &Arc<Distribution1D>,
    camera: &'a Box<Camera + Send + Sync>,
    sampler: &mut Box<Sampler + Send + Sync>,
    p_raster: &mut Point2f,
) -> Spectrum {
    // TODO: ProfilePhase _(Prof::BDPTConnectSubpaths);
    let mut l: Spectrum = Spectrum::default();
    // ignore invalid connections related to infinite area lights
    if t > 1 && s != 0 && camera_vertices[t - 1].vertex_type == VertexType::Light {
        return Spectrum::default();
    }
    // perform connection and write contribution to _L_
    // Vertex sampled;
    if s == 0 {
        // interpret the camera subpath as a complete path
        if camera_vertices[t - 1].is_light() {
            l = camera_vertices[t - 1].le(scene, &camera_vertices[t - 2])
                * camera_vertices[t - 1].beta;
        }
        assert!(!l.has_nans());
    } else if t == 1 {
        // sample a point on the camera and connect it to the light subpath
        if light_vertices[s - 1].is_connectible() {
            let mut iref: InteractionCommon = InteractionCommon::default();
            // qs.GetInteraction()
            match light_vertices[s - 1].vertex_type {
                VertexType::Medium => {}
                VertexType::Surface => {
                    if let Some(ref si) = light_vertices[s - 1].si {
                        iref.p = si.p;
                        iref.time = si.time;
                        iref.p_error = si.p_error;
                        iref.wo = si.wo;
                        iref.n = si.n;
                    } else {
                    }
                }
                _ => {}
            }
            let mut wi: Vector3f = Vector3f::default();
            let mut pdf: Float = 0.0 as Float;
            let mut vis: VisibilityTester = VisibilityTester::default();
            let wi_color: Spectrum = camera.sample_wi(
                &iref,
                &sampler.get_2d(),
                &mut wi,
                &mut pdf,
                p_raster,
                &mut vis,
            );
            if pdf > 0.0 as Float && !wi_color.is_black() {
                // initialize dynamically sampled vertex and _L_ for $t=1$ case
                let sampled: Vertex =
                    Vertex::create_camera_from_interaction(camera, &vis.p1, &(wi_color / pdf));
                l = light_vertices[s - 1].beta
                    * light_vertices[s - 1].f(&sampled, TransportMode::Importance)
                    * sampled.beta;
                println!("l = {:?}", l);
                if light_vertices[s - 1].is_on_surface() {
                    l *= Spectrum::new(vec3_abs_dot_nrm(&wi, &light_vertices[s - 1].ns()));
                }
                assert!(!l.has_nans());
                // only check visibility after we know that the path
                // would make a non-zero contribution.
                if !l.is_black() {
                    l *= vis.tr(scene, sampler);
                }
            }
        }
    } else if s == 1 {
        // sample a point on a light and connect it to the camera subpath
        if camera_vertices[t - 1].is_connectible() {
            let mut wi: Vector3f = Vector3f::default();
            let mut pdf: Float = 0.0 as Float;
            let mut light_pdf: Option<Float> = Some(0.0 as Float);
            let mut vis: VisibilityTester = VisibilityTester::default();
            let light_num: usize =
                light_distr.sample_discrete(sampler.get_1d(), light_pdf.as_mut());
            //         const std::shared_ptr<Light> &light = scene.lights[light_num];
            let mut iref: InteractionCommon = InteractionCommon::default();
            // qs.GetInteraction()
            match camera_vertices[t - 1].vertex_type {
                VertexType::Medium => {}
                VertexType::Surface => {
                    if let Some(ref si) = camera_vertices[t - 1].si {
                        iref.p = si.p;
                        iref.time = si.time;
                        iref.p_error = si.p_error;
                        iref.wo = si.wo;
                        iref.n = si.n;
                    } else {
                    }
                }
                _ => {}
            }
            let light_weight: Spectrum = scene.lights[light_num].sample_li(
                &iref,
                &sampler.get_2d(),
                &mut wi,
                &mut pdf,
                &mut vis,
            );
            if pdf > 0.0 as Float && !light_weight.is_black() {
                let ei: EndpointInteraction = EndpointInteraction::new_interaction_from_light(
                    &vis.p1,
                    &scene.lights[light_num],
                );
                let sampled: Vertex = Vertex::create_light_interaction(
                    ei,
                    &(light_weight / (pdf * light_pdf.unwrap())),
                    0.0 as Float,
                );
                // sampled.pdfFwd =
                //     sampled.PdfLightOrigin(scene, camera_vertices[t - 1], light_distr, lightToIndex);
                sampled.pdf_light_origin(scene, &camera_vertices[t - 1], light_distr);
                // L = camera_vertices[t - 1].beta * camera_vertices[t - 1].f(sampled, TransportMode::Radiance) * sampled.beta;
                // if (camera_vertices[t - 1].IsOnSurface()) L *= AbsDot(wi, camera_vertices[t - 1].ns());
                // // Only check visibility if the path would carry radiance.
                // if (!L.IsBlack()) L *= vis.Tr(scene, sampler);
            }
        }
    } else {
        // handle all other bidirectional connection cases
        // const Vertex &qs = lightVertices[s - 1], &pt = cameraVertices[t - 1];
        // if (qs.IsConnectible() && pt.IsConnectible()) {
        //     L = qs.beta * qs.f(pt, TransportMode::Importance) *
        //     pt.f(qs, TransportMode::Radiance) * pt.beta;
        //     VLOG(2) << "General connect s: " << s << ", t: " <<
        //         t << " qs: " << qs << ", pt: " << pt << ",
        //         qs.f(pt): " << qs.f(pt,
        //         TransportMode::Importance) << ", pt.f(qs): " <<
        //         pt.f(qs, TransportMode::Radiance) << ", G: " <<
        //         G(scene, sampler, qs, pt) << ", dist^2: " <<
        //         DistanceSquared(qs.p(), pt.p());
        //     if (!L.IsBlack()) L *= G(scene, sampler, qs, pt);
        // }
    }

    // ++totalPaths;
    // if (L.IsBlack()) ++zeroRadiancePaths;
    // ReportValue(pathLength, s + t - 2);

    // // Compute MIS weight for connection strategy
    // Float misWeight =
    //     L.IsBlack() ? 0.f : MISWeight(scene, lightVertices, cameraVertices,
    //                                   sampled, s, t, light_distr, lightToIndex);
    // VLOG(2) << "MIS weight for (s,t) = (" << s << ", " << t << ") connection: "
    //         << misWeight;
    // DCHECK(!std::isnan(misWeight));
    // L *= misWeight;
    // if (misWeightPtr) *misWeightPtr = misWeight;
    // WORK
    l
}

pub fn infinite_light_density<'a>(
    scene: &'a Scene,
    light_distr: &Arc<Distribution1D>,
    // const std::unordered_map<const Light *, size_t> &lightToDistrIndex,
    w: &Vector3f,
) -> Float {
    let mut pdf: Float = 0.0 as Float;
    println!("TODO: infinite_light_density()");
    for light in &scene.infinite_lights {
        // for i in 0..scene.infinite_lights.len() {
        //     CHECK(lightToDistrIndex.find(light.get()) != lightToDistrIndex.end());
        //     size_t index = lightToDistrIndex.find(light.get())->second;
        let index: usize = 0; // TODO: calculate index (see above)
        pdf += light.pdf_li(&SurfaceInteraction::default(), -(*w)) * light_distr.func[index];
    }
    // TODO: Old loop (without cache) !!!
    // for (size_t i = 0; i < scene.lights.size(); ++i)
    //     if (scene.lights[i]->flags & (int)LightFlags::Infinite)
    //         pdf +=
    //             scene.lights[i]->Pdf_Li(Interaction(), -w) * light_distr.func[i];
    pdf / (light_distr.func_int * light_distr.count() as Float)
}
