//! As the scene file is parsed, objects are created that represent
//! the lights and geometric primitives in the scene. These are all
//! stored in the **Scene** object.
//!

// std
use std::sync::Arc;
// pbrt
use crate::core::geometry::{Bounds3f, Ray, Vector3f};
use crate::core::interaction::{Interaction, SurfaceInteraction};
use crate::core::light::{Light, LightFlags};
use crate::core::pbrt::{Float, Spectrum};
use crate::core::primitive::Primitive;
use crate::core::sampler::Sampler;

// see scene.h

#[derive(Clone)]
pub struct Scene {
    pub lights: Vec<Arc<dyn Light + Sync + Send>>,
    pub infinite_lights: Vec<Arc<dyn Light + Sync + Send>>,
    pub aggregate: Arc<Primitive>,
    pub world_bound: Bounds3f,
}

impl Scene {
    pub fn new(
        aggregate: Arc<Primitive>,
        lights: Vec<Arc<dyn Light + Sync + Send>>,
    ) -> Self {
        let world_bound: Bounds3f = aggregate.world_bound();
        let scene: Scene = Scene {
            lights: Vec::new(),
            infinite_lights: Vec::new(),
            aggregate: aggregate.clone(),
            world_bound,
        };
        let mut changed_lights = Vec::new();
        let mut infinite_lights = Vec::new();
        for light in lights {
            light.preprocess(&scene);
            changed_lights.push(light.clone());
            let check: u8 = light.get_flags() & LightFlags::Infinite as u8;
            if check == LightFlags::Infinite as u8 {
                infinite_lights.push(light);
            }
        }
        Scene {
            lights: changed_lights,
            infinite_lights,
            aggregate,
            world_bound,
        }
    }
    pub fn world_bound(&self) -> Bounds3f {
        self.world_bound
    }
    pub fn intersect(&self, ray: &mut Ray) -> Option<SurfaceInteraction> {
        // TODO: ++nIntersectionTests;
        assert_ne!(
            ray.d,
            Vector3f {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            }
        );
        self.aggregate.intersect(ray)
    }
    pub fn intersect_p(&self, ray: &mut Ray) -> bool {
        // TODO: ++nShadowTests;
        assert_ne!(
            ray.d,
            Vector3f {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            }
        );
        self.aggregate.intersect_p(ray)
    }
    pub fn intersect_tr(
        &self,
        ray: &mut Ray,
        sampler: &mut Box<dyn Sampler + Send + Sync>,
    ) -> (Option<SurfaceInteraction>, Spectrum) {
        let mut tr: Spectrum = Spectrum::new(1.0 as Float);
        loop {
            // bool hit_surface = Intersect(ray, isect);
            if let Some(isect) = self.intersect(ray) {
                // accumulate beam transmittance for ray segment
                if let Some(ref medium_arc) = ray.medium {
                    tr *= medium_arc.tr(&ray, sampler);
                }
                // initialize next ray segment or terminate transmittance computation
                if let Some(ref primitive) = isect.primitive {
                    if let Some(_material) = primitive.get_material() {
                        return (Some(isect), tr);
                    }
                }
                *ray = isect.spawn_ray(&ray.d);
            } else {
                // accumulate beam transmittance for ray segment
                if let Some(ref medium_arc) = ray.medium {
                    tr *= medium_arc.tr(&ray, sampler);
                }
                return (None, tr);
            }
        }
    }
}
