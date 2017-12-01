//! # Scene
//!
//! As the scene file is parsed, objects are created that represent
//! the lights and geometric primitives in the scene. These are all
//! stored in the **Scene** object.
//!

// std
use std::sync::Arc;
// pbrt
use accelerators::BVHAccel;
use core::interaction::SurfaceInteraction;
use core::light::{Light, LightFlags};
use core::primitive::Primitive;
use geometry::{Bounds3f, Ray, Vector3f};

// see scene.h

#[derive(Clone)]
pub struct Scene {
    pub lights: Vec<Arc<Light + Sync + Send>>,
    pub infinite_lights: Vec<Arc<Light + Sync + Send>>,
    pub aggregate: Arc<BVHAccel>, // TODO: Primitive,
    pub world_bound: Bounds3f,
}

impl Scene {
    pub fn new(aggregate: Arc<BVHAccel>,
               lights: Vec<Arc<Light + Sync + Send>>)
               -> Self {
        let world_bound: Bounds3f = aggregate.world_bound();
        let scene: Scene = Scene {
            lights: Vec::new(),
            infinite_lights: Vec::new(),
            aggregate: aggregate.clone(),
            world_bound: world_bound,
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
            infinite_lights: infinite_lights,
            aggregate: aggregate,
            world_bound: world_bound,
        }
    }
    pub fn world_bound(&self) -> Bounds3f {
        self.world_bound
    }
    pub fn intersect(&self, ray: &mut Ray) -> Option<SurfaceInteraction> {
        // TODO: ++nIntersectionTests;
        assert_ne!(ray.d,
                   Vector3f {
                       x: 0.0,
                       y: 0.0,
                       z: 0.0,
                   });
        self.aggregate.intersect(ray)
    }
    pub fn intersect_p(&self, ray: &mut Ray) -> bool {
        // TODO: ++nShadowTests;
        assert_ne!(ray.d,
                   Vector3f {
                       x: 0.0,
                       y: 0.0,
                       z: 0.0,
                   });
        self.aggregate.intersect_p(ray)
    }
}
