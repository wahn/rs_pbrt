//! The abstract **Primitive** base class is the bridge between the
//! geometry processing and shading subsystems of pbrt.

// std
use std::sync::Arc;
// pbrt
use crate::accelerators::bvh::BVHAccel;
use crate::accelerators::kdtreeaccel::KdTreeAccel;
use crate::core::geometry::nrm_dot_nrmf;
use crate::core::geometry::{Bounds3f, Ray};
use crate::core::interaction::SurfaceInteraction;
use crate::core::light::Light;
use crate::core::material::{Material, TransportMode};
use crate::core::medium::{Medium, MediumInterface};
use crate::core::pbrt::Float;
use crate::core::shape::Shape;
use crate::core::transform::{AnimatedTransform, Transform};

// see primitive.h

pub enum Primitive {
    Geometric(Box<GeometricPrimitive>),
    Transformed(Box<TransformedPrimitive>),
    BVH(Box<BVHAccel>),
    KdTree(Box<KdTreeAccel>),
}

impl Primitive {
    pub fn world_bound(&self) -> Bounds3f {
        match self {
            Primitive::Geometric(primitive) => primitive.world_bound(),
            Primitive::Transformed(primitive) => primitive.world_bound(),
            Primitive::BVH(primitive) => primitive.world_bound(),
            Primitive::KdTree(primitive) => primitive.world_bound(),
        }
    }
    pub fn intersect(&self, ray: &mut Ray, isect: &mut SurfaceInteraction) -> bool {
        match self {
            Primitive::Geometric(primitive) => {
                let hit_surface: bool = primitive.intersect(ray, isect);
                if hit_surface {
                    isect.primitive = Some(self);
                }
                hit_surface
            }
            Primitive::Transformed(primitive) => primitive.intersect(ray, isect),
            Primitive::BVH(primitive) => primitive.intersect(ray, isect),
            Primitive::KdTree(primitive) => primitive.intersect(ray, isect),
        }
    }
    pub fn intersect_p(&self, ray: &Ray) -> bool {
        match self {
            Primitive::Geometric(primitive) => primitive.intersect_p(ray),
            Primitive::Transformed(primitive) => primitive.intersect_p(ray),
            Primitive::BVH(primitive) => primitive.intersect_p(ray),
            Primitive::KdTree(primitive) => primitive.intersect_p(ray),
        }
    }
    pub fn get_area_light(&self) -> Option<Arc<Light>> {
        match self {
            Primitive::Geometric(primitive) => primitive.get_area_light(),
            Primitive::Transformed(primitive) => primitive.get_area_light(),
            Primitive::BVH(primitive) => primitive.get_area_light(),
            Primitive::KdTree(primitive) => primitive.get_area_light(),
        }
    }
    pub fn get_material(&self) -> Option<Arc<Material>> {
        match self {
            Primitive::Geometric(primitive) => primitive.get_material(),
            Primitive::Transformed(primitive) => primitive.get_material(),
            Primitive::BVH(primitive) => primitive.get_material(),
            Primitive::KdTree(primitive) => primitive.get_material(),
        }
    }
    pub fn compute_scattering_functions(
        &self,
        isect: &mut SurfaceInteraction,
        mode: TransportMode,
        allow_multiple_lobes: bool,
    ) {
        if let Some(ref material) = self.get_material() {
            material.compute_scattering_functions(
                isect,
                mode,
                allow_multiple_lobes,
                self.get_material(),
                None,
            );
        }
        assert!(
            nrm_dot_nrmf(&isect.common.n, &isect.shading.n) >= 0.0,
            "n: {:?} dot shading.n: {:?}",
            isect.common.n,
            isect.shading.n
        );
    }
}

#[derive(Clone)]
pub struct GeometricPrimitive {
    pub shape: Arc<Shape>,
    pub material: Option<Arc<Material>>,
    pub area_light: Option<Arc<Light>>,
    pub medium_interface: Option<Arc<MediumInterface>>,
}

impl GeometricPrimitive {
    pub fn new(
        shape: Arc<Shape>,
        material: Option<Arc<Material>>,
        area_light: Option<Arc<Light>>,
        medium_interface: Option<Arc<MediumInterface>>,
    ) -> Self {
        if let Some(area_light) = area_light {
            if let Some(medium_interface) = medium_interface {
                GeometricPrimitive {
                    shape,
                    material,
                    area_light: Some(area_light),
                    medium_interface: Some(medium_interface),
                }
            } else {
                GeometricPrimitive {
                    shape,
                    material,
                    area_light: Some(area_light),
                    medium_interface: None,
                }
            }
        } else if let Some(medium_interface) = medium_interface {
            GeometricPrimitive {
                shape,
                material,
                area_light: None,
                medium_interface: Some(medium_interface),
            }
        } else {
            GeometricPrimitive {
                shape,
                material,
                area_light: None,
                medium_interface: None,
            }
        }
    }
    // Primitive
    pub fn world_bound(&self) -> Bounds3f {
        self.shape.world_bound()
    }
    pub fn intersect(&self, ray: &mut Ray, isect: &mut SurfaceInteraction) -> bool {
        let mut t_hit: Float = 0.0;
        if self.shape.intersect(ray, &mut t_hit, isect) {
            // TODO: isect.primitive
            ray.t_max = t_hit;
            // let it: &SurfaceInteraction = isect_rc.borrow();
            assert!(nrm_dot_nrmf(&isect.common.n, &isect.shading.n) >= 0.0 as Float);
            // initialize _SurfaceInteraction::mediumInterface_ after
            // _Shape_ intersection
            if let Some(ref medium_interface) = self.medium_interface {
                if medium_interface.is_medium_transition() {
                    isect.common.medium_interface = Some(medium_interface.clone());
                } else if let Some(ref medium_arc) = ray.medium {
                    let inside: Option<Arc<Medium>> = Some(medium_arc.clone());
                    let outside: Option<Arc<Medium>> = Some(medium_arc.clone());
                    isect.common.medium_interface =
                        Some(Arc::new(MediumInterface::new(inside, outside)));
                }
                // print!("medium_interface = {{inside = ");
                // if let Some(ref inside) = medium_interface.inside {
                //     print!("{:p} , outside = ", inside);
                // } else {
                //     print!("0x0 , outside = ")
                // }
                // if let Some(ref outside) = medium_interface.outside {
                //     println!("{:p}}}", outside);
                // } else {
                //     println!("0x0}}")
                // }
            }
            true
        } else {
            false
        }
    }
    pub fn intersect_p(&self, r: &Ray) -> bool {
        self.shape.intersect_p(r)
    }
    pub fn get_material(&self) -> Option<Arc<Material>> {
        if let Some(ref material) = self.material {
            Some(material.clone())
        } else {
            None
        }
    }
    pub fn get_area_light(&self) -> Option<Arc<Light>> {
        if let Some(ref area_light) = self.area_light {
            Some(area_light.clone())
        } else {
            None
        }
    }
}

pub struct TransformedPrimitive {
    pub primitive: Arc<Primitive>,
    pub primitive_to_world: AnimatedTransform,
}

impl TransformedPrimitive {
    pub fn new(primitive: Arc<Primitive>, primitive_to_world: AnimatedTransform) -> Self {
        TransformedPrimitive {
            primitive,
            primitive_to_world,
        }
    }
    // Primitive
    pub fn world_bound(&self) -> Bounds3f {
        self.primitive_to_world
            .motion_bounds(&self.primitive.world_bound())
    }
    pub fn intersect(&self, r: &mut Ray, isect: &mut SurfaceInteraction) -> bool {
        // compute _ray_ after transformation by _self.primitive_to_world_
        let mut interpolated_prim_to_world: Transform = Transform::default();
        self.primitive_to_world
            .interpolate(r.time, &mut interpolated_prim_to_world);
        let mut ray: Ray = Transform::inverse(&interpolated_prim_to_world).transform_ray(&*r);
        if self.primitive.intersect(&mut ray, isect) {
            r.t_max = ray.t_max;
            // transform instance's intersection data to world space
            if !interpolated_prim_to_world.is_identity() {
                interpolated_prim_to_world.transform_surface_interaction(isect);
                // let new_isect = interpolated_prim_to_world.transform_surface_interaction(isect);
                // assert!(nrm_dot_nrmf(&new_isect.n, &new_isect.shading.n) >= 0.0 as Float);
                // let mut is: SurfaceInteraction = SurfaceInteraction::new(
                //     &new_isect.p,
                //     &new_isect.p_error,
                //     new_isect.uv,
                //     &new_isect.wo,
                //     &new_isect.dpdu,
                //     &new_isect.dpdv,
                //     &new_isect.dndu,
                //     &new_isect.dndv,
                //     new_isect.time,
                //     None,
                // );
                // // we need to preserve the primitive pointer
                // if let Some(primitive) = isect.primitive {
                //     is.primitive = Some(primitive);
                // }
                // // keep shading (and normal)
                // is.n = new_isect.n;
                // is.shading.n = new_isect.shading.n;
                // is.shading.dpdu = new_isect.shading.dpdu;
                // is.shading.dpdv = new_isect.shading.dpdv;
                // is.shading.dndu = new_isect.shading.dndu;
                // is.shading.dndv = new_isect.shading.dndv;
                return true;
            }
            false
        } else {
            false
        }
    }
    pub fn intersect_p(&self, r: &Ray) -> bool {
        let mut interpolated_prim_to_world: Transform = Transform::default();
        self.primitive_to_world
            .interpolate(r.time, &mut interpolated_prim_to_world);
        interpolated_prim_to_world = Transform::inverse(&interpolated_prim_to_world);
        self.primitive
            .intersect_p(&interpolated_prim_to_world.transform_ray(&*r))
    }
    pub fn get_material(&self) -> Option<Arc<Material>> {
        None
    }
    pub fn get_area_light(&self) -> Option<Arc<Light>> {
        None
    }
}
