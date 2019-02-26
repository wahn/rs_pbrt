//! The abstract **Primitive** base class is the bridge between the
//! geometry processing and shading subsystems of pbrt.

// std
use std::sync::Arc;
// pbrt
use core::geometry::nrm_dot_nrm;
use core::geometry::{Bounds3f, Ray};
use core::interaction::SurfaceInteraction;
use core::light::AreaLight;
use core::material::{Material, TransportMode};
use core::medium::{Medium, MediumInterface};
use core::pbrt::Float;
use core::shape::Shape;
use core::transform::{AnimatedTransform, Transform};

// see primitive.h

pub trait Primitive {
    fn world_bound(&self) -> Bounds3f;
    fn intersect(&self, ray: &mut Ray) -> Option<SurfaceInteraction>;
    fn intersect_p(&self, r: &Ray) -> bool;
    fn get_area_light(&self) -> Option<Arc<AreaLight + Send + Sync>>;
    fn get_material(&self) -> Option<Arc<Material + Send + Sync>>;
    fn compute_scattering_functions(
        &self,
        isect: &mut SurfaceInteraction,
        // arena: &mut Arena,
        mode: TransportMode,
        allow_multiple_lobes: bool,
    ) {
        if let Some(ref material) = self.get_material() {
            material.compute_scattering_functions(
                isect, // arena,
                mode,
                allow_multiple_lobes,
                self.get_material(),
            );
        }
        if isect.shading.n.x.is_nan() {
            println!("DEBUG ME");
        }
        assert!(
            nrm_dot_nrm(&isect.n, &isect.shading.n) > 0.0,
            "n: {:?} dot shading.n: {:?}",
            isect.n,
            isect.shading.n
        );
    }
}

pub struct GeometricPrimitive {
    pub shape: Arc<Shape + Send + Sync>,
    pub material: Option<Arc<Material + Send + Sync>>,
    pub area_light: Option<Arc<AreaLight + Send + Sync>>,
    pub medium_interface: Option<Arc<MediumInterface>>,
}

impl GeometricPrimitive {
    pub fn new(
        shape: Arc<Shape + Send + Sync>,
        material: Option<Arc<Material + Send + Sync>>,
        area_light: Option<Arc<AreaLight + Send + Sync>>,
        medium_interface: Option<Arc<MediumInterface>>,
    ) -> Self {
        if let Some(area_light) = area_light {
            if let Some(medium_interface) = medium_interface {
                GeometricPrimitive {
                    shape: shape,
                    material: material,
                    area_light: Some(area_light),
                    medium_interface: Some(medium_interface),
                }
            } else {
                GeometricPrimitive {
                    shape: shape,
                    material: material,
                    area_light: Some(area_light),
                    medium_interface: None,
                }
            }
        } else {
            if let Some(medium_interface) = medium_interface {
                GeometricPrimitive {
                    shape: shape,
                    material: material,
                    area_light: None,
                    medium_interface: Some(medium_interface),
                }
            } else {
                GeometricPrimitive {
                    shape: shape,
                    material: material,
                    area_light: None,
                    medium_interface: None,
                }
            }
        }
    }
}

impl Primitive for GeometricPrimitive {
    fn world_bound(&self) -> Bounds3f {
        self.shape.world_bound()
    }
    fn intersect(&self, ray: &mut Ray) -> Option<SurfaceInteraction> {
        if let Some((mut isect, t_hit)) = self.shape.intersect(ray) {
            isect.primitive = Some(self.clone());
            ray.t_max = t_hit;
            assert!(nrm_dot_nrm(&isect.n, &isect.shading.n) >= 0.0 as Float);
            // initialize _SurfaceInteraction::mediumInterface_ after
            // _Shape_ intersection
            if let Some(ref medium_interface) = self.medium_interface {
                if medium_interface.is_medium_transition() {
                    isect.medium_interface = Some(medium_interface.clone());
                } else {
                    if let Some(ref medium_arc) = ray.medium {
                        let inside: Option<Arc<Medium + Send + Sync>> = Some(medium_arc.clone());
                        let outside: Option<Arc<Medium + Send + Sync>> = Some(medium_arc.clone());
                        isect.medium_interface =
                            Some(Arc::new(MediumInterface::new(inside, outside)));
                    }
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
            Some(isect)
        } else {
            None
        }
    }
    fn intersect_p(&self, r: &Ray) -> bool {
        self.shape.intersect_p(r)
    }
    fn get_material(&self) -> Option<Arc<Material + Send + Sync>> {
        if let Some(ref material) = self.material {
            Some(material.clone())
        } else {
            None
        }
    }
    fn get_area_light(&self) -> Option<Arc<AreaLight + Send + Sync>> {
        if let Some(ref area_light) = self.area_light {
            Some(area_light.clone())
        } else {
            None
        }
    }
}

pub struct TransformedPrimitive {
    pub primitive: Arc<Primitive + Sync + Send>,
    pub primitive_to_world: AnimatedTransform,
}

impl TransformedPrimitive {
    pub fn new(
        primitive: Arc<Primitive + Sync + Send>,
        primitive_to_world: AnimatedTransform,
    ) -> Self {
        TransformedPrimitive {
            primitive: primitive,
            primitive_to_world: primitive_to_world,
        }
    }
}

impl Primitive for TransformedPrimitive {
    fn world_bound(&self) -> Bounds3f {
        self.primitive_to_world
            .motion_bounds(&self.primitive.world_bound())
    }
    fn intersect(&self, r: &mut Ray) -> Option<SurfaceInteraction> {
        // compute _ray_ after transformation by _self.primitive_to_world_
        let mut interpolated_prim_to_world: Transform = Transform::default();
        self.primitive_to_world
            .interpolate(r.time, &mut interpolated_prim_to_world);
        let mut ray: Ray = Transform::inverse(&interpolated_prim_to_world).transform_ray(&*r);
        if let Some(isect) = self.primitive.intersect(&mut ray) {
            r.t_max = ray.t_max;
            // transform instance's intersection data to world space
            if !interpolated_prim_to_world.is_identity() {
                let new_isect = interpolated_prim_to_world.transform_surface_interaction(&isect);
                assert!(nrm_dot_nrm(&new_isect.n, &new_isect.shading.n) >= 0.0 as Float);
                let mut is: SurfaceInteraction = SurfaceInteraction::new(
                    &new_isect.p,
                    &new_isect.p_error,
                    &new_isect.uv,
                    &new_isect.wo,
                    &new_isect.dpdu,
                    &new_isect.dpdv,
                    &new_isect.dndu,
                    &new_isect.dndv,
                    new_isect.time,
                    None,
                );
                // we need to preserve the primitive pointer
                if let Some(primitive) = isect.primitive {
                    is.primitive = Some(primitive);
                }
                return Some(is);
            }
            None
        } else {
            None
        }
    }
    fn intersect_p(&self, r: &Ray) -> bool {
        let mut interpolated_prim_to_world: Transform = Transform::default();
        self.primitive_to_world
            .interpolate(r.time, &mut interpolated_prim_to_world);
        interpolated_prim_to_world = Transform::inverse(&interpolated_prim_to_world);
        self.primitive
            .intersect_p(&interpolated_prim_to_world.transform_ray(&*r))
    }
    fn get_material(&self) -> Option<Arc<Material + Send + Sync>> {
        None
    }
    fn get_area_light(&self) -> Option<Arc<AreaLight + Send + Sync>> {
        None
    }
}
