// std
use std::f32::consts::PI;
use std::sync::Arc;
// pbrt
use core::geometry::{nrm_abs_dot_vec3, pnt3_distance_squared};
use core::geometry::{Bounds3f, Normal3f, Point2f, Point3f, Ray, Vector3f};
use core::interaction::{Interaction, InteractionCommon, SurfaceInteraction};
use core::material::Material;
use core::pbrt::Float;
use core::pbrt::{clamp_t, radians};
use core::sampling::concentric_sample_disk;
use core::shape::Shape;
use core::transform::Transform;

// see disk.h

#[derive(Clone)]
pub struct Disk {
    pub height: Float,
    pub radius: Float,
    pub inner_radius: Float,
    pub phi_max: Float,
    // inherited from class Shape (see shape.h)
    object_to_world: Transform,
    world_to_object: Transform,
    reverse_orientation: bool,
    transform_swaps_handedness: bool,
    pub material: Option<Arc<Material + Send + Sync>>,
}

impl Default for Disk {
    fn default() -> Self {
        let object_to_world: Transform = Transform::default();
        Disk {
            // Shape
            object_to_world: object_to_world,
            world_to_object: Transform::default(),
            reverse_orientation: false,
            transform_swaps_handedness: object_to_world.swaps_handedness(),
            // Disk
            height: 0.0,
            radius: 1.0,
            inner_radius: 0.0,
            phi_max: radians(360.0),
            material: None,
        }
    }
}

impl Disk {
    pub fn new(
        object_to_world: Transform,
        world_to_object: Transform,
        reverse_orientation: bool,
        transform_swaps_handedness: bool,
        height: Float,
        radius: Float,
        inner_radius: Float,
        phi_max: Float,
    ) -> Self {
        Disk {
            // Shape
            object_to_world: object_to_world,
            world_to_object: world_to_object,
            reverse_orientation: reverse_orientation,
            transform_swaps_handedness: transform_swaps_handedness,
            // Disk
            height: height,
            radius: radius,
            inner_radius: inner_radius,
            phi_max: radians(clamp_t(phi_max, 0.0, 360.0)),
            material: None,
        }
    }
}

impl Shape for Disk {
    fn object_bound(&self) -> Bounds3f {
        Bounds3f {
            p_min: Point3f {
                x: -self.radius,
                y: -self.radius,
                z: self.height,
            },
            p_max: Point3f {
                x: self.radius,
                y: self.radius,
                z: self.height,
            },
        }
    }
    fn world_bound(&self) -> Bounds3f {
        // in C++: Bounds3f Shape::WorldBound() const { return (*ObjectToWorld)(ObjectBound()); }
        self.object_to_world.transform_bounds(&self.object_bound())
    }
    fn intersect(&self, r: &Ray) -> Option<(SurfaceInteraction, Float)> {
        // TODO: ProfilePhase p(Prof::ShapeIntersect);
        // transform _Ray_ to object space
        let mut o_err: Vector3f = Vector3f::default();
        let mut d_err: Vector3f = Vector3f::default();
        let ray: Ray = self
            .world_to_object
            .transform_ray_with_error(r, &mut o_err, &mut d_err);

        // compute plane intersection for disk

        // reject disk intersections for rays parallel to the disk's plane
        if ray.d.z == 0.0 {
            return None;
        }
        let t_shape_hit: Float = (self.height - ray.o.z) / ray.d.z;
        if t_shape_hit <= 0.0 || t_shape_hit >= ray.t_max {
            return None;
        }
        // see if hit point is inside disk radii and $\phimax$
        let mut p_hit: Point3f = ray.position(t_shape_hit);
        let dist2: Float = p_hit.x * p_hit.x + p_hit.y * p_hit.y;
        if dist2 > self.radius * self.radius || dist2 < self.inner_radius * self.inner_radius {
            return None;
        }
        // test disk $\phi$ value against $\phimax$
        let mut phi: Float = p_hit.y.atan2(p_hit.x);
        if phi < 0.0 {
            phi += 2.0_f32 * PI;
        }
        if phi > self.phi_max {
            return None;
        }
        // find parametric representation of disk hit
        let u: Float = phi / self.phi_max;
        let r_hit: Float = dist2.sqrt();
        let one_minus_v: Float = (r_hit - self.inner_radius) / (self.radius - self.inner_radius);
        let v: Float = 1.0 - one_minus_v;
        let dpdu: Vector3f = Vector3f {
            x: -self.phi_max * p_hit.y,
            y: self.phi_max * p_hit.x,
            z: 0.0,
        };
        let dpdv: Vector3f = Vector3f {
            x: p_hit.x,
            y: p_hit.y,
            z: 0.0,
        } * (self.inner_radius - self.radius)
            / r_hit;
        let dndu: Normal3f = Normal3f::default();
        let dndv: Normal3f = Normal3f::default();
        // refine disk intersection point
        p_hit.z = self.height;
        // compute error bounds for disk intersection
        let p_error: Vector3f = Vector3f::default();
        // initialize _SurfaceInteraction_ from parametric information
        let uv_hit: Point2f = Point2f { x: u, y: v };
        let wo: Vector3f = -ray.d;
        let si: SurfaceInteraction = SurfaceInteraction::new(
            &p_hit, &p_error, &uv_hit, &wo, &dpdu, &dpdv, &dndu, &dndv, ray.time, None,
        );
        let mut isect: SurfaceInteraction = self.object_to_world.transform_surface_interaction(&si);
        if let Some(_shape) = si.shape {
            isect.shape = si.shape;
        }
        if let Some(_primitive) = si.primitive {
            isect.primitive = si.primitive;
        }
        Some((isect, t_shape_hit))
    }
    fn intersect_p(&self, r: &Ray) -> bool {
        // TODO: ProfilePhase p(Prof::ShapeIntersectP);
        // transform _Ray_ to object space
        let mut o_err: Vector3f = Vector3f::default();
        let mut d_err: Vector3f = Vector3f::default();
        let ray: Ray = self
            .world_to_object
            .transform_ray_with_error(r, &mut o_err, &mut d_err);

        // compute plane intersection for disk

        // reject disk intersections for rays parallel to the disk's plane
        if ray.d.z == 0.0 {
            return false;
        }
        let t_shape_hit: Float = (self.height - ray.o.z) / ray.d.z;
        if t_shape_hit <= 0.0 || t_shape_hit >= ray.t_max {
            return false;
        }
        // see if hit point is inside disk radii and $\phimax$
        let p_hit: Point3f = ray.position(t_shape_hit);
        let dist2: Float = p_hit.x * p_hit.x + p_hit.y * p_hit.y;
        if dist2 > self.radius * self.radius || dist2 < self.inner_radius * self.inner_radius {
            return false;
        }
        // test disk $\phi$ value against $\phimax$
        let mut phi: Float = p_hit.y.atan2(p_hit.x);
        if phi < 0.0 {
            phi += 2.0_f32 * PI;
        }
        if phi > self.phi_max {
            return false;
        }
        true
    }
    fn get_reverse_orientation(&self) -> bool {
        self.reverse_orientation
    }
    fn get_transform_swaps_handedness(&self) -> bool {
        self.transform_swaps_handedness
    }
    fn area(&self) -> Float {
        self.phi_max
            * 0.5 as Float
            * (self.radius * self.radius - self.inner_radius * self.inner_radius)
    }
    fn sample(&self, u: &Point2f, pdf: &mut Float) -> InteractionCommon {
        let pd: Point2f = concentric_sample_disk(u);
        let p_obj: Point3f = Point3f {
            x: pd.x * self.radius,
            y: pd.y * self.radius,
            z: self.height,
        };
        let mut it: InteractionCommon = InteractionCommon::default();
        it.n = self
            .object_to_world
            .transform_normal(&Normal3f {
                x: 0.0 as Float,
                y: 0.0 as Float,
                z: 1.0 as Float,
            }).normalize();
        if self.reverse_orientation {
            it.n *= -1.0 as Float;
        }
        let pt_error: Vector3f = Vector3f::default();
        it.p =
            self.object_to_world
                .transform_point_with_abs_error(&p_obj, &pt_error, &mut it.p_error);
        *pdf = 1.0 as Float / self.area();
        it
    }
    fn sample_with_ref_point(
        &self,
        iref: &InteractionCommon,
        u: &Point2f,
        pdf: &mut Float,
    ) -> InteractionCommon {
        let intr: InteractionCommon = self.sample(u, pdf);
        let mut wi: Vector3f = intr.p - iref.p;
        if wi.length_squared() == 0.0 as Float {
            *pdf = 0.0 as Float;
        } else {
            wi = wi.normalize();
            // convert from area measure, as returned by the Sample()
            // call above, to solid angle measure.
            *pdf *= pnt3_distance_squared(&iref.p, &intr.p) / nrm_abs_dot_vec3(&intr.n, &-wi);
            if (*pdf).is_infinite() {
                *pdf = 0.0 as Float;
            }
        }
        intr
    }
    fn pdf_with_ref_point(&self, iref: &Interaction, wi: &Vector3f) -> Float {
        // intersect sample ray with area light geometry
        let ray: Ray = iref.spawn_ray(wi);
        // ignore any alpha textures used for trimming the shape when
        // performing this intersection. Hack for the "San Miguel"
        // scene, where this is used to make an invisible area light.
        if let Some((isect_light, _t_hit)) = self.intersect(&ray) {
            // convert light sample weight to solid angle measure
            let mut pdf: Float = pnt3_distance_squared(&iref.get_p(), &isect_light.p)
                / (nrm_abs_dot_vec3(&isect_light.n, &-(*wi)) * self.area());
            if pdf.is_infinite() {
                pdf = 0.0 as Float;
            }
            pdf
        } else {
            0.0 as Float
        }
    }
}
