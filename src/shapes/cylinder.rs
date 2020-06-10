// std
use std::f32::consts::PI;
use std::sync::Arc;
// pbrt
use crate::core::efloat::quadratic_efloat;
use crate::core::efloat::EFloat;
use crate::core::geometry::{
    nrm_abs_dot_vec3, pnt3_distance_squared, vec3_cross_vec3, vec3_dot_vec3,
};
use crate::core::geometry::{Bounds3f, Normal3f, Point2f, Point3f, Ray, Vector3f, XYEnum};
use crate::core::interaction::{Interaction, InteractionCommon, SurfaceInteraction};
use crate::core::material::Material;
use crate::core::pbrt::Float;
use crate::core::pbrt::{clamp_t, gamma, lerp, radians};
use crate::core::transform::Transform;

// see cylinder.h

#[derive(Clone)]
pub struct Cylinder {
    pub radius: Float,
    pub z_min: Float,
    pub z_max: Float,
    pub phi_max: Float,
    // inherited from class Shape (see shape.h)
    pub object_to_world: Transform,
    pub world_to_object: Transform,
    pub reverse_orientation: bool,
    pub transform_swaps_handedness: bool,
    pub material: Option<Arc<Material>>,
}

impl Default for Cylinder {
    fn default() -> Self {
        let object_to_world: Transform = Transform::default();
        Cylinder {
            // Shape
            object_to_world,
            world_to_object: Transform::default(),
            reverse_orientation: false,
            transform_swaps_handedness: object_to_world.swaps_handedness(),
            // Cylinder
            radius: 1.0,
            z_min: -1.0,
            z_max: 1.0,
            phi_max: radians(360.0),
            material: None,
        }
    }
}

impl Cylinder {
    pub fn new(
        object_to_world: Transform,
        world_to_object: Transform,
        reverse_orientation: bool,
        radius: Float,
        z_min: Float,
        z_max: Float,
        phi_max: Float,
    ) -> Self {
        Cylinder {
            // Shape
            object_to_world,
            world_to_object,
            reverse_orientation,
            transform_swaps_handedness: object_to_world.swaps_handedness(),
            // Cylinder
            radius,
            z_min: z_min.min(z_max),
            z_max: z_min.max(z_max),
            phi_max: radians(clamp_t(phi_max, 0.0, 360.0)),
            material: None,
        }
    }
    // Shape
    pub fn object_bound(&self) -> Bounds3f {
        Bounds3f {
            p_min: Point3f {
                x: -self.radius,
                y: -self.radius,
                z: self.z_min,
            },
            p_max: Point3f {
                x: self.radius,
                y: self.radius,
                z: self.z_max,
            },
        }
    }
    pub fn world_bound(&self) -> Bounds3f {
        // in C++: Bounds3f Shape::WorldBound() const { return (*ObjectToWorld)(ObjectBound()); }
        self.object_to_world.transform_bounds(&self.object_bound())
    }
    pub fn intersect(&self, r: &Ray, t_hit: &mut Float, isect: &mut SurfaceInteraction) -> bool {
        // TODO: ProfilePhase p(Prof::ShapeIntersect);
        // transform _Ray_ to object space
        let mut o_err: Vector3f = Vector3f::default();
        let mut d_err: Vector3f = Vector3f::default();
        let ray: Ray = self
            .world_to_object
            .transform_ray_with_error(r, &mut o_err, &mut d_err);

        // compute quadratic cylinder coefficients

        // initialize _EFloat_ ray coordinate values
        let ox = EFloat::new(ray.o.x as f32, o_err.x as f32);
        let oy = EFloat::new(ray.o.y as f32, o_err.y as f32);
        // let oz = EFloat::new(ray.o.z as f32, o_err.z as f32);
        let dx = EFloat::new(ray.d.x as f32, d_err.x as f32);
        let dy = EFloat::new(ray.d.y as f32, d_err.y as f32);
        // let dz = EFloat::new(ray.d.z as f32, d_err.z as f32);
        let a: EFloat = dx * dx + dy * dy;
        let b: EFloat = (dx * ox + dy * oy) * 2.0f32;
        let c: EFloat = ox * ox + oy * oy
            - EFloat::new(self.radius as f32, 0.0) * EFloat::new(self.radius as f32, 0.0);

        // Solve quadratic equation for _t_ values
        let mut t0: EFloat = EFloat::default();
        let mut t1: EFloat = EFloat::default();
        if !quadratic_efloat(a, b, c, &mut t0, &mut t1) {
            return false;
        }
        // check quadric shape _t0_ and _t1_ for nearest intersection
        if t0.upper_bound() > ray.t_max as f32 || t1.lower_bound() <= 0.0f32 {
            return false;
        }
        let mut t_shape_hit: EFloat = t0;
        if t_shape_hit.lower_bound() <= 0.0f32 {
            t_shape_hit = t1;
            if t_shape_hit.upper_bound() > ray.t_max as f32 {
                return false;
            }
        }
        // compute cylinder hit point and $\phi$
        let mut p_hit: Point3f = ray.position(t_shape_hit.v);
        // refine cylinder intersection point
        let hit_rad: Float = (p_hit.x * p_hit.x + p_hit.y * p_hit.y).sqrt();
        p_hit.x *= self.radius / hit_rad;
        p_hit.y *= self.radius / hit_rad;
        let mut phi: Float = p_hit.y.atan2(p_hit.x);
        if phi < 0.0 as Float {
            phi += 2.0 as Float * PI;
        }
        // test cylinder intersection against clipping parameters
        if p_hit.z < self.z_min || p_hit.z > self.z_max || phi > self.phi_max {
            if t_shape_hit == t1 {
                return false;
            }
            t_shape_hit = t1;
            if t1.upper_bound() > ray.t_max {
                return false;
            }
            // compute cylinder hit point and $\phi$
            p_hit = ray.position(t_shape_hit.v);

            // refine cylinder intersection point
            let hit_rad: Float = (p_hit.x * p_hit.x + p_hit.y * p_hit.y).sqrt();
            p_hit.x *= self.radius / hit_rad;
            p_hit.y *= self.radius / hit_rad;
            phi = p_hit.y.atan2(p_hit.x);
            if phi < 0.0 as Float {
                phi += 2.0 as Float * PI;
            }
            if p_hit.z < self.z_min || p_hit.z > self.z_max || phi > self.phi_max {
                return false;
            }
        }
        // find parametric representation of cylinder hit
        let u: Float = phi / self.phi_max;
        let v: Float = (p_hit.z - self.z_min) / (self.z_max - self.z_min);
        // Compute cylinder $\dpdu$ and $\dpdv$
        let dpdu: Vector3f = Vector3f {
            x: -self.phi_max * p_hit.y,
            y: self.phi_max * p_hit.x,
            z: 0.0,
        };
        let dpdv: Vector3f = Vector3f {
            x: 0.0,
            y: 0.0,
            z: self.z_max - self.z_min,
        };
        // compute cylinder $\dndu$ and $\dndv$
        let d2_p_duu: Vector3f = Vector3f {
            x: p_hit.x,
            y: p_hit.y,
            z: 0.0,
        } * -self.phi_max
            * self.phi_max;
        let d2_p_duv: Vector3f = Vector3f {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };
        let d2_p_dvv: Vector3f = Vector3f {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        };
        // compute coefficients for fundamental forms
        let ec: Float = vec3_dot_vec3(&dpdu, &dpdu);
        let fc: Float = vec3_dot_vec3(&dpdu, &dpdv);
        let gc: Float = vec3_dot_vec3(&dpdv, &dpdv);
        let nc: Vector3f = vec3_cross_vec3(&dpdu, &dpdv).normalize();
        let el: Float = vec3_dot_vec3(&nc, &d2_p_duu);
        let fl: Float = vec3_dot_vec3(&nc, &d2_p_duv);
        let gl: Float = vec3_dot_vec3(&nc, &d2_p_dvv);
        // compute $\dndu$ and $\dndv$ from fundamental form coefficients
        let inv_egf2: Float = 1.0 / (ec * gc - fc * fc);
        let dndu = dpdu * (fl * fc - el * gc) * inv_egf2 + dpdv * (el * fc - fl * ec) * inv_egf2;
        let dndu = Normal3f {
            x: dndu.x,
            y: dndu.y,
            z: dndu.z,
        };
        let dndv = dpdu * (gl * fc - fl * gc) * inv_egf2 + dpdv * (fl * fc - gl * ec) * inv_egf2;
        let dndv = Normal3f {
            x: dndv.x,
            y: dndv.y,
            z: dndv.z,
        };
        // compute error bounds for cylinder intersection
        let p_error: Vector3f = Vector3f {
            x: p_hit.x,
            y: p_hit.y,
            z: 0.0,
        }
        .abs()
            * gamma(3_i32);
        // initialize _SurfaceInteraction_ from parametric information
        let uv_hit: Point2f = Point2f { x: u, y: v };
        let wo: Vector3f = -ray.d;
        *isect = SurfaceInteraction::new(
            &p_hit, &p_error, uv_hit, &wo, &dpdu, &dpdv, &dndu, &dndv, ray.time, None,
        );
        self.object_to_world.transform_surface_interaction(isect);
        // if let Some(ref shape) = si.shape {
        //     isect.shape = Some(shape.clone());
        // }
        // if let Some(primitive) = si.primitive {
        //     isect.primitive = Some(primitive.clone());
        // }
        *t_hit = t_shape_hit.v as Float;
        true
    }
    pub fn intersect_p(&self, r: &Ray) -> bool {
        // TODO: ProfilePhase p(Prof::ShapeIntersect);
        // transform _Ray_ to object space
        let mut o_err: Vector3f = Vector3f::default();
        let mut d_err: Vector3f = Vector3f::default();
        let ray: Ray = self
            .world_to_object
            .transform_ray_with_error(r, &mut o_err, &mut d_err);

        // compute quadratic cylinder coefficients

        // initialize _EFloat_ ray coordinate values
        let ox = EFloat::new(ray.o.x as f32, o_err.x as f32);
        let oy = EFloat::new(ray.o.y as f32, o_err.y as f32);
        // let oz = EFloat::new(ray.o.z as f32, o_err.z as f32);
        let dx = EFloat::new(ray.d.x as f32, d_err.x as f32);
        let dy = EFloat::new(ray.d.y as f32, d_err.y as f32);
        // let dz = EFloat::new(ray.d.z as f32, d_err.z as f32);
        let a: EFloat = dx * dx + dy * dy;
        let b: EFloat = (dx * ox + dy * oy) * 2.0f32;
        let c: EFloat = ox * ox + oy * oy
            - EFloat::new(self.radius as f32, 0.0) * EFloat::new(self.radius as f32, 0.0);

        // Solve quadratic equation for _t_ values
        let mut t0: EFloat = EFloat::default();
        let mut t1: EFloat = EFloat::default();
        if !quadratic_efloat(a, b, c, &mut t0, &mut t1) {
            return false;
        }
        // check quadric shape _t0_ and _t1_ for nearest intersection
        if t0.upper_bound() > ray.t_max as f32 || t1.lower_bound() <= 0.0f32 {
            return false;
        }
        let mut t_shape_hit: EFloat = t0;
        if t_shape_hit.lower_bound() <= 0.0f32 {
            t_shape_hit = t1;
            if t_shape_hit.upper_bound() > ray.t_max as f32 {
                return false;
            }
        }
        // compute cylinder hit point and $\phi$
        let mut p_hit: Point3f = ray.position(t_shape_hit.v);
        // refine cylinder intersection point
        let hit_rad: Float = (p_hit.x * p_hit.x + p_hit.y * p_hit.y).sqrt();
        p_hit.x *= self.radius / hit_rad;
        p_hit.y *= self.radius / hit_rad;
        let mut phi: Float = p_hit.y.atan2(p_hit.x);
        if phi < 0.0 as Float {
            phi += 2.0 as Float * PI;
        }
        // test cylinder intersection against clipping parameters
        if p_hit.z < self.z_min || p_hit.z > self.z_max || phi > self.phi_max {
            if t_shape_hit == t1 {
                return false;
            }
            t_shape_hit = t1;
            if t1.upper_bound() > ray.t_max {
                return false;
            }
            // compute cylinder hit point and $\phi$
            p_hit = ray.position(t_shape_hit.v);

            // refine cylinder intersection point
            let hit_rad: Float = (p_hit.x * p_hit.x + p_hit.y * p_hit.y).sqrt();
            p_hit.x *= self.radius / hit_rad;
            p_hit.y *= self.radius / hit_rad;
            phi = p_hit.y.atan2(p_hit.x);
            if phi < 0.0 as Float {
                phi += 2.0 as Float * PI;
            }
            if p_hit.z < self.z_min || p_hit.z > self.z_max || phi > self.phi_max {
                return false;
            }
        }
        true
    }
    pub fn get_reverse_orientation(&self) -> bool {
        self.reverse_orientation
    }
    pub fn get_transform_swaps_handedness(&self) -> bool {
        self.transform_swaps_handedness
    }
    pub fn get_object_to_world(&self) -> Transform {
        self.object_to_world
    }
    pub fn area(&self) -> Float {
        (self.z_max - self.z_min) * self.radius * self.phi_max
    }
    pub fn sample(&self, u: Point2f, pdf: &mut Float) -> InteractionCommon {
        let z: Float = lerp(u[XYEnum::X], self.z_min, self.z_max);
        let phi: Float = u[XYEnum::Y] * self.phi_max;
        let mut p_obj: Point3f = Point3f {
            x: self.radius * phi.cos(),
            y: self.radius * phi.sin(),
            z,
        };
        let mut it: InteractionCommon = InteractionCommon::default();
        it.n = self
            .object_to_world
            .transform_normal(&Normal3f {
                x: p_obj.x,
                y: p_obj.y,
                z: 0.0,
            })
            .normalize();
        if self.reverse_orientation {
            it.n *= -1.0 as Float;
        }
        // reproject _p_obj_ to cylinder surface and compute _p_obj_error_
        let hit_rad: Float = (p_obj.x * p_obj.x + p_obj.y * p_obj.y).sqrt();
        p_obj.x *= self.radius / hit_rad;
        p_obj.y *= self.radius / hit_rad;
        let p_obj_error: Vector3f = Vector3f {
            x: p_obj.x,
            y: p_obj.y,
            z: 0.0,
        }
        .abs()
            * gamma(3_i32);
        it.p = self.object_to_world.transform_point_with_abs_error(
            &p_obj,
            &p_obj_error,
            &mut it.p_error,
        );
        *pdf = 1.0 as Float / self.area();
        it
    }
    pub fn sample_with_ref_point(
        &self,
        iref: &InteractionCommon,
        u: Point2f,
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
    pub fn pdf_with_ref_point(&self, iref: &dyn Interaction, wi: &Vector3f) -> Float {
        // intersect sample ray with area light geometry
        let ray: Ray = iref.spawn_ray(wi);
        // ignore any alpha textures used for trimming the shape when
        // performing this intersection. Hack for the "San Miguel"
        // scene, where this is used to make an invisible area light.
        let mut t_hit: Float = 0.0;
        let mut isect_light: SurfaceInteraction = SurfaceInteraction::default();
        if self.intersect(&ray, &mut t_hit, &mut isect_light) {
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
