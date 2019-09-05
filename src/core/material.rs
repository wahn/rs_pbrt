//! The abstract **Material** class defines the interface that
//! material implementations must provide.

//std
use std::sync::{Arc, RwLock};
// pbrt
use crate::core::geometry::vec3_cross_vec3;
use crate::core::geometry::{Normal3f, Vector2f, Vector3f};
use crate::core::interaction::SurfaceInteraction;
use crate::core::pbrt::Float;
use crate::core::texture::Texture;

// see material.h

/// Is used to inform non-symetric BSDFs about the transported
/// quantity so that they can correctly switch between the adjoint and
/// non-adjoint forms.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TransportMode {
    Radiance,
    Importance,
}

/// **Material** defines the interface that material implementations
/// must provide.
pub trait Material {
    /// The method is given a **SurfaceInteraction** object that
    /// contains geometric properties at an intersection point on the
    /// surface of a shape and is responsible for determining the
    /// reflective properties at the point and initializing some
    /// member variables.
    fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        // arena: &mut Arena,
        mode: TransportMode,
        allow_multiple_lobes: bool,
        material: Option<Arc<dyn Material + Send + Sync>>,
    );
    /// Computing the effect of bump mapping at the point being shaded
    /// given a particular displacement texture.
    fn bump(d: &Arc<dyn Texture<Float> + Send + Sync>, si: &mut SurfaceInteraction)
    where
        Self: Sized,
    {
        // compute offset positions and evaluate displacement texture
        let mut si_eval: SurfaceInteraction = SurfaceInteraction::default();
        si_eval.p = si.p.clone();
        si_eval.time = si.time;
        si_eval.p_error = si.p_error.clone();
        si_eval.wo = si.wo.clone();
        si_eval.n = si.n.clone();
        if let Some(medium_interface) = &si.medium_interface {
            Some(Arc::new(medium_interface.clone()));
        } else {
            si_eval.medium_interface = None
        }
        si_eval.uv = si.uv.clone();
        si_eval.dpdu = si.dpdu.clone();
        si_eval.dpdv = si.dpdv.clone();
        si_eval.dndu = si.dndu.clone();
        si_eval.dndv = si.dndv.clone();
        let dudx: Float = *si.dudx.read().unwrap();
        si_eval.dudx = RwLock::new(dudx);
        let dvdx: Float = *si.dvdx.read().unwrap();
        si_eval.dvdx = RwLock::new(dvdx);
        let dudy: Float = *si.dudy.read().unwrap();
        si_eval.dudy = RwLock::new(dudy);
        let dvdy: Float = *si.dvdy.read().unwrap();
        si_eval.dvdy = RwLock::new(dvdy);
        let dpdx: Vector3f = *si.dpdx.read().unwrap();
        si_eval.dpdx = RwLock::new(dpdx);
        let dpdy: Vector3f = *si.dpdy.read().unwrap();
        si_eval.dpdy = RwLock::new(dpdy);
        if let Some(primitive) = &si.primitive {
            Some(Arc::new(primitive.clone()));
        } else {
            si_eval.primitive = None
        }
        si_eval.shading.n = si.shading.n.clone();
        si_eval.shading.dpdu = si.shading.dpdu.clone();
        si_eval.shading.dpdv = si.shading.dpdv.clone();
        si_eval.shading.dndu = si.shading.dndu.clone();
        si_eval.shading.dndv = si.shading.dndv.clone();
        if let Some(bsdf) = &si.bsdf {
            Some(Arc::new(bsdf.clone()));
        } else {
            si_eval.bsdf = None
        }
        if let Some(bssrdf) = &si.bssrdf {
            Some(Arc::new(bssrdf.clone()));
        } else {
            si_eval.bssrdf = None
        }
        if let Some(shape) = &si.shape {
            Some(Arc::new(shape.clone()));
        } else {
            si_eval.shape = None
        }
        // shift _si_eval_ _du_ in the $u$ direction
        let dudx: Float = *si.dudx.read().unwrap();
        let dudy: Float = *si.dudy.read().unwrap();
        let mut du: Float = 0.5 as Float * (dudx.abs() + dudy.abs());
        // The most common reason for du to be zero is for ray that start from
        // light sources, where no differentials are available. In this case,
        // we try to choose a small enough du so that we still get a decently
        // accurate bump value.
        if du == 0.0 as Float {
            du = 0.0005 as Float;
        }
        si_eval.p = si.p + si.shading.dpdu * du;
        si_eval.uv = si.uv
            + Vector2f {
                x: du,
                y: 0.0 as Float,
            };
        si_eval.n = (Normal3f::from(vec3_cross_vec3(&si.shading.dpdu, &si.shading.dpdv))
            + si.dndu * du)
            .normalize();
        let u_displace: Float = d.evaluate(&si_eval);
        // shift _si_eval_ _dv_ in the $v$ direction
        let dvdx: Float = *si.dvdx.read().unwrap();
        let dvdy: Float = *si.dvdy.read().unwrap();
        let mut dv: Float = 0.5 as Float * (dvdx.abs() + dvdy.abs());
        if dv == 00 as Float {
            dv = 0.0005 as Float;
        }
        si_eval.p = si.p + si.shading.dpdv * dv;
        si_eval.uv = si.uv
            + Vector2f {
                x: 0.0 as Float,
                y: dv,
            };
        si_eval.n = (Normal3f::from(vec3_cross_vec3(&si.shading.dpdu, &si.shading.dpdv))
            + si.dndv * dv)
            .normalize();
        let v_displace: Float = d.evaluate(&si_eval);
        let displace: Float = d.evaluate(&si);
        // compute bump-mapped differential geometry
        let dpdu: Vector3f = si.shading.dpdu
            + Vector3f::from(si.shading.n) * ((u_displace - displace) / du)
            + Vector3f::from(si.shading.dndu) * displace;
        let dpdv: Vector3f = si.shading.dpdv
            + Vector3f::from(si.shading.n) * ((v_displace - displace) / dv)
            + Vector3f::from(si.shading.dndv) * displace;
        let dndu = si.shading.dndu;
        let dndv = si.shading.dndv;
        si.set_shading_geometry(&dpdu, &dpdv, &dndu, &dndv, false);
    }
}
