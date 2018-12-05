//! The abstract **Material** class defines the interface that
//! material implementations must provide.

//std
use std::sync::Arc;
// pbrt
use core::geometry::vec3_cross_vec3;
use core::geometry::{Normal3f, Vector2f, Vector3f};
use core::interaction::SurfaceInteraction;
use core::pbrt::Float;
use core::texture::Texture;

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
    );
    /// Computing the effect of bump mapping at the point being shaded
    /// given a particular displacement texture.
    fn bump(d: &Arc<Texture<Float> + Send + Sync>, si: &mut SurfaceInteraction)
    where
        Self: Sized,
    {
        // compute offset positions and evaluate displacement texture
        let mut si_eval: SurfaceInteraction = si.clone();
        // shift _si_eval_ _du_ in the $u$ direction
        let mut du: Float = 0.5 as Float * (si.dudx.abs() + si.dudy.abs());
        // The most common reason for du to be zero is for ray that start from
        // light sources, where no differentials are available. In this case,
        // we try to choose a small enough du so that we still get a decently
        // accurate bump value.
        if du == 0.0 as Float {
            du = 0.0005 as Float;
        }
        si_eval.p = si.p + si.shading.dpdu * du;
        si_eval.uv = si.uv + Vector2f {
            x: du,
            y: 0.0 as Float,
        };
        si_eval.n = (Normal3f::from(vec3_cross_vec3(&si.shading.dpdu, &si.shading.dpdv))
            + si.dndu * du)
            .normalize();
        let u_displace: Float = d.evaluate(&si_eval);
        // shift _si_eval_ _dv_ in the $v$ direction
        let mut dv: Float = 0.5 as Float * (si.dvdx.abs() + si.dvdy.abs());
        if dv == 00 as Float {
            dv = 0.0005 as Float;
        }
        si_eval.p = si.p + si.shading.dpdv * dv;
        si_eval.uv = si.uv + Vector2f {
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
