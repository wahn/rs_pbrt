// pbrt
use core::interaction::SurfaceInteraction;

// see material.h

/// Is used to inform non-symetric BSDFs about the transported
/// quantity so that they can correctly switch between the adjoint and
/// non-adjoint forms.
#[derive(PartialEq)]
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
    fn compute_scattering_functions(&self,
                                    si: &mut SurfaceInteraction,
                                    // arena: &mut Arena,
                                    mode: TransportMode,
                                    allow_multiple_lobes: bool);
}
