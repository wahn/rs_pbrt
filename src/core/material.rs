// see material.h

/// Is used to inform non-symetric BSDFs about the transported
/// quantity so that they can correctly switch between the adjoint and
/// non-adjoint forms.
#[derive(PartialEq)]
pub enum TransportMode {
    Radiance,
    Importance,
}
