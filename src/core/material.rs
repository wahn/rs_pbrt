//! The abstract **Material** class defines the interface that
//! material implementations must provide.

//std
use std::cell::Cell;
use std::sync::Arc;
// pbrt
use crate::core::geometry::vec3_cross_vec3;
use crate::core::geometry::{Normal3f, Vector2f, Vector3f};
use crate::core::interaction::SurfaceInteraction;
use crate::core::pbrt::{Float, Spectrum};
use crate::core::texture::Texture;
use crate::materials::disney::DisneyMaterial;
use crate::materials::fourier::FourierMaterial;
use crate::materials::glass::GlassMaterial;
use crate::materials::hair::HairMaterial;
use crate::materials::matte::MatteMaterial;
use crate::materials::metal::MetalMaterial;
use crate::materials::mirror::MirrorMaterial;
use crate::materials::mixmat::MixMaterial;
use crate::materials::plastic::PlasticMaterial;
use crate::materials::substrate::SubstrateMaterial;
use crate::materials::subsurface::SubsurfaceMaterial;
use crate::materials::translucent::TranslucentMaterial;
use crate::materials::uber::UberMaterial;

// see material.h

/// Is used to inform non-symetric BSDFs about the transported
/// quantity so that they can correctly switch between the adjoint and
/// non-adjoint forms.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TransportMode {
    Radiance,
    Importance,
}

pub enum Material {
    Disney(Box<DisneyMaterial>),
    Fourier(Box<FourierMaterial>),
    Glass(Box<GlassMaterial>),
    Hair(Box<HairMaterial>),
    Matte(Box<MatteMaterial>),
    Metal(Box<MetalMaterial>),
    Mirror(Box<MirrorMaterial>),
    Mix(Box<MixMaterial>),
    Plastic(Box<PlasticMaterial>),
    Substrate(Box<SubstrateMaterial>),
    Subsurface(Box<SubsurfaceMaterial>),
    Translucent(Box<TranslucentMaterial>),
    Uber(Box<UberMaterial>),
}

/// **Material** defines the interface that material implementations
/// must provide.
impl Material {
    /// The method is given a **SurfaceInteraction** object that
    /// contains geometric properties at an intersection point on the
    /// surface of a shape and is responsible for determining the
    /// reflective properties at the point and initializing some
    /// member variables.
    pub fn compute_scattering_functions(
        &self,
        si: &mut SurfaceInteraction,
        // arena: &mut Arena,
        mode: TransportMode,
        allow_multiple_lobes: bool,
        mat: Option<Arc<Material>>,
        scale: Option<Spectrum>,
    ) {
        match self {
            Material::Disney(material) => {
                material.compute_scattering_functions(si, mode, allow_multiple_lobes, mat, scale)
            }
            Material::Fourier(material) => {
                material.compute_scattering_functions(si, mode, allow_multiple_lobes, mat, scale)
            }
            Material::Glass(material) => {
                material.compute_scattering_functions(si, mode, allow_multiple_lobes, mat, scale)
            }
            Material::Hair(material) => {
                material.compute_scattering_functions(si, mode, allow_multiple_lobes, mat, scale)
            }
            Material::Matte(material) => {
                material.compute_scattering_functions(si, mode, allow_multiple_lobes, mat, scale)
            }
            Material::Metal(material) => {
                material.compute_scattering_functions(si, mode, allow_multiple_lobes, mat, scale)
            }
            Material::Mirror(material) => {
                material.compute_scattering_functions(si, mode, allow_multiple_lobes, mat, scale)
            }
            Material::Mix(material) => {
                material.compute_scattering_functions(si, mode, allow_multiple_lobes, mat, scale)
            }
            Material::Plastic(material) => {
                material.compute_scattering_functions(si, mode, allow_multiple_lobes, mat, scale)
            }
            Material::Substrate(material) => {
                material.compute_scattering_functions(si, mode, allow_multiple_lobes, mat, scale)
            }
            Material::Subsurface(material) => {
                material.compute_scattering_functions(si, mode, allow_multiple_lobes, mat, scale)
            }
            Material::Translucent(material) => {
                material.compute_scattering_functions(si, mode, allow_multiple_lobes, mat, scale)
            }
            Material::Uber(material) => {
                material.compute_scattering_functions(si, mode, allow_multiple_lobes, mat, scale)
            }
        }
    }
    /// Computing the effect of bump mapping at the point being shaded
    /// given a particular displacement texture.
    pub fn bump(d: &Arc<dyn Texture<Float> + Send + Sync>, si: &mut SurfaceInteraction)
    where
        Self: Sized,
    {
        // compute offset positions and evaluate displacement texture
        let mut si_eval: SurfaceInteraction = SurfaceInteraction::default();
        si_eval.common.p = si.common.p;
        si_eval.common.time = si.common.time;
        si_eval.common.p_error = si.common.p_error;
        si_eval.common.wo = si.common.wo;
        si_eval.common.n = si.common.n;
        if let Some(medium_interface) = &si.common.medium_interface {
            Arc::new(medium_interface.clone());
        } else {
            si_eval.common.medium_interface = None
        }
        si_eval.uv = si.uv;
        si_eval.dpdu = si.dpdu;
        si_eval.dpdv = si.dpdv;
        si_eval.dndu = si.dndu;
        si_eval.dndv = si.dndv;
        si_eval.dudx = Cell::new(si.dudx.get());
        si_eval.dvdx = Cell::new(si.dvdx.get());
        si_eval.dudy = Cell::new(si.dudy.get());
        si_eval.dvdy = Cell::new(si.dvdy.get());
        si_eval.dpdx = Cell::new(si.dpdx.get());
        si_eval.dpdy = Cell::new(si.dpdy.get());
        if let Some(primitive) = &si.primitive {
            Arc::new(*primitive);
        } else {
            si_eval.primitive = None
        }
        si_eval.shading.n = si.shading.n;
        si_eval.shading.dpdu = si.shading.dpdu;
        si_eval.shading.dpdv = si.shading.dpdv;
        si_eval.shading.dndu = si.shading.dndu;
        si_eval.shading.dndv = si.shading.dndv;
        if let Some(bsdf) = &si.bsdf {
            Arc::new(bsdf.clone());
        } else {
            si_eval.bsdf = None
        }
        // if let Some(bssrdf) = &si.bssrdf {
        //     Some(Arc::new(bssrdf.clone()));
        // } else {
        //     si_eval.bssrdf = None
        // }
        if let Some(shape) = &si.shape {
            Arc::new(shape);
        } else {
            si_eval.shape = None
        }
        // shift _si_eval_ _du_ in the $u$ direction
        let mut du: Float = 0.5 as Float * (si.dudx.get().abs() + si.dudy.get().abs());
        // The most common reason for du to be zero is for ray that start from
        // light sources, where no differentials are available. In this case,
        // we try to choose a small enough du so that we still get a decently
        // accurate bump value.
        if du == 0.0 as Float {
            du = 0.0005 as Float;
        }
        si_eval.common.p = si.common.p + si.shading.dpdu * du;
        si_eval.uv = si.uv
            + Vector2f {
                x: du,
                y: 0.0 as Float,
            };
        si_eval.common.n = (Normal3f::from(vec3_cross_vec3(&si.shading.dpdu, &si.shading.dpdv))
            + si.dndu * du)
            .normalize();
        let u_displace: Float = d.evaluate(&si_eval);
        // shift _si_eval_ _dv_ in the $v$ direction
        let mut dv: Float = 0.5 as Float * (si.dvdx.get().abs() + si.dvdy.get().abs());
        if dv == 00 as Float {
            dv = 0.0005 as Float;
        }
        si_eval.common.p = si.common.p + si.shading.dpdv * dv;
        si_eval.uv = si.uv
            + Vector2f {
                x: 0.0 as Float,
                y: dv,
            };
        si_eval.common.n = (Normal3f::from(vec3_cross_vec3(&si.shading.dpdu, &si.shading.dpdv))
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
