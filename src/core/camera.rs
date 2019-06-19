//! The abstract **Camera** base class holds generic camera options
//! and defines the interface that all camera implementations must
//! provide.

// std
use std::sync::Arc;
// pbrt
use crate::core::film::Film;
use crate::core::geometry::{Point2f, Ray, Vector3f};
use crate::core::interaction::InteractionCommon;
use crate::core::light::VisibilityTester;
use crate::core::pbrt::{Float, Spectrum};

// see camera.h

pub trait Camera {
    fn generate_ray_differential(&self, sample: &CameraSample, ray: &mut Ray) -> Float;
    fn we(&self, ray: &Ray, p_raster2: Option<&mut Point2f>) -> Spectrum;
    fn pdf_we(&self, ray: &Ray) -> (Float, Float);
    fn sample_wi(
        &self,
        iref: &InteractionCommon,
        u: &Point2f,
        wi: &mut Vector3f,
        pdf: &mut Float,
        p_raster: &mut Point2f,
        vis: &mut VisibilityTester,
    ) -> Spectrum;
    fn get_shutter_open(&self) -> Float;
    fn get_shutter_close(&self) -> Float;
    fn get_film(&self) -> Arc<Film>;
}

#[derive(Debug, Default, Copy, Clone)]
pub struct CameraSample {
    pub p_film: Point2f,
    pub p_lens: Point2f,
    pub time: Float,
}
