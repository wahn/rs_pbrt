//! The abstract **Camera** base class holds generic camera options
//! and defines the interface that all camera implementations must
//! provide.

// std
use std::sync::Arc;
// pbrt
use crate::cameras::environment::EnvironmentCamera;
use crate::cameras::orthographic::OrthographicCamera;
use crate::cameras::perspective::PerspectiveCamera;
use crate::cameras::realistic::RealisticCamera;
use crate::core::film::Film;
use crate::core::geometry::{Point2f, Ray, Vector3f};
use crate::core::interaction::InteractionCommon;
use crate::core::light::VisibilityTester;
use crate::core::pbrt::{Float, Spectrum};

// see camera.h

pub enum Camera {
    Environment(Box<EnvironmentCamera>),
    Orthographic(Box<OrthographicCamera>),
    Perspective(Box<PerspectiveCamera>),
    Realistic(Box<RealisticCamera>),
}

impl Camera {
    pub fn generate_ray_differential(&self, sample: &CameraSample, ray: &mut Ray) -> Float {
        match self {
            Camera::Environment(camera) => camera.generate_ray_differential(sample, ray),
            Camera::Orthographic(camera) => camera.generate_ray_differential(sample, ray),
            Camera::Perspective(camera) => camera.generate_ray_differential(sample, ray),
            Camera::Realistic(camera) => camera.generate_ray_differential(sample, ray),
        }
    }
    pub fn we(&self, ray: &Ray, p_raster2: Option<&mut Point2f>) -> Spectrum {
        match self {
            Camera::Environment(camera) => camera.we(ray, p_raster2),
            Camera::Orthographic(camera) => camera.we(ray, p_raster2),
            Camera::Perspective(camera) => camera.we(ray, p_raster2),
            Camera::Realistic(camera) => camera.we(ray, p_raster2),
        }
    }
    pub fn pdf_we(&self, ray: &Ray) -> (Float, Float) {
        match self {
            Camera::Environment(camera) => camera.pdf_we(ray),
            Camera::Orthographic(camera) => camera.pdf_we(ray),
            Camera::Perspective(camera) => camera.pdf_we(ray),
            Camera::Realistic(camera) => camera.pdf_we(ray),
        }
    }
    pub fn sample_wi<'a, 'b>(
        &self,
        iref: &'a InteractionCommon,
        lens_intr: &'b mut InteractionCommon,
        u: Point2f,
        wi: &mut Vector3f,
        pdf: &mut Float,
        p_raster: &mut Point2f,
        vis: &mut VisibilityTester<'a, 'b>,
    ) -> Spectrum {
        match self {
            Camera::Environment(camera) => {
                camera.sample_wi(iref, lens_intr, u, wi, pdf, p_raster, vis)
            }
            Camera::Orthographic(camera) => {
                camera.sample_wi(iref, lens_intr, u, wi, pdf, p_raster, vis)
            }
            Camera::Perspective(camera) => {
                camera.sample_wi(iref, lens_intr, u, wi, pdf, p_raster, vis)
            }
            Camera::Realistic(camera) => {
                camera.sample_wi(iref, lens_intr, u, wi, pdf, p_raster, vis)
            }
        }
    }
    pub fn get_shutter_open(&self) -> Float {
        match self {
            Camera::Environment(camera) => camera.get_shutter_open(),
            Camera::Orthographic(camera) => camera.get_shutter_open(),
            Camera::Perspective(camera) => camera.get_shutter_open(),
            Camera::Realistic(camera) => camera.get_shutter_open(),
        }
    }
    pub fn get_shutter_close(&self) -> Float {
        match self {
            Camera::Environment(camera) => camera.get_shutter_close(),
            Camera::Orthographic(camera) => camera.get_shutter_close(),
            Camera::Perspective(camera) => camera.get_shutter_close(),
            Camera::Realistic(camera) => camera.get_shutter_close(),
        }
    }
    pub fn get_film(&self) -> Arc<Film> {
        match self {
            Camera::Environment(camera) => camera.get_film(),
            Camera::Orthographic(camera) => camera.get_film(),
            Camera::Perspective(camera) => camera.get_film(),
            Camera::Realistic(camera) => camera.get_film(),
        }
    }
    // extra
    pub fn get_clipping_start(&self) -> Float {
        match self {
            Camera::Perspective(camera) => camera.get_clipping_start(),
            _ => 0.0 as Float,
        }
    }
    pub fn adjust_to_clipping_start(&self, sample: &CameraSample, ray: &mut Ray) {
        if let Camera::Perspective(camera) = self {
            camera.adjust_to_clipping_start(sample, ray)
        }
    }
}

#[derive(Debug, Default, Copy, Clone)]
pub struct CameraSample {
    pub p_film: Point2f,
    pub p_lens: Point2f,
    pub time: Float,
}
