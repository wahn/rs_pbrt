//! The type of film or sensor in a camera has a dramatic effect on
//! the way that incident light is transformed into colors in an
//! image. In **pbrt**, the **Film** class models the sensing device
//! in the simulated camera. After the radiance is found for each
//! camera ray, the **Film** implementation determines the sample's
//! contribution to the pixel around the point on the film plane where
//! the camera ray began and updates its representation of the
//! image. When the main rendering loop exits, the **Film** writes the
//! final image to file.
//!

// std
#[cfg(feature = "openexr")]
use std;
use std::ops::DerefMut;
use std::path::Path;
use std::sync::{Arc, RwLock, RwLockWriteGuard};
// others
use image;
#[cfg(feature = "openexr")]
use openexr::{FrameBuffer, Header, PixelType, ScanlineOutputFile};
// pbrt
use core::filter::Filter;
use core::geometry::{Bounds2f, Bounds2i, Point2f, Point2i, Vector2f};
use core::geometry::{bnd2_intersect_bnd2, pnt2_ceil, pnt2_floor, pnt2_inside_exclusive,
                     pnt2_max_pnt2, pnt2_min_pnt2};
use core::parallel::AtomicFloat;
use core::pbrt::{Float, Spectrum};
use core::pbrt::{clamp_t, gamma_correct};
use core::spectrum::xyz_to_rgb;

// see film.h

const FILTER_TABLE_WIDTH: usize = 16;

#[derive(Debug, Default, Clone)]
pub struct Pixel {
    xyz: [Float; 3],
    filter_weight_sum: Float,
    splat_xyz: [AtomicFloat; 3],
    pad: Float,
}

#[derive(Debug, Default, Copy, Clone)]
pub struct FilmTilePixel {
    contrib_sum: Spectrum,
    filter_weight_sum: Float,
}

pub struct FilmTile<'a> {
    pub pixel_bounds: Bounds2i,
    filter_radius: Vector2f,
    inv_filter_radius: Vector2f,
    filter_table: &'a [Float; FILTER_TABLE_WIDTH * FILTER_TABLE_WIDTH],
    filter_table_size: usize,
    pixels: Vec<FilmTilePixel>,
    max_sample_luminance: Float,
}

impl<'a> FilmTile<'a> {
    pub fn new(
        pixel_bounds: Bounds2i,
        filter_radius: Vector2f,
        filter_table: &'a [Float; FILTER_TABLE_WIDTH * FILTER_TABLE_WIDTH],
        filter_table_size: usize,
        max_sample_luminance: Float,
    ) -> Self {
        FilmTile {
            pixel_bounds: pixel_bounds,
            filter_radius: filter_radius,
            inv_filter_radius: Vector2f {
                x: 1.0 / filter_radius.x,
                y: 1.0 / filter_radius.y,
            },
            filter_table: filter_table,
            filter_table_size: filter_table_size,
            // TODO: pixels = std::vector<FilmTilePixel>(std::max(0, pixelBounds.Area()));
            pixels: vec![FilmTilePixel::default(); pixel_bounds.area() as usize],
            max_sample_luminance: max_sample_luminance,
        }
    }
    pub fn add_sample(&mut self, p_film: &Point2f, l: &mut Spectrum, sample_weight: Float) {
        // TODO: ProfilePhase _(Prof::AddFilmSample);
        if l.y() > self.max_sample_luminance {
            *l *= Spectrum::new(self.max_sample_luminance / l.y());
        }
        // compute sample's raster bounds
        let p_film_discrete: Point2f = *p_film - Vector2f { x: 0.5, y: 0.5 };
        let p0f: Point2f = pnt2_ceil(&(p_film_discrete - self.filter_radius));
        let mut p0: Point2i = Point2i {
            x: p0f.x as i32,
            y: p0f.y as i32,
        };
        let p1f: Point2f = pnt2_floor(&(p_film_discrete + self.filter_radius));
        let mut p1: Point2i = Point2i {
            x: p1f.x as i32 + 1,
            y: p1f.y as i32 + 1,
        };
        p0 = pnt2_max_pnt2(p0, self.pixel_bounds.p_min);
        p1 = pnt2_min_pnt2(p1, self.pixel_bounds.p_max);

        // loop over filter support and add sample to pixel arrays

        // precompute $x$ and $y$ filter table offsets
        let mut ifx: Vec<usize> = Vec::with_capacity(p1.x as usize - p0.x as usize);
        for x in p0.x..p1.x {
            let fx: Float = ((x as Float - p_film_discrete.x) * self.inv_filter_radius.x
                * self.filter_table_size as Float)
                .abs();
            ifx.push(fx.floor().min(self.filter_table_size as Float - 1.0) as usize);
        }
        let mut ify: Vec<usize> = Vec::with_capacity(p1.y as usize - p0.y as usize);
        for y in p0.y..p1.y {
            let fy: Float = ((y as Float - p_film_discrete.y) * self.inv_filter_radius.y
                * self.filter_table_size as Float)
                .abs();
            ify.push(fy.floor().min(self.filter_table_size as Float - 1.0) as usize);
        }
        for y in p0.y..p1.y {
            for x in p0.x..p1.x {
                // evaluate filter value at $(x,y)$ pixel
                let offset: usize =
                    ify[(y - p0.y) as usize] * self.filter_table_size + ifx[(x - p0.x) as usize];
                let filter_weight: Float = self.filter_table[offset];
                // update pixel values with filtered sample contribution
                let idx = self.get_pixel_index(x, y);
                let ref mut pixel = self.pixels[idx];
                pixel.contrib_sum +=
                    *l * Spectrum::new(sample_weight) * Spectrum::new(filter_weight);
                pixel.filter_weight_sum += filter_weight;
            }
        }
    }
    fn get_pixel_index(&self, x: i32, y: i32) -> usize {
        let width: i32 = self.pixel_bounds.p_max.x - self.pixel_bounds.p_min.x;
        let pidx = (y - self.pixel_bounds.p_min.y) * width + (x - self.pixel_bounds.p_min.x);
        pidx as usize
    }
}

pub struct Film {
    // Film Public Data
    /// The overall resolution of the image in pixels
    pub full_resolution: Point2i,
    /// The length of the diagonal of the film's physical area (specified in mm, stored in meters)
    pub diagonal: Float,
    /// A filter function
    pub filter: Arc<Filter + Sync + Send>,
    /// The filename of the output image
    pub filename: String,
    /// A crop window that may specify a subset of the image to render
    pub cropped_pixel_bounds: Bounds2i,

    // Film Private Data
    pub pixels: RwLock<Vec<Pixel>>,
    filter_table: [Float; FILTER_TABLE_WIDTH * FILTER_TABLE_WIDTH],
    scale: Float,
    max_sample_luminance: Float,
}

impl Film {
    pub fn new(
        resolution: Point2i,
        crop_window: Bounds2f,
        filter: Arc<Filter + Sync + Send>,
        diagonal: Float,
        filename: String,
        scale: Float,
        max_sample_luminance: Float,
    ) -> Self {
        let cropped_pixel_bounds: Bounds2i = Bounds2i {
            p_min: Point2i {
                x: (resolution.x as Float * crop_window.p_min.x).ceil() as i32,
                y: (resolution.y as Float * crop_window.p_min.y).ceil() as i32,
            },
            p_max: Point2i {
                x: (resolution.x as Float * crop_window.p_max.x).ceil() as i32,
                y: (resolution.y as Float * crop_window.p_max.y).ceil() as i32,
            },
        };
        // allocate film image storage
        // let pixels: Vec<Pixel> = vec![Pixel::default(); cropped_pixel_bounds.area() as usize];
        // precompute filter weight table
        let mut filter_table: [Float; FILTER_TABLE_WIDTH * FILTER_TABLE_WIDTH] =
            [0.0; FILTER_TABLE_WIDTH * FILTER_TABLE_WIDTH];
        let mut offset: usize = 0;
        let filter_radius: Vector2f = filter.get_radius();
        for y in 0..FILTER_TABLE_WIDTH {
            for x in 0..FILTER_TABLE_WIDTH {
                let p: Point2f = Point2f {
                    x: (x as Float + 0.5) * filter_radius.x / FILTER_TABLE_WIDTH as Float,
                    y: (y as Float + 0.5) * filter_radius.y / FILTER_TABLE_WIDTH as Float,
                };
                filter_table[offset] = filter.evaluate(p);
                offset += 1;
            }
        }
        Film {
            full_resolution: resolution,
            diagonal: diagonal * 0.001,
            filter: filter,
            filename: filename,
            cropped_pixel_bounds: cropped_pixel_bounds,
            pixels: RwLock::new(vec![Pixel::default(); cropped_pixel_bounds.area() as usize]),
            filter_table: filter_table,
            scale: scale,
            max_sample_luminance: max_sample_luminance,
        }
    }
    pub fn get_sample_bounds(&self) -> Bounds2i {
        let f: Point2f = pnt2_floor(
            &(Point2f {
                x: self.cropped_pixel_bounds.p_min.x as Float,
                y: self.cropped_pixel_bounds.p_min.y as Float,
            } + Vector2f { x: 0.5, y: 0.5 } - self.filter.get_radius()),
        );
        let c: Point2f = pnt2_ceil(
            &(Point2f {
                x: self.cropped_pixel_bounds.p_max.x as Float,
                y: self.cropped_pixel_bounds.p_max.y as Float,
            } - Vector2f { x: 0.5, y: 0.5 } + self.filter.get_radius()),
        );
        let float_bounds: Bounds2f = Bounds2f { p_min: f, p_max: c };
        Bounds2i {
            p_min: Point2i {
                x: float_bounds.p_min.x as i32,
                y: float_bounds.p_min.y as i32,
            },
            p_max: Point2i {
                x: float_bounds.p_max.x as i32,
                y: float_bounds.p_max.y as i32,
            },
        }
    }
    pub fn get_film_tile(&self, sample_bounds: &Bounds2i) -> FilmTile {
        // bound image pixels that samples in _sample_bounds_ contribute to
        let half_pixel: Vector2f = Vector2f { x: 0.5, y: 0.5 };
        let float_bounds: Bounds2f = Bounds2f {
            p_min: Point2f {
                x: sample_bounds.p_min.x as Float,
                y: sample_bounds.p_min.y as Float,
            },
            p_max: Point2f {
                x: sample_bounds.p_max.x as Float,
                y: sample_bounds.p_max.y as Float,
            },
        };
        let p_min: Point2f = float_bounds.p_min - half_pixel - self.filter.get_radius();
        let p0: Point2i = Point2i {
            x: p_min.x.ceil() as i32,
            y: p_min.y.ceil() as i32,
        };
        let p_max: Point2f = float_bounds.p_max - half_pixel + self.filter.get_radius();
        let p1: Point2i = Point2i {
            x: p_max.x.floor() as i32,
            y: p_max.y.floor() as i32,
        } + Point2i { x: 1, y: 1 };
        let tile_pixel_bounds: Bounds2i = bnd2_intersect_bnd2(
            &Bounds2i {
                p_min: p0,
                p_max: p1,
            },
            &self.cropped_pixel_bounds,
        );
        FilmTile::new(
            tile_pixel_bounds,
            self.filter.get_radius(),
            &self.filter_table,
            FILTER_TABLE_WIDTH,
            self.max_sample_luminance,
        )
    }
    pub fn merge_film_tile(&self, tile: &FilmTile) {
        // TODO: ProfilePhase p(Prof::MergeFilmTile);
        // println!("Merging film tile {:?}", tile.pixel_bounds);
        // TODO: std::lock_guard<std::mutex> lock(mutex);
        for pixel in &tile.pixel_bounds {
            // merge _pixel_ into _Film::pixels_
            let idx = tile.get_pixel_index(pixel.x, pixel.y);
            let ref tile_pixel = tile.pixels[idx];
            // START let mut merge_pixel: &mut Pixel = self.get_pixel_mut(pixel);
            assert!(pnt2_inside_exclusive(&pixel, &self.cropped_pixel_bounds));
            let width: i32 = self.cropped_pixel_bounds.p_max.x - self.cropped_pixel_bounds.p_min.x;
            let offset: i32 = (pixel.x - self.cropped_pixel_bounds.p_min.x)
                + (pixel.y - self.cropped_pixel_bounds.p_min.y) * width;
            let mut pixels_write = self.pixels.write().unwrap();
            let mut merge_pixel = pixels_write[offset as usize].clone();
            // END let mut merge_pixel: &mut Pixel = self.get_pixel_mut(pixel);
            let mut xyz: [Float; 3] = [0.0; 3];
            tile_pixel.contrib_sum.to_xyz(&mut xyz);
            for i in 0..3 {
                merge_pixel.xyz[i] += xyz[i];
            }
            merge_pixel.filter_weight_sum += tile_pixel.filter_weight_sum;
            // write pixel back
            pixels_write[offset as usize] = merge_pixel;
        }
    }
    pub fn add_splat(&self, p: &Point2f, v: &Spectrum) {
        let mut v: Spectrum = *v;
        // TODO: ProfilePhase pp(Prof::SplatFilm);
        if v.has_nans() {
            println!(
                "ERROR: Ignoring splatted spectrum with NaN values at ({:?}, {:?})",
                p.x, p.y
            );
            return;
        } else if v.y() < 0.0 as Float {
            println!(
                "ERROR: Ignoring splatted spectrum with negative luminance {:?} at ({:?}, {:?})",
                v.y(),
                p.x,
                p.y
            );
            return;
        } else if v.y().is_infinite() {
            println!(
                "ERROR: Ignoring splatted spectrum with infinite luminance at ({:?}, {:?})",
                p.x, p.y
            );
            return;
        }

        let pi: Point2i = Point2i {
            x: p.x as i32,
            y: p.y as i32,
        };
        if !pnt2_inside_exclusive(&pi, &self.cropped_pixel_bounds) {
            return;
        }
        if v.y() > self.max_sample_luminance {
            v = v * self.max_sample_luminance / v.y();
        }
        let mut xyz: [Float; 3] = [Float::default(); 3];
        println!("v = {:?}", v);
        v.to_xyz(&mut xyz);
        println!("x = {:?}", xyz);
        // let pixel: Pixel = self.get_pixel(&pi);

        let width: i32 = self.cropped_pixel_bounds.p_max.x - self.cropped_pixel_bounds.p_min.x;
        let offset: i32 = (pi.x - self.cropped_pixel_bounds.p_min.x)
            + (pi.y - self.cropped_pixel_bounds.p_min.y) * width;
        let mut pixels_write: RwLockWriteGuard<Vec<Pixel>> = self.pixels.write().unwrap();
        println!("offset: {:?}", offset);
        let pixel_vec: &mut Vec<Pixel> = pixels_write.deref_mut();
        let pixel: &mut Pixel = &mut pixel_vec[offset as usize];

        pixel.splat_xyz[0].add(xyz[0]);
        pixel.splat_xyz[1].add(xyz[1]);
        pixel.splat_xyz[2].add(xyz[2]);
        // pixel_vec[offset as usize] = *pixel;
    }
    #[cfg(not(feature = "openexr"))]
    pub fn write_image(&self, splat_scale: Float) {
        println!("Converting image to RGB and computing final weighted pixel values");
        let mut rgb: Vec<Float> =
            vec![0.0 as Float; (3 * self.cropped_pixel_bounds.area()) as usize];
        let mut offset: usize = 0;
        for p in &self.cropped_pixel_bounds {
            // convert pixel XYZ color to RGB
            let pixel: Pixel = self.get_pixel(&p);
            let start = 3 * offset;
            let mut rgb_array: [Float; 3] = [0.0 as Float; 3];
            xyz_to_rgb(&pixel.xyz, &mut rgb_array); // TODO: Use 'rgb' directly.
            rgb[start + 0] = rgb_array[0];
            rgb[start + 1] = rgb_array[1];
            rgb[start + 2] = rgb_array[2];
            // normalize pixel with weight sum
            let filter_weight_sum: Float = pixel.filter_weight_sum;
            if filter_weight_sum != 0.0 as Float {
                let inv_wt: Float = 1.0 as Float / filter_weight_sum;
                rgb[start + 0] = (rgb[start + 0] * inv_wt).max(0.0 as Float);
                rgb[start + 1] = (rgb[start + 1] * inv_wt).max(0.0 as Float);
                rgb[start + 2] = (rgb[start + 2] * inv_wt).max(0.0 as Float);
            }
            // add splat value at pixel
            let mut splat_rgb: [Float; 3] = [0.0 as Float; 3];
            let splat_xyz: [Float; 3] = [
                Float::from(pixel.splat_xyz[0].clone()),
                Float::from(pixel.splat_xyz[1].clone()),
                Float::from(pixel.splat_xyz[2].clone()),
            ];
            println!("splat_xyz = {:?}", splat_xyz);
            xyz_to_rgb(&splat_xyz, &mut splat_rgb);
            rgb[start + 0] += splat_scale * splat_rgb[0];
            rgb[start + 1] += splat_scale * splat_rgb[1];
            rgb[start + 2] += splat_scale * splat_rgb[2];
            // scale pixel value by _scale_
            rgb[start + 0] *= self.scale;
            rgb[start + 1] *= self.scale;
            rgb[start + 2] *= self.scale;
            offset += 1;
        }
        let filename = "pbrt.png";
        println!(
            "Writing image {:?} with bounds {:?}",
            filename, // TODO: self.filename,
            self.cropped_pixel_bounds
        );
        // TODO: pbrt::WriteImage(filename, &rgb[0], croppedPixelBounds, fullResolution);
        let mut buffer: Vec<u8> = vec![0.0 as u8; (3 * self.cropped_pixel_bounds.area()) as usize];
        // 8-bit format; apply gamma (see WriteImage(...) in imageio.cpp)
        let width: u32 =
            (self.cropped_pixel_bounds.p_max.x - self.cropped_pixel_bounds.p_min.x) as u32;
        let height: u32 =
            (self.cropped_pixel_bounds.p_max.y - self.cropped_pixel_bounds.p_min.y) as u32;
        for y in 0..height {
            for x in 0..width {
                // red
                let index: usize = (3 * (y * width + x) + 0) as usize;
                buffer[index] = clamp_t(
                    255.0 as Float * gamma_correct(rgb[index]) + 0.5,
                    0.0 as Float,
                    255.0 as Float,
                ) as u8;
                // green
                let index: usize = (3 * (y * width + x) + 1) as usize;
                buffer[index] = clamp_t(
                    255.0 as Float * gamma_correct(rgb[index]) + 0.5,
                    0.0 as Float,
                    255.0 as Float,
                ) as u8;
                // blue
                let index: usize = (3 * (y * width + x) + 2) as usize;
                buffer[index] = clamp_t(
                    255.0 as Float * gamma_correct(rgb[index]) + 0.5,
                    0.0 as Float,
                    255.0 as Float,
                ) as u8;
            }
        }
        // write "pbrt.png" to disk
        image::save_buffer(
            &Path::new("pbrt.png"),
            &buffer,
            width,
            height,
            image::RGB(8),
        ).unwrap();
    }
    #[cfg(feature = "openexr")]
    pub fn write_image(&self, splat_scale: Float) {
        println!("Converting image to RGB and computing final weighted pixel values");
        let mut rgb: Vec<Float> =
            vec![0.0 as Float; (3 * self.cropped_pixel_bounds.area()) as usize];
        let mut exr: Vec<(Float, Float, Float)> = // copy data for OpenEXR image
            vec![(0.0_f32, 0.0_f32, 0.0_f32); self.cropped_pixel_bounds.area() as usize];
        let mut offset: usize = 0;
        for p in &self.cropped_pixel_bounds {
            // convert pixel XYZ color to RGB
            let pixel: Pixel = self.get_pixel(&p);
            let start = 3 * offset;
            let mut rgb_array: [Float; 3] = [0.0 as Float; 3];
            xyz_to_rgb(&pixel.xyz, &mut rgb_array); // TODO: Use 'rgb' directly.
            rgb[start + 0] = rgb_array[0];
            rgb[start + 1] = rgb_array[1];
            rgb[start + 2] = rgb_array[2];
            // normalize pixel with weight sum
            let filter_weight_sum: Float = pixel.filter_weight_sum;
            if filter_weight_sum != 0.0 as Float {
                let inv_wt: Float = 1.0 as Float / filter_weight_sum;
                rgb[start + 0] = (rgb[start + 0] * inv_wt).max(0.0 as Float);
                rgb[start + 1] = (rgb[start + 1] * inv_wt).max(0.0 as Float);
                rgb[start + 2] = (rgb[start + 2] * inv_wt).max(0.0 as Float);
            }
            // add splat value at pixel
            let mut splat_rgb: [Float; 3] = [0.0 as Float; 3];
            let splat_xyz: [Float; 3] = [
                Float::from(pixel.splat_xyz[0].clone()),
                Float::from(pixel.splat_xyz[1].clone()),
                Float::from(pixel.splat_xyz[2].clone()),
            ];
            xyz_to_rgb(&splat_xyz, &mut splat_rgb);
            rgb[start + 0] += splat_scale * splat_rgb[0];
            rgb[start + 1] += splat_scale * splat_rgb[1];
            rgb[start + 2] += splat_scale * splat_rgb[2];
            // scale pixel value by _scale_
            rgb[start + 0] *= self.scale;
            rgb[start + 1] *= self.scale;
            rgb[start + 2] *= self.scale;
            // copy data for OpenEXR image
            exr[offset].0 = rgb[start + 0];
            exr[offset].1 = rgb[start + 1];
            exr[offset].2 = rgb[start + 2];
            offset += 1;
        }
        let filename = "pbrt.png";
        println!(
            "Writing image {:?} with bounds {:?}",
            filename, // TODO: self.filename,
            self.cropped_pixel_bounds
        );
        // TODO: pbrt::WriteImage(filename, &rgb[0], croppedPixelBounds, fullResolution);
        let mut buffer: Vec<u8> = vec![0.0 as u8; (3 * self.cropped_pixel_bounds.area()) as usize];
        // 8-bit format; apply gamma (see WriteImage(...) in imageio.cpp)
        let width: u32 =
            (self.cropped_pixel_bounds.p_max.x - self.cropped_pixel_bounds.p_min.x) as u32;
        let height: u32 =
            (self.cropped_pixel_bounds.p_max.y - self.cropped_pixel_bounds.p_min.y) as u32;
        // OpenEXR
        let filename = "pbrt_rust.exr";
        println!(
            "Writing image {:?} with bounds {:?}",
            filename, // TODO: self.filename,
            self.cropped_pixel_bounds
        );
        let mut file = std::fs::File::create("pbrt_rust.exr").unwrap();
        let mut output_file = ScanlineOutputFile::new(
            &mut file,
            Header::new()
                .set_resolution(width, height)
                .add_channel("R", PixelType::FLOAT)
                .add_channel("G", PixelType::FLOAT)
                .add_channel("B", PixelType::FLOAT),
        ).unwrap();
        let mut fb = FrameBuffer::new(width as u32, height as u32);
        fb.insert_channels(&["R", "G", "B"], &exr);
        output_file.write_pixels(&fb).unwrap();

        // OpenEXR
        for y in 0..height {
            for x in 0..width {
                // red
                let index: usize = (3 * (y * width + x) + 0) as usize;
                buffer[index] = clamp_t(
                    255.0 as Float * gamma_correct(rgb[index]) + 0.5,
                    0.0 as Float,
                    255.0 as Float,
                ) as u8;
                // green
                let index: usize = (3 * (y * width + x) + 1) as usize;
                buffer[index] = clamp_t(
                    255.0 as Float * gamma_correct(rgb[index]) + 0.5,
                    0.0 as Float,
                    255.0 as Float,
                ) as u8;
                // blue
                let index: usize = (3 * (y * width + x) + 2) as usize;
                buffer[index] = clamp_t(
                    255.0 as Float * gamma_correct(rgb[index]) + 0.5,
                    0.0 as Float,
                    255.0 as Float,
                ) as u8;
            }
        }
        // write "pbrt.png" to disk
        image::save_buffer(
            &Path::new("pbrt.png"),
            &buffer,
            width,
            height,
            image::RGB(8),
        ).unwrap();
    }
    pub fn get_pixel(&self, p: &Point2i) -> Pixel {
        assert!(pnt2_inside_exclusive(p, &self.cropped_pixel_bounds));
        let width: i32 = self.cropped_pixel_bounds.p_max.x - self.cropped_pixel_bounds.p_min.x;
        let offset: i32 = (p.x - self.cropped_pixel_bounds.p_min.x)
            + (p.y - self.cropped_pixel_bounds.p_min.y) * width;
        self.pixels.read().unwrap()[offset as usize].clone()
    }
}
