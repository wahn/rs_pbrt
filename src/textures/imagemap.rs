extern crate image;

// std
use std::path::Path;
use std::sync::Arc;
// others
use image::{DynamicImage, ImageResult};
// pbrt
use core::geometry::{Point2f, Point2i, Vector2f};
use core::interaction::SurfaceInteraction;
use core::mipmap::{ImageWrap, MipMap};
use core::pbrt::{Float, Spectrum};
use core::texture::{Texture, TextureMapping2D};

// see imagemap.h

pub struct ImageTexture {
    pub mapping: Box<TextureMapping2D + Send + Sync>,
    pub mipmap: Arc<MipMap>,
}

impl ImageTexture {
    pub fn new(
        mapping: Box<TextureMapping2D + Send + Sync>,
        filename: String,
        do_trilinear: bool,
        max_aniso: Float,
        wrap_mode: ImageWrap,
        _scale: Float,
        _gamma: bool,
    ) -> Self {
        let path = Path::new(&filename);
        let img_result: ImageResult<DynamicImage> = image::open(path);
        if !img_result.is_ok() {
            panic!("Error reading \"{}\"", filename);
        }
        let buf = img_result.unwrap();
        let rgb = buf.to_rgb();
        let res = Point2i {
            x: rgb.width() as i32,
            y: rgb.height() as i32,
        };
        // instead of convertIn(texels[i], &convertedTexels[i], scale, gamma);
        let mut texels: Vec<Spectrum> =
            rgb.pixels().map(|p| Spectrum::from_srgb(&p.data)).collect();
        // flip image in y; texture coordinate space has (0,0) at the
        // lower left corner.
        for y in 0..res.y / 2 {
            for x in 0..res.x {
                let o1 = (y * res.x + x) as usize;
                let o2 = ((res.y - 1 - y) * res.x + x) as usize;
                texels.swap(o1, o2);
            }
        }
        // create _MipMap_ from converted texels (see above)
        let mipmap = Arc::new(MipMap::new(
            &res,
            &texels[..],
            do_trilinear,
            max_aniso,
            wrap_mode,
        ));
        ImageTexture {
            mapping: mapping,
            mipmap: mipmap,
        }
    }
    pub fn convert_out(from: &Spectrum, to: &mut Spectrum) {
        let mut rgb: [Float; 3] = [0.0 as Float; 3];
        from.to_rgb(&mut rgb);
        *to = Spectrum::from_rgb(&rgb);
    }
}

impl Texture<Float> for ImageTexture {
    fn evaluate(&self, si: &SurfaceInteraction) -> Float {
        // Vector2f dstdx, dstdy;
        // Point2f st = mapping->Map(si, &dstdx, &dstdy);
        // Tmemory mem = mipmap->Lookup(st, dstdx, dstdy);
        // Treturn ret;
        // convertOut(mem, &ret);
        // return ret;
        // let mut dstdx: Vector2f = Vector2f::default();
        // let mut dstdy: Vector2f = Vector2f::default();
        // let st: Point2f = self.mapping.map(si, &mut dstdx, &mut dstdy);
        // let mem: Spectrum = self.mipmap.lookup_pnt_vec_vec(&st, &mut dstdx, &mut dstdy);
        let mut ret: Float = 0.0 as Float;
        // ImageTexture::convert_out(&mem, &mut ret);
        ret
    }
}

impl Texture<Spectrum> for ImageTexture {
    fn evaluate(&self, si: &SurfaceInteraction) -> Spectrum {
        // Vector2f dstdx, dstdy;
        // Point2f st = mapping->Map(si, &dstdx, &dstdy);
        // Tmemory mem = mipmap->Lookup(st, dstdx, dstdy);
        // Treturn ret;
        // convertOut(mem, &ret);
        // return ret;
        let mut dstdx: Vector2f = Vector2f::default();
        let mut dstdy: Vector2f = Vector2f::default();
        let st: Point2f = self.mapping.map(si, &mut dstdx, &mut dstdy);
        let mem: Spectrum = self.mipmap.lookup_pnt_vec_vec(&st, &mut dstdx, &mut dstdy);
        let mut ret: Spectrum = Spectrum::new(0.0);
        ImageTexture::convert_out(&mem, &mut ret);
        ret
    }
}
