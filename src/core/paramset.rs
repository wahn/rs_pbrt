//! Bundle up parameters and their values in a generic way.

// std
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Arc;
// pbrt
use crate::core::floatfile::read_float_file;
use crate::core::geometry::{Normal3f, Point2f, Point3f, Vector2f, Vector3f};
use crate::core::pbrt::{Float, Spectrum};
use crate::core::spectrum::blackbody_normalized;
use crate::core::spectrum::{CIE_LAMBDA, N_CIE_SAMPLES};
use crate::core::texture::Texture;
use crate::textures::constant::ConstantTexture;

// see paramset.h

pub struct ParamSetItem<T> {
    pub name: String,
    pub values: Vec<T>,
    pub n_values: usize,
    pub looked_up: bool, // false
}

#[derive(Default)]
pub struct ParamSet {
    pub key_word: String,
    pub name: String,
    pub tex_type: String,
    pub tex_name: String,
    pub bools: Vec<ParamSetItem<bool>>,
    pub ints: Vec<ParamSetItem<i32>>,
    pub floats: Vec<ParamSetItem<Float>>,
    pub point2fs: Vec<ParamSetItem<Point2f>>,
    pub vector2fs: Vec<ParamSetItem<Vector2f>>,
    pub point3fs: Vec<ParamSetItem<Point3f>>,
    pub vector3fs: Vec<ParamSetItem<Vector3f>>,
    pub normals: Vec<ParamSetItem<Normal3f>>,
    pub spectra: Vec<ParamSetItem<Spectrum>>,
    pub strings: Vec<ParamSetItem<String>>,
    pub textures: Vec<ParamSetItem<String>>,
}

impl ParamSet {
    pub fn reset(&mut self, key_word: String, name: String, tex_type: String, tex_name: String) {
        self.key_word = key_word;
        self.name = name;
        self.tex_type = tex_type;
        self.tex_name = tex_name;
        self.bools.clear();
        self.ints.clear();
        self.floats.clear();
        self.point2fs.clear();
        self.vector2fs.clear();
        self.point3fs.clear();
        self.vector3fs.clear();
        self.normals.clear();
        self.spectra.clear();
        self.strings.clear();
        self.textures.clear();
    }
    pub fn add_float(&mut self, name: String, value: Float) {
        self.floats.push(ParamSetItem::<Float> {
            name,
            values: vec![value],
            n_values: 1_usize,
            looked_up: false,
        });
    }
    pub fn add_floats(&mut self, name: String, values: Vec<Float>) {
        let n_values: usize = values.len();
        self.floats.push(ParamSetItem::<Float> {
            name,
            values,
            n_values,
            looked_up: false,
        });
    }
    pub fn add_int(&mut self, name: String, value: i32) {
        self.ints.push(ParamSetItem::<i32> {
            name,
            values: vec![value],
            n_values: 1_usize,
            looked_up: false,
        });
    }
    pub fn add_ints(&mut self, name: String, values: Vec<i32>) {
        let n_values: usize = values.len();
        self.ints.push(ParamSetItem::<i32> {
            name,
            values,
            n_values,
            looked_up: false,
        });
    }
    pub fn add_bool(&mut self, name: String, value: bool) {
        self.bools.push(ParamSetItem::<bool> {
            name,
            values: vec![value],
            n_values: 1_usize,
            looked_up: false,
        });
    }
    pub fn add_point2f(&mut self, name: String, value: Point2f) {
        self.point2fs.push(ParamSetItem::<Point2f> {
            name,
            values: vec![value],
            n_values: 1_usize,
            looked_up: false,
        });
    }
    pub fn add_point2fs(&mut self, name: String, values: Vec<Float>) {
        let n_values: usize = values.len();
        let mut p_values: Vec<Point2f> = Vec::new();
        let n_points: usize = values.len() / 2_usize;
        assert!(
            n_values % 2 == 0,
            "point parameters need 2 coordinates ({} found for {:?})",
            n_values,
            name
        );
        for i in 0..n_points {
            let x: Float = values[i * 2 + 0];
            let y: Float = values[i * 2 + 1];
            p_values.push(Point2f { x, y });
        }
        self.point2fs.push(ParamSetItem::<Point2f> {
            name,
            values: p_values,
            n_values: n_points,
            looked_up: false,
        });
    }
    pub fn add_point3f(&mut self, name: String, value: Point3f) {
        self.point3fs.push(ParamSetItem::<Point3f> {
            name,
            values: vec![value],
            n_values: 1_usize,
            looked_up: false,
        });
    }
    pub fn add_point3fs(&mut self, name: String, values: Vec<Float>) {
        let n_values: usize = values.len();
        let mut p_values: Vec<Point3f> = Vec::new();
        let n_points: usize = values.len() / 3_usize;
        assert!(
            n_values % 3 == 0,
            "point parameters need 3 coordinates ({} found for {:?})",
            n_values,
            name
        );
        for i in 0..n_points {
            let x: Float = values[i * 3 + 0];
            let y: Float = values[i * 3 + 1];
            let z: Float = values[i * 3 + 2];
            p_values.push(Point3f { x, y, z });
        }
        self.point3fs.push(ParamSetItem::<Point3f> {
            name,
            values: p_values,
            n_values: n_points,
            looked_up: false,
        });
    }
    pub fn add_sampled_spectrum_files(&mut self, name: String, names: Vec<String>) {
        // TODO: cachedSpectra
        self.erase_spectrum(name.clone());
        let mut s: Vec<Spectrum> = Vec::with_capacity(names.len());
        for i in 0..names.len() {
            // std::string filename = AbsolutePath(ResolveFilename(names[i]));
            let fn_str: &String = &names[i];
            let _f = File::open(fn_str.clone()).unwrap();
            let ip: &Path = Path::new(fn_str.as_str());
            if ip.is_relative() {
                let cp: PathBuf = env::current_dir().unwrap();
                let pb: PathBuf = cp.join(ip);
                let search_directory: &Path = pb.as_path().parent().unwrap();
                let mut path_buf: PathBuf = PathBuf::from("/");
                path_buf.push(search_directory);
                path_buf.push(ip.file_name().unwrap());
                let filename = String::from(path_buf.to_str().unwrap());
                let mut vals: Vec<Float> = Vec::new();
                if !read_float_file(&filename, &mut vals) {
                    println!(
                        "WARNING: Unable to read SPD file {:?}. Using black distribution.",
                        filename
                    );
                    s.push(Spectrum::default());
                } else {
                    if vals.len() % 2 == 1_usize {
                        println!(
                            "WARNING: Extra value found in spectrum file {:?}. Ignoring it.",
                            filename
                        );
                    }
                    let mut wls: Vec<Float> = Vec::new();
                    let mut v: Vec<Float> = Vec::new();
                    for j in 0..(vals.len() / 2_usize) {
                        wls.push(vals[2 * j]);
                        v.push(vals[2 * j + 1]);
                    }
                    s.push(Spectrum::from_sampled(&wls[..], &v[..], wls.len() as i32));
                }
            }
        }
        let n_values: usize = s.len();
        self.spectra.push(ParamSetItem::<Spectrum> {
            name: name.clone(),
            values: s,
            n_values,
            looked_up: false,
        });
    }
    pub fn add_string(&mut self, name: String, value: String) {
        self.strings.push(ParamSetItem::<String> {
            name,
            values: vec![value],
            n_values: 1_usize,
            looked_up: false,
        });
    }
    pub fn add_texture(&mut self, name: String, value: String) {
        self.textures.push(ParamSetItem::<String> {
            name,
            values: vec![value],
            n_values: 1_usize,
            looked_up: false,
        });
    }
    pub fn add_vector3f(&mut self, name: String, value: Vector3f) {
        self.vector3fs.push(ParamSetItem::<Vector3f> {
            name,
            values: vec![value],
            n_values: 1_usize,
            looked_up: false,
        });
    }
    pub fn add_vector3fs(&mut self, name: String, values: Vec<Float>) {
        let n_values: usize = values.len();
        let mut p_values: Vec<Vector3f> = Vec::new();
        let n_vectors: usize = values.len() / 3_usize;
        assert!(n_values % 3 == 0, "vector parameters need 3 coordinates");
        for i in 0..n_vectors {
            let x: Float = values[i * 3 + 0];
            let y: Float = values[i * 3 + 1];
            let z: Float = values[i * 3 + 2];
            p_values.push(Vector3f { x, y, z });
        }
        self.vector3fs.push(ParamSetItem::<Vector3f> {
            name,
            values: p_values,
            n_values: n_vectors,
            looked_up: false,
        });
    }
    pub fn add_normal3f(&mut self, name: String, value: Normal3f) {
        self.normals.push(ParamSetItem::<Normal3f> {
            name,
            values: vec![value],
            n_values: 1_usize,
            looked_up: false,
        });
    }
    pub fn add_normal3fs(&mut self, name: String, values: Vec<Float>) {
        let n_values: usize = values.len();
        let mut p_values: Vec<Normal3f> = Vec::new();
        let n_normals: usize = values.len() / 3_usize;
        assert!(n_values % 3 == 0, "normal parameters need 3 coordinates");
        for i in 0..n_normals {
            let x: Float = values[i * 3 + 0];
            let y: Float = values[i * 3 + 1];
            let z: Float = values[i * 3 + 2];
            p_values.push(Normal3f { x, y, z });
        }
        self.normals.push(ParamSetItem::<Normal3f> {
            name,
            values: p_values,
            n_values: n_normals,
            looked_up: false,
        });
    }
    pub fn add_rgb_spectrum(&mut self, name: String, value: Spectrum) {
        self.spectra.push(ParamSetItem::<Spectrum> {
            name,
            values: vec![value],
            n_values: 1_usize,
            looked_up: false,
        });
    }
    pub fn add_blackbody_spectrum(&mut self, name: String, values: Vec<Float>) {
        assert!(values.len() % 2 == 0);
        // temperature (K), scale, ...
        let n_values: usize = values.len() / 2_usize;
        let mut s: Vec<Spectrum> = Vec::with_capacity(n_values);
        let mut v: Vec<Float> = Vec::with_capacity(N_CIE_SAMPLES as usize);
        for i in 0..n_values {
            blackbody_normalized(&CIE_LAMBDA, N_CIE_SAMPLES as usize, values[2 * i], &mut v);
            s.push(
                Spectrum::from_sampled(&CIE_LAMBDA, &v, N_CIE_SAMPLES as i32) * values[2 * i + 1],
            );
        }
        self.spectra.push(ParamSetItem::<Spectrum> {
            name,
            values: s,
            n_values,
            looked_up: false,
        });
    }
    pub fn copy_from(&mut self, param_set: &ParamSet) {
        self.key_word = param_set.key_word.clone();
        // self.name = param_set.name.clone();
        self.bools.clear();
        for b in &param_set.bools {
            let mut values: Vec<bool> = Vec::new();
            for ix in 0..b.n_values {
                values.push(b.values[ix]);
            }
            self.bools.push(ParamSetItem::<bool> {
                name: b.name.clone(),
                values,
                n_values: b.n_values,
                looked_up: false,
            });
        }
        self.ints.clear();
        for i in &param_set.ints {
            let mut values: Vec<i32> = Vec::new();
            for ix in 0..i.n_values {
                values.push(i.values[ix]);
            }
            self.ints.push(ParamSetItem::<i32> {
                name: i.name.clone(),
                values,
                n_values: i.n_values,
                looked_up: false,
            });
        }
        self.floats.clear();
        for f in &param_set.floats {
            let mut values: Vec<Float> = Vec::new();
            for ix in 0..f.n_values {
                values.push(f.values[ix]);
            }
            self.floats.push(ParamSetItem::<Float> {
                name: f.name.clone(),
                values,
                n_values: f.n_values,
                looked_up: false,
            });
        }
        self.point2fs.clear();
        self.vector2fs.clear();
        self.point3fs.clear();
        for p in &param_set.point3fs {
            let mut values: Vec<Point3f> = Vec::new();
            for ix in 0..p.n_values {
                values.push(p.values[ix].clone());
            }
            self.point3fs.push(ParamSetItem::<Point3f> {
                name: p.name.clone(),
                values,
                n_values: p.n_values,
                looked_up: false,
            });
        }
        self.vector3fs.clear();
        self.normals.clear();
        self.spectra.clear();
        for s in &param_set.spectra {
            let mut values: Vec<Spectrum> = Vec::new();
            for ix in 0..s.n_values {
                values.push(s.values[ix].clone());
            }
            self.spectra.push(ParamSetItem::<Spectrum> {
                name: s.name.clone(),
                values,
                n_values: s.n_values,
                looked_up: false,
            });
        }
        self.strings.clear();
        for s in &param_set.strings {
            let mut values: Vec<String> = Vec::new();
            for ix in 0..s.n_values {
                values.push(s.values[ix].clone());
            }
            self.strings.push(ParamSetItem::<String> {
                name: s.name.clone(),
                values,
                n_values: s.n_values,
                looked_up: false,
            });
        }
        self.textures.clear();
        for s in &param_set.textures {
            let mut values: Vec<String> = Vec::new();
            for ix in 0..s.n_values {
                values.push(s.values[ix].clone());
            }
            self.textures.push(ParamSetItem::<String> {
                name: s.name.clone(),
                values,
                n_values: s.n_values,
                looked_up: false,
            });
        }
    }
    pub fn erase_spectrum(&mut self, name: String) -> bool {
        for i in 0..self.spectra.len() {
            if self.spectra[i].name == name {
                self.spectra.remove(i);
                return true;
            }
        }
        false
    }
    pub fn find_one_float(&self, name: &str, d: Float) -> Float {
        for v in &self.floats {
            if v.name == name && v.n_values == 1 {
                // v.looked_up = true;
                return v.values[0];
            }
        }
        d
    }
    pub fn find_one_int(&self, name: &str, d: i32) -> i32 {
        for v in &self.ints {
            if v.name == name && v.n_values == 1 {
                // v.looked_up = true;
                return v.values[0];
            }
        }
        d
    }
    pub fn find_one_bool(&self, name: &str, d: bool) -> bool {
        for v in &self.bools {
            if v.name == name && v.n_values == 1 {
                // v.looked_up = true;
                return v.values[0];
            }
        }
        d
    }
    pub fn find_one_point3f(&self, name: &str, d: Point3f) -> Point3f {
        for v in &self.point3fs {
            if v.name == name && v.n_values == 1 {
                // v.looked_up = true;
                return v.values[0];
            }
        }
        d
    }
    pub fn find_one_vector3f(&self, name: &str, d: Vector3f) -> Vector3f {
        for v in &self.vector3fs {
            if v.name == name && v.n_values == 1 {
                // v.looked_up = true;
                return v.values[0];
            }
        }
        d
    }
    pub fn find_one_spectrum(&self, name: &str, d: Spectrum) -> Spectrum {
        for v in &self.spectra {
            if v.name == name && v.n_values == 1 {
                // v.looked_up = true;
                return v.values[0];
            }
        }
        d
    }
    pub fn find_one_string(&self, name: &str, d: String) -> String {
        for v in &self.strings {
            if v.name == name && v.n_values == 1 {
                // v.looked_up = true;
                return v.values[0].clone();
            }
        }
        d
    }
    pub fn find_one_filename(&self, name: &str, d: String) -> String {
        let filename: String = self.find_one_string(name, String::new());
        if filename == "" {
            return d;
        }
        // TODO: filename = AbsolutePath(ResolveFilename(filename));
        filename
    }
    pub fn find_texture(&self, name: &str) -> String {
        let d: String = String::new();
        lookup_one(&self.textures, name, d)
    }
    pub fn find_int(&self, name: &str) -> Vec<i32> {
        let mut values: Vec<i32> = Vec::new();
        for v in &self.ints {
            if v.name == name {
                let n_values = v.n_values;
                // v.looked_up = true;
                for i in 0..n_values {
                    values.push(v.values[i]);
                }
            }
        }
        values
    }
    pub fn find_float(&self, name: &str) -> Vec<Float> {
        let mut values: Vec<Float> = Vec::new();
        for v in &self.floats {
            if v.name == name {
                let n_values = v.n_values;
                // v.looked_up = true;
                for i in 0..n_values {
                    values.push(v.values[i]);
                }
            }
        }
        values
    }
    pub fn find_point2f(&self, name: &str) -> Vec<Point2f> {
        let mut values: Vec<Point2f> = Vec::new();
        for v in &self.point2fs {
            if v.name == name {
                let n_values = v.n_values;
                // v.looked_up = true;
                for i in 0..n_values {
                    values.push(v.values[i]);
                }
            }
        }
        values
    }
    pub fn find_vector2f(&self, name: &str) -> Vec<Vector2f> {
        let mut values: Vec<Vector2f> = Vec::new();
        for v in &self.vector2fs {
            if v.name == name {
                let n_values = v.n_values;
                // v.looked_up = true;
                for i in 0..n_values {
                    values.push(v.values[i]);
                }
            }
        }
        values
    }
    pub fn find_point3f(&self, name: &str) -> Vec<Point3f> {
        let mut values: Vec<Point3f> = Vec::new();
        for v in &self.point3fs {
            if v.name == name {
                let n_values = v.n_values;
                // v.looked_up = true;
                for i in 0..n_values {
                    values.push(v.values[i]);
                }
            }
        }
        values
    }
    pub fn find_vector3f(&self, name: &str) -> Vec<Vector3f> {
        let mut values: Vec<Vector3f> = Vec::new();
        for v in &self.vector3fs {
            if v.name == name {
                let n_values = v.n_values;
                // v.looked_up = true;
                for i in 0..n_values {
                    values.push(v.values[i]);
                }
            }
        }
        values
    }
    pub fn find_normal3f(&self, name: &str) -> Vec<Normal3f> {
        let mut values: Vec<Normal3f> = Vec::new();
        for v in &self.normals {
            if v.name == name {
                let n_values = v.n_values;
                // v.looked_up = true;
                for i in 0..n_values {
                    values.push(v.values[i]);
                }
            }
        }
        values
    }
    pub fn find_spectrum(&self, name: &str) -> Vec<Spectrum> {
        let mut values: Vec<Spectrum> = Vec::new();
        for v in &self.spectra {
            if v.name == name {
                let n_values = v.n_values;
                // v.looked_up = true;
                for i in 0..n_values {
                    values.push(v.values[i]);
                }
            }
        }
        values
    }
}

#[derive(Default)]
pub struct TextureParams {
    pub float_textures: Arc<HashMap<String, Arc<dyn Texture<Float> + Send + Sync>>>,
    pub spectrum_textures: Arc<HashMap<String, Arc<dyn Texture<Spectrum> + Send + Sync>>>,
    pub geom_params: ParamSet,
    pub material_params: ParamSet,
}

impl TextureParams {
    pub fn new(
        geom_params: ParamSet,
        material_params: ParamSet,
        f_tex: Arc<HashMap<String, Arc<dyn Texture<Float> + Send + Sync>>>,
        s_tex: Arc<HashMap<String, Arc<dyn Texture<Spectrum> + Send + Sync>>>,
    ) -> Self {
        TextureParams {
            float_textures: f_tex,
            spectrum_textures: s_tex,
            geom_params,
            material_params,
        }
    }
    pub fn get_spectrum_texture(
        &mut self,
        n: &str,
        def: Spectrum,
    ) -> Arc<dyn Texture<Spectrum> + Send + Sync> {
        let mut name: String = self.geom_params.find_texture(n);
        if name == "" {
            name = self.material_params.find_texture(n);
        }
        if name != "" {
            match self.spectrum_textures.get(name.as_str()) {
                Some(spectrum_texture) => {
                    return spectrum_texture.clone();
                }
                None => {
                    panic!(
                        "Couldn't find spectrum texture named \"{}\" for parameter \"{}\"",
                        name, n
                    );
                }
            }
        }
        let mut val: Spectrum = self.material_params.find_one_spectrum(n.clone(), def);
        val = self.geom_params.find_one_spectrum(n.clone(), val);
        Arc::new(ConstantTexture { value: val })
    }
    pub fn get_spectrum_texture_or_null(
        &mut self,
        n: &str,
    ) -> Option<Arc<dyn Texture<Spectrum> + Send + Sync>> {
        let mut name: String = self.geom_params.find_texture(n);
        if name == "" {
            name = self.material_params.find_texture(n);
        }
        if name != String::new() {
            match self.spectrum_textures.get(name.as_str()) {
                Some(spectrum_texture) => return Some(spectrum_texture.clone()),
                None => {
                    println!(
                        "Couldn't find spectrum texture named \"{}\" for parameter \"{}\"",
                        name, n
                    );
                    return None;
                }
            }
        }
        let mut val: Vec<Spectrum> = self.material_params.find_spectrum(n);
        if val.len() == 0_usize {
            val = self.geom_params.find_spectrum(n);
        }
        if val.len() == 0_usize {
            None
        } else {
            Some(Arc::new(ConstantTexture { value: val[0] }))
        }
    }
    pub fn get_float_texture(&mut self, n: &str, def: Float) -> Arc<dyn Texture<Float> + Send + Sync> {
        let tex_option = self.get_float_texture_or_null(n);
        if let Some(tex) = tex_option {
            tex
        } else {
            let mut val: Float = self.material_params.find_one_float(n, def);
            val = self.geom_params.find_one_float(n, val);
            Arc::new(ConstantTexture { value: val })
        }
    }
    pub fn get_float_texture_or_null(
        &mut self,
        n: &str,
    ) -> Option<Arc<dyn Texture<Float> + Send + Sync>> {
        let mut name: String = self.geom_params.find_texture(n);
        if name == "" {
            let s: Vec<Float> = self.geom_params.find_float(n);
            if s.len() > 1 {
                println!(
                    "Ignoring excess values provided with parameter \"{}\"",
                    n.clone()
                );
            } else if s.len() != 0 {
                return Some(Arc::new(ConstantTexture { value: s[0] }));
            }
            name = self.material_params.find_texture(n);
        }
        if name != String::new() {
            match self.float_textures.get(name.as_str()) {
                Some(float_texture) => {
                    return Some(float_texture.clone());
                }
                None => {
                    println!(
                        "Couldn't find float texture named \"{}\" for parameter \"{}\"",
                        name, n
                    );
                    return None;
                }
            }
        }
        let mut val: Vec<Float> = self.material_params.find_float(n);
        if val.len() == 0_usize {
            val = self.geom_params.find_float(n);
        }
        if val.len() == 0_usize {
            None
        } else {
            Some(Arc::new(ConstantTexture { value: val[0] }))
        }
    }
    pub fn find_float(&mut self, name: &str, d: Float) -> Float {
        self.geom_params
            .find_one_float(name, self.material_params.find_one_float(name, d))
    }
    pub fn find_string(&mut self, name: &str, d: String) -> String {
        self.geom_params
            .find_one_string(name, self.material_params.find_one_string(name, d))
    }
    pub fn find_filename(&mut self, name: &str, d: String) -> String {
        self.geom_params
            .find_one_filename(name, self.material_params.find_one_filename(name, d))
    }
    pub fn find_int(&mut self, name: &str, d: i32) -> i32 {
        self.geom_params
            .find_one_int(name, self.material_params.find_one_int(name, d))
    }
    pub fn find_bool(&mut self, name: &str, d: bool) -> bool {
        self.geom_params
            .find_one_bool(name, self.material_params.find_one_bool(name, d))
    }
    pub fn find_vector3f(&mut self, name: &str, d: Vector3f) -> Vector3f {
        self.geom_params
            .find_one_vector3f(name, self.material_params.find_one_vector3f(name, d))
    }
    pub fn find_spectrum(&mut self, name: &str, d: Spectrum) -> Spectrum {
        self.geom_params
            .find_one_spectrum(name, self.material_params.find_one_spectrum(name, d))
    }
}

/// Replaces a macro on the C++ side.
pub fn lookup_one<T>(vec: &Vec<ParamSetItem<T>>, name: &str, d: T) -> T
where
    T: Clone,
{
    for v in vec {
        if v.name == name && v.n_values == 1_usize {
            // v.looked_up = true;
            return v.values[0].clone();
        }
    }
    d
}
