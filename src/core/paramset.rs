// std
use std::collections::HashMap;
use std::sync::Arc;
// pbrt
use core::geometry::{Normal3f, Point2f, Point3f, Vector2f, Vector3f};
use core::pbrt::{Float, Spectrum};
use core::spectrum::blackbody_normalized;
use core::spectrum::{CIE_LAMBDA, N_CIE_SAMPLES};
use core::texture::Texture;
use textures::constant::ConstantTexture;

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
            name: name,
            values: vec![value],
            n_values: 1_usize,
            looked_up: false,
        });
    }
    pub fn add_floats(&mut self, name: String, values: Vec<Float>) {
        let n_values: usize = values.len();
        self.floats.push(ParamSetItem::<Float> {
            name: name,
            values: values,
            n_values: n_values,
            looked_up: false,
        });
    }
    pub fn add_int(&mut self, name: String, value: i32) {
        self.ints.push(ParamSetItem::<i32> {
            name: name,
            values: vec![value],
            n_values: 1_usize,
            looked_up: false,
        });
    }
    pub fn add_ints(&mut self, name: String, values: Vec<i32>) {
        let n_values: usize = values.len();
        self.ints.push(ParamSetItem::<i32> {
            name: name,
            values: values,
            n_values: n_values,
            looked_up: false,
        });
    }
    pub fn add_bool(&mut self, name: String, value: bool) {
        self.bools.push(ParamSetItem::<bool> {
            name: name,
            values: vec![value],
            n_values: 1_usize,
            looked_up: false,
        });
    }
    pub fn add_point3f(&mut self, name: String, value: Point3f) {
        self.point3fs.push(ParamSetItem::<Point3f> {
            name: name,
            values: vec![value],
            n_values: 1_usize,
            looked_up: false,
        });
    }
    pub fn add_point3fs(&mut self, name: String, values: Vec<Float>) {
        let n_values: usize = values.len();
        let mut p_values: Vec<Point3f> = Vec::new();
        let n_points: usize = values.len() / 3_usize;
        assert!(n_values % 3 == 0, "point parameters need 3 coordinates");
        for i in 0..n_points {
            let x: Float = values[i * 3 + 0];
            let y: Float = values[i * 3 + 1];
            let z: Float = values[i * 3 + 2];
            p_values.push(Point3f { x: x, y: y, z: z });
        }
        self.point3fs.push(ParamSetItem::<Point3f> {
            name: name,
            values: p_values,
            n_values: n_points,
            looked_up: false,
        });
    }
    pub fn add_string(&mut self, name: String, value: String) {
        self.strings.push(ParamSetItem::<String> {
            name: name,
            values: vec![value],
            n_values: 1_usize,
            looked_up: false,
        });
    }
    pub fn add_texture(&mut self, name: String, value: String) {
        self.textures.push(ParamSetItem::<String> {
            name: name,
            values: vec![value],
            n_values: 1_usize,
            looked_up: false,
        });
    }
    pub fn add_vector3f(&mut self, name: String, value: Vector3f) {
        self.vector3fs.push(ParamSetItem::<Vector3f> {
            name: name,
            values: vec![value],
            n_values: 1_usize,
            looked_up: false,
        });
    }
    pub fn add_normal3f(&mut self, name: String, value: Normal3f) {
        self.normals.push(ParamSetItem::<Normal3f> {
            name: name,
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
            p_values.push(Normal3f { x: x, y: y, z: z });
        }
        self.normals.push(ParamSetItem::<Normal3f> {
            name: name,
            values: p_values,
            n_values: n_normals,
            looked_up: false,
        });
    }
    pub fn add_rgb_spectrum(&mut self, name: String, value: Spectrum) {
        self.spectra.push(ParamSetItem::<Spectrum> {
            name: name,
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
            name: name,
            values: s,
            n_values: n_values,
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
                values: values,
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
                values: values,
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
                values: values,
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
                values: values,
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
                values: values,
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
                values: values,
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
                values: values,
                n_values: s.n_values,
                looked_up: false,
            });
        }
    }
    pub fn find_one_float(&self, name: String, d: Float) -> Float {
        for v in &self.floats {
            if v.name == name && v.n_values == 1 {
                // v.looked_up = true;
                return v.values[0];
            }
        }
        d
    }
    pub fn find_one_int(&self, name: String, d: i32) -> i32 {
        for v in &self.ints {
            if v.name == name && v.n_values == 1 {
                // v.looked_up = true;
                return v.values[0];
            }
        }
        d
    }
    pub fn find_one_bool(&self, name: String, d: bool) -> bool {
        for v in &self.bools {
            if v.name == name && v.n_values == 1 {
                // v.looked_up = true;
                return v.values[0];
            }
        }
        d
    }
    pub fn find_one_point3f(&self, name: String, d: Point3f) -> Point3f {
        for v in &self.point3fs {
            if v.name == name && v.n_values == 1 {
                // v.looked_up = true;
                return v.values[0];
            }
        }
        d
    }
    pub fn find_one_vector3f(&self, name: String, d: Vector3f) -> Vector3f {
        for v in &self.vector3fs {
            if v.name == name && v.n_values == 1 {
                // v.looked_up = true;
                return v.values[0];
            }
        }
        d
    }
    pub fn find_one_spectrum(&self, name: String, d: Spectrum) -> Spectrum {
        for v in &self.spectra {
            if v.name == name && v.n_values == 1 {
                // v.looked_up = true;
                return v.values[0];
            }
        }
        d
    }
    pub fn find_one_string(&self, name: String, d: String) -> String {
        for v in &self.strings {
            if v.name == name && v.n_values == 1 {
                // v.looked_up = true;
                return v.values[0].clone();
            }
        }
        d
    }
    pub fn find_one_filename(&self, name: String, d: String) -> String {
        let filename: String = self.find_one_string(name, String::new());
        if filename == String::new() {
            return d;
        }
        // TODO: filename = AbsolutePath(ResolveFilename(filename));
        filename
    }
    pub fn find_texture(&self, name: String) -> String {
        let d: String = String::new();
        lookup_one(&self.textures, name, d)
    }
    pub fn find_int(&self, name: String) -> Vec<i32> {
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
    pub fn find_float(&self, name: String) -> Vec<Float> {
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
    pub fn find_point2f(&self, name: String) -> Vec<Point2f> {
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
    pub fn find_vector2f(&self, name: String) -> Vec<Vector2f> {
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
    pub fn find_point3f(&self, name: String) -> Vec<Point3f> {
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
    pub fn find_vector3f(&self, name: String) -> Vec<Vector3f> {
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
    pub fn find_normal3f(&self, name: String) -> Vec<Normal3f> {
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
    pub fn find_spectrum(&self, name: String) -> Vec<Spectrum> {
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
    pub float_textures: HashMap<String, Arc<Texture<Float> + Send + Sync>>,
    pub spectrum_textures: HashMap<String, Arc<Texture<Spectrum> + Send + Sync>>,
    pub geom_params: ParamSet,
    pub material_params: ParamSet,
}

impl TextureParams {
    pub fn new(
        geom_params: ParamSet,
        material_params: ParamSet,
        f_tex: HashMap<String, Arc<Texture<Float> + Send + Sync>>,
        s_tex: HashMap<String, Arc<Texture<Spectrum> + Send + Sync>>,
    ) -> Self {
        TextureParams {
            float_textures: f_tex,
            spectrum_textures: s_tex,
            geom_params: geom_params,
            material_params: material_params,
        }
    }
    pub fn get_spectrum_texture(
        &mut self,
        n: String,
        def: Spectrum,
    ) -> Arc<Texture<Spectrum> + Send + Sync> {
        let mut name: String = self.geom_params.find_texture(n.clone());
        if name == String::new() {
            name = self.material_params.find_texture(n.clone());
        }
        if name != String::new() {
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
        n: String,
    ) -> Option<Arc<Texture<Spectrum> + Send + Sync>> {
        let mut name: String = self.geom_params.find_texture(n.clone());
        if name == String::new() {
            name = self.material_params.find_texture(n.clone());
        }
        if name != String::new() {
            match self.spectrum_textures.get(name.as_str()) {
                Some(_spectrum_texture) => {}
                None => {
                    println!(
                        "Couldn't find spectrum texture named \"{}\" for parameter \"{}\"",
                        name, n
                    );
                    return None;
                }
            }
        }
        let mut val: Vec<Spectrum> = self.material_params.find_spectrum(n.clone());
        if val.len() == 0_usize {
            val = self.geom_params.find_spectrum(n.clone());
        }
        if val.len() == 0_usize {
            None
        } else {
            Some(Arc::new(ConstantTexture { value: val[0] }))
        }
    }
    pub fn get_float_texture(
        &mut self,
        n: String,
        def: Float,
    ) -> Arc<Texture<Float> + Send + Sync> {
        let mut name: String = self.geom_params.find_texture(n.clone());
        if name == String::new() {
            name = self.material_params.find_texture(n.clone());
        }
        if name != String::new() {
            match self.float_textures.get(name.as_str()) {
                Some(_float_texture) => {}
                None => {
                    panic!(
                        "Couldn't find float texture named \"{}\" for parameter \"{}\"",
                        name, n
                    );
                }
            }
        }
        let mut val: Float = self.material_params.find_one_float(n.clone(), def);
        val = self.geom_params.find_one_float(n.clone(), val);
        Arc::new(ConstantTexture { value: val })
    }
    pub fn get_float_texture_or_null(
        &mut self,
        n: String,
    ) -> Option<Arc<Texture<Float> + Send + Sync>> {
        let mut name: String = self.geom_params.find_texture(n.clone());
        if name == String::new() {
            name = self.material_params.find_texture(n.clone());
        }
        if name != String::new() {
            match self.float_textures.get(name.as_str()) {
                Some(_float_texture) => {}
                None => {
                    println!(
                        "Couldn't find float texture named \"{}\" for parameter \"{}\"",
                        name, n
                    );
                    return None;
                }
            }
        }
        let mut val: Vec<Float> = self.material_params.find_float(n.clone());
        if val.len() == 0_usize {
            val = self.geom_params.find_float(n.clone());
        }
        if val.len() == 0_usize {
            None
        } else {
            Some(Arc::new(ConstantTexture { value: val[0] }))
        }
    }
    pub fn find_float(&mut self, name: String, d: Float) -> Float {
        self.geom_params.find_one_float(
            name.clone(),
            self.material_params.find_one_float(name.clone(), d),
        )
    }
    pub fn find_string(&mut self, name: String, d: String) -> String {
        self.geom_params.find_one_string(
            name.clone(),
            self.material_params.find_one_string(name.clone(), d),
        )
    }
    pub fn find_filename(&mut self, name: String, d: String) -> String {
        self.geom_params.find_one_filename(
            name.clone(),
            self.material_params.find_one_filename(name.clone(), d),
        )
    }
    pub fn find_int(&mut self, name: String, d: i32) -> i32 {
        self.geom_params.find_one_int(
            name.clone(),
            self.material_params.find_one_int(name.clone(), d),
        )
    }
    pub fn find_bool(&mut self, name: String, d: bool) -> bool {
        self.geom_params.find_one_bool(
            name.clone(),
            self.material_params.find_one_bool(name.clone(), d),
        )
    }
    pub fn find_vector3f(&mut self, name: String, d: Vector3f) -> Vector3f {
        self.geom_params.find_one_vector3f(
            name.clone(),
            self.material_params.find_one_vector3f(name.clone(), d),
        )
    }
}

/// Replaces a macro on the C++ side.
pub fn lookup_one<T>(vec: &Vec<ParamSetItem<T>>, name: String, d: T) -> T
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
