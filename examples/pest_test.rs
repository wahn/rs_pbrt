#![recursion_limit="2000"]
#![feature(drop_types_in_const)]

#[macro_use]
extern crate pest;
extern crate getopts;
extern crate pbrt;

use pbrt::{AnimatedTransform, Bounds2f, BoxFilter, Camera, Checkerboard2DTexture, ConstantTexture,
           DistantLight, Film, Filter, Float, GeometricPrimitive, GraphicsState, ImageTexture,
           ImageWrap, Material, MatteMaterial, Matrix4x4, MirrorMaterial, ParamSet,
           PerspectiveCamera, PlanarMapping2D, Point2i, Point2f, Point3f, RenderOptions, Sampler,
           Spectrum, Sphere, Texture, TextureMapping2D, TextureParams, Transform, TransformSet,
           UVMapping2D, Vector2f, Vector3f, ZeroTwoSequenceSampler};
// parser
use pest::prelude::*;
// getopts
use getopts::Options;
// std
use std::str::FromStr;
use std::collections::{HashMap, LinkedList};
use std::env;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::sync::Arc;

pub const VERSION: &'static str = env!("CARGO_PKG_VERSION");

static mut SEARCH_DIRECTORY: Option<Box<PathBuf>> = None;
static mut CUR_TRANSFORM: TransformSet = TransformSet {
    t: [Transform {
        m: Matrix4x4 {
            m: [[1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0]],
        },
        m_inv: Matrix4x4 {
            m: [[1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0]],
        },
    }; 2],
};
static mut NAMED_COORDINATE_SYSTEMS: Option<Box<HashMap<&str, TransformSet>>> = None;
static mut RENDER_OPTIONS: Option<Box<RenderOptions>> = None;
static mut GRAPHICS_STATE: Option<Box<GraphicsState>> = None;
static mut PUSHED_GRAPHICS_STATES: Option<Box<Vec<GraphicsState>>> = None;
static mut PUSHED_TRANSFORMS: Option<Box<Vec<TransformSet>>> = None;
// not used in original C++ code:
static mut PARAM_SET: Option<Box<ParamSet>> = None;

#[derive(Debug, PartialEq)]
pub enum IntNode {
    Ints(LinkedList<IntNode>),
    OneInt(i32),
}

#[derive(Debug, PartialEq)]
pub enum FloatNode {
    Floats(LinkedList<FloatNode>),
    OneFloat(Float),
}

impl_rdp! {
    grammar! {
        pbrt = _{ whitespace? ~ statement* ~ last_statement }
        statement = { look_at | translate | rotate | named_statement | keyword }
        named_statement = { camera |
                            pixel_filter |
                            sampler |
                            film |
                            integrator |
                            coord_sys_transform |
                            light_source |
                            texture |
                            material |
                            shape }
        parameter = { float_param |
                      string_param |
                      integer_param |
                      point_param |
                      vector_param |
                      rgb_param |
                      spectrum_param |
                      texture_param }
        float_param = { ["\"float"] ~ ident ~ ["\""] ~ lbrack ~ number+ ~ rbrack }
        string_param = { (["\"string"] ~ ident ~ ["\""] ~ lbrack ~ string ~ rbrack) |
                         (["\"string"] ~ ident ~ ["\""] ~ string) }
        integer_param = { ["\"integer"] ~ ident ~ ["\""] ~ lbrack ~ integer+ ~ rbrack }
        point_param = { ["\"point"] ~ ident ~ ["\""] ~ lbrack ~ number+ ~ rbrack }
        vector_param = { ["\"vector"] ~ ident ~ ["\""] ~ lbrack ~ number ~ number ~ number ~ rbrack }
        rgb_param = { (["\"rgb"] ~ ident ~ ["\""] ~ lbrack ~ number ~ number ~ number ~ rbrack) |
                      (["\"color"] ~ ident ~ ["\""] ~ lbrack ~ number ~ number ~ number ~ rbrack) }
        spectrum_param = { ["\"spectrum\""] ~ string }
        texture_param = { ["\"texture"] ~ ident ~ ["\""] ~ string }
        // Translate x y z
        translate = { ["Translate"] ~
                   // followed by 3 numbers:
                   number ~ number ~ number
        }
        // Rotate angle x y z
        rotate = { ["Rotate"] ~
                   // followed by 4 numbers:
                   number ~ number ~ number ~ number
        }
        // LookAt eye_x eye_y eye_z look_x look_y look_z up_x up_y up_z
        look_at = { ["LookAt"] ~
                    // followed by 9 numbers:

                    // eye_x eye_y eye_z
                    number ~ number ~ number ~
                    // look_x look_y look_z
                    number ~ number ~ number ~
                    // up_x up_y up_z
                    number ~ number ~ number
        }
        // Camera "perspective" "float fov" [90] ...
        camera = { ["Camera"] ~ string ~ parameter* }
        // PixelFilter "mitchell" "float xwidth" [2] "float ywidth" [2]
        pixel_filter = { ["PixelFilter"] ~ string ~ parameter* }
        // Sampler "halton"
        sampler = { ["Sampler"] ~ string ~ parameter* }
        // Film "image" "string filename" ["..."] ...
        film = { ["Film"] ~ string ~ parameter* }
        // Integrator "path" "integer maxdepth" [5]
        integrator = { ["Integrator"] ~ string ~ parameter* }
        // CoordSysTransform "camera"
        coord_sys_transform = { ["CoordSysTransform"] ~ string }
        // LightSource "point" "rgb I" [ .5 .5 .5 ]
        light_source = { ["LightSource"] ~ string ~ parameter* }
        // Texture "mydiffuse" "spectrum" "imagemap" "string filename" "image.tga"
        texture = { ["Texture"] ~ string ~ string ~ string ~ parameter* }
        // Material "matte" "texture Kd" "mydiffuse"
        material = { ["Material"] ~ string ~ parameter* }
        // Shape "sphere" "float radius" [0.25]
        shape = { ["Shape"] ~ string ~ parameter* }
        // keywords
        keyword = {
            (["Accelerator"] |
             ["ActiveTransform"] |
             ["All"] |
             ["AreaLightSource"] |
             attribute_begin |
             attribute_end |
             ["ConcatTransform"] |
             ["CoordinateSystem"] |
             ["EndTime"] |
             ["Identity"] |
             ["Include"] |
             ["MakeNamedMedium"] |
             ["MakeNamedMaterial"] |
             ["MediumInterface"] |
             ["NamedMaterial"] |
             ["ObjectBegin"] |
             ["ObjectEnd"] |
             ["ObjectInstance"] |
             ["ReverseOrientation"] |
             ["Scale"] |
             ["StartTime"] |
             ["TransformBegin"] |
             ["TransformEnd"] |
             ["TransformTimes"] |
             ["Transform"] |
             world_begin
            )
        }
        attribute_begin = { ["AttributeBegin"] }
        attribute_end = { ["AttributeEnd"] }
        world_begin = { ["WorldBegin"] }
        // IDENT [a-zA-Z_][a-zA-Z_0-9]*
        ident =  { (['a'..'z'] | ['A'..'Z'] | ["_"]) ~
                   (['a'..'z'] | ['A'..'Z'] | ["_"] | ["-"] | ['0'..'9'])* }
        string = { (["\""] ~ ident ~ ["\""]) | (["\""] ~ filename ~ ["\""]) }
        filename = { (['a'..'z'] | ['A'..'Z'] | ["_"]) ~ // TODO: can be a full path
                     (['a'..'z'] | ['A'..'Z'] | ["_"] | ["-"] | ["."] | ["/"] | ['0'..'9'])* }
        // "[" { return LBRACK; }
        lbrack = { ["["] }
        // "]" { return RBRACK; }
        rbrack = { ["]"] }
        // NUMBER [-+]?([0-9]+|(([0-9]+\.[0-9]*)|(\.[0-9]+)))([eE][-+]?[0-9]+)?
        number = @{
            (["-"] | ["+"])? ~ // optional sign, followed by
            (
                (
                    (["."] ~ ['0'..'9']+) // dot and digits
                        | // or
                    (['0'..'9']+ ~ ["."] ~ ['0'..'9']*) // digits, dot, and (optional digits)
                )
                    | // or
                ['0'..'9']+ // just digits
            ) ~ ( // followed by (optional)
                (["e"] | ["E"]) ~ // 'e' or 'E', followed by
                (["-"] | ["+"])? ~ // optional sign, followed by
                ['0'..'9']+ // digits
            )?
        }
        integer = @{
            (["-"] | ["+"])? ~ // optional sign, followed by
                (
                    ['1'..'9'] ~ // at least one non-zero digit, followed by
                    ['0'..'9']* // just digits
                )
                    | // or
                ['0'..'9'] // single digit
        }
        last_statement = @{ whitespace? ~ ["WorldEnd"] ~ whitespace? }
        whitespace = _{ ([" "] | ["\t"] | ["\r"] | ["\n"]) }
    }
    process! {
        main(&self) -> () {
            (_list: _pbrt()) => {
            }
        }
        _pbrt(&self) -> () {
            (_head: statement, _tail: _statement()) => {},
            (_l: last_statement) => { println!("WorldEnd"); },
        }
        // statements
        _statement(&self) -> () {
            (_head: look_at, _tail: _look_at()) => {},
            (_head: translate, _tail: _translate()) => {},
            (_head: rotate, _tail: _rotate()) => {},
            (_head: named_statement, _tail: _named_statement()) => {},
            (_head: keyword, _tail: _keyword()) => {},
        }
        _look_at(&self) -> () {
            (eye_x: _number(), eye_y: _number(), eye_z: _number(),
             look_x: _number(), look_y: _number(), look_z: _number(),
             up_x: _number(), up_y: _number(), up_z: _number()) => {
                println!("LookAt {} {} {} {} {} {} {} {} {}",
                         eye_x, eye_y, eye_z,
                         look_x, look_y, look_z,
                         up_x, up_y, up_z,);
                let pos: Point3f = Point3f { x: eye_x, y: eye_y, z: eye_z, };
                let look: Point3f = Point3f { x: look_x, y: look_y, z: look_z, };
                let up: Vector3f = Vector3f { x: up_x, y: up_y, z: up_z, };
                let look_at: Transform = Transform::look_at(pos, look, up);
                unsafe {
                    CUR_TRANSFORM.t[0] = CUR_TRANSFORM.t[0] * look_at;
                    CUR_TRANSFORM.t[1] = CUR_TRANSFORM.t[1] * look_at;
                    // println!("CUR_TRANSFORM: {:?}", CUR_TRANSFORM);
                }
                self._pbrt();
            }
        }
        _translate(&self) -> () {
            (x: _number(), y: _number(), z: _number()) => {
                println!("Translate {} {} {}", x, y, z);
                let translate: Transform = Transform::translate(Vector3f { x: x, y: y, z: z, });
                unsafe {
                    CUR_TRANSFORM.t[0] = CUR_TRANSFORM.t[0] * translate;
                    CUR_TRANSFORM.t[1] = CUR_TRANSFORM.t[1] * translate;
                    // println!("CUR_TRANSFORM: {:?}", CUR_TRANSFORM);
                }
                self._pbrt();
            }
        }
        _rotate(&self) -> () {
            (angle: _number(), x: _number(), y: _number(), z: _number()) => {
                println!("Rotate {} {} {} {}",
                         angle, x, y, z);
                let rotate: Transform = Transform::rotate(angle, Vector3f { x: x, y: y, z: z, });
                unsafe {
                    CUR_TRANSFORM.t[0] = CUR_TRANSFORM.t[0] * rotate;
                    CUR_TRANSFORM.t[1] = CUR_TRANSFORM.t[1] * rotate;
                    // println!("CUR_TRANSFORM: {:?}", CUR_TRANSFORM);
                }
                self._pbrt();
            }
        }
        // named statements
        _named_statement(&self) -> () {
            (_head: camera, _tail: _camera()) => {},
            (_head: pixel_filter, _tail: _pixel_filter()) => {},
            (_head: sampler, _tail: _sampler()) => {},
            (_head: film, _tail: _film()) => {},
            (_head: integrator, _tail: _integrator()) => {},
            (_head: coord_sys_transform, _tail: _coord_sys_transform()) => {},
            (_head: light_source, _tail: _light_source()) => {},
            (_head: texture, _tail: _texture()) => {},
            (_head: material, _tail: _material()) => {},
            (_head: shape, _tail: _shape()) => {},
        }
        _camera(&self) -> () {
            (name: _string(), optional_parameters) => {
                unsafe {
                    if let Some(ref mut ro) = RENDER_OPTIONS {
                        ro.camera_name = name;
                        ro.camera_to_world.t[0] =
                            Transform::inverse(CUR_TRANSFORM.t[0]);
                        ro.camera_to_world.t[1] =
                            Transform::inverse(CUR_TRANSFORM.t[1]);
                        if let Some(ref mut named_coordinate_systems) = NAMED_COORDINATE_SYSTEMS {
                            named_coordinate_systems.insert("camera",
                                                            TransformSet {
                                                                t: [ro.camera_to_world.t[0],
                                                                    ro.camera_to_world.t[1]]
                                                            });
                            // println!("NAMED_COORDINATE_SYSTEMS: {:?}", named_coordinate_systems);
                        }
                        // println!("ro.camera_to_world: {:?}",
                        //          ro.camera_to_world);
                    }
                    if let Some(ref mut param_set) = PARAM_SET {
                        param_set.reset(String::from("Camera"),
                                        String::from(""),
                                        String::from(""),
                                        String::from(""));
                    }
                }
                if optional_parameters.rule == Rule::parameter {
                    self._parameter();
                } else if optional_parameters.rule == Rule::statement {
                    self._statement();
                } else if optional_parameters.rule == Rule::last_statement {
                    world_end();
                } else {
                    println!("ERROR: parameter expected, {:?} found ...", optional_parameters);
                }
            },
        }
        _pixel_filter(&self) -> () {
            (name: _string(), optional_parameters) => {
                unsafe {
                    if let Some(ref mut ro) = RENDER_OPTIONS {
                        ro.filter_name = name;
                        if optional_parameters.rule == Rule::statement ||
                            optional_parameters.rule == Rule::last_statement {
                            println!("PixelFilter \"{}\" ", ro.filter_name);
                        }
                    }
                    if let Some(ref mut param_set) = PARAM_SET {
                        param_set.reset(String::from("PixelFilter"),
                                        String::from(""),
                                        String::from(""),
                                        String::from(""));
                    }
                }
                if optional_parameters.rule == Rule::parameter {
                    self._parameter();
                } else if optional_parameters.rule == Rule::statement {
                    self._statement();
                } else if optional_parameters.rule == Rule::last_statement {
                    world_end();
                } else {
                    println!("ERROR: parameter expected, {:?} found ...", optional_parameters);
                }
            },
        }
        _sampler(&self) -> () {
            (name: _string(), optional_parameters) => {
                unsafe {
                    if let Some(ref mut ro) = RENDER_OPTIONS {
                        ro.sampler_name = name;
                        if optional_parameters.rule == Rule::statement ||
                            optional_parameters.rule == Rule::last_statement {
                            println!("Sampler \"{}\" ", ro.sampler_name);
                        }
                    }
                    if let Some(ref mut param_set) = PARAM_SET {
                        param_set.reset(String::from("Sampler"),
                                        String::from(""),
                                        String::from(""),
                                        String::from(""));
                    }
                }
                if optional_parameters.rule == Rule::parameter {
                    self._parameter();
                } else if optional_parameters.rule == Rule::statement {
                    self._statement();
                } else if optional_parameters.rule == Rule::last_statement {
                    world_end();
                } else {
                    println!("ERROR: parameter expected, {:?} found ...", optional_parameters);
                }
            },
        }
        _film(&self) -> () {
            (name: _string(), optional_parameters) => {
                unsafe {
                    if let Some(ref mut ro) = RENDER_OPTIONS {
                        ro.film_name = name;
                    }
                    if let Some(ref mut param_set) = PARAM_SET {
                        param_set.reset(String::from("Film"),
                                        String::from(""),
                                        String::from(""),
                                        String::from(""));
                    }
                }
                if optional_parameters.rule == Rule::parameter {
                    self._parameter();
                } else if optional_parameters.rule == Rule::statement {
                    self._statement();
                } else if optional_parameters.rule == Rule::last_statement {
                    world_end();
                } else {
                    println!("ERROR: parameter expected, {:?} found ...", optional_parameters);
                }
            },
        }
        _integrator(&self) -> () {
            (name: _string(), optional_parameters) => {
                unsafe {
                    if let Some(ref mut ro) = RENDER_OPTIONS {
                        ro.integrator_name = name;
                    }
                    if let Some(ref mut param_set) = PARAM_SET {
                        param_set.reset(String::from("Integrator"),
                                        String::from(""),
                                        String::from(""),
                                        String::from(""));
                    }
                }
                if optional_parameters.rule == Rule::parameter {
                    self._parameter();
                } else if optional_parameters.rule == Rule::statement {
                    self._statement();
                } else if optional_parameters.rule == Rule::last_statement {
                    world_end();
                } else {
                    println!("ERROR: parameter expected, {:?} found ...", optional_parameters);
                }
            },
        }
        _coord_sys_transform(&self) -> () {
            (name: _string()) => {
                println!("CoordSysTransform \"{}\" ", name);
                unsafe {
                    if let Some(ref mut named_coordinate_systems) = NAMED_COORDINATE_SYSTEMS {
                        match named_coordinate_systems.get(name.as_str()) {
                            Some(transform_set) => {
                                CUR_TRANSFORM.t[0] = transform_set.t[0];
                                CUR_TRANSFORM.t[1] = transform_set.t[1];
                                // println!("CUR_TRANSFORM: {:?}", CUR_TRANSFORM);
                            },
                            None => {
                                println!("Couldn't find named coordinate system \"{}\"", name);
                            },
                        };
                    }
                }
                self._pbrt();
            },
        }
        _light_source(&self) -> () {
            (name: _string(), optional_parameters) => {
                unsafe {
                    if let Some(ref mut param_set) = PARAM_SET {
                        param_set.reset(String::from("LightSource"),
                                        String::from(name),
                                        String::from(""),
                                        String::from(""));
                    }
                }
                if optional_parameters.rule == Rule::parameter {
                    self._parameter();
                } else if optional_parameters.rule == Rule::statement {
                    self._statement();
                } else if optional_parameters.rule == Rule::last_statement {
                    world_end();
                } else {
                    println!("ERROR: parameter expected, {:?} found ...", optional_parameters);
                }
            },
        }
        _texture(&self) -> () {
            (name: _string(), tex_type: _string(), tex_name: _string(), optional_parameters) => {
                unsafe {
                    if optional_parameters.rule == Rule::statement ||
                        optional_parameters.rule == Rule::last_statement
                    {
                        println!("Texture \"{}\" \"{}\" \"{}\" ",
                                 name,
                                 tex_type,
                                 tex_name);
                    }
                    if let Some(ref mut param_set) = PARAM_SET {
                        param_set.reset(String::from("Texture"),
                                        String::from(name),
                                        String::from(tex_type),
                                        String::from(tex_name));
                    }
                }
                if optional_parameters.rule == Rule::parameter {
                    self._parameter();
                } else if optional_parameters.rule == Rule::statement {
                    self._statement();
                } else if optional_parameters.rule == Rule::last_statement {
                    world_end();
                } else {
                    println!("ERROR: parameter expected, {:?} found ...", optional_parameters);
                }
            },
        }
        _material(&self) -> () {
            (name: _string(), optional_parameters) => {
                unsafe {
                    if let Some(ref mut param_set) = PARAM_SET {
                        if optional_parameters.rule == Rule::statement ||
                            optional_parameters.rule == Rule::last_statement
                        {
                            println!("Material \"{}\" ", name.clone());
                        }
                        // pbrtMaterial (api.cpp:1082)
                        unsafe {
                            if let Some(ref mut graphics_state) = GRAPHICS_STATE {
                                graphics_state.material = name.clone();
                                graphics_state.current_named_material = String::new();
                            }
                        }
                        param_set.reset(String::from("Material"),
                                        String::from(name),
                                        String::from(""),
                                        String::from(""));
                    }
                }
                if optional_parameters.rule == Rule::parameter {
                    self._parameter();
                } else if optional_parameters.rule == Rule::statement {
                    self._statement();
                } else if optional_parameters.rule == Rule::last_statement {
                    world_end();
                } else {
                    println!("ERROR: parameter expected, {:?} found ...", optional_parameters);
                }
            },
        }
        _shape(&self) -> () {
            (name: _string(), optional_parameters) => {
                unsafe {
                    if let Some(ref mut param_set) = PARAM_SET {
                        if optional_parameters.rule == Rule::statement ||
                            optional_parameters.rule == Rule::last_statement
                        {
                            println!("Shape \"{}\" ", name);
                            // pbrtShape (api.cpp:1153)
                            // TODO: if (!curTransform.IsAnimated()) { ... }
                            // TODO: transformCache.Lookup(curTransform[0], &ObjToWorld, &WorldToObj);
                            let obj_to_world: Transform = Transform {
                                m: CUR_TRANSFORM.t[0].m,
                                m_inv: CUR_TRANSFORM.t[0].m_inv,
                            };
                            let world_to_obj: Transform =Transform {
                                m: CUR_TRANSFORM.t[0].m_inv,
                                m_inv: CUR_TRANSFORM.t[0].m,
                            };
                            // MakeShapes (api.cpp:296)
                            if name == String::from("sphere") {
                                // CreateSphereShape
                                let radius: Float = 1.0; // default
                                let z_min: Float = -radius; // default
                                let z_max: Float = radius; // default
                                let phi_max: Float = 360.0; // default
                                let sphere = Arc::new(Sphere::new(obj_to_world,
                                                                  world_to_obj,
                                                                  false,
                                                                  false,
                                                                  radius,
                                                                  z_min,
                                                                  z_max,
                                                                  phi_max));
                                print!("Sphere {{ object_to_world: {:?}, world_to_object: {:?}, ",
                                       obj_to_world,
                                       world_to_obj);
                                println!("radius: {}, z_min: {}, z_max: {}, phi_max: {} }}",
                                         radius,
                                         z_min,
                                         z_max,
                                         phi_max);
                                let mtl: Arc<Material + Send + Sync> = create_material();
                                let geo_prim = Arc::new(GeometricPrimitive::new(sphere, mtl.clone()));
                                if let Some(ref mut ro) = RENDER_OPTIONS {
                                    ro.primitives.push(geo_prim.clone());
                                }
                            } else if name == String::from("cylinder") {
                                println!("TODO: CreateCylinderShape");
                            } else if name == String::from("disk") {
                                println!("TODO: CreateDiskShape");
                            } else if name == String::from("cone") {
                                println!("TODO: CreateConeShape");
                            } else if name == String::from("paraboloid") {
                                println!("TODO: CreateParaboloidShape");
                            } else if name == String::from("hyperboloid") {
                                    println!("TODO: CreateHyperboloidShape");
                            } else if name == String::from("curve") {
                                println!("TODO: CreateCurveShape");
                            } else if name == String::from("trianglemesh") {
                                println!("TODO: CreateTriangleMeshShape");
                            } else if name == String::from("plymesh") {
                                println!("TODO: CreatePLYMesh");
                            } else if name == String::from("heightfield") {
                                println!("TODO: CreateHeightfield");
                            } else if name == String::from("loopsubdiv") {
                                println!("TODO: CreateLoopSubdiv");
                            } else if name == String::from("nurbs") {
                                println!("TODO: CreateNURBS");
                            } else {
                                panic!("Shape \"{}\" unknown.", name);
                            }
                        }
                        param_set.reset(String::from("Shape"),
                                        String::from(name),
                                        String::from(""),
                                        String::from(""));
                    }
                }
                if optional_parameters.rule == Rule::parameter {
                    self._parameter();
                } else if optional_parameters.rule == Rule::statement {
                    self._statement();
                } else if optional_parameters.rule == Rule::last_statement {
                    world_end();
                } else {
                    println!("ERROR: parameter expected, {:?} found ...", optional_parameters);
                }
            },
        }
        // parameters
        _parameter(&self) -> () {
            (_head: float_param, tail: _float_param()) => {
                let (string, numbers) = tail;
                unsafe {
                    if let Some(ref mut param_set) = PARAM_SET {
                        if numbers.len() == 1 {
                            param_set.add_float(string, numbers[0]);
                        } else {
                            param_set.add_floats(string, numbers);
                        }
                    }
                }
                self._parameter();
            },
            (_head: string_param, tail: _string_param()) => {
                let (string1, string2) = tail;
                unsafe {
                    if let Some(ref mut param_set) = PARAM_SET {
                        param_set.add_string(string1, string2);
                    }
                }
                self._parameter();
            },
            (_head: integer_param, tail: _integer_param()) => {
                let (string, numbers) = tail;
                unsafe {
                    if let Some(ref mut param_set) = PARAM_SET {
                        if numbers.len() == 1 {
                            param_set.add_int(string, numbers[0]);
                        } else {
                            param_set.add_ints(string, numbers);
                        }
                    }
                }
                self._parameter();
            },
            (_head: point_param, tail: _point_param()) => {
                let (string, numbers) = tail;
                unsafe {
                    if let Some(ref mut param_set) = PARAM_SET {
                        if numbers.len() == 3 {
                            param_set.add_point3f(string,
                                                  Point3f {
                                                      x: numbers[0],
                                                      y: numbers[1],
                                                      z: numbers[2],
                                                  });
                        } else {
                            param_set.add_point3fs(string, numbers);
                        }
                    }
                }
                self._parameter();
            },
            (_head: vector_param, tail: _vector_param()) => {
                let (string, number1, number2, number3) = tail;
                unsafe {
                    if let Some(ref mut param_set) = PARAM_SET {
                        param_set.add_vector3f(string,
                                              Vector3f {
                                                  x: number1,
                                                  y: number2,
                                                  z: number3,
                                              });
                    }
                }
                self._parameter();
            },
            (_head: rgb_param, tail: _rgb_param()) => {
                let (string, number1, number2, number3) = tail;
                unsafe {
                    if let Some(ref mut param_set) = PARAM_SET {
                        param_set.add_rgb_spectrum(string,
                                                   Spectrum {
                                                       c: [number1, number2, number3],
                                                   });
                    }
                }
                self._parameter();
            },
            (_head: spectrum_param, tail: _spectrum_param()) => {
                let string = tail;
                print!("\"spectrum\" {} ", string);
                // unsafe {
                //     if let Some(ref mut param_set) = PARAM_SET {
                //         param_set.add_string(string1, string2);
                //     }
                // }
                self._parameter();
            },
            (_head: texture_param, tail: _texture_param()) => {
                let (string1, string2) = tail;
                unsafe {
                    if let Some(ref mut param_set) = PARAM_SET {
                        param_set.add_texture(string1, string2);
                    }
                }
                self._parameter();
            },
            (optional_parameters) => {
                if optional_parameters.rule == Rule::statement ||
                    optional_parameters.rule == Rule::last_statement {
                    unsafe {
                        if let Some(ref mut param_set) = PARAM_SET {
                            let mut name: String = String::new();
                            name.push_str(param_set.name.as_str());
                            if param_set.key_word == String::from("Camera") {
                                if let Some(ref mut ro) = RENDER_OPTIONS {
                                    println!("Camera \"{}\" ", ro.camera_name);
                                    ro.camera_params.copy_from(param_set);
                                    print_params(&ro.camera_params);
                                }
                            } else if param_set.key_word == String::from("PixelFilter") {
                                if let Some(ref mut ro) = RENDER_OPTIONS {
                                    println!("PixelFilter \"{}\" ", ro.filter_name);
                                    ro.filter_params.copy_from(param_set);
                                    print_params(&ro.filter_params);
                                }
                            } else if param_set.key_word == String::from("Sampler") {
                                if let Some(ref mut ro) = RENDER_OPTIONS {
                                    println!("Sampler \"{}\" ", ro.sampler_name);
                                    ro.sampler_params.copy_from(param_set);
                                    print_params(&ro.sampler_params);
                                }
                            } else if param_set.key_word == String::from("Film") {
                                if let Some(ref mut ro) = RENDER_OPTIONS {
                                    println!("Film \"{}\" ", ro.film_name);
                                    ro.film_params.copy_from(param_set);
                                    print_params(&ro.film_params);
                                }
                            } else if param_set.key_word == String::from("Integrator") {
                                if let Some(ref mut ro) = RENDER_OPTIONS {
                                    println!("Integrator \"{}\" ", ro.integrator_name);
                                    ro.integrator_params.copy_from(param_set);
                                    print_params(&ro.integrator_params);
                                }
                            } else if param_set.key_word == String::from("LightSource") {
                                if let Some(ref mut ro) = RENDER_OPTIONS {
                                    println!("LightSource \"{}\" ", param_set.name);
                                    print_params(&param_set);
                                    // MakeLight (api.cpp:591)
                                    if param_set.name == String::from("point") {
                                        println!("TODO: CreatePointLight");
                                    } else if param_set.name == String::from("spot") {
                                        println!("TODO: CreateSpotLight");
                                    } else if param_set.name == String::from("goniometric") {
                                        println!("TODO: CreateGoniometricLight");
                                    } else if param_set.name == String::from("projection") {
                                        println!("TODO: CreateProjectionLight");
                                    } else if param_set.name == String::from("distant") {
                                        // CreateDistantLight
                                        let l: Spectrum = param_set.find_one_spectrum(String::from("L"),
                                                                                      Spectrum::new(1.0
                                                                                                    as Float));
                                        let sc: Spectrum = param_set.find_one_spectrum(String::from("scale"),
                                                                                       Spectrum::new(1.0
                                                                                                     as Float));
                                        let from: Point3f = param_set.find_one_point3f(String::from("from"),
                                                                                       Point3f { x: 0.0,
                                                                                                 y: 0.0,
                                                                                                 z: 0.0 });
                                        let to: Point3f = param_set.find_one_point3f(String::from("to"),
                                                                                       Point3f { x: 0.0,
                                                                                                 y: 0.0,
                                                                                                 z: 0.0 });
                                        let dir: Vector3f = from - to;
                                        // return std::make_shared<DistantLight>(light2world, L * sc, dir);
                                        let distant_light =
                                            Arc::new(DistantLight::new(&CUR_TRANSFORM.t[0], &(l * sc), &dir));
                                        println!("{:?}", distant_light);
                                        ro.lights.push(distant_light);
                                    } else if param_set.name == String::from("infinite") {
                                        println!("TODO: CreateInfiniteLight");
                                    } else if param_set.name == String::from("exinfinite") {
                                        println!("TODO: CreateInfiniteLight");
                                    } else {
                                        panic!("MakeLight: unknown name {}", param_set.name);
                                    }
                                }
                                // let lt = make_light(name, params, CUR_TRANSFORM.t[0]);
                            } else if param_set.key_word == String::from("Texture") {
                                if let Some(ref mut ro) = RENDER_OPTIONS {
                                    println!("Texture \"{}\" \"{}\" \"{}\" ",
                                             param_set.name,
                                             param_set.tex_type,
                                             param_set.tex_name);
                                    print_params(&param_set);
                                    // pbrtTexture (api.cpp:1049)
                                    unsafe {
                                        if let Some(ref mut graphics_state) = GRAPHICS_STATE {
                                            let mut geom_params: ParamSet = ParamSet::default();
                                            let mut material_params: ParamSet = ParamSet::default();
                                            geom_params.copy_from(param_set);
                                            material_params.copy_from(param_set);
                                            let mut tp: TextureParams = TextureParams {
                                                float_textures: graphics_state.float_textures.clone(),
                                                spectrum_textures: graphics_state.spectrum_textures.clone(),
                                                geom_params: geom_params,
                                                material_params: material_params,
                                            };
                                            if param_set.tex_type == String::from("float") {
                                                println!("TODO: MakeFloatTexture");
                                            } else if param_set.tex_type == String::from("color") ||
                                                param_set.tex_type == String::from("spectrum") {
                                                    match graphics_state.spectrum_textures.get(param_set.name.as_str()) {
                                                        Some(_spectrum_texture) => {
                                                            println!("Texture \"{}\" being redefined",
                                                                     param_set.name);
                                                        },
                                                        None => {},
                                                    }
                                                    // MakeSpectrumTexture(texname, curTransform[0], tp);
                                                    if param_set.tex_name == String::from("constant") {
                                                        println!("TODO: CreateConstantSpectrumTexture");
                                                    } else if param_set.tex_name == String::from("scale") {
                                                        println!("TODO: CreateScaleSpectrumTexture");
                                                    } else if param_set.tex_name == String::from("mix") {
                                                        println!("TODO: CreateMixSpectrumTexture");
                                                    } else if param_set.tex_name == String::from("bilerp") {
                                                        println!("TODO: CreateBilerpSpectrumTexture");
                                                    } else if param_set.tex_name == String::from("imagemap") {
                                                        // CreateImageSpectrumTexture
                                                        let mut map: Option<Box<TextureMapping2D + Send + Sync>> = None;
                                                        let mapping: String =
                                                            tp.find_string(String::from("mapping"), String::from("uv"));
                                                        if mapping == String::from("uv") {
                                                            let su: Float = tp.find_float(String::from("uscale"), 1.0);
                                                            let sv: Float = tp.find_float(String::from("vscale"), 1.0);
                                                            let du: Float = tp.find_float(String::from("udelta"), 0.0);
                                                            let dv: Float = tp.find_float(String::from("vdelta"), 0.0);
                                                            map = Some(Box::new(UVMapping2D {
                                                                su: su,
                                                                sv: sv,
                                                                du: du,
                                                                dv: dv,
                                                            }));
                                                        } else if mapping == String::from("spherical") {
                                                            println!("TODO: SphericalMapping2D");
                                                        } else if mapping == String::from("cylindrical") {
                                                            println!("TODO: CylindricalMapping2D");
                                                        } else if mapping == String::from("planar") {
                                                            map = Some(Box::new(PlanarMapping2D {
                                                                vs: tp.find_vector3f(String::from("v1"),
                                                                                     Vector3f {
                                                                                         x: 1.0,
                                                                                         y: 0.0,
                                                                                         z: 0.0
                                                                                     }),
                                                                vt: tp.find_vector3f(String::from("v2"),
                                                                                     Vector3f {
                                                                                         x: 0.0,
                                                                                         y: 1.0,
                                                                                         z: 0.0
                                                                                     }),
                                                                ds: tp.find_float(String::from("udelta"),
                                                                                      0.0),
                                                                dt: tp.find_float(String::from("vdelta"),
                                                                                  0.0),
                                                            }));
                                                        } else {
                                                            panic!("2D texture mapping \"{}\" unknown",
                                                                   mapping);
                                                        }
                                                        // initialize _ImageTexture_ parameters
                                                        let max_aniso: Float =
                                                            tp.find_float(String::from("maxanisotropy"), 8.0);
                                                        let do_trilinear: bool =
                                                            tp.find_bool(String::from("trilinear"), false);
                                                        let wrap: String =
                                                            tp.find_string(String::from("wrap"),
                                                                           String::from("repeat"));
                                                        let mut wrap_mode: ImageWrap = ImageWrap::Repeat;
                                                        if wrap == String::from("black") {
                                                            wrap_mode = ImageWrap::Black;
                                                        }
                                                        else if wrap == String::from("clamp") {
                                                            wrap_mode = ImageWrap::Clamp;
                                                        }
                                                        let scale: Float =
                                                            tp.find_float(String::from("scale"), 1.0);
                                                        let mut filename: String =
                                                            tp.find_filename(String::from("filename"),
                                                                             String::new());
                                                        if let Some(ref search_directory) = SEARCH_DIRECTORY {
                                                            // filename = AbsolutePath(ResolveFilename(filename));
                                                            let mut path_buf: PathBuf = PathBuf::from("/");
                                                            path_buf.push(search_directory.as_ref());
                                                            path_buf.push(filename);
                                                            filename = String::from(path_buf.to_str().unwrap());
                                                        }
                                                        // TODO: default depends on:
                                                        // HasExtension(filename,
                                                        // ".tga") ||
                                                        // HasExtension(filename,
                                                        // ".png"));
                                                        let gamma: bool = tp.find_bool(String::from("gamma"), true);
                                                        
                                                        if let Some(mapping) = map {
                                                            let st = Arc::new(ImageTexture::new(mapping,
                                                                                                filename,
                                                                                                do_trilinear,
                                                                                                max_aniso,
                                                                                                wrap_mode,
                                                                                                scale,
                                                                                                gamma));
                                                            graphics_state.spectrum_textures.insert(
                                                                name, st);
                                                        }
                                                    } else if param_set.tex_name == String::from("uv") {
                                                        println!("TODO: CreateUVSpectrumTexture");
                                                    } else if param_set.tex_name == String::from("checkerboard") {
                                                        // CreateCheckerboardSpectrumTexture
                                                        let dim: i32 = tp.find_int(String::from("dimension"), 2);
                                                        if dim != 2 && dim != 3 {
                                                            panic!("{} dimensional checkerboard texture not supported",
                                                                   dim);
                                                        }
                                                        let tex1: Arc<Texture<Spectrum> + Send + Sync> =
                                                            tp.get_spectrum_texture(String::from("tex1"),
                                                                                    Spectrum::new(1.0));
                                                        let tex2: Arc<Texture<Spectrum> + Send + Sync> =
                                                            tp.get_spectrum_texture(String::from("tex2"),
                                                                                    Spectrum::new(0.0));
                                                        if dim == 2 {
                                                            let mut map: Option<Box<TextureMapping2D + Send + Sync>> = None;
                                                            let mapping: String =
                                                                tp.find_string(String::from("mapping"), String::from("uv"));
                                                            if mapping == String::from("uv") {
                                                                println!("TODO: UVMapping2D");
                                                            } else if mapping == String::from("spherical") {
                                                                println!("TODO: SphericalMapping2D");
                                                            } else if mapping == String::from("cylindrical") {
                                                                println!("TODO: CylindricalMapping2D");
                                                            } else if mapping == String::from("planar") {
                                                                map = Some(Box::new(PlanarMapping2D {
                                                                    vs: tp.find_vector3f(String::from("v1"),
                                                                                         Vector3f {
                                                                                             x: 1.0,
                                                                                             y: 0.0,
                                                                                             z: 0.0
                                                                                         }),
                                                                    vt: tp.find_vector3f(String::from("v2"),
                                                                                         Vector3f {
                                                                                             x: 0.0,
                                                                                             y: 1.0,
                                                                                             z: 0.0
                                                                                         }),
                                                                    ds: tp.find_float(String::from("udelta"),
                                                                                      0.0),
                                                                    dt: tp.find_float(String::from("vdelta"),
                                                                                      0.0),
                                                                }));
                                                            } else {
                                                                panic!("2D texture mapping \"{}\" unknown",
                                                                       mapping);
                                                            }
                                                            // TODO: aamode
                                                            if let Some(mapping) = map {
                                                                Arc::new(Checkerboard2DTexture::new(mapping,
                                                                                                    tex1,
                                                                                                    tex2));
                                                            }
                                                        } else { // dim == 3
                                                            println!("TODO: TextureMapping3D");
                                                        }
                                                    } else if param_set.tex_name == String::from("dots") {
                                                        println!("TODO: CreateDotsSpectrumTexture");
                                                    } else if param_set.tex_name == String::from("fbm") {
                                                        println!("TODO: CreateFBmSpectrumTexture");
                                                    } else if param_set.tex_name == String::from("wrinkled") {
                                                        println!("TODO: CreateWrinkledSpectrumTexture");
                                                    } else if param_set.tex_name == String::from("marble") {
                                                        println!("TODO: CreateMarbleSpectrumTexture");
                                                    } else if param_set.tex_name == String::from("windy") {
                                                        println!("TODO: CreateWindySpectrumTexture");
                                                    } else {
                                                        println!("Spectrum texture \"{}\" unknown.",
                                                        param_set.tex_name);
                                                    }
                                                } else {
                                                    panic!("Texture type \"{}\" unknown.", param_set.tex_type);
                                                }
                                        }
                                    }
                                }
                                // MakeFloatTexture(texname, curTransform[0], tp);
                                // or
                                // MakeSpectrumTexture(texname, curTransform[0], tp);
                            } else if param_set.key_word == String::from("Material") {
                                if let Some(ref mut ro) = RENDER_OPTIONS {
                                    println!("Material \"{}\" ", param_set.name);
                                    print_params(&param_set);
                                    // pbrtMaterial (api.cpp:1082)
                                    unsafe {
                                        if let Some(ref mut graphics_state) = GRAPHICS_STATE {
                                            graphics_state.material = name;
                                            graphics_state.material_params.copy_from(&param_set);
                                            graphics_state.current_named_material = String::new();
                                        }
                                    }
                                }
                            } else if param_set.key_word == String::from("Shape") {
                                if let Some(ref mut ro) = RENDER_OPTIONS {
                                    println!("Shape \"{}\" ", param_set.name);
                                    print_params(&param_set);
                                    unsafe {
                                        // pbrtShape (api.cpp:1153)
                                        // TODO: if (!curTransform.IsAnimated()) { ... }
                                        // TODO: transformCache.Lookup(curTransform[0], &ObjToWorld, &WorldToObj);
                                        let obj_to_world: Transform = Transform {
                                            m: CUR_TRANSFORM.t[0].m,
                                            m_inv: CUR_TRANSFORM.t[0].m_inv,
                                        };
                                        let world_to_obj: Transform =Transform {
                                            m: CUR_TRANSFORM.t[0].m_inv,
                                            m_inv: CUR_TRANSFORM.t[0].m,
                                        };
                                        // MakeShapes (api.cpp:296)
                                        if param_set.name == String::from("sphere") {
                                            // CreateSphereShape
                                            let radius: Float = param_set.find_one_float(String::from("radius"),
                                                                                         1.0 as Float);
                                            let z_min: Float = param_set.find_one_float(String::from("zmin"),
                                                                                        -radius);
                                            let z_max: Float = param_set.find_one_float(String::from("zmin"),
                                                                                        radius);
                                            let phi_max: Float = param_set.find_one_float(String::from("phimax"),
                                                                                          360.0 as Float);
                                            let sphere = Arc::new(Sphere::new(obj_to_world,
                                                                              world_to_obj,
                                                                              false,
                                                                              false,
                                                                              radius,
                                                                              z_min,
                                                                              z_max,
                                                                              phi_max));
                                            print!("Sphere {{ object_to_world: {:?}, world_to_object: {:?}, ",
                                                   obj_to_world,
                                                   world_to_obj);
                                            println!("radius: {}, z_min: {}, z_max: {}, phi_max: {} }}",
                                                     radius,
                                                     z_min,
                                                     z_max,
                                                     phi_max)
                                        } else if param_set.name == String::from("cylinder") {
                                            println!("TODO: CreateCylinderShape");
                                        } else if param_set.name == String::from("disk") {
                                            println!("TODO: CreateDiskShape");
                                        } else if param_set.name == String::from("cone") {
                                            println!("TODO: CreateConeShape");
                                        } else if param_set.name == String::from("paraboloid") {
                                            println!("TODO: CreateParaboloidShape");
                                        } else if param_set.name == String::from("hyperboloid") {
                                            println!("TODO: CreateHyperboloidShape");
                                        } else if param_set.name == String::from("curve") {
                                            println!("TODO: CreateCurveShape");
                                        } else if param_set.name == String::from("trianglemesh") {
                                            println!("TODO: CreateTriangleMeshShape");
                                        } else if param_set.name == String::from("plymesh") {
                                            println!("TODO: CreatePLYMesh");
                                        } else if param_set.name == String::from("heightfield") {
                                            println!("TODO: CreateHeightfield");
                                        } else if param_set.name == String::from("loopsubdiv") {
                                            println!("TODO: CreateLoopSubdiv");
                                        } else if param_set.name == String::from("nurbs") {
                                            println!("TODO: CreateNURBS");
                                        } else {
                                            panic!("Shape \"{}\" unknown.", param_set.name);
                                        }
                                    }
                                }
                            } else {
                                println!("PARAM_SET: {}", param_set.key_word);
                            }
                            param_set.reset(String::from(""),
                                            String::from(""),
                                            String::from(""),
                                            String::from(""));
                        }
                    }
                    if optional_parameters.rule == Rule::last_statement {
                        world_end();
                    } else { // statement
                        self._statement();
                    }
                } else if optional_parameters.rule == Rule::parameter {
                    self._parameter();
                } else {
                    println!("ERROR: statement or parameter expected, {:?} found ...",
                             optional_parameters);
                }
            }
        }
        _float_param(&self) -> (String, Vec<Float>) {
            // single float
            (&i: ident, _l: lbrack, &n: number, _r: rbrack) => {
                let string: String = String::from_str(i).unwrap();
                let number: Float = f32::from_str(n).unwrap();
                let fvec: Vec<Float> = vec![number];
                (string, fvec)
            },
            // more than one float
            (&i: ident, _l: lbrack, mut list: _list_of_floats(), _r: rbrack) => {
                let string: String = String::from_str(i).unwrap();
                let mut fvec: Vec<Float> = Vec::new();
                if let Some(floats) = list.pop_front() {
                    match floats {
                        FloatNode::Floats(list) => {
                            for element in list.iter() {
                                match *element {
                                    FloatNode::Floats(ref _l) => {},
                                    FloatNode::OneFloat(ref f) => {
                                        fvec.push(*f);
                                    },
                                };
                            }
                        },
                        FloatNode::OneFloat(_f) => {},
                    };
                }
                (string, fvec)
            },
        }
        _list_of_floats(&self) -> LinkedList<FloatNode> {
            (&n: number, mut head: _float(), mut tail: _list_of_floats()) => {
                let number: Float = Float::from_str(n).unwrap();
                head.push_front(FloatNode::OneFloat(number));
                tail.push_front(FloatNode::Floats(head));
                tail
            },
            () => {
                LinkedList::new()
            }
        }
        _float(&self) -> LinkedList<FloatNode> {
            (&head: number, mut tail: _float()) => {
                let number: Float = Float::from_str(head).unwrap();
                tail.push_front(FloatNode::OneFloat(number));
                tail
            },
            () => {
                LinkedList::new()
            }
        }
        _string_param(&self) -> (String, String) {
            (&i: ident, _l: lbrack, _s: string, &f: filename, _r: rbrack) => {
                let string1: String = String::from_str(i).unwrap();
                let string2: String = String::from_str(f).unwrap();
                (string1, string2)
            },
            (&i1: ident, _l: lbrack, _s: string, &i2: ident, _r: rbrack) => {
                let string1: String = String::from_str(i1).unwrap();
                let string2: String = String::from_str(i2).unwrap();
                (string1, string2)
            },
            (&i: ident, _s: string, &f: filename) => {
                let string1: String = String::from_str(i).unwrap();
                let string2: String = String::from_str(f).unwrap();
                (string1, string2)
            },
            (&i1: ident, _s: string, &i2: ident) => {
                let string1: String = String::from_str(i1).unwrap();
                let string2: String = String::from_str(i2).unwrap();
                (string1, string2)
            },
        }
        _integer_param(&self) -> (String, Vec<i32>) {
            // single integer
            (&i: ident, _l: lbrack, &n: integer, _r: rbrack) => {
                let string: String = String::from_str(i).unwrap();
                let number: i32 = i32::from_str(n).unwrap();
                let ivec: Vec<i32> = vec![number];
                (string, ivec)
            },
            // more than one integer
            (&i: ident, _l: lbrack, mut list: _list_of_integers(), _r: rbrack) => {
                let string: String = String::from_str(i).unwrap();
                let mut ivec: Vec<i32> = Vec::new();
                if let Some(integers) = list.pop_front() {
                    match integers {
                        IntNode::Ints(list) => {
                            for element in list.iter() {
                                match *element {
                                    IntNode::Ints(ref _l) => {},
                                    IntNode::OneInt(ref i) => {
                                        ivec.push(*i);
                                    },
                                };
                            }
                        },
                        IntNode::OneInt(_i) => {},
                    };
                }
                (string, ivec)
            },
        }
        _list_of_integers(&self) -> LinkedList<IntNode> {
            (&i: integer, mut head: _integer(), mut tail: _list_of_integers()) => {
                let number: i32 = i32::from_str(i).unwrap();
                head.push_front(IntNode::OneInt(number));
                tail.push_front(IntNode::Ints(head));
                tail
            },
            () => {
                LinkedList::new()
            }
        }
        _integer(&self) -> LinkedList<IntNode> {
            (&head: integer, mut tail: _integer()) => {
                let number: i32 = i32::from_str(head).unwrap();
                tail.push_front(IntNode::OneInt(number));
                tail
            },
            () => {
                LinkedList::new()
            }
        }
        _point_param(&self) -> (String, Vec<Float>) {
            // single point
            (&i: ident, _l: lbrack, &n1: number, &n2: number, &n3: number, _r: rbrack) => {
                let string: String = String::from_str(i).unwrap();
                let number1: Float = f32::from_str(n1).unwrap();
                let number2: Float = f32::from_str(n2).unwrap();
                let number3: Float = f32::from_str(n3).unwrap();
                let fvec: Vec<Float> = vec![number1, number2, number3];
                (string, fvec)
            },
            // more than one point (3 floats define one point)
            (&i: ident, _l: lbrack, mut list: _list_of_floats(), _r: rbrack) => {
                let string: String = String::from_str(i).unwrap();
                let mut fvec: Vec<Float> = Vec::new();
                if let Some(floats) = list.pop_front() {
                    match floats {
                        FloatNode::Floats(list) => {
                            for element in list.iter() {
                                match *element {
                                    FloatNode::Floats(ref _l) => {},
                                    FloatNode::OneFloat(ref f) => {
                                        fvec.push(*f);
                                    },
                                };
                            }
                        },
                        FloatNode::OneFloat(_f) => {},
                    };
                }
                assert!(fvec.len() % 3 == 0, "point parameters need 3 coordinates");
                (string, fvec)
            },
        }
        _vector_param(&self) -> (String, Float, Float, Float) {
            (&i: ident, _l: lbrack, &n1: number, &n2: number, &n3: number, _r: rbrack) => {
                let string: String = String::from_str(i).unwrap();
                let number1: Float = f32::from_str(n1).unwrap();
                let number2: Float = f32::from_str(n2).unwrap();
                let number3: Float = f32::from_str(n3).unwrap();
                (string, number1, number2, number3)
            },
        }
        _rgb_param(&self) -> (String, Float, Float, Float) {
            (&i: ident, _l: lbrack, &n1: number, &n2: number, &n3: number, _r: rbrack) => {
                let string: String = String::from_str(i).unwrap();
                let number1: Float = f32::from_str(n1).unwrap();
                let number2: Float = f32::from_str(n2).unwrap();
                let number3: Float = f32::from_str(n3).unwrap();
                (string, number1, number2, number3)
            },
        }
        _spectrum_param(&self) -> String {
            (&s: string, _i: ident) => {
                let string: String = String::from_str(s).unwrap();
                string
            },
        }
        _texture_param(&self) -> (String, String) {
            (&i1: ident, _s: string, &i2: ident) => {
                let string1: String = String::from_str(i1).unwrap();
                let string2: String = String::from_str(i2).unwrap();
                (string1, string2)
            },
        }
        // keywords
        _keyword(&self) -> () {
            (_ab: attribute_begin) => {
                println!("AttributeBegin");
                unsafe {
                    if let Some(ref mut graphics_state) = GRAPHICS_STATE {
                        if let Some(ref mut pushed_graphics_states) = PUSHED_GRAPHICS_STATES {
                            let mut param_set: ParamSet = ParamSet::default();
                            param_set.copy_from(&graphics_state.material_params);
                            pushed_graphics_states.push(GraphicsState {
                                float_textures: graphics_state.float_textures.clone(),
                                spectrum_textures: graphics_state.spectrum_textures.clone(),
                                material_params: param_set,
                                material: String::from(graphics_state.material.as_ref()),
                                current_named_material: String::from(graphics_state.current_named_material.as_ref()),
                            });
                        }
                        if let Some(ref mut pt) = PUSHED_TRANSFORMS {
                            pt.push(TransformSet {
                                t: [
                                    Transform {
                                        m: CUR_TRANSFORM.t[0].m,
                                        m_inv: CUR_TRANSFORM.t[0].m_inv,},
                                    Transform {
                                        m: CUR_TRANSFORM.t[1].m,
                                        m_inv: CUR_TRANSFORM.t[1].m_inv,},
                                    ]
                            });
                        }
                        // TODO? pushedActiveTransformBits.push_back(activeTransformBits);
                    }
                }
                self._pbrt();
            },
            (_ae: attribute_end) => {
                println!("AttributeEnd");
                unsafe {
                    if let Some(ref mut graphics_state) = GRAPHICS_STATE {
                        if let Some(ref mut pushed_graphics_states) = PUSHED_GRAPHICS_STATES {
                            if !(pushed_graphics_states.len() >= 1_usize) {
                                panic!("Unmatched pbrtAttributeEnd() encountered.")
                            }
                            let pgs: GraphicsState = pushed_graphics_states.pop().unwrap();
                            // material_params
                            graphics_state.material_params.reset(String::new(),
                                                                 String::from(""),
                                                                 String::from(""),
                                                                 String::new());
                            graphics_state.material_params.copy_from(&pgs.material_params);
                            // material
                            graphics_state.material = String::from(pgs.material.as_ref());
                        }
                        if let Some(ref mut pt) = PUSHED_TRANSFORMS {
                            let popped_transform_set: TransformSet = pt.pop().unwrap();
                            CUR_TRANSFORM.t[0] = popped_transform_set.t[0];
                            CUR_TRANSFORM.t[1] = popped_transform_set.t[1];
                            // println!("CUR_TRANSFORM: {:?}", CUR_TRANSFORM);
                        }
                        // TODO? pushedActiveTransformBits.push_back(activeTransformBits);
                    }
                }
                self._pbrt();
            },
            (_wb: world_begin) => {
                println!("WorldBegin");
                unsafe {
                    CUR_TRANSFORM.t[0] = Transform::default();
                    CUR_TRANSFORM.t[1] = Transform::default();
                    // println!("CUR_TRANSFORM: {:?}", CUR_TRANSFORM);
                    if let Some(ref mut named_coordinate_systems) = NAMED_COORDINATE_SYSTEMS {
                        named_coordinate_systems.insert("world",
                                                        TransformSet {
                                                            t: [Transform::default(); 2]
                                                        });
                        // println!("NAMED_COORDINATE_SYSTEMS: {:?}", named_coordinate_systems);
                    }
                }
                self._pbrt();
            },
            (t) => { println!("TODO: {:?}", t); },
        }
        // numbers
        _number(&self) -> Float {
            (&n: number) => {
                let number: Float = f32::from_str(n).unwrap();
                number
            },
        }
        // strings
        _string(&self) -> String {
            (_s: string) => {
                self._ident()
            },
        }
        // identifiers
        _ident(&self) -> String {
            (&i: ident) => {
                let string: String = String::from_str(i).unwrap();
                string
            },
        }
    }
}

fn print_usage(program: &str, opts: Options) {
    let brief = format!("Usage: {} [options]", program);
    print!("{}", opts.usage(&brief));
}

fn print_version(program: &str) {
    println!("{} {}", program, VERSION);
}

fn print_params(params: &ParamSet) {
    for p in &params.ints {
        if p.n_values == 1_usize {
            println!("  \"integer {}\" [{}]", p.name, p.values[0]);
        } else {
            print!("  \"integer {}\" [ ", p.name);
            for i in 0..p.n_values {
                print!("{} ", p.values[i]);
            }
            println!("]");
        }
    }
    for p in &params.floats {
        if p.n_values == 1_usize {
            println!("  \"float {}\" [{}]", p.name, p.values[0]);
        } else {
            print!("  \"float {}\" [ ", p.name);
            for i in 0..p.n_values {
                print!("{} ", p.values[i]);
            }
            println!("]");
        }
    }
    for p in &params.point3fs {
        if p.n_values == 1_usize {
            println!("  \"point {}\" [{} {} {}]",
                     p.name,
                     p.values[0].x,
                     p.values[0].y,
                     p.values[0].z);
        } else {
            println!("  \"point {}\" [", p.name);
            for i in 0..p.n_values {
                println!("    {} {} {} ", p.values[i].x, p.values[i].y, p.values[i].z);
            }
            println!("  ]");
        }
    }
    for p in &params.vector3fs {
        if p.n_values == 1_usize {
            println!("  \"vector {}\" [{} {} {}]",
                     p.name,
                     p.values[0].x,
                     p.values[0].y,
                     p.values[0].z);
        }
    }
    for p in &params.spectra {
        if p.n_values == 1_usize {
            println!("  \"rgb {}\" [{} {} {}]",
                     p.name,
                     p.values[0].c[0],
                     p.values[0].c[1],
                     p.values[0].c[2]);
        }
    }
    for p in &params.strings {
        if p.n_values == 1_usize {
            println!("  \"string {}\" [\"{}\"]", p.name, p.values[0]);
        }
    }
    for p in &params.textures {
        if p.n_values == 1_usize {
            println!("  \"texture {}\" \"{}\"", p.name, p.values[0]);
        }
    }
}

fn create_material() -> Arc<Material + Send + Sync> {
    unsafe {
        if let Some(ref mut graphics_state) = GRAPHICS_STATE {
            // CreateMaterial
            let mut material_params = ParamSet::default();
            material_params.copy_from(&graphics_state.material_params);
            let mut mp: TextureParams = TextureParams {
                float_textures: graphics_state.float_textures.clone(),
                spectrum_textures: graphics_state.spectrum_textures.clone(),
                geom_params: ParamSet::default(),
                material_params: material_params,
            };
            if graphics_state.current_named_material != String::new() {
                println!("TODO: CreateMaterial, if (currentNamedMaterial != \"\")");
            } else {
                // MakeMaterial
                assert_ne!(graphics_state.material, String::new());
                assert_ne!(graphics_state.material, String::from("none"));
                if graphics_state.material == String::from("matte") {
                    println!("TODO: CreateMatteMaterial");
                } else if graphics_state.material == String::from("plastic") {
                    println!("TODO: CreatePlasticMaterial");
                } else if graphics_state.material == String::from("translucent") {
                    println!("TODO: CreateTranslucentMaterial");
                } else if graphics_state.material == String::from("glass") {
                    println!("TODO: CreateGlassMaterial");
                } else if graphics_state.material == String::from("mirror") {
                    println!("TODO: CreateMirrorMaterial");
                    let kr = mp.get_spectrum_texture(String::from("Kr"),
                                                     Spectrum::new(0.9 as Float));
                    // TODO: std::shared_ptr<Texture<Float>> bumpMap = mp.GetFloatTextureOrNull("bumpmap");
                    let mirror = Arc::new(MirrorMaterial { kr: kr });
                    return mirror;
                } else if graphics_state.material == String::from("hair") {
                    println!("TODO: CreateHairMaterial");
                } else if graphics_state.material == String::from("mix") {
                    println!("TODO: CreateMixMaterial");
                } else if graphics_state.material == String::from("metal") {
                    println!("TODO: CreateMetalMaterial");
                } else if graphics_state.material == String::from("substrate") {
                    println!("TODO: CreateSubstrateMaterial");
                } else if graphics_state.material == String::from("uber") {
                    println!("TODO: CreateUberMaterial");
                } else if graphics_state.material == String::from("subsurface") {
                    println!("TODO: CreateSubsurfaceMaterial");
                } else if graphics_state.material == String::from("kdsubsurface") {
                    println!("TODO: CreateKdsubsurfaceMaterial");
                } else if graphics_state.material == String::from("fourier") {
                    println!("TODO: CreateFourierMaterial");
                } else {
                    panic!("Material \"{}\" unknown.", graphics_state.material);
                }
            }
        }
    }
    let kd = Arc::new(ConstantTexture::new(Spectrum::new(0.5)));
    Arc::new(MatteMaterial::new(kd, 0.0 as Float))
}

fn world_end() {
    println!("WorldEnd");
    unsafe {
        if let Some(ref mut pushed_graphics_states) = PUSHED_GRAPHICS_STATES {
            assert!(pushed_graphics_states.len() == 0_usize,
                    "Missing end to pbrtAttributeBegin()");
            if let Some(ref mut pt) = PUSHED_TRANSFORMS {
                assert!(pt.len() == 0_usize, "Missing end to pbrtTransformBegin()");
                if let Some(ref mut ro) = RENDER_OPTIONS {
                    // MakeFilter
                    let mut some_filter: Option<Arc<Filter + Sync + Send>> = None;
                    if ro.filter_name == String::from("box") {
                        let xw: Float = ro.filter_params.find_one_float(String::from("xwidth"),
                                                                        0.5);
                        let yw: Float = ro.filter_params.find_one_float(String::from("ywidth"),
                                                                        0.5);
                        let box_filter: Arc<Filter + Sync + Send> =
                            Arc::new(BoxFilter {
                                         radius: Vector2f { x: xw, y: yw },
                                         inv_radius: Vector2f {
                                             x: 1.0 / xw,
                                             y: 1.0 / yw,
                                         },
                                     });
                        some_filter = Some(box_filter);
                    } else if ro.filter_name == String::from("gaussian") {
                        println!("TODO: CreateGaussianFilter");
                    } else if ro.filter_name == String::from("mitchell") {
                        println!("TODO: CreateMitchellFilter");
                    } else if ro.filter_name == String::from("sinc") {
                        println!("TODO: CreateSincFilter");
                    } else if ro.filter_name == String::from("triangle") {
                        println!("TODO: CreateTriangleFilter");
                    } else {
                        panic!("Filter \"{}\" unknown.", ro.filter_name);
                    }
                    // MakeFilm
                    if ro.film_name == String::from("image") {
                        println!("TODO: CreateFilm");
                        let filename: String =
                            ro.film_params
                                .find_one_string(String::from("filename"), String::new());
                        println!("filename = {:?}", filename);
                        let xres: i32 = ro.film_params
                            .find_one_int(String::from("xresolution"), 1280);
                        let yres: i32 = ro.film_params
                            .find_one_int(String::from("yresolution"), 720);
                        // TODO: if (PbrtOptions.quickRender) xres = std::max(1, xres / 4);
                        // TODO: if (PbrtOptions.quickRender) yres = std::max(1, yres / 4);
                        println!("xres = {:?}", xres);
                        println!("yres = {:?}", yres);
                        let crop: Bounds2f = Bounds2f {
                            p_min: Point2f { x: 0.0, y: 0.0 },
                            p_max: Point2f { x: 1.0, y: 1.0 },
                        };
                        // TODO: const Float *cr = params.FindFloat("cropwindow", &cwi);
                        let scale: Float = ro.film_params.find_one_float(String::from("scale"),
                                                                         1.0);
                        println!("scale = {:?}", scale);
                        let diagonal: Float = ro.film_params
                            .find_one_float(String::from("diagonal"), 35.0);
                        println!("diagonal = {:?}", diagonal);
                        let max_sample_luminance: Float =
                            ro.film_params
                                .find_one_float(String::from("maxsampleluminance"),
                                                std::f32::INFINITY);
                        if let Some(filter) = some_filter {
                            let film: Film = Film::new(Point2i { x: xres, y: yres },
                                                       crop,
                                                       filter,
                                                       diagonal,
                                                       filename,
                                                       scale,
                                                       max_sample_luminance);
                            // MakeCamera
                            // TODO: MediumInterface mediumInterface = graphicsState.CreateMediumInterface();
                            println!("camera_to_world: {:?}", ro.camera_to_world);
                            let animated_cam_to_world: AnimatedTransform =
                                AnimatedTransform::new(&ro.camera_to_world.t[0],
                                                       ro.transform_start_time,
                                                       &ro.camera_to_world.t[1],
                                                       ro.transform_end_time);
                            if ro.camera_name == String::from("perspective") {
                                let shutteropen: Float =
                                    ro.camera_params
                                        .find_one_float(String::from("shutteropen"), 0.0);
                                let shutterclose: Float =
                                    ro.camera_params
                                        .find_one_float(String::from("shutterclose"), 1.0);
                                // TODO: std::swap(shutterclose, shutteropen);
                                assert!(shutterclose >= shutteropen);
                                let lensradius: Float =
                                    ro.camera_params
                                        .find_one_float(String::from("lensradius"), 0.0);
                                let focaldistance: Float =
                                    ro.camera_params
                                        .find_one_float(String::from("focaldistance"), 1e6);
                                let frame: Float =
                                    ro.camera_params
                                        .find_one_float(String::from("frameaspectratio"),
                                                        (film.full_resolution.x as Float) /
                                                        (film.full_resolution.y as Float));
                                let mut screen: Bounds2f = Bounds2f::default();
                                if frame > 1.0 {
                                    screen.p_min.x = -frame;
                                    screen.p_max.x = frame;
                                    screen.p_min.y = -1.0;
                                    screen.p_max.y = 1.0;
                                } else {
                                    screen.p_min.x = -1.0;
                                    screen.p_max.x = 1.0;
                                    screen.p_min.y = -1.0 / frame;
                                    screen.p_max.y = 1.0 / frame;
                                }
                                // TODO: const Float *sw = params.FindFloat("screenwindow", &swi);
                                let fov: Float =
                                    ro.camera_params.find_one_float(String::from("fov"), 90.0);
                                // let halffov: Float =
                                //     ro.camera_params.find_one_float(String::from("halffov"), -1.0);
                                // TODO: if (halffov > 0.f)
                                let perspective_camera: PerspectiveCamera =
                                    PerspectiveCamera::new(animated_cam_to_world,
                                                           screen,
                                                           shutteropen,
                                                           shutterclose,
                                                           lensradius,
                                                           focaldistance,
                                                           fov,
                                                           film);
                            } else if ro.camera_name == String::from("orthographic") {
                                println!("TODO: CreateOrthographicCamera");
                            } else if ro.camera_name == String::from("realistic") {
                                println!("TODO: CreateRealisticCamera");
                            } else if ro.camera_name == String::from("environment") {
                                println!("TODO: CreateEnvironmentCamera");
                            } else {
                                panic!("Camera \"{}\" unknown.", ro.camera_name);
                            }
                            // MakeSampler
                            let mut some_sampler: Option<Arc<Sampler + Sync + Send>> = None;
                            if ro.sampler_name == String::from("lowdiscrepancy") ||
                                ro.sampler_name == String::from("02sequence")
                            {
                                let nsamp: i32 =
                                    ro.sampler_params
                                    .find_one_int(String::from("pixelsamples"), 16);
                                let sd: i32 = ro.sampler_params
                                    .find_one_int(String::from("dimensions"), 4);
                                // TODO: if (PbrtOptions.quickRender) nsamp = 1;
                                let sampler =
                                    Arc::new(ZeroTwoSequenceSampler::new(nsamp as i64,
                                                                         sd as i64));
                                some_sampler = Some(sampler);
                            } else if ro.sampler_name == String::from("maxmindist") {
                                println!("TODO: CreateMaxMinDistSampler");
                            } else if ro.sampler_name == String::from("halton") {
                                println!("TODO: CreateHaltonSampler");
                            } else if ro.sampler_name == String::from("sobol") {
                                println!("TODO: CreateSobolSampler");
                            } else if ro.sampler_name == String::from("random") {
                                println!("TODO: CreateRandomSampler");
                            } else if ro.sampler_name == String::from("stratified") {
                                println!("TODO: CreateStratifiedSampler");
                            } else {
                                panic!("Sampler \"{}\" unknown.", ro.sampler_name);
                            }
                            // MakeIntegrator
                            // MakeScene
                        } else {
                            panic!("Unable to create film.");
                        }
                    } else {
                        panic!("Film \"{}\" unknown.", ro.film_name);
                    }
                }
            }
        }
    }
}

fn main() {
    // handle command line options
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();
    let mut opts = Options::new();
    opts.optflag("h", "help", "print this help menu");
    opts.optopt("i", "", "parse an input file", "FILE");
    opts.optflag("v", "version", "print version number");
    let matches = match opts.parse(&args[1..]) {
        Ok(m) => m,
        Err(f) => panic!(f.to_string()),
    };
    if matches.opt_present("h") {
        print_usage(&program, opts);
        return;
    } else if matches.opt_present("i") {
        let infile = matches.opt_str("i");
        match infile {
            Some(x) => {
                println!("FILE = {}", x);
                let f = File::open(x.clone()).unwrap();
                let ip: &Path = Path::new(x.as_str());
                if ip.is_relative() {
                    let cp: PathBuf = env::current_dir().unwrap();
                    let pb: PathBuf = cp.join(ip);
                    let search_directory: &Path = pb.as_path().parent().unwrap();
                    println!("search_directory is {}", search_directory.display());
                    unsafe {
                        SEARCH_DIRECTORY = Some(Box::new(PathBuf::from(search_directory)));
                    }
                }
                let mut reader = BufReader::new(f);
                let mut str_buf: String = String::default();
                reader.read_to_string(&mut str_buf);
                unsafe {
                    // render options
                    NAMED_COORDINATE_SYSTEMS = Some(Box::new(HashMap::new()));
                    RENDER_OPTIONS = Some(Box::new(RenderOptions::default()));
                    GRAPHICS_STATE = Some(Box::new(GraphicsState::default()));
                    PUSHED_GRAPHICS_STATES = Some(Box::new(Vec::new()));
                    PUSHED_TRANSFORMS = Some(Box::new(Vec::new()));
                    PARAM_SET = Some(Box::new(ParamSet::default()));
                    // parser
                    let mut parser = Rdp::new(StringInput::new(&str_buf));
                    assert!(parser.pbrt());
                    assert!(parser.end());
                    println!("{:?}", parser.queue());
                    println!("do something with created tokens ...");
                    parser.main();
                    println!("done.");
                }
            }
            None => panic!("no input file name"),
        }
        return;
    } else if matches.opt_present("v") {
        print_version(&program);
        return;
    } else {
        print_usage(&program, opts);
        return;
    }
}
