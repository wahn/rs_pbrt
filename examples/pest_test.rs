#![recursion_limit="2000"]
#![feature(drop_types_in_const)]

extern crate pest;
#[macro_use]
extern crate pest_derive;
extern crate getopts;
extern crate pbrt;

use pbrt::{AnimatedTransform, AOIntegrator, Bounds2f, Bounds2i, BoxFilter, BVHAccel,
           Checkerboard2DTexture, ConstantTexture, Cylinder, DiffuseAreaLight,
           DirectLightingIntegrator, Disk, DistantLight, Film, Filter, Float, GaussianFilter,
           GeometricPrimitive, GlassMaterial, GraphicsState, ImageTexture, ImageWrap, Light,
           LightStrategy, Material, MatteMaterial, Matrix4x4, MirrorMaterial, Normal3f, ParamSet,
           PathIntegrator, PerspectiveCamera, PlanarMapping2D, PlasticMaterial, Point2f, Point2i,
           Point3f, PointLight, RenderOptions, SamplerIntegrator, Scene, Shape, Spectrum, Sphere,
           SplitMethod, Texture, TextureMapping2D, TextureParams, Transform, TransformSet,
           Triangle, TriangleMesh, UVMapping2D, Vector2f, Vector3f, ZeroTwoSequenceSampler};
// parser
use pest::Parser;
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
use std::thread;

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

#[derive(Parser)]
#[grammar = "../examples/pbrt.pest"]
struct PbrtParser;

// impl_rdop! {
//     grammar! {
//         pbrt = _{ whitespace? ~ (statement | comment)* ~ last_statement }
//         statement = { look_at | translate | rotate | scale | transform | concat_transform | named_statement | keyword }
//         named_statement = { accelerator |
//                             camera |
//                             pixel_filter |
//                             sampler |
//                             film |
//                             integrator |
//                             coord_sys_transform |
//                             area_light_source |
//                             light_source |
//                             texture |
//                             material |
//                             make_named_material |
//                             named_material |
//                             shape }
//         parameter = { bool_param |
//                       float_param |
//                       string_param |
//                       integer_param |
//                       point_param |
//                       vector_param |
//                       normal_param |
//                       rgb_param |
//                       spectrum_param |
//                       texture_param }
//         bool_param = { (["\"bool"] ~ ident ~ ["\""] ~ lbrack ~ string ~ rbrack) }
//         float_param = { (["\"float"] ~ ident ~ ["\""] ~ lbrack ~ number+ ~ rbrack) |
//                         (["\"float"] ~ ident ~ ["\""] ~ number) }
//         string_param = { (["\"string"] ~ ident ~ ["\""] ~ lbrack ~ string ~ rbrack) |
//                          (["\"string"] ~ ident ~ ["\""] ~ string) }
//         integer_param = { ["\"integer"] ~ ident ~ ["\""] ~ lbrack ~ integer+ ~ rbrack }
//         point_param = { ["\"point"] ~ ident ~ ["\""] ~ lbrack ~ number+ ~ rbrack }
//         vector_param = { ["\"vector"] ~ ident ~ ["\""] ~ lbrack ~ number ~ number ~ number ~ rbrack }
//         normal_param = { ["\"normal"] ~ ident ~ ["\""] ~ lbrack ~ number+ ~ rbrack }
//         rgb_param = { (["\"rgb"] ~ ident ~ ["\""] ~ lbrack ~ number ~ number ~ number ~ rbrack) |
//                       (["\"color"] ~ ident ~ ["\""] ~ lbrack ~ number ~ number ~ number ~ rbrack) }
//         spectrum_param = { ["\"spectrum\""] ~ string }
//         texture_param = { ["\"texture"] ~ ident ~ ["\""] ~ string }
//         // Translate x y z
//         translate = { ["Translate"] ~
//                    // followed by 3 numbers:
//                    number ~ number ~ number
//         }
//         // Rotate angle x y z
//         rotate = { ["Rotate"] ~
//                    // followed by 4 numbers:
//                    number ~ number ~ number ~ number
//         }
//         // Scale x y z
//         scale = { ["Scale"] ~
//                    // followed by 3 numbers:
//                    number ~ number ~ number
//         }
//         // Transform m00 .. m33
//         transform = { (["Transform"] ~ lbrack ~
//                        // followed by 16 numbers:
//                        number ~ number ~ number ~ number ~
//                        number ~ number ~ number ~ number ~
//                        number ~ number ~ number ~ number ~
//                        number ~ number ~ number ~ number ~ rbrack) |
//                       (["Transform"] ~
//                        // followed by 16 numbers:
//                        number ~ number ~ number ~ number ~
//                        number ~ number ~ number ~ number ~
//                        number ~ number ~ number ~ number ~
//                        number ~ number ~ number ~ number)
//         }
//         // ConcatTransform m00 .. m33
//         concat_transform = { (["ConcatTransform"] ~ lbrack ~
//                               // followed by 16 numbers:
//                               number ~ number ~ number ~ number ~
//                               number ~ number ~ number ~ number ~
//                               number ~ number ~ number ~ number ~
//                               number ~ number ~ number ~ number ~ rbrack) |
//                              (["ConcatTransform"] ~
//                               // followed by 16 numbers:
//                               number ~ number ~ number ~ number ~
//                               number ~ number ~ number ~ number ~
//                               number ~ number ~ number ~ number ~
//                               number ~ number ~ number ~ number)
//         }
//         // LookAt eye_x eye_y eye_z look_x look_y look_z up_x up_y up_z
//         look_at = { ["LookAt"] ~
//                     // followed by 9 numbers:

//                     // eye_x eye_y eye_z
//                     number ~ number ~ number ~
//                     // look_x look_y look_z
//                     number ~ number ~ number ~
//                     // up_x up_y up_z
//                     number ~ number ~ number
//         }
//         // Accelerator "kdtree" "float emptybonus" [0.1]
//         accelerator = { ["Accelerator"] ~ string ~ parameter* }
//         // Camera "perspective" "float fov" [90] ...
//         camera = { ["Camera"] ~ string ~ parameter* }
//         // PixelFilter "mitchell" "float xwidth" [2] "float ywidth" [2]
//         pixel_filter = { ["PixelFilter"] ~ string ~ parameter* }
//         // Sampler "halton"
//         sampler = { ["Sampler"] ~ string ~ parameter* }
//         // Film "image" "string filename" ["..."] ...
//         film = { ["Film"] ~ string ~ parameter* }
//         // Integrator "path" "integer maxdepth" [5]
//         integrator = { ["Integrator"] ~ string ~ parameter* }
//         // CoordSysTransform "camera"
//         coord_sys_transform = { ["CoordSysTransform"] ~ string }
//         // AreaLightSource "diffuse" "rgb L" [ .5 .5 .5 ]
//         area_light_source = { ["AreaLightSource"] ~ string ~ parameter* }
//         // LightSource "point" "rgb I" [ .5 .5 .5 ]
//         light_source = { ["LightSource"] ~ string ~ parameter* }
//         // Texture "mydiffuse" "spectrum" "imagemap" "string filename" "image.tga"
//         texture = { ["Texture"] ~ string ~ string ~ string ~ parameter* }
//         // Material "matte" "texture Kd" "mydiffuse"
//         material = { ["Material"] ~ string ~ parameter* }
//         // MakeNamedMaterial "myplastic" "string type" "plastic" "float roughness" [0.1]
//         make_named_material = { ["MakeNamedMaterial"] ~ string ~ parameter* }
//         // NamedMaterial "myplastic"
//         named_material = { ["NamedMaterial"] ~ string ~ parameter* }
//         // Shape "sphere" "float radius" [0.25]
//         shape = { ["Shape"] ~ string ~ parameter* }
//         // keywords
//         keyword = {
//             (["ActiveTransform"] |
//              ["All"] |
//              attribute_begin |
//              attribute_end |
//              ["CoordinateSystem"] |
//              ["EndTime"] |
//              ["Identity"] |
//              ["Include"] |
//              ["MakeNamedMedium"] |
//              ["MediumInterface"] |
//              ["ObjectBegin"] |
//              ["ObjectEnd"] |
//              ["ObjectInstance"] |
//              ["ReverseOrientation"] |
//              ["StartTime"] |
//              ["TransformBegin"] |
//              ["TransformEnd"] |
//              ["TransformTimes"] |
//              world_begin
//             )
//         }
//         attribute_begin = { ["AttributeBegin"] }
//         attribute_end = { ["AttributeEnd"] }
//         world_begin = { ["WorldBegin"] }
//         // IDENT [a-zA-Z_][a-zA-Z_0-9]*
//         ident =  { (['a'..'z'] | ['A'..'Z'] | ["_"]) ~
//                    (['a'..'z'] | ['A'..'Z'] | ["_"] | ["-"] | ['0'..'9'])* }
//         string = { (["\""] ~ ident ~ ["\""]) | (["\""] ~ filename ~ ["\""]) }
//         filename = { (['a'..'z'] | ['A'..'Z'] | ["_"]) ~ // TODO: can be a full path
//                      (['a'..'z'] | ['A'..'Z'] | ["_"] | ["-"] | ["."] | ["/"] | ['0'..'9'])* }
//         // "[" { return LBRACK; }
//         lbrack = { ["["] }
//         // "]" { return RBRACK; }
//         rbrack = { ["]"] }
//         // NUMBER [-+]?([0-9]+|(([0-9]+\.[0-9]*)|(\.[0-9]+)))([eE][-+]?[0-9]+)?
//         number = @{
//             (["-"] | ["+"])? ~ // optional sign, followed by
//             (
//                 (
//                     (["."] ~ ['0'..'9']+) // dot and digits
//                         | // or
//                     (['0'..'9']+ ~ ["."] ~ ['0'..'9']*) // digits, dot, and (optional digits)
//                 )
//                     | // or
//                 ['0'..'9']+ // just digits
//             ) ~ ( // followed by (optional)
//                 (["e"] | ["E"]) ~ // 'e' or 'E', followed by
//                 (["-"] | ["+"])? ~ // optional sign, followed by
//                 ['0'..'9']+ // digits
//             )?
//         }
//         integer = @{
//             (["-"] | ["+"])? ~ // optional sign, followed by
//                 (
//                     ['1'..'9'] ~ // at least one non-zero digit, followed by
//                     ['0'..'9']* // just digits
//                 )
//                     | // or
//                 ['0'..'9'] // single digit
//         }
//         last_statement = @{ whitespace? ~ ["WorldEnd"] ~ (whitespace | comment)* }
//         whitespace = _{ ([" "] | ["\t"] | ["\r"] | ["\n"]) }
//         comment = _{ ( ["#"] ~ (!(["\r"] | ["\n"]) ~ any)* ~ (["\n"] | ["\r\n"] | ["\r"] | eoi) ) }
//     }
//     process! {
//         main(&self) -> () {
//             (_list: _pbrt()) => {
//             }
//         }
//         _pbrt(&self) -> () {
//             (_head: statement, _tail: _statement()) => {},
//             (_l: last_statement) => { pbrt_world_end(); },
//         }
//         // statements
//         _statement(&self) -> () {
//             (_head: look_at, _tail: _look_at()) => {},
//             (_head: translate, _tail: _translate()) => {},
//             (_head: rotate, _tail: _rotate()) => {},
//             (_head: scale, _tail: _scale()) => {},
//             (_head: transform, _tail: _transform()) => {},
//             (_head: concat_transform, _tail: _concat_transform()) => {},
//             (_head: named_statement, _tail: _named_statement()) => {},
//             (_head: keyword, _tail: _keyword()) => {},
//         }
//         _look_at(&self) -> () {
//             (eye_x: _number(), eye_y: _number(), eye_z: _number(),
//              look_x: _number(), look_y: _number(), look_z: _number(),
//              up_x: _number(), up_y: _number(), up_z: _number()) => {
//                 // println!("LookAt {} {} {} {} {} {} {} {} {}",
//                 //          eye_x, eye_y, eye_z,
//                 //          look_x, look_y, look_z,
//                 //          up_x, up_y, up_z,);
//                 let pos: Point3f = Point3f { x: eye_x, y: eye_y, z: eye_z, };
//                 let look: Point3f = Point3f { x: look_x, y: look_y, z: look_z, };
//                 let up: Vector3f = Vector3f { x: up_x, y: up_y, z: up_z, };
//                 let look_at: Transform = Transform::look_at(pos, look, up);
//                 unsafe {
//                     CUR_TRANSFORM.t[0] = CUR_TRANSFORM.t[0] * look_at;
//                     CUR_TRANSFORM.t[1] = CUR_TRANSFORM.t[1] * look_at;
//                     // println!("CUR_TRANSFORM: {:?}", CUR_TRANSFORM);
//                 }
//                 self._pbrt();
//             }
//         }
//         _translate(&self) -> () {
//             (x: _number(), y: _number(), z: _number()) => {
//                 // println!("Translate {} {} {}", x, y, z);
//                 let translate: Transform = Transform::translate(Vector3f { x: x, y: y, z: z, });
//                 unsafe {
//                     CUR_TRANSFORM.t[0] = CUR_TRANSFORM.t[0] * translate;
//                     CUR_TRANSFORM.t[1] = CUR_TRANSFORM.t[1] * translate;
//                     // println!("CUR_TRANSFORM: {:?}", CUR_TRANSFORM);
//                 }
//                 self._pbrt();
//             }
//         }
//         _rotate(&self) -> () {
//             (angle: _number(), x: _number(), y: _number(), z: _number()) => {
//                 // println!("Rotate {} {} {} {}",
//                 //          angle, x, y, z);
//                 let rotate: Transform = Transform::rotate(angle, Vector3f { x: x, y: y, z: z, });
//                 unsafe {
//                     CUR_TRANSFORM.t[0] = CUR_TRANSFORM.t[0] * rotate;
//                     CUR_TRANSFORM.t[1] = CUR_TRANSFORM.t[1] * rotate;
//                     // println!("CUR_TRANSFORM: {:?}", CUR_TRANSFORM);
//                 }
//                 self._pbrt();
//             }
//         }
//         _scale(&self) -> () {
//             (x: _number(), y: _number(), z: _number()) => {
//                 // println!("Scale {} {} {}",
//                 //          x, y, z);
//                 let scale: Transform = Transform::scale(x, y, z);
//                 unsafe {
//                     CUR_TRANSFORM.t[0] = CUR_TRANSFORM.t[0] * scale;
//                     CUR_TRANSFORM.t[1] = CUR_TRANSFORM.t[1] * scale;
//                     // println!("CUR_TRANSFORM: {:?}", CUR_TRANSFORM);
//                 }
//                 self._pbrt();
//             }
//         }
//         _transform(&self) -> () {
//             (_l: lbrack,
//              m00: _number(), m01: _number(), m02: _number(), m03: _number(),
//              m10: _number(), m11: _number(), m12: _number(), m13: _number(),
//              m20: _number(), m21: _number(), m22: _number(), m23: _number(),
//              m30: _number(), m31: _number(), m32: _number(), m33: _number(),
//              _r: rbrack) => {
//                 // println!("Transform [");
//                 // println!("  {} {} {} {}", m00, m01, m02, m03);
//                 // println!("  {} {} {} {}", m10, m11, m12, m13);
//                 // println!("  {} {} {} {}", m20, m21, m22, m23);
//                 // println!("  {} {} {} {}", m30, m31, m32, m33);
//                 // println!("]");
//                 // INFO: The order in PBRT file is different !!!
//                 let transform: Transform = Transform::new(m00, m10, m20, m30,
//                                                           m01, m11, m21, m31,
//                                                           m02, m12, m22, m32,
//                                                           m03, m13, m23, m33);
//                 unsafe {
//                     CUR_TRANSFORM.t[0] = CUR_TRANSFORM.t[0] * transform;
//                     CUR_TRANSFORM.t[1] = CUR_TRANSFORM.t[1] * transform;
//                     // println!("CUR_TRANSFORM: {:?}", CUR_TRANSFORM);
//                 }
//                 self._pbrt();
//             },
//             (m00: _number(), m01: _number(), m02: _number(), m03: _number(),
//              m10: _number(), m11: _number(), m12: _number(), m13: _number(),
//              m20: _number(), m21: _number(), m22: _number(), m23: _number(),
//              m30: _number(), m31: _number(), m32: _number(), m33: _number()) => {
//                 // println!("Transform [");
//                 // println!("  {} {} {} {}", m00, m01, m02, m03);
//                 // println!("  {} {} {} {}", m10, m11, m12, m13);
//                 // println!("  {} {} {} {}", m20, m21, m22, m23);
//                 // println!("  {} {} {} {}", m30, m31, m32, m33);
//                 // println!("]");
//                 // INFO: The order in PBRT file is different !!!
//                 let transform: Transform = Transform::new(m00, m10, m20, m30,
//                                                           m01, m11, m21, m31,
//                                                           m02, m12, m22, m32,
//                                                           m03, m13, m23, m33);
//                 unsafe {
//                     CUR_TRANSFORM.t[0] = CUR_TRANSFORM.t[0] * transform;
//                     CUR_TRANSFORM.t[1] = CUR_TRANSFORM.t[1] * transform;
//                     // println!("CUR_TRANSFORM: {:?}", CUR_TRANSFORM);
//                 }
//                 self._pbrt();
//             },
//         }
//         _concat_transform(&self) -> () {
//             (_l: lbrack,
//              m00: _number(), m01: _number(), m02: _number(), m03: _number(),
//              m10: _number(), m11: _number(), m12: _number(), m13: _number(),
//              m20: _number(), m21: _number(), m22: _number(), m23: _number(),
//              m30: _number(), m31: _number(), m32: _number(), m33: _number(),
//              _r: rbrack) => {
//                 // println!("ConcatTransform [");
//                 // println!("  {} {} {} {}", m00, m01, m02, m03);
//                 // println!("  {} {} {} {}", m10, m11, m12, m13);
//                 // println!("  {} {} {} {}", m20, m21, m22, m23);
//                 // println!("  {} {} {} {}", m30, m31, m32, m33);
//                 // println!("]");
//                 // INFO: The order in PBRT file is different !!!
//                 let transform: Transform = Transform::new(m00, m10, m20, m30,
//                                                           m01, m11, m21, m31,
//                                                           m02, m12, m22, m32,
//                                                           m03, m13, m23, m33);
//                 unsafe {
//                     CUR_TRANSFORM.t[0] = CUR_TRANSFORM.t[0] * transform;
//                     CUR_TRANSFORM.t[1] = CUR_TRANSFORM.t[1] * transform;
//                     // println!("CUR_TRANSFORM: {:?}", CUR_TRANSFORM);
//                 }
//                 self._pbrt();
//             },
//             (m00: _number(), m01: _number(), m02: _number(), m03: _number(),
//              m10: _number(), m11: _number(), m12: _number(), m13: _number(),
//              m20: _number(), m21: _number(), m22: _number(), m23: _number(),
//              m30: _number(), m31: _number(), m32: _number(), m33: _number()) => {
//                 // println!("ConcatTransform [");
//                 // println!("  {} {} {} {}", m00, m01, m02, m03);
//                 // println!("  {} {} {} {}", m10, m11, m12, m13);
//                 // println!("  {} {} {} {}", m20, m21, m22, m23);
//                 // println!("  {} {} {} {}", m30, m31, m32, m33);
//                 // println!("]");
//                 // INFO: The order in PBRT file is different !!!
//                 let transform: Transform = Transform::new(m00, m10, m20, m30,
//                                                           m01, m11, m21, m31,
//                                                           m02, m12, m22, m32,
//                                                           m03, m13, m23, m33);
//                 unsafe {
//                     CUR_TRANSFORM.t[0] = CUR_TRANSFORM.t[0] * transform;
//                     CUR_TRANSFORM.t[1] = CUR_TRANSFORM.t[1] * transform;
//                     // println!("CUR_TRANSFORM: {:?}", CUR_TRANSFORM);
//                 }
//                 self._pbrt();
//             },
//         }
//         // named statements
//         _named_statement(&self) -> () {
//             (_head: accelerator, _tail: _accelerator()) => {},
//             (_head: camera, _tail: _camera()) => {},
//             (_head: pixel_filter, _tail: _pixel_filter()) => {},
//             (_head: sampler, _tail: _sampler()) => {},
//             (_head: film, _tail: _film()) => {},
//             (_head: integrator, _tail: _integrator()) => {},
//             (_head: coord_sys_transform, _tail: _coord_sys_transform()) => {},
//             (_head: area_light_source, _tail: _area_light_source()) => {},
//             (_head: light_source, _tail: _light_source()) => {},
//             (_head: texture, _tail: _texture()) => {},
//             (_head: material, _tail: _material()) => {},
//             (_head: make_named_material, _tail: _make_named_material()) => {},
//             (_head: named_material, _tail: _named_material()) => {},
//             (_head: shape, _tail: _shape()) => {},
//         }
//         _accelerator(&self) -> () {
//             (name: _string(), optional_parameters) => {
//                 unsafe {
//                     if let Some(ref mut ro) = RENDER_OPTIONS {
//                         ro.accelerator_name = name;
//                         if optional_parameters.rule == Rule::statement ||
//                             optional_parameters.rule == Rule::last_statement {
//                             // println!("Accelerator \"{}\" ", ro.accelerator_name);
//                         }
//                     }
//                     if let Some(ref mut param_set) = PARAM_SET {
//                         param_set.reset(String::from("Accelerator"),
//                                         String::from(""),
//                                         String::from(""),
//                                         String::from(""));
//                     }
//                 }
//                 if optional_parameters.rule == Rule::parameter {
//                     self._parameter();
//                 } else if optional_parameters.rule == Rule::statement {
//                     self._statement();
//                 } else if optional_parameters.rule == Rule::last_statement {
//                     pbrt_world_end();
//                 } else {
//                     println!("ERROR: parameter expected, {:?} found ...", optional_parameters);
//                 }
//             },
//         }
//         _camera(&self) -> () {
//             (name: _string(), optional_parameters) => {
//                 unsafe {
//                     if let Some(ref mut ro) = RENDER_OPTIONS {
//                         ro.camera_name = name;
//                         ro.camera_to_world.t[0] =
//                             Transform::inverse(CUR_TRANSFORM.t[0]);
//                         ro.camera_to_world.t[1] =
//                             Transform::inverse(CUR_TRANSFORM.t[1]);
//                         if let Some(ref mut named_coordinate_systems) = NAMED_COORDINATE_SYSTEMS {
//                             named_coordinate_systems.insert("camera",
//                                                             TransformSet {
//                                                                 t: [ro.camera_to_world.t[0],
//                                                                     ro.camera_to_world.t[1]]
//                                                             });
//                             // println!("NAMED_COORDINATE_SYSTEMS: {:?}", named_coordinate_systems);
//                         }
//                         // println!("ro.camera_to_world: {:?}",
//                         //          ro.camera_to_world);
//                     }
//                     if let Some(ref mut param_set) = PARAM_SET {
//                         param_set.reset(String::from("Camera"),
//                                         String::from(""),
//                                         String::from(""),
//                                         String::from(""));
//                     }
//                 }
//                 if optional_parameters.rule == Rule::parameter {
//                     self._parameter();
//                 } else if optional_parameters.rule == Rule::statement {
//                     self._statement();
//                 } else if optional_parameters.rule == Rule::last_statement {
//                     pbrt_world_end();
//                 } else {
//                     println!("ERROR: parameter expected, {:?} found ...", optional_parameters);
//                 }
//             },
//         }
//         _pixel_filter(&self) -> () {
//             (name: _string(), optional_parameters) => {
//                 unsafe {
//                     if let Some(ref mut ro) = RENDER_OPTIONS {
//                         ro.filter_name = name;
//                         if optional_parameters.rule == Rule::statement ||
//                             optional_parameters.rule == Rule::last_statement {
//                             // println!("PixelFilter \"{}\" ", ro.filter_name);
//                         }
//                     }
//                     if let Some(ref mut param_set) = PARAM_SET {
//                         param_set.reset(String::from("PixelFilter"),
//                                         String::from(""),
//                                         String::from(""),
//                                         String::from(""));
//                     }
//                 }
//                 if optional_parameters.rule == Rule::parameter {
//                     self._parameter();
//                 } else if optional_parameters.rule == Rule::statement {
//                     self._statement();
//                 } else if optional_parameters.rule == Rule::last_statement {
//                     pbrt_world_end();
//                 } else {
//                     println!("ERROR: parameter expected, {:?} found ...", optional_parameters);
//                 }
//             },
//         }
//         _sampler(&self) -> () {
//             (name: _string(), optional_parameters) => {
//                 unsafe {
//                     if let Some(ref mut ro) = RENDER_OPTIONS {
//                         ro.sampler_name = name;
//                         if optional_parameters.rule == Rule::statement ||
//                             optional_parameters.rule == Rule::last_statement {
//                             // println!("Sampler \"{}\" ", ro.sampler_name);
//                         }
//                     }
//                     if let Some(ref mut param_set) = PARAM_SET {
//                         param_set.reset(String::from("Sampler"),
//                                         String::from(""),
//                                         String::from(""),
//                                         String::from(""));
//                     }
//                 }
//                 if optional_parameters.rule == Rule::parameter {
//                     self._parameter();
//                 } else if optional_parameters.rule == Rule::statement {
//                     self._statement();
//                 } else if optional_parameters.rule == Rule::last_statement {
//                     pbrt_world_end();
//                 } else {
//                     println!("ERROR: parameter expected, {:?} found ...", optional_parameters);
//                 }
//             },
//         }
//         _film(&self) -> () {
//             (name: _string(), optional_parameters) => {
//                 unsafe {
//                     if let Some(ref mut ro) = RENDER_OPTIONS {
//                         ro.film_name = name;
//                     }
//                     if let Some(ref mut param_set) = PARAM_SET {
//                         param_set.reset(String::from("Film"),
//                                         String::from(""),
//                                         String::from(""),
//                                         String::from(""));
//                     }
//                 }
//                 if optional_parameters.rule == Rule::parameter {
//                     self._parameter();
//                 } else if optional_parameters.rule == Rule::statement {
//                     self._statement();
//                 } else if optional_parameters.rule == Rule::last_statement {
//                     pbrt_world_end();
//                 } else {
//                     println!("ERROR: parameter expected, {:?} found ...", optional_parameters);
//                 }
//             },
//         }
//         _integrator(&self) -> () {
//             (name: _string(), optional_parameters) => {
//                 unsafe {
//                     if let Some(ref mut ro) = RENDER_OPTIONS {
//                         ro.integrator_name = name;
//                     }
//                     if let Some(ref mut param_set) = PARAM_SET {
//                         param_set.reset(String::from("Integrator"),
//                                         String::from(""),
//                                         String::from(""),
//                                         String::from(""));
//                     }
//                 }
//                 if optional_parameters.rule == Rule::parameter {
//                     self._parameter();
//                 } else if optional_parameters.rule == Rule::statement {
//                     self._statement();
//                 } else if optional_parameters.rule == Rule::last_statement {
//                     pbrt_world_end();
//                 } else {
//                     println!("ERROR: parameter expected, {:?} found ...", optional_parameters);
//                 }
//             },
//         }
//         _coord_sys_transform(&self) -> () {
//             (name: _string()) => {
//                 // println!("CoordSysTransform \"{}\" ", name);
//                 unsafe {
//                     if let Some(ref mut named_coordinate_systems) = NAMED_COORDINATE_SYSTEMS {
//                         match named_coordinate_systems.get(name.as_str()) {
//                             Some(transform_set) => {
//                                 CUR_TRANSFORM.t[0] = transform_set.t[0];
//                                 CUR_TRANSFORM.t[1] = transform_set.t[1];
//                                 // println!("CUR_TRANSFORM: {:?}", CUR_TRANSFORM);
//                             },
//                             None => {
//                                 println!("Couldn't find named coordinate system \"{}\"", name);
//                             },
//                         };
//                     }
//                 }
//                 self._pbrt();
//             },
//         }
//         _area_light_source(&self) -> () {
//             (name: _string(), optional_parameters) => {
//                 unsafe {
//                     if let Some(ref mut param_set) = PARAM_SET {
//                         param_set.reset(String::from("AreaLightSource"),
//                                         String::from(name),
//                                         String::from(""),
//                                         String::from(""));
//                     }
//                 }
//                 if optional_parameters.rule == Rule::parameter {
//                     self._parameter();
//                 } else if optional_parameters.rule == Rule::statement {
//                     self._statement();
//                 } else if optional_parameters.rule == Rule::last_statement {
//                     pbrt_world_end();
//                 } else {
//                     println!("ERROR: parameter expected, {:?} found ...", optional_parameters);
//                 }
//             },
//         }
//         _light_source(&self) -> () {
//             (name: _string(), optional_parameters) => {
//                 unsafe {
//                     if let Some(ref mut param_set) = PARAM_SET {
//                         param_set.reset(String::from("LightSource"),
//                                         String::from(name),
//                                         String::from(""),
//                                         String::from(""));
//                     }
//                 }
//                 if optional_parameters.rule == Rule::parameter {
//                     self._parameter();
//                 } else if optional_parameters.rule == Rule::statement {
//                     self._statement();
//                 } else if optional_parameters.rule == Rule::last_statement {
//                     pbrt_world_end();
//                 } else {
//                     println!("ERROR: parameter expected, {:?} found ...", optional_parameters);
//                 }
//             },
//         }
//         _texture(&self) -> () {
//             (name: _string(), tex_type: _string(), tex_name: _string(), optional_parameters) => {
//                 unsafe {
//                     if optional_parameters.rule == Rule::statement ||
//                         optional_parameters.rule == Rule::last_statement
//                     {
//                         // println!("Texture \"{}\" \"{}\" \"{}\" ",
//                         //          name,
//                         //          tex_type,
//                         //          tex_name);
//                     }
//                     if let Some(ref mut param_set) = PARAM_SET {
//                         param_set.reset(String::from("Texture"),
//                                         String::from(name),
//                                         String::from(tex_type),
//                                         String::from(tex_name));
//                     }
//                 }
//                 if optional_parameters.rule == Rule::parameter {
//                     self._parameter();
//                 } else if optional_parameters.rule == Rule::statement {
//                     self._statement();
//                 } else if optional_parameters.rule == Rule::last_statement {
//                     pbrt_world_end();
//                 } else {
//                     println!("ERROR: parameter expected, {:?} found ...", optional_parameters);
//                 }
//             },
//         }
//         _material(&self) -> () {
//             (name: _string(), optional_parameters) => {
//                 unsafe {
//                     if let Some(ref mut param_set) = PARAM_SET {
//                         if optional_parameters.rule == Rule::statement ||
//                             optional_parameters.rule == Rule::last_statement
//                         {
//                             // println!("Material \"{}\" ", name.clone());
//                         }
//                         // pbrtMaterial (api.cpp:1082)
//                         if let Some(ref mut graphics_state) = GRAPHICS_STATE {
//                             graphics_state.material = name.clone();
//                             graphics_state.current_named_material = String::new();
//                         }
//                         param_set.reset(String::from("Material"),
//                                         String::from(name),
//                                         String::from(""),
//                                         String::from(""));
//                     }
//                 }
//                 if optional_parameters.rule == Rule::parameter {
//                     self._parameter();
//                 } else if optional_parameters.rule == Rule::statement {
//                     self._statement();
//                 } else if optional_parameters.rule == Rule::last_statement {
//                     pbrt_world_end();
//                 } else {
//                     println!("ERROR: parameter expected, {:?} found ...", optional_parameters);
//                 }
//             },
//         }
//         _make_named_material(&self) -> () {
//             (name: _string(), optional_parameters) => {
//                 unsafe {
//                     if let Some(ref mut param_set) = PARAM_SET {
//                         if optional_parameters.rule == Rule::statement ||
//                             optional_parameters.rule == Rule::last_statement
//                         {
//                             // println!("MakeNamedMaterial \"{}\" ", name.clone());
//                         }
//                         // pbrtMakeNamedMaterial (api.cpp:1094)
//                         if let Some(ref mut graphics_state) = GRAPHICS_STATE {
//                             graphics_state.material = name.clone();
//                             graphics_state.current_named_material = String::new();
//                         }
//                         param_set.reset(String::from("MakeNamedMaterial"),
//                                         String::from(name),
//                                         String::from(""),
//                                         String::from(""));
//                     }
//                 }
//                 if optional_parameters.rule == Rule::parameter {
//                     self._parameter();
//                 } else if optional_parameters.rule == Rule::statement {
//                     self._statement();
//                 } else if optional_parameters.rule == Rule::last_statement {
//                     pbrt_world_end();
//                 } else {
//                     println!("ERROR: parameter expected, {:?} found ...", optional_parameters);
//                 }
//             },
//         }
//         _named_material(&self) -> () {
//             (name: _string(), optional_parameters) => {
//                 unsafe {
//                     if let Some(ref mut param_set) = PARAM_SET {
//                         if optional_parameters.rule == Rule::statement ||
//                             optional_parameters.rule == Rule::last_statement
//                         {
//                             // println!("NamedMaterial \"{}\" ", name.clone());
//                         }
//                         // pbrtNamedMaterial (api.cpp:1119)
//                         if let Some(ref mut graphics_state) = GRAPHICS_STATE {
//                             graphics_state.current_named_material = name.clone();
//                         }
//                         param_set.reset(String::from("NamedMaterial"),
//                                         String::from(name),
//                                         String::from(""),
//                                         String::from(""));
//                     }
//                 }
//                 if optional_parameters.rule == Rule::parameter {
//                     self._parameter();
//                 } else if optional_parameters.rule == Rule::statement {
//                     self._statement();
//                 } else if optional_parameters.rule == Rule::last_statement {
//                     pbrt_world_end();
//                 } else {
//                     println!("ERROR: parameter expected, {:?} found ...", optional_parameters);
//                 }
//             },
//         }
//         _shape(&self) -> () {
//             (name: _string(), optional_parameters) => {
//                 unsafe {
//                     if let Some(ref mut param_set) = PARAM_SET {
//                         if optional_parameters.rule == Rule::statement ||
//                             optional_parameters.rule == Rule::last_statement
//                         {
//                             // println!("Shape \"{}\" ", name);
//                             // WARNING: Reset BEFORE calling pbrt_shape() !
//                             param_set.reset(String::from("Shape"),
//                                             String::from(name),
//                                             String::from(""),
//                                             String::from(""));
//                             let (shapes, materials) = pbrt_shape(&param_set);
//                             assert_eq!(shapes.len(), materials.len());
//                             for i in 0..shapes.len() {
//                                 let shape = &shapes[i];
//                                 let material = &materials[i];
//                                 let geo_prim = Arc::new(GeometricPrimitive::new(shape.clone(),
//                                                                                 material.clone(),
//                                                                                 None));
//                                 if let Some(ref mut ro) = RENDER_OPTIONS {
//                                     ro.primitives.push(geo_prim.clone());
//                                 }
//                             }
//                         } else {
//                             param_set.reset(String::from("Shape"),
//                                             String::from(name),
//                                             String::from(""),
//                                             String::from(""));
//                         }
//                     }
//                 }
//                 if optional_parameters.rule == Rule::parameter {
//                     self._parameter();
//                 } else if optional_parameters.rule == Rule::statement {
//                     self._statement();
//                 } else if optional_parameters.rule == Rule::last_statement {
//                     pbrt_world_end();
//                 } else {
//                     println!("ERROR: parameter expected, {:?} found ...", optional_parameters);
//                 }
//             },
//         }
//         // parameters
//         _parameter(&self) -> () {
//             (_head: bool_param, tail: _bool_param()) => {
//                 let (string1, string2) = tail;
//                 unsafe {
//                     if let Some(ref mut param_set) = PARAM_SET {
//                         param_set.add_bool(string1, string2);
//                     }
//                 }
//                 self._parameter();
//             },
//             (_head: float_param, tail: _float_param()) => {
//                 let (string, numbers) = tail;
//                 unsafe {
//                     if let Some(ref mut param_set) = PARAM_SET {
//                         if numbers.len() == 1 {
//                             param_set.add_float(string, numbers[0]);
//                         } else {
//                             param_set.add_floats(string, numbers);
//                         }
//                     }
//                 }
//                 self._parameter();
//             },
//             (_head: string_param, tail: _string_param()) => {
//                 let (string1, string2) = tail;
//                 unsafe {
//                     if let Some(ref mut param_set) = PARAM_SET {
//                         param_set.add_string(string1, string2);
//                     }
//                 }
//                 self._parameter();
//             },
//             (_head: integer_param, tail: _integer_param()) => {
//                 let (string, numbers) = tail;
//                 unsafe {
//                     if let Some(ref mut param_set) = PARAM_SET {
//                         if numbers.len() == 1 {
//                             param_set.add_int(string, numbers[0]);
//                         } else {
//                             param_set.add_ints(string, numbers);
//                         }
//                     }
//                 }
//                 self._parameter();
//             },
//             (_head: point_param, tail: _point_param()) => {
//                 let (string, numbers) = tail;
//                 unsafe {
//                     if let Some(ref mut param_set) = PARAM_SET {
//                         if numbers.len() == 3 {
//                             param_set.add_point3f(string,
//                                                   Point3f {
//                                                       x: numbers[0],
//                                                       y: numbers[1],
//                                                       z: numbers[2],
//                                                   });
//                         } else {
//                             param_set.add_point3fs(string, numbers);
//                         }
//                     }
//                 }
//                 self._parameter();
//             },
//             (_head: vector_param, tail: _vector_param()) => {
//                 let (string, number1, number2, number3) = tail;
//                 unsafe {
//                     if let Some(ref mut param_set) = PARAM_SET {
//                         param_set.add_vector3f(string,
//                                               Vector3f {
//                                                   x: number1,
//                                                   y: number2,
//                                                   z: number3,
//                                               });
//                     }
//                 }
//                 self._parameter();
//             },
//             (_head: normal_param, tail: _normal_param()) => {
//                 let (string, numbers) = tail;
//                 unsafe {
//                     if let Some(ref mut param_set) = PARAM_SET {
//                         if numbers.len() == 3 {
//                             param_set.add_normal3f(string,
//                                                    Normal3f {
//                                                        x: numbers[0],
//                                                        y: numbers[1],
//                                                        z: numbers[2],
//                                                    });
//                         } else {
//                             param_set.add_normal3fs(string, numbers);
//                         }
//                     }
//                 }
//                 self._parameter();
//             },
//             (_head: rgb_param, tail: _rgb_param()) => {
//                 let (string, number1, number2, number3) = tail;
//                 unsafe {
//                     if let Some(ref mut param_set) = PARAM_SET {
//                         param_set.add_rgb_spectrum(string,
//                                                    Spectrum {
//                                                        c: [number1, number2, number3],
//                                                    });
//                     }
//                 }
//                 self._parameter();
//             },
//             (_head: spectrum_param, tail: _spectrum_param()) => {
//                 let string = tail;
//                 print!("\"spectrum\" {} ", string);
//                 // unsafe {
//                 //     if let Some(ref mut param_set) = PARAM_SET {
//                 //         param_set.add_string(string1, string2);
//                 //     }
//                 // }
//                 self._parameter();
//             },
//             (_head: texture_param, tail: _texture_param()) => {
//                 let (string1, string2) = tail;
//                 unsafe {
//                     if let Some(ref mut param_set) = PARAM_SET {
//                         param_set.add_texture(string1, string2);
//                     }
//                 }
//                 self._parameter();
//             },
//             (optional_parameters) => {
//                 if optional_parameters.rule == Rule::statement ||
//                     optional_parameters.rule == Rule::last_statement {
//                     unsafe {
//                         if let Some(ref mut param_set) = PARAM_SET {
//                             let mut name: String = String::new();
//                             name.push_str(param_set.name.as_str());
//                             if param_set.key_word == String::from("Accelerator") {
//                                 if let Some(ref mut ro) = RENDER_OPTIONS {
//                                     // println!("Accelerator \"{}\" ", ro.accelerator_name);
//                                     ro.accelerator_params.copy_from(param_set);
//                                     print_params(&ro.accelerator_params);
//                                 }
//                             } else if param_set.key_word == String::from("Camera") {
//                                 if let Some(ref mut ro) = RENDER_OPTIONS {
//                                     // println!("Camera \"{}\" ", ro.camera_name);
//                                     ro.camera_params.copy_from(param_set);
//                                     print_params(&ro.camera_params);
//                                 }
//                             } else if param_set.key_word == String::from("PixelFilter") {
//                                 if let Some(ref mut ro) = RENDER_OPTIONS {
//                                     // println!("PixelFilter \"{}\" ", ro.filter_name);
//                                     ro.filter_params.copy_from(param_set);
//                                     print_params(&ro.filter_params);
//                                 }
//                             } else if param_set.key_word == String::from("Sampler") {
//                                 if let Some(ref mut ro) = RENDER_OPTIONS {
//                                     // println!("Sampler \"{}\" ", ro.sampler_name);
//                                     ro.sampler_params.copy_from(param_set);
//                                     print_params(&ro.sampler_params);
//                                 }
//                             } else if param_set.key_word == String::from("Film") {
//                                 if let Some(ref mut ro) = RENDER_OPTIONS {
//                                     // println!("Film \"{}\" ", ro.film_name);
//                                     ro.film_params.copy_from(param_set);
//                                     print_params(&ro.film_params);
//                                 }
//                             } else if param_set.key_word == String::from("Integrator") {
//                                 if let Some(ref mut ro) = RENDER_OPTIONS {
//                                     // println!("Integrator \"{}\" ", ro.integrator_name);
//                                     ro.integrator_params.copy_from(param_set);
//                                     print_params(&ro.integrator_params);
//                                 }
//                             } else if param_set.key_word == String::from("AreaLightSource") {
//                                 // println!("AreaLightSource \"{}\" ", param_set.name);
//                                 print_params(&param_set);
//                                 // pbrtAreaLightSource (api.cpp:1142)
//                                 if let Some(ref mut graphics_state) = GRAPHICS_STATE {
//                                     graphics_state.area_light = name;
//                                     graphics_state.area_light_params.copy_from(&param_set);
//                                 }
//                             } else if param_set.key_word == String::from("LightSource") {
//                                 if let Some(ref mut ro) = RENDER_OPTIONS {
//                                     // println!("LightSource \"{}\" ", param_set.name);
//                                     print_params(&param_set);
//                                     make_light(&param_set, ro);
//                                 }
//                                 // let lt = make_light(name, params, CUR_TRANSFORM.t[0]);
//                             } else if param_set.key_word == String::from("Texture") {
//                                 // println!("Texture \"{}\" \"{}\" \"{}\" ",
//                                 //          param_set.name,
//                                 //          param_set.tex_type,
//                                 //          param_set.tex_name);
//                                 print_params(&param_set);
//                                 // pbrtTexture (api.cpp:1049)
//                                 if let Some(ref mut graphics_state) = GRAPHICS_STATE {
//                                     let mut geom_params: ParamSet = ParamSet::default();
//                                     let mut material_params: ParamSet = ParamSet::default();
//                                     geom_params.copy_from(param_set);
//                                     material_params.copy_from(param_set);
//                                     let mut tp: TextureParams = TextureParams {
//                                         float_textures: graphics_state.float_textures.clone(),
//                                         spectrum_textures: graphics_state.spectrum_textures.clone(),
//                                         geom_params: geom_params,
//                                         material_params: material_params,
//                                     };
//                                     if param_set.tex_type == String::from("float") {
//                                         println!("TODO: MakeFloatTexture");
//                                     } else if param_set.tex_type == String::from("color") ||
//                                         param_set.tex_type == String::from("spectrum") {
//                                             match graphics_state.spectrum_textures.get(param_set.name.as_str()) {
//                                                 Some(_spectrum_texture) => {
//                                                     println!("Texture \"{}\" being redefined",
//                                                              param_set.name);
//                                                 },
//                                                 None => {},
//                                             }
//                                             // MakeSpectrumTexture(texname, curTransform[0], tp);
//                                             if param_set.tex_name == String::from("constant") {
//                                                 println!("TODO: CreateConstantSpectrumTexture");
//                                             } else if param_set.tex_name == String::from("scale") {
//                                                 println!("TODO: CreateScaleSpectrumTexture");
//                                             } else if param_set.tex_name == String::from("mix") {
//                                                 println!("TODO: CreateMixSpectrumTexture");
//                                             } else if param_set.tex_name == String::from("bilerp") {
//                                                 println!("TODO: CreateBilerpSpectrumTexture");
//                                             } else if param_set.tex_name == String::from("imagemap") {
//                                                 // CreateImageSpectrumTexture
//                                                 let mut map: Option<Box<TextureMapping2D + Send + Sync>> = None;
//                                                 let mapping: String =
//                                                     tp.find_string(String::from("mapping"), String::from("uv"));
//                                                 if mapping == String::from("uv") {
//                                                     let su: Float = tp.find_float(String::from("uscale"), 1.0);
//                                                     let sv: Float = tp.find_float(String::from("vscale"), 1.0);
//                                                     let du: Float = tp.find_float(String::from("udelta"), 0.0);
//                                                     let dv: Float = tp.find_float(String::from("vdelta"), 0.0);
//                                                     map = Some(Box::new(UVMapping2D {
//                                                         su: su,
//                                                         sv: sv,
//                                                         du: du,
//                                                         dv: dv,
//                                                     }));
//                                                 } else if mapping == String::from("spherical") {
//                                                     println!("TODO: SphericalMapping2D");
//                                                 } else if mapping == String::from("cylindrical") {
//                                                     println!("TODO: CylindricalMapping2D");
//                                                 } else if mapping == String::from("planar") {
//                                                     map = Some(Box::new(PlanarMapping2D {
//                                                         vs: tp.find_vector3f(String::from("v1"),
//                                                                              Vector3f {
//                                                                                  x: 1.0,
//                                                                                  y: 0.0,
//                                                                                  z: 0.0
//                                                                              }),
//                                                         vt: tp.find_vector3f(String::from("v2"),
//                                                                              Vector3f {
//                                                                                  x: 0.0,
//                                                                                  y: 1.0,
//                                                                                  z: 0.0
//                                                                              }),
//                                                         ds: tp.find_float(String::from("udelta"),
//                                                                           0.0),
//                                                         dt: tp.find_float(String::from("vdelta"),
//                                                                           0.0),
//                                                     }));
//                                                 } else {
//                                                     panic!("2D texture mapping \"{}\" unknown",
//                                                            mapping);
//                                                 }
//                                                 // initialize _ImageTexture_ parameters
//                                                 let max_aniso: Float =
//                                                     tp.find_float(String::from("maxanisotropy"), 8.0);
//                                                 let do_trilinear: bool =
//                                                     tp.find_bool(String::from("trilinear"), false);
//                                                 let wrap: String =
//                                                     tp.find_string(String::from("wrap"),
//                                                                    String::from("repeat"));
//                                                 let mut wrap_mode: ImageWrap = ImageWrap::Repeat;
//                                                 if wrap == String::from("black") {
//                                                     wrap_mode = ImageWrap::Black;
//                                                 }
//                                                 else if wrap == String::from("clamp") {
//                                                     wrap_mode = ImageWrap::Clamp;
//                                                 }
//                                                 let scale: Float =
//                                                     tp.find_float(String::from("scale"), 1.0);
//                                                 let mut filename: String =
//                                                     tp.find_filename(String::from("filename"),
//                                                                      String::new());
//                                                 if let Some(ref search_directory) = SEARCH_DIRECTORY {
//                                                     // filename = AbsolutePath(ResolveFilename(filename));
//                                                     let mut path_buf: PathBuf = PathBuf::from("/");
//                                                     path_buf.push(search_directory.as_ref());
//                                                     path_buf.push(filename);
//                                                     filename = String::from(path_buf.to_str().unwrap());
//                                                 }
//                                                 // TODO: default depends on:
//                                                 // HasExtension(filename,
//                                                 // ".tga") ||
//                                                 // HasExtension(filename,
//                                                 // ".png"));
//                                                 let gamma: bool = tp.find_bool(String::from("gamma"), true);

//                                                 if let Some(mapping) = map {
//                                                     let st = Arc::new(ImageTexture::new(mapping,
//                                                                                         filename,
//                                                                                         do_trilinear,
//                                                                                         max_aniso,
//                                                                                         wrap_mode,
//                                                                                         scale,
//                                                                                         gamma));
//                                                     graphics_state.spectrum_textures.insert(
//                                                         name, st);
//                                                 }
//                                             } else if param_set.tex_name == String::from("uv") {
//                                                 println!("TODO: CreateUVSpectrumTexture");
//                                             } else if param_set.tex_name == String::from("checkerboard") {
//                                                 // CreateCheckerboardSpectrumTexture
//                                                 let dim: i32 = tp.find_int(String::from("dimension"), 2);
//                                                 if dim != 2 && dim != 3 {
//                                                     panic!("{} dimensional checkerboard texture not supported",
//                                                            dim);
//                                                 }
//                                                 let tex1: Arc<Texture<Spectrum> + Send + Sync> =
//                                                     tp.get_spectrum_texture(String::from("tex1"),
//                                                                             Spectrum::new(1.0));
//                                                 let tex2: Arc<Texture<Spectrum> + Send + Sync> =
//                                                     tp.get_spectrum_texture(String::from("tex2"),
//                                                                             Spectrum::new(0.0));
//                                                 if dim == 2 {
//                                                     let mut map: Option<Box<TextureMapping2D + Send + Sync>> = None;
//                                                     let mapping: String =
//                                                         tp.find_string(String::from("mapping"), String::from("uv"));
//                                                     if mapping == String::from("uv") {
//                                                         println!("TODO: UVMapping2D");
//                                                     } else if mapping == String::from("spherical") {
//                                                         println!("TODO: SphericalMapping2D");
//                                                     } else if mapping == String::from("cylindrical") {
//                                                         println!("TODO: CylindricalMapping2D");
//                                                     } else if mapping == String::from("planar") {
//                                                         map = Some(Box::new(PlanarMapping2D {
//                                                             vs: tp.find_vector3f(String::from("v1"),
//                                                                                  Vector3f {
//                                                                                      x: 1.0,
//                                                                                      y: 0.0,
//                                                                                      z: 0.0
//                                                                                  }),
//                                                             vt: tp.find_vector3f(String::from("v2"),
//                                                                                  Vector3f {
//                                                                                      x: 0.0,
//                                                                                      y: 1.0,
//                                                                                      z: 0.0
//                                                                                  }),
//                                                             ds: tp.find_float(String::from("udelta"),
//                                                                               0.0),
//                                                             dt: tp.find_float(String::from("vdelta"),
//                                                                               0.0),
//                                                         }));
//                                                     } else {
//                                                         panic!("2D texture mapping \"{}\" unknown",
//                                                                mapping);
//                                                     }
//                                                     // TODO: aamode
//                                                     if let Some(mapping) = map {
//                                                         Arc::new(Checkerboard2DTexture::new(mapping,
//                                                                                             tex1,
//                                                                                             tex2));
//                                                     }
//                                                 } else { // dim == 3
//                                                     println!("TODO: TextureMapping3D");
//                                                 }
//                                             } else if param_set.tex_name == String::from("dots") {
//                                                 println!("TODO: CreateDotsSpectrumTexture");
//                                             } else if param_set.tex_name == String::from("fbm") {
//                                                 println!("TODO: CreateFBmSpectrumTexture");
//                                             } else if param_set.tex_name == String::from("wrinkled") {
//                                                 println!("TODO: CreateWrinkledSpectrumTexture");
//                                             } else if param_set.tex_name == String::from("marble") {
//                                                 println!("TODO: CreateMarbleSpectrumTexture");
//                                             } else if param_set.tex_name == String::from("windy") {
//                                                 println!("TODO: CreateWindySpectrumTexture");
//                                             } else {
//                                                 println!("Spectrum texture \"{}\" unknown.",
//                                                          param_set.tex_name);
//                                             }
//                                         } else {
//                                             panic!("Texture type \"{}\" unknown.", param_set.tex_type);
//                                         }
//                                 }
//                                 // MakeFloatTexture(texname, curTransform[0], tp);
//                                 // or
//                                 // MakeSpectrumTexture(texname, curTransform[0], tp);
//                             } else if param_set.key_word == String::from("Material") {
//                                 // println!("Material \"{}\" ", param_set.name);
//                                 print_params(&param_set);
//                                 // pbrtMaterial (api.cpp:1082)
//                                 if let Some(ref mut graphics_state) = GRAPHICS_STATE {
//                                     graphics_state.material = name;
//                                     graphics_state.material_params.copy_from(&param_set);
//                                     graphics_state.current_named_material = String::new();
//                                 }
//                             } else if param_set.key_word == String::from("MakeNamedMaterial") {
//                                 // println!("MakeNamedMaterial \"{}\" ", param_set.name);
//                                 print_params(&param_set);
//                                 // pbrtMakeNamedMaterial (api.cpp:1094)
//                                 let mat_type: String = param_set.find_one_string(String::from("type"),
//                                                                                  String::new());
//                                 if mat_type == String::new() {
//                                     panic!("No parameter string \"type\" found in MakeNamedMaterial");
//                                 }
//                                 if let Some(ref mut graphics_state) = GRAPHICS_STATE {
//                                     graphics_state.material = mat_type.clone();
//                                     graphics_state.material_params.copy_from(&param_set);
//                                     let mtl: Arc<Material + Send + Sync> = create_material();
//                                     match graphics_state.named_materials.get(mat_type.as_str()) {
//                                         Some(_named_material) => {
//                                             // println!("Named material \"{}\" redefined",
//                                             //          mat_type);
//                                         },
//                                         None => {},
//                                     }
//                                     graphics_state.named_materials.insert(param_set.name.clone(), mtl);
//                                 }
//                             } else if param_set.key_word == String::from("NamedMaterial") {
//                                 // println!("NamedMaterial \"{}\" ", param_set.name);
//                                 print_params(&param_set);
//                                 // pbrtNamedMaterial (api.cpp:1119)
//                                 if let Some(ref mut graphics_state) = GRAPHICS_STATE {
//                                     graphics_state.current_named_material = param_set.name.clone();
//                                 }
//                             } else if param_set.key_word == String::from("Shape") {
//                                 // println!("Shape \"{}\" ", param_set.name);
//                                 print_params(&param_set);
//                                 // collect area lights
//                                 let mut area_lights: Vec<Arc<Light + Send + Sync>> = Vec::new();
//                                 // possibly create area light for shape (see pbrtShape())
//                                 if let Some(ref mut graphics_state) = GRAPHICS_STATE {
//                                     if graphics_state.area_light != String::new() {
//                                         // MakeAreaLight
//                                         if graphics_state.area_light == String::from("area") ||
//                                             graphics_state.area_light == String::from("diffuse")
//                                         {
//                                             // first create the shape
//                                             let (shapes, materials) = pbrt_shape(&param_set);
//                                             assert_eq!(shapes.len(), materials.len());
//                                             for i in 0..shapes.len() {
//                                                 let shape = &shapes[i];
//                                                 let material = &materials[i];
//                                                 // CreateDiffuseAreaLight
//                                                 let light_to_world: Transform = CUR_TRANSFORM.t[0];
//                                                 let l: Spectrum =
//                                                     graphics_state.area_light_params.find_one_spectrum(String::from("L"),
//                                                                                                        Spectrum::new(1.0));
//                                                 let sc: Spectrum =
//                                                     graphics_state.area_light_params.find_one_spectrum(String::from("scale"),
//                                                                                                        Spectrum::new(1.0));
//                                                 let n_samples: i32 = // try "nsamples" first
//                                                     graphics_state.area_light_params.find_one_int(String::from("nsamples"),
//                                                                                                   1);
//                                                 let n_samples: i32 = // try "samples"next
//                                                     graphics_state.area_light_params.find_one_int(String::from("samples"),
//                                                                                                   n_samples);
//                                                 let two_sided: bool =
//                                                     graphics_state.area_light_params.find_one_bool(String::from("twosided"),
//                                                                                                    false);
//                                                 // TODO: if (PbrtOptions.quickRender) nSamples = std::max(1, nSamples / 4);
//                                                 let l_emit: Spectrum = l * sc;
//                                                 let area_light: Arc<DiffuseAreaLight> =
//                                                     Arc::new(DiffuseAreaLight::new(
//                                                         &light_to_world,
//                                                         &l_emit,
//                                                         n_samples,
//                                                         shape.clone(),
//                                                         two_sided
//                                                     ));
//                                                 area_lights.push(area_light.clone());
//                                                 let geo_prim = Arc::new(GeometricPrimitive::new(shape.clone(),
//                                                                                                 material.clone(),
//                                                                                                 Some(area_light.clone())));
//                                                 if let Some(ref mut ro) = RENDER_OPTIONS {
//                                                     ro.primitives.push(geo_prim.clone());
//                                                 }
//                                             }
//                                         }
//                                     } else {
//                                         // continue with shape itself
//                                         let (shapes, materials) = pbrt_shape(&param_set);
//                                         assert_eq!(shapes.len(), materials.len());
//                                         for i in 0..shapes.len() {
//                                             let shape = &shapes[i];
//                                             let material = &materials[i];
//                                             let geo_prim = Arc::new(GeometricPrimitive::new(shape.clone(),
//                                                                                             material.clone(),
//                                                                                             None));
//                                             if let Some(ref mut ro) = RENDER_OPTIONS {
//                                                 ro.primitives.push(geo_prim.clone());
//                                             }
//                                         }
//                                     }
//                                 } else {
//                                     // continue with shape itself
//                                     let (shapes, materials) = pbrt_shape(&param_set);
//                                     assert_eq!(shapes.len(), materials.len());
//                                     for i in 0..shapes.len() {
//                                         let shape = &shapes[i];
//                                         let material = &materials[i];
//                                         let geo_prim = Arc::new(GeometricPrimitive::new(shape.clone(),
//                                                                                         material.clone(),
//                                                                                         None));
//                                         if let Some(ref mut ro) = RENDER_OPTIONS {
//                                             ro.primitives.push(geo_prim.clone());
//                                         }
//                                     }
//                                 }
//                                 // add _prims_ and _areaLights_ to scene or current instance
//                                 // if (renderOptions->currentInstance) {
//                                 //     if (areaLights.size())
//                                 //         Warning("Area lights not supported with object instancing");
//                                 //     renderOptions->currentInstance->insert(
//                                 //         renderOptions->currentInstance->end(), prims.begin(), prims.end());
//                                 // } else {
//                                 if let Some(ref mut ro) = RENDER_OPTIONS {
//                                     // ro.primitives.insert(ro.primitives.end(),
//                                     //                      prims.begin(), prims.end());
//                                     if area_lights.len() > 0 {
//                                         for area_light in area_lights {
//                                             ro.lights.push(area_light);
//                                         }
//                                     }
//                                 }
//                             } else {
//                                 println!("PARAM_SET: {}", param_set.key_word);
//                             }
//                             param_set.reset(String::from(""),
//                                             String::from(""),
//                                             String::from(""),
//                                             String::from(""));
//                         }
//                     }
//                     if optional_parameters.rule == Rule::last_statement {
//                         pbrt_world_end();
//                     } else { // statement
//                         self._statement();
//                     }
//                 } else if optional_parameters.rule == Rule::parameter {
//                     self._parameter();
//                 } else {
//                     println!("ERROR: statement or parameter expected, {:?} found ...",
//                              optional_parameters);
//                 }
//             }
//         }
//         _bool_param(&self) -> (String, bool) {
//             (&i1: ident, _l: lbrack, _s: string, &i2: ident, _r: rbrack) => {
//                 let string1: String = String::from_str(i1).unwrap();
//                 let string2: String = String::from_str(i2).unwrap();
//                 let b: bool;
//                 if string2 == String::from("true") {
//                     b = true;
//                 } else if string2 == String::from("false") {
//                     b = false
//                 } else {
//                     println!("WARNING: parameter {:?} not well defined, defaulting to false", string1);
//                     b = false
//                 }
//                 (string1, b)
//             },
//         }
//         _float_param(&self) -> (String, Vec<Float>) {
//             // single float without brackets
//             (&i: ident, &n: number) => {
//                 let string: String = String::from_str(i).unwrap();
//                 let number: Float = f32::from_str(n).unwrap();
//                 let fvec: Vec<Float> = vec![number];
//                 (string, fvec)
//             },
//             // single float
//             (&i: ident, _l: lbrack, &n: number, _r: rbrack) => {
//                 let string: String = String::from_str(i).unwrap();
//                 let number: Float = f32::from_str(n).unwrap();
//                 let fvec: Vec<Float> = vec![number];
//                 (string, fvec)
//             },
//             // more than one float
//             (&i: ident, _l: lbrack, mut list: _list_of_floats(), _r: rbrack) => {
//                 let string: String = String::from_str(i).unwrap();
//                 let mut fvec: Vec<Float> = Vec::new();
//                 if let Some(floats) = list.pop_front() {
//                     match floats {
//                         FloatNode::Floats(list) => {
//                             for element in list.iter() {
//                                 match *element {
//                                     FloatNode::Floats(ref _l) => {},
//                                     FloatNode::OneFloat(ref f) => {
//                                         fvec.push(*f);
//                                     },
//                                 };
//                             }
//                         },
//                         FloatNode::OneFloat(_f) => {},
//                     };
//                 }
//                 (string, fvec)
//             },
//         }
//         _list_of_floats(&self) -> LinkedList<FloatNode> {
//             (&n: number, mut head: _float(), mut tail: _list_of_floats()) => {
//                 let number: Float = Float::from_str(n).unwrap();
//                 head.push_front(FloatNode::OneFloat(number));
//                 tail.push_front(FloatNode::Floats(head));
//                 tail
//             },
//             () => {
//                 LinkedList::new()
//             }
//         }
//         _float(&self) -> LinkedList<FloatNode> {
//             (&head: number, mut tail: _float()) => {
//                 let number: Float = Float::from_str(head).unwrap();
//                 tail.push_front(FloatNode::OneFloat(number));
//                 tail
//             },
//             () => {
//                 LinkedList::new()
//             }
//         }
//         _string_param(&self) -> (String, String) {
//             (&i: ident, _l: lbrack, _s: string, &f: filename, _r: rbrack) => {
//                 let string1: String = String::from_str(i).unwrap();
//                 let string2: String = String::from_str(f).unwrap();
//                 (string1, string2)
//             },
//             (&i1: ident, _l: lbrack, _s: string, &i2: ident, _r: rbrack) => {
//                 let string1: String = String::from_str(i1).unwrap();
//                 let string2: String = String::from_str(i2).unwrap();
//                 (string1, string2)
//             },
//             (&i: ident, _s: string, &f: filename) => {
//                 let string1: String = String::from_str(i).unwrap();
//                 let string2: String = String::from_str(f).unwrap();
//                 (string1, string2)
//             },
//             (&i1: ident, _s: string, &i2: ident) => {
//                 let string1: String = String::from_str(i1).unwrap();
//                 let string2: String = String::from_str(i2).unwrap();
//                 (string1, string2)
//             },
//         }
//         _integer_param(&self) -> (String, Vec<i32>) {
//             // single integer
//             (&i: ident, _l: lbrack, &n: integer, _r: rbrack) => {
//                 let string: String = String::from_str(i).unwrap();
//                 let number: i32 = i32::from_str(n).unwrap();
//                 let ivec: Vec<i32> = vec![number];
//                 (string, ivec)
//             },
//             // more than one integer
//             (&i: ident, _l: lbrack, mut list: _list_of_integers(), _r: rbrack) => {
//                 let string: String = String::from_str(i).unwrap();
//                 let mut ivec: Vec<i32> = Vec::new();
//                 if let Some(integers) = list.pop_front() {
//                     match integers {
//                         IntNode::Ints(list) => {
//                             for element in list.iter() {
//                                 match *element {
//                                     IntNode::Ints(ref _l) => {},
//                                     IntNode::OneInt(ref i) => {
//                                         ivec.push(*i);
//                                     },
//                                 };
//                             }
//                         },
//                         IntNode::OneInt(_i) => {},
//                     };
//                 }
//                 (string, ivec)
//             },
//         }
//         _list_of_integers(&self) -> LinkedList<IntNode> {
//             (&i: integer, mut head: _integer(), mut tail: _list_of_integers()) => {
//                 let number: i32 = i32::from_str(i).unwrap();
//                 head.push_front(IntNode::OneInt(number));
//                 tail.push_front(IntNode::Ints(head));
//                 tail
//             },
//             () => {
//                 LinkedList::new()
//             }
//         }
//         _integer(&self) -> LinkedList<IntNode> {
//             (&head: integer, mut tail: _integer()) => {
//                 let number: i32 = i32::from_str(head).unwrap();
//                 tail.push_front(IntNode::OneInt(number));
//                 tail
//             },
//             () => {
//                 LinkedList::new()
//             }
//         }
//         _point_param(&self) -> (String, Vec<Float>) {
//             // single point
//             (&i: ident, _l: lbrack, &n1: number, &n2: number, &n3: number, _r: rbrack) => {
//                 let string: String = String::from_str(i).unwrap();
//                 let number1: Float = f32::from_str(n1).unwrap();
//                 let number2: Float = f32::from_str(n2).unwrap();
//                 let number3: Float = f32::from_str(n3).unwrap();
//                 let fvec: Vec<Float> = vec![number1, number2, number3];
//                 (string, fvec)
//             },
//             // more than one point (3 floats define one point)
//             (&i: ident, _l: lbrack, mut list: _list_of_floats(), _r: rbrack) => {
//                 let string: String = String::from_str(i).unwrap();
//                 let mut fvec: Vec<Float> = Vec::new();
//                 if let Some(floats) = list.pop_front() {
//                     match floats {
//                         FloatNode::Floats(list) => {
//                             for element in list.iter() {
//                                 match *element {
//                                     FloatNode::Floats(ref _l) => {},
//                                     FloatNode::OneFloat(ref f) => {
//                                         fvec.push(*f);
//                                     },
//                                 };
//                             }
//                         },
//                         FloatNode::OneFloat(_f) => {},
//                     };
//                 }
//                 assert!(fvec.len() % 3 == 0, "point parameters need 3 coordinates");
//                 (string, fvec)
//             },
//         }
//         _vector_param(&self) -> (String, Float, Float, Float) {
//             (&i: ident, _l: lbrack, &n1: number, &n2: number, &n3: number, _r: rbrack) => {
//                 let string: String = String::from_str(i).unwrap();
//                 let number1: Float = f32::from_str(n1).unwrap();
//                 let number2: Float = f32::from_str(n2).unwrap();
//                 let number3: Float = f32::from_str(n3).unwrap();
//                 (string, number1, number2, number3)
//             },
//         }
//         _normal_param(&self) -> (String, Vec<Float>) {
//             // single normal
//             (&i: ident, _l: lbrack, &n1: number, &n2: number, &n3: number, _r: rbrack) => {
//                 let string: String = String::from_str(i).unwrap();
//                 let number1: Float = f32::from_str(n1).unwrap();
//                 let number2: Float = f32::from_str(n2).unwrap();
//                 let number3: Float = f32::from_str(n3).unwrap();
//                 let fvec: Vec<Float> = vec![number1, number2, number3];
//                 (string, fvec)
//             },
//             // more than one normal (3 floats define one normal)
//             (&i: ident, _l: lbrack, mut list: _list_of_floats(), _r: rbrack) => {
//                 let string: String = String::from_str(i).unwrap();
//                 let mut fvec: Vec<Float> = Vec::new();
//                 if let Some(floats) = list.pop_front() {
//                     match floats {
//                         FloatNode::Floats(list) => {
//                             for element in list.iter() {
//                                 match *element {
//                                     FloatNode::Floats(ref _l) => {},
//                                     FloatNode::OneFloat(ref f) => {
//                                         fvec.push(*f);
//                                     },
//                                 };
//                             }
//                         },
//                         FloatNode::OneFloat(_f) => {},
//                     };
//                 }
//                 assert!(fvec.len() % 3 == 0, "normal parameters need 3 coordinates");
//                 (string, fvec)
//             },
//         }
//         _rgb_param(&self) -> (String, Float, Float, Float) {
//             (&i: ident, _l: lbrack, &n1: number, &n2: number, &n3: number, _r: rbrack) => {
//                 let string: String = String::from_str(i).unwrap();
//                 let number1: Float = f32::from_str(n1).unwrap();
//                 let number2: Float = f32::from_str(n2).unwrap();
//                 let number3: Float = f32::from_str(n3).unwrap();
//                 (string, number1, number2, number3)
//             },
//         }
//         _spectrum_param(&self) -> String {
//             (&s: string, _i: ident) => {
//                 let string: String = String::from_str(s).unwrap();
//                 string
//             },
//         }
//         _texture_param(&self) -> (String, String) {
//             (&i1: ident, _s: string, &i2: ident) => {
//                 let string1: String = String::from_str(i1).unwrap();
//                 let string2: String = String::from_str(i2).unwrap();
//                 (string1, string2)
//             },
//         }
//         // keywords
//         _keyword(&self) -> () {
//             (_ab: attribute_begin) => {
//                 // println!("AttributeBegin");
//                 unsafe {
//                     if let Some(ref mut graphics_state) = GRAPHICS_STATE {
//                         if let Some(ref mut pushed_graphics_states) = PUSHED_GRAPHICS_STATES {
//                             let mut material_param_set: ParamSet = ParamSet::default();
//                             material_param_set.copy_from(&graphics_state.material_params);
//                             let mut area_light_param_set: ParamSet = ParamSet::default();
//                             area_light_param_set.copy_from(&graphics_state.area_light_params);
//                             pushed_graphics_states.push(GraphicsState {
//                                 float_textures: graphics_state.float_textures.clone(),
//                                 spectrum_textures: graphics_state.spectrum_textures.clone(),
//                                 material_params: material_param_set,
//                                 material: String::from(graphics_state.material.as_ref()),
//                                 named_materials: graphics_state.named_materials.clone(),
//                                 current_named_material: String::from(graphics_state.current_named_material.as_ref()),
//                                 area_light_params: area_light_param_set,
//                                 area_light: String::from(graphics_state.area_light.as_ref()),
//                             });
//                         }
//                         if let Some(ref mut pt) = PUSHED_TRANSFORMS {
//                             pt.push(TransformSet {
//                                 t: [
//                                     Transform {
//                                         m: CUR_TRANSFORM.t[0].m,
//                                         m_inv: CUR_TRANSFORM.t[0].m_inv,},
//                                     Transform {
//                                         m: CUR_TRANSFORM.t[1].m,
//                                         m_inv: CUR_TRANSFORM.t[1].m_inv,},
//                                     ]
//                             });
//                         }
//                         // TODO? pushedActiveTransformBits.push_back(activeTransformBits);
//                     }
//                 }
//                 self._pbrt();
//             },
//             (_ae: attribute_end) => {
//                 // println!("AttributeEnd");
//                 unsafe {
//                     if let Some(ref mut graphics_state) = GRAPHICS_STATE {
//                         if let Some(ref mut pushed_graphics_states) = PUSHED_GRAPHICS_STATES {
//                             if !(pushed_graphics_states.len() >= 1_usize) {
//                                 panic!("Unmatched pbrtAttributeEnd() encountered.")
//                             }
//                             let pgs: GraphicsState = pushed_graphics_states.pop().unwrap();
//                             // material_params
//                             graphics_state.material_params.reset(String::new(),
//                                                                  String::from(""),
//                                                                  String::from(""),
//                                                                  String::new());
//                             graphics_state.material_params.copy_from(&pgs.material_params);
//                             // material
//                             graphics_state.material = String::from(pgs.material.as_ref());
//                             // area_light_params
//                             graphics_state.area_light_params.reset(String::new(),
//                                                                  String::from(""),
//                                                                  String::from(""),
//                                                                  String::new());
//                             graphics_state.area_light_params.copy_from(&pgs.area_light_params);
//                             // area_light
//                             graphics_state.area_light = String::from(pgs.area_light.as_ref());
//                         }
//                         if let Some(ref mut pt) = PUSHED_TRANSFORMS {
//                             let popped_transform_set: TransformSet = pt.pop().unwrap();
//                             CUR_TRANSFORM.t[0] = popped_transform_set.t[0];
//                             CUR_TRANSFORM.t[1] = popped_transform_set.t[1];
//                             // println!("CUR_TRANSFORM: {:?}", CUR_TRANSFORM);
//                         }
//                         // TODO? pushedActiveTransformBits.push_back(activeTransformBits);
//                     }
//                 }
//                 self._pbrt();
//             },
//             (_wb: world_begin) => {
//                 // println!("WorldBegin");
//                 unsafe {
//                     CUR_TRANSFORM.t[0] = Transform::default();
//                     CUR_TRANSFORM.t[1] = Transform::default();
//                     // println!("CUR_TRANSFORM: {:?}", CUR_TRANSFORM);
//                     if let Some(ref mut named_coordinate_systems) = NAMED_COORDINATE_SYSTEMS {
//                         named_coordinate_systems.insert("world",
//                                                         TransformSet {
//                                                             t: [Transform::default(); 2]
//                                                         });
//                         // println!("NAMED_COORDINATE_SYSTEMS: {:?}", named_coordinate_systems);
//                     }
//                 }
//                 self._pbrt();
//             },
//             (t) => { println!("TODO: {:?}", t); },
//         }
//         // numbers
//         _number(&self) -> Float {
//             (&n: number) => {
//                 let number: Float = f32::from_str(n).unwrap();
//                 number
//             },
//         }
//         // strings
//         _string(&self) -> String {
//             (_s: string) => {
//                 self._ident()
//             },
//         }
//         // identifiers
//         _ident(&self) -> String {
//             (&i: ident) => {
//                 let string: String = String::from_str(i).unwrap();
//                 string
//             },
//         }
//     }
// }

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
    for p in &params.normals {
        if p.n_values == 1_usize {
            println!("  \"normal {}\" [{} {} {}]",
                     p.name,
                     p.values[0].x,
                     p.values[0].y,
                     p.values[0].z);
        } else {
            println!("  \"normal {}\" [", p.name);
            for i in 0..p.n_values {
                println!("    {} {} {} ", p.values[i].x, p.values[i].y, p.values[i].z);
            }
            println!("  ]");
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
                match graphics_state
                          .named_materials
                          .get(graphics_state.current_named_material.as_str()) {
                    Some(named_material) => {
                        return named_material.clone();
                    }
                    None => {
                        println!("WARNING: Named material \"{}\" not defined. Using \"matte\".",
                                 graphics_state.current_named_material);
                    }
                }
            } else {
                // MakeMaterial
                assert_ne!(graphics_state.material, String::new());
                assert_ne!(graphics_state.material, String::from("none"));
                if graphics_state.material == String::from("matte") {
                    let kd = mp.get_spectrum_texture(String::from("Kd"),
                                                     Spectrum::new(0.5 as Float));
                    // TODO: std::shared_ptr<Texture<Float>> sigma = mp.GetFloatTexture("sigma", 0.f);
                    let matte = Arc::new(MatteMaterial::new(kd, 0.0 as Float));
                    return matte;
                } else if graphics_state.material == String::from("plastic") {
                    let kd = mp.get_spectrum_texture(String::from("Kd"),
                                                     Spectrum::new(0.25 as Float));
                    let ks = mp.get_spectrum_texture(String::from("Ks"),
                                                     Spectrum::new(0.25 as Float));
                    let roughness = mp.get_float_texture(String::from("roughness"), 0.1 as Float);
                    // TODO: std::shared_ptr<Texture<Float>> bumpMap = mp.GetFloatTextureOrNull("bumpmap");
                    let remap_roughness = mp.find_bool(String::from("remaproughness"), true);
                    let plastic =
                        Arc::new(PlasticMaterial::new(kd, ks, roughness, remap_roughness));
                    return plastic;
                } else if graphics_state.material == String::from("translucent") {
                    println!("TODO: CreateTranslucentMaterial");
                } else if graphics_state.material == String::from("glass") {
                    let kr = mp.get_spectrum_texture(String::from("Kr"),
                                                     Spectrum::new(1.0 as Float));
                    let kt = mp.get_spectrum_texture(String::from("Kt"),
                                                     Spectrum::new(1.0 as Float));
                    // let some_eta = mp.get_float_texture_or_null(String::from("eta"));
                    // if let Some(eta) = some_eta {
                    //     println!("some eta");
                    // } else {
                    let eta = mp.get_float_texture(String::from("index"), 1.5);
                    // }
                    // std::shared_ptr<Texture<Float>> roughu =
                    //     mp.GetFloatTexture("uroughness", 0.f);
                    // std::shared_ptr<Texture<Float>> roughv =
                    //     mp.GetFloatTexture("vroughness", 0.f);
                    // std::shared_ptr<Texture<Float>> bumpMap =
                    //     mp.GetFloatTextureOrNull("bumpmap");
                    let remap_roughness: bool = mp.find_bool(String::from("remaproughness"), true);
                    let glass = Arc::new(GlassMaterial {
                                             kr: kr,
                                             kt: kt,
                                             u_roughness: 0.0 as Float,
                                             v_roughness: 0.0 as Float,
                                             index: eta,
                                             remap_roughness: remap_roughness,
                                         });
                    return glass;
                } else if graphics_state.material == String::from("mirror") {
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

fn make_light(param_set: &ParamSet, ro: &mut Box<RenderOptions>) {
    // MakeLight (api.cpp:591)
    if param_set.name == String::from("point") {
        let i: Spectrum = param_set
            .find_one_spectrum(String::from("I"), Spectrum::new(1.0 as Float));
        // Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
        // Point3f P = paramSet.FindOnePoint3f("from", Point3f(0, 0, 0));
        // Transform l2w = Translate(Vector3f(P.x, P.y, P.z)) * light2world;
        // return std::make_shared<PointLight>(l2w, medium, I * sc);
        unsafe {
            let point_light = Arc::new(PointLight::new(&CUR_TRANSFORM.t[0], &i));
            ro.lights.push(point_light);
        }
    } else if param_set.name == String::from("spot") {
        println!("TODO: CreateSpotLight");
    } else if param_set.name == String::from("goniometric") {
        println!("TODO: CreateGoniometricLight");
    } else if param_set.name == String::from("projection") {
        println!("TODO: CreateProjectionLight");
    } else if param_set.name == String::from("distant") {
        // CreateDistantLight
        let l: Spectrum = param_set
            .find_one_spectrum(String::from("L"), Spectrum::new(1.0 as Float));
        let sc: Spectrum =
            param_set.find_one_spectrum(String::from("scale"), Spectrum::new(1.0 as Float));
        let from: Point3f = param_set.find_one_point3f(String::from("from"),
                                                       Point3f {
                                                           x: 0.0,
                                                           y: 0.0,
                                                           z: 0.0,
                                                       });
        let to: Point3f = param_set.find_one_point3f(String::from("to"),
                                                     Point3f {
                                                         x: 0.0,
                                                         y: 0.0,
                                                         z: 0.0,
                                                     });
        let dir: Vector3f = from - to;
        // return std::make_shared<DistantLight>(light2world, L * sc, dir);
        unsafe {
            let distant_light = Arc::new(DistantLight::new(&CUR_TRANSFORM.t[0], &(l * sc), &dir));
            ro.lights.push(distant_light);
        }
    } else if param_set.name == String::from("infinite") {
        println!("TODO: CreateInfiniteLight");
    } else if param_set.name == String::from("exinfinite") {
        println!("TODO: CreateInfiniteLight");
    } else {
        panic!("MakeLight: unknown name {}", param_set.name);
    }
}

fn make_texture(param_set: &ParamSet) {
    unsafe {
        // pbrtTexture (api.cpp:1049)
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
                                param_set.name.clone(), st);
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
        // MakeFloatTexture(texname, curTransform[0], tp);
        // or
        // MakeSpectrumTexture(texname, curTransform[0], tp);
    }
}

fn pbrt_float_parameter<R, I>(pairs: &mut pest::iterators::Pairs<R, I>) -> (String, Vec<Float>)
    where I: pest::inputs::Input, R: pest::RuleType
{
    let mut floats: Vec<Float> = Vec::new();
    // single float or several floats using brackets
    let ident = pairs.next();
    let string: String =
        String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
    let _lbrack = pairs.next(); // assume opening bracket
    let mut number = pairs.next();
    while number.is_some() {
        let pair = number.unwrap().clone();
        if pair.as_str() == String::from("]") {
            // closing bracket found
            break;
        } else {
            let float: Float =
                f32::from_str(pair.into_span().as_str()).unwrap();
            floats.push(float);
        }
        number = pairs.next();
    }
    (string, floats)
}

fn pbrt_integer_parameter<R, I>(pairs: &mut pest::iterators::Pairs<R, I>) -> (String, Vec<i32>)
    where I: pest::inputs::Input, R: pest::RuleType
{
    let mut integers: Vec<i32> = Vec::new();
    // single integer with brackets
    let ident = pairs.next();
    let string: String =
        String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
    let _lbrack = pairs.next();
    let number = pairs.next();
    let number: i32 =
        i32::from_str(number.unwrap().clone().into_span().as_str()).unwrap();
    integers.push(number);
    (string, integers)
}

fn pbrt_string_parameter<R, I>(pairs: &mut pest::iterators::Pairs<R, I>) -> (String, String)
    where I: pest::inputs::Input, R: pest::RuleType
{
    // single string without brackets
    let ident = pairs.next();
    let string1: String =
        String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
    let option = pairs.next();
    let lbrack = option.clone().unwrap();
    let string2: String;
    if lbrack.as_str() == String::from("[") {
        let string = pairs.next();
        let pair = string.unwrap().clone();
        let ident = pair.into_inner().next();
        string2 = String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
    } else {
        let string = option.clone();
        let pair = string.unwrap().clone();
        let ident = pair.into_inner().next();
        string2 = String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
    }
    (string1, string2)
}

fn pbrt_shape(param_set: &ParamSet)
              -> (Vec<Arc<Shape + Send + Sync>>, Vec<Arc<Material + Send + Sync>>) {
    let mut shapes: Vec<Arc<Shape + Send + Sync>> = Vec::new();
    let mut materials: Vec<Arc<Material + Send + Sync>> = Vec::new();
    // pbrtShape (api.cpp:1153)
    // TODO: if (!curTransform.IsAnimated()) { ... }
    // TODO: transformCache.Lookup(curTransform[0], &ObjToWorld, &WorldToObj);
    unsafe {
        let obj_to_world: Transform = Transform {
            m: CUR_TRANSFORM.t[0].m,
            m_inv: CUR_TRANSFORM.t[0].m_inv,
        };
        let world_to_obj: Transform = Transform {
            m: CUR_TRANSFORM.t[0].m_inv,
            m_inv: CUR_TRANSFORM.t[0].m,
        };
        // MakeShapes (api.cpp:296)
        if param_set.name == String::from("sphere") {
            // CreateSphereShape
            let radius: Float = param_set.find_one_float(String::from("radius"), 1.0 as Float);
            let z_min: Float = param_set.find_one_float(String::from("zmin"), -radius);
            let z_max: Float = param_set.find_one_float(String::from("zmax"), radius);
            let phi_max: Float = param_set.find_one_float(String::from("phimax"), 360.0 as Float);
            let sphere = Arc::new(Sphere::new(obj_to_world,
                                              world_to_obj,
                                              false,
                                              false,
                                              radius,
                                              z_min,
                                              z_max,
                                              phi_max));
            let mtl: Arc<Material + Send + Sync> = create_material();
            shapes.push(sphere.clone());
            materials.push(mtl.clone());
        } else if param_set.name == String::from("cylinder") {
            let radius: Float = param_set.find_one_float(String::from("radius"), 1.0);
            let z_min: Float = param_set.find_one_float(String::from("zmin"), -radius);
            let z_max: Float = param_set.find_one_float(String::from("zmax"), radius);
            let phi_max: Float = param_set.find_one_float(String::from("phimax"), 360.0 as Float);
            let cylinder = Arc::new(Cylinder::new(obj_to_world,
                                                  world_to_obj,
                                                  false,
                                                  radius,
                                                  z_min,
                                                  z_max,
                                                  phi_max));
            let mtl: Arc<Material + Send + Sync> = create_material();
            shapes.push(cylinder.clone());
            materials.push(mtl.clone());
        } else if param_set.name == String::from("disk") {
            let height: Float = param_set.find_one_float(String::from("height"), 0.0);
            let radius: Float = param_set.find_one_float(String::from("radius"), 1.0);
            let inner_radius: Float = param_set.find_one_float(String::from("innerradius"), 0.0);
            let phi_max: Float = param_set.find_one_float(String::from("phimax"), 360.0);
            let disk = Arc::new(Disk::new(obj_to_world,
                                          world_to_obj,
                                          false,
                                          false,
                                          height,
                                          radius,
                                          inner_radius,
                                          phi_max));
            let mtl: Arc<Material + Send + Sync> = create_material();
            shapes.push(disk.clone());
            materials.push(mtl.clone());
        } else if param_set.name == String::from("cone") {
            println!("TODO: CreateConeShape");
        } else if param_set.name == String::from("paraboloid") {
            println!("TODO: CreateParaboloidShape");
        } else if param_set.name == String::from("hyperboloid") {
            println!("TODO: CreateHyperboloidShape");
        } else if param_set.name == String::from("curve") {
            println!("TODO: CreateCurveShape");
        } else if param_set.name == String::from("trianglemesh") {
            let vi = param_set.find_int(String::from("indices"));
            let p = param_set.find_point3f(String::from("P"));
            // try "uv" with Point2f
            let mut uvs = param_set.find_point2f(String::from("uv"));
            if uvs.is_empty() {
                // try "st" with Point2f
                uvs = param_set.find_point2f(String::from("st"));
            }
            if uvs.is_empty() {
                // try "uv" with float
                let mut fuv = param_set.find_float(String::from("uv"));
                if fuv.is_empty() {
                    // try "st" with float
                    fuv = param_set.find_float(String::from("st"));
                }
                if !fuv.is_empty() {
                    // found some float UVs
                    for i in 0..(fuv.len() / 2) {
                        uvs.push(Point2f {
                                     x: fuv[2 * i],
                                     y: fuv[2 * i + 1],
                                 });
                    }
                }
            }
            if !uvs.is_empty() {
                // TODO: if (nuvi < npi) {...} else if (nuvi > npi) ...
                assert!(uvs.len() == p.len());
            }
            assert!(vi.len() > 0_usize);
            assert!(p.len() > 0_usize);
            let s = param_set.find_vector3f(String::from("S"));
            if !s.is_empty() {
                assert!(s.len() == p.len());
            }
            let n = param_set.find_normal3f(String::from("N"));
            if !n.is_empty() {
                assert!(n.len() == p.len());
            }
            for i in 0..vi.len() {
                if vi[i] as usize >= p.len() {
                    panic!("trianglemesh has out of-bounds vertex index {} ({} \"P\" values were given)",
                           vi[i],
                           p.len());
                }
            }
            // TODO: alpha
            // CreateTriangleMesh
            // transform mesh vertices to world space
            let mut p_ws: Vec<Point3f> = Vec::new();
            let n_vertices: usize = p.len();
            for i in 0..n_vertices {
                p_ws.push(obj_to_world.transform_point(p[i]));
            }
            // vertex indices are expected as usize, not i32
            let mut vertex_indices: Vec<usize> = Vec::new();
            for i in 0..vi.len() {
                vertex_indices.push(vi[i] as usize);
            }
            let mesh = Arc::new(TriangleMesh::new(obj_to_world,
                                                  world_to_obj,
                                                  false, // reverse_orientation
                                                  false, // transform_swaps_handedness
                                                  vi.len() / 3, // n_triangles
                                                  vertex_indices,
                                                  n_vertices,
                                                  p_ws, // in world space
                                                  s,
                                                  n,
                                                  uvs));
            let mtl: Arc<Material + Send + Sync> = create_material();
            for id in 0..mesh.n_triangles {
                let triangle = Arc::new(Triangle::new(mesh.object_to_world,
                                                      mesh.world_to_object,
                                                      mesh.transform_swaps_handedness,
                                                      mesh.clone(),
                                                      id));
                shapes.push(triangle.clone());
                materials.push(mtl.clone());
            }
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
    (shapes, materials)
}

fn pbrt_world_end() {
    // println!("WorldEnd");
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
                        // println!("TODO: CreateGaussianFilter");
                        let xw: Float = ro.filter_params.find_one_float(String::from("xwidth"),
                                                                        2.0);
                        let yw: Float = ro.filter_params.find_one_float(String::from("ywidth"),
                                                                        2.0);
                        let alpha: Float = ro.filter_params.find_one_float(String::from("alpha"),
                                                                           2.0);
                        // see gaussian.h (GaussianFilter constructor)
                        let exp_x: Float = (-alpha * xw * xw).exp();
                        let exp_y: Float = (-alpha * yw * yw).exp();
                        let gaussian_filter: Arc<Filter + Sync + Send> =
                            Arc::new(GaussianFilter {
                                         alpha: alpha,
                                         exp_x: exp_x,
                                         exp_y: exp_y,
                                         radius: Vector2f { x: xw, y: yw },
                                         inv_radius: Vector2f {
                                             x: 1.0 / xw,
                                             y: 1.0 / yw,
                                         },
                                     });
                        some_filter = Some(gaussian_filter);
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
                        let filename: String =
                            ro.film_params
                                .find_one_string(String::from("filename"), String::new());
                        let xres: i32 = ro.film_params
                            .find_one_int(String::from("xresolution"), 1280);
                        let yres: i32 = ro.film_params
                            .find_one_int(String::from("yresolution"), 720);
                        // TODO: if (PbrtOptions.quickRender) xres = std::max(1, xres / 4);
                        // TODO: if (PbrtOptions.quickRender) yres = std::max(1, yres / 4);
                        let crop: Bounds2f = Bounds2f {
                            p_min: Point2f { x: 0.0, y: 0.0 },
                            p_max: Point2f { x: 1.0, y: 1.0 },
                        };
                        // TODO: const Float *cr = params.FindFloat("cropwindow", &cwi);
                        let scale: Float = ro.film_params.find_one_float(String::from("scale"),
                                                                         1.0);
                        let diagonal: Float = ro.film_params
                            .find_one_float(String::from("diagonal"), 35.0);
                        let max_sample_luminance: Float =
                            ro.film_params
                                .find_one_float(String::from("maxsampleluminance"),
                                                std::f32::INFINITY);
                        if let Some(filter) = some_filter {
                            let film: Arc<Film> = Arc::new(Film::new(Point2i { x: xres, y: yres },
                                                                     crop,
                                                                     filter,
                                                                     diagonal,
                                                                     filename,
                                                                     scale,
                                                                     max_sample_luminance));
                            // MakeCamera
                            // TODO: let mut some_camera: Option<Arc<Camera + Sync + Send>> = None;
                            let mut some_camera: Option<Arc<PerspectiveCamera>> = None;
                            // TODO: MediumInterface mediumInterface = graphicsState.CreateMediumInterface();
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
                                // TODO: let perspective_camera: Arc<Camera + Sync + Send> =
                                let perspective_camera: Arc<PerspectiveCamera> =
                                    Arc::new(PerspectiveCamera::new(animated_cam_to_world,
                                                                    screen,
                                                                    shutteropen,
                                                                    shutterclose,
                                                                    lensradius,
                                                                    focaldistance,
                                                                    fov,
                                                                    film));
                                some_camera = Some(perspective_camera);
                            } else if ro.camera_name == String::from("orthographic") {
                                println!("TODO: CreateOrthographicCamera");
                            } else if ro.camera_name == String::from("realistic") {
                                println!("TODO: CreateRealisticCamera");
                            } else if ro.camera_name == String::from("environment") {
                                println!("TODO: CreateEnvironmentCamera");
                            } else {
                                panic!("Camera \"{}\" unknown.", ro.camera_name);
                            }
                            if let Some(camera) = some_camera {
                                // MakeSampler
                                let mut sampler: ZeroTwoSequenceSampler = ZeroTwoSequenceSampler::default();
                                if ro.sampler_name == String::from("lowdiscrepancy") ||
                                   ro.sampler_name == String::from("02sequence") {
                                    let nsamp: i32 =
                                        ro.sampler_params
                                            .find_one_int(String::from("pixelsamples"), 16);
                                    let sd: i32 = ro.sampler_params
                                        .find_one_int(String::from("dimensions"), 4);
                                    // TODO: if (PbrtOptions.quickRender) nsamp = 1;
                                    let new_sampler = ZeroTwoSequenceSampler::new(nsamp as i64,
                                                                                  sd as i64);
                                    sampler = new_sampler.clone()
                                } else if ro.sampler_name == String::from("maxmindist") {
                                    println!("TODO: CreateMaxMinDistSampler");
                                } else if ro.sampler_name == String::from("halton") {
                                    // println!("TODO: CreateHaltonSampler");
                                    // int nsamp = params.FindOneInt("pixelsamples", 16);
                                    // if (PbrtOptions.quickRender) nsamp = 1;
                                    // bool sampleAtCenter = params.FindOneBool("samplepixelcenter", false);
                                    // return new HaltonSampler(nsamp, sampleBounds, sampleAtCenter);
                                    // WARNING: Use ZeroTwoSequenceSampler for now !!!
                                    let nsamp: i64 = 16;
                                    let sd: i64 = 4;
                                    let new_sampler = ZeroTwoSequenceSampler::new(nsamp, sd);
                                    sampler = new_sampler.clone()
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
                                // if let Some(mut sampler) = some_sampler {
                                let mut some_integrator: Option<Arc<SamplerIntegrator + Sync + Send>> = None;
                                if ro.integrator_name == String::from("whitted") {
                                    println!("TODO: CreateWhittedIntegrator");
                                } else if ro.integrator_name == String::from("directlighting") {
                                    // CreateDirectLightingIntegrator
                                    let max_depth: i32 =
                                        ro.integrator_params
                                            .find_one_int(String::from("maxdepth"), 5);
                                    let st: String = ro.integrator_params
                                        .find_one_string(String::from("strategy"),
                                                         String::from("all"));
                                    let strategy: LightStrategy;
                                    if st == String::from("one") {
                                        strategy = LightStrategy::UniformSampleOne;
                                    } else if st == String::from("all") {
                                        strategy = LightStrategy::UniformSampleAll;
                                    } else {
                                        panic!("Strategy \"{}\" for direct lighting unknown.", st);
                                    }
                                    // TODO: const int *pb = params.FindInt("pixelbounds", &np);
                                    let pixel_bounds: Bounds2i = Bounds2i {
                                        p_min: Point2i { x: 0, y: 0 },
                                        p_max: Point2i { x: xres, y: yres },
                                    };
                                    let integrator =
                                        Arc::new(DirectLightingIntegrator::new(strategy,
                                                                               max_depth as i64,
                                                                               pixel_bounds));
                                    some_integrator = Some(integrator);
                                } else if ro.integrator_name == String::from("path") {
                                    // CreatePathIntegrator
                                    let max_depth: i32 =
                                        ro.integrator_params
                                            .find_one_int(String::from("maxdepth"), 5);
                                    let pb: Vec<i32> = ro.integrator_params
                                        .find_int(String::from("pixelbounds"));
                                    let np: usize = pb.len();
                                    let pixel_bounds: Bounds2i = camera.film.get_sample_bounds();
                                    if np > 0 as usize {
                                        if np != 4 as usize {
                                            panic!("Expected four values for \"pixelbounds\" parameter. Got {}.",
                                                   np);
                                        } else {
                                            println!("TODO: pixelBounds = Intersect(...)");
                                            // pixelBounds = Intersect(pixelBounds,
                                            //                         Bounds2i{{pb[0], pb[2]}, {pb[1], pb[3]}});
                                            // if (pixelBounds.Area() == 0)
                                            //     Error("Degenerate \"pixelbounds\" specified.");
                                        }
                                    }
                                    let rr_threshold: Float =
                                        ro.integrator_params
                                            .find_one_float(String::from("rrthreshold"),
                                                            1.0 as Float);
                                    println!("DEBUG: rr_threshold = {:?}", rr_threshold);
                                    // std::string lightStrategy =
                                    //     params.FindOneString("lightsamplestrategy", "spatial");
                                    let light_strategy: String =
                                        ro.integrator_params
                                            .find_one_string(String::from("lightsamplestrategy"),
                                                             String::from("spatial"));
                                    let integrator = Arc::new(PathIntegrator::new(max_depth as
                                                                                  u32,
                                                                                  &camera,
                                                                                  &sampler,
                                                                                  pixel_bounds,
                                                                                  rr_threshold,
                                                                                  light_strategy));
                                    some_integrator = Some(integrator);
                                } else if ro.integrator_name == String::from("volpath") {
                                    println!("TODO: CreateVolPathIntegrator");
                                } else if ro.integrator_name == String::from("bdpt") {
                                    println!("TODO: CreateBDPTIntegrator");
                                } else if ro.integrator_name == String::from("mlt") {
                                    println!("TODO: CreateMLTIntegrator");
                                } else if ro.integrator_name == String::from("ambientocclusion") {
                                    // CreateAOIntegrator
                                    let pb: Vec<i32> = ro.integrator_params
                                        .find_int(String::from("pixelbounds"));
                                    let np: usize = pb.len();
                                    let pixel_bounds: Bounds2i = camera.film.get_sample_bounds();
                                    if np > 0 as usize {
                                        if np != 4 as usize {
                                            panic!("Expected four values for \"pixelbounds\" parameter. Got {}.",
                                                   np);
                                        } else {
                                            println!("TODO: pixelBounds = Intersect(...)");
                                            // pixelBounds = Intersect(pixelBounds,
                                            //                         Bounds2i{{pb[0], pb[2]}, {pb[1], pb[3]}});
                                            // if (pixelBounds.Area() == 0)
                                            //     Error("Degenerate \"pixelbounds\" specified.");
                                        }
                                    }
                                    let rr_threshold: Float =
                                        ro.integrator_params
                                            .find_one_float(String::from("rrthreshold"),
                                                            1.0 as Float);
                                    println!("DEBUG: rr_threshold = {:?}", rr_threshold);
                                    let cos_sample: bool =
                                        ro.integrator_params
                                            .find_one_bool(String::from("cossample"), true);
                                    println!("DEBUG: cos_sample = {:?}", cos_sample);
                                    // int nSamples = params.Find_One_Int("nsamples", 64);
                                    let n_samples: i32 =
                                        ro.integrator_params
                                            .find_one_int(String::from("nsamples"), 64 as i32);
                                    // return new AOIntegrator(cosSample, nSamples, camera, sampler, pixelBounds);

                                    let integrator = Arc::new(AOIntegrator::new(cos_sample,
                                                                                n_samples,
                                                                                &camera,
                                                                                &sampler,
                                                                                pixel_bounds));
                                    some_integrator = Some(integrator);
                                } else if ro.integrator_name == String::from("sppm") {
                                    println!("TODO: CreateSPPMIntegrator");
                                } else {
                                    panic!("Integrator \"{}\" unknown.", ro.integrator_name);
                                }
                                if let Some(mut integrator) = some_integrator {
                                    // MakeIntegrator
                                    // TODO: if (renderOptions->haveScatteringMedia && ...)
                                    if ro.lights.is_empty() {
                                        // warn if no light sources are defined
                                        println!("WARNING: No light sources defined in scene; rendering a black image.",);
                                    }
                                    // MakeAccelerator
                                    if ro.accelerator_name == String::from("bvh") {
                                        //  CreateBVHAccelerator
                                        let split_method_name: String =
                                                ro.accelerator_params.find_one_string(String::from("splitmethod"),
                                                                                      String::from("sah"));
                                        let split_method;
                                        if split_method_name == String::from("sah") {
                                            split_method = SplitMethod::SAH;
                                        } else if split_method_name == String::from("hlbvh") {
                                            split_method = SplitMethod::HLBVH;
                                        } else if split_method_name == String::from("middle") {
                                            split_method = SplitMethod::Middle;
                                        } else if split_method_name == String::from("equal") {
                                            split_method = SplitMethod::EqualCounts;
                                        } else {
                                            println!("WARNING: BVH split method \"{}\" unknown.  Using \"sah\".",
                                                     split_method_name);
                                            split_method = SplitMethod::SAH;
                                        }
                                        let max_prims_in_node: i32 =
                                                ro.accelerator_params.find_one_int(String::from("maxnodeprims"), 4);
                                        let accelerator =
                                            Arc::new(BVHAccel::new(ro.primitives.clone(),
                                                                   max_prims_in_node as usize,
                                                                   split_method));
                                        // MakeScene
                                        let scene: Scene = Scene::new(accelerator.clone(),
                                                                      ro.lights.clone());
                                        // TODO: primitives.erase(primitives.begin(), primitives.end());
                                        // TODO: lights.erase(lights.begin(), lights.end());
                                        pbrt::render(&scene, camera, &mut sampler, &mut integrator);
                                    } else if ro.accelerator_name == String::from("kdtree") {
                                        // println!("TODO: CreateKdTreeAccelerator");
                                        // WARNING: Use BVHAccel for now !!!
                                        let accelerator = Arc::new(BVHAccel::new(ro.primitives
                                                                                     .clone(),
                                                                                 4,
                                                                                 SplitMethod::SAH));
                                        // MakeScene
                                        let scene: Scene = Scene::new(accelerator.clone(),
                                                                      ro.lights.clone());
                                        // TODO: primitives.erase(primitives.begin(), primitives.end());
                                        // TODO: lights.erase(lights.begin(), lights.end());
                                        pbrt::render(&scene, camera, &mut sampler, &mut integrator);
                                    } else {
                                        panic!("Accelerator \"{}\" unknown.", ro.accelerator_name);
                                    }
                                } else {
                                    panic!("Unable to create integrator.");
                                }
                                // } else {
                                //     panic!("Unable to create sampler.");
                                // }
                            } else {
                                panic!("Unable to create camera.");
                            }
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
    let child = thread::Builder::new()
        .stack_size(32 * 1024 * 1024)
        .spawn(move || {
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
                        let num_bytes = reader.read_to_string(&mut str_buf);
                        if num_bytes.is_ok() {
                            let n_bytes = num_bytes.unwrap();
                            println!("{} bytes read", n_bytes);
                        }
                        unsafe {
                            // render options
                            NAMED_COORDINATE_SYSTEMS = Some(Box::new(HashMap::new()));
                            RENDER_OPTIONS = Some(Box::new(RenderOptions::default()));
                            GRAPHICS_STATE = Some(Box::new(GraphicsState::default()));
                            PUSHED_GRAPHICS_STATES = Some(Box::new(Vec::new()));
                            PUSHED_TRANSFORMS = Some(Box::new(Vec::new()));
                            PARAM_SET = Some(Box::new(ParamSet::default()));
                            // parser
                            let pairs = PbrtParser::parse_str(Rule::pbrt, &str_buf).unwrap_or_else(|e| panic!("{}", e));
                            // assert!(parser.pbrt());
                            // assert!(parser.end());
                            // println!("{:?}", parser.queue());
                            println!("do something with created tokens ...");
                            for pair in pairs {
                                // a pair is a combination of the rule which matched and a span of input
                                match pair.as_rule() {
                                    Rule::statement => {
                                        for statement_pair in pair.into_inner() {
                                            match statement_pair.as_rule() {
                                                Rule::concat_transform => println!("TODO: Rule::concat_transform"),
                                                Rule::keyword => {
                                                    for keyword_pair in statement_pair.into_inner() {
                                                        match keyword_pair.as_rule() {
                                                            Rule::attribute_begin => {
                                                                if let Some(ref mut graphics_state) = GRAPHICS_STATE {
                                                                    if let Some(ref mut pushed_graphics_states) = PUSHED_GRAPHICS_STATES {
                                                                        let mut material_param_set: ParamSet = ParamSet::default();
                                                                        material_param_set.copy_from(&graphics_state.material_params);
                                                                        let mut area_light_param_set: ParamSet = ParamSet::default();
                                                                        area_light_param_set.copy_from(&graphics_state.area_light_params);
                                                                        pushed_graphics_states.push(GraphicsState {
                                                                            float_textures: graphics_state.float_textures.clone(),
                                                                            spectrum_textures: graphics_state.spectrum_textures.clone(),
                                                                            material_params: material_param_set,
                                                                            material: String::from(graphics_state.material.as_ref()),
                                                                            named_materials: graphics_state.named_materials.clone(),
                                                                            current_named_material: String::from(graphics_state.current_named_material.as_ref()),
                                                                            area_light_params: area_light_param_set,
                                                                            area_light: String::from(graphics_state.area_light.as_ref()),
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
                                                            },
                                                            Rule::attribute_end => {
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
                                                                        // area_light_params
                                                                        graphics_state.area_light_params.reset(String::new(),
                                                                                                               String::from(""),
                                                                                                               String::from(""),
                                                                                                               String::new());
                                                                        graphics_state.area_light_params.copy_from(&pgs.area_light_params);
                                                                        // area_light
                                                                        graphics_state.area_light = String::from(pgs.area_light.as_ref());
                                                                    }
                                                                    if let Some(ref mut pt) = PUSHED_TRANSFORMS {
                                                                        let popped_transform_set: TransformSet = pt.pop().unwrap();
                                                                        CUR_TRANSFORM.t[0] = popped_transform_set.t[0];
                                                                        CUR_TRANSFORM.t[1] = popped_transform_set.t[1];
                                                                    }
                                                                    // TODO? pushedActiveTransformBits.push_back(activeTransformBits);
                                                                }
                                                            },
                                                            Rule::world_begin => {
                                                                println!("WorldBegin");
                                                                CUR_TRANSFORM.t[0] = Transform::default();
                                                                CUR_TRANSFORM.t[1] = Transform::default();
                                                                if let Some(ref mut named_coordinate_systems) = NAMED_COORDINATE_SYSTEMS {
                                                                    named_coordinate_systems.insert("world",
                                                                                                    TransformSet {
                                                                                                        t: [Transform::default(); 2]
                                                                                                    });
                                                                }
                                                            },
                                                            _ => unreachable!()
                                                        }
                                                    }
                                                },
                                                Rule::look_at => {
                                                    let mut numbers: Vec<Float> = Vec::new();
                                                    for look_at_pair in statement_pair.into_inner() {
                                                        let number: Float =
                                                            f32::from_str(look_at_pair.clone().into_span().as_str()).unwrap();
                                                        numbers.push(number);
                                                    }
                                                    let pos: Point3f = Point3f { x: numbers[0], y: numbers[1], z: numbers[2], };
                                                    let look: Point3f = Point3f { x: numbers[3], y: numbers[4], z: numbers[5], };
                                                    let up: Vector3f = Vector3f { x: numbers[6], y: numbers[7], z: numbers[8], };
                                                    let look_at: Transform = Transform::look_at(pos, look, up);
                                                    CUR_TRANSFORM.t[0] = CUR_TRANSFORM.t[0] * look_at;
                                                    CUR_TRANSFORM.t[1] = CUR_TRANSFORM.t[1] * look_at;
                                                },
                                                Rule::named_statement => {
                                                    for named_statement_pair in statement_pair.into_inner() {
                                                        match named_statement_pair.as_rule() {
                                                            Rule::accelerator => println!("TODO: Rule::accelerator"),
                                                            Rule::area_light_source => println!("TODO: Rule::area_light_source"),
                                                            Rule::camera => {
                                                                for camera_pair in named_statement_pair.into_inner() {
                                                                    match camera_pair.as_rule() {
                                                                        Rule::string => {
                                                                            let mut string_pairs = camera_pair.into_inner();
                                                                            let ident = string_pairs.next();
                                                                            let name: String =
                                                                                String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
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
                                                                                }
                                                                            } else {
                                                                                panic!("Can't get render options.");
                                                                            }
                                                                            if let Some(ref mut param_set) = PARAM_SET {
                                                                                param_set.reset(String::from("Camera"),
                                                                                                String::from(""),
                                                                                                String::from(""),
                                                                                                String::from(""));
                                                                            } else {
                                                                                panic!("Can't get parameter set.");
                                                                            }
                                                                        },
                                                                        Rule::parameter => {
                                                                            for parameter_pair in camera_pair.into_inner() {
                                                                                match parameter_pair.as_rule() {
                                                                                    Rule::bool_param => println!("TODO: Rule::bool_param"),
                                                                                    Rule::float_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if floats.len() == 1 {
                                                                                                param_set.add_float(string, floats[0]);
                                                                                            } else {
                                                                                                param_set.add_floats(string, floats);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::integer_param => {
                                                                                        let tuple: (String, Vec<i32>) = pbrt_integer_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let integers: Vec<i32> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if integers.len() == 1 {
                                                                                                param_set.add_int(string, integers[0]);
                                                                                            } else {
                                                                                                param_set.add_ints(string, integers);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::point_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if floats.len() == 3 {
                                                                                                param_set.add_point3f(string,
                                                                                                                      Point3f {
                                                                                                                          x: floats[0],
                                                                                                                          y: floats[1],
                                                                                                                          z: floats[2],
                                                                                                                      });
                                                                                            } else {
                                                                                                param_set.add_point3fs(string, floats);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::normal_param => println!("TODO: Rule::normal_param"),
                                                                                    Rule::rgb_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            param_set.add_rgb_spectrum(string,
                                                                                                                       Spectrum {
                                                                                                                           c: [floats[0], floats[1], floats[2]],
                                                                                                                       });
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::spectrum_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            param_set.add_rgb_spectrum(string,
                                                                                                                       Spectrum {
                                                                                                                           c: [floats[0], floats[1], floats[2]],
                                                                                                                       });
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::string_param => {
                                                                                        let tuple: (String, String) = pbrt_string_parameter(&mut parameter_pair.into_inner());
                                                                                        let string1: String = tuple.0;
                                                                                        let string2: String = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                                param_set.add_string(string1, string2);
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::texture_param => println!("TODO: Rule::texture_param"),
                                                                                    Rule::vector_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if floats.len() == 3 {
                                                                                                param_set.add_vector3f(string,
                                                                                                                       Vector3f {
                                                                                                                           x: floats[0],
                                                                                                                           y: floats[1],
                                                                                                                           z: floats[2],
                                                                                                                       });
                                                                                            } else {
                                                                                                param_set.add_point3fs(string, floats);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    _ => unreachable!()
                                                                                };
                                                                            }
                                                                        },
                                                                        _ => unreachable!()
                                                                    }
                                                                }
                                                                // we should have the camera parameters by now
                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                    if let Some(ref mut ro) = RENDER_OPTIONS {
                                                                        println!("Camera \"{}\" ", ro.camera_name);
                                                                        ro.camera_params.copy_from(param_set);
                                                                        print_params(&ro.camera_params);
                                                                    } else {
                                                                        panic!("Can't get render options.");
                                                                    }
                                                                } else {
                                                                    panic!("Can't get parameter set.");
                                                                }
                                                            },
                                                            Rule::coord_sys_transform => println!("TODO: Rule::coord_sys_transform"),
                                                            Rule::film => {
                                                                for film_pair in named_statement_pair.into_inner() {
                                                                    match film_pair.as_rule() {
                                                                        Rule::string => {
                                                                            let mut string_pairs = film_pair.into_inner();
                                                                            let ident = string_pairs.next();
                                                                            let name: String =
                                                                                String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
                                                                            if let Some(ref mut ro) = RENDER_OPTIONS {
                                                                                ro.film_name = name;
                                                                            } else {
                                                                                panic!("Can't get render options.");
                                                                            }
                                                                            if let Some(ref mut param_set) = PARAM_SET {
                                                                                param_set.reset(String::from("Film"),
                                                                                                String::from(""),
                                                                                                String::from(""),
                                                                                                String::from(""));
                                                                            } else {
                                                                                panic!("Can't get parameter set.");
                                                                            }
                                                                        },
                                                                        Rule::parameter => {
                                                                            for parameter_pair in film_pair.into_inner() {
                                                                                match parameter_pair.as_rule() {
                                                                                    Rule::bool_param => println!("TODO: Rule::bool_param"),
                                                                                    Rule::float_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if floats.len() == 1 {
                                                                                                param_set.add_float(string, floats[0]);
                                                                                            } else {
                                                                                                param_set.add_floats(string, floats);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::integer_param => {
                                                                                        let tuple: (String, Vec<i32>) = pbrt_integer_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let integers: Vec<i32> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if integers.len() == 1 {
                                                                                                param_set.add_int(string, integers[0]);
                                                                                            } else {
                                                                                                param_set.add_ints(string, integers);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::point_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if floats.len() == 3 {
                                                                                                param_set.add_point3f(string,
                                                                                                                      Point3f {
                                                                                                                          x: floats[0],
                                                                                                                          y: floats[1],
                                                                                                                          z: floats[2],
                                                                                                                      });
                                                                                            } else {
                                                                                                param_set.add_point3fs(string, floats);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::normal_param => println!("TODO: Rule::normal_param"),
                                                                                    Rule::rgb_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            param_set.add_rgb_spectrum(string,
                                                                                                                       Spectrum {
                                                                                                                           c: [floats[0], floats[1], floats[2]],
                                                                                                                       });
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::spectrum_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            param_set.add_rgb_spectrum(string,
                                                                                                                       Spectrum {
                                                                                                                           c: [floats[0], floats[1], floats[2]],
                                                                                                                       });
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::string_param => {
                                                                                        let tuple: (String, String) = pbrt_string_parameter(&mut parameter_pair.into_inner());
                                                                                        let string1: String = tuple.0;
                                                                                        let string2: String = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                                param_set.add_string(string1, string2);
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::texture_param => println!("TODO: Rule::texture_param"),
                                                                                    Rule::vector_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if floats.len() == 3 {
                                                                                                param_set.add_vector3f(string,
                                                                                                                       Vector3f {
                                                                                                                           x: floats[0],
                                                                                                                           y: floats[1],
                                                                                                                           z: floats[2],
                                                                                                                       });
                                                                                            } else {
                                                                                                param_set.add_point3fs(string, floats);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    _ => unreachable!()
                                                                                };
                                                                            }
                                                                        },
                                                                        _ => unreachable!()
                                                                    }
                                                                }
                                                                // we should have the film parameters by now
                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                    if let Some(ref mut ro) = RENDER_OPTIONS {
                                                                        println!("Film \"{}\" ", ro.film_name);
                                                                        ro.film_params.copy_from(param_set);
                                                                        print_params(&ro.film_params);
                                                                    } else {
                                                                        panic!("Can't get render options.");
                                                                    }
                                                                } else {
                                                                    panic!("Can't get parameter set.");
                                                                }
                                                            },
                                                            Rule::integrator => {
                                                                for integrator_pair in named_statement_pair.into_inner() {
                                                                    match integrator_pair.as_rule() {
                                                                        Rule::string => {
                                                                            let mut string_pairs = integrator_pair.into_inner();
                                                                            let ident = string_pairs.next();
                                                                            let name: String =
                                                                                String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
                                                                            if let Some(ref mut ro) = RENDER_OPTIONS {
                                                                                ro.integrator_name = name;
                                                                            } else {
                                                                                panic!("Can't get render options.");
                                                                            }
                                                                            if let Some(ref mut param_set) = PARAM_SET {
                                                                                param_set.reset(String::from("Integrator"),
                                                                                                String::from(""),
                                                                                                String::from(""),
                                                                                                String::from(""));
                                                                            } else {
                                                                                panic!("Can't get parameter set.");
                                                                            }
                                                                        },
                                                                        Rule::parameter => {
                                                                            for parameter_pair in integrator_pair.into_inner() {
                                                                                match parameter_pair.as_rule() {
                                                                                    Rule::bool_param => println!("TODO: Rule::bool_param"),
                                                                                    Rule::float_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if floats.len() == 1 {
                                                                                                param_set.add_float(string, floats[0]);
                                                                                            } else {
                                                                                                param_set.add_floats(string, floats);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::integer_param => {
                                                                                        let tuple: (String, Vec<i32>) = pbrt_integer_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let integers: Vec<i32> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if integers.len() == 1 {
                                                                                                param_set.add_int(string, integers[0]);
                                                                                            } else {
                                                                                                param_set.add_ints(string, integers);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::point_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if floats.len() == 3 {
                                                                                                param_set.add_point3f(string,
                                                                                                                      Point3f {
                                                                                                                          x: floats[0],
                                                                                                                          y: floats[1],
                                                                                                                          z: floats[2],
                                                                                                                      });
                                                                                            } else {
                                                                                                param_set.add_point3fs(string, floats);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::normal_param => println!("TODO: Rule::normal_param"),
                                                                                    Rule::rgb_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            param_set.add_rgb_spectrum(string,
                                                                                                                       Spectrum {
                                                                                                                           c: [floats[0], floats[1], floats[2]],
                                                                                                                       });
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::spectrum_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            param_set.add_rgb_spectrum(string,
                                                                                                                       Spectrum {
                                                                                                                           c: [floats[0], floats[1], floats[2]],
                                                                                                                       });
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::string_param => {
                                                                                        let tuple: (String, String) = pbrt_string_parameter(&mut parameter_pair.into_inner());
                                                                                        let string1: String = tuple.0;
                                                                                        let string2: String = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                                param_set.add_string(string1, string2);
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::texture_param => println!("TODO: Rule::texture_param"),
                                                                                    Rule::vector_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if floats.len() == 3 {
                                                                                                param_set.add_vector3f(string,
                                                                                                                       Vector3f {
                                                                                                                           x: floats[0],
                                                                                                                           y: floats[1],
                                                                                                                           z: floats[2],
                                                                                                                       });
                                                                                            } else {
                                                                                                param_set.add_point3fs(string, floats);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    _ => unreachable!()
                                                                                };
                                                                            }
                                                                        },
                                                                        _ => unreachable!()
                                                                    }
                                                                }
                                                                // we should have the integrator parameters by now
                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                    if let Some(ref mut ro) = RENDER_OPTIONS {
                                                                        println!("Integrator \"{}\" ", ro.integrator_name);
                                                                        ro.integrator_params.copy_from(param_set);
                                                                        print_params(&ro.integrator_params);
                                                                    } else {
                                                                        panic!("Can't get render options.");
                                                                    }
                                                                } else {
                                                                    panic!("Can't get parameter set.");
                                                                }
                                                            },
                                                            Rule::light_source => {
                                                                for light_source_pair in named_statement_pair.into_inner() {
                                                                    match light_source_pair.as_rule() {
                                                                        Rule::string => {
                                                                            let mut string_pairs = light_source_pair.into_inner();
                                                                            let ident = string_pairs.next();
                                                                            let name: String =
                                                                                String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
                                                                            if let Some(ref mut param_set) = PARAM_SET {
                                                                                param_set.reset(String::from("Light_Source"),
                                                                                                String::from(name),
                                                                                                String::from(""),
                                                                                                String::from(""));
                                                                            } else {
                                                                                panic!("Can't get parameter set.");
                                                                            }
                                                                        },
                                                                        Rule::parameter => {
                                                                            for parameter_pair in light_source_pair.into_inner() {
                                                                                match parameter_pair.as_rule() {
                                                                                    Rule::bool_param => println!("TODO: Rule::bool_param"),
                                                                                    Rule::float_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if floats.len() == 1 {
                                                                                                param_set.add_float(string, floats[0]);
                                                                                            } else {
                                                                                                param_set.add_floats(string, floats);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::integer_param => {
                                                                                        let tuple: (String, Vec<i32>) = pbrt_integer_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let integers: Vec<i32> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if integers.len() == 1 {
                                                                                                param_set.add_int(string, integers[0]);
                                                                                            } else {
                                                                                                param_set.add_ints(string, integers);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::point_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if floats.len() == 3 {
                                                                                                param_set.add_point3f(string,
                                                                                                                      Point3f {
                                                                                                                          x: floats[0],
                                                                                                                          y: floats[1],
                                                                                                                          z: floats[2],
                                                                                                                      });
                                                                                            } else {
                                                                                                param_set.add_point3fs(string, floats);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::normal_param => println!("TODO: Rule::normal_param"),
                                                                                    Rule::rgb_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            param_set.add_rgb_spectrum(string,
                                                                                                                       Spectrum {
                                                                                                                           c: [floats[0], floats[1], floats[2]],
                                                                                                                       });
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::spectrum_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            param_set.add_rgb_spectrum(string,
                                                                                                                       Spectrum {
                                                                                                                           c: [floats[0], floats[1], floats[2]],
                                                                                                                       });
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::string_param => {
                                                                                        let tuple: (String, String) = pbrt_string_parameter(&mut parameter_pair.into_inner());
                                                                                        let string1: String = tuple.0;
                                                                                        let string2: String = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                                param_set.add_string(string1, string2);
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::texture_param => println!("TODO: Rule::texture_param"),
                                                                                    Rule::vector_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if floats.len() == 3 {
                                                                                                param_set.add_vector3f(string,
                                                                                                                       Vector3f {
                                                                                                                           x: floats[0],
                                                                                                                           y: floats[1],
                                                                                                                           z: floats[2],
                                                                                                                       });
                                                                                            } else {
                                                                                                param_set.add_point3fs(string, floats);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    _ => unreachable!()
                                                                                };
                                                                            }
                                                                        },
                                                                        _ => unreachable!()
                                                                    }
                                                                }
                                                                // we should have the light_source parameters by now
                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                    if let Some(ref mut ro) = RENDER_OPTIONS {
                                                                        println!("LightSource \"{}\" ", param_set.name);
                                                                        print_params(&param_set);
                                                                        make_light(&param_set, ro);
                                                                    } else {
                                                                        panic!("Can't get render options.");
                                                                    }
                                                                } else {
                                                                    panic!("Can't get parameter set.");
                                                                }
                                                            },
                                                            Rule::make_named_material => println!("TODO: Rule::make_named_material"),
                                                            Rule::material => println!("TODO: Rule::material"),
                                                            Rule::named_material => println!("TODO: Rule::named_material"),
                                                            Rule::pixel_filter => {
                                                                for pixel_filter_pair in named_statement_pair.into_inner() {
                                                                    match pixel_filter_pair.as_rule() {
                                                                        Rule::string => {
                                                                            let mut string_pairs = pixel_filter_pair.into_inner();
                                                                            let ident = string_pairs.next();
                                                                            let name: String =
                                                                                String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
                                                                            if let Some(ref mut ro) = RENDER_OPTIONS {
                                                                                ro.filter_name = name;
                                                                            } else {
                                                                                panic!("Can't get render options.");
                                                                            }
                                                                            if let Some(ref mut param_set) = PARAM_SET {
                                                                                param_set.reset(String::from("Pixel_Filter"),
                                                                                                String::from(""),
                                                                                                String::from(""),
                                                                                                String::from(""));
                                                                            } else {
                                                                                panic!("Can't get parameter set.");
                                                                            }
                                                                        },
                                                                        Rule::parameter => {
                                                                            for parameter_pair in pixel_filter_pair.into_inner() {
                                                                                match parameter_pair.as_rule() {
                                                                                    Rule::bool_param => println!("TODO: Rule::bool_param"),
                                                                                    Rule::float_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if floats.len() == 1 {
                                                                                                param_set.add_float(string, floats[0]);
                                                                                            } else {
                                                                                                param_set.add_floats(string, floats);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::integer_param => {
                                                                                        let tuple: (String, Vec<i32>) = pbrt_integer_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let integers: Vec<i32> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if integers.len() == 1 {
                                                                                                param_set.add_int(string, integers[0]);
                                                                                            } else {
                                                                                                param_set.add_ints(string, integers);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::point_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if floats.len() == 3 {
                                                                                                param_set.add_point3f(string,
                                                                                                                      Point3f {
                                                                                                                          x: floats[0],
                                                                                                                          y: floats[1],
                                                                                                                          z: floats[2],
                                                                                                                      });
                                                                                            } else {
                                                                                                param_set.add_point3fs(string, floats);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::normal_param => println!("TODO: Rule::normal_param"),
                                                                                    Rule::rgb_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            param_set.add_rgb_spectrum(string,
                                                                                                                       Spectrum {
                                                                                                                           c: [floats[0], floats[1], floats[2]],
                                                                                                                       });
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::spectrum_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            param_set.add_rgb_spectrum(string,
                                                                                                                       Spectrum {
                                                                                                                           c: [floats[0], floats[1], floats[2]],
                                                                                                                       });
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::string_param => {
                                                                                        let tuple: (String, String) = pbrt_string_parameter(&mut parameter_pair.into_inner());
                                                                                        let string1: String = tuple.0;
                                                                                        let string2: String = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                                param_set.add_string(string1, string2);
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::texture_param => println!("TODO: Rule::texture_param"),
                                                                                    Rule::vector_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if floats.len() == 3 {
                                                                                                param_set.add_vector3f(string,
                                                                                                                       Vector3f {
                                                                                                                           x: floats[0],
                                                                                                                           y: floats[1],
                                                                                                                           z: floats[2],
                                                                                                                       });
                                                                                            } else {
                                                                                                param_set.add_point3fs(string, floats);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    _ => unreachable!()
                                                                                };
                                                                            }
                                                                        },
                                                                        _ => unreachable!()
                                                                    }
                                                                }
                                                                // we should have the pixel_filter parameters by now
                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                    if let Some(ref mut ro) = RENDER_OPTIONS {
                                                                        println!("Pixel_Filter \"{}\" ", ro.filter_name);
                                                                        ro.filter_params.copy_from(param_set);
                                                                        print_params(&ro.filter_params);
                                                                    } else {
                                                                        panic!("Can't get render options.");
                                                                    }
                                                                } else {
                                                                    panic!("Can't get parameter set.");
                                                                }
                                                            },
                                                            Rule::sampler => {
                                                                for sampler_pair in named_statement_pair.into_inner() {
                                                                    match sampler_pair.as_rule() {
                                                                        Rule::string => {
                                                                            let mut string_pairs = sampler_pair.into_inner();
                                                                            let ident = string_pairs.next();
                                                                            let name: String =
                                                                                String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
                                                                            if let Some(ref mut ro) = RENDER_OPTIONS {
                                                                                ro.sampler_name = name;
                                                                            } else {
                                                                                panic!("Can't get render options.");
                                                                            }
                                                                            if let Some(ref mut param_set) = PARAM_SET {
                                                                                param_set.reset(String::from("Sampler"),
                                                                                                String::from(""),
                                                                                                String::from(""),
                                                                                                String::from(""));
                                                                            } else {
                                                                                panic!("Can't get parameter set.");
                                                                            }
                                                                        },
                                                                        Rule::parameter => {
                                                                            for parameter_pair in sampler_pair.into_inner() {
                                                                                match parameter_pair.as_rule() {
                                                                                    Rule::bool_param => println!("TODO: Rule::bool_param"),
                                                                                    Rule::float_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if floats.len() == 1 {
                                                                                                param_set.add_float(string, floats[0]);
                                                                                            } else {
                                                                                                param_set.add_floats(string, floats);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::integer_param => {
                                                                                        let tuple: (String, Vec<i32>) = pbrt_integer_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let integers: Vec<i32> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if integers.len() == 1 {
                                                                                                param_set.add_int(string, integers[0]);
                                                                                            } else {
                                                                                                param_set.add_ints(string, integers);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::point_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if floats.len() == 3 {
                                                                                                param_set.add_point3f(string,
                                                                                                                      Point3f {
                                                                                                                          x: floats[0],
                                                                                                                          y: floats[1],
                                                                                                                          z: floats[2],
                                                                                                                      });
                                                                                            } else {
                                                                                                param_set.add_point3fs(string, floats);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::normal_param => println!("TODO: Rule::normal_param"),
                                                                                    Rule::rgb_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            param_set.add_rgb_spectrum(string,
                                                                                                                       Spectrum {
                                                                                                                           c: [floats[0], floats[1], floats[2]],
                                                                                                                       });
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::spectrum_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            param_set.add_rgb_spectrum(string,
                                                                                                                       Spectrum {
                                                                                                                           c: [floats[0], floats[1], floats[2]],
                                                                                                                       });
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::string_param => {
                                                                                        let tuple: (String, String) = pbrt_string_parameter(&mut parameter_pair.into_inner());
                                                                                        let string1: String = tuple.0;
                                                                                        let string2: String = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                                param_set.add_string(string1, string2);
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::texture_param => println!("TODO: Rule::texture_param"),
                                                                                    Rule::vector_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if floats.len() == 3 {
                                                                                                param_set.add_vector3f(string,
                                                                                                                       Vector3f {
                                                                                                                           x: floats[0],
                                                                                                                           y: floats[1],
                                                                                                                           z: floats[2],
                                                                                                                       });
                                                                                            } else {
                                                                                                param_set.add_point3fs(string, floats);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    _ => unreachable!()
                                                                                };
                                                                            }
                                                                        },
                                                                        _ => unreachable!()
                                                                    }
                                                                }
                                                                // we should have the sampler parameters by now
                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                    if let Some(ref mut ro) = RENDER_OPTIONS {
                                                                        println!("Sampler \"{}\" ", ro.sampler_name);
                                                                        ro.sampler_params.copy_from(param_set);
                                                                        print_params(&ro.sampler_params);
                                                                    } else {
                                                                        panic!("Can't get render options.");
                                                                    }
                                                                } else {
                                                                    panic!("Can't get parameter set.");
                                                                }
                                                            },
                                                            Rule::shape => println!("TODO: Rule::shape"),
                                                            Rule::texture => {
                                                                let mut counter: u8 = 0_u8;
                                                                let mut name: String = String::from("undefined");
                                                                let mut tex_type: String = String::from("undefined");
                                                                let mut tex_name;
                                                                for texture_pair in named_statement_pair.into_inner() {
                                                                    match texture_pair.as_rule() {
                                                                        Rule::string => {
                                                                            match counter {
                                                                                0 => {
                                                                                    // name
                                                                                    let mut string_pairs = texture_pair.into_inner();
                                                                                    let ident = string_pairs.next();
                                                                                    name = String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
                                                                                },
                                                                                1 => {
                                                                                    // tex_type
                                                                                    let mut string_pairs = texture_pair.into_inner();
                                                                                    let ident = string_pairs.next();
                                                                                    tex_type = String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
                                                                                },
                                                                                2 => {
                                                                                    // tex_name
                                                                                    let mut string_pairs = texture_pair.into_inner();
                                                                                    let ident = string_pairs.next();
                                                                                    tex_name = String::from_str(ident.unwrap().clone().into_span().as_str()).unwrap();
                                                                                    if let Some(ref mut param_set) = PARAM_SET {
                                                                                        param_set.reset(String::from("Texture"),
                                                                                                        String::from(name.clone()),
                                                                                                        String::from(tex_type.clone()),
                                                                                                        String::from(tex_name.clone()));
                                                                                    } else {
                                                                                        panic!("Can't get parameter set.");
                                                                                    }
                                                                                },
                                                                                _ => unreachable!()
                                                                            };
                                                                            counter += 1_u8;
                                                                        },
                                                                        Rule::parameter => {
                                                                            for parameter_pair in texture_pair.into_inner() {
                                                                                match parameter_pair.as_rule() {
                                                                                    Rule::bool_param => println!("TODO: Rule::bool_param"),
                                                                                    Rule::float_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if floats.len() == 1 {
                                                                                                param_set.add_float(string, floats[0]);
                                                                                            } else {
                                                                                                param_set.add_floats(string, floats);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::integer_param => {
                                                                                        let tuple: (String, Vec<i32>) = pbrt_integer_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let integers: Vec<i32> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if integers.len() == 1 {
                                                                                                param_set.add_int(string, integers[0]);
                                                                                            } else {
                                                                                                param_set.add_ints(string, integers);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::point_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if floats.len() == 3 {
                                                                                                param_set.add_point3f(string,
                                                                                                                      Point3f {
                                                                                                                          x: floats[0],
                                                                                                                          y: floats[1],
                                                                                                                          z: floats[2],
                                                                                                                      });
                                                                                            } else {
                                                                                                param_set.add_point3fs(string, floats);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::normal_param => println!("TODO: Rule::normal_param"),
                                                                                    Rule::rgb_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            param_set.add_rgb_spectrum(string,
                                                                                                                       Spectrum {
                                                                                                                           c: [floats[0], floats[1], floats[2]],
                                                                                                                       });
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::spectrum_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            param_set.add_rgb_spectrum(string,
                                                                                                                       Spectrum {
                                                                                                                           c: [floats[0], floats[1], floats[2]],
                                                                                                                       });
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::string_param => {
                                                                                        let tuple: (String, String) = pbrt_string_parameter(&mut parameter_pair.into_inner());
                                                                                        let string1: String = tuple.0;
                                                                                        let string2: String = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                                param_set.add_string(string1, string2);
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    Rule::texture_param => println!("TODO: Rule::texture_param"),
                                                                                    Rule::vector_param => {
                                                                                        let tuple: (String, Vec<Float>) = pbrt_float_parameter(&mut parameter_pair.into_inner());
                                                                                        let string: String = tuple.0;
                                                                                        let floats: Vec<Float> = tuple.1;
                                                                                        if let Some(ref mut param_set) = PARAM_SET {
                                                                                            if floats.len() == 3 {
                                                                                                param_set.add_vector3f(string,
                                                                                                                       Vector3f {
                                                                                                                           x: floats[0],
                                                                                                                           y: floats[1],
                                                                                                                           z: floats[2],
                                                                                                                       });
                                                                                            } else {
                                                                                                param_set.add_point3fs(string, floats);
                                                                                            }
                                                                                        } else {
                                                                                            panic!("Can't get parameter set.");
                                                                                        }
                                                                                    },
                                                                                    _ => unreachable!()
                                                                                };
                                                                            }
                                                                        },
                                                                        _ => unreachable!()
                                                                    }
                                                                }
                                                                // we should have the texture parameters by now
                                                                if let Some(ref mut param_set) = PARAM_SET {
                                                                    println!("Texture \"{}\" \"{}\" \"{}\" ",
                                                                             param_set.name,
                                                                             param_set.tex_type,
                                                                             param_set.tex_name);
                                                                    print_params(&param_set);
                                                                    make_texture(&param_set);
                                                                } else {
                                                                    panic!("Can't get parameter set.");
                                                                }
                                                            },
                                                            _ => unreachable!()
                                                        }
                                                    }
                                                },
                                                Rule::rotate => println!("TODO: Rule::rotate"),
                                                Rule::scale => println!("TODO: Rule::scale"),
                                                Rule::transform => println!("TODO: Rule::transform"),
                                                Rule::translate => {
                                                    let mut numbers: Vec<Float> = Vec::new();
                                                    for translate_pair in statement_pair.into_inner() {
                                                        let number: Float =
                                                            f32::from_str(translate_pair.clone().into_span().as_str()).unwrap();
                                                        numbers.push(number);
                                                    }
                                                    let x: Float = numbers[0];
                                                    let y: Float = numbers[1];
                                                    let z: Float = numbers[2];
                                                    let translate: Transform = Transform::translate(Vector3f { x: x, y: y, z: z, });
                                                    CUR_TRANSFORM.t[0] = CUR_TRANSFORM.t[0] * translate;
                                                    CUR_TRANSFORM.t[1] = CUR_TRANSFORM.t[1] * translate;
                                                },
                                                _ => unreachable!()
                                            };
                                        }
                                    },
                                    Rule::last_statement => {
                                        println!("WorldEnd");
                                        println!("TODO: Rule::last_statement");
                                        // pbrt_world_end();
                                    },
                                    _ => unreachable!()
                                }
                            }
                            // parser.main();
                            println!("done.");
                        }
                    }
                    None => panic!("No input file name."),
                }
                return;
            } else if matches.opt_present("v") {
                print_version(&program);
                return;
            } else {
                print_usage(&program, opts);
                return;
            }
        })
        .unwrap();
    let _res = child.join().unwrap();
}
