use pest_derive::*;

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(Parser)]
#[grammar = "../examples/rs_pbrt.pest"]
struct PbrtParser;

// parser
use pest::Parser;

// command line options
use structopt::StructOpt;
// pbrt
use rs_pbrt::core::api::{
    pbrt_accelerator, pbrt_active_transform_all, pbrt_active_transform_end_time,
    pbrt_active_transform_start_time, pbrt_area_light_source, pbrt_attribute_begin,
    pbrt_attribute_end, pbrt_camera, pbrt_cleanup, pbrt_concat_transform, pbrt_coord_sys_transform,
    pbrt_film, pbrt_init, pbrt_integrator, pbrt_light_source, pbrt_look_at,
    pbrt_make_named_material, pbrt_make_named_medium, pbrt_material, pbrt_medium_interface,
    pbrt_named_material, pbrt_object_begin, pbrt_object_end, pbrt_object_instance,
    pbrt_pixel_filter, pbrt_reverse_orientation, pbrt_rotate, pbrt_sampler, pbrt_scale, pbrt_shape,
    pbrt_texture, pbrt_transform, pbrt_transform_begin, pbrt_transform_end, pbrt_translate,
    pbrt_world_begin,
};
use rs_pbrt::core::api::{ApiState, BsdfState};
use rs_pbrt::core::geometry::{Normal3f, Point2f, Point3f, Vector3f};
use rs_pbrt::core::paramset::ParamSet;
use rs_pbrt::core::pbrt::{Float, Spectrum};
use rs_pbrt::core::transform::Transform;
// std
use std::env;
use std::fs::File;
use std::io::BufReader;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::str::FromStr;

/// Parse a PBRT scene file (extension .pbrt) and render it.
#[derive(StructOpt)]
struct Cli {
    /// use specified number of threads for rendering
    #[structopt(short = "t", long = "nthreads", default_value = "0")]
    nthreads: u8,
    /// The path to the file to read
    #[structopt(parse(from_os_str))]
    path: std::path::PathBuf,
}

// Accelerator
// CoordinateSystem
// Identity
// TransformTimes

fn pbrt_bool_parameter(pairs: &mut pest::iterators::Pairs<Rule>) -> (String, bool) {
    // single string with or without brackets
    let ident = pairs.next();
    let string: String = String::from_str(ident.unwrap().clone().as_span().as_str()).unwrap();
    let option = pairs.next();
    let lbrack = option.clone().unwrap();
    let string2 = if lbrack.as_str() == "[" {
        // check for brackets
        let string = pairs.next();
        let pair = string.unwrap().clone();
        let ident = pair.into_inner().next();
        String::from_str(ident.unwrap().clone().as_span().as_str()).unwrap()
    } else {
        // no brackets
        let string = option.clone();
        let pair = string.unwrap().clone();
        let ident = pair.into_inner().next();
        String::from_str(ident.unwrap().clone().as_span().as_str()).unwrap()
    };
    // return boolean (instead of string)
    let b: bool;
    if string2 == "true" {
        b = true;
    } else if string2 == "false" {
        b = false
    } else {
        println!(
            "WARNING: parameter {:?} not well defined, defaulting to false",
            string
        );
        b = false
    }
    (string, b)
}

fn pbrt_float_parameter(pairs: &mut pest::iterators::Pairs<Rule>) -> (String, Vec<Float>) {
    let mut floats: Vec<Float> = Vec::new();
    // single float or several floats using brackets
    let ident = pairs.next();
    let string: String = String::from_str(ident.unwrap().clone().as_span().as_str()).unwrap();
    let option = pairs.next();
    let lbrack = option.clone().unwrap();
    if lbrack.as_str() == "[" {
        // check for brackets
        let mut number = pairs.next();
        while number.is_some() {
            let pair = number.unwrap().clone();
            if pair.as_str() == "]" {
                // closing bracket found
                break;
            } else {
                let float: Float = f32::from_str(pair.as_span().as_str()).unwrap();
                floats.push(float);
            }
            number = pairs.next();
        }
    } else {
        // no brackets
        let mut number = option.clone();
        while number.is_some() {
            let pair = number.unwrap().clone();
            let float: Float = f32::from_str(pair.as_span().as_str()).unwrap();
            floats.push(float);
            number = pairs.next();
        }
    }
    (string, floats)
}

fn pbrt_integer_parameter(pairs: &mut pest::iterators::Pairs<Rule>) -> (String, Vec<i32>) {
    let mut integers: Vec<i32> = Vec::new();
    // single integer or several integers using brackets
    let ident = pairs.next();
    let string: String = String::from_str(ident.unwrap().clone().as_span().as_str()).unwrap();
    let option = pairs.next();
    let lbrack = option.clone().unwrap();
    if lbrack.as_str() == "[" {
        // check for brackets
        let mut number = pairs.next();
        while number.is_some() {
            let pair = number.unwrap().clone();
            if pair.as_str() == "]" {
                // closing bracket found
                break;
            } else {
                let integer: i32 = i32::from_str(pair.as_span().as_str()).unwrap();
                integers.push(integer);
            }
            number = pairs.next();
        }
    } else {
        // no brackets
        let mut number = option.clone();
        while number.is_some() {
            let pair = number.unwrap().clone();
            let integer: i32 = i32::from_str(pair.as_span().as_str()).unwrap();
            integers.push(integer);
            number = pairs.next();
        }
    }
    (string, integers)
}

fn pbrt_string_parameter(pairs: &mut pest::iterators::Pairs<Rule>) -> (String, String) {
    // single string with or without brackets
    let ident = pairs.next();
    let string1: String = String::from_str(ident.unwrap().clone().as_span().as_str()).unwrap();
    let option = pairs.next();
    let lbrack = option.clone().unwrap();
    let string2 = if lbrack.as_str() == "[" {
        // check for brackets
        let string = pairs.next();
        let pair = string.unwrap().clone();
        let ident = pair.into_inner().next();
        String::from_str(ident.unwrap().clone().as_span().as_str()).unwrap()
    } else {
        // no brackets
        let string = option.clone();
        let pair = string.unwrap().clone();
        let ident = pair.into_inner().next();
        String::from_str(ident.unwrap().clone().as_span().as_str()).unwrap()
    };
    (string1, string2)
}

fn pbrt_texture_parameter(pairs: &mut pest::iterators::Pairs<Rule>) -> (String, String) {
    // single string with or without brackets
    let ident = pairs.next();
    let string1: String = String::from_str(ident.unwrap().clone().as_span().as_str()).unwrap();
    let option = pairs.next();
    let lbrack = option.clone().unwrap();
    let string2 = if lbrack.as_str() == "[" {
        // check for brackets
        let string = pairs.next();
        let pair = string.unwrap().clone();
        let ident = pair.into_inner().next();
        String::from_str(ident.unwrap().clone().as_span().as_str()).unwrap()
    } else {
        // no brackets
        let string = option.clone();
        let pair = string.unwrap().clone();
        let ident = pair.into_inner().next();
        String::from_str(ident.unwrap().clone().as_span().as_str()).unwrap()
    };
    (string1, string2)
}

fn extract_params(key_word: String, pairs: pest::iterators::Pair<Rule>) -> ParamSet {
    let mut params: ParamSet = ParamSet::default();
    params.key_word = key_word;
    let mut counter: u8 = 0_u8;
    for pair in pairs.into_inner() {
        // let span = pair.clone().as_span();
        // println!("Rule:    {:?}", pair.as_rule());
        // println!("Span:    {:?}", span);
        // println!("Text:    {}", span.as_str());
        match pair.as_rule() {
            Rule::identifier => {
                // ignore (was added above)
            }
            Rule::empty_string => {}
            Rule::string => {
                match counter {
                    0 => {
                        // name
                        let mut string_pairs = pair.into_inner();
                        let ident = string_pairs.next();
                        params.name =
                            String::from_str(ident.unwrap().clone().as_span().as_str()).unwrap();
                    }
                    1 => {
                        // tex_type
                        let mut string_pairs = pair.into_inner();
                        let ident = string_pairs.next();
                        params.tex_type =
                            String::from_str(ident.unwrap().clone().as_span().as_str()).unwrap();
                    }
                    2 => {
                        // tex_name
                        let mut string_pairs = pair.into_inner();
                        let ident = string_pairs.next();
                        params.tex_name =
                            String::from_str(ident.unwrap().clone().as_span().as_str()).unwrap();
                    }
                    _ => unreachable!(),
                };
                counter += 1_u8;
            }
            Rule::type_name => {
                // name
                let mut string_pairs = pair.into_inner();
                let ident = string_pairs.next();
                params.name = String::from_str(ident.unwrap().clone().as_span().as_str()).unwrap();
            }
            Rule::file_name => {
                // name
                let mut string_pairs = pair.into_inner();
                let ident = string_pairs.next();
                params.name = String::from_str(ident.unwrap().clone().as_span().as_str()).unwrap();
            }
            Rule::parameter => {
                for parameter_pair in pair.into_inner() {
                    // println!("DEBUG: {:?}", parameter_pair.as_rule());
                    match parameter_pair.as_rule() {
                        Rule::bool_param => {
                            let tuple: (String, bool) =
                                pbrt_bool_parameter(&mut parameter_pair.into_inner());
                            let string: String = tuple.0;
                            let b: bool = tuple.1;
                            params.add_bool(string, b);
                        }
                        Rule::blackbody_param => {
                            let tuple: (String, Vec<Float>) =
                                pbrt_float_parameter(&mut parameter_pair.into_inner());
                            let string: String = tuple.0;
                            let floats: Vec<Float> = tuple.1;
                            params.add_blackbody_spectrum(string, floats);
                        }
                        Rule::float_param => {
                            let tuple: (String, Vec<Float>) =
                                pbrt_float_parameter(&mut parameter_pair.into_inner());
                            let string: String = tuple.0;
                            let floats: Vec<Float> = tuple.1;
                            if floats.len() == 1 {
                                params.add_float(string, floats[0]);
                            } else {
                                params.add_floats(string, floats);
                            }
                        }
                        Rule::integer_param => {
                            let tuple: (String, Vec<i32>) =
                                pbrt_integer_parameter(&mut parameter_pair.into_inner());
                            let string: String = tuple.0;
                            let integers: Vec<i32> = tuple.1;
                            if integers.len() == 1 {
                                params.add_int(string, integers[0]);
                            } else {
                                params.add_ints(string, integers);
                            }
                        }
                        Rule::point_param => {
                            let tuple: (String, Vec<Float>) =
                                pbrt_float_parameter(&mut parameter_pair.into_inner());
                            let string: String = tuple.0;
                            let floats: Vec<Float> = tuple.1;
                            if floats.len() == 3 {
                                params.add_point3f(
                                    string,
                                    Point3f {
                                        x: floats[0],
                                        y: floats[1],
                                        z: floats[2],
                                    },
                                );
                            } else {
                                params.add_point3fs(string, floats);
                            }
                        }
                        Rule::point2_param => {
                            let tuple: (String, Vec<Float>) =
                                pbrt_float_parameter(&mut parameter_pair.into_inner());
                            let string: String = tuple.0;
                            let floats: Vec<Float> = tuple.1;
                            if floats.len() == 2 {
                                params.add_point2f(
                                    string,
                                    Point2f {
                                        x: floats[0],
                                        y: floats[1],
                                    },
                                );
                            } else {
                                params.add_point2fs(string, floats);
                            }
                        }
                        Rule::normal_param => {
                            let tuple: (String, Vec<Float>) =
                                pbrt_float_parameter(&mut parameter_pair.into_inner());
                            let string: String = tuple.0;
                            let floats: Vec<Float> = tuple.1;
                            if floats.len() == 3 {
                                params.add_normal3f(
                                    string,
                                    Normal3f {
                                        x: floats[0],
                                        y: floats[1],
                                        z: floats[2],
                                    },
                                );
                            } else {
                                params.add_normal3fs(string, floats);
                            }
                        }
                        Rule::rgb_param => {
                            let tuple: (String, Vec<Float>) =
                                pbrt_float_parameter(&mut parameter_pair.into_inner());
                            let string: String = tuple.0;
                            let floats: Vec<Float> = tuple.1;
                            params.add_rgb_spectrum(
                                string,
                                Spectrum {
                                    c: [floats[0], floats[1], floats[2]],
                                },
                            );
                        }
                        Rule::spectrum_param => {
                            // TODO: "spectrum Kd" [ 300 .3  400 .6   410 .65  415 .8  500 .2  600 .1 ]
                            // let tuple: (String, Vec<Float>) =
                            //     pbrt_float_parameter(&mut parameter_pair.into_inner());
                            // let string: String = tuple.0;
                            // let floats: Vec<Float> = tuple.1;
                            // params.add_rgb_spectrum(
                            //     string,
                            //     Spectrum {
                            //         c: [floats[0], floats[1], floats[2]],
                            //     },
                            // );
                            // or
                            // "spectrum Kd" "filename"
                            let tuple: (String, String) =
                                pbrt_string_parameter(&mut parameter_pair.into_inner());
                            let string1: String = tuple.0;
                            let string2: String = tuple.1;
                            let mut strings: Vec<String> = Vec::with_capacity(1_usize);
                            strings.push(string2);
                            params.add_sampled_spectrum_files(string1, strings);
                        }
                        Rule::string_param => {
                            let tuple: (String, String) =
                                pbrt_string_parameter(&mut parameter_pair.into_inner());
                            let string1: String = tuple.0;
                            let string2: String = tuple.1;
                            params.add_string(string1, string2);
                        }
                        Rule::texture_param => {
                            let tuple: (String, String) =
                                pbrt_texture_parameter(&mut parameter_pair.into_inner());
                            let string1: String = tuple.0;
                            let string2: String = tuple.1;
                            params.add_texture(string1, string2);
                        }
                        Rule::vector_param => {
                            let tuple: (String, Vec<Float>) =
                                pbrt_float_parameter(&mut parameter_pair.into_inner());
                            let string: String = tuple.0;
                            let floats: Vec<Float> = tuple.1;
                            if floats.len() == 3 {
                                params.add_vector3f(
                                    string,
                                    Vector3f {
                                        x: floats[0],
                                        y: floats[1],
                                        z: floats[2],
                                    },
                                );
                            } else {
                                params.add_vector3fs(string, floats);
                            }
                        }
                        // TODO: more rules
                        _ => println!("TODO: {:?}", parameter_pair.as_rule()),
                    }
                }
            }
            _ => println!("TODO: {:?}", pair.as_rule()),
        }
    }
    params
}

fn parse_line(
    api_state: &mut ApiState,
    bsdf_state: &mut BsdfState,
    identifier: &str,
    str_buf: String,
) {
    if str_buf == "" {
        // no additional arguments
        match identifier {
            "AttributeBegin" => {
                // AttributeBegin
                // println!("{} {}", identifier, str_buf);
                pbrt_attribute_begin(api_state);
            }
            "AttributeEnd" => {
                // AttributeEnd
                // println!("{} {}", identifier, str_buf);
                pbrt_attribute_end(api_state);
            }
            "ObjectEnd" => {
                // ObjectEnd
                // println!("{} {}", identifier, str_buf);
                pbrt_object_end(api_state);
            }
            "ReverseOrientation" => {
                // ReverseOrientation
                // println!("{} {}", identifier, str_buf);
                pbrt_reverse_orientation(api_state);
            }
            "TransformBegin" => {
                // TransformBegin
                pbrt_transform_begin(api_state);
            }
            "TransformEnd" => {
                // TransformEnd
                pbrt_transform_end(api_state);
            }
            "WorldBegin" => {
                // WorldBegin
                // println!("{} {}", identifier, str_buf);
                pbrt_world_begin(api_state);
            }
            "WorldEnd" => {
                // WorldEnd
                // println!("{} {}", identifier, str_buf);
                pbrt_cleanup(api_state);
            }
            _ => println!("{} {:?}", identifier, str_buf),
        }
    } else {
        let statement = String::from(identifier) + " " + &str_buf;
        // println!("DEBUG: {:?}", &statement);
        let pairs = PbrtParser::parse(Rule::name_and_or_params, &statement)
            .expect("unsuccessful parse")
            .next()
            .unwrap();
        for inner_pair in pairs.into_inner() {
            // println!("DEBUG: {:?}", inner_pair.as_rule());
            match inner_pair.as_rule() {
                Rule::type_params => {
                    // identifier "type" parameter-list
                    let for_printing = inner_pair.as_str();
                    // println!("DEBUG: {}", for_printing);
                    let params = extract_params(String::from(identifier), inner_pair);
                    match identifier {
                        "Accelerator" => {
                            // Accelerator
                            pbrt_accelerator(api_state, params);
                        }
                        "AreaLightSource" => {
                            // AreaLightSource
                            pbrt_area_light_source(api_state, params);
                        }
                        "Camera" => {
                            // Camera
                            pbrt_camera(api_state, params);
                        }
                        "CoordSysTransform" => {
                            // CoordSysTransform
                            pbrt_coord_sys_transform(api_state, params);
                        }
                        "Film" => {
                            // Film
                            pbrt_film(api_state, params);
                        }
                        "Include" => {
                            // Include
                            let mut include_file: String = params.name.clone();
                            if let Some(ref search_directory) = api_state.search_directory {
                                let mut path_buf: PathBuf = PathBuf::from("/");
                                path_buf.push(search_directory.as_ref());
                                path_buf.push(params.name);
                                include_file = String::from(path_buf.to_str().unwrap());
                                // println!("DEBUG: {:?}", include_file);
                            }
                            let todo: Vec<&str> = for_printing.splitn(3, '"').collect();
                            println!("Include {:?}", include_file);
                            parse_file(include_file, api_state, bsdf_state, todo[2]);
                        }
                        "Integrator" => {
                            // Integrator
                            pbrt_integrator(api_state, params);
                        }
                        "LightSource" => {
                            // LightSource
                            pbrt_light_source(api_state, params);
                        }
                        "MakeNamedMaterial" => {
                            // MakeNamedMaterial
                            pbrt_make_named_material(api_state, bsdf_state, params);
                        }
                        "MakeNamedMedium" => {
                            // MakeNamedMedium
                            pbrt_make_named_medium(api_state, params);
                        }
                        "Material" => {
                            // Material
                            pbrt_material(api_state, params);
                        }
                        "NamedMaterial" => {
                            // NamedMaterial
                            pbrt_named_material(api_state, params);
                        }
                        "ObjectBegin" => {
                            // ObjectBegin
                            pbrt_object_begin(api_state, params);
                        }
                        "ObjectInstance" => {
                            // ObjectInstance
                            pbrt_object_instance(api_state, params);
                        }
                        "PixelFilter" => {
                            // PixelFilter
                            pbrt_pixel_filter(api_state, params);
                        }
                        "Sampler" => {
                            // Sampler
                            pbrt_sampler(api_state, params);
                        }
                        "Shape" => {
                            // Shape
                            pbrt_shape(api_state, bsdf_state, params);
                        }
                        "Texture" => {
                            // Texture
                            pbrt_texture(api_state, params);
                        }
                        _ => println!("> {}", for_printing),
                    }
                }
                Rule::active_transform => {
                    // ActiveTransform
                    for rule_pair in inner_pair.into_inner() {
                        match rule_pair.as_rule() {
                            Rule::all => {
                                pbrt_active_transform_all(api_state);
                            }
                            Rule::start_time => {
                                pbrt_active_transform_start_time(api_state);
                            }
                            Rule::end_time => {
                                pbrt_active_transform_end_time(api_state);
                            }
                            _ => unreachable!(),
                        }
                    }
                }
                Rule::concat_transform => {
                    // ConcatTransform m00 .. m33
                    let mut m: Vec<Float> = Vec::new();
                    for rule_pair in inner_pair.into_inner() {
                        // ignore brackets
                        let not_opening: bool = rule_pair.as_str() != "[";
                        let not_closing: bool = rule_pair.as_str() != "]";
                        if not_opening && not_closing {
                            let number: Float =
                                f32::from_str(rule_pair.clone().as_span().as_str()).unwrap();
                            m.push(number);
                        }
                    }
                    let m00: Float = m[0];
                    let m01: Float = m[1];
                    let m02: Float = m[2];
                    let m03: Float = m[3];
                    let m10: Float = m[4];
                    let m11: Float = m[5];
                    let m12: Float = m[6];
                    let m13: Float = m[7];
                    let m20: Float = m[8];
                    let m21: Float = m[9];
                    let m22: Float = m[10];
                    let m23: Float = m[11];
                    let m30: Float = m[12];
                    let m31: Float = m[13];
                    let m32: Float = m[14];
                    let m33: Float = m[15];
                    let tr: Transform = Transform::new(
                        m00, m10, m20, m30, m01, m11, m21, m31, m02, m12, m22, m32, m03, m13, m23,
                        m33,
                    );
                    pbrt_concat_transform(api_state, &tr);
                }
                Rule::look_at => {
                    // LookAt eye_x eye_y eye_z look_x look_y look_z up_x up_y up_z
                    let mut v: Vec<Float> = Vec::new();
                    for rule_pair in inner_pair.into_inner() {
                        let number: Float =
                            f32::from_str(rule_pair.clone().as_span().as_str()).unwrap();
                        v.push(number);
                    }
                    // println!(
                    //     "LookAt {} {} {} {} {} {} {} {} {}",
                    //     v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8],
                    // );
                    pbrt_look_at(
                        api_state, v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8],
                    );
                }
                Rule::medium_interface => {
                    // MediumInterface
                    let mut strings: Vec<String> = Vec::new();
                    for rule_pair in inner_pair.into_inner() {
                        match rule_pair.as_rule() {
                            Rule::empty_string => {
                                strings.push(String::from(""));
                            }
                            Rule::string => {
                                let ident = rule_pair.into_inner().next();
                                let string: String =
                                    String::from_str(ident.unwrap().clone().as_span().as_str())
                                        .unwrap();
                                strings.push(string);
                            }
                            _ => unreachable!(),
                        }
                    }
                    assert!(
                        strings.len() == 2_usize,
                        "ERROR: expected two strings, found {:?}",
                        strings.len()
                    );
                    pbrt_medium_interface(api_state, &strings[0], &strings[1]);
                }
                Rule::rotate => {
                    // Rotate angle x y z
                    let mut v: Vec<Float> = Vec::new();
                    for rule_pair in inner_pair.into_inner() {
                        let number: Float =
                            f32::from_str(rule_pair.clone().as_span().as_str()).unwrap();
                        v.push(number);
                    }
                    // println!("Rotate {} {} {} {}", v[0], v[1], v[2], v[3]);
                    pbrt_rotate(api_state, v[0], v[1], v[2], v[3]);
                }
                Rule::scale => {
                    // Scale x y z
                    let mut v: Vec<Float> = Vec::new();
                    for rule_pair in inner_pair.into_inner() {
                        let number: Float =
                            f32::from_str(rule_pair.clone().as_span().as_str()).unwrap();
                        v.push(number);
                    }
                    // println!("Scale {} {} {}", v[0], v[1], v[2]);
                    pbrt_scale(api_state, v[0], v[1], v[2]);
                }
                Rule::transform => {
                    // Transform m00 .. m33
                    let mut m: Vec<Float> = Vec::new();
                    for rule_pair in inner_pair.into_inner() {
                        // ignore brackets
                        let not_opening: bool = rule_pair.as_str() != "[";
                        let not_closing: bool = rule_pair.as_str() != "]";
                        if not_opening && not_closing {
                            let number: Float =
                                f32::from_str(rule_pair.clone().as_span().as_str()).unwrap();
                            m.push(number);
                        }
                    }
                    let m00: Float = m[0];
                    let m01: Float = m[1];
                    let m02: Float = m[2];
                    let m03: Float = m[3];
                    let m10: Float = m[4];
                    let m11: Float = m[5];
                    let m12: Float = m[6];
                    let m13: Float = m[7];
                    let m20: Float = m[8];
                    let m21: Float = m[9];
                    let m22: Float = m[10];
                    let m23: Float = m[11];
                    let m30: Float = m[12];
                    let m31: Float = m[13];
                    let m32: Float = m[14];
                    let m33: Float = m[15];
                    let tr: Transform = Transform::new(
                        m00, m10, m20, m30, m01, m11, m21, m31, m02, m12, m22, m32, m03, m13, m23,
                        m33,
                    );
                    pbrt_transform(api_state, &tr);
                }
                Rule::translate => {
                    // Translate x y z
                    let mut v: Vec<Float> = Vec::new();
                    for rule_pair in inner_pair.into_inner() {
                        let number: Float =
                            f32::from_str(rule_pair.clone().as_span().as_str()).unwrap();
                        v.push(number);
                    }
                    // println!("Translate {} {} {}", v[0], v[1], v[2]);
                    pbrt_translate(api_state, v[0], v[1], v[2]);
                }
                Rule::remaining_line => {
                    // predetermined number of arguments of predetermined type
                    println!("< {}", inner_pair.as_str());
                }
                // _ => unreachable!(),
                _ => println!("TODO: {:?}", inner_pair.as_rule()),
            }
        }
    }
}

fn parse_file(
    filename: String,
    api_state: &mut ApiState,
    bsdf_state: &mut BsdfState,
    append: &str,
) {
    // println!("FILE = {}", x);
    let f = File::open(filename.clone()).unwrap();
    let ip: &Path = Path::new(filename.as_str());
    if ip.is_relative() {
        let cp: PathBuf = env::current_dir().unwrap();
        let pb: PathBuf = cp.join(ip);
        let search_directory: &Path = pb.as_path().parent().unwrap();
        // println!("search_directory is {}", search_directory.display());
        api_state.search_directory = Some(Box::new(PathBuf::from(search_directory)));
    }
    let mut reader = BufReader::new(f);
    let mut str_buf: String = String::default();
    let _num_bytes = reader.read_to_string(&mut str_buf);
    // if num_bytes.is_ok() {
    //     let n_bytes = num_bytes.unwrap();
    //     println!("{} bytes read", n_bytes);
    // }
    if append != "" {
        str_buf += append;
        str_buf += "\n";
    }
    let pairs = PbrtParser::parse(Rule::pbrt, &str_buf)
        .expect("unsuccessful parse")
        .next()
        .unwrap();
    let mut identifier: &str = "";
    // let mut comment_count: u64 = 0;
    // let mut empty_count: u64 = 0;
    // let mut todo_count: u64 = 0;
    let mut parse_again: String = String::default();
    // first parse file line by line
    for inner_pair in pairs.into_inner() {
        match inner_pair.as_rule() {
            // comment lines (starting with '#')
            Rule::comment_line => {
                // comment_count += 1;
            }
            Rule::statement_line => {
                for statement_pair in inner_pair.into_inner() {
                    match statement_pair.as_rule() {
                        Rule::identifier => {
                            if identifier != "" {
                                parse_line(api_state, bsdf_state, identifier, parse_again.clone());
                            }
                            identifier = statement_pair.as_str();
                            parse_again = String::default();
                        }
                        Rule::remaining_line => {
                            if parse_again != "" {
                                parse_again = parse_again + " " + statement_pair.as_str();
                            } else {
                                parse_again += statement_pair.as_str();
                            }
                        }
                        Rule::trailing_comment => {
                            // ignore (only if there are no '"' chars)
                            if statement_pair.as_str().contains('\"') {
                                if parse_again != "" {
                                    parse_again = parse_again + " " + statement_pair.as_str();
                                } else {
                                    parse_again += statement_pair.as_str();
                                }
                            }
                        }
                        _ => println!("TODO: {:?}", statement_pair.as_rule()),
                    }
                }
            }
            Rule::empty_line => {
                // empty_count += 1;
            }
            Rule::todo_line => {
                // todo_count += 1;
                for params_pair in inner_pair.into_inner() {
                    match params_pair.as_rule() {
                        Rule::remaining_params => {
                            if parse_again != "" {
                                parse_again = parse_again + " " + params_pair.as_str();
                            } else {
                                parse_again += params_pair.as_str();
                            }
                        }
                        Rule::trailing_comment => {
                            // ignore
                        }
                        _ => println!("TODO: {:?}", params_pair.as_rule()),
                    }
                }
            }
            Rule::EOI => parse_line(api_state, bsdf_state, identifier, parse_again.clone()),
            _ => unreachable!(),
        }
    }
    // println!("Number of comment line(s):   {}", comment_count);
    // println!("Number of parameter line(s): {}", todo_count);
    // println!("Number of empty line(s):     {}", empty_count);
}

fn main() {
    // handle command line options
    let args = Cli::from_args();
    let number_of_threads: u8 = args.nthreads;
    let num_cores = num_cpus::get();
    let git_describe = option_env!("GIT_DESCRIBE").unwrap_or("unknown");
    println!(
        "pbrt version {} ({}) [Detected {} cores]",
        VERSION, git_describe, num_cores
    );
    println!("Copyright (c) 2016-2021 Jan Douglas Bert Walter.");
    println!("Rust code based on C++ code by Matt Pharr, Greg Humphreys, and Wenzel Jakob.");
    let (mut api_state, mut bsdf_state) = pbrt_init(number_of_threads);
    parse_file(
        args.path.into_os_string().into_string().unwrap(),
        &mut api_state,
        &mut bsdf_state,
        "",
    );
}
