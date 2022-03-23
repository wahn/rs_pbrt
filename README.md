# rs_pbrt

[![crates.io](https://img.shields.io/crates/v/rs_pbrt.svg)](https://crates.io/crates/rs_pbrt)
[![Documentation Status](https://readthedocs.org/projects/rs-pbrt/badge/?version=latest)](https://rs-pbrt.readthedocs.io/en/latest/?badge=latest)
[![dependency status](https://deps.rs/repo/github/wahn/rs_pbrt/status.svg)](https://deps.rs/repo/github/wahn/rs_pbrt)
[![builds.sr.ht status](https://builds.sr.ht/~wahn/rs-pbrt.svg)](https://builds.sr.ht/~wahn/rs-pbrt?)
<!-- [![](https://tokei.rs/b1/github/wahn/rs_pbrt?category=code)](https://github.com/wahn/rs_pbrt) -->
<!-- [![](https://img.shields.io/github/release-date/wahn/rs_pbrt.svg)](https://github.com/wahn/rs_pbrt/releases) -->
<!-- [![](https://img.shields.io/github/issues-raw/wahn/rs_pbrt.svg)](https://github.com/wahn/rs_pbrt/issues) -->
<!-- [![Build Status](https://travis-ci.com/wahn/rs_pbrt.svg?branch=master)](https://travis-ci.com/wahn/rs_pbrt) -->

You can find more information about `rs_pbrt` at https://www.rs-pbrt.org/about ...

**Rust** crate to implement a counterpart to the PBRT book's (3rd edition) C++ code:

http://www.pbrt.org

Current **Rust** (development) documentation:

https://www.janwalter.org/doc/rust/rs_pbrt/index.html
or
https://www.rs-pbrt.org/doc/crates/rs_pbrt/index.html

## Usage

```shell
> cargo build --release --no-default-features
> ./target/release/rs_pbrt --help
rs_pbrt 0.9.6
Parse a PBRT scene file (extension .pbrt) and render it

USAGE:
    rs_pbrt [OPTIONS] <path>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
        --cropx0 <cropx0>            Specify an image crop window <x0 x1 y0 y1> [default: 0.0]
        --cropx1 <cropx1>            Specify an image crop window <x0 x1 y0 y1> [default: 1.0]
        --cropy0 <cropy0>            Specify an image crop window <x0 x1 y0 y1> [default: 0.0]
        --cropy1 <cropy1>            Specify an image crop window <x0 x1 y0 y1> [default: 1.0]
    -i, --integrator <integrator>    ao, directlighting, whitted, path, bdpt, mlt, sppm, volpath
    -t, --nthreads <nthreads>        use specified number of threads for rendering [default: 0]
    -s, --samples <samples>          pixel samples [default: 0]

ARGS:
    <path>    The path to the file to read
```

## Test Scenes

Some images of the test scenes are shown below, but you can find more
test scenes on [GitLab][test-scenes].

## Ganesha Statue

Very detailed scan of a small statue with over 4.3 million triangles,
illuminated by a few area light sources.

![Ganesha Statue](https://www.janwalter.org/assets/ganesha.png)

The scene can be found within the [GitLab repository][gitlab_pbrt]
(`pbrt_ganesha.tar.gz`).

## Subsurface Scattering (SSS)

![SSS Dragon](https://www.janwalter.org/assets/sss_dragon.png)

The scene can be found within the [GitLab repository][gitlab_pbrt]
(`sss_dragon.tar.gz`).

## Stochastic Progressive Photon Mapping (SPPM)

![SPPM Caustic
Glass](https://www.janwalter.org/assets/caustic_glass_pbrt_rust_sppm.png)

The scene can be found within the [GitLab repository][gitlab_pbrt]
(`caustic_glass.tar.gz`).

## Ecosystem (Cover image for the first edition of the PBRT book)

![Ecosystem](https://www.janwalter.org/assets/ecosys.png)

The scene can be found within the [GitLab repository][gitlab_pbrt]
(`pbrt_ecosys.tar.gz`).

## Landscape (Cover image for the third edition of the PBRT book)

![Landscape](https://www.janwalter.org/assets/landscape_rust_pbrt_view_0.png)

## Hair

The [hair scattering][hair-scattering] model in action:

![Curly and straight hair rendered by Rust version of
PBRT](https://www.janwalter.org/assets/hair_rust_pbrt.png)

The scene can be found within the [GitLab repository][gitlab_pbrt]
(`pbrt_hair.tar.gz`).

## Japanes Classroom by NovaZeeke

![Classroom room rendered by
rs_pbrt](https://www.janwalter.org/assets/classroom_pbrt_rust.png)

The scene can be found within the [GitLab repository][gitlab_pbrt]
(`bitterli_classroom.tar.gz`).

## The White Room by [Jay-Artist][jay-artist]

![The White Room rendered by
rs_pbrt](https://www.janwalter.org/assets/living-room-2_pbrt_rust_mlt.png)

The scene can be found within the [GitLab repository][gitlab_pbrt]
(`bitterli_living-room-2.tar.gz`).

## Country Kitchen by [Jay-Artist][jay-artist]

![Kitchen rendered by
rs_pbrt](https://www.janwalter.org/assets/kitchen_pbrt_rust.png)

The scene can be found within the [GitLab repository][gitlab_pbrt]
(`bitterli_kitchen.tar.gz`).

## The Wooden Staircase by [Wig42][wig42]

![Staircase rendered by
rs_pbrt](https://www.janwalter.org/assets/staircase_pbrt_rust.png)

The scene can be found within the [GitLab repository][gitlab_pbrt]
(`bitterli_staircase.tar.gz`).

## Conference Room by Anat Grynberg and Greg Ward

![Conference room rendered by
rs_pbrt](https://www.janwalter.org/assets/conference_room_pbrt_rust_current.png)

The scene can be found within the [GitLab repository][gitlab_pbrt]
(`radiance_conference_room` folder).

## Theater by Charles Ehrlich and Greg Ward

![Theater rendered by
rs_pbrt](https://www.janwalter.org/assets/theater_pbrt_rust_corner.png)

![Theater rendered by
rs_pbrt](https://www.janwalter.org/assets/theater_pbrt_rust_stage.png)

The scene can be found within the [GitLab repository][gitlab_pbrt]
(`radiance_theater` folder).

For more information check out the `rs-pbrt` [web page][rs-pbrt].

<!-- For more info look at the [Wiki][wiki] page or the [release notes][release-notes]. -->

<!-- Here you find another Rust implementation: -->

<!-- https://bitbucket.org/abusch/rustracer -->

## License

Licensed under either of

* Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE) or
  http://www.apache.org/licenses/LICENSE-2.0)
* MIT license ([LICENSE-MIT](LICENSE-MIT) or
  http://opensource.org/licenses/MIT)

at your option.

### Contribution

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in the work by you, as defined in the
Apache-2.0 license, shall be dual licensed as above, without any
additional terms or conditions.

[test-scenes]:          https://gitlab.com/jdb-walter/rs-pbrt-test-scenes/wikis/home
[gitlab_pbrt]:          https://gitlab.com/jdb-walter/rs-pbrt-test-scenes/-/tree/master/pbrt
[jay-artist]:           https://www.blendswap.com/user/Jay-Artist
[wig42]:                https://www.blendswap.com/user/Wig42
[hair-scattering]:      http://www.pbrt.org/hair.pdf
[rs-pbrt]:              https://www.rs-pbrt.org
