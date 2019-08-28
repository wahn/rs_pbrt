# rs_pbrt

[![Build Status](https://travis-ci.org/wahn/rs_pbrt.svg?branch=master)](https://travis-ci.org/wahn/rs_pbrt)
[![dependency status](https://deps.rs/repo/github/wahn/rs_pbrt/status.svg)](https://deps.rs/repo/github/wahn/rs_pbrt)
[![](https://img.shields.io/github/issues-raw/XAMPPRocky/tokei.svg)](https://github.com/wahn/rs_pbrt/issues)
[![](https://tokei.rs/b1/github/wahn/rs_pbrt?category=code)](https://github.com/Aaronepower/tokei)

**Rust** crate to implement at least parts of the PBRT book's C++ code:

http://www.pbrt.org

Current **Rust** documentation:

https://www.janwalter.org/doc/rust/pbrt/index.html

## Usage

```shell
> cargo run --release -- -h
Usage: target/release/rs_pbrt [options]

Options:
    -h, --help          print this help menu
    -i FILE             parse an input file
    -t, --nthreads NUM  use specified number of threads for rendering
    -v, --version       print version number
```

## Test Scenes

Some images of the test scenes are shown below, but you can find more
test scenes on [GitLab][test-scenes].

## Ganesha Statue

Very detailed scan of a small statue with over 4.3 million triangles,
illuminated by a few area light sources.

![Ganesha Statue](https://www.janwalter.org/assets/ganesha.png)

The scene can be found within the repository
(`assets/scenes/ganesha.tar.gz`).

## Subsurface Scattering (SSS)

![SSS Dragon](https://www.janwalter.org/assets/sss_dragon.png)

[sss_dragon_pbrt.tar.gz][sss_dragon_pbrt]

## Stochastic Progressive Photon Mapping (SPPM)

![SPPM Caustic
Glass](https://www.janwalter.org/assets/caustic_glass_pbrt_rust_sppm.png)

[caustic_glass.tar.gz][caustic_glass_pbrt]

## Ecosystem (Cover image for the first edition of the PBRT book)

![Ecosystem](https://www.janwalter.org/assets/ecosys.png)

[pbrt_ecosys.tar.gz][ecosys_pbrt]

## Hair

The [hair scattering][hair-scattering] model in action:

![Curly and straight hair rendered by Rust version of
PBRT](https://www.janwalter.org/assets/hair_rust_pbrt.png)

## Japanes Classroom by NovaZeeke

![Classroom room rendered by
rs_pbrt](https://www.janwalter.org/assets/classroom_pbrt_rust.png)

[classroom_pbrt.tar.gz][classroom_pbrt]

## The White Room by [Jay-Artist][jay-artist]

![The White Room rendered by
rs_pbrt](https://www.janwalter.org/assets/living-room-2_pbrt_rust_mlt.png)

[living-room-2_pbrt.tar.gz][living-room-2_pbrt]

## Country Kitchen by [Jay-Artist][jay-artist]

![Kitchen rendered by
rs_pbrt](https://www.janwalter.org/assets/kitchen_pbrt_rust.png)

[kitchen_pbrt.tar.gz][kitchen_pbrt]

## The Wooden Staircase by [Wig42][wig42]

![Staircase rendered by
rs_pbrt](https://www.janwalter.org/assets/staircase_pbrt_rust.png)

[staircase_pbrt.tar.gz][staircase_pbrt]

## Conference Room by Anat Grynberg and Greg Ward

![Conference room rendered by
rs_pbrt](https://www.janwalter.org/assets/conference_room_pbrt_rust_current.png)

[conference_room_pbrt.tar.gz][conference_room_pbrt]

## Theater by Charles Ehrlich and Greg Ward

![Theater rendered by
rs_pbrt](https://www.janwalter.org/assets/theater_pbrt_rust_corner.png)

![Theater rendered by
rs_pbrt](https://www.janwalter.org/assets/theater_pbrt_rust_stage.png)

[theater_pbrt.tar.gz][theater_pbrt]

For more info look at the [Wiki][wiki] page or the [release notes][release-notes].

Here you find another Rust implementation:

https://bitbucket.org/abusch/rustracer

[test-scenes]:          https://gitlab.com/jdb-walter/rs-pbrt-test-scenes/wikis/home
[wiki]:                 https://github.com/wahn/rs_pbrt/wiki
[release-notes]:        https://github.com/wahn/rs_pbrt/wiki/Release-Notes
[novazeeke]:            https://www.blendswap.com/user/NovaZeeke
[jay-artist]:           https://www.blendswap.com/user/Jay-Artist
[wig42]:                https://www.blendswap.com/user/Wig42
[classroom_pbrt]:       https://www.janwalter.org/Download/Scenes/PBRT/classroom_pbrt.tar.gz
[living-room-2_pbrt]:   https://www.janwalter.org/Download/Scenes/PBRT/living-room-2_pbrt.tar.gz
[kitchen_pbrt]:         https://www.janwalter.org/Download/Scenes/PBRT/kitchen_pbrt.tar.gz
[staircase_pbrt]:       https://www.janwalter.org/Download/Scenes/PBRT/staircase_pbrt.tar.gz
[conference_room_pbrt]: https://www.janwalter.org/Download/Scenes/conference_room_pbrt.tar.gz
[theater_pbrt]:         https://www.janwalter.org/Download/Scenes/theater_pbrt.tar.gz
[hair-scattering]:      http://www.pbrt.org/hair.pdf
[sss_dragon_pbrt]:      https://www.janwalter.org/Download/Scenes/sss_dragon_pbrt.tar.gz
[caustic_glass_pbrt]:   https://www.janwalter.org/Download/Scenes/caustic_glass.tar.gz
[ecosys_pbrt]:          https://www.janwalter.org/Download/Scenes/pbrt_ecosys.tar.gz
