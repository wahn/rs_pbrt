# rs_pbrt

[![Build Status](https://travis-ci.org/wahn/rs_pbrt.svg?branch=master)](https://travis-ci.org/wahn/rs_pbrt)

**Rust** crate to implement at least parts of the PBRT book's C++ code:

http://www.pbrt.org

Current **Rust** documentation:

https://www.janwalter.org/doc/rust/pbrt/index.html

Scene with a **glass** material on the first sphere and a **mirror**
material on the second sphere. The ground triangles use a procedural
**checker** texture on a **matte** material. Rendered via the **Rust**
version of **PBRT**:

![Rendered via Rust version of PBRT](https://www.janwalter.org/assets/spheres-differentials-texfilt_v0_1_5.png)

Same scene with a non-procedural texture:

![Rendered via Rust version of PBRT](https://www.janwalter.org/assets/spheres-differentials-texfilt_v0_1_6.png)

The scene above can be rendered with one of the compiled example programs:

```shell
> ./target/release/examples/pbrt_spheres_differentials_texfilt --help
Usage: ./target/release/examples/pbrt_spheres_differentials_texfilt [options]

Options:
    -h, --help          print this help menu
    -c, --checker       use procedural texture
    -i, --image         use image texture
    -n, --none          use no texture
    -m, --matte         use only matte materials
    -v, --version       print version number
```

There are also a couple of scenes which can be rendered from the
command line by parsing a scene description file:

```shell
> ./target/release/examples/pest_test -i assets/scenes/teapot-area-light.pbrt
```

![Teapot scene rendered via Rust version of PBRT](https://www.janwalter.org/assets/teapot-area-light_v0.1.12.png)

```
> ./target/release/examples/pest_test -i assets/scenes/cornell_box.pbrt
```

With direct lighting:

```
Integrator "directlighting" "integer maxdepth" [10]
```

![Cornell Box scene rendered via Rust version of PBRT](https://www.janwalter.org/assets/cornell_box_v0.1.13.png)

With path tracing:

```
Sampler "lowdiscrepancy" "integer pixelsamples" [512]
PixelFilter "gaussian" "float xwidth" [2.000000 ] "float ywidth" [2.000000 ]
Integrator "path"
```

![Cornell Box scene rendered via Rust version of PBRT](https://www.janwalter.org/assets/cornell_box_v0.2.0_high.png)

With ambient occlusion:

```
Integrator "ambientocclusion"
```

![Cornell Box scene rendered via Rust using
 AO](https://www.janwalter.org/assets/cornell_box_pbrt_rust_ao.png)

If you look for a more complete Rust implementation:

https://bitbucket.org/abusch/rustracer
