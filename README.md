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

If you look for a more complete Rust implementation:

https://bitbucket.org/abusch/rustracer
