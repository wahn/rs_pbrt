===============
Getting Started
===============

Get Rust
========

First `install Rust`_. To
keep it up to date use `rustup`_:

.. code:: shell

          rustup update

Compiling rs_pbrt
=================

Cloning the repository
----------------------

There are two repositories you can get the **Rust source code** from:

1. `GitHub repository`_
2. `Codeberg repository`_

Both should be exactly the same, but reported issues etc. might differ.

GitHub
......

.. code:: shell

          # using SSH
          git clone git@github.com:wahn/rs_pbrt.git
          # or using HTTPS
          git clone https://github.com/wahn/rs_pbrt.git

Codeberg
........

.. code:: shell

          # using SSH
          git clone codeberg@codeberg.org:wahn/rs_pbrt.git
          # or using HTTPS
          git clone https://codeberg.org/wahn/rs_pbrt.git

Use Cargo to compile executables
--------------------------------

.. code:: shell

          # enter repository
          cd rs_pbrt
          # compile without OpenEXR support
          cargo test --release --no-default-features
          # compile executable rs_pbrt (and run it to see options)
          cargo run --release --no-default-features

For a **debug** version compile without the `--release` option.

.. code:: shell

          # compile without OpenEXR support
          cargo test --no-default-features
          # compile executable rs_pbrt (and run it to see options)
          cargo run --no-default-features

For more information about **Cargo**, check out `its
documentation`_.

The executables can be found in either the **release** or the
**debug** target directory:

.. code:: shell

          # release
          ls ./target/release
          # debug
          ls ./target/debug

Create a local copy of the documentation
========================================

.. code:: shell

          # no OpenEXR support
          cargo doc --no-default-features

Use your favourite `web browser`_ to open the local (Rust source code)
documentation:

.. code:: shell

          firefox target/doc/rs_pbrt/index.html

You can also find the `official documentation`_ (of the latest
release) on the `rs_pbrt` web site.

Running the renderer
====================

.. _install Rust: https://www.rust-lang.org/tools/install
.. _rustup: https://github.com/rust-lang-nursery/rustup.rs
.. _GitHub repository: https://github.com/wahn/rs_pbrt
.. _Codeberg repository: https://codeberg.org/wahn/rs_pbrt
.. _its documentation: https://doc.rust-lang.org/cargo
.. _web browser: https://en.wikipedia.org/wiki/Web_browser
.. _official documentation: https://www.rs-pbrt.org/doc/crates/rs_pbrt/index.html
