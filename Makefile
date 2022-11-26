all: message without-exr

clean:
	-rm -f *~ examples/*~ src/*~ examples/*.rs.bk src/*.rs.bk *.exr *.png

clobber: clean
	-rm -fr target

distclean: clobber
	-rm -f *.exr *.png master.zip
	-rm -fr openexr-rs-master

message:
	$(info # if you want to have a git commit tag in your executables reported:")
	$(info # bash: export GIT_DESCRIBE=something")
	$(info # tcsh: setenv GIT_DESCRIBE something")
	$(info # replace 'something' above with the output of:")
	$(info git describe --tags)

doc:
	cargo doc --no-default-features

browse: doc
	firefox ./target/doc/pbrt/index.html

debug: # master.zip
	cargo test --no-default-features
	cargo run --bin parse_blend_file --no-default-features -- --help
	cargo run --bin rs_pbrt --no-default-features -- --help

release:
	cargo test --release
	cargo run --bin parse_blend_file --release -- --help
	cargo run --bin rs_pbrt --release -- --help

without-exr:
	cargo test --release --no-default-features
	cargo run --bin parse_blend_file --release --no-default-features -- --help
	cargo run --bin rs_pbrt --release --no-default-features -- --help
