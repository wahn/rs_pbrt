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
	$(info # replace 'something' abobe with the output of:")
	$(info git describe --tags)

doc:
	cargo doc --no-default-features

browse: doc
	firefox ./target/doc/pbrt/index.html

debug: # master.zip
	cargo test --no-default-features
	cargo run --no-default-features -- --help

release: # master.zip
	cargo test --release
	cargo run --release -- --help

without-exr:
	cargo test --release --no-default-features
	cargo run --release --no-default-features -- --help

# master.zip:
# 	wget https://github.com/cessen/openexr-rs/archive/master.zip
# 	unzip master.zip

examples: without-exr
	./target/release/examples/parse_ass_file -i ./assets/ass/cornell_box.ass
	./target/release/examples/parse_blend_file ./assets/blend/suzanne_integrator_test_2_79.blend
