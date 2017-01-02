all: release

clean:
	-rm -f *~ examples/*~ src/*~ examples/*.rs.bk src/*.rs.bk

clobber: clean
	-rm -fr target

doc:
	cargo doc

browse: doc
	firefox ./target/doc/pbrt/index.html

debug:
	cargo test

release:
	cargo test --release

examples: release
	./target/release/examples/geometry_bounds2_unit_cube
	./target/release/examples/geometry_bounds3_unit_cube
	./target/release/examples/geometry_length
	./target/release/examples/geometry_length_squared
	./target/release/examples/geometry_normal3_null
	./target/release/examples/geometry_point2_origin
	./target/release/examples/geometry_point3_origin
	./target/release/examples/geometry_ray_creation
	./target/release/examples/geometry_vector2_null
	./target/release/examples/geometry_vector3_null
	./target/release/examples/pbrt
	./target/release/examples/shapes_sphere_default
	./target/release/examples/shapes_sphere_new
	./target/release/examples/transform_matrix4x4_identity
	./target/release/examples/transform_matrix4x4_new
	./target/release/examples/transform_matrix4x4_transpose
	./target/release/examples/transform_transform_identity
	./target/release/examples/transform_transform_look_at
	./target/release/examples/transform_transform_new
	./target/release/examples/transform_transform_rotate
	./target/release/examples/transform_transform_rotate_x
	./target/release/examples/transform_transform_rotate_y
	./target/release/examples/transform_transform_rotate_z
	./target/release/examples/transform_transform_scale
	./target/release/examples/transform_transform_translate
