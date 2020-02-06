all: without-exr

clean:
	-rm -f *~ examples/*~ src/*~ examples/*.rs.bk src/*.rs.bk *.exr *.png

clobber: clean
	-rm -fr target

distclean: clobber
	-rm -f *.exr *.png master.zip
	-rm -fr openexr-rs-master

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
	./target/release/examples/api_make_camera
	./target/release/examples/cameras_perspective_generate_ray_differential
	./target/release/examples/core_lowdiscrepancy_radical_inverse
	./target/release/examples/core_next_float_down
	./target/release/examples/core_next_float_up
	./target/release/examples/core_quadratic
	./target/release/examples/core_read_float_file
	./target/release/examples/core_rng_set_sequence
	./target/release/examples/filters_create_box_filter
	./target/release/examples/filters_create_gaussian_filter
	./target/release/examples/filters_create_mitchell_filter
	./target/release/examples/filters_create_sinc_filter
	./target/release/examples/filters_create_triangle_filter
	./target/release/examples/geometry_bounds2_unit_cube
	./target/release/examples/geometry_bounds3_unit_cube
	./target/release/examples/geometry_length
	./target/release/examples/geometry_length_squared
	./target/release/examples/geometry_normal3_null
	./target/release/examples/geometry_point2_origin
	./target/release/examples/geometry_point3_origin
	./target/release/examples/geometry_ray_creation
	./target/release/examples/geometry_spherical_direction_vec3
	./target/release/examples/geometry_vector2_null
	./target/release/examples/geometry_vector3_null
	./target/release/examples/integrators_sspm_hash
	./target/release/examples/lights_diffuse_area_light_new
	./target/release/examples/lights_distant_light_new
	./target/release/examples/lights_infinite_area_light_new
	./target/release/examples/lights_point_light_new
	./target/release/examples/parse_ass_file -i ./assets/ass/cornell_box.ass
	./target/release/examples/parse_blend_file ./assets/blend/suzanne_integrator_test_2_79.blend
	./target/release/examples/pbrt_spheres_differentials_texfilt
	./target/release/examples/pbrt_teapot_area_light
	./target/release/examples/shapes_cylinder_create_cylinder_shape
	./target/release/examples/shapes_disk_create_disk_shape
	./target/release/examples/shapes_sphere_create_sphere_shape
	./target/release/examples/shapes_sphere_intersect
	./target/release/examples/shapes_sphere_world_bound
	./target/release/examples/shapes_triangle_create_triangle_mesh
	./target/release/examples/shapes_triangle_intersect
	./target/release/examples/shapes_triangle_world_bound
	./target/release/examples/transform_matrix4x4_identity
	./target/release/examples/transform_matrix4x4_new
	./target/release/examples/transform_matrix4x4_transpose
	./target/release/examples/transform_quaternion_default
	./target/release/examples/transform_transform_identity
	./target/release/examples/transform_transform_look_at
	./target/release/examples/transform_transform_new
	./target/release/examples/transform_transform_point_with_error
	./target/release/examples/transform_transform_ray_with_error
	./target/release/examples/transform_transform_rotate
	./target/release/examples/transform_transform_rotate_x
	./target/release/examples/transform_transform_rotate_y
	./target/release/examples/transform_transform_rotate_z
	./target/release/examples/transform_transform_scale
	./target/release/examples/transform_transform_translate
	./target/release/examples/transform_transform_vector_with_error
