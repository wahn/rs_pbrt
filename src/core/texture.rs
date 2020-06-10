//! **Texture** is a template class parameterized by return type of
//! its evaluation function. This design makes it possible to reuse
//! almost all of the code among textures that return different
//! types. PBRT currently uses only **Float** and **Spectrum**
//! textures.

// std
use std::f32::consts::PI;
// pbrt
use crate::core::geometry::{spherical_phi, spherical_theta, vec3_dot_vec3};
use crate::core::geometry::{Point2f, Point3f, Vector2f, Vector3f, XYEnum};
use crate::core::interaction::SurfaceInteraction;
use crate::core::pbrt::Float;
use crate::core::pbrt::{clamp_t, lerp, log_2};
use crate::core::pbrt::{INV_2_PI, INV_PI};
use crate::core::transform::Transform;

// see texture.h

// Perlin Noise Data
pub const NOISE_PERM_SIZE: usize = 256;
pub const NOISE_PERM: [u8; 2 * NOISE_PERM_SIZE] = [
    151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, 69,
    142, // remainder of the noise permutation table
    8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203,
    117, 35, 11, 32, 57, 177, 33, 88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74,
    165, 71, 134, 139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220,
    105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132,
    187, 208, 89, 18, 169, 200, 196, 135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3,
    64, 52, 217, 226, 250, 124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59,
    227, 47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70,
    221, 153, 101, 155, 167, 43, 172, 9, 129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232,
    178, 185, 112, 104, 218, 246, 97, 228, 251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162,
    241, 81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204,
    176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141,
    128, 195, 78, 66, 215, 61, 156, 180, 151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194,
    233, 7, 225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 247, 120, 234,
    75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33, 88, 237, 149, 56, 87, 174,
    20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166, 77, 146, 158, 231, 83,
    111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54, 65, 25,
    63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196, 135, 130, 116, 188,
    159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123, 5, 202, 38, 147,
    118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 170,
    213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9, 129, 22, 39, 253,
    19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97, 228, 251, 34, 242, 193,
    238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31,
    181, 199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 93,
    222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180,
];

pub enum TextureMapping2D {
    UV(UVMapping2D),
    Spherical(SphericalMapping2D),
    Cylindrical(CylindricalMapping2D),
    Planar(PlanarMapping2D),
}

impl TextureMapping2D {
    pub fn map(
        &self,
        si: &SurfaceInteraction,
        dstdx: &mut Vector2f,
        dstdy: &mut Vector2f,
    ) -> Point2f {
        match self {
            TextureMapping2D::UV(texturemapping2d) => texturemapping2d.map(si, dstdx, dstdy),
            TextureMapping2D::Spherical(texturemapping2d) => texturemapping2d.map(si, dstdx, dstdy),
            TextureMapping2D::Cylindrical(texturemapping2d) => {
                texturemapping2d.map(si, dstdx, dstdy)
            }
            TextureMapping2D::Planar(texturemapping2d) => texturemapping2d.map(si, dstdx, dstdy),
        }
    }
}

pub enum TextureMapping3D {
    Identity(IdentityMapping3D),
}

impl TextureMapping3D {
    pub fn map(
        &self,
        si: &SurfaceInteraction,
        dpdx: &mut Vector3f,
        dpdy: &mut Vector3f,
    ) -> Point3f {
        match self {
            TextureMapping3D::Identity(texturemapping3d) => texturemapping3d.map(si, dpdx, dpdy),
        }
    }
}

#[derive(Debug, Default, Copy, Clone)]
pub struct UVMapping2D {
    pub su: Float,
    pub sv: Float,
    pub du: Float,
    pub dv: Float,
}

impl UVMapping2D {
    pub fn map(
        &self,
        si: &SurfaceInteraction,
        dstdx: &mut Vector2f,
        dstdy: &mut Vector2f,
    ) -> Point2f {
        // compute texture differentials for 2D identity mapping
        *dstdx = Vector2f {
            x: si.dudx.get() * self.su,
            y: si.dvdx.get() * self.sv,
        };
        *dstdy = Vector2f {
            x: si.dudy.get() * self.su,
            y: si.dvdy.get() * self.sv,
        };
        Point2f {
            x: si.uv[XYEnum::X] * self.su + self.du,
            y: si.uv[XYEnum::Y] * self.sv + self.dv,
        }
    }
}

#[derive(Debug, Default, Copy, Clone)]
pub struct SphericalMapping2D {
    sphere: Point2f,
    pub world_to_texture: Transform,
}

impl SphericalMapping2D {
    pub fn new(world_to_texture: Transform) -> Self {
        SphericalMapping2D {
            sphere: Point2f::default(),
            world_to_texture,
        }
    }
    pub fn sphere(&self, p: &Point3f) -> Point2f {
        let vec3f: Vector3f =
            (self.world_to_texture.transform_point(p) - Point3f::default()).normalize();
        let theta: Float = spherical_theta(&vec3f);
        let phi: Float = spherical_phi(&vec3f);
        Point2f {
            x: theta * INV_PI,
            y: phi * INV_2_PI,
        }
    }
}

impl SphericalMapping2D {
    pub fn map(
        &self,
        si: &SurfaceInteraction,
        dstdx: &mut Vector2f,
        dstdy: &mut Vector2f,
    ) -> Point2f {
        let st: Point2f = self.sphere(&si.p);
        // compute texture coordinate differentials for sphere $(u,v)$ mapping
        let delta: Float = 0.1;
        let st_delta_x: Point2f = self.sphere(&(si.p + si.dpdx.get() * delta));
        *dstdx = (st_delta_x - st) / delta;
        let st_delta_y: Point2f = self.sphere(&(si.p + si.dpdy.get() * delta));
        *dstdy = (st_delta_y - st) / delta;
        // handle sphere mapping discontinuity for coordinate differentials
        if (*dstdx)[XYEnum::Y] > 0.5 as Float {
            (*dstdx)[XYEnum::Y] = 1.0 as Float - (*dstdx)[XYEnum::Y];
        } else if (*dstdx)[XYEnum::Y] < -0.5 as Float {
            (*dstdx)[XYEnum::Y] = -((*dstdx)[XYEnum::Y] + 1.0 as Float);
        }
        if (*dstdy)[XYEnum::Y] > 0.5 as Float {
            (*dstdy)[XYEnum::Y] = 1.0 as Float - (*dstdy)[XYEnum::Y];
        } else if (*dstdy)[XYEnum::Y] < -0.5 as Float {
            (*dstdy)[XYEnum::Y] = -((*dstdy)[XYEnum::Y] + 1.0 as Float);
        }
        st
    }
}

#[derive(Debug, Default, Copy, Clone)]
pub struct CylindricalMapping2D {
    pub world_to_texture: Transform,
}

impl CylindricalMapping2D {
    pub fn new(world_to_texture: Transform) -> Self {
        CylindricalMapping2D { world_to_texture }
    }
    pub fn cylinder(&self, p: &Point3f) -> Point2f {
        let vec3f: Vector3f =
            (self.world_to_texture.transform_point(p) - Point3f::default()).normalize();
        Point2f {
            x: PI + vec3f.y.atan2(vec3f.x) * INV_2_PI,
            y: vec3f.z,
        }
    }
}

impl CylindricalMapping2D {
    pub fn map(
        &self,
        si: &SurfaceInteraction,
        dstdx: &mut Vector2f,
        dstdy: &mut Vector2f,
    ) -> Point2f {
        let st: Point2f = self.cylinder(&si.p);
        // compute texture coordinate differentials for cylinder $(u,v)$ mapping
        let delta: Float = 0.01;
        let st_delta_x: Point2f = self.cylinder(&(si.p + si.dpdx.get() * delta));
        *dstdx = (st_delta_x - st) / delta;
        if (*dstdx)[XYEnum::Y] > 0.5 as Float {
            (*dstdx)[XYEnum::Y] = 1.0 as Float - (*dstdx)[XYEnum::Y];
        } else if (*dstdx)[XYEnum::Y] < -0.5 as Float {
            (*dstdx)[XYEnum::Y] = -((*dstdx)[XYEnum::Y] + 1.0 as Float);
        }
        let st_delta_y: Point2f = self.cylinder(&(si.p + si.dpdy.get() * delta));
        *dstdy = (st_delta_y - st) / delta;
        if (*dstdy)[XYEnum::Y] > 0.5 as Float {
            (*dstdy)[XYEnum::Y] = 1.0 as Float - (*dstdy)[XYEnum::Y];
        } else if (*dstdy)[XYEnum::Y] < -0.5 as Float {
            (*dstdy)[XYEnum::Y] = -((*dstdy)[XYEnum::Y] + 1.0 as Float);
        }
        st
    }
}

#[derive(Debug, Default, Copy, Clone)]
pub struct PlanarMapping2D {
    pub vs: Vector3f,
    pub vt: Vector3f,
    pub ds: Float,
    pub dt: Float,
}

impl PlanarMapping2D {
    pub fn map(
        &self,
        si: &SurfaceInteraction,
        dstdx: &mut Vector2f,
        dstdy: &mut Vector2f,
    ) -> Point2f {
        let vec: Vector3f = Vector3f {
            x: si.p.x,
            y: si.p.y,
            z: si.p.z,
        };
        *dstdx = Vector2f {
            x: vec3_dot_vec3(&si.dpdx.get(), &self.vs),
            y: vec3_dot_vec3(&si.dpdx.get(), &self.vt),
        };
        *dstdy = Vector2f {
            x: vec3_dot_vec3(&si.dpdy.get(), &self.vs),
            y: vec3_dot_vec3(&si.dpdy.get(), &self.vt),
        };
        Point2f {
            x: self.ds + vec3_dot_vec3(&vec, &self.vs),
            y: self.dt + vec3_dot_vec3(&vec, &self.vt),
        }
    }
}

#[derive(Debug, Default, Copy, Clone)]
pub struct IdentityMapping3D {
    pub world_to_texture: Transform,
}

impl IdentityMapping3D {
    pub fn new(world_to_texture: Transform) -> Self {
        IdentityMapping3D { world_to_texture }
    }
    pub fn get_world_to_texture(&self) -> Transform {
        self.world_to_texture
    }
    // TextureMapping3D
    pub fn map(
        &self,
        si: &SurfaceInteraction,
        dpdx: &mut Vector3f,
        dpdy: &mut Vector3f,
    ) -> Point3f {
        let world_to_texture = self.get_world_to_texture();
        *dpdx = world_to_texture.transform_vector(&si.dpdx.get());
        *dpdy = world_to_texture.transform_vector(&si.dpdy.get());
        world_to_texture.transform_point(&si.p)
    }
}

pub trait Texture<T> {
    fn evaluate(&self, si: &SurfaceInteraction) -> T;
}

pub fn smooth_step(min: Float, max: Float, value: Float) -> Float {
    let v: Float = clamp_t((value - min) / (max - min), 0.0 as Float, 1.0 as Float);
    v * v * (-2.0 as Float * v + 3.0 as Float)
}

pub fn noise_flt(x: Float, y: Float, z: Float) -> Float {
    // compute noise cell coordinates and offsets
    let mut ix: i32 = x.floor() as i32;
    let mut iy: i32 = y.floor() as i32;
    let mut iz: i32 = z.floor() as i32;
    let dx: Float = x - ix as Float;
    let dy: Float = y - iy as Float;
    let dz: Float = z - iz as Float;
    // compute gradient weights
    ix &= NOISE_PERM_SIZE as i32 - 1;
    iy &= NOISE_PERM_SIZE as i32 - 1;
    iz &= NOISE_PERM_SIZE as i32 - 1;
    let w000: Float = grad(ix, iy, iz, dx, dy, dz);
    let w100: Float = grad(ix + 1, iy, iz, dx - 1.0 as Float, dy, dz);
    let w010: Float = grad(ix, iy + 1, iz, dx, dy - 1.0 as Float, dz);
    let w110: Float = grad(ix + 1, iy + 1, iz, dx - 1.0 as Float, dy - 1.0 as Float, dz);
    let w001: Float = grad(ix, iy, iz + 1, dx, dy, dz - 1.0 as Float);
    let w101: Float = grad(ix + 1, iy, iz + 1, dx - 1.0 as Float, dy, dz - 1.0 as Float);
    let w011: Float = grad(ix, iy + 1, iz + 1, dx, dy - 1.0 as Float, dz - 1.0 as Float);
    let w111: Float = grad(
        ix + 1,
        iy + 1,
        iz + 1,
        dx - 1.0 as Float,
        dy - 1.0 as Float,
        dz - 1.0 as Float,
    );
    // compute trilinear interpolation of weights
    let wx: Float = noise_weight(dx);
    let wy: Float = noise_weight(dy);
    let wz: Float = noise_weight(dz);
    let x00: Float = lerp(wx, w000, w100);
    let x10: Float = lerp(wx, w010, w110);
    let x01: Float = lerp(wx, w001, w101);
    let x11: Float = lerp(wx, w011, w111);
    let y0: Float = lerp(wy, x00, x10);
    let y1: Float = lerp(wy, x01, x11);
    let ret: Float = lerp(wz, y0, y1);
    ret
}

pub fn noise_pnt3(p: &Point3f) -> Float {
    noise_flt(p.x, p.y, p.z)
}

pub fn grad(x: i32, y: i32, z: i32, dx: Float, dy: Float, dz: Float) -> Float {
    let mut h: u8 =
        NOISE_PERM[NOISE_PERM[NOISE_PERM[x as usize] as usize + y as usize] as usize + z as usize];
    h &= 15_u8;
    let u = if h < 8_u8 || h == 12_u8 || h == 13_u8 {
        dx
    } else {
        dy
    };
    let v = if h < 4_u8 || h == 12_u8 || h == 13_u8 {
        dy
    } else {
        dz
    };
    let ret_u = if h & 1_u8 > 0_u8 { -u } else { u };
    let ret_v = if h & 2_u8 > 0_u8 { -v } else { v };
    ret_u + ret_v
}

pub fn noise_weight(t: Float) -> Float {
    let t3: Float = t * t * t;
    let t4: Float = t3 * t;
    6.0 as Float * t4 * t - 15.0 as Float * t4 + 10.0 as Float * t3
}

pub fn fbm(p: &Point3f, dpdx: &Vector3f, dpdy: &Vector3f, omega: Float, max_octaves: i32) -> Float {
    // compute number of octaves for antialiased FBm
    let len2: Float = dpdx.length_squared().max(dpdy.length_squared());
    let n: Float = clamp_t(
        -1.0 as Float - 0.5 as Float * log_2(len2),
        0.0 as Float,
        max_octaves as Float,
    );
    let n_int: i32 = n.floor() as i32;
    // compute sum of octaves of noise for FBm
    let mut sum: Float = 0.0;
    let mut lambda: Float = 1.0;
    let mut o: Float = 1.0;
    for _i in 0..n_int {
        sum += o * noise_pnt3(&(*p * lambda));
        lambda *= 1.99 as Float;
        o *= omega;
    }
    let n_partial: Float = n - n_int as Float;
    sum += o * smooth_step(0.3 as Float, 0.7 as Float, n_partial) * noise_pnt3(&(*p * lambda));
    sum
}

pub fn turbulence(
    p: &Point3f,
    dpdx: &Vector3f,
    dpdy: &Vector3f,
    omega: Float,
    max_octaves: i32,
) -> Float {
    // compute number of octaves for antialiased FBm
    let len2: Float = dpdx.length_squared().max(dpdy.length_squared());
    let n: Float = clamp_t(
        -1.0 as Float - 0.5 as Float * log_2(len2),
        0.0 as Float,
        max_octaves as Float,
    );
    let n_int: usize = n.floor() as usize;
    // compute sum of octaves of noise for turbulence
    let mut sum: Float = 0.0;
    let mut lambda: Float = 1.0;
    let mut o: Float = 1.0;
    for _i in 0..n_int {
        sum += o * noise_pnt3(&(*p * lambda)).abs();
        lambda *= 1.99 as Float;
        o *= omega;
    }
    // account for contributions of clamped octaves in turbulence
    let n_partial: Float = n - n_int as Float;
    sum += o * lerp(
        smooth_step(0.3 as Float, 0.7 as Float, n_partial),
        0.2,
        noise_pnt3(&(*p * lambda)).abs(),
    );
    for _i in n_int..max_octaves as usize {
        sum += o * 0.2 as Float;
        o *= omega;
    }
    sum
}

pub fn lanczos(x: Float, tau: Float) -> Float {
    let mut x: Float = x;
    x = x.abs();
    if x < 1e-5 as Float {
        return 1.0 as Float;
    }
    if x > 1.0 as Float {
        return 0.0 as Float;
    }
    x *= PI;
    let s: Float = (x * tau).sin() / (x * tau);
    let lanczos: Float = x.sin() / x;
    s * lanczos
}
