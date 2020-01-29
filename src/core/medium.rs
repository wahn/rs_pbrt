//! A key operation that **Medium** implementations must perform is to
//! compute the beam transmittance along a given ray.

// std
use std::f32::consts::PI;
use std::sync::Arc;
// pbrt
use crate::core::geometry::{spherical_direction_vec3, vec3_coordinate_system, vec3_dot_vec3};
use crate::core::geometry::{Point2f, Ray, Vector3f};
use crate::core::interaction::MediumInteraction;
use crate::core::pbrt::INV_4_PI;
use crate::core::pbrt::{Float, Spectrum};
use crate::core::sampler::Sampler;
use crate::media::grid::GridDensityMedium;
use crate::media::homogeneous::HomogeneousMedium;

pub const SUBSURFACE_PARAMETER_TABLE: [MeasuredSS; 47] = [
    // From "A Practical Model for Subsurface Light Transport"
    // Jensen, Marschner, Levoy, Hanrahan
    // Proc SIGGRAPH 2001
    MeasuredSS {
        name: "Apple",
        sigma_prime_s: [2.29, 2.39, 1.97],
        sigma_a: [0.0030, 0.0034, 0.046],
    },
    MeasuredSS {
        name: "Chicken1",
        sigma_prime_s: [0.15, 0.21, 0.38],
        sigma_a: [0.015, 0.077, 0.19],
    },
    MeasuredSS {
        name: "Chicken2",
        sigma_prime_s: [0.19, 0.25, 0.32],
        sigma_a: [0.018, 0.088, 0.20],
    },
    MeasuredSS {
        name: "Cream",
        sigma_prime_s: [7.38, 5.47, 3.15],
        sigma_a: [0.0002, 0.0028, 0.0163],
    },
    MeasuredSS {
        name: "Ketchup",
        sigma_prime_s: [0.18, 0.07, 0.03],
        sigma_a: [0.061, 0.97, 1.45],
    },
    MeasuredSS {
        name: "Marble",
        sigma_prime_s: [2.19, 2.62, 3.00],
        sigma_a: [0.0021, 0.0041, 0.0071],
    },
    MeasuredSS {
        name: "Potato",
        sigma_prime_s: [0.68, 0.70, 0.55],
        sigma_a: [0.0024, 0.0090, 0.12],
    },
    MeasuredSS {
        name: "Skimmilk",
        sigma_prime_s: [0.70, 1.22, 1.90],
        sigma_a: [0.0014, 0.0025, 0.0142],
    },
    MeasuredSS {
        name: "Skin1",
        sigma_prime_s: [0.74, 0.88, 1.01],
        sigma_a: [0.032, 0.17, 0.48],
    },
    MeasuredSS {
        name: "Skin2",
        sigma_prime_s: [1.09, 1.59, 1.79],
        sigma_a: [0.013, 0.070, 0.145],
    },
    MeasuredSS {
        name: "Spectralon",
        sigma_prime_s: [11.6, 20.4, 14.9],
        sigma_a: [0.00, 0.00, 0.00],
    },
    MeasuredSS {
        name: "Wholemilk",
        sigma_prime_s: [2.55, 3.21, 3.77],
        sigma_a: [0.0011, 0.0024, 0.014],
    },
    // From "Acquiring Scattering Properties of Participating Media by
    // Dilution",
    // Narasimhan, Gupta, Donner, Ramamoorthi, Nayar, Jensen
    // Proc SIGGRAPH 2006
    MeasuredSS {
        name: "Lowfat Milk",
        sigma_prime_s: [0.89187, 1.5136, 2.532],
        sigma_a: [0.002_875, 0.00575, 0.0115],
    },
    MeasuredSS {
        name: "Reduced Milk",
        sigma_prime_s: [2.4858, 3.1669, 4.5214],
        sigma_a: [0.002_555_6, 0.005_111_1, 0.012_778],
    },
    MeasuredSS {
        name: "Regular Milk",
        sigma_prime_s: [4.5513, 5.8294, 7.136],
        sigma_a: [0.001_533_3, 0.0046, 0.019_933],
    },
    MeasuredSS {
        name: "Espresso",
        sigma_prime_s: [0.72378, 0.84557, 1.0247],
        sigma_a: [4.7984, 6.5751, 8.8493],
    },
    MeasuredSS {
        name: "Mint Mocha Coffee",
        sigma_prime_s: [0.31602, 0.38538, 0.48131],
        sigma_a: [3.772, 5.8228, 7.82],
    },
    MeasuredSS {
        name: "Lowfat Soy Milk",
        sigma_prime_s: [0.30576, 0.34233, 0.61664],
        sigma_a: [0.001_437_5, 0.007_187_5, 0.035_937],
    },
    MeasuredSS {
        name: "Regular Soy Milk",
        sigma_prime_s: [0.59223, 0.73866, 1.4693],
        sigma_a: [0.001_916_7, 0.009_583_3, 0.065_167],
    },
    MeasuredSS {
        name: "Lowfat Chocolate Milk",
        sigma_prime_s: [0.64925, 0.83916, 1.1057],
        sigma_a: [0.0115, 0.0368, 0.1564],
    },
    MeasuredSS {
        name: "Regular Chocolate Milk",
        sigma_prime_s: [1.4585, 2.1289, 2.9527],
        sigma_a: [0.010_063, 0.043_125, 0.14375],
    },
    MeasuredSS {
        name: "Coke",
        sigma_prime_s: [8.9053e-05, 8.372e-05, 0.0],
        sigma_a: [0.10014, 0.16503, 0.2468],
    },
    MeasuredSS {
        name: "Pepsi",
        sigma_prime_s: [6.1697e-05, 4.2564e-05, 0.0],
        sigma_a: [0.091_641, 0.14158, 0.20729],
    },
    MeasuredSS {
        name: "Sprite",
        sigma_prime_s: [6.0306e-06, 6.4139e-06, 6.5504e-06],
        sigma_a: [0.001_886, 0.001_830_8, 0.002_002_5],
    },
    MeasuredSS {
        name: "Gatorade",
        sigma_prime_s: [0.002_457_4, 0.003_007, 0.003_732_5],
        sigma_a: [0.024_794, 0.019_289, 0.008_878],
    },
    MeasuredSS {
        name: "Chardonnay",
        sigma_prime_s: [1.7982e-05, 1.3758e-05, 1.2023e-05],
        sigma_a: [0.010_782, 0.011_855, 0.023_997],
    },
    MeasuredSS {
        name: "White Zinfandel",
        sigma_prime_s: [1.7501e-05, 1.9069e-05, 1.288e-05],
        sigma_a: [0.012_072, 0.016_184, 0.019_843],
    },
    MeasuredSS {
        name: "Merlot",
        sigma_prime_s: [2.1129e-05, 0.0, 0.0],
        sigma_a: [0.11632, 0.25191, 0.29434],
    },
    MeasuredSS {
        name: "Budweiser Beer",
        sigma_prime_s: [2.4356e-05, 2.4079e-05, 1.0564e-05],
        sigma_a: [0.011_492, 0.024_911, 0.057_786],
    },
    MeasuredSS {
        name: "Coors Light Beer",
        sigma_prime_s: [5.0922e-05, 4.301e-05, 0.0],
        sigma_a: [0.006_164, 0.013_984, 0.034_983],
    },
    MeasuredSS {
        name: "Clorox",
        sigma_prime_s: [0.002_403_5, 0.003_137_3, 0.003_991],
        sigma_a: [0.003_354_2, 0.014_892, 0.026_297],
    },
    MeasuredSS {
        name: "Apple Juice",
        sigma_prime_s: [0.000_136_12, 0.000_158_36, 0.000_227],
        sigma_a: [0.012_957, 0.023_741, 0.052_184],
    },
    MeasuredSS {
        name: "Cranberry Juice",
        sigma_prime_s: [0.000_104_02, 0.000_116_46, 7.8139e-05],
        sigma_a: [0.039_437, 0.094_223, 0.12426],
    },
    MeasuredSS {
        name: "Grape Juice",
        sigma_prime_s: [5.382e-05, 0.0, 0.0],
        sigma_a: [0.10404, 0.23958, 0.29325],
    },
    MeasuredSS {
        name: "Ruby Grapefruit Juice",
        sigma_prime_s: [0.011_002, 0.010_927, 0.011_036],
        sigma_a: [0.085_867, 0.18314, 0.25262],
    },
    MeasuredSS {
        name: "White Grapefruit Juice",
        sigma_prime_s: [0.22826, 0.23998, 0.32748],
        sigma_a: [0.0138, 0.018_831, 0.056_781],
    },
    MeasuredSS {
        name: "Shampoo",
        sigma_prime_s: [0.000_717_6, 0.000_830_3, 0.000_901_6],
        sigma_a: [0.014_107, 0.045_693, 0.061_717],
    },
    MeasuredSS {
        name: "Strawberry Shampoo",
        sigma_prime_s: [0.000_156_71, 0.000_159_47, 1.518e-05],
        sigma_a: [0.01449, 0.05796, 0.075_823],
    },
    MeasuredSS {
        name: "Head & Shoulders Shampoo",
        sigma_prime_s: [0.023_805, 0.028_804, 0.034_306],
        sigma_a: [0.084_621, 0.15688, 0.20365],
    },
    MeasuredSS {
        name: "Lemon Tea Powder",
        sigma_prime_s: [0.040_224, 0.045_264, 0.051_081],
        sigma_a: [2.4288, 4.5757, 7.2127],
    },
    MeasuredSS {
        name: "Orange Powder",
        sigma_prime_s: [0.000_156_17, 0.000_174_82, 0.000_176_2],
        sigma_a: [0.001_449, 0.003_441, 0.007_863],
    },
    MeasuredSS {
        name: "Pink Lemonade Powder",
        sigma_prime_s: [0.000_121_03, 0.000_130_73, 0.000_125_28],
        sigma_a: [0.001_165, 0.002_366, 0.003_195],
    },
    MeasuredSS {
        name: "Cappuccino Powder",
        sigma_prime_s: [1.8436, 2.5851, 2.1662],
        sigma_a: [35.844, 49.547, 61.084],
    },
    MeasuredSS {
        name: "Salt Powder",
        sigma_prime_s: [0.027_333, 0.032_451, 0.031_979],
        sigma_a: [0.28415, 0.3257, 0.34148],
    },
    MeasuredSS {
        name: "Sugar Powder",
        sigma_prime_s: [0.000_222_72, 0.000_255_13, 0.000_271],
        sigma_a: [0.012_638, 0.031_051, 0.050_124],
    },
    MeasuredSS {
        name: "Suisse Mocha Powder",
        sigma_prime_s: [2.7979, 3.5452, 4.3365],
        sigma_a: [17.502, 27.004, 35.433],
    },
    MeasuredSS {
        name: "Pacific Ocean Surface Water",
        sigma_prime_s: [0.000_176_4, 0.000_320_95, 0.000_196_17],
        sigma_a: [0.031_845, 0.031_324, 0.030_147],
    },
];

pub struct MeasuredSS {
    pub name: &'static str,
    pub sigma_prime_s: [Float; 3],
    pub sigma_a: [Float; 3],
}

pub struct NoMedium {}

pub enum Medium {
    Empty(NoMedium),
    GridDensity(GridDensityMedium),
    Homogeneous(HomogeneousMedium),
}

impl Medium {
    pub fn tr(&self, r_world: &Ray, sampler: &mut Box<Sampler>) -> Spectrum {
        match self {
            Medium::Empty(_medium) => Spectrum::default(),
            Medium::GridDensity(medium) => medium.tr(r_world, sampler),
            Medium::Homogeneous(medium) => medium.tr(r_world, sampler),
        }
    }
    pub fn sample(
        &self,
        r_world: &Ray,
        sampler: &mut Box<Sampler>,
    ) -> (Spectrum, Option<MediumInteraction>) {
        match self {
            Medium::Empty(_medium) => (Spectrum::default(), None),
            Medium::GridDensity(medium) => medium.sample(r_world, sampler),
            Medium::Homogeneous(medium) => medium.sample(r_world, sampler),
        }
    }
}

pub struct HenyeyGreenstein {
    pub g: Float,
}

impl HenyeyGreenstein {
    pub fn p(&self, wo: &Vector3f, wi: &Vector3f) -> Float {
        // TODO: ProfilePhase _(Prof::PhaseFuncEvaluation);
        phase_hg(vec3_dot_vec3(wo, wi), self.g)
    }
    pub fn sample_p(&self, wo: &Vector3f, wi: &mut Vector3f, u: &Point2f) -> Float {
        // TODO: ProfilePhase _(Prof::PhaseFuncSampling);
        // compute $\cos \theta$ for Henyey--Greenstein sample
        let cos_theta: Float;
        if self.g.abs() < 1e-3 as Float {
            cos_theta = 1.0 as Float - 2.0 as Float * u[0];
        } else {
            let sqr_term: Float = (1.0 as Float - self.g * self.g)
                / (1.0 as Float - self.g + 2.0 as Float * self.g * u[0]);
            cos_theta =
                (1.0 as Float + self.g * self.g - sqr_term * sqr_term) / (2.0 as Float * self.g);
        }
        // compute direction _wi_ for Henyey--Greenstein sample
        let sin_theta: Float = (0.0 as Float)
            .max(1.0 as Float - cos_theta * cos_theta)
            .sqrt();
        let phi: Float = 2.0 as Float * PI * u[1];
        let mut v1: Vector3f = Vector3f::default();
        let mut v2: Vector3f = Vector3f::default();
        vec3_coordinate_system(wo, &mut v1, &mut v2);
        *wi = spherical_direction_vec3(sin_theta, cos_theta, phi, &v1, &v2, &(-*wo));
        phase_hg(-cos_theta, self.g)
    }
}

#[derive(Default, Clone)]
pub struct MediumInterface {
    pub inside: Option<Arc<Medium>>,
    pub outside: Option<Arc<Medium>>,
}

impl MediumInterface {
    pub fn new(inside: Option<Arc<Medium>>, outside: Option<Arc<Medium>>) -> Self {
        MediumInterface { inside, outside }
    }
    pub fn is_medium_transition(&self) -> bool {
        if let Some(ref inside) = self.inside {
            // self.inside == Some
            if let Some(ref outside) = self.outside {
                // self.outside == Some
                let pi = &*inside as *const _ as *const usize;
                let po = &*outside as *const _ as *const usize;
                pi != po
            } else {
                // self.outside == None
                true
            }
        } else {
            // self.inside == None
            if let Some(ref _outside) = self.outside {
                // self.outside == Some
                true
            } else {
                // self.outside == None
                false
            }
        }
    }
    pub fn get_inside(&self) -> Option<Arc<Medium>> {
        if let Some(ref inside) = self.inside {
            Some(inside.clone())
        } else {
            None
        }
    }
    pub fn get_outside(&self) -> Option<Arc<Medium>> {
        if let Some(ref outside) = self.outside {
            Some(outside.clone())
        } else {
            None
        }
    }
}

pub fn get_medium_scattering_properties(
    name: &String,
    sigma_a: &mut Spectrum,
    sigma_prime_s: &mut Spectrum,
) -> bool {
    if *name == "" {
        return false;
    }
    for mss in SUBSURFACE_PARAMETER_TABLE.iter() {
        if name == mss.name {
            *sigma_a = Spectrum::from_rgb(&mss.sigma_a);
            *sigma_prime_s = Spectrum::from_rgb(&mss.sigma_prime_s);
            return true;
        }
    }
    false
}

pub fn phase_hg(cos_theta: Float, g: Float) -> Float {
    let denom: Float = 1.0 as Float + g * g + 2.0 as Float * g * cos_theta;
    INV_4_PI * (1.0 as Float - g * g) / (denom * denom.sqrt())
}
