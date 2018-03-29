// pbrt
use core::pbrt::{Float, Spectrum};

pub const SUBSURFACE_PARAMETER_TABLE: [MeasuredSS; 12] = [
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
];

pub struct MeasuredSS {
    pub name: &'static str,
    pub sigma_prime_s: [Float; 3],
    pub sigma_a: [Float; 3],
}

pub fn get_medium_scattering_properties(
    name: &String,
    sigma_a: &mut Spectrum,
    sigma_prime_s: &mut Spectrum,
) -> bool {
    // WORK
    false
}
