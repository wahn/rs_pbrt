// pbrt
use core::pbrt::{Float, Spectrum};

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
        sigma_a: [0.002875, 0.00575, 0.0115],
    },
    MeasuredSS {
        name: "Reduced Milk",
        sigma_prime_s: [2.4858, 3.1669, 4.5214],
        sigma_a: [0.0025556, 0.0051111, 0.012778],
    },
    MeasuredSS {
        name: "Regular Milk",
        sigma_prime_s: [4.5513, 5.8294, 7.136],
        sigma_a: [0.0015333, 0.0046, 0.019933],
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
        sigma_a: [0.0014375, 0.0071875, 0.035937],
    },
    MeasuredSS {
        name: "Regular Soy Milk",
        sigma_prime_s: [0.59223, 0.73866, 1.4693],
        sigma_a: [0.0019167, 0.0095833, 0.065167],
    },
    MeasuredSS {
        name: "Lowfat Chocolate Milk",
        sigma_prime_s: [0.64925, 0.83916, 1.1057],
        sigma_a: [0.0115, 0.0368, 0.1564],
    },
    MeasuredSS {
        name: "Regular Chocolate Milk",
        sigma_prime_s: [1.4585, 2.1289, 2.9527],
        sigma_a: [0.010063, 0.043125, 0.14375],
    },
    MeasuredSS {
        name: "Coke",
        sigma_prime_s: [8.9053e-05, 8.372e-05, 0.0],
        sigma_a: [0.10014, 0.16503, 0.2468],
    },
    MeasuredSS {
        name: "Pepsi",
        sigma_prime_s: [6.1697e-05, 4.2564e-05, 0.0],
        sigma_a: [0.091641, 0.14158, 0.20729],
    },
    MeasuredSS {
        name: "Sprite",
        sigma_prime_s: [6.0306e-06, 6.4139e-06, 6.5504e-06],
        sigma_a: [0.001886, 0.0018308, 0.0020025],
    },
    MeasuredSS {
        name: "Gatorade",
        sigma_prime_s: [0.0024574, 0.003007, 0.0037325],
        sigma_a: [0.024794, 0.019289, 0.008878],
    },
    MeasuredSS {
        name: "Chardonnay",
        sigma_prime_s: [1.7982e-05, 1.3758e-05, 1.2023e-05],
        sigma_a: [0.010782, 0.011855, 0.023997],
    },
    MeasuredSS {
        name: "White Zinfandel",
        sigma_prime_s: [1.7501e-05, 1.9069e-05, 1.288e-05],
        sigma_a: [0.012072, 0.016184, 0.019843],
    },
    MeasuredSS {
        name: "Merlot",
        sigma_prime_s: [2.1129e-05, 0.0, 0.0],
        sigma_a: [0.11632, 0.25191, 0.29434],
    },
    MeasuredSS {
        name: "Budweiser Beer",
        sigma_prime_s: [2.4356e-05, 2.4079e-05, 1.0564e-05],
        sigma_a: [0.011492, 0.024911, 0.057786],
    },
    MeasuredSS {
        name: "Coors Light Beer",
        sigma_prime_s: [5.0922e-05, 4.301e-05, 0.0],
        sigma_a: [0.006164, 0.013984, 0.034983],
    },
    MeasuredSS {
        name: "Clorox",
        sigma_prime_s: [0.0024035, 0.0031373, 0.003991],
        sigma_a: [0.0033542, 0.014892, 0.026297],
    },
    MeasuredSS {
        name: "Apple Juice",
        sigma_prime_s: [0.00013612, 0.00015836, 0.000227],
        sigma_a: [0.012957, 0.023741, 0.052184],
    },
    MeasuredSS {
        name: "Cranberry Juice",
        sigma_prime_s: [0.00010402, 0.00011646, 7.8139e-05],
        sigma_a: [0.039437, 0.094223, 0.12426],
    },
    MeasuredSS {
        name: "Grape Juice",
        sigma_prime_s: [5.382e-05, 0.0, 0.0],
        sigma_a: [0.10404, 0.23958, 0.29325],
    },
    MeasuredSS {
        name: "Ruby Grapefruit Juice",
        sigma_prime_s: [0.011002, 0.010927, 0.011036],
        sigma_a: [0.085867, 0.18314, 0.25262],
    },
    MeasuredSS {
        name: "White Grapefruit Juice",
        sigma_prime_s: [0.22826, 0.23998, 0.32748],
        sigma_a: [0.0138, 0.018831, 0.056781],
    },
    MeasuredSS {
        name: "Shampoo",
        sigma_prime_s: [0.0007176, 0.0008303, 0.0009016],
        sigma_a: [0.014107, 0.045693, 0.061717],
    },
    MeasuredSS {
        name: "Strawberry Shampoo",
        sigma_prime_s: [0.00015671, 0.00015947, 1.518e-05],
        sigma_a: [0.01449, 0.05796, 0.075823],
    },
    MeasuredSS {
        name: "Head & Shoulders Shampoo",
        sigma_prime_s: [0.023805, 0.028804, 0.034306],
        sigma_a: [0.084621, 0.15688, 0.20365],
    },
    MeasuredSS {
        name: "Lemon Tea Powder",
        sigma_prime_s: [0.040224, 0.045264, 0.051081],
        sigma_a: [2.4288, 4.5757, 7.2127],
    },
    MeasuredSS {
        name: "Orange Powder",
        sigma_prime_s: [0.00015617, 0.00017482, 0.0001762],
        sigma_a: [0.001449, 0.003441, 0.007863],
    },
    MeasuredSS {
        name: "Pink Lemonade Powder",
        sigma_prime_s: [0.00012103, 0.00013073, 0.00012528],
        sigma_a: [0.001165, 0.002366, 0.003195],
    },
    MeasuredSS {
        name: "Cappuccino Powder",
        sigma_prime_s: [1.8436, 2.5851, 2.1662],
        sigma_a: [35.844, 49.547, 61.084],
    },
    MeasuredSS {
        name: "Salt Powder",
        sigma_prime_s: [0.027333, 0.032451, 0.031979],
        sigma_a: [0.28415, 0.3257, 0.34148],
    },
    MeasuredSS {
        name: "Sugar Powder",
        sigma_prime_s: [0.00022272, 0.00025513, 0.000271],
        sigma_a: [0.012638, 0.031051, 0.050124],
    },
    MeasuredSS {
        name: "Suisse Mocha Powder",
        sigma_prime_s: [2.7979, 3.5452, 4.3365],
        sigma_a: [17.502, 27.004, 35.433],
    },
    MeasuredSS {
        name: "Pacific Ocean Surface Water",
        sigma_prime_s: [0.0001764, 0.00032095, 0.00019617],
        sigma_a: [0.031845, 0.031324, 0.030147],
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
