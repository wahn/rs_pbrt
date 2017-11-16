// pbrt
use core::pbrt::Float;
use core::rng::ONE_MINUS_EPSILON;
use core::rng::Rng;
use core::sampling::shuffle;
use geometry::{Point2i, Point2f};

// see lowdiscrepancy.h

pub const PRIME_TABLE_SIZE: u16 = 1000_u16;

pub const PRIME_SUMS: [u32; PRIME_TABLE_SIZE as usize] =
    [0, 2, 5, 10, 17, 28, 41, 58, 77, 100, 129, 160, 197, 238, 281, 328, 381, 440, 501, 568, 639,
     712, 791, 874, 963, 1060, 1161, 1264, 1371, 1480, 1593, 1720, 1851, 1988, 2127, 2276, 2427,
     2584, 2747, 2914, 3087, 3266, 3447, 3638, 3831, 4028, 4227, 4438, 4661, 4888, 5117, 5350,
     5589, 5830, 6081, 6338, 6601, 6870, 7141, 7418, 7699, 7982, 8275, 8582, 8893, 9206, 9523,
     9854, 10191, 10538, 10887, 11240, 11599, 11966, 12339, 12718, 13101, 13490, 13887, 14288,
     14697, 15116, 15537, 15968, 16401, 16840, 17283, 17732, 18189, 18650, 19113, 19580, 20059,
     20546, 21037, 21536, 22039, 22548, 23069, 23592, 24133, 24680, 25237, 25800, 26369, 26940,
     27517, 28104, 28697, 29296, 29897, 30504, 31117, 31734, 32353, 32984, 33625, 34268, 34915,
     35568, 36227, 36888, 37561, 38238, 38921, 39612, 40313, 41022, 41741, 42468, 43201, 43940,
     44683, 45434, 46191, 46952, 47721, 48494, 49281, 50078, 50887, 51698, 52519, 53342, 54169,
     54998, 55837, 56690, 57547, 58406, 59269, 60146, 61027, 61910, 62797, 63704, 64615, 65534,
     66463, 67400, 68341, 69288, 70241, 71208, 72179, 73156, 74139, 75130, 76127, 77136, 78149,
     79168, 80189, 81220, 82253, 83292, 84341, 85392, 86453, 87516, 88585, 89672, 90763, 91856,
     92953, 94056, 95165, 96282, 97405, 98534, 99685, 100838, 102001, 103172, 104353, 105540,
     106733, 107934, 109147, 110364, 111587, 112816, 114047, 115284, 116533, 117792, 119069,
     120348, 121631, 122920, 124211, 125508, 126809, 128112, 129419, 130738, 132059, 133386,
     134747, 136114, 137487, 138868, 140267, 141676, 143099, 144526, 145955, 147388, 148827,
     150274, 151725, 153178, 154637, 156108, 157589, 159072, 160559, 162048, 163541, 165040,
     166551, 168074, 169605, 171148, 172697, 174250, 175809, 177376, 178947, 180526, 182109,
     183706, 185307, 186914, 188523, 190136, 191755, 193376, 195003, 196640, 198297, 199960,
     201627, 203296, 204989, 206686, 208385, 210094, 211815, 213538, 215271, 217012, 218759,
     220512, 222271, 224048, 225831, 227618, 229407, 231208, 233019, 234842, 236673, 238520,
     240381, 242248, 244119, 245992, 247869, 249748, 251637, 253538, 255445, 257358, 259289,
     261222, 263171, 265122, 267095, 269074, 271061, 273054, 275051, 277050, 279053, 281064,
     283081, 285108, 287137, 289176, 291229, 293292, 295361, 297442, 299525, 301612, 303701,
     305800, 307911, 310024, 312153, 314284, 316421, 318562, 320705, 322858, 325019, 327198,
     329401, 331608, 333821, 336042, 338279, 340518, 342761, 345012, 347279, 349548, 351821,
     354102, 356389, 358682, 360979, 363288, 365599, 367932, 370271, 372612, 374959, 377310,
     379667, 382038, 384415, 386796, 389179, 391568, 393961, 396360, 398771, 401188, 403611,
     406048, 408489, 410936, 413395, 415862, 418335, 420812, 423315, 425836, 428367, 430906,
     433449, 435998, 438549, 441106, 443685, 446276, 448869, 451478, 454095, 456716, 459349,
     461996, 464653, 467312, 469975, 472646, 475323, 478006, 480693, 483382, 486075, 488774,
     491481, 494192, 496905, 499624, 502353, 505084, 507825, 510574, 513327, 516094, 518871,
     521660, 524451, 527248, 530049, 532852, 535671, 538504, 541341, 544184, 547035, 549892,
     552753, 555632, 558519, 561416, 564319, 567228, 570145, 573072, 576011, 578964, 581921,
     584884, 587853, 590824, 593823, 596824, 599835, 602854, 605877, 608914, 611955, 615004,
     618065, 621132, 624211, 627294, 630383, 633492, 636611, 639732, 642869, 646032, 649199,
     652368, 655549, 658736, 661927, 665130, 668339, 671556, 674777, 678006, 681257, 684510,
     687767, 691026, 694297, 697596, 700897, 704204, 707517, 710836, 714159, 717488, 720819,
     724162, 727509, 730868, 734229, 737600, 740973, 744362, 747753, 751160, 754573, 758006,
     761455, 764912, 768373, 771836, 775303, 778772, 782263, 785762, 789273, 792790, 796317,
     799846, 803379, 806918, 810459, 814006, 817563, 821122, 824693, 828274, 831857, 835450,
     839057, 842670, 846287, 849910, 853541, 857178, 860821, 864480, 868151, 871824, 875501,
     879192, 882889, 886590, 890299, 894018, 897745, 901478, 905217, 908978, 912745, 916514,
     920293, 924086, 927883, 931686, 935507, 939330, 943163, 947010, 950861, 954714, 958577,
     962454, 966335, 970224, 974131, 978042, 981959, 985878, 989801, 993730, 997661, 1001604,
     1005551, 1009518, 1013507, 1017508, 1021511, 1025518, 1029531, 1033550, 1037571, 1041598,
     1045647, 1049698, 1053755, 1057828, 1061907, 1065998, 1070091, 1074190, 1078301, 1082428,
     1086557, 1090690, 1094829, 1098982, 1103139, 1107298, 1111475, 1115676, 1119887, 1124104,
     1128323, 1132552, 1136783, 1141024, 1145267, 1149520, 1153779, 1158040, 1162311, 1166584,
     1170867, 1175156, 1179453, 1183780, 1188117, 1192456, 1196805, 1201162, 1205525, 1209898,
     1214289, 1218686, 1223095, 1227516, 1231939, 1236380, 1240827, 1245278, 1249735, 1254198,
     1258679, 1263162, 1267655, 1272162, 1276675, 1281192, 1285711, 1290234, 1294781, 1299330,
     1303891, 1308458, 1313041, 1317632, 1322229, 1326832, 1331453, 1336090, 1340729, 1345372,
     1350021, 1354672, 1359329, 1363992, 1368665, 1373344, 1378035, 1382738, 1387459, 1392182,
     1396911, 1401644, 1406395, 1411154, 1415937, 1420724, 1425513, 1430306, 1435105, 1439906,
     1444719, 1449536, 1454367, 1459228, 1464099, 1468976, 1473865, 1478768, 1483677, 1488596,
     1493527, 1498460, 1503397, 1508340, 1513291, 1518248, 1523215, 1528184, 1533157, 1538144,
     1543137, 1548136, 1553139, 1558148, 1563159, 1568180, 1573203, 1578242, 1583293, 1588352,
     1593429, 1598510, 1603597, 1608696, 1613797, 1618904, 1624017, 1629136, 1634283, 1639436,
     1644603, 1649774, 1654953, 1660142, 1665339, 1670548, 1675775, 1681006, 1686239, 1691476,
     1696737, 1702010, 1707289, 1712570, 1717867, 1723170, 1728479, 1733802, 1739135, 1744482,
     1749833, 1755214, 1760601, 1765994, 1771393, 1776800, 1782213, 1787630, 1793049, 1798480,
     1803917, 1809358, 1814801, 1820250, 1825721, 1831198, 1836677, 1842160, 1847661, 1853164,
     1858671, 1864190, 1869711, 1875238, 1880769, 1886326, 1891889, 1897458, 1903031, 1908612,
     1914203, 1919826, 1925465, 1931106, 1936753, 1942404, 1948057, 1953714, 1959373, 1965042,
     1970725, 1976414, 1982107, 1987808, 1993519, 1999236, 2004973, 2010714, 2016457, 2022206,
     2027985, 2033768, 2039559, 2045360, 2051167, 2056980, 2062801, 2068628, 2074467, 2080310,
     2086159, 2092010, 2097867, 2103728, 2109595, 2115464, 2121343, 2127224, 2133121, 2139024,
     2144947, 2150874, 2156813, 2162766, 2168747, 2174734, 2180741, 2186752, 2192781, 2198818,
     2204861, 2210908, 2216961, 2223028, 2229101, 2235180, 2241269, 2247360, 2253461, 2259574,
     2265695, 2271826, 2277959, 2284102, 2290253, 2296416, 2302589, 2308786, 2314985, 2321188,
     2327399, 2333616, 2339837, 2346066, 2352313, 2358570, 2364833, 2371102, 2377373, 2383650,
     2389937, 2396236, 2402537, 2408848, 2415165, 2421488, 2427817, 2434154, 2440497, 2446850,
     2453209, 2459570, 2465937, 2472310, 2478689, 2485078, 2491475, 2497896, 2504323, 2510772,
     2517223, 2523692, 2530165, 2536646, 2543137, 2549658, 2556187, 2562734, 2569285, 2575838,
     2582401, 2588970, 2595541, 2602118, 2608699, 2615298, 2621905, 2628524, 2635161, 2641814,
     2648473, 2655134, 2661807, 2668486, 2675175, 2681866, 2688567, 2695270, 2701979, 2708698,
     2715431, 2722168, 2728929, 2735692, 2742471, 2749252, 2756043, 2762836, 2769639, 2776462,
     2783289, 2790118, 2796951, 2803792, 2810649, 2817512, 2824381, 2831252, 2838135, 2845034,
     2851941, 2858852, 2865769, 2872716, 2879665, 2886624, 2893585, 2900552, 2907523, 2914500,
     2921483, 2928474, 2935471, 2942472, 2949485, 2956504, 2963531, 2970570, 2977613, 2984670,
     2991739, 2998818, 3005921, 3013030, 3020151, 3027278, 3034407, 3041558, 3048717, 3055894,
     3063081, 3070274, 3077481, 3084692, 3091905, 3099124, 3106353, 3113590, 3120833, 3128080,
     3135333, 3142616, 3149913, 3157220, 3164529, 3171850, 3179181, 3186514, 3193863, 3201214,
     3208583, 3215976, 3223387, 3230804, 3238237, 3245688, 3253145, 3260604, 3268081, 3275562,
     3283049, 3290538, 3298037, 3305544, 3313061, 3320584, 3328113, 3335650, 3343191, 3350738,
     3358287, 3365846, 3373407, 3380980, 3388557, 3396140, 3403729, 3411320, 3418923, 3426530,
     3434151, 3441790, 3449433, 3457082, 3464751, 3472424, 3480105, 3487792, 3495483, 3503182,
     3510885, 3518602, 3526325, 3534052, 3541793, 3549546, 3557303, 3565062, 3572851, 3580644,
     3588461, 3596284, 3604113, 3611954, 3619807, 3627674, 3635547, 3643424, 3651303, 3659186,
     3667087, 3674994];

/// The bits of an integer quantity can be efficiently reversed with a
/// series of logical bit operations.
pub fn reverse_bits_32(n: u32) -> u32 {
    let mut n = (n << 16) | (n >> 16);
    n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8);
    n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4);
    n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2);
    n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1);
    n
}

/// The bits of a 64-bit value can be reversed by reversing the two
/// 32-bit components individually and then interchanging them.
pub fn reverse_bits_64(n: u64) -> u64 {
    let n0: u64 = reverse_bits_32(n as u32) as u64;
    let n1: u64 = reverse_bits_32((n >> 32) as u32) as u64;
    (n0 << 32) | n1
}


/// Compute the inverse of the radical inverse function.
pub fn inverse_radical_inverse(base: u8, inverse: u64, n_digits: u64) -> u64 {
    let mut inverse: u64 = inverse;
    let mut index: u64 = 0_64;
    for _i in 0..n_digits {
        let digit: u64 = inverse % base as u64;
        inverse /= base as u64;
        index = index * base as u64 + digit;
    }
    index
}

/// Takes a generator matrix *c*, a number of 1D samples to generate
/// *n*, and stores the corresponding samples in memory at the
/// location pointed to by *p*.
pub fn gray_code_sample_1d(c: [u32; 32], n: u32, scramble: u32, p: &mut [Float]) {
    let mut v: u32 = scramble;
    for i in 0..n as usize {
        // 1/2^32
        p[i] = (v as Float * 2.3283064365386963e-10 as Float).min(ONE_MINUS_EPSILON);
        v ^= c[(i + 1).trailing_zeros() as usize];
    }
}

/// Takes two generator matrices *c0* and *c1*, a number of 2D samples
/// to generate *n*, and stores the corresponding samples in memory at
/// the location pointed to by *p*.
pub fn gray_code_sample_2d(c0: &[u32], c1: &[u32], n: u32, scramble: &Point2i, p: &mut [Point2f]) {
    let mut v: [u32; 2] = [scramble.x as u32, scramble.y as u32];
    for i in 0..n as usize {
        p[i].x = (v[0] as Float * 2.3283064365386963e-10 as Float).min(ONE_MINUS_EPSILON);
        p[i].y = (v[1] as Float * 2.3283064365386963e-10 as Float).min(ONE_MINUS_EPSILON);
        v[0] ^= c0[(i + 1).trailing_zeros() as usize];
        v[1] ^= c1[(i + 1).trailing_zeros() as usize];
    }
}

/// Generates a number of scrambled 1D sample values using the Gray
/// code-based sampling machinery.
pub fn van_der_corput(n_samples_per_pixel_sample: i32,
                      n_pixel_samples: i32,
                      samples: &mut [Float],
                      rng: &mut Rng) {
    let scramble: u32 = rng.uniform_uint32();
    let c_van_der_corput: [u32; 32] = [0x80000000, 0x40000000, 0x20000000, 0x10000000, 0x8000000,
                                       0x4000000, 0x2000000, 0x1000000, 0x800000, 0x400000,
                                       0x200000, 0x100000, 0x80000, 0x40000, 0x20000, 0x10000,
                                       0x8000, 0x4000, 0x2000, 0x1000, 0x800, 0x400, 0x200, 0x100,
                                       0x80, 0x40, 0x20, 0x10, 0x8, 0x4, 0x2, 0x1];
    let total_samples: i32 = n_samples_per_pixel_sample * n_pixel_samples;
    gray_code_sample_1d(c_van_der_corput, total_samples as u32, scramble, samples);
    // randomly shuffle 1D sample points
    for i in 0..n_pixel_samples as usize {
        shuffle(&mut samples[(i * n_samples_per_pixel_sample as usize)..],
                n_samples_per_pixel_sample,
                1,
                rng);
    }
    shuffle(&mut samples[..],
            n_pixel_samples,
            n_samples_per_pixel_sample,
            rng);
}

/// Similar to *van_der_corput()*, but uses two generator matrices to
/// generate the first two dimensions of Sobol' points.
pub fn sobol_2d(n_samples_per_pixel_sample: i32,
                n_pixel_samples: i32,
                samples: &mut [Point2f],
                rng: &mut Rng) {
    let x: i32 = rng.uniform_uint32() as i32;
    let y: i32 = rng.uniform_uint32() as i32;
    let scramble: Point2i = Point2i { x: x, y: y };
    // define 2D Sobol$'$ generator matrices _c_sobol[2]_
    let c_sobol: [[u32; 32]; 2] =
        [[0x80000000_u32,
          0x40000000,
          0x20000000,
          0x10000000,
          0x8000000,
          0x4000000,
          0x2000000,
          0x1000000,
          0x800000,
          0x400000,
          0x200000,
          0x100000,
          0x80000,
          0x40000,
          0x20000,
          0x10000,
          0x8000,
          0x4000,
          0x2000,
          0x1000,
          0x800,
          0x400,
          0x200,
          0x100,
          0x80,
          0x40,
          0x20,
          0x10,
          0x8,
          0x4,
          0x2,
          0x1],
         [0x80000000, 0xc0000000, 0xa0000000, 0xf0000000, 0x88000000, 0xcc000000, 0xaa000000,
          0xff000000, 0x80800000, 0xc0c00000, 0xa0a00000, 0xf0f00000, 0x88880000, 0xcccc0000,
          0xaaaa0000, 0xffff0000, 0x80008000, 0xc000c000, 0xa000a000, 0xf000f000, 0x88008800,
          0xcc00cc00, 0xaa00aa00, 0xff00ff00, 0x80808080, 0xc0c0c0c0, 0xa0a0a0a0, 0xf0f0f0f0,
          0x88888888, 0xcccccccc, 0xaaaaaaaa, 0xffffffff]];
    gray_code_sample_2d(&c_sobol[0],
                        &c_sobol[1],
                        (n_samples_per_pixel_sample * n_pixel_samples) as u32,
                        &scramble,
                        &mut samples[..]);
    for i in 0..n_pixel_samples as usize {
        shuffle(&mut samples[(i * n_samples_per_pixel_sample as usize)..],
                n_samples_per_pixel_sample,
                1,
                rng);
    }
    shuffle(&mut samples[..],
            n_pixel_samples,
            n_samples_per_pixel_sample,
            rng);
}

// see lowdiscrepancy.cpp

/// Once we have an appropriate prime number, use it to compute the
/// radical inverse.
pub fn radical_inverse_specialized(base: u16, a: u64) -> Float {
    let inv_base: Float = 1.0 as Float / base as Float;
    let mut reversed_digits: u64 = 0_u64;
    let mut inv_base_n: Float = 1.0 as Float;
    let mut a: u64 = a; // shadowing input parameter
    while a != 0_u64 {
        let next: u64 = a / base as u64;
        let digit: u64 = a - next * base as u64;
        reversed_digits = reversed_digits * base as u64 + digit;
        inv_base_n *= inv_base;
        a = next;
    }
    assert!(reversed_digits as Float * inv_base_n < 1.00001 as Float);
    (reversed_digits as Float * inv_base_n).min(ONE_MINUS_EPSILON)
}

/// Compute the radical inverse, but put each pixel through the
/// permutation table for the given base. Deal with a special case
/// that can arise.
pub fn scrambled_radical_inverse_specialized(base: u16, perm: Vec<u16>, a: u64) -> Float {
    let inv_base: Float = 1.0 as Float / base as Float;
    let mut reversed_digits: u64 = 0_u64;
    let mut inv_base_n: Float = 1.0 as Float;
    let mut a: u64 = a; // shadowing input parameter
    while a != 0_u64 {
        let next: u64 = a / base as u64;
        let digit: u64 = a - next * base as u64;
        assert!(perm[digit as usize] < base);
        reversed_digits = reversed_digits * base as u64 + perm[digit as usize] as u64;
        inv_base_n *= inv_base;
        a = next;
    }
    assert!((inv_base_n *
             (reversed_digits as Float * inv_base_n * perm[0] as Float /
              (1.0 as Float - inv_base))) < 1.00001 as Float);
    (inv_base_n *
     (reversed_digits as Float * inv_base_n * perm[0] as Float / (1.0 as Float - inv_base)))
            .min(ONE_MINUS_EPSILON)
}

/// Map to an appropriate prime number and delegate to another
/// function to compute the radical inverse.
pub fn radical_inverse(base_index: u16, a: u64) -> Float {
    match base_index {
        0 => {
            // TODO: #ifndef PBRT_HAVE_HEX_FP_CONSTANTS
            // 0x1p-64 = (2.0 as Float).powi(-64 as i32)
            return reverse_bits_64(a) as Float * (2.0 as Float).powi(-64 as i32);
        }
        1 => {
            return radical_inverse_specialized(3_u16, a);
        }
        2 => {
            return radical_inverse_specialized(5_u16, a);
        }
        3 => {
            return radical_inverse_specialized(7_u16, a);
        }
        4 => {
            return radical_inverse_specialized(11_u16, a);
        }
        5 => {
            return radical_inverse_specialized(13_u16, a);
        }
        6 => {
            return radical_inverse_specialized(17_u16, a);
        }
        // WORK
        _ => {
            panic!("TODO: radical_inverse({:?}, {:?})", base_index, a);
        }
    };
}

/// Compute the radical inverse, but put each pixel through the
/// permutation table for the given base.
pub fn scrambled_radical_inverse(base_index: u16, a: u64, perm: Vec<u16>) -> Float {
    match base_index {
        0 => {
            return scrambled_radical_inverse_specialized(2_u16, perm, a);
        }
        1 => {
            return scrambled_radical_inverse_specialized(3_u16, perm, a);
        }
        2 => {
            return scrambled_radical_inverse_specialized(5_u16, perm, a);
        }
        3 => {
            return scrambled_radical_inverse_specialized(7_u16, perm, a);
        }
        // WORK
        _ => {
            panic!("TODO: scrambled_radical_inverse({:?}, {:?})", base_index, a);
        }
    };
}
