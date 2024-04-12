// first one clean
const alpha_pt: [[f64; 13]; 2] = [
    [
        0., -1.5045, -1.0045, -0.5045, -0.1870, -0.1860, -0.1850, -0.1840, -0.1830, -0.1820,
        -0.1220, -0.0340, -0.0330,
    ],
    [
        0., -1.1843, -0.8843, -0.1843, -0.2389, -0.2379, -0.2369, -0.2359, -0.2349, -0.2339,
        -0.2139, -0.0862, -0.0852,
    ],
];

//second one clean
const alpha_pd: [[f64; 13]; 2] = [
    [
        0., -1.3012, -0.8012, -0.3012, -0.1102, -0.1092, -0.1082, -0.1072, -0.1062, -0.1052,
        -0.0436, -0.0068, -0.0058,
    ],
    [
        0., -1.0094, -1.0094, -0.5094, -0.1497, -0.1487, -0.1403, -0.1393, -0.1352, -0.1342,
        -0.0906, -0.0896, -0.0871,
    ],
];

// fn morse_pot(tot_e: f64, dist: f64) -> f64 {
//     let a = 1.;
//     let x: f64 = (-a * dist);
//     // let x: f64 = 0.7;
//     let res = tot_e * (1. - x.exp()).powi(2) - tot_e;
//     res
// }

pub fn e_barrier(prev_e: f64, future_e: f64) -> f64 {
    const offset: f64 = 0.7;
    // println!("prev_e: {}", prev_e);
    assert!(!prev_e.is_nan());
    // println!("future_e: {}", future_e);
    assert!(!future_e.is_nan());
    let e_barr_correction = -(offset * (1. - offset) * (prev_e.abs() * future_e.abs()).sqrt());
    // println!("res: {}", res);

    assert!(!e_barr_correction.is_nan());
    (prev_e - e_barr_correction).abs()
}

fn sum_alphas(
    atom_type: usize,
    nn_atom_type_count: [u8; super::NUM_ATOM_TYPES],
    cn_metal: usize,
) -> f64 {
    // if could be taken out for summing over central atom
    if cn_metal == 0 {
        return 0.;
    }
    let e = if atom_type == 1 {
        alpha_pt[0][cn_metal]
            + (alpha_pt[1][cn_metal] - alpha_pt[0][cn_metal]) / cn_metal as f64
                * nn_atom_type_count[1] as f64
        // + alpha_pt[1][cn_metal]
        // + (alpha_pt[0][cn_metal] - alpha_pt[1][cn_metal]) / cn_metal as f64
        //     * nn_atom_type_count[0] as f64
    } else if atom_type == 2 {
        alpha_pd[1][cn_metal]
            + (alpha_pd[0][cn_metal] - alpha_pd[1][cn_metal]) / cn_metal as f64
                * nn_atom_type_count[0] as f64
        // + alpha_pd[1][cn_metal]
        // + (alpha_pd[0][cn_metal] - alpha_pd[1][cn_metal]) / cn_metal as f64
        //     * nn_atom_type_count[0] as f64
        // alpha_pd[1][cn_metal] * nn_atom_type_count[1] as f64
        //     + alpha_pd[0][cn_metal] * nn_atom_type_count[0] as f64
    } else {
        println!("atom type: {}", atom_type);
        panic!("wtf");
    };
    // println!("{} {}", e, cn_metal);
    e
}

pub fn e_one_atom<I>(
    cn_metal_range: (usize, usize),
    nn_atom_type_count: [u8; super::NUM_ATOM_TYPES],
    nn_nn_atom_type_count: I,
    atom_type: usize,
    at_supp: u8,
    supp_ee: i64,
) -> f64
where
    I: Iterator<Item = (usize, [u8; super::NUM_ATOM_TYPES], usize)>,
{
    assert!(nn_atom_type_count.iter().sum::<u8>() as usize == cn_metal_range.1 - cn_metal_range.0);
    let mut energy = 0.;

    for cn_i in cn_metal_range.0..=cn_metal_range.1 {
        energy += sum_alphas(atom_type, nn_atom_type_count, cn_i)
    }
    for nn_atom_type_count in nn_nn_atom_type_count {
        // println!("atom neig type: {}", nn_atom_type_count.2);

        assert!(nn_atom_type_count.1.iter().sum::<u8>() as usize == nn_atom_type_count.0);
        energy += sum_alphas(
            nn_atom_type_count.2,
            nn_atom_type_count.1,
            nn_atom_type_count.0,
        )
    }
    assert!(!energy.is_nan());
    energy
}

// pub fn simd_sum(values: &[f32]) -> f32 {
//     let chunks = values.chunks_exact(LANES);
//     let remainder = chunks.remainder();
//
//     let sum = chunks.fold([0.0f32; LANES], |mut acc, chunk| {
//         let chunk: [f32; LANES] = chunk.try_into().unwrap();
//         for i in 0..LANES {
//             acc[i] += chunk[i];
//         }
//         acc
//     });
//
//     let remainder: f32 = remainder.iter().copied().sum();
//
//     let mut reduced = 0.0f32;
//     for i in 0..LANES {
//         reduced += sum[i];
//     }
//     reduced + remainder
// }
