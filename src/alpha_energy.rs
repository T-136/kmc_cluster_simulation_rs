// first one clean
const alpha_pt: [[f64; 12]; 2] = [
    [
        -1.5045, -1.0045, -0.5045, -0.1870, -0.1860, -0.1850, -0.1840, -0.1830, -0.1820, -0.1220,
        -0.0340, -0.0330,
    ],
    [
        -1.1843, -0.8843, -0.1843, -0.2389, -0.2379, -0.2369, -0.2359, -0.2349, -0.2339, -0.2139,
        -0.0862, -0.0852,
    ],
];

//second one clean
const alpha_pd: [[f64; 12]; 2] = [
    [
        -1.3012, -0.8012, -0.3012, -0.1102, -0.1092, -0.1082, -0.1072, -0.1062, -0.1052, -0.0436,
        -0.0068, -0.0058,
    ],
    [
        -1.0094, -1.0094, -0.5094, -0.1497, -0.1487, -0.1403, -0.1393, -0.1352, -0.1342, -0.0906,
        -0.0896, -0.0871,
    ],
];

const energy_const: [[[f64; 12]; super::NUM_ATOM_TYPES]; super::NUM_ATOM_TYPES] =
    [alpha_pt, alpha_pd];
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

fn get_alpha_vector(
    atom_type: usize,
    nn_atom_type_count: [u8; super::NUM_ATOM_TYPES],
    cn_metal: usize,
    // energy_const: [[[f64; 13]; METALS_N]; METALS_N],
) -> f64 {
    // if could be taken out for summing over central atom
    if cn_metal == 0 {
        return 0.;
    }
    let atom_type_index = atom_type - 1;
    energy_const[atom_type_index][1][cn_metal] / cn_metal as f64 * nn_atom_type_count[1] as f64
}

pub fn e_one_atom_tst(
    cn_metal_range: (usize, usize),
    nn_atom_type_count: [u8; super::NUM_ATOM_TYPES],
    atom_type: usize,
    at_supp: u8,
    supp_ee: i64,
) -> f64 {
    // assert!(nn_atom_type_count.iter().sum::<u8>() as usize == cn_metal_range.1 - cn_metal_range.0);
    let mut energy = 0.;

    nn_atom_type_count
        .iter()
        .enumerate()
        .for_each(|(metal_type, nn_atom_type_count_num)| {
            energy += (0..(cn_metal_range.1))
                // .filter(|cn_i| *cn_i != 0)
                .map(|cn_i| {
                    // energy += get_alpha_vector(atom_type, nn_atom_type_count, cn_i)
                    energy_const[atom_type - 1][metal_type][cn_i] / cn_metal_range.1 as f64
                        * *nn_atom_type_count_num as f64
                })
                .sum::<f64>();
        });
    assert!(!energy.is_nan());
    energy
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
    // assert!(nn_atom_type_count.iter().sum::<u8>() as usize == cn_metal_range.1 - cn_metal_range.0);
    let mut energy = 0.;

    nn_atom_type_count
        .into_iter()
        .enumerate()
        .for_each(|(metal_type, nn_atom_type_count_num)| {
            energy += (0..(cn_metal_range.1))
                // .filter(|cn_i| *cn_i != 0)
                .map(|cn_i| {
                    // energy += get_alpha_vector(atom_type, nn_atom_type_count, cn_i)
                    energy_const[atom_type - 1][metal_type][cn_i] / cn_metal_range.1 as f64
                        * nn_atom_type_count_num as f64
                })
                .sum::<f64>();
        });
    for nn_atom_type_counts in nn_nn_atom_type_count {
        assert!(nn_atom_type_counts.1.iter().sum::<u8>() as usize == nn_atom_type_counts.0);
        for (metal_type, nn_atom_type_count_num) in nn_atom_type_counts.1.iter().enumerate() {
            energy += energy_const[nn_atom_type_counts.2 - 1][metal_type][nn_atom_type_counts.0 - 1]
                / nn_atom_type_counts.0 as f64
                * *nn_atom_type_count_num as f64
        }
        // energy += get_alpha_vector(
        //     nn_atom_type_counts.2,
        //     nn_atom_type_counts.1,
        //     nn_atom_type_counts.0,
        // )
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
