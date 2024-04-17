// first one clean
const alpha_pt: [[f64; 12]; 2] = [
    [
        // -2.74574, -0.86615, -0.68165, -0.18156, -0.18146, -0.18136, -0.18126, -0.15201, -0.15191,
        // -0.03632, -0.03622,
        // -0.019765,
        -2.74574, -0.86615, -0.68165, -0.18156, -0.18146, -0.18136, -0.18126, -0.15201, -0.15191,
        -0.03632, -0.03622,
        -0.03612,
        //old:
        // -1.5045, -1.0045, -0.5045, -0.1870, -0.1860, -0.1850, -0.1840, -0.1830, -0.1820, -0.1220,
        // -0.0340, -0.0330,
    ],
    [
        // -1.70894, -0.65027, -0.60655, -0.08254, -0.08244, -0.08234, -0.08224, -0.08214, -0.08204,
        // -0.08194, -0.08184,
        // -0.074097,
        -1.70894, -0.65027, -0.60655, -0.08254, -0.08244, -0.08234, -0.08224, -0.08214, -0.08204,
        -0.08194, -0.08184,
        -0.05881,
        //old:
        // -1.1843, -0.8843, -0.1843, -0.2389, -0.2379, -0.2369, -0.2359, -0.2349, -0.2339, -0.2139,
        // -0.0862, -0.0852,
    ],
];

//second one clean
const alpha_pd: [[f64; 12]; 2] = [
    [
        // -1.61277, -0.76781, -0.76267, -0.33190, -0.33180, -0.33170, -0.33160, -0.27197, -0.27187,
        // -0.02403, -0.02393,
        // 0.41626,
        -1.61277, -0.76781, -0.76267, -0.33190, -0.33180, -0.33170, -0.33160, -0.27197, -0.27187,
        -0.02403, -0.02393,
        -0.02383,
        //old:
        // -0.416261.3012, -0.8012, -0.3012, -0.1102, -0.1092, -0.1082, -0.1072, -0.1062, -0.1052, -0.0436,
        // -0.0068, -0.0058,
    ],
    [
        // -0.79154, -0.62089, -0.62079, -0.16854, -0.13978, -0.13968, -0.13958, -0.13948, -0.13938,
        // -0.10637, -0.09711,
        // -0.08424,
        -0.79154, -0.62089, -0.62079, -0.16854, -0.13978, -0.13968, -0.13958, -0.13948, -0.13938,
        -0.10637, -0.09711,
        -0.09701,
        //old:
        // -1.0094, -1.0094, -0.5094, -0.1497, -0.1487, -0.1403, -0.1393, -0.1352, -0.1342, -0.0906,
        // -0.0896, -0.0871,
    ],
];

pub const energy_const: [[[f64; 12]; super::NUM_ATOM_TYPES]; super::NUM_ATOM_TYPES] =
    [alpha_pt, alpha_pd];
// fn morse_pot(tot_e: f64, dist: f64) -> f64 {
//     let a = 1.;
//     let x: f64 = (-a * dist);
//     // let x: f64 = 0.7;
//     let res = tot_e * (1. - x.exp()).powi(2) - tot_e;
//     res
// }

#[derive(Clone)]
pub struct Alphas {
    //alphas and alphas_summed_to_x are div by cn
    //atom_type_index;in_atom_type_index;cn-1
    pub cn: [[[f64; 12]; super::NUM_ATOM_TYPES]; super::NUM_ATOM_TYPES],
    pub summed_to_x_div_cn: [[[f64; 12]; super::NUM_ATOM_TYPES]; super::NUM_ATOM_TYPES],
}

impl Alphas {
    pub fn new(
        mut alphas_input: [[[f64; 12]; super::NUM_ATOM_TYPES]; super::NUM_ATOM_TYPES],
    ) -> Alphas {
        let alphas_summed_to_x = Alphas::summ_alphas_to_x(&alphas_input);
        //maybe there is a way to increase performance with div_by_cn but understand floating
        //better first
        // Alphas::alphas_div_by_cn(&mut alphas_input);
        Alphas {
            cn: alphas_input,
            summed_to_x_div_cn: alphas_summed_to_x,
        }
    }

    fn sum_up_to_cn(atom_type_index: usize, in_a_type_index: usize, cn: usize) -> f64 {
        (0..(cn))
            // .filter(|cn_i| *cn_i != 0)
            .map(|cn_i| {
                // energy += get_alpha_vector(atom_type, nn_atom_type_count, cn_i)
                energy_const[atom_type_index][in_a_type_index][cn_i] / cn as f64
            })
            .sum::<f64>()
    }

    //in_metal;cn-1
    fn map_clean_diluted_to_matrix_sum(
        clean_and_dilluted: &[[f64; 12]],
        metal_i: usize,
    ) -> [[f64; 12]; super::NUM_ATOM_TYPES] {
        let mut metal = [[0_f64; 12]; super::NUM_ATOM_TYPES];

        for (nn_m_i, metal_e_list) in clean_and_dilluted.iter().enumerate() {
            for cn_index in 0..metal_e_list.len() {
                metal[nn_m_i][cn_index] = Alphas::sum_up_to_cn(metal_i, nn_m_i, cn_index + 1);
            }
        }

        metal
    }

    fn summ_alphas_to_x(
        alphas_input: &[[[f64; 12]; super::NUM_ATOM_TYPES]; super::NUM_ATOM_TYPES],
    ) -> [[[f64; 12]; super::NUM_ATOM_TYPES]; super::NUM_ATOM_TYPES] {
        let mut alphas_summed_to_x: [[[f64; 12]; super::NUM_ATOM_TYPES]; super::NUM_ATOM_TYPES] =
            [[[0_f64; 12]; super::NUM_ATOM_TYPES]; super::NUM_ATOM_TYPES];

        for (i, alpha_metal) in alphas_input.iter().enumerate() {
            alphas_summed_to_x[i] = Alphas::map_clean_diluted_to_matrix_sum(alpha_metal, i);
        }
        alphas_summed_to_x
    }

    fn alphas_div_by_cn(
        alphas_input: &mut [[[f64; 12]; super::NUM_ATOM_TYPES]; super::NUM_ATOM_TYPES],
    ) {
        for i1 in 0..alphas_input.len() {
            for i2 in 0..alphas_input[i1].len() {
                for i3 in 0..alphas_input[i1][i2].len() {
                    alphas_input[i1][i2][i3] = (alphas_input[i1][i2][i3] / (i3 + 1) as f64);
                }
            }
        }
    }

    pub fn e_one_atom_tst(
        &self,
        cn_metal: usize,
        nn_atom_type_count: [u8; super::NUM_ATOM_TYPES],
        atom_type: usize,
        at_supp: u8,
        supp_ee: i64,
    ) -> f64 {
        if cn_metal == 0 {
            return 0.;
        }

        let mut energy = 0.;

        nn_atom_type_count
            .iter()
            .enumerate()
            .for_each(|(metal_type, nn_atom_type_count_num)| {
                energy += self.summed_to_x_div_cn[atom_type - 1][metal_type][cn_metal - 1]
                    * *nn_atom_type_count_num as f64
                // / cn_metal as f64
            });
        energy
    }

    pub fn e_one_atom<I>(
        &self,
        cn_metal: usize,
        nn_atom_type_count: [u8; super::NUM_ATOM_TYPES],
        nn_nn_atom_type_count: I,
        atom_type: usize,
        at_supp: u8,
        supp_ee: i64,
    ) -> f64
    where
        I: Iterator<Item = NnData>,
    {
        let mut energy = 0.;

        if cn_metal != 0 {
            nn_atom_type_count.iter().enumerate().for_each(
                |(metal_type, nn_atom_type_count_num)| {
                    energy += self.summed_to_x_div_cn[atom_type - 1][metal_type][cn_metal - 1]
                        * *nn_atom_type_count_num as f64
                    // / cn_metal as f64
                },
            );
        }
        nn_nn_atom_type_count
            .into_iter()
            .for_each(|nn_atom_type_counts| {
                assert!(
                    nn_atom_type_counts
                        .nn_atom_type_count_num_list
                        .iter()
                        .sum::<u8>() as usize
                        == nn_atom_type_counts.cn_metal
                );
                nn_atom_type_counts
                    .nn_atom_type_count_num_list
                    .iter()
                    .enumerate()
                    .for_each(|(metal_type_index, nn_atom_type_count_num)| {
                        energy += self.cn[nn_atom_type_counts.atom_type - 1][metal_type_index]
                            [nn_atom_type_counts.cn_metal - 1]
                            / nn_atom_type_counts.cn_metal as f64
                            * *nn_atom_type_count_num as f64
                    })
            });
        assert!(!energy.is_nan());
        energy
    }
}

pub struct NnData {
    pub cn_metal: usize,
    pub nn_atom_type_count_num_list: [u8; super::NUM_ATOM_TYPES],
    pub atom_type: usize,
}

pub fn e_barrier(prev_e: f64, future_e: f64) -> f64 {
    const offset: f64 = 0.7;
    // println!("prev_e: {}", prev_e);
    // assert!(!prev_e.is_nan());
    // println!("future_e: {}", future_e);
    // assert!(!future_e.is_nan());
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_e_summed() {
        let alphas_inp = [alpha_pt, alpha_pd];
        let alphas = Alphas::new(alphas_inp);

        let mut summed_alphas = 0.;
        let cn = 3;
        for cn_i in 0..cn {
            summed_alphas += (alpha_pt[0][cn_i]) / 4. * 2.;
            summed_alphas += alpha_pt[1][cn_i] / 4. * 2.;
        }
        assert!(
            alphas.summed_to_x_div_cn[0][0][cn - 1] * 2.
                + alphas.summed_to_x_div_cn[1][0][cn - 1] * 2.
                - summed_alphas
                < 0.00001
        )
    }
}
