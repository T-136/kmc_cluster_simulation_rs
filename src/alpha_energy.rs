///```let binding_energy = alphas.cn[atom_type_index][in_atom_type_index][cn-1]```
///summed_to_x contains the alphas pre-summed to each coordination number
#[derive(Clone, Debug)]
pub struct Alphas {
    pub cn: [[[f64; 12]; super::NUM_ATOM_TYPES]; super::NUM_ATOM_TYPES],
    pub summed_to_x: [[[f64; 12]; super::NUM_ATOM_TYPES]; super::NUM_ATOM_TYPES],
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
            summed_to_x: alphas_summed_to_x,
        }
    }

    fn sum_up_to_cn(clean_and_dilluted: &[f64; 12], cn: usize) -> f64 {
        clean_and_dilluted[..cn].iter().sum::<f64>()
    }

    //in_metal;cn-1
    fn map_clean_diluted_to_matrix_sum(
        clean_and_dilluted: &[[f64; 12]],
        metal_i: usize,
    ) -> [[f64; 12]; super::NUM_ATOM_TYPES] {
        let mut metal = [[0_f64; 12]; super::NUM_ATOM_TYPES];

        for (nn_m_i, metal_e_list) in clean_and_dilluted.iter().enumerate() {
            for cn_index in 0..metal_e_list.len() {
                metal[nn_m_i][cn_index] =
                    // Alphas::sum_up_to_cn(&clean_and_dilluted[nn_m_i], cn_index + 1);
                clean_and_dilluted[nn_m_i][..(cn_index + 1)]
                    .iter()
                    .sum::<f64>();
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

        if cn_metal == 3 {
            energy += 0.3 ;
        } else if cn_metal == 4 {
            energy += 0.4 ;
        }

        nn_atom_type_count
            .iter()
            .enumerate()
            .for_each(|(metal_type, nn_atom_type_count_num)| {
                energy += self.summed_to_x[atom_type][metal_type][cn_metal - 1]
                    * *nn_atom_type_count_num as f64
                    / cn_metal as f64
            });

        // if cn_metal == 3 || cn_metal == 4 {
        //     return morse_potential(energy);
        // }
        energy

    }

    pub fn e_one_atom<I>(
        &self,
        cn_metal: u8,
        nn_atom_type_count: [u8; super::NUM_ATOM_TYPES],
        nn_nn_atom_type_count: I,
        atom_type: u8,
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
                    energy += self.summed_to_x[atom_type as usize][metal_type]
                        [(cn_metal - 1) as usize]
                        * *nn_atom_type_count_num as f64
                        / cn_metal as f64
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
                        .sum::<u8>()
                        == nn_atom_type_counts.cn_metal
                );
                nn_atom_type_counts
                    .nn_atom_type_count_num_list
                    .iter()
                    .enumerate()
                    .for_each(|(metal_type_index, nn_atom_type_count_num)| {
                        energy += self.cn[nn_atom_type_counts.atom_type][metal_type_index]
                            [(nn_atom_type_counts.cn_metal - 1) as usize]
                            / nn_atom_type_counts.cn_metal as f64
                            * *nn_atom_type_count_num as f64
                    })
            });
        assert!(!energy.is_nan());
        energy
    }
}


#[derive(Clone, Debug)]
pub struct NnData {
    pub cn_metal: u8,
    pub nn_atom_type_count_num_list: [u8; super::NUM_ATOM_TYPES],
    pub atom_type: usize,
}

pub fn e_barrier(prev_e: f64, future_e: f64) -> f64 {
    const offset: f64 = 0.7;
    let e_barr_correction = -(offset * (1. - offset) * (prev_e.abs() * future_e.abs()).sqrt());

    assert!(!e_barr_correction.is_nan());
    (prev_e - e_barr_correction).abs()
    // (prev_e - future_e).abs()
}

fn morse_potential(e_well: f64) -> f64 {
    // https://www.sciencedirect.com/science/article/pii/S0927025622000180
    //          E       A           r_min
    // pd: 7.559e−20	1.576e＋10	2.907e−10
    // pt: 	1.136e−19	1.581e＋10	2.929e−10
    
    // const D_NORMAL: f64 = 1./(2_f64.powf(0.5_f64));
    const D_CUBE_LEN: f64 = std::f64::consts::SQRT_2;
     
    const A: f64 = 1.576e+10;
    const A_EV: f64 = 1.576e+10 / D_MIN;
    const D_MIN: f64 = 2.907e-10;
    const D_E: f64 = 7.559e-20;
    const E: f64 = std::f64::consts::E;


    (-e_well * (1.-E.powf(-A * (-D_MIN * 0.15)).powi(2)) + e_well) 
}

fn lennard_potential() -> f64 {
    // Pd	7.383e−20	2.515e−10	7.559e−20	
    // Pt	1.098e−19	2.541e−10	1.136e−19
    const D_E: f64 = 1.098e-19;
    const D: f64 = 0.9;

    4. * D_E * ((1. /D).powi(12)-(1./D).powi(6))
}


#[cfg(test)]
mod tests {
    use super::*;

    // first one clean
    const alpha_pt: [[f64; 12]; 2] = [
        [
            -2.74574, -0.86615, -0.68165, -0.18156, -0.18146, -0.18136, -0.18126, -0.15201,
            -0.15191, -0.03632, -0.03622, -0.03612,
        ],
        [
            -1.70894, -0.65027, -0.60655, -0.08254, -0.08244, -0.08234, -0.08224, -0.08214,
            -0.08204, -0.08194, -0.08184, -0.05881,
        ],
    ];

    //second one clean
    const alpha_pd: [[f64; 12]; 2] = [
        [
            -1.61277, -0.76781, -0.76267, -0.33190, -0.33180, -0.33170, -0.33160, -0.27197,
            -0.27187, -0.02403, -0.02393, -0.02383,
        ],
        [
            -0.79154, -0.62089, -0.62079, -0.16854, -0.13978, -0.13968, -0.13958, -0.13948,
            -0.13938, -0.10637, -0.09711, -0.09701,
        ],
    ];

    pub const energy_const: [[[f64; 12]; 2]; 2] = [alpha_pt, alpha_pd];

    #[test]
    fn test_e_summed() {

        let mut atom_names: std::collections::HashMap<String, u8> = std::collections::HashMap::new();
        let alpha_file: String = "./Pt_Pd.6.bat".to_string();
        let alphas_inp = super::super::read_alphas(alpha_file, &mut atom_names);
        // let alphas_inp = [alpha_pt, alpha_pd];
        let alphas = Alphas::new(alphas_inp);

        let mut summed_alphas = 0.;
        let cn = 3;
        for cn_i in 0..cn {
            summed_alphas += (alpha_pt[0][cn_i]) / 4. * 2.;
            summed_alphas += alpha_pt[1][cn_i] / 4. * 2.;
        }
        println!(
            "alphas {:?} \n alphas sumed: {:?}",
            alphas.cn, alphas.summed_to_x
        );
        assert!(
            alphas.summed_to_x[0][0][cn - 1] * 2. + alphas.summed_to_x[1][0][cn - 1] * 2.
                - summed_alphas
                < 0.00001
        )
    }

    #[test]
    fn test_show_alphas() {
        let mut atom_names: std::collections::HashMap<String, u8> = std::collections::HashMap::new();
        let alpha_file: String = "./Pt_Pd.6.bat".to_string();
        let alphas_inp = super::super::read_alphas(alpha_file, &mut atom_names);
        // let alphas_inp = [alpha_pt, alpha_pd];
        let alphas = Alphas::new(alphas_inp);

        struct TestPositin  {
        cn_metal: u8,
        nn_atom_type_count: [u8; super::super::NUM_ATOM_TYPES],
        nn_nn_atom_type_count: [NnData],
        }

        
        let nn_nn_atom_bulk = NnData {
            cn_metal: 12,
            nn_atom_type_count_num_list: [12,0],
            atom_type: 0,
        };

        let mut nn_nn_atom_type_count: Vec<NnData> = Vec::new();
        for _ in 0..=7 {
            nn_nn_atom_type_count.push(nn_nn_atom_bulk.clone());
        }

        let start_e = alphas.e_one_atom(11, [11,0], nn_nn_atom_type_count.into_iter(), 0, 0, 0);
        let tst_e = alphas.e_one_atom_tst(4, [4,0], 0, 0, 0);
        let barr = super::e_barrier(start_e-tst_e, start_e-tst_e);

        

        let mut atom_0 = "";
        for ( atom_name, index) in atom_names.iter(){
            if index == &0 {
                atom_0 =  atom_name;
            }

        }


        println!("distance {}", std::f64::consts::SQRT_2 / 2.  );
        println!("morse pot {}", super::morse_potential(barr) );
        println!("move in bulk of {}: {}", atom_0, barr);

    }
}
