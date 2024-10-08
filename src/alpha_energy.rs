use core::panic;
use serde_json;
use std::collections::HashMap;
use std::fmt;
use std::fs;
use std::io::BufRead;
use std::io::BufReader;

use crate::CN;

///```let binding_energy = alphas.cn[atom_type_index][in_atom_type_index][cn-1]```
///summed_to_x contains the alphas pre-summed to each coordination number
#[derive(Clone, Debug)]
pub struct Alphas {
    pub cn: [[[f64; 12]; super::NUM_ATOM_TYPES]; super::NUM_ATOM_TYPES],
    pub summed_to_x: [[[f64; 12]; super::NUM_ATOM_TYPES]; super::NUM_ATOM_TYPES],
}

#[derive(Debug, Clone)]
struct AlphaJsonError;

impl std::error::Error for AlphaJsonError {}

impl fmt::Display for AlphaJsonError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Bad alphas JSON file")
    }
}

#[derive(Debug, Clone)]
struct AlphaEnergyValueError;

impl std::error::Error for AlphaEnergyValueError {}

impl fmt::Display for AlphaEnergyValueError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Energy value is not a float")
    }
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

    pub fn new_from_json(file: String, atom_names: &mut HashMap<String, u8>) -> Alphas {
        let mut alphas_input: [[[f64; 12]; super::NUM_ATOM_TYPES]; super::NUM_ATOM_TYPES] =
            [[[0_f64; 12]; super::NUM_ATOM_TYPES]; super::NUM_ATOM_TYPES];
        let file = fs::File::open(file).expect("file should open read only");
        let json: serde_json::Value =
            serde_json::from_reader(file).expect("file should be proper JSON");

        let mut unique_atom_index = 0;
        for key_value in json.as_object().expect("root should be object").iter() {
            let diluted_in_atom_index: usize;
            let atom_index: usize;
            let key_collection: Vec<&str> = key_value.0.split('_').collect::<Vec<&str>>();
            if key_collection[0] == "dilute" {
                assert_eq!(key_collection.len(), 4, "wrong key in alphas json file");

                atom_index = *atom_names
                    .get(key_collection[1])
                    .expect("diluted energy key needs format 'diluted_Metall_on_Metall' ")
                    as usize;
                diluted_in_atom_index = *atom_names
                    .get(key_collection[3])
                    .ok_or(AlphaJsonError)
                    .expect("diluted energy key needs format 'diluted_Metall_on_Metall' ")
                    as usize;
            } else if key_collection[0] == "clean" {
                atom_index = unique_atom_index;
                diluted_in_atom_index = unique_atom_index;
                atom_names.insert(key_collection[1].to_string(), atom_index as u8);
            } else {
                panic!("bad key in alphas JSON file");
            }

            assert!(
                key_value
                    .1
                    .as_object()
                    .ok_or(AlphaJsonError)
                    .expect("bad nested object")
                    .len()
                    >= CN,
                "one of the energy lists has a wrong format"
            );
            let mut average_cn_12: f64 = 0.;
            for (cn, energy) in key_value
                .1
                .as_object()
                .ok_or(AlphaJsonError)
                .expect("bad nested object")
                .iter()
            {
                let cn_key_collection: Vec<&str> = cn.split('-').collect::<Vec<&str>>();
                if cn_key_collection.len() == 1 {
                    let cn_value = cn.parse::<usize>().expect("bad JSON");
                    alphas_input[atom_index][diluted_in_atom_index][cn_value - 1] = energy
                        .as_f64()
                        .ok_or(AlphaEnergyValueError)
                        .expect("energy value is not a flaot");
                } else {
                    if cn_key_collection[0] != "12" {
                        panic!("no 12");
                    }
                    average_cn_12 += energy
                        .as_f64()
                        .ok_or(AlphaEnergyValueError)
                        .expect("energy value is not a flaot");
                }
            }
            alphas_input[atom_index][diluted_in_atom_index][11] = average_cn_12 / 3.;

            if key_collection[0] == "clean" {
                unique_atom_index += 1;
            }
        }

        Alphas::new(alphas_input)
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

    fn _alphas_div_by_cn(
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

        // if cn_metal == 3 {
        //     energy += 0.3 ;
        // } else if cn_metal == 4 {
        //     energy += 0.4 ;
        // }

        // if cn_metal == 3 || cn_metal == 4 {
        //     return morse_potential(energy);
        // }

        nn_atom_type_count
            .iter()
            .enumerate()
            .for_each(|(metal_type, nn_atom_type_count_num)| {
                energy += self.summed_to_x[atom_type][metal_type][cn_metal - 1]
                    * *nn_atom_type_count_num as f64
                    / cn_metal as f64
            });

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

    (-e_well * (1. - E.powf(-A * (-D_MIN * 0.15)).powi(2)) + e_well)
}

///returns alphas where alhpas[x][x] is pure atom_x and alhpas[x][y] is atom_x in atom_y
pub fn read_alphas(
    alphas_file: String,
    atom_names: &mut HashMap<String, u8>,
) -> [[[f64; 12]; 2]; 2] {
    const LINE_COUNT: usize = 14;

    let path = std::path::Path::new(&alphas_file);
    let file_name = path.file_name().unwrap();
    // let mut x = alphas_file.split('.');
    let atom_names_string = file_name.to_str().unwrap().split('.').next().unwrap();
    for (i, metal) in atom_names_string.split('_').enumerate() {
        atom_names.insert(metal.to_string(), i as u8);
    }
    let pairlist = fs::File::open(alphas_file).expect("Should have been able to read the file");

    let lines = BufReader::new(pairlist);
    let mut alphas: [[[f64; 12]; 2]; 2] = [[[0.; 12]; 2]; 2];

    for (i, line) in lines.lines().enumerate() {
        let r = line.unwrap();
        let num = r.parse::<f64>().unwrap();
        if i < LINE_COUNT {
            if i >= 12 {
                continue;
            }
            alphas[0][0][i] = num;
        } else if i < LINE_COUNT * 2 {
            if i >= LINE_COUNT * 2 - 2 {
                continue;
            }
            alphas[1][1][i - LINE_COUNT * 1] = num;
        } else if i < LINE_COUNT * 3 {
            if i >= LINE_COUNT * 3 - 2 {
                continue;
            }
            alphas[1][0][i - LINE_COUNT * 2] = num;
        } else {
            if i >= LINE_COUNT * 4 - 2 {
                continue;
            }
            alphas[0][1][i - LINE_COUNT * 3] = num;
        }
    }
    alphas
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
        let mut atom_names: std::collections::HashMap<String, u8> =
            std::collections::HashMap::new();
        let alpha_file: String = "./Pt_Pd.6.bat".to_string();
        let alphas_inp = super::read_alphas(alpha_file, &mut atom_names);
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
        let mut atom_names: std::collections::HashMap<String, u8> =
            std::collections::HashMap::new();
        let alpha_file: String = "./Pt_Pd.6.bat".to_string();
        let alphas_inp = super::read_alphas(alpha_file, &mut atom_names);
        let alphas = Alphas::new(alphas_inp);

        struct TestPositin {
            cn_metal: u8,
            nn_atom_type_count: [u8; super::super::NUM_ATOM_TYPES],
            nn_nn_atom_type_count: [NnData],
        }

        let nn_nn_atom_bulk = NnData {
            cn_metal: 12,
            nn_atom_type_count_num_list: [12, 0],
            atom_type: 0,
        };

        let mut nn_nn_atom_type_count: Vec<NnData> = Vec::new();
        for _ in 0..=7 {
            nn_nn_atom_type_count.push(nn_nn_atom_bulk.clone());
        }

        let start_e = alphas.e_one_atom(11, [11, 0], nn_nn_atom_type_count.into_iter(), 0, 0, 0);
        let tst_e = alphas.e_one_atom_tst(4, [4, 0], 0, 0, 0);
        let barr = super::e_barrier(start_e - tst_e, start_e - tst_e);

        let mut atom_0 = "";
        for (atom_name, index) in atom_names.iter() {
            if index == &0 {
                atom_0 = atom_name;
            }
        }

        println!("distance {}", std::f64::consts::SQRT_2 / 2.);
        println!("morse pot {}", super::morse_potential(barr));
        println!("move in bulk of {}: {}", atom_0, barr);
    }

    #[test]
    fn test_alphas_from_json() {
        let file_path = "./alpha_alt.json".to_string();
        let mut atom_names: std::collections::HashMap<String, u8> =
            std::collections::HashMap::new();
        let mut old_atom_names: std::collections::HashMap<String, u8> =
            std::collections::HashMap::new();

        let alphas = Alphas::new_from_json(file_path, &mut atom_names);

        let old_alpha_file: String = "./alphas/Pt_Pd.6.bat".to_string();
        let alphas_inp = super::read_alphas(old_alpha_file, &mut old_atom_names);
        let old_alphas = Alphas::new(alphas_inp);

        println!("old atom_names: {:?}", old_atom_names);
        println!("new atom_names: {:?}", atom_names);

        for atom_name in atom_names.keys() {
            for in_atom_name in atom_names.keys() {
                for cn_i in 0..CN {
                    println!(
                        "{}",
                        alphas.cn[*atom_names.get(atom_name).unwrap() as usize]
                            [*atom_names.get(in_atom_name).unwrap() as usize][cn_i]
                    );
                    // assert_eq!(alphas.cn[*atom_names.get(atom_name).unwrap() as usize][*atom_names.get(in_atom_name).unwrap() as usize][cn_i]
                    //     , old_alphas.cn[*old_atom_names.get(atom_name).unwrap() as usize][*old_atom_names.get(in_atom_name).unwrap() as usize][cn_i]);
                }
            }
        }
    }
}
