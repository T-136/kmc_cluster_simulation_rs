use anyhow;
use atom_change::AtomChangeHow;
use csv::Writer;
use energy::EnergyInput;
use rand;
use rand::distributions::{Distribution, Uniform};
use rand::prelude::*;
use rand::rngs::SmallRng;
use rayon::prelude::*;
use std::collections::hash_map::Entry;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::sync::Arc;
use std::{cmp, eprint, fs, println, usize};

// mod add_remove_list;
pub mod alpha_energy;
mod atom_change;
mod buckets_linear;
pub mod energy;
mod grid_structure;
// mod main;
mod moves;
mod read_and_write;
mod setup;
mod sim;

pub use alpha_energy::Alphas;
pub use buckets_linear::ItemEnum;
pub use grid_structure::GridStructure;
pub use sim::Results;

const CN: usize = 12;
// const GCN: usize = 54;
const NUM_ATOM_TYPES: usize = 2;
const GCN: usize = 145;
const NN_PAIR_NUMBER: usize = 20;
const AMOUNT_SECTIONS: usize = 1000;

const GRID_SIZE: [u32; 3] = [30, 30, 30];

const SAVE_ENTIRE_SIM: bool = true;

// const how: add_remove::AddRemoveHow = add_remove::AddRemoveHow::RemoveAndAdd(1, 0);
// const how: add_remove::AtomChangeHow = add_remove::AtomChangeHow::Remove;
const how: atom_change::AtomChangeHow = atom_change::AtomChangeHow::Remove;

#[derive(Clone, Default)]
pub struct AtomPosition {
    occ: u8,
    cn_metal: u8,
    energy: f64,
    nn_support: u8,
    nn_atom_type_count: [u8; NUM_ATOM_TYPES],
    frozen: bool,
}

fn min_max_xyz(xyz: &Vec<[f64; 3]>, i: usize) -> (f64, f64) {
    let min = xyz
        .iter()
        .map(|p| p[i])
        .min_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    let max = xyz
        .iter()
        .map(|p| p[i])
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    (min, max)
}

fn freeze_atoms(pos: &mut [AtomPosition], freez: &[String], xyz: &Vec<[f64; 3]>) {
    let (min_x, max_x) = min_max_xyz(xyz, 0);
    let (min_y, max_y) = min_max_xyz(xyz, 1);
    let (min_z, max_z) = min_max_xyz(xyz, 2);

    xyz.iter().enumerate().for_each(|(i, &p)| {
        if (p[0] == min_x && freez.contains(&"x".to_string()))
            || (p[0] == max_x && freez.contains(&"-x".to_string()))
            || (p[1] == min_y && freez.contains(&"y".to_string()))
            || (p[1] == max_y && freez.contains(&"-y".to_string()))
            || (p[2] == min_z && freez.contains(&"z".to_string()))
            || (p[2] == max_z && freez.contains(&"-z".to_string()))
        {
            if pos[i].occ != 100 && pos[i].occ != 255 {
                println!("freezing atom: {:?}", xyz[i]);
                pos[i].frozen = true;
            }
        }
    });
}

#[derive(Clone)]
pub struct Simulation {
    atom_pos: Vec<AtomPosition>,
    atom_names: HashMap<String, u8>,
    niter: u64,
    composition: f64,
    number_all_atoms: u32,
    onlyocc: HashSet<u32, fnv::FnvBuildHasher>,
    possible_moves: buckets_linear::Buckets,
    sim_time: f64,
    total_energy: f64,
    cn_dict: [u32; CN + 1],
    cn_dict_at_supp: [u32; CN + 1],
    surface_count: Vec<i32>,
    time_per_section: Vec<f64>,
    save_folder: String,
    temperature: f64,
    cn_dict_sections: Vec<HashMap<u8, f64>>,
    energy_sections_list: Vec<f64>,
    surface_composition: Vec<f64>,
    optimization_cut_off_fraction: Vec<u64>,
    heat_map: Option<Vec<u64>>,
    snap_shot_sections: Option<Vec<Vec<u8>>>,
    heat_map_sections: Vec<Vec<u64>>,
    alphas: Arc<Alphas>,
    gridstructure: Arc<GridStructure>,
    support_e: i64,
    is_supported: bool,
    // add_or_remove: add_remove_list::AddOrRemove,
}

impl Simulation {
    pub fn new(
        atom_names: HashMap<String, u8>,
        niter: u64,
        input_file: Option<String>,
        atoms_input: Option<u32>,
        temperature: f64,
        save_folder_name: String,
        write_snap_shots: bool,
        is_heat_map: bool,
        repetition: usize,
        optimization_cut_off_fraction: Vec<u64>,
        alphas: Arc<Alphas>,
        support_indices: Option<Vec<u32>>,
        gridstructure: Arc<GridStructure>,
        coating: Option<String>,
        support_e: i64,
        freez: Option<Vec<String>>,
    ) -> Simulation {
        let nsites: u32 = GRID_SIZE[0] * GRID_SIZE[1] * GRID_SIZE[2] * 4;
        let mut atom_pos: Vec<AtomPosition> = vec![
            AtomPosition {
                occ: 255,
                ..Default::default()
            };
            nsites as usize
        ];
        let mut cn_dict: [u32; CN + 1] = [0; CN + 1];
        let mut cn_dict_at_supp: [u32; CN + 1] = [0; CN + 1];
        let is_supported = if support_e != 0 { true } else { false };
        let (onlyocc, number_all_atoms) = if input_file.is_some() {
            let xyz = read_and_write::read_sample(&input_file.unwrap());
            let onlyocc = setup::occ_onlyocc_from_xyz(
                &mut atom_pos,
                &xyz,
                &gridstructure.xsites_positions,
                &atom_names,
                coating,
                &gridstructure.nn,
            );
            // let mut nn_support = vec![0_u8; gridstructure.xsites_positions.len()];
            setup::nn_support_from_supprt(&mut atom_pos, &gridstructure.nn);
            let number_of_atoms: u32 = onlyocc.len() as u32;
            // let mut opt_nn_support: Option<Vec<u8>> = None;
            // if nn_support.iter().sum::<u8>() > 0 {
            //     opt_nn_support = Some(nn_support);
            // }
            (onlyocc, number_of_atoms)
        } else if atoms_input.is_some() {
            let onlyocc = setup::create_input_cluster(
                &mut atom_pos,
                atoms_input.unwrap(),
                &gridstructure.xsites_positions,
                &gridstructure.nn,
                // nsites,
                support_indices,
                coating,
                &atom_names,
            );
            let number_of_atom: u32 = onlyocc.len() as u32;
            (onlyocc, number_of_atom)
        } else {
            panic!("gib input atoms or input file");
        };
        // let mut cn_metal: Vec<usize> = Vec::with_capacity(nsites as usize);
        let mut surface_count = vec![0; 3];
        // let mut neighboring_atom_type_count: Vec<[u8; NUM_ATOM_TYPES]> =
        //     vec![[0; NUM_ATOM_TYPES]; nsites as usize];

        for o in 0..nsites as usize {
            let mut neighbors: u8 = 0;
            for o1 in gridstructure.nn[&(o as u32)].iter() {
                if atom_pos[*o1 as usize].occ != 255 && atom_pos[*o1 as usize].occ != 100 {
                    // cn.entry(o).and_modify(|x| *x += 1).or_insert(1);
                    let atom_type = atom_pos[*o1 as usize].occ as usize;
                    atom_pos[o].nn_atom_type_count[atom_type] += 1;
                    neighbors += 1;
                }
            }
            atom_pos[o].cn_metal = neighbors;
            if atom_pos[o].occ != 255 && atom_pos[o].occ != 100 {
                cn_dict[atom_pos[o].cn_metal as usize] += 1;
                if atom_pos[o].cn_metal != 12 {
                    surface_count[atom_pos[o].occ as usize] += 1;
                }
                if is_supported {
                    if atom_pos[o as usize].nn_support == 1 {
                        cn_dict_at_supp[atom_pos[o as usize].cn_metal as usize] += 1;
                    }
                }
            };
        }
        let mut gcn_metal: Vec<usize> = Vec::with_capacity(nsites as usize);
        for o in 0..nsites {
            let mut gcn: usize = 0;
            for o1 in gridstructure.nn[&o].iter() {
                if atom_pos[*o1 as usize].occ != 255 {
                    gcn += atom_pos[*o1 as usize].cn_metal as usize;
                }
            }
            gcn_metal.push(gcn);
        }

        let mut total_energy: f64 = 0.;
        let possible_moves: buckets_linear::Buckets = buckets_linear::Buckets::new();

        for o in onlyocc.iter() {
            let at_support = atom_pos[*o as usize].nn_support;

            let energy = alphas.e_one_atom(
                atom_pos[*o as usize].cn_metal,
                atom_pos[*o as usize].nn_atom_type_count,
                gridstructure.nn[o].iter().filter_map(|x| {
                    if atom_pos[*x as usize].occ != 255 {
                        Some(alpha_energy::NnData {
                            cn_metal: atom_pos[*x as usize].cn_metal,
                            nn_atom_type_count_num_list: atom_pos[*x as usize].nn_atom_type_count,
                            atom_type: atom_pos[*x as usize].occ as usize,
                        })
                    } else {
                        None
                    }
                }),
                atom_pos[*o as usize].occ,
                at_support,
                support_e,
            );
            atom_pos[*o as usize].energy = energy;
            total_energy += energy;
        }

        let simulation_folder_name = std::format!("{}K_{}I_{}A", temperature, niter, onlyocc.len());

        let mut sub_folder = save_folder_name + &simulation_folder_name;

        sub_folder = std::format!("{}_{}", sub_folder, repetition);

        fs::create_dir(&sub_folder).unwrap_or_else(|error| {
            println!(
                "could not create folder: {} with error: {}",
                sub_folder, error
            )
        });

        let cn_dict_sections = Vec::with_capacity(AMOUNT_SECTIONS);
        let energy_sections_list = Vec::with_capacity(AMOUNT_SECTIONS);

        let snap_shot_sections: Option<Vec<Vec<u8>>> = if write_snap_shots {
            Some(Vec::new())
        } else {
            None
        };

        let heat_map: Option<Vec<u64>> = if is_heat_map {
            Some(vec![0; nsites as usize])
        } else {
            None
        };

        let heat_map_sections: Vec<Vec<u64>> = Vec::new();

        let mut time_per_section: Vec<f64> = Vec::new();
        let mut surface_composition: Vec<f64> = Vec::new();
        let mut atom_1 = 0;
        let mut atom_2 = 0;
        for o in atom_pos.iter() {
            if o.occ == 0 {
                atom_1 += 1;
            }
            if o.occ == 1 {
                atom_2 += 1;
            }
        }
        let composition = atom_1 as f64 / (atom_1 + atom_2) as f64;
        surface_composition
            .push(surface_count[0] as f64 / (surface_count[0] + surface_count[1]) as f64);
        time_per_section.push(0.);

        if let Some(freez) = freez {
            freeze_atoms(&mut atom_pos, &freez, &gridstructure.xsites_positions);
        }

        Simulation {
            atom_pos,
            atom_names,
            niter,
            number_all_atoms,
            composition,
            onlyocc,
            possible_moves,
            sim_time: 0.,
            total_energy,
            cn_dict,
            cn_dict_at_supp,
            save_folder: sub_folder,
            temperature,
            cn_dict_sections,
            energy_sections_list,
            surface_composition,
            time_per_section,
            surface_count,
            optimization_cut_off_fraction,
            snap_shot_sections,
            heat_map,
            heat_map_sections,
            alphas,
            gridstructure,
            // neighboring_atom_type_count,
            support_e,
            is_supported,
            // add_or_remove,
        }
    }

    fn check_moves_not_start_255(&self, bucket: &buckets_linear::Bucket) {
        for item in bucket.items.iter() {
            if let buckets_linear::ItemEnum::Move(mmove) = item {
                if self.atom_pos[mmove.from as usize].occ == 255 {
                    println!("atom from {}", mmove.from);
                    println!("atomt from {}", self.atom_pos[mmove.from as usize].occ);
                    println!("atomt to {}", self.atom_pos[mmove.to as usize].occ);
                    println!("panicked item {:?}", item);
                    panic!("checked for 255");
                }
            }
        }
    }

    pub fn run(&mut self) -> Results {
        let save_every_nth: u64 = 100;
        let mut rng_choose = SmallRng::from_entropy();
        let mut bucket_pick = SmallRng::from_entropy();
        let mut coin_toss = SmallRng::from_entropy();

        let cut_off_perc = self.optimization_cut_off_fraction[0] as f64
            / self.optimization_cut_off_fraction[1] as f64;

        let mut lowest_energy_struct = sim::LowestEnergy::default();

        let mut temp_energy_section: f64 = 0.;
        let mut temp_surface_composition: f64 = 0.;
        let mut temp_cn_dict_section: [u64; CN + 1] = [0; CN + 1];

        let start = sim::Start::new(self.total_energy, &self.cn_dict);

        let mut lowest_e_onlyocc: HashSet<u32, fnv::FnvBuildHasher> =
            fnv::FnvHashSet::with_capacity_and_hasher(
                self.number_all_atoms as usize,
                Default::default(),
            );

        for (i, o) in self.atom_pos.iter().enumerate() {
            let pos_change = atom_change::AtomPosChange::new(
                i as u32,
                o.cn_metal as u8,
                o.occ,
                self.temperature,
                how.clone(),
            );
            if let Some(pos_change) = pos_change {
                self.possible_moves
                    .cond_add_item(ItemEnum::AddOrRemove(pos_change));
            }
            if o.occ != 255 && o.occ != 100 {
                for u in &self.gridstructure.nn[&(i as u32)] {
                    if self.atom_pos[*u as usize].occ == 255 {
                        // >1 so that atoms cant leave the cluster
                        // <x cant move if all neighbors are occupied
                        if self.atom_pos[*u as usize].cn_metal > 1
                            && !self.atom_pos[i as usize].frozen
                        {
                            let (prev_e, future_e) =
                                self.calc_energy_change_by_move(i as u32, *u, o.occ);
                            let e_barr = alpha_energy::e_barrier(prev_e, future_e);
                            let mmove = moves::Move::new(
                                i as u32,
                                *u,
                                e_barr,
                                future_e - prev_e,
                                self.temperature,
                            );
                            assert!(self.atom_pos[i as usize].occ != 255);
                            self.possible_moves.cond_add_item(ItemEnum::Move(mmove))
                        }
                    }
                }
            }
        }

        let section_size: u64 = self.niter / AMOUNT_SECTIONS as u64;
        println!("section_size: {}", section_size);
        println!("SAVE_TH: {}", save_every_nth);
        println!("niter: {}", self.niter);

        // print!("poss moves: {:?}", self.possible_moves.moves);

        for iiter in 0..self.niter {
            if iiter % section_size == 0 {
                println!(
                    "iteration {}; {}%",
                    iiter,
                    (iiter as f64 / self.niter as f64 * 100.)
                );
                println!(
                    "possible_moves len {} k {}",
                    self.possible_moves
                        .iter()
                        .map(|x| x.items.len())
                        .sum::<usize>(),
                    self.possible_moves.total_k
                );
                // println!("poss addremove: {:?}", self.add_or_remove.atoms);
                // println!(
                //     "buckets: {} {:?} ",
                //     self.possible_moves.buckets_list.len(),
                //     self.possible_moves
                //         .buckets_list
                //         .iter()
                //         .map(|bucket| (bucket.bucket_power, bucket.own_k))
                //         .collect::<Vec<(i32, f64)>>()
                // );
            }
            let is_recording_sections = iiter * self.optimization_cut_off_fraction[1]
                >= self.niter * self.optimization_cut_off_fraction[0];

            if !SAVE_ENTIRE_SIM
                && iiter * self.optimization_cut_off_fraction[1]
                    == self.niter * self.optimization_cut_off_fraction[0]
            {
                self.cn_dict.iter_mut().for_each(|x| {
                    *x = 0;
                });
                for o in 0..self.atom_pos.len() {
                    if self.atom_pos[o].occ != 255 && self.atom_pos[o].occ != 100 {
                        self.update_cn_dict(
                            self.atom_pos[o as usize].nn_support,
                            self.atom_pos[o as usize].cn_metal,
                            true,
                        );
                        // self.cn_dict[self.atom_pos[o].cn_metal] += 1;
                    };
                }
            };
            self.cond_snap_and_heat_map(&iiter);

            // println!(
            //     "total cn: {}, ",
            //     self.atom_pos
            //         .iter()
            //         .filter(|x| x.occ != 255)
            //         .map(|x| x.cn_metal)
            //         .sum::<usize>()
            // );

            let (item, k_tot) = self
                .possible_moves
                .choose_ramdom_move_kmc(
                    &mut rng_choose,
                    &mut bucket_pick,
                    &mut coin_toss,
                    self.temperature,
                )
                .expect("kmc pick move failed");
            match item.clone() {
                ItemEnum::Move(mmove) => {
                    // println!(
                    //     "item xyz: {:?}",
                    //     self.gridstructure.xsites_positions[mmove.from as usize]
                    // );
                    // println!("tot energy {}", self.total_energy);
                    // println!("energy {}", mmove.e_diff);
                    // println!(
                    //     "from to {:?}, {:?}",
                    //     self.atom_pos[mmove.from as usize].nn_atom_type_count,
                    //     self.atom_pos[mmove.to as usize].nn_atom_type_count
                    // );
                    self.increment_time(k_tot, &mut rng_choose);

                    self.perform_move(mmove.from, mmove.to, mmove.e_diff, is_recording_sections);
                    self.update_moves(mmove.from, mmove.to);
                    self.update_possible_moves(mmove.from, mmove.to);
                    self.update_add_remove(mmove.from, mmove.to, &how);
                    if let Some(map) = &mut self.heat_map {
                        map[mmove.to as usize] += 1;
                        map[mmove.from as usize] += 1;
                    }
                }
                ItemEnum::AddOrRemove(change) => {
                    println!("change something");
                    self.change_item(change.pos, &how, is_recording_sections);
                    self.redox_update_total_k(change.pos);
                    self.redox_update_possibel_moves(change.pos, &how, self.temperature);
                    self.redox_update_add_remove(change.pos, &how);
                }
            }

            if SAVE_ENTIRE_SIM || is_recording_sections {
                (temp_energy_section, temp_surface_composition) = self.save_sections(
                    &iiter,
                    temp_energy_section,
                    temp_surface_composition,
                    &mut temp_cn_dict_section,
                    section_size,
                    save_every_nth,
                );
            }
        }

        println!("heatmap section len: {:?}", self.heat_map_sections.len());

        read_and_write::write_occ_as_xyz(
            self.save_folder.clone(),
            lowest_e_onlyocc,
            &self.gridstructure.xsites_positions,
            &self.gridstructure.unit_cell,
            &self.atom_pos,
            &self.atom_names,
        );

        if self.heat_map.is_some() {
            let mut wtr = Writer::from_path(self.save_folder.clone() + "/heat_map.csv").unwrap();
            for heat_section in &self.heat_map_sections {
                wtr.write_record(heat_section.iter().map(|x| x.to_string()))
                    .unwrap();
            }
            wtr.flush().unwrap();
        }

        if self.snap_shot_sections.is_some() {
            let mut wtr =
                Writer::from_path(self.save_folder.clone() + "/snap_shot_sections.csv").unwrap();
            if let Some(snap_shot_sections) = &self.snap_shot_sections {
                for heat_section in snap_shot_sections {
                    wtr.write_record(heat_section.iter().map(|x| x.to_string()))
                        .unwrap();
                }
            }
            wtr.flush().unwrap();
        }

        let duration = sim::Duration {
            sec: self.sim_time,
            minutes: self.sim_time / 60.,
            hours: self.sim_time / 60. / 60.,
        };

        Results {
            start,
            lowest_energy_struct,
            composition: self.composition,
            number_all_atoms: self.number_all_atoms,
            surface_composition: self.surface_composition.clone(),
            energy_section_list: self.energy_sections_list.clone(),
            time_per_section: self.time_per_section.clone(),
            // cn_dict_sections: self.cn_dict_sections.clone(),
            duration,
        }
    }

    fn increment_time(&mut self, k_tot: f64, rng_e_number: &mut SmallRng) {
        let between = Uniform::new_inclusive(0., 1.);
        let rand_value: f64 = between.sample(rng_e_number);

        self.sim_time -= rand_value.ln() / k_tot;
        // rand_value.ln() / k_tot
        // self.sim_time += 1. / k_tot
    }

    pub fn write_exp_file(&self, exp: &Results) {
        let mut file = File::create(self.save_folder.clone() + "/exp_file.json").unwrap();
        file.write_all(serde_json::to_string_pretty(exp).unwrap().as_bytes())
            .unwrap();
    }

    fn save_sections(
        &mut self,
        iiter: &u64,
        mut temp_energy_section: f64,
        mut temp_surface_composition: f64,
        temp_cn_dict_section: &mut [u64; CN + 1],
        section_size: u64,
        SAVE_TH: u64,
    ) -> (f64, f64) {
        if (iiter + 1) % SAVE_TH == 0 {
            temp_energy_section += self.total_energy;
            temp_surface_composition += self.surface_count[0] as f64
                / (self.surface_count[0] + self.surface_count[1]) as f64;

            temp_cn_dict_section
                .iter_mut()
                .enumerate()
                .for_each(|(i, v)| *v += self.cn_dict[i] as u64);
        }

        if (iiter + 1) % section_size == 0 {
            self.energy_sections_list
                .push(temp_energy_section as f64 / (section_size / SAVE_TH) as f64);

            let mut section: HashMap<u8, f64> = HashMap::new();
            for (k, list) in temp_cn_dict_section.iter_mut().enumerate() {
                section.insert(k as u8, *list as f64 / (section_size / SAVE_TH) as f64);
                *list = 0;
            }
            assert_eq!(temp_cn_dict_section, &mut [0_u64; CN + 1]);
            self.cn_dict_sections.push(section.clone());

            self.surface_composition
                .push(temp_surface_composition as f64 / (section_size / SAVE_TH) as f64);

            self.time_per_section.push(self.sim_time);

            return (0., 0.);
        }
        (temp_energy_section, temp_surface_composition)
    }

    // fn opt_save_lowest_energy(
    //     &mut self,
    //     iiter: &u64,
    //     lowest_energy_struct: &mut sim::LowestEnergy,
    // ) -> Option<HashSet<u32, fnv::FnvBuildHasher>> {
    //     if lowest_energy_struct.energy > (self.total_energy as f64 / 1000.) {
    //         let mut empty_neighbor_cn: HashMap<u8, u32> = HashMap::new();
    //         let empty_set: HashSet<&u32> =
    //             HashSet::from_iter(self.possible_moves.iter().map(|mmove| &mmove.to));
    //         for empty in empty_set {
    //             if self.atom_pos[*empty as usize].cn_metal > 3 {
    //                 empty_neighbor_cn
    //                     .entry(self.atom_pos[*empty as usize].cn_metal as u8)
    //                     .and_modify(|x| *x += 1)
    //                     .or_insert(1);
    //             }
    //         }
    //         lowest_energy_struct.empty_cn = empty_neighbor_cn;
    //         lowest_energy_struct.energy = self.total_energy as f64 / 1000.;
    //         lowest_energy_struct.iiter = *iiter;
    //
    //         let mut cn_hash_map: HashMap<u8, u32> = HashMap::new();
    //         for (i, v) in self.cn_dict.into_iter().enumerate() {
    //             cn_hash_map.insert(i as u8, v);
    //         }
    //         lowest_energy_struct.cn_total = cn_hash_map;
    //
    //         let mut cn_hash_map_at_supp: HashMap<u8, u32> = HashMap::new();
    //         for (i, v) in self.cn_dict_at_supp.into_iter().enumerate() {
    //             cn_hash_map_at_supp.insert(i as u8, v);
    //         }
    //         lowest_energy_struct.cn_dict_at_supp = cn_hash_map_at_supp;
    //
    //         Some(self.onlyocc.clone())
    //     } else {
    //         None
    //     }
    // }

    #[inline]
    fn update_cn_dict(&mut self, nn_supp: u8, cn: u8, change_is_positiv: bool) {
        match change_is_positiv {
            true => {
                if self.is_supported {
                    if nn_supp == 1 {
                        self.cn_dict_at_supp[cn as usize] += 1;
                    }
                }
                self.cn_dict[cn as usize] += 1;
            }
            false => {
                if self.is_supported {
                    if nn_supp == 1 {
                        self.cn_dict_at_supp[cn as usize] -= 1;
                    }
                }
                self.cn_dict[cn as usize] -= 1;
            }
        }
    }

    fn filter_map_energy_diff_cn(self, x: &u32) -> Option<(usize, u8)> {
        let atom_typ_index = self.atom_pos[*x as usize].occ;
        if atom_typ_index != 0 {
            Some((
                self.atom_pos[*x as usize].cn_metal as usize,
                self.atom_pos[*x as usize].occ,
            ))
        } else {
            None
        }
    }

    fn cond_snap_and_heat_map(&mut self, iiter: &u64) {
        const NUMBER_HEAT_MAP_SECTIONS: u64 = 200;

        if let Some(snap_shot_sections) = &mut self.snap_shot_sections {
            if (iiter) % (self.niter / NUMBER_HEAT_MAP_SECTIONS) == 0 {
                let mut t_vec = vec![0; self.atom_pos.len()];
                t_vec
                    .iter_mut()
                    .enumerate()
                    .for_each(|(i, x)| *x = self.atom_pos[i].occ);
                snap_shot_sections.push(t_vec);
            }
        }

        if let Some(heat_map) = &mut self.heat_map {
            if (iiter) % (self.niter / NUMBER_HEAT_MAP_SECTIONS) == 0 {
                let mut t_vec = vec![0; heat_map.len()];
                t_vec
                    .iter_mut()
                    .enumerate()
                    .for_each(|(i, x)| *x = heat_map[i]);
                self.heat_map_sections.push(t_vec);
            }
        }
    }
}

pub fn find_simulation_with_lowest_energy(folder: String) -> anyhow::Result<()> {
    // let mut folder_with_lowest_e: PathBuf = PathBuf::new();
    let mut lowest_e: f64 = f64::INFINITY;

    for _ in 0..2 {
        let paths = fs::read_dir(&folder).unwrap();
        for path in paths {
            let ok_path = match path {
                Ok(ok_path) => ok_path,
                Err(e) => {
                    eprintln!("{:?}", e);
                    continue;
                }
            };

            if !ok_path.path().is_dir() {
                continue;
            }
            let folder = fs::read_dir(ok_path.path().as_path())?;
            for folder_entry in folder {
                let file = match folder_entry {
                    Ok(path2) => {
                        if !path2.path().is_dir() {
                            path2
                        } else {
                            println!("unexpected folder");
                            continue;
                        }
                    }
                    Err(e) => {
                        eprintln!("{:?}", e);
                        continue;
                    }
                };
                if !file.path().to_str().unwrap().ends_with(".json") {
                    continue;
                }

                let file = fs::File::open(file.path()).unwrap();
                let reader = BufReader::new(file);
                let res: Result<Results, serde_json::Error> = serde_json::from_reader(reader);
                match res {
                    Ok(res) => {
                        if res.lowest_energy_struct.energy <= lowest_e {
                            // folder_with_lowest_e = path2.path();
                            lowest_e = res.lowest_energy_struct.energy;
                        } else {
                            println!("{:?}", res.lowest_energy_struct.energy);
                            println!("{:?}", ok_path.path())
                            // fs::remove_dir_all(ok_path.path())?;
                        }
                    }
                    Err(e) => {
                        eprintln!("{:?}", format!("{:?} in folder {:?}", e, ok_path.path()));
                        // fs::remove_dir_all(ok_path.path())?;
                    }
                }
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_perform_move() {
        fn file_paths(
            grid_folder: String,
        ) -> (
            String,
            String,
            String,
            String,
            // String,
            String,
            String,
            String,
            String,
        ) {
            (
                format!("{}bulk.poscar", grid_folder),
                format!("{}nearest_neighbor", grid_folder),
                format!("{}next_nearest_neighbor", grid_folder),
                format!("{}nn_pairlist", grid_folder),
                // format!("{}nnn_pairlist", grid_folder),
                format!("{}atom_sites", grid_folder),
                format!("{}nn_pair_no_intersec", grid_folder),
                format!("{}nnn_gcn_no_intersec.json", grid_folder),
                format!("{}surrounding_moves.json", grid_folder),
            )
        }
        #[allow(unused_variables)]
        let (
            bulk_file_name,
            nn_file,
            nnn_file,
            nn_pairlist_file,
            // nnn_pairlist_file,
            atom_sites,
            nn_pair_no_int_file,
            nnn_pair_no_int_file,
            surrounding_moves_file,
        ) = file_paths("../303030-pair_kmc/".to_string());

        let gridstructure = GridStructure::new(
            nn_file,
            // nnn_file,
            nn_pair_no_int_file,
            // nnn_pair_no_int_file,
            atom_sites,
            bulk_file_name,
            surrounding_moves_file,
        );
        let gridstructure = Arc::new(gridstructure);

        let mut atom_names: HashMap<String, u8> = HashMap::new();
        // atom_names.insert("Pt".to_string(), 0);
        // atom_names.insert("Pd".to_string(), 1);
        atom_names.insert("Al".to_string(), 100);

        let alphas_arr = read_alphas("./Pt_Pd.6.bat".to_string(), &mut atom_names);
        let alphas = alpha_energy::Alphas::new(alphas_arr);

        let mut sim = Simulation::new(
            atom_names,
            1000000,
            Some("./711A.xyz".to_string()),
            None,
            900.,
            String::from("./sim/"),
            false,
            false,
            0_usize,
            vec![4, 4],
            Arc::new(alphas),
            None,
            gridstructure,
            None,
            0,
            None,
        );
        let mut rng_choose = SmallRng::from_entropy();
        let mut bucket_pick = SmallRng::from_entropy();
        let mut coin_toss = SmallRng::from_entropy();

        for (i, o) in sim.atom_pos.iter().enumerate() {
            let pos_change = atom_change::AtomPosChange::new(
                i as u32,
                o.cn_metal as u8,
                o.occ,
                sim.temperature,
                how.clone(),
            );
            if let Some(pos_change) = pos_change {
                sim.possible_moves
                    .cond_add_item(ItemEnum::AddOrRemove(pos_change));
            }
            if o.occ != 255 && o.occ != 100 {
                for u in &sim.gridstructure.nn[&(i as u32)] {
                    if sim.atom_pos[*u as usize].occ == 255 {
                        // >1 so that atoms cant leave the cluster
                        // <x cant move if all neighbors are occupied
                        if sim.atom_pos[*u as usize].cn_metal > 1
                            && !sim.atom_pos[i as usize].frozen
                        {
                            let (prev_e, future_e) =
                                sim.calc_energy_change_by_move(i as u32, *u, o.occ);
                            let e_barr = alpha_energy::e_barrier(prev_e, future_e);
                            let mmove = moves::Move::new(
                                i as u32,
                                *u,
                                e_barr,
                                future_e - prev_e,
                                sim.temperature,
                            );
                            assert!(sim.atom_pos[i as usize].occ != 255);
                            sim.possible_moves.cond_add_item(ItemEnum::Move(mmove))
                        }
                    }
                }
            }
        }

        let (item, k_tot) = sim
            .possible_moves
            .choose_ramdom_move_kmc(
                &mut rng_choose,
                &mut bucket_pick,
                &mut coin_toss,
                sim.temperature,
            )
            .expect("kmc pick move failed");
        let item = if let ItemEnum::Move(item) = item {
            item
        } else {
            panic!("p");
        };
        let (from, to, to2) = 'bar: {
            for x in sim.gridstructure.nn[&item.to] {
                if sim.gridstructure.nn[&item.from].contains(&x)
                    && sim.atom_pos[x as usize].occ == 255
                {
                    break 'bar (item.from, item.to, x);
                }
            }
            panic!("p");
            (255, 255, 255)
        };
        println!("{},{},{}", from, to, to2);
        // let (from, to) = sim.possible_moves.moves[1];

        println!(
            "nn_from: {:?} nn_to: {:?} nn_to2: {:?}",
            sim.gridstructure.nn[&from], sim.gridstructure.nn[&to], sim.gridstructure.nn[&to2]
        );
        // println!("gcn before: {:?}", sim.gcn_metal);
        assert!(
            sim.atom_pos[to2 as usize].occ == 255,
            "{}",
            sim.atom_pos[to2 as usize].occ
        );
        assert!(
            sim.atom_pos[to as usize].occ == 255,
            "{}",
            sim.atom_pos[to as usize].occ
        );

        // println!(
        //     "occ {:?}",
        //     sim.atom_pos.occ
        //         .iter()
        //         .enumerate()
        //         .filter(|(_, x)| **x == 1)
        //         .map(|(i, _)| i)
        //         .collect::<Vec<usize>>()
        // );
        let e_before_move = sim.total_energy;
        println!("total_e {}", sim.total_energy);
        println!("e_diff 1: {}\n", item.e_diff);
        sim.perform_move(from, to, item.e_diff, false);
        sim.update_moves(from, to);
        sim.update_possible_moves(from, to);
        // println!("gcn between: {:?}", sim.gcn_metal);
        // println!("gcn between: ");
        let item = sim.possible_moves.get(to, to2).unwrap();
        let item = if let ItemEnum::Move(item) = item {
            item
        } else {
            panic!("p")
        };
        println!("total_e {}", sim.total_energy);
        println!("e_diff 2: {}\n", item.e_diff);

        sim.perform_move(to, to2, item.e_diff, false);
        sim.update_moves(to, to2);
        sim.update_possible_moves(to, to2);
        // println!("gcn between2: ");
        let item = sim.possible_moves.get(to2, from).unwrap();
        let item = if let ItemEnum::Move(item) = item {
            item
        } else {
            panic!("p")
        };
        println!("total_e {}", sim.total_energy);
        println!("e_diff 2: {}\n", item.e_diff);
        sim.perform_move(to2, from, item.e_diff, false);
        // println!("gcn after: {:?}", sim.gcn_metal);
        println!("from: {}, to: {}, to2: {}\n", from, to, to2);
        println!("total_e {}", sim.total_energy);
        assert_eq!(e_before_move, sim.total_energy);
    }
}

fn read_alphas(alphas_file: String, atom_names: &mut HashMap<String, u8>) -> [[[f64; 12]; 2]; 2] {
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
        println!("{}", i);
        let r = line.unwrap();
        let num = r.parse::<f64>().unwrap();
        println!("{}", num);
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

enum FromOrTo {
    From,
    To,
}
