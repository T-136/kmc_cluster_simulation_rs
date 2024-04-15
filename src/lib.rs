use anyhow;
use chemfiles::{Atom, Frame, Trajectory, UnitCell};
use core::panic;
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

pub mod alpha_energy;
pub mod energy;
mod grid_structure;
mod listdict;
mod read_and_write;
mod setup;
mod sim;

pub use alpha_energy::AlphasTable;
pub use grid_structure::GridStructure;
pub use sim::Results;

const CN: usize = 12;
// const GCN: usize = 54;
const NUM_ATOM_TYPES: usize = 2;
const GCN: usize = 145;
const NN_PAIR_NUMBER: usize = 20;
const NN_PAIR_NO_INTERSEC_NUMBER: usize = 7;
const NN_PAIR_ONLY_INTERSEC_NUMBER: usize = 4;
const NNN_PAIR_NO_INTERSEC_NUMBER: usize = 20;
const AMOUNT_SECTIONS: usize = 1000;

const GRID_SIZE: [u32; 3] = [20, 20, 20];

const SAVE_ENTIRE_SIM: bool = true;

#[derive(Clone, Default)]
pub struct AtomPosition {
    occ: u8,
    cn_metal: usize,
    energy: f64,
    nn_support: u8,
    nn_atom_type_count: [u8; NUM_ATOM_TYPES],
    nn: [u32; CN],
}

#[derive(Clone)]
pub struct Simulation {
    atom_pos: Vec<AtomPosition>,
    atom_names: HashMap<String, u8>,
    niter: u64,
    composition: f64,
    number_all_atoms: u32,
    onlyocc: HashSet<u32, fnv::FnvBuildHasher>,
    possible_moves: listdict::ListDict,
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
    unique_levels: HashMap<BTreeMap<u8, u32>, (f64, u64)>,
    heat_map: Option<Vec<u64>>,
    snap_shot_sections: Option<Vec<Vec<u8>>>,
    heat_map_sections: Vec<Vec<u64>>,
    energy: EnergyInput,
    alphas: Arc<AlphasTable>,
    gridstructure: Arc<GridStructure>,
    support_e: i64,
    is_supported: bool,
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
        energy: EnergyInput,
        alphas: Arc<AlphasTable>,
        support_indices: Option<Vec<u32>>,
        gridstructure: Arc<GridStructure>,
        coating: Option<String>,
        support_e: i64,
    ) -> Simulation {
        let nsites: u32 = GRID_SIZE[0] * GRID_SIZE[1] * GRID_SIZE[2] * 4;
        let mut atom_pos: Vec<AtomPosition> = vec![AtomPosition::default(); nsites as usize];
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
                if atom_pos[*o1 as usize].occ != 0 && atom_pos[*o1 as usize].occ != 100 {
                    // cn.entry(o).and_modify(|x| *x += 1).or_insert(1);
                    let atom_type = atom_pos[*o1 as usize].occ as usize - 1;
                    atom_pos[o].nn_atom_type_count[atom_type] += 1;
                    neighbors += 1;
                }
            }
            atom_pos[o].cn_metal = neighbors as usize;
            if atom_pos[o].occ != 0 && atom_pos[o].occ != 100 {
                cn_dict[atom_pos[o].cn_metal] += 1;
                if atom_pos[o].cn_metal != 12 {
                    surface_count[atom_pos[o].occ as usize] += 1;
                }
                if is_supported {
                    if atom_pos[o as usize].nn_support == 1 {
                        cn_dict_at_supp[atom_pos[o as usize].cn_metal] += 1;
                    }
                }
            };
        }
        let mut gcn_metal: Vec<usize> = Vec::with_capacity(nsites as usize);
        for o in 0..nsites {
            let mut gcn: usize = 0;
            for o1 in gridstructure.nn[&o].iter() {
                if atom_pos[*o1 as usize].occ != 0 {
                    gcn += atom_pos[*o1 as usize].cn_metal;
                }
            }
            gcn_metal.push(gcn);
        }

        let mut total_energy: f64 = 0.;
        let possible_moves: listdict::ListDict = listdict::ListDict::new(GRID_SIZE);
        for o in onlyocc.iter() {
            let at_support = atom_pos[*o as usize].nn_support;
            // match energy {
            //     EnergyInput::LinearCn(_) | EnergyInput::Cn(_) => {
            //         let energy_1000: i64 = energy::energy_1000_calculation(
            //             &energy,
            //             atom_pos[*o as usize].cn_metal,
            //             atom_pos[*o as usize].occ as usize,
            //             at_support,
            //             support_e,
            //         );
            //         total_energy += energy_1000;
            //     }

            // EnergyInput::LinearGcn(_) | EnergyInput::Gcn(_) => {
            // println!("occ check{}", atom_pos[*o as usize].occ);
            let energy = alpha_energy::e_one_atom(
                (0, atom_pos[*o as usize].cn_metal),
                atom_pos[*o as usize].nn_atom_type_count,
                gridstructure.nn[o].iter().filter_map(|x| {
                    if atom_pos[*x as usize].occ != 0 {
                        Some((
                            atom_pos[*x as usize].cn_metal,
                            atom_pos[*x as usize].nn_atom_type_count,
                            atom_pos[*x as usize].occ as usize,
                        ))
                    } else {
                        None
                    }
                }),
                atom_pos[*o as usize].occ as usize,
                at_support,
                support_e,
            );
            atom_pos[*o as usize].energy = energy;
            total_energy += energy;
            // }
            // };

            // for u in &gridstructure.nn[o] {
            //     if occ[*u as usize] == 0 {
            //         // >1 so that atoms cant leave the cluster
            //         // <x cant move if all neighbors are occupied
            //         if cn_metal[*u as usize] > 1 {
            //             possible_moves.add_item(*o, *u, Some(10_i64))
            //         }
            //     }
            // }
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
        let unique_levels = HashMap::new();

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
            if o.occ == 1 {
                atom_1 += 1;
            }
            if o.occ == 2 {
                atom_2 += 1;
            }
        }
        let composition = atom_1 as f64 / (atom_1 + atom_2) as f64;
        surface_composition
            .push(surface_count[1] as f64 / (surface_count[1] + surface_count[2]) as f64);
        time_per_section.push(0.);

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
            unique_levels,
            snap_shot_sections,
            heat_map,
            heat_map_sections,
            energy,
            alphas,
            gridstructure,
            // neighboring_atom_type_count,
            support_e,
            is_supported,
        }
    }

    pub fn run(&mut self, mut amount_unique_levels: i32) -> Results {
        let save_every_nth: u64 = 100;
        let mut rng_choose = SmallRng::from_entropy();

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
        if self.niter == 0 {
            if let Some(x) = self.opt_save_lowest_energy(&0, &mut lowest_energy_struct) {
                lowest_e_onlyocc = x;
            };
        }

        for (i, o) in self.atom_pos.iter().enumerate() {
            if o.occ != 0 && o.occ != 100 {
                for u in &self.gridstructure.nn[&(i as u32)] {
                    if self.atom_pos[*u as usize].occ == 0 {
                        // >1 so that atoms cant leave the cluster
                        // <x cant move if all neighbors are occupied
                        if self.atom_pos[*u as usize].cn_metal > 1 {
                            let (prev_e, future_e) =
                                self.calc_energy_change_by_move(i as u32, *u, o.occ);
                            let e_barr = alpha_energy::e_barrier(prev_e, future_e);
                            self.possible_moves.add_item(
                                i as u32,
                                *u,
                                future_e - prev_e,
                                e_barr,
                                self.temperature,
                            )
                        }
                    }
                }
            }
        }

        self.possible_moves.calc_total_k_change(self.temperature);

        let section_size: u64 = self.niter / AMOUNT_SECTIONS as u64;
        println!("section_size: {}", section_size);
        println!("SAVE_TH: {}", save_every_nth);
        println!("niter: {}", self.niter);

        print!("poss moves: {:?}", self.possible_moves.moves);

        for iiter in 0..self.niter {
            if iiter % section_size == 0 {
                println!(
                    "iteration {}; {}%",
                    iiter,
                    (iiter as f64 / self.niter as f64 * 100.)
                );
                println!("possible_moves len {}", self.possible_moves.moves.len())
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
                    if self.atom_pos[o].occ != 0 && self.atom_pos[o].occ != 100 {
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

            let (move_from, move_to, e_diff, e_barr, k_tot, move_k) = self
                .possible_moves
                .choose_ramdom_move_kmc(&mut rng_choose, self.temperature)
                .expect("kmc pick move failed");
            // println!(
            //     "{}-{}: ktot:move_K {:.5}:{:.5} {:.2} {:.2}",
            //     move_from, move_to, k_tot, move_k, e_diff, e_barr
            // );
            self.increment_time(k_tot, &mut rng_choose);
            self.perform_move(move_from, move_to, e_diff, is_recording_sections);
            self.update_total_k(move_from, move_to);
            self.update_possible_moves(move_from, move_to);
            if let Some(map) = &mut self.heat_map {
                map[move_to as usize] += 1;
                map[move_from as usize] += 1;
            }
            if k_tot * 0.2 <= move_k || iiter % 10000 == 0 {
                self.possible_moves.calc_total_k_change(self.temperature);
            }

            if iiter * self.optimization_cut_off_fraction[1]
                >= self.niter * self.optimization_cut_off_fraction[0]
            {
                if let Some(x) = self.opt_save_lowest_energy(&iiter, &mut lowest_energy_struct) {
                    lowest_e_onlyocc = x;
                };
            }

            if SAVE_ENTIRE_SIM || is_recording_sections {
                (temp_energy_section, temp_surface_composition) = self.save_sections(
                    &iiter,
                    temp_energy_section,
                    temp_surface_composition,
                    &mut temp_cn_dict_section,
                    &mut amount_unique_levels,
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
            hour: self.sim_time / 60. / 60.,
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
            unique_levels: self.unique_levels.clone(),
            duration,
        }
    }
    fn increment_time(&mut self, k_tot: f64, rng_e_number: &mut SmallRng) {
        let between = Uniform::new_inclusive(0., 1.);
        let rand_value: f64 = between.sample(rng_e_number);
        // println!("k_tot: {}, rand_value ln: {}", k_tot, rand_value.ln());
        self.sim_time -= rand_value.ln() / k_tot
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
        amount_unique_levels: &mut i32,
        section_size: u64,
        SAVE_TH: u64,
    ) -> (f64, f64) {
        if (iiter + 1) % SAVE_TH == 0 {
            temp_energy_section += self.total_energy;
            temp_surface_composition += self.surface_count[1] as f64
                / (self.surface_count[1] + self.surface_count[2]) as f64;

            temp_cn_dict_section
                .iter_mut()
                .enumerate()
                .for_each(|(i, v)| *v += self.cn_dict[i] as u64);
        }

        if *amount_unique_levels != 0 {
            let mut cn_hash_map = HashMap::new();
            for (i, v) in self.cn_dict.into_iter().enumerate() {
                cn_hash_map.insert(i as u8, v);
            }

            if *iiter * self.optimization_cut_off_fraction[1]
                >= self.niter * self.optimization_cut_off_fraction[0]
            {
                let cn_btree: BTreeMap<_, _> = cn_hash_map.into_iter().collect();
                match self.unique_levels.entry(cn_btree) {
                    Entry::Occupied(mut entry) => {
                        let (_, x) = entry.get_mut();
                        *x += 1;
                    }
                    Entry::Vacant(entry) => {
                        if *amount_unique_levels == 1 {
                            eprint!("amount_unique_levels reached");
                        }
                        *amount_unique_levels -= 1;
                        entry.insert((self.total_energy, 1));
                    }
                }
            }
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

    fn opt_save_lowest_energy(
        &mut self,
        iiter: &u64,
        lowest_energy_struct: &mut sim::LowestEnergy,
    ) -> Option<HashSet<u32, fnv::FnvBuildHasher>> {
        if lowest_energy_struct.energy > (self.total_energy as f64 / 1000.) {
            let mut empty_neighbor_cn: HashMap<u8, u32> = HashMap::new();
            let empty_set: HashSet<&u32> =
                HashSet::from_iter(self.possible_moves.iter().map(|mmove| &mmove.to));
            for empty in empty_set {
                if self.atom_pos[*empty as usize].cn_metal > 3 {
                    empty_neighbor_cn
                        .entry(self.atom_pos[*empty as usize].cn_metal as u8)
                        .and_modify(|x| *x += 1)
                        .or_insert(1);
                }
            }
            lowest_energy_struct.empty_cn = empty_neighbor_cn;
            lowest_energy_struct.energy = self.total_energy as f64 / 1000.;
            lowest_energy_struct.iiter = *iiter;

            let mut cn_hash_map: HashMap<u8, u32> = HashMap::new();
            for (i, v) in self.cn_dict.into_iter().enumerate() {
                cn_hash_map.insert(i as u8, v);
            }
            lowest_energy_struct.cn_total = cn_hash_map;

            let mut cn_hash_map_at_supp: HashMap<u8, u32> = HashMap::new();
            for (i, v) in self.cn_dict_at_supp.into_iter().enumerate() {
                cn_hash_map_at_supp.insert(i as u8, v);
            }
            lowest_energy_struct.cn_dict_at_supp = cn_hash_map_at_supp;

            Some(self.onlyocc.clone())
        } else {
            None
        }
    }

    #[inline]
    fn update_cn_dict(&mut self, nn_supp: u8, cn: usize, change_is_positiv: bool) {
        match change_is_positiv {
            true => {
                if self.is_supported {
                    if nn_supp == 1 {
                        self.cn_dict_at_supp[cn] += 1;
                    }
                }
                self.cn_dict[cn] += 1;
            }
            false => {
                if self.is_supported {
                    if nn_supp == 1 {
                        self.cn_dict_at_supp[cn] -= 1;
                    }
                }
                self.cn_dict[cn] -= 1;
            }
        }
    }

    fn perform_move(
        &mut self,
        move_from: u32,
        move_to: u32,
        e_diff: f64,
        is_recording_sections: bool,
    ) {
        let (from_change, to_change, inter) =
            no_int_nn_from_move(move_from, move_to, &self.gridstructure.nn_pair_no_intersec);

        self.atom_pos[move_to as usize].occ = self.atom_pos[move_from as usize].occ; // covers different alloys also
        self.atom_pos[move_from as usize].occ = 0;

        self.onlyocc.remove(&move_from);
        self.onlyocc.insert(move_to);

        if SAVE_ENTIRE_SIM || is_recording_sections {
            self.update_cn_dict(
                self.atom_pos[move_from as usize].nn_support,
                self.atom_pos[move_from as usize].cn_metal,
                false,
            );
            // self.cn_dict[self.atom_pos[move_from as usize].cn_metal] -= 1;
        }
        // println!("possible moves: {:?}", self.possible_moves.moves);
        // for o in self.gridstructure.nn[&move_from] {
        let atom_type = self.atom_pos[move_to as usize].occ as usize - 1;
        for o in from_change {
            if (SAVE_ENTIRE_SIM || is_recording_sections)
                && self.atom_pos[o as usize].occ != 0
                && self.atom_pos[o as usize].occ != 100
                && o != move_to
            {
                self.update_cn_dict(
                    self.atom_pos[o as usize].nn_support,
                    self.atom_pos[o as usize].cn_metal,
                    false,
                );
                self.update_cn_dict(
                    self.atom_pos[o as usize].nn_support,
                    self.atom_pos[o as usize].cn_metal - 1,
                    true,
                );
                // self.cn_dict[self.atom_pos[o as usize].cn_metal] -= 1;
                // self.cn_dict[self.atom_pos[o as usize].cn_metal - 1] += 1;

                // self.surface_count[self.atom_pos.occ[o as usize] as usize] -=
                //     one_if_12(self.atom_pos.cn_metal[o as usize]);
                if self.atom_pos[o as usize].cn_metal == 12 {
                    self.surface_count[(self.atom_pos[o as usize].occ) as usize] += 1;
                }
            }
            self.atom_pos[o as usize].cn_metal -= 1;
            self.atom_pos[o as usize].nn_atom_type_count[atom_type] -= 1;
        }

        for o in to_change {
            // for o in self.gridstructure.nn[&move_to] {
            if (SAVE_ENTIRE_SIM || is_recording_sections)
                && self.atom_pos[o as usize].occ != 0
                && self.atom_pos[o as usize].occ != 100
                && o != move_from
            {
                self.update_cn_dict(
                    self.atom_pos[o as usize].nn_support,
                    self.atom_pos[o as usize].cn_metal,
                    false,
                );
                self.update_cn_dict(
                    self.atom_pos[o as usize].nn_support,
                    self.atom_pos[o as usize].cn_metal + 1,
                    true,
                );
                // self.cn_dict[self.atom_pos[o as usize].cn_metal] -= 1;
                // self.cn_dict[self.atom_pos[o as usize].cn_metal + 1] += 1;

                if self.atom_pos[o as usize].cn_metal == 11 {
                    self.surface_count[(self.atom_pos[o as usize].occ) as usize] -= 1;
                }
            }
            self.atom_pos[o as usize].cn_metal += 1;
            self.atom_pos[o as usize].nn_atom_type_count[atom_type] += 1;
        }

        self.atom_pos[move_to as usize].cn_metal -= 1;
        self.atom_pos[move_from as usize].cn_metal += 1;
        self.atom_pos[move_to as usize].nn_atom_type_count[atom_type] -= 1;
        self.atom_pos[move_from as usize].nn_atom_type_count[atom_type] += 1;

        // if SAVE_ENTIRE_SIM || is_recording_sections {
        //     self.cn_dict[self.atom_pos.cn_metal[move_to as usize]] += 1;
        // }

        // match self.energy {
        //     EnergyInput::LinearCn(_) | EnergyInput::Cn(_) => {}
        //     EnergyInput::LinearGcn(_) | EnergyInput::Gcn(_) => {
        //         let (from_change, to_change, intersect, is_reverse) = no_int_nnn_from_move(
        //             move_from,
        //             move_to,
        //             &self.gridstructure.nnn_pair_no_intersec,
        //         );
        //         let (from_change_nn, to_change_nn) = no_int_nn_from_move(
        //             move_from,
        //             move_to,
        //             &self.gridstructure.nn_pair_no_intersec,
        //         );
        //         for o in from_change_nn {
        //             if o == move_to {
        //                 continue;
        //             }
        //             if self.atom_pos[o as usize].occ == 1 {
        //                 self.atom_pos[move_from as usize].gcn_metal -= 1;
        //             }
        //         }
        //         self.atom_pos[move_from as usize].gcn_metal +=
        //             self.atom_pos[move_to as usize].cn_metal;
        //         for o in to_change_nn {
        //             if o == move_from {
        //                 continue;
        //             }
        //             if self.atom_pos[o as usize].occ == 1 {
        //                 self.atom_pos[move_to as usize].gcn_metal += 1;
        //             }
        //         }
        //         self.atom_pos[move_to as usize].gcn_metal -=
        //             self.atom_pos[move_from as usize].cn_metal - 1;
        //         for atom_and_neighbors in to_change {
        //             for n in atom_and_neighbors.iter().skip(1) {
        //                 if n == &move_to {
        //                     self.atom_pos[atom_and_neighbors[0] as usize].gcn_metal +=
        //                         self.atom_pos[*n as usize].cn_metal;
        //                     continue;
        //                 }
        //                 #[cfg(debug_assertions)]
        //                 if n == &move_from {
        //                     panic!(
        //                         "found move from {:?}, move_from {}",
        //                         atom_and_neighbors, move_from
        //                     );
        //                 }
        //                 if self.atom_pos[*n as usize].occ == 1 {
        //                     self.atom_pos[atom_and_neighbors[0] as usize].gcn_metal += 1;
        //                 }
        //             }
        //         }
        //         for atom_and_neighbors in from_change {
        //             // println!("gcn bef{:?}", self.atom_pos.gcn_metal[*atom as usize]);
        //             for n in atom_and_neighbors.iter().skip(1) {
        //                 #[cfg(debug_assertions)]
        //                 if n == &move_to {
        //                     panic!("found move to");
        //                 }
        //                 if n == &move_from {
        //                     // println!("cn move from {:?}", self.atom_pos.cn_metal[*n as usize] - 1);
        //                     self.atom_pos[atom_and_neighbors[0] as usize].gcn_metal -=
        //                         self.atom_pos[*n as usize].cn_metal - 1;
        //                     continue;
        //                 }
        //                 if self.atom_pos[*n as usize].occ == 1 {
        //                     // println!("-1",);
        //                     self.atom_pos[atom_and_neighbors[0] as usize].gcn_metal -= 1;
        //                 }
        //             }
        //         }
        //         for (atom, first_neighbors, second_neighbors, to_from_atoms) in intersect {
        //             if !is_reverse {
        //                 for n in first_neighbors {
        //                     if self.atom_pos[*n as usize].occ == 1 {
        //                         self.atom_pos[*atom as usize].gcn_metal -= 1;
        //                     }
        //                 }
        //                 for n in second_neighbors {
        //                     if self.atom_pos[*n as usize].occ == 1 {
        //                         self.atom_pos[*atom as usize].gcn_metal += 1;
        //                     }
        //                 }
        //             } else if is_reverse {
        //                 for n in second_neighbors {
        //                     if self.atom_pos[*n as usize].occ == 1 {
        //                         self.atom_pos[*atom as usize].gcn_metal -= 1;
        //                     }
        //                 }
        //                 for n in first_neighbors {
        //                     if self.atom_pos[*n as usize].occ == 1 {
        //                         self.atom_pos[*atom as usize].gcn_metal += 1;
        //                     }
        //                 }
        //             }
        //             for n in to_from_atoms {
        //                 if n == &move_to {
        //                     self.atom_pos[*atom as usize].gcn_metal +=
        //                         self.atom_pos[*n as usize].cn_metal;
        //                     continue;
        //                 }
        //                 if n == &move_from {
        //                     self.atom_pos[*atom as usize].gcn_metal -=
        //                         self.atom_pos[*n as usize].cn_metal - 1;
        //                     continue;
        //                 }
        //                 panic!("neither start nor end found");
        //             }
        //         }
        //     }
        // }
        if SAVE_ENTIRE_SIM || is_recording_sections {
            self.update_cn_dict(
                self.atom_pos[move_to as usize].nn_support,
                self.atom_pos[move_to as usize].cn_metal,
                true,
            );
            // self.cn_dict[self.atom_pos[move_to as usize].cn_metal] += 1;
        }

        self.total_energy += e_diff;
    }

    fn filter_map_energy_diff_cn(self, x: &u32) -> Option<(usize, u8)> {
        let atom_typ_index = self.atom_pos[*x as usize].occ;
        if atom_typ_index != 0 {
            Some((
                self.atom_pos[*x as usize].cn_metal,
                self.atom_pos[*x as usize].occ,
            ))
        } else {
            None
        }
    }

    fn calc_energy_change_by_move(&self, move_from: u32, move_to: u32, atom_typ: u8) -> (f64, f64) {
        let atom_typ_index = atom_typ - 1;
        // println!("{}", atom_typ_index);
        // println!("{:?}", self.atom_pos.neigboring_atom_type_count);
        //
        let (from_change_nn, to_change_nn, inter_nn) =
            no_int_nn_from_move(move_from, move_to, &self.gridstructure.nn_pair_no_intersec);

        let mut nn_atom_type_count_tst = [0; NUM_ATOM_TYPES];
        let mut cn_tst = 0;
        let mut from_nn_atom_type_no_tst = self.atom_pos[move_from as usize].nn_atom_type_count;
        let mut to_nn_atom_type_count = self.atom_pos[move_to as usize].nn_atom_type_count;
        to_nn_atom_type_count[atom_typ_index as usize] -= 1;

        // inter_nn.iter().for_each(|x| {
        for x in inter_nn {
            if self.atom_pos[x as usize].occ != 0 && x != move_from {
                nn_atom_type_count_tst[self.atom_pos[x as usize].occ as usize - 1] += 1;
                // from_nn_atom_type_no_tst[self.atom_pos[x as usize].occ as usize - 1] -= 1;
                // to_nn_atom_type_count[self.atom_pos[x as usize].occ as usize - 1] -= 1;
                cn_tst += 1;
            }
        }

        let tst_e = self.alphas.e_one_atom_tst(
            cn_tst,
            nn_atom_type_count_tst,
            // inter_nn.iter().filter_map(|x| {
            //     if self.atom_pos[*x as usize].occ != 0 {
            //         let mut cn_metal_future = self.atom_pos[*x as usize].cn_metal;
            //         let mut nn_atom_type_count_future =
            //             self.atom_pos[*x as usize].nn_atom_type_count;
            //         nn_atom_type_count_future[atom_typ_index as usize] += 1;
            //         cn_metal_future += 1;
            //
            //         Some((
            //             cn_metal_future,
            //             nn_atom_type_count_future,
            //             self.atom_pos[*x as usize].occ as usize,
            //         ))
            //     } else {
            //         None
            //     }
            // }),
            self.atom_pos[move_from as usize].occ as usize,
            0,
            0,
        );

        let prev_e = self.alphas.e_one_atom(
            self.atom_pos[move_from as usize].cn_metal,
            // self.atom_pos[move_from as usize].nn_atom_type_count,
            from_nn_atom_type_no_tst,
            from_change_nn.iter().filter_map(|x| {
                if self.atom_pos[*x as usize].occ != 0 {
                    Some(
                        alpha_energy::NnData {
                            cn_metal: self.atom_pos[*x as usize].cn_metal,
                            nn_atom_type_count_num_list: self.atom_pos[*x as usize]
                                .nn_atom_type_count,
                            atom_type: self.atom_pos[*x as usize].occ as usize,
                        }, //     (
                           //     self.atom_pos[*x as usize].cn_metal,
                           //     self.atom_pos[*x as usize].nn_atom_type_count,
                           //     self.atom_pos[*x as usize].occ as usize,
                           // )
                    )
                } else {
                    None
                }
            }),
            self.atom_pos[move_from as usize].occ as usize,
            0,
            0,
        );

        let future_e = self.alphas.e_one_atom(
            self.atom_pos[move_to as usize].cn_metal - 1,
            to_nn_atom_type_count,
            to_change_nn.iter().filter_map(|x| {
                if self.atom_pos[*x as usize].occ != 0 && *x != move_from {
                    let mut cn_metal_future = self.atom_pos[*x as usize].cn_metal;
                    let mut nn_atom_type_count_future =
                        self.atom_pos[*x as usize].nn_atom_type_count;
                    nn_atom_type_count_future[atom_typ_index as usize] += 1;
                    cn_metal_future += 1;

                    Some(
                        alpha_energy::NnData {
                            cn_metal: cn_metal_future,
                            nn_atom_type_count_num_list: nn_atom_type_count_future,
                            atom_type: self.atom_pos[*x as usize].occ as usize,
                        },
                        //     (
                        //     cn_metal_future,
                        //     nn_atom_type_count_future,
                        //     self.atom_pos[*x as usize].occ as usize,
                        // )
                    )
                } else {
                    None
                }
            }),
            self.atom_pos[move_from as usize].occ as usize,
            0,
            0,
        );

        // println!("prev_e: {}, future_e: {} tst_e {}", prev_e, future_e, tst_e);
        (prev_e - tst_e, future_e - tst_e)
    }

    fn update_total_k(&mut self, move_from: u32, move_to: u32) {
        // match self.energy {
        //     EnergyInput::LinearCn(_) | EnergyInput::Cn(_) => {
        // let i_slice = &self.gridstructure.nn_pair_no_intersec
        //     [&(move_from as u64 + ((move_to as u64) << 32))];
        for pos in self.gridstructure.surrounding_moves[&(std::cmp::min(move_from, move_to) as u64
            + ((std::cmp::max(move_to, move_from) as u64) << 32))]
            .iter()
        {
            if self.atom_pos[pos.0 as usize].occ != 0
                && self.atom_pos[pos.0 as usize].occ != 100
                && self.atom_pos[pos.1 as usize].occ == 0
            {
                let (prev_e, future_e) = self.calc_energy_change_by_move(
                    pos.0,
                    pos.1,
                    self.atom_pos[pos.0 as usize].occ,
                );
                let e_barr = alpha_energy::e_barrier(prev_e, future_e);
                self.possible_moves.update_k_if_move_exists(
                    pos.0,
                    pos.1,
                    future_e - prev_e,
                    e_barr,
                    self.temperature,
                );
            }
            if self.atom_pos[pos.1 as usize].occ != 100
                && self.atom_pos[pos.1 as usize].occ != 0
                && self.atom_pos[pos.0 as usize].occ == 0
            {
                let (prev_e, future_e) = self.calc_energy_change_by_move(
                    pos.1,
                    pos.0,
                    self.atom_pos[pos.1 as usize].occ,
                );
                let e_barr = alpha_energy::e_barrier(prev_e, future_e);
                self.possible_moves.update_k_if_move_exists(
                    pos.1,
                    pos.0,
                    future_e - prev_e,
                    e_barr,
                    self.temperature,
                );
            }
        }
        //     }
        //     EnergyInput::LinearGcn(_) => todo!(),
        //     EnergyInput::Gcn(_) => todo!(),
        // }
    }

    fn update_possible_moves(&mut self, move_from: u32, move_to: u32) {
        self.possible_moves.remove_item(move_from, move_to);
        for nn_to_from in self.gridstructure.nn[&move_from] {
            if self.atom_pos[nn_to_from as usize].occ == 0 {
                self.possible_moves.remove_item(move_from, nn_to_from);
            }
            if self.atom_pos[nn_to_from as usize].occ != 0
                && self.atom_pos[nn_to_from as usize].occ != 100
            {
                // greater than one because of neighbor moving in this spot
                if self.atom_pos[move_from as usize].cn_metal > 1 {
                    let (prev_e, future_e) = self.calc_energy_change_by_move(
                        nn_to_from,
                        move_from,
                        self.atom_pos[nn_to_from as usize].occ,
                    );
                    let e_barr = alpha_energy::e_barrier(prev_e, future_e);
                    assert!(self.atom_pos[nn_to_from as usize].occ != 100);
                    self.possible_moves.add_item(
                        nn_to_from,
                        move_from,
                        future_e - prev_e,
                        e_barr,
                        self.temperature,
                    );
                }
            }
        }

        for nn_to_to in self.gridstructure.nn[&move_to] {
            if self.atom_pos[nn_to_to as usize].occ != 0
                && self.atom_pos[nn_to_to as usize].occ != 100
            {
                self.possible_moves.remove_item(nn_to_to, move_to);
            }
            if self.atom_pos[nn_to_to as usize].occ == 0 {
                // greater than one because of neighbor moving in this spot
                if self.atom_pos[nn_to_to as usize].cn_metal > 1 {
                    let (prev_e, future_e) = self.calc_energy_change_by_move(
                        move_to,
                        nn_to_to,
                        self.atom_pos[move_to as usize].occ,
                    );
                    let e_barr = alpha_energy::e_barrier(prev_e, future_e);
                    assert!(self.atom_pos[move_to as usize].occ != 100);
                    self.possible_moves.add_item(
                        move_to,
                        nn_to_to,
                        future_e - prev_e,
                        e_barr,
                        self.temperature,
                    );
                }
            }
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
fn no_int_nn_from_move(
    move_from: u32,
    move_to: u32,
    nn_pair_no_intersec: &HashMap<
        u64,
        (
            [u32; NN_PAIR_NO_INTERSEC_NUMBER],
            [u32; NN_PAIR_NO_INTERSEC_NUMBER],
            [u32; NN_PAIR_ONLY_INTERSEC_NUMBER],
        ),
        fnv::FnvBuildHasher,
    >,
) -> ([u32; 7], [u32; 7], [u32; NN_PAIR_ONLY_INTERSEC_NUMBER]) {
    let no_int = nn_pair_no_intersec[&(std::cmp::min(move_from, move_to) as u64
        + ((std::cmp::max(move_to, move_from) as u64) << 32))];
    if move_to > move_from {
        (no_int.0, no_int.1, no_int.2)
    } else {
        (no_int.1, no_int.0, no_int.2)
    }
}

fn no_int_nnn_from_move(
    move_from: u32,
    move_to: u32,
    nnn_pair_no_intersec: &HashMap<
        u64,
        (
            Vec<Vec<u32>>,
            Vec<Vec<u32>>,
            Vec<(u32, Vec<u32>, Vec<u32>, Vec<u32>)>,
        ),
        fnv::FnvBuildHasher,
    >,
) -> (
    &Vec<Vec<u32>>,
    &Vec<Vec<u32>>,
    &Vec<(u32, Vec<u32>, Vec<u32>, Vec<u32>)>,
    bool,
) {
    let (min, max, inter) = &nnn_pair_no_intersec[&(std::cmp::min(move_from, move_to) as u64
        + ((std::cmp::max(move_to, move_from) as u64) << 32))];
    if move_to > move_from {
        (&min, &max, &inter, false)
    } else {
        (&max, &min, &inter, true)
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

// #[cfg(test)]
// mod tests {
//     // Note this useful idiom: importing names from outer (for mod tests) scope.
//     use super::*;
//
//     #[test]
//     fn test_perform_move() {
//         fn file_paths(
//             grid_folder: String,
//         ) -> (String, String, String, String, String, String, String) {
//             (
//                 format!("{}/nearest_neighbor", grid_folder),
//                 format!("{}/next_nearest_neighbor", grid_folder),
//                 format!("{}/nn_pairlist", grid_folder),
//                 format!("{}/nnn_pairlist", grid_folder),
//                 format!("{}/atom_sites", grid_folder),
//                 format!("{}/nn_pair_no_intersec", grid_folder),
//                 format!("{}/nnn_gcn_no_intersec.json", grid_folder),
//                 format!("{}surrounding_moves.json", grid_folder),
//             )
//         }
//         let (
//             pairlist_file,
//             n_pairlist_file,
//             nn_pairlist_file,
//             nnn_pairlist_file,
//             atom_sites,
//             nn_pair_no_int_file,
//             nnn_pair_no_int_file,
//             surrounding_moves_file,
//         ) = file_paths("../999-pair_inline".to_string());
//
//         let energy = EnergyInput::Gcn([
//             4752, 4719, 4686, 4653, 4620, 4587, 4554, 4521, 4488, 4455, 4422, 4389, 4356, 4323,
//             4290, 4257, 4224, 4191, 4158, 4125, 4092, 4059, 4026, 3993, 3960, 3927, 3894, 3861,
//             3828, 3795, 3762, 3729, 3696, 3663, 3630, 3597, 3564, 3531, 3498, 3465, 3432, 3399,
//             3366, 3333, 3300, 3267, 3234, 3201, 3168, 3135, 3102, 3069, 3036, 3003, 2970, 2937,
//             2904, 2871, 2838, 2805, 2772, 2739, 2706, 2673, 2640, 2607, 2574, 2541, 2508, 2475,
//             2442, 2409, 2376, 2343, 2310, 2277, 2244, 2211, 2178, 2145, 2112, 2079, 2046, 2013,
//             1980, 1947, 1914, 1881, 1848, 1815, 1782, 1749, 1716, 1683, 1650, 1617, 1584, 1551,
//             1518, 1485, 1452, 1419, 1386, 1353, 1320, 1287, 1254, 1221, 1188, 1155, 1122, 1089,
//             1056, 1023, 990, 957, 924, 891, 858, 825, 792, 759, 726, 693, 660, 627, 594, 561, 528,
//             495, 462, 429, 396, 363, 330, 297, 264, 231, 198, 165, 132, 99, 66, 33, 0,
//         ]);
//
//         let gridstructure = GridStructure::new(
//             pairlist_file,
//             n_pairlist_file,
//             nn_pair_no_int_file,
//             nnn_pair_no_int_file,
//             atom_sites,
//             String::from("../input_cluster/bulk.poscar"),
//             surrounding_moves_file,
//         );
//
//         let mut sim = Simulation::new(
//             1000000000,
//             None,
//             Some(20),
//             300.,
//             String::from("./sim/"),
//             false,
//             false,
//             0_usize,
//             vec![3, 4],
//             energy,
//             None,
//             Arc::new(gridstructure),
//         );
//         let (from, to, to2) = 'bar: {
//             for (from, to, _) in &sim.possible_moves.moves {
//                 for x in sim.gridstructure.nn[to] {
//                     if sim.gridstructure.nn[from].contains(&x) && sim.occ[x as usize] == 0 {
//                         break 'bar (*from, *to, x);
//                     }
//                 }
//             }
//             (0, 0, 0)
//         };
//         println!("{},{},{}", from, to, to2);
//         // let (from, to) = sim.possible_moves.moves[1];
//         let mut effected_gcn = Vec::new();
//         for o in sim.gridstructure.nn[&from] {
//             for x in sim.gridstructure.nn[&o] {
//                 if !effected_gcn.contains(&x) && sim.cn_metal[x as usize] != 0 {
//                     effected_gcn.push(x)
//                 }
//             }
//         }
//         for o in sim.gridstructure.nn[&to] {
//             for x in sim.gridstructure.nn[&o] {
//                 if !effected_gcn.contains(&x) && sim.cn_metal[x as usize] != 0 {
//                     effected_gcn.push(x)
//                 }
//             }
//         }
//         for o in sim.gridstructure.nn[&to2] {
//             for x in sim.gridstructure.nn[&o] {
//                 if !effected_gcn.contains(&x) && sim.cn_metal[x as usize] != 0 {
//                     effected_gcn.push(x)
//                 }
//             }
//         }
//         let old_gcn = sim
//             .gcn_metal
//             .clone()
//             .into_iter()
//             .enumerate()
//             .filter(|(i, _)| effected_gcn.contains(&(*i as u32)))
//             .map(|(i, x)| (i, x))
//             .collect::<Vec<(usize, usize)>>();
//
//         println!(
//             "nn_from: {:?} nn_to: {:?} nn_to2: {:?}",
//             sim.gridstructure.nn[&from], sim.gridstructure.nn[&to], sim.gridstructure.nn[&to2]
//         );
//         // println!("gcn before: {:?}", sim.gcn_metal);
//         assert!(sim.occ[to2 as usize] == 0, "{}", sim.occ[to2 as usize]);
//         assert!(sim.occ[to as usize] == 0, "{}", sim.occ[to as usize]);
//
//         println!(
//             "occ {:?}",
//             sim.occ
//                 .iter()
//                 .enumerate()
//                 .filter(|(_, x)| **x == 1)
//                 .map(|(i, _)| i)
//                 .collect::<Vec<usize>>()
//         );
//         sim.perform_move(from, to, 0, false);
//         // println!("gcn between: {:?}", sim.gcn_metal);
//         // println!("gcn between: ");
//         sim.perform_move(to, to2, 0, false);
//         // println!("gcn between2: ");
//         println!(
//             "{:?}",
//             sim.gcn_metal
//                 .clone()
//                 .into_iter()
//                 .enumerate()
//                 .filter(|(i, _)| effected_gcn.contains(&(*i as u32)))
//                 .map(|(i, x)| (i, x))
//                 .collect::<Vec<(usize, usize)>>()
//         );
//         sim.perform_move(to2, from, 0, false);
//         // println!("gcn after: {:?}", sim.gcn_metal);
//         let (from_change, to_change, intercet, is_reverse) =
//             no_int_nnn_from_move(from, to, &sim.gridstructure.nnn_pair_no_intersec);
//         let (from_change2, to_change2, intercet2, _) =
//             no_int_nnn_from_move(to, to2, &sim.gridstructure.nnn_pair_no_intersec);
//         let (from_change3, to_change3, intercet3, _) =
//             no_int_nnn_from_move(to2, from, &sim.gridstructure.nnn_pair_no_intersec);
//         println!("from: {}, to: {}, to2: {}\n", from, to, to2);
//         println!(
//             "from: {:?},\n to: {:?},\n inter: {:?} \n",
//             from_change, to_change, intercet
//         );
//         println!(
//             "from: {:?},\n to: {:?},\n inter: {:?} \n",
//             from_change2, to_change2, intercet2
//         );
//         println!(
//             "from: {:?},\n to: {:?},\n inter: {:?} \n",
//             from_change3, to_change3, intercet3
//         );
//         let new_gcn = sim
//             .gcn_metal
//             .clone()
//             .into_iter()
//             .enumerate()
//             .filter(|(i, _)| effected_gcn.contains(&(*i as u32)))
//             .map(|(i, x)| (i, x))
//             .collect::<Vec<(usize, usize)>>();
//         assert_eq!(old_gcn, new_gcn);
//     }
// }
enum FromOrTo {
    From,
    To,
}
