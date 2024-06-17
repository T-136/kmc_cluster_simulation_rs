use crate::{atom_change, buckets_linear, moves};

pub use super::alpha_energy;

// pub const ADD_ATOM: bool = true;
// pub const REMOVE_ATOM: bool = false;
// pub const EXCHANGE_ATOM: bool = false;
// pub const CHANGE_TYPES: [bool; 3] = [ADD_ATOM, REMOVE_ATOM, EXCHANGE_ATOM];
pub const ADD_ATOM_TYPE: u8 = 0;
pub const REMOVE_ATOM_TYPE: u8 = 1;

const CN_FOR_INV: u8 = 12;
// const E_RATIO: f64 = 0.20; 400K
const E_RATIO_BARR: f64 = 0.1100000;

const CN_E: [f64; 13] = [0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.];
// const CN_E_ADD: [f64; 13] = [12., 11., 10., 6., 9., 8., 6., 5., 4., 3., 2., 1., 0.];
const CN_E_ADD: [f64; 13] = [12., 11., 9., 9., 4., 9., 9., 9., 9., 9., 90., 10000., 1000.];

#[derive(Clone, Debug)]
pub struct AtomPosChange {
    pub pos: u32,
    pub k: f64,
    pub how: AtomChangeHow,
}

impl AtomPosChange {
    pub fn new(
        pos: u32,
        cn: u8,
        atom_type: u8,
        temperature: f64,
        how: AtomChangeHow,
    ) -> Option<AtomPosChange> {
        return None;
        match how {
            AtomChangeHow::Remove | AtomChangeHow::Exchange => {
                if atom_type == REMOVE_ATOM_TYPE {
                    let k = tst_rate_calculation(cn as f64 * E_RATIO_BARR, temperature);
                    return Some(AtomPosChange { pos, k, how });
                }
            }
            AtomChangeHow::Add => {
                if atom_type == 255 {
                    let k = if cn == 0 {
                        tst_rate_calculation((100) as f64 * E_RATIO_BARR, temperature)
                    } else {
                        tst_rate_calculation(CN_E_ADD[cn as usize] * E_RATIO_BARR, temperature)
                    };
                    // println!("k of change {}", k);
                    return Some(AtomPosChange { pos, k, how });
                }
            }
            AtomChangeHow::RemoveAndAdd => todo!(),
        }

        None
    }

    fn all_enum_ids(pos: u32) -> impl Iterator<Item = u64> {
        (1..5)
            .into_iter()
            .map(move |x| pos as u64 + ((x as u64) << 32))
    }
}

fn tst_rate_calculation(e_barr: f64, temperature: f64) -> f64 {
    // let e_use = if e_diff.is_negative() { 0. } else { e_diff };
    // println!("barr {}", e_barr);
    // println!("diff {}", e_diff);
    const KB_joul: f64 = 1.380649e-23;
    const h_joul: f64 = 6.62607015e-34;
    const KB_DIV_H: f64 = KB_joul / h_joul;
    const KB_eV: f64 = 8.6173324e-5;
    (KB_joul * temperature / h_joul) * ((-e_barr) / (KB_eV * temperature)).exp()
}

#[derive(Clone, Debug)]
pub enum AtomChangeHow {
    // atom_type_index
    Add = 1,
    Remove = 2,
    Exchange = 3,
    RemoveAndAdd = 4,
}

impl AtomChangeHow {
    fn get_hash_pos_enum_iter(pos: u32) -> impl std::iter::Iterator<Item = u64> {
        (0..5)
            .into_iter()
            .map(move |i| pos as u64 + ((i as u64) << 32))
    }
}

impl crate::Simulation {
    pub fn change_item(&mut self, pos: u32, how: &AtomChangeHow, is_recording_sections: bool) {
        match how {
            AtomChangeHow::Remove => {
                self.atom_pos[pos as usize].occ = 255;
                self.onlyocc.remove(&pos);

                for o in self.gridstructure.nn[&pos] {
                    if (super::SAVE_ENTIRE_SIM || is_recording_sections)
                        && self.atom_pos[o as usize].occ != 255
                        && self.atom_pos[o as usize].occ != 100
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
                        if self.atom_pos[o as usize].cn_metal == 12 {
                            self.surface_count[(self.atom_pos[o as usize].occ) as usize] += 1;
                        }
                    }
                    self.atom_pos[o as usize].cn_metal -= 1;
                    self.atom_pos[o as usize].nn_atom_type_count[REMOVE_ATOM_TYPE as usize] -= 1;
                }

                if super::SAVE_ENTIRE_SIM || is_recording_sections {
                    self.update_cn_dict(
                        self.atom_pos[pos as usize].nn_support,
                        self.atom_pos[pos as usize].cn_metal,
                        false,
                    );
                    // self.cn_dict[self.atom_pos[move_from as usize].cn_metal] -= 1;
                }
            }
            AtomChangeHow::Exchange => {
                self.atom_pos[pos as usize].occ = ADD_ATOM_TYPE;

                for o in self.gridstructure.nn[&pos] {
                    if (super::SAVE_ENTIRE_SIM || is_recording_sections)
                        && self.atom_pos[o as usize].occ != 255
                        && self.atom_pos[o as usize].occ != 100
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
                        if self.atom_pos[o as usize].cn_metal == 12 {
                            self.surface_count[(self.atom_pos[o as usize].occ) as usize] += 1;
                        }
                    }
                    // self.atom_pos[o as usize].cn_metal -= 1;
                    self.atom_pos[o as usize].nn_atom_type_count[REMOVE_ATOM_TYPE as usize] -= 1;
                    self.atom_pos[o as usize].nn_atom_type_count[ADD_ATOM_TYPE as usize] += 1;
                }
            }
            AtomChangeHow::Add => {
                self.atom_pos[pos as usize].occ = ADD_ATOM_TYPE;

                self.onlyocc.insert(pos);

                if super::SAVE_ENTIRE_SIM || is_recording_sections {
                    self.update_cn_dict(
                        self.atom_pos[pos as usize].nn_support,
                        self.atom_pos[pos as usize].cn_metal,
                        false,
                    );
                    // self.cn_dict[self.atom_pos[move_from as usize].cn_metal] -= 1;
                }

                let mut cn_changed = 0;
                for o in self.gridstructure.nn[&pos] {
                    if (super::SAVE_ENTIRE_SIM || is_recording_sections)
                        && self.atom_pos[o as usize].occ != 255
                        && self.atom_pos[o as usize].occ != 100
                        && o != pos
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
                    self.atom_pos[o as usize].nn_atom_type_count[ADD_ATOM_TYPE as usize] += 1;
                }

                if super::SAVE_ENTIRE_SIM || is_recording_sections {
                    self.update_cn_dict(
                        self.atom_pos[pos as usize].nn_support,
                        self.atom_pos[pos as usize].cn_metal,
                        true,
                    );
                    // self.cn_dict[self.atom_pos[move_to as usize].cn_metal] += 1;
                }

                // self.total_energy += self.e_change_add_move(pos, *add_atom_type);
            }
            AtomChangeHow::RemoveAndAdd => todo!(),
        }

        // println!("possible moves: {:?}", self.possible_moves.moves);

        self.total_energy += self.e_change(pos, how);
    }

    fn e_change(&self, pos: u32, how: &AtomChangeHow) -> f64 {
        let mut pos_nn_atom_type_no_tst = self.atom_pos[pos as usize].nn_atom_type_count;
        let mut energy = 0.;
        match how {
            AtomChangeHow::Remove => {
                energy -= self.alphas.e_one_atom(
                    self.atom_pos[pos as usize].cn_metal,
                    // self.atom_pos[move_from as usize].nn_atom_type_count,
                    pos_nn_atom_type_no_tst,
                    self.gridstructure.nn[&pos].iter().filter_map(|x| {
                        if self.atom_pos[*x as usize].occ != 255 {
                            Some(alpha_energy::NnData {
                                cn_metal: self.atom_pos[*x as usize].cn_metal,
                                nn_atom_type_count_num_list: self.atom_pos[*x as usize]
                                    .nn_atom_type_count,
                                atom_type: self.atom_pos[*x as usize].occ as usize,
                            })
                        } else {
                            None
                        }
                    }),
                    REMOVE_ATOM_TYPE,
                    0,
                    0,
                );
            }
            AtomChangeHow::Exchange => {
                energy -= self.alphas.e_one_atom(
                    self.atom_pos[pos as usize].cn_metal,
                    // self.atom_pos[move_from as usize].nn_atom_type_count,
                    pos_nn_atom_type_no_tst,
                    self.gridstructure.nn[&pos].iter().filter_map(|x| {
                        if self.atom_pos[*x as usize].occ != 255 {
                            Some(alpha_energy::NnData {
                                cn_metal: self.atom_pos[*x as usize].cn_metal,
                                nn_atom_type_count_num_list: self.atom_pos[*x as usize]
                                    .nn_atom_type_count,
                                atom_type: self.atom_pos[*x as usize].occ as usize,
                            })
                        } else {
                            None
                        }
                    }),
                    REMOVE_ATOM_TYPE,
                    0,
                    0,
                );

                energy += self.alphas.e_one_atom(
                    self.atom_pos[pos as usize].cn_metal,
                    pos_nn_atom_type_no_tst,
                    self.gridstructure.nn[&pos].iter().filter_map(|x| {
                        if self.atom_pos[*x as usize].occ != 255 {
                            Some(alpha_energy::NnData {
                                cn_metal: self.atom_pos[*x as usize].cn_metal,
                                nn_atom_type_count_num_list: self.atom_pos[*x as usize]
                                    .nn_atom_type_count,
                                atom_type: self.atom_pos[*x as usize].occ as usize,
                            })
                        } else {
                            None
                        }
                    }),
                    ADD_ATOM_TYPE,
                    0,
                    0,
                );
            }
            AtomChangeHow::Add => {
                let pos_nn_atom_type_ = self.atom_pos[pos as usize].nn_atom_type_count;
                energy += self.alphas.e_one_atom(
                    self.atom_pos[pos as usize].cn_metal,
                    pos_nn_atom_type_,
                    self.gridstructure.nn[&pos].iter().filter_map(|x| {
                        if self.atom_pos[*x as usize].occ != 255 {
                            Some(alpha_energy::NnData {
                                cn_metal: self.atom_pos[*x as usize].cn_metal,
                                nn_atom_type_count_num_list: self.atom_pos[*x as usize]
                                    .nn_atom_type_count,
                                atom_type: self.atom_pos[*x as usize].occ as usize,
                            })
                        } else {
                            None
                        }
                    }),
                    ADD_ATOM_TYPE,
                    0,
                    0,
                )
            }
            AtomChangeHow::RemoveAndAdd => todo!(),
        }

        energy
    }

    pub fn redox_update_possibel_moves(&mut self, pos: u32, how: &AtomChangeHow, temperature: f64) {
        match how {
            AtomChangeHow::Remove => {
                for nn_to_pos in self.gridstructure.nn[&pos] {
                    if self.atom_pos[nn_to_pos as usize].occ == 255 {
                        self.possible_moves.remove_move(pos, nn_to_pos);
                    }
                    if self.atom_pos[nn_to_pos as usize].occ != 255
                        && self.atom_pos[nn_to_pos as usize].occ != 100
                    {
                        // greater than one because of neighbor moving in this spot
                        if self.atom_pos[pos as usize].cn_metal > 1
                            && !self.atom_pos[nn_to_pos as usize].frozen
                        {
                            let (prev_e, future_e) = self.calc_energy_change_by_move(
                                nn_to_pos,
                                pos,
                                self.atom_pos[nn_to_pos as usize].occ,
                            );
                            let e_barr = alpha_energy::e_barrier(prev_e, future_e);
                            assert!(self.atom_pos[nn_to_pos as usize].occ != 100);

                            let mmove = moves::Move::new(
                                nn_to_pos,
                                pos,
                                e_barr,
                                future_e - prev_e,
                                self.temperature,
                            );
                            self.possible_moves
                                .cond_add_item(crate::ItemEnum::Move(mmove));
                        }
                    }
                }
            }
            AtomChangeHow::Exchange => todo!(),

            AtomChangeHow::Add => {
                for nn_to_atom in self.gridstructure.nn[&pos] {
                    if self.atom_pos[nn_to_atom as usize].occ != 255
                        && self.atom_pos[nn_to_atom as usize].occ != 100
                    {
                        self.possible_moves.remove_move(nn_to_atom, pos);
                    }
                    if self.atom_pos[nn_to_atom as usize].occ == 255 {
                        // greater than one because of neighbor moving in this spot
                        if self.atom_pos[nn_to_atom as usize].cn_metal > 1
                            && !self.atom_pos[pos as usize].frozen
                        {
                            let (prev_e, future_e) = self.calc_energy_change_by_move(
                                pos,
                                nn_to_atom,
                                self.atom_pos[pos as usize].occ,
                            );
                            let e_barr = alpha_energy::e_barrier(prev_e, future_e);
                            assert!(self.atom_pos[pos as usize].occ != 100);
                            let mmove = moves::Move::new(
                                pos,
                                nn_to_atom,
                                e_barr,
                                future_e - prev_e,
                                self.temperature,
                            );
                            self.possible_moves
                                .cond_add_item(crate::ItemEnum::Move(mmove));
                        }
                    }
                }
            }
            AtomChangeHow::RemoveAndAdd => todo!(),
        }
    }

    pub fn redox_update_total_k(&mut self, pos: u32) {
        for nn_to_pos in self.gridstructure.nn[&pos] {
            for nnn_to_pos in self.gridstructure.nn[&nn_to_pos] {
                let (full, empty) = if self.atom_pos[nn_to_pos as usize].occ == 255
                    && self.atom_pos[nnn_to_pos as usize].occ != 255
                    && self.atom_pos[nnn_to_pos as usize].occ != 100
                {
                    (nnn_to_pos, nn_to_pos)
                } else if self.atom_pos[nnn_to_pos as usize].occ == 255
                    && self.atom_pos[nn_to_pos as usize].occ != 255
                    && self.atom_pos[nn_to_pos as usize].occ != 100
                {
                    (nn_to_pos, nnn_to_pos)
                } else {
                    continue;
                };
                // println!(
                //     "{} {} {}",
                //     self.atom_pos[full as usize].cn_metal,
                //     self.atom_pos[empty as usize].cn_metal,
                //     self.atom_pos[full as usize].occ
                // );
                let (prev_e, future_e) =
                    self.calc_energy_change_by_move(full, empty, self.atom_pos[full as usize].occ);
                let e_barr = alpha_energy::e_barrier(prev_e, future_e);
                let mmove =
                    moves::Move::new(full, empty, e_barr, future_e - prev_e, self.temperature);
                self.possible_moves
                    .update_k_if_item_exists(crate::ItemEnum::Move(mmove));
            }
        }
    }

    pub fn redox_update_add_remove(&mut self, pos: u32, how: &AtomChangeHow) {
        match how {
            AtomChangeHow::Remove => {
                self.possible_moves.remove_add_remove(pos);
                for x in self.gridstructure.nn[&pos] {
                    let pos = atom_change::AtomPosChange::new(
                        x,
                        self.atom_pos[x as usize].cn_metal as u8,
                        self.atom_pos[x as usize].occ as u8,
                        self.temperature,
                        how.clone(),
                    );
                    if let Some(pos) = pos {
                        self.possible_moves
                            .update_k_if_item_exists(crate::ItemEnum::AddOrRemove(pos));
                    }
                }
            }
            AtomChangeHow::Exchange => {
                self.possible_moves.remove_add_remove(pos);
            }
            AtomChangeHow::Add => {
                self.possible_moves.remove_add_remove(pos);
                for neighbor in self.gridstructure.nn[&pos] {
                    let pos = atom_change::AtomPosChange::new(
                        neighbor,
                        self.atom_pos[neighbor as usize].cn_metal as u8,
                        self.atom_pos[neighbor as usize].occ,
                        self.temperature,
                        how.clone(),
                    );
                    if let Some(pos) = pos {
                        self.possible_moves
                            .cond_add_item(crate::ItemEnum::AddOrRemove(pos));
                    }

                    // maybe not necessary if prev does it
                    // self.add_or_remove.cond_update_cn(
                    //     neighbor,
                    //     self.atom_pos[neighbor as usize].cn_metal as u8,
                    //     self.temperature,
                    // );
                }
            }
            AtomChangeHow::RemoveAndAdd => todo!(),
        }
    }
}
