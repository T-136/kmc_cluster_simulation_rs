pub use super::alpha_energy;

pub enum AddRemoveHow {
    // atom_type_index
    // Add(u8),
    Remove(u8),
    RemoveAndAdd(u8, u8),
}

impl crate::Simulation {
    pub fn change_item(&mut self, pos: u32, how: &AddRemoveHow, is_recording_sections: bool) {
        match how {
            AddRemoveHow::Remove(remove_atom_index) => {
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
                    self.atom_pos[o as usize].nn_atom_type_count[*remove_atom_index as usize] -= 1;
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
            AddRemoveHow::RemoveAndAdd(remove_atom_index, add_atom_index) => {
                self.atom_pos[pos as usize].occ = *add_atom_index + 1;

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
                    self.atom_pos[o as usize].nn_atom_type_count[*remove_atom_index as usize] -= 1;
                    self.atom_pos[o as usize].nn_atom_type_count[*add_atom_index as usize] += 1;
                }
            }
        }

        // println!("possible moves: {:?}", self.possible_moves.moves);

        self.total_energy += self.e_change_add_remove(pos, how);
    }

    fn e_change_add_remove(&self, pos: u32, how: &AddRemoveHow) -> f64 {
        let mut pos_nn_atom_type_no_tst = self.atom_pos[pos as usize].nn_atom_type_count;
        let mut energy = 0.;
        match how {
            AddRemoveHow::Remove(remove_atom_index) => {
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
                    *remove_atom_index as usize + 1,
                    0,
                    0,
                );
            }
            AddRemoveHow::RemoveAndAdd(add_atom_index, remove_atom_index) => {
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
                    *remove_atom_index as usize + 1,
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
                    *add_atom_index as usize + 1,
                    0,
                    0,
                );
            }
        }

        energy
    }

    pub fn redox_update_possibel_moves(&mut self, pos: u32, how: &AddRemoveHow, temperature: f64) {
        match how {
            AddRemoveHow::Remove(remove_atom_index) => {
                for nn_to_pos in self.gridstructure.nn[&pos] {
                    if self.atom_pos[nn_to_pos as usize].occ == 255 {
                        self.possible_moves.remove_item(pos, nn_to_pos);
                    }
                    if self.atom_pos[nn_to_pos as usize].occ != 255
                        && self.atom_pos[nn_to_pos as usize].occ != 100
                    {
                        // greater than one because of neighbor moving in this spot
                        if self.atom_pos[pos as usize].cn_metal > 1 {
                            let (prev_e, future_e) = self.calc_energy_change_by_move(
                                nn_to_pos,
                                pos,
                                self.atom_pos[nn_to_pos as usize].occ,
                            );
                            let e_barr = alpha_energy::e_barrier(prev_e, future_e);
                            assert!(self.atom_pos[nn_to_pos as usize].occ != 100);
                            self.possible_moves.add_item(
                                nn_to_pos,
                                pos,
                                future_e - prev_e,
                                e_barr,
                                self.temperature,
                            );
                        }
                    }
                }
            }
            AddRemoveHow::RemoveAndAdd(add_atom_index, remove_atom_index) => {
                for nn_to_pos in self.gridstructure.nn[&pos] {
                    if self.atom_pos[nn_to_pos as usize].occ == 255 {
                        let (prev_e, future_e) = self.calc_energy_change_by_move(
                            pos,
                            nn_to_pos,
                            self.atom_pos[pos as usize].occ,
                        );
                        let e_barr = alpha_energy::e_barrier(prev_e, future_e);
                        self.possible_moves.update_k_if_move_exists(
                            pos,
                            nn_to_pos,
                            future_e - prev_e,
                            e_barr,
                            temperature,
                        );
                    }
                }
            }
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
                let (prev_e, future_e) =
                    self.calc_energy_change_by_move(full, empty, self.atom_pos[full as usize].occ);
                let e_barr = alpha_energy::e_barrier(prev_e, future_e);
                self.possible_moves.update_k_if_move_exists(
                    full,
                    empty,
                    future_e - prev_e,
                    e_barr,
                    self.temperature,
                );
            }
        }
    }

    pub fn redox_update_add_remove(&mut self, pos: u32, how: &AddRemoveHow) {
        match how {
            AddRemoveHow::Remove(_) => {
                self.add_or_remove.remove_item(pos);
                // self.add_or_remove
                //     .add_item(move_to, self.atom_pos[move_to as usize].cn_metal as u8);
                for x in self.gridstructure.nn[&pos] {
                    self.add_or_remove
                        .cond_update_cn(x, self.atom_pos[x as usize].cn_metal as u8);
                }
            }
            AddRemoveHow::RemoveAndAdd(_, _) => {
                self.add_or_remove.remove_item(pos);
                // self.add_or_remove
                //     .add_item(move_to, self.atom_pos[move_to as usize].cn_metal as u8);
                // for x in self.gridstructure.nn[&pos] {
                //     self.add_or_remove
                //         .cond_update_cn(x, self.atom_pos[x as usize].cn_metal as u8);
                // }
            }
        }
    }
}
