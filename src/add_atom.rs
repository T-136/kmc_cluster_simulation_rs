pub use super::alpha_energy;

impl crate::Simulation {
    pub fn add_atom(&mut self, atom: u32, is_recording_sections: bool, add_atom_type: u8) {
        self.atom_pos[atom as usize].occ = add_atom_type;

        self.onlyocc.insert(atom);

        if super::SAVE_ENTIRE_SIM || is_recording_sections {
            self.update_cn_dict(
                self.atom_pos[atom as usize].nn_support,
                self.atom_pos[atom as usize].cn_metal,
                false,
            );
            // self.cn_dict[self.atom_pos[move_from as usize].cn_metal] -= 1;
        }
        // println!("possible moves: {:?}", self.possible_moves.moves);
        // for o in self.gridstructure.nn[&move_from] {
        // let atom_type = self.atom_pos[atom as usize].occ as usize;

        for o in self.gridstructure.nn[&atom] {
            // for o in self.gridstructure.nn[&move_to] {
            if (super::SAVE_ENTIRE_SIM || is_recording_sections)
                && self.atom_pos[o as usize].occ != 255
                && self.atom_pos[o as usize].occ != 100
                && o != atom
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
            self.atom_pos[o as usize].nn_atom_type_count[add_atom_type as usize] += 1;
        }

        if super::SAVE_ENTIRE_SIM || is_recording_sections {
            self.update_cn_dict(
                self.atom_pos[atom as usize].nn_support,
                self.atom_pos[atom as usize].cn_metal,
                true,
            );
            // self.cn_dict[self.atom_pos[move_to as usize].cn_metal] += 1;
        }

        self.total_energy += self.e_change_add_move(atom, add_atom_type);
    }

    pub fn e_change_add_move(&mut self, atom: u32, add_atom_type: u8) -> f64 {
        let mut pos_nn_atom_type_ = self.atom_pos[atom as usize].nn_atom_type_count;
        self.alphas.e_one_atom(
            self.atom_pos[atom as usize].cn_metal,
            pos_nn_atom_type_,
            self.gridstructure.nn[&atom].iter().filter_map(|x| {
                if self.atom_pos[*x as usize].occ != 255 {
                    Some(alpha_energy::NnData {
                        cn_metal: self.atom_pos[*x as usize].cn_metal,
                        nn_atom_type_count_num_list: self.atom_pos[*x as usize].nn_atom_type_count,
                        atom_type: self.atom_pos[*x as usize].occ as usize,
                    })
                } else {
                    None
                }
            }),
            add_atom_type as usize,
            0,
            0,
        )
    }

    pub fn add_atom_update_total_k(&mut self, atom: u32, atom_type: u8) {
        // match self.energy {
        //     EnergyInput::LinearCn(_) | EnergyInput::Cn(_) => {
        // let i_slice = &self.gridstructure.nn_pair_no_intersec
        //     [&(move_from as u64 + ((move_to as u64) << 32))];
        for nn_to_pos in self.gridstructure.nn[&atom] {
            if self.atom_pos[nn_to_pos as usize].occ != 255
                && self.atom_pos[nn_to_pos as usize].occ != 100
            {
                self.possible_moves.remove_item(nn_to_pos, atom);
            }
            if self.atom_pos[nn_to_pos as usize].occ == 255 {
                // greater than one because of neighbor moving in this spot
                // if self.atom_pos[atom as usize].cn_metal > 1 {
                let (prev_e, future_e) =
                    self.calc_energy_change_by_move(atom, nn_to_pos, atom_type);
                let e_barr = alpha_energy::e_barrier(prev_e, future_e);
                assert!(self.atom_pos[nn_to_pos as usize].occ != 100);
                self.possible_moves.add_item(
                    atom,
                    nn_to_pos,
                    future_e - prev_e,
                    e_barr,
                    self.temperature,
                );
                // }
            }
        }
        //     }
        //     EnergyInput::LinearGcn(_) => todo!(),
        //     EnergyInput::Gcn(_) => todo!(),
        // }
    }

    pub fn add_atom_update_possible_moves(&mut self, atom: u32) {
        for nn_to_atom in self.gridstructure.nn[&atom] {
            if self.atom_pos[nn_to_atom as usize].occ != 255
                && self.atom_pos[nn_to_atom as usize].occ != 100
            {
                self.possible_moves.remove_item(nn_to_atom, atom);
            }
            if self.atom_pos[nn_to_atom as usize].occ == 255 {
                // greater than one because of neighbor moving in this spot
                if self.atom_pos[atom as usize].cn_metal > 1 {
                    let (prev_e, future_e) = self.calc_energy_change_by_move(
                        atom,
                        nn_to_atom,
                        self.atom_pos[atom as usize].occ,
                    );
                    let e_barr = alpha_energy::e_barrier(prev_e, future_e);
                    assert!(self.atom_pos[atom as usize].occ != 100);
                    self.possible_moves.add_item(
                        atom,
                        nn_to_atom,
                        future_e - prev_e,
                        e_barr,
                        self.temperature,
                    );
                }
            }
        }
    }

    pub fn add_atom_update_add_remove(
        &mut self,
        atom: u32,
        atom_type: u8,
        how: &super::add_remove::AddRemoveHow,
    ) {
        self.add_or_remove.cond_add_item(
            atom,
            self.atom_pos[atom as usize].cn_metal as u8,
            atom_type,
            self.temperature,
            how,
        );
        // self.add_or_remove
        //     .add_item(move_to, self.atom_pos[move_to as usize].cn_metal as u8);
        for x in self.gridstructure.nn[&atom] {
            self.add_or_remove.cond_update_cn(
                x,
                self.atom_pos[x as usize].cn_metal as u8,
                self.temperature,
            );
        }
    }
}
