use super::SAVE_ENTIRE_SIM;
use crate::alpha_energy;
use crate::atom_change;
use crate::atom_change::AtomChangeHow;
use crate::buckets_linear::ItemEnum;
use crate::grid_structure::{
    NNN_PAIR_NO_INTERSEC_NUMBER, NN_PAIR_NO_INTERSEC_NUMBER, NN_PAIR_ONLY_INTERSEC_NUMBER,
};
use crate::moves;
use crate::NUM_ATOM_TYPES;

#[derive(Clone, Debug)]
pub struct Move {
    pub from: u32,
    pub to: u32,
    pub e_diff: f64,
    pub e_barr: f64,
    pub k: f64,
}

impl Move {
    pub fn new(from: u32, to: u32, e_barr: f64, e_diff: f64, temperature: f64) -> Move {
        Move {
            from,
            to,
            e_diff,
            e_barr,
            k: tst_rate_calculation(e_diff, e_barr, temperature),
        }
    }
}

pub fn tst_rate_calculation(e_diff: f64, e_barr: f64, temperature: f64) -> f64 {
    // let e_use = if e_diff.is_negative() { 0. } else { e_diff };
    // println!("barr {}", e_barr);
    // println!("diff {}", e_diff);
    const KB_joul: f64 = 1.380649e-23;
    const h_joul: f64 = 6.62607015e-34;
    const KB_DIV_H: f64 = KB_joul / h_joul;
    const KB_eV: f64 = 8.6173324e-5;
    (KB_joul * temperature / h_joul) * ((-e_barr) / (KB_eV * temperature)).exp()
}

impl crate::Simulation {
    pub fn perform_move(
        &mut self,
        move_from: u32,
        move_to: u32,
        e_diff: f64,
        is_recording_sections: bool,
    ) {
        let (from_change, to_change, inter) =
            no_int_nn_from_move(move_from, move_to, &self.gridstructure.nn_pair_no_intersec);

        self.atom_pos[move_to as usize].occ = self.atom_pos[move_from as usize].occ; // covers different alloys also
        self.atom_pos[move_from as usize].occ = 255;

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
        let atom_type = self.atom_pos[move_to as usize].occ as usize;
        for o in from_change {
            if (SAVE_ENTIRE_SIM || is_recording_sections)
                && self.atom_pos[o as usize].occ != 255
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
                && self.atom_pos[o as usize].occ != 255
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

    pub fn update_moves(&mut self, move_from: u32, move_to: u32) {
        for pos in self.gridstructure.surrounding_moves[&(std::cmp::min(move_from, move_to) as u64
            + ((std::cmp::max(move_to, move_from) as u64) << 32))]
            .iter()
        {
            if self.atom_pos[pos.0 as usize].occ != 255
                && self.atom_pos[pos.0 as usize].occ != 100
                && self.atom_pos[pos.1 as usize].occ == 255
            {
                let (prev_e, future_e) = self.calc_energy_change_by_move(
                    pos.0,
                    pos.1,
                    self.atom_pos[pos.0 as usize].occ,
                );
                let e_barr = alpha_energy::e_barrier(prev_e, future_e);
                let mmove = Move::new(pos.0, pos.1, future_e - prev_e, e_barr, self.temperature);
                assert!(self.atom_pos[pos.0 as usize].occ != 255);
                self.possible_moves
                    .update_k_if_item_exists(ItemEnum::Move(mmove));
            }
            if self.atom_pos[pos.1 as usize].occ != 100
                && self.atom_pos[pos.1 as usize].occ != 255
                && self.atom_pos[pos.0 as usize].occ == 255
            {
                let (prev_e, future_e) = self.calc_energy_change_by_move(
                    pos.1,
                    pos.0,
                    self.atom_pos[pos.1 as usize].occ,
                );
                let e_barr = alpha_energy::e_barrier(prev_e, future_e);
                let mmove = Move::new(pos.1, pos.0, future_e - prev_e, e_barr, self.temperature);
                assert!(self.atom_pos[pos.1 as usize].occ != 255);
                self.possible_moves
                    .update_k_if_item_exists(ItemEnum::Move(mmove));
            }
        }
    }

    pub fn update_possible_moves(&mut self, move_from: u32, move_to: u32) {
        self.possible_moves.remove_move(move_from, move_to);
        for nn_to_from in self.gridstructure.nn[&move_from] {
            if self.atom_pos[nn_to_from as usize].occ == 255 {
                self.possible_moves.remove_move(move_from, nn_to_from);
            }
            if self.atom_pos[nn_to_from as usize].occ != 255
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
                    assert!(self.atom_pos[move_to as usize].occ != 255);
                    let mmove = moves::Move::new(
                        nn_to_from,
                        move_from,
                        future_e - prev_e,
                        e_barr,
                        self.temperature,
                    );
                    assert!(self.atom_pos[nn_to_from as usize].occ != 255);
                    self.possible_moves.cond_add_item(ItemEnum::Move(mmove));
                }
            }
        }

        for nn_to_to in self.gridstructure.nn[&move_to] {
            if self.atom_pos[nn_to_to as usize].occ != 255
                && self.atom_pos[nn_to_to as usize].occ != 100
            {
                self.possible_moves.remove_move(nn_to_to, move_to);
            }
            if self.atom_pos[nn_to_to as usize].occ == 255 {
                // greater than one because of neighbor moving in this spot
                if self.atom_pos[nn_to_to as usize].cn_metal > 1 {
                    let (prev_e, future_e) = self.calc_energy_change_by_move(
                        move_to,
                        nn_to_to,
                        self.atom_pos[move_to as usize].occ,
                    );
                    let e_barr = alpha_energy::e_barrier(prev_e, future_e);
                    assert!(self.atom_pos[move_to as usize].occ != 255);
                    let mmove = moves::Move::new(
                        move_to,
                        nn_to_to,
                        e_barr,
                        future_e - prev_e,
                        self.temperature,
                    );
                    assert!(self.atom_pos[move_to as usize].occ != 255);
                    self.possible_moves.cond_add_item(ItemEnum::Move(mmove));
                }
            }
        }
    }

    pub fn update_add_remove(&mut self, move_from: u32, move_to: u32, how: &AtomChangeHow) {
        let (from_change, to_change, inter) =
            no_int_nn_from_move(move_from, move_to, &self.gridstructure.nn_pair_no_intersec);
        // if add_remove::REMOVE_ATOM || add_remove::EXCHANGE_ATOM {
        let pos_change = match how {
            AtomChangeHow::Add => {
                self.possible_moves.remove_add_remove(move_to);
                atom_change::AtomPosChange::new(
                    move_from,
                    self.atom_pos[move_to as usize].cn_metal as u8,
                    self.atom_pos[move_to as usize].occ as u8,
                    self.temperature,
                    how.clone(),
                )
            }
            AtomChangeHow::Remove => {
                self.possible_moves.remove_add_remove(move_from);
                atom_change::AtomPosChange::new(
                    move_to,
                    self.atom_pos[move_to as usize].cn_metal as u8,
                    self.atom_pos[move_to as usize].occ as u8,
                    self.temperature,
                    how.clone(),
                )
            }
            AtomChangeHow::Exchange => todo!(),
            AtomChangeHow::RemoveAndAdd => todo!(),
        };
        if let Some(pos_change) = pos_change {
            self.possible_moves
                .cond_add_item(ItemEnum::AddOrRemove(pos_change));
        }
        for x in from_change {
            let pos_change = atom_change::AtomPosChange::new(
                x,
                self.atom_pos[x as usize].cn_metal as u8,
                self.atom_pos[x as usize].occ as u8,
                self.temperature,
                how.clone(),
            );
            if let Some(pos_change) = pos_change {
                self.possible_moves
                    .update_k_if_item_exists(ItemEnum::AddOrRemove(pos_change));
            }
        }
        for x in to_change {
            let pos_change = atom_change::AtomPosChange::new(
                x,
                self.atom_pos[x as usize].cn_metal as u8,
                self.atom_pos[x as usize].occ as u8,
                self.temperature,
                how.clone(),
            );
            if let Some(pos_change) = pos_change {
                self.possible_moves
                    .update_k_if_item_exists(ItemEnum::AddOrRemove(pos_change));
            }
        }
    }

    pub fn calc_energy_change_by_move(
        &self,
        move_from: u32,
        move_to: u32,
        atom_typ: u8,
    ) -> (f64, f64) {
        // println!("{}", atom_typ_index);
        // println!("{:?}", self.atom_pos.neigboring_atom_type_count);
        //
        let (from_change_nn, to_change_nn, inter_nn) =
            no_int_nn_from_move(move_from, move_to, &self.gridstructure.nn_pair_no_intersec);

        let mut nn_atom_type_count_tst = [0; NUM_ATOM_TYPES];
        let mut cn_tst = 0;
        let mut from_nn_atom_type_no_tst = self.atom_pos[move_from as usize].nn_atom_type_count;
        let mut to_nn_atom_type_count = self.atom_pos[move_to as usize].nn_atom_type_count;
        to_nn_atom_type_count[atom_typ as usize] -= 1;

        // inter_nn.iter().for_each(|x| {
        for x in inter_nn {
            if self.atom_pos[x as usize].occ != 255 && x != move_from {
                nn_atom_type_count_tst[self.atom_pos[x as usize].occ as usize] += 1;
                // from_nn_atom_type_no_tst[self.atom_pos[x as usize].occ as usize - 1] -= 1;
                // to_nn_atom_type_count[self.atom_pos[x as usize].occ as usize - 1] -= 1;
                cn_tst += 1;
            }
        }

        let tst_e = self.alphas.e_one_atom_tst(
            cn_tst,
            nn_atom_type_count_tst,
            self.atom_pos[move_from as usize].occ as usize,
            0,
            0,
        );

        let prev_e = self.alphas.e_one_atom(
            self.atom_pos[move_from as usize].cn_metal,
            // self.atom_pos[move_from as usize].nn_atom_type_count,
            from_nn_atom_type_no_tst,
            from_change_nn.iter().filter_map(|x| {
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
            self.atom_pos[move_from as usize].occ,
            0,
            0,
        );

        let future_e = self.alphas.e_one_atom(
            self.atom_pos[move_to as usize].cn_metal - 1,
            to_nn_atom_type_count,
            to_change_nn.iter().filter_map(|x| {
                if self.atom_pos[*x as usize].occ != 255 && *x != move_from {
                    let mut cn_metal_future = self.atom_pos[*x as usize].cn_metal;
                    let mut nn_atom_type_count_future =
                        self.atom_pos[*x as usize].nn_atom_type_count;
                    nn_atom_type_count_future[atom_typ as usize] += 1;
                    cn_metal_future += 1;

                    Some(alpha_energy::NnData {
                        cn_metal: cn_metal_future,
                        nn_atom_type_count_num_list: nn_atom_type_count_future,
                        atom_type: self.atom_pos[*x as usize].occ as usize,
                    })
                } else {
                    None
                }
            }),
            self.atom_pos[move_from as usize].occ,
            0,
            0,
        );

        // println!("prev_e: {}, future_e: {} tst_e {}", prev_e, future_e, tst_e);
        (prev_e - tst_e, future_e - tst_e)
    }
}

fn no_int_nn_from_move(
    move_from: u32,
    move_to: u32,
    nn_pair_no_intersec: &std::collections::HashMap<
        u64,
        (
            [u32; NN_PAIR_NO_INTERSEC_NUMBER],
            [u32; NN_PAIR_NO_INTERSEC_NUMBER],
            [u32; NN_PAIR_ONLY_INTERSEC_NUMBER],
        ),
        fnv::FnvBuildHasher,
    >,
) -> (
    [u32; NN_PAIR_NO_INTERSEC_NUMBER],
    [u32; NN_PAIR_NO_INTERSEC_NUMBER],
    [u32; NN_PAIR_ONLY_INTERSEC_NUMBER],
) {
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
    nnn_pair_no_intersec: &std::collections::HashMap<
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

// pub unsafe fn update_k_if_move_exists_par(
//     ptr: *mut ListDict,
//     move_from: u32,
//     move_to: u32,
//     new_energy: i64,
//     temperature: f64,
// ) -> Option<f64> {
//     if let Some(position) = (*ptr)
//         .move_to_position
//         .get(&(move_from as u64 + ((move_to as u64) << 32)))
//     {
//         // let old_energy = std::mem::replace(&mut self.moves[*position].energy, new_energy);
//
//         let new_k = tst_rate_calculation(new_energy, temperature);
//         // println!("new_k: {} e: {}", new_k, new_energy);
//         let total_k_change = new_k - (*ptr).moves[*position].k;
//         // (*ptr).total_k += new_k;
//         // (*ptr).total_k -= (*ptr).moves[*position].k;
//         (*ptr).moves[*position].k = new_k;
//         (*ptr).moves[*position].energy = new_energy;
//         Some(total_k_change)
//     } else {
//         None
//     }
// }
// pub struct WrapperListDict(pub *mut ListDict);
// unsafe impl Send for ListDict {}
// unsafe impl Sync for ListDict {}
// unsafe impl Send for WrapperListDict {}
// unsafe impl Sync for WrapperListDict {}

use std::convert::TryInto;

const LANES: usize = 16;

pub fn simd_sum(values: &[f64]) -> f64 {
    let chunks = values.chunks_exact(LANES);
    let remainder = chunks.remainder();

    let sum = chunks.fold([0.0_f64; LANES], |mut acc, chunk| {
        let chunk: [f64; LANES] = chunk.try_into().unwrap();
        for i in 0..LANES {
            acc[i] += chunk[i];
        }
        acc
    });

    let remainder: f64 = remainder.iter().copied().sum();
    let reduced: f64 = sum.iter().copied().take(LANES).sum();
    reduced + remainder
}

pub fn simd_sum_to(values: &[f64], cond: f64) -> f64 {
    let chunks = values.chunks_exact(LANES);
    let remainder = chunks.remainder();

    let mut i = 0;
    let mut sum = 0.;

    let sum = chunks.fold([0.0_f64; LANES], |mut acc, chunk| {
        let chunk: [f64; LANES] = chunk.try_into().unwrap();
        for i in 0..LANES {
            acc[i] += chunk[i];
        }
        acc
    });

    let remainder: f64 = remainder.iter().copied().sum();
    let reduced: f64 = sum.iter().copied().take(LANES).sum();
    reduced + remainder
}
