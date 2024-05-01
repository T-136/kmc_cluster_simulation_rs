use core::panic;
use rand::distributions::{Distribution, Uniform};
use rand::rngs::SmallRng;
use rand::seq::SliceRandom;
use rand_distr::num_traits::Float;
use std::collections::HashMap;
use std::usize;

#[derive(Clone, Debug)]
pub struct ListDict {
    move_to_position: HashMap<u64, usize, ahash::RandomState>,
    pub moves: Vec<Move>, // [(from, to, energy_change)]
    pub moves_k: Vec<f64>,
    pub total_k: f64,
}

#[derive(Clone, Debug)]
pub struct Move {
    from: u32,
    pub to: u32,
    e_diff: f64,
    e_barr: f64,
}

impl ListDict {
    pub fn new(grid_size: [u32; 3]) -> ListDict {
        let largest_atom_position = grid_size[0] * grid_size[1] * grid_size[2] * 4;
        let item_to_position: HashMap<u64, usize, ahash::RandomState> = HashMap::default();
        // let item_to_position: HashMap<u64, usize, fnv::FnvBuildHasher> =
        //     fnv::FnvHashMap::with_capacity_and_hasher(32000, Default::default());
        ListDict {
            move_to_position: item_to_position,
            moves: Vec::with_capacity((largest_atom_position * 3) as usize),
            moves_k: Vec::new(),
            total_k: f64::INFINITY,
        }
    }

    pub fn calc_total_k_change(&mut self, temp: f64) {
        self.total_k = simd_sum(&self.moves_k);
        // self.total_k = self.moves_k.iter().sum::<f64>()
    }

    pub fn add_item(
        &mut self,
        move_from: u32,
        move_to: u32,
        e_diff: f64,
        e_barr: f64,
        temperature: f64,
    ) {
        match self
            .move_to_position
            .entry((move_from as u64 + ((move_to as u64) << 32)))
        {
            std::collections::hash_map::Entry::Vacant(e) => {
                let k = tst_rate_calculation(e_diff, e_barr, temperature);
                self.moves.push(Move {
                    from: move_from,
                    to: move_to,
                    e_diff,
                    e_barr,
                });
                self.moves_k.push(k);
                e.insert(self.moves.len() - 1);
                self.total_k += k;
            }
            _ => return,
        }
    }

    pub fn remove_item(&mut self, move_from: u32, move_to: u32) {
        if let Some(position) = self
            .move_to_position
            .remove(&(move_from as u64 + ((move_to as u64) << 32)))
        {
            let mmove = self.moves.pop().unwrap();
            let mmove_k = self.moves_k.pop().unwrap();
            if position != self.moves.len() {
                let old_move = std::mem::replace(&mut self.moves[position], mmove);
                let old_k = std::mem::replace(&mut self.moves_k[position], mmove_k);
                // self.moves[position] = mmove;
                self.move_to_position.insert(
                    (self.moves[position].from as u64 + ((self.moves[position].to as u64) << 32)),
                    position,
                );
                self.total_k -= old_k;
            } else {
                self.total_k -= mmove_k;
            }
        }
    }

    pub fn update_k_if_move_exists(
        &mut self,
        move_from: u32,
        move_to: u32,
        e_diff: f64,
        e_barr: f64,
        temperature: f64,
    ) {
        if let Some(position) = self
            .move_to_position
            .get(&(move_from as u64 + ((move_to as u64) << 32)))
        {
            // let old_energy = std::mem::replace(&mut self.moves[*position].energy, new_energy);

            let new_k = tst_rate_calculation(e_diff, e_barr, temperature);
            // println!("new_k: {} e: {}", new_k, new_energy);
            self.total_k += new_k;
            // let old_k = std::mem::replace(&mut self.moves[*position].k, new_k);

            // println!(
            //     "old_k: {} e: {}",
            //     self.moves[*position].k, self.moves[*position].energy
            // );
            self.total_k -= self.moves_k[*position];
            self.moves_k[*position] = new_k;
            self.moves[*position].e_diff = e_diff;
            self.moves[*position].e_barr = e_barr;
        }
    }

    pub fn choose_ramdom_move_kmc(
        &self,
        rng_choose: &mut SmallRng,
        temp: f64,
    ) -> Option<(u32, u32, f64, f64, f64, f64)> {
        // self.calc_total_k_change(temp);
        let between = Uniform::new_inclusive(0., 1.);
        let mut k_time_rng = between.sample(rng_choose) * self.total_k;
        // println!(
        //     "ktot: {} krng: {}",
        //     format!("{:e}", self.total_k),
        //     format!("{:e}", k_time_rng),
        // );
        let mut cur_k = 0_f64;
        let mut res: Option<(u32, u32, f64, f64, f64, f64)> = None;

        for (iter, k) in self.moves_k.iter().enumerate() {
            cur_k += k;
            if k_time_rng <= cur_k {
                return Some((
                    self.moves[iter].from,
                    self.moves[iter].to,
                    self.moves[iter].e_diff,
                    self.moves[iter].e_barr,
                    self.total_k,
                    *k,
                ));
            }
        }
        return None;
        // for (i, mmove) in self.iter().enumerate() {

        // return res;
        // }
        // let calc_tot_k = self.iter().map(|mmove| mmove.k).sum::<f64>();
        // println!(
        //     "clac_tot_k/tot_k: {}, # moves: {}",
        //     format!("{}", calc_tot_k / self.total_k),
        //     format!("{}", self.iter().count()),
        // );
        // // println!("tot_k: {}", self.total_k);
        // res
    }

    pub fn iter(&self) -> std::slice::Iter<'_, Move> {
        self.moves.iter()
    }

    // pub fn contains(&self, move_from: u32, move_to: u32) -> bool {
    //     self.item_to_position
    //         .contains_key(&(move_from as u64 + ((move_to as u64) << 32)))
    // }

    // pub fn remove_by_index(&mut self, index: usize) {
    //     self.item_to_position.remove(&self.items.swap_remove(index));
    // }

    // pub fn drain_filter(&mut self, cn: &Vec<usize>, move_from: &u32, move_to: &u32) {
    //     let mut i = 0;
    //     while i < self.items.len() {
    //         let (o, u) = self.items[i];
    //         if (cn[u as usize] == 0 || &o == move_from || &u == move_to) {
    //             let (move_from, move_to) = self.items.remove(i);
    //             self.item_to_position
    //                 .remove(&(move_from as u64 + ((move_to as u64) << 32)));
    //         } else {
    //             i += 1;
    //         }
    //     }
    // }
    // pub fn filter(&self) {
    //     self.items.iter().filter()
    // }

    // pub fn iter_mut(self) -> std::vec::IntoIter<(u64, u64)> {
    //     self.items.into_iter()
    // }

    pub fn _len(&self) -> usize {
        self.moves.len()
    }
}

fn tst_rate_calculation(e_diff: f64, e_barr: f64, temperature: f64) -> f64 {
    // let e_use = if e_diff.is_negative() { 0. } else { e_diff };
    // println!("barr {}", e_barr);
    // println!("diff {}", e_diff);
    const KB_joul: f64 = 1.380649e-23;
    const h_joul: f64 = 6.62607015e-34;
    const KB_DIV_H: f64 = KB_joul / h_joul;
    const KB_eV: f64 = 8.6173324e-5;
    (KB_joul * temperature / h_joul) * ((-e_barr) / (KB_eV * temperature)).exp()
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

// fn slow_sum(x: &[f32]) -> f32 {
//     assert!(x.len() % 4 == 0);
//     let mut sum: f32 = 0.;
//     for i in (0..x.len()).step_by(4) {
//         sum += f32x4::from_slice_unaligned(&x[i..]).sum();
//     }
//     sum
// }
