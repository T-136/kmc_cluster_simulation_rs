use rand::distributions::{Distribution, Uniform};
use rand::rngs::SmallRng;
use rand::seq::SliceRandom;
use std::collections::HashMap;

#[derive(Clone)]
pub struct ListDict {
    move_to_position: HashMap<u64, usize, ahash::RandomState>,
    pub moves: Vec<(u32, u32, Option<i64>)>, // [(from, to, energy_change)]
    total_k: Option<f64>,
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
            total_k: None,
        }
    }

    pub fn calc_total_k_change(&mut self, temp: f64) {
        self.total_k = Some(
            self.iter()
                .map(|(_, _, energy_change)| {
                    // println!(
                    //     "energy_change: {} K: {}",
                    //     energy_change.unwrap(),
                    //     arrhenius_equation(energy_change.unwrap(), temp)
                    // );
                    arrhenius_equation(energy_change.unwrap(), temp)
                })
                .sum::<f64>(),
        )
    }

    pub fn add_item(&mut self, move_from: u32, move_to: u32, energy_change: Option<i64>) {
        match self
            .move_to_position
            .entry((move_from as u64 + ((move_to as u64) << 32)))
        {
            std::collections::hash_map::Entry::Vacant(e) => {
                self.moves.push((move_from, move_to, energy_change));
                e.insert(self.moves.len() - 1);
            }
            _ => return,
        }
    }

    pub fn remove_item(&mut self, move_from: u32, move_to: u32) {
        if let Some(position) = self
            .move_to_position
            .remove(&(move_from as u64 + ((move_to as u64) << 32)))
        {
            let (move_from, move_to, energy_change) = self.moves.pop().unwrap();
            if position != self.moves.len() {
                self.moves[position] = (move_from, move_to, energy_change);
                self.move_to_position
                    .insert((move_from as u64 + ((move_to as u64) << 32)), position);
            }
        }
    }

    pub fn choose_random_move_mc(&self, rng_choose: &mut SmallRng) -> (u32, u32, Option<i64>) {
        self.moves.choose(rng_choose).unwrap().clone()
    }

    pub fn choose_ramdom_move_kmc(
        &mut self,
        rng_choose: &mut SmallRng,
        temp: f64,
    ) -> Option<(u32, u32, i64, f64)> {
        self.calc_total_k_change(temp);
        let between = Uniform::new_inclusive(0., 1.);
        let k_time_rng = between.sample(rng_choose) * self.total_k.unwrap();
        let mut cur_k = 0_f64;
        let mut res: Option<(u32, u32, i64, f64)> = None;
        self.iter().for_each(|(from, to, energy_change)| {
            cur_k += arrhenius_equation(energy_change.expect("no energy change in ListDict"), temp);
            if cur_k >= k_time_rng {
                res = Some((*from, *to, energy_change.unwrap(), self.total_k.unwrap()))
            }
        });
        res
    }

    pub fn iter(&self) -> std::slice::Iter<'_, (u32, u32, Option<i64>)> {
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

fn arrhenius_equation(energy: i64, temperature: f64) -> f64 {
    let e_use = if energy.is_negative() { 0 } else { energy };
    const e_barrier: i64 = 1;
    const KB_joul: f64 = 1.380649e-23;
    // const R: f64 = 8.31446261815324;
    const h_joul: f64 = 6.62607015e-34;
    const KB_eV: f64 = 8.6173324e-5;
    (KB_joul * temperature / h_joul)
        * (-(e_use + e_barrier) as f64 / (KB_eV * temperature * 1000.)).exp()
}

// #[cfg(test)]
// mod tests {
//     use super::*;
//
//     #[test]
//     fn test_arrhenius_eq() {
//         let k = arrhenius_equation(-360, 300.);
//         println!("{}", k);
//         assert_eq!(k, -100.);
//     }
// }
