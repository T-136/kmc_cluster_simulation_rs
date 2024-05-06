use rand::distributions::{Distribution, Uniform};
use rand::prelude::*;
use rand::rngs::SmallRng;
use std::collections::HashMap;

use crate::add_remove;

const CN_FOR_INV: u8 = 12;
// const E_RATIO: f64 = 0.20; 400K
const E_RATIO_BARR: f64 = 0.10;

#[derive(Clone, Debug)]
pub struct AddOrRemove {
    atom_to_position: HashMap<u32, usize, ahash::RandomState>,
    pub atoms: Vec<Removable>, // [(from, to, energy_change)]
    pub total_inv_cn: u64,
    pub total_k: f64,
    // potential: f64,
}

#[derive(Clone, Debug)]
pub struct Removable {
    pos: u32,
    inv_cn: u8,
    k: f64,
}

impl AddOrRemove {
    pub fn new() -> AddOrRemove {
        let item_to_position: HashMap<u32, usize, ahash::RandomState> = HashMap::default();
        // let item_to_position: HashMap<u64, usize, fnv::FnvBuildHasher> =
        //     fnv::FnvHashMap::with_capacity_and_hasher(32000, Default::default());
        AddOrRemove {
            atom_to_position: item_to_position,
            atoms: Vec::new(),
            total_inv_cn: 0_u64,
            total_k: 0.,
        }
    }

    pub fn calc_total_cn_change(&mut self) {
        self.total_inv_cn = self.atoms.iter().map(|x| x.inv_cn as u64).sum::<u64>();
        self.total_k = self.atoms.iter().map(|x| x.k).sum::<f64>();
    }

    pub fn cond_add_item(
        &mut self,
        pos: u32,
        cn: u8,
        atom_type: u8,
        temperature: f64,
        how: &add_remove::AddRemoveHow,
    ) {
        if cn >= 11 {
            return;
        }
        match how {
            add_remove::AddRemoveHow::Remove(remove_atom_type)
            | add_remove::AddRemoveHow::RemoveAndAdd(remove_atom_type, _) => {
                if atom_type == *remove_atom_type {
                    match self.atom_to_position.entry(pos) {
                        std::collections::hash_map::Entry::Vacant(e) => {
                            let k = tst_rate_calculation(cn as f64 * E_RATIO_BARR, temperature);
                            self.atoms.push(Removable {
                                pos,
                                inv_cn: CN_FOR_INV - cn,
                                k,
                            });
                            e.insert(self.atoms.len() - 1);
                            self.total_k += k;
                            // self.total_inv_cn += (CN_FOR_INV - cn) as u64;
                        }
                        _ => return,
                    }
                }
            }
        }
    }

    pub fn remove_item(&mut self, pos: u32) {
        if let Some(position) = self.atom_to_position.remove(&pos) {
            let atom = self.atoms.pop().unwrap();
            if position != self.atoms.len() {
                let old_move = std::mem::replace(&mut self.atoms[position], atom);
                // self.moves[position] = mmove;
                self.atom_to_position
                    .insert(self.atoms[position].pos, position);
                // self.total_inv_cn -= old_move.inv_cn as u64;
                self.total_k -= old_move.k;
            } else {
                self.total_k -= atom.k;
                // self.total_inv_cn -= atom.inv_cn as u64;
            }
        }
    }

    pub fn cond_update_cn(&mut self, pos: u32, cn: u8, temperature: f64) {
        let inv_cn = CN_FOR_INV - cn;
        if let Some(position) = self.atom_to_position.get(&pos) {
            let k = tst_rate_calculation((cn) as f64 * E_RATIO_BARR, temperature);
            // self.total_inv_cn += inv_cn as u64;
            // self.total_inv_cn -= self.atoms[*position].inv_cn as u64;
            // self.atoms[*position].inv_cn = inv_cn;
            self.total_k -= self.atoms[*position].k;
            self.atoms[*position].k = k;
            self.total_k += k;
        }
    }
    // pub fn cond_update_type(&mut self, pos: u32, cn: u8) {
    //     let inv_cn = CN_FOR_INV - cn;
    //     if let Some(position) = self.atom_to_position.get(&pos) {
    //         self.total_inv_cn += inv_cn as u64;
    //         self.total_inv_cn -= self.atoms[*position].inv_cn as u64;
    //         self.atoms[*position].inv_cn = inv_cn;
    //     }
    // }

    pub fn choose_ramdom_atom_to_remove(
        &self,
        rng_choose: &mut SmallRng,
    ) -> Option<(u32, f64, f64)> {
        // self.calc_total_k_change(temp);
        let between = Uniform::new_inclusive(0., self.total_k);
        let mut k_time_rng = between.sample(rng_choose);
        // println!(
        //     "ktot: {} krng: {}",
        //     format!("{:e}", self.total_k),
        //     format!("{:e}", k_time_rng),
        // );
        let mut cur_cn = 0.;
        let mut res: Option<(u32, u32, f64, f64, f64, f64)> = None;

        for atom in self.atoms.iter() {
            cur_cn += atom.k;
            // println!("{}", cur_cn);
            if k_time_rng <= cur_cn {
                return Some((atom.pos, atom.k, self.total_k));
            }
        }
        return None;
    }

    pub fn iter(&self) -> std::slice::Iter<'_, Removable> {
        self.atoms.iter()
    }

    pub fn _len(&self) -> usize {
        self.atoms.len()
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
