use core::panic;
use rand::distributions::{Distribution, Uniform};
use rand::rngs::SmallRng;
use rand::seq::SliceRandom;
use rand_distr::num_traits::Float;
use std::collections::HashMap;
use std::sync::Arc;
use std::usize;

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
