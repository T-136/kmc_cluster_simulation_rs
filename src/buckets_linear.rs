use core::panic;
use rand::distributions::{Distribution, Uniform};
use rand::rngs::SmallRng;
use std::collections::HashMap;

use super::atom_change;
use super::moves;

const MAX_EDIT_COUNTER: i32 = 1000;

#[derive(Clone, Debug)]
pub enum ItemEnum {
    Move(moves::Move),
    AddOrRemove(atom_change::AtomPosChange),
}

impl ItemEnum {
    fn get_id(&self) -> u64 {
        match self {
            ItemEnum::Move(mmove) => (mmove.from as u64 + ((mmove.to as u64) << 32)),
            ItemEnum::AddOrRemove(add_or_remove) => {
                // add_or_remove.pos as u64 + ((add_or_remove.how as u64) << 32)
                add_or_remove.pos as u64
            }
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct Bucket {
    edit_counter: i32,
    pub bucket_power: i32,
    pub own_k: f64,
    pub items: Vec<ItemEnum>,
    pub k_s: Vec<f64>,
}

impl Bucket {
    fn cond_update_bucket_k(&mut self, new_k: f64) {
        self.edit_counter += 1;
        if self.own_k * 0.2 <= new_k || self.own_k < 0.001 || self.edit_counter >= MAX_EDIT_COUNTER
        {
            self.edit_counter = 0;
            self.own_k = self.k_s.iter().sum::<f64>();
        }
    }
}

#[derive(Clone, Debug)]
pub struct ItemIndexes {
    pub power_of_k: i32,
    pub vec_index: usize,
}

#[derive(Clone, Debug)]
pub struct Buckets {
    pub edit_counter: i32,
    pub bucket_power_to_pos: HashMap<i32, usize, ahash::RandomState>,
    // key = from + (to << 32)
    pub move_to_position: HashMap<u64, ItemIndexes, ahash::RandomState>,
    // key = pos + (enum discrimenent << 32)
    pub add_remove_to_position: HashMap<u64, ItemIndexes, ahash::RandomState>,
    pub total_k: f64,
    pub buckets_list: Vec<Bucket>,
}

impl Buckets {
    pub fn new() -> Buckets {
        let bucket_power_to_pos: HashMap<i32, usize, ahash::RandomState> = HashMap::default();
        let move_to_position: HashMap<u64, ItemIndexes, ahash::RandomState> = HashMap::default();
        let add_remove_to_position: HashMap<u64, ItemIndexes, ahash::RandomState> =
            HashMap::default();
        Buckets {
            edit_counter: 0,
            bucket_power_to_pos,
            move_to_position,
            add_remove_to_position,
            total_k: 0.,
            buckets_list: Vec::new(),
        }
    }

    pub fn get(&mut self, from: u32, to: u32) -> Option<ItemEnum> {
        if let Some(item_indexes) = self
            .move_to_position
            .get(&(from as u64 + ((to as u64) << 32)))
        {
            let bucket_index = self
                .bucket_power_to_pos
                .get(&item_indexes.power_of_k)
                .unwrap();

            let bucket = &mut self.buckets_list[*bucket_index];

            let item = bucket.items[item_indexes.vec_index].clone();
            return Some(item);
        }
        None
    }

    fn cond_update_ks(&mut self, new_k: f64, bucket_index: usize) {
        self.buckets_list[bucket_index].cond_update_bucket_k(new_k);
        self.edit_counter += 1;
        if self.total_k * 0.2 <= new_k
            || self.total_k < 0.001
            || self.edit_counter >= MAX_EDIT_COUNTER
        {
            self.edit_counter = 0;
            self.total_k = self.buckets_list.iter().map(|item| item.own_k).sum::<f64>();
        }
    }

    pub fn cond_add_item(&mut self, item: ItemEnum, k: f64) {
        let power_of_k = k.log2().ceil() as i32;
        let bucket_index = if let Some(x) = self.bucket_power_to_pos.get(&power_of_k) {
            *x
        } else {
            self.cond_add_bucket(power_of_k)
        };

        match item {
            ItemEnum::Move(_) => match self.move_to_position.entry(item.get_id()) {
                std::collections::hash_map::Entry::Occupied(_) => return,
                std::collections::hash_map::Entry::Vacant(e) => e.insert(ItemIndexes {
                    power_of_k,
                    vec_index: self.buckets_list[bucket_index].items.len(),
                }),
            },
            ItemEnum::AddOrRemove(_) => match self.add_remove_to_position.entry(item.get_id()) {
                std::collections::hash_map::Entry::Occupied(_) => return,
                std::collections::hash_map::Entry::Vacant(e) => e.insert(ItemIndexes {
                    power_of_k,
                    vec_index: self.buckets_list[bucket_index].items.len(),
                }),
            },
        };
        self.total_k += k;
        self.buckets_list[bucket_index].own_k += k;
        self.buckets_list[bucket_index].items.push(item);
        self.buckets_list[bucket_index].k_s.push(k);
        self.cond_update_ks(k, bucket_index);
    }

    pub fn cond_add_bucket(&mut self, power_of_k: i32) -> usize {
        if let Some(bucket_index) = self.bucket_power_to_pos.get(&power_of_k) {
            return *bucket_index;
        }
        println!("bucket count {}", self.buckets_list.len());
        let bucket = Bucket {
            bucket_power: power_of_k,
            ..Default::default()
        };
        self.buckets_list.push(bucket);
        self.bucket_power_to_pos
            .insert(power_of_k, self.buckets_list.len() - 1);
        self.buckets_list.len() - 1
    }

    pub fn remove_move(&mut self, from: u32, to: u32) -> bool {
        if let Some(item_indexes) = self
            .move_to_position
            .remove(&(from as u64 + ((to as u64) << 32)))
        {
            let bucket_index = self
                .bucket_power_to_pos
                .get(&item_indexes.power_of_k)
                .unwrap();

            let bucket = &mut self.buckets_list[*bucket_index];

            let poped_item = bucket.items.pop().unwrap();
            let poped_k = bucket.k_s.pop().unwrap();
            let k = if bucket.items.len() != item_indexes.vec_index {
                let id = poped_item.get_id();
                match poped_item {
                    ItemEnum::Move(_) => {
                        self.move_to_position.insert(id, item_indexes.clone());
                    }
                    ItemEnum::AddOrRemove(_) => {
                        self.add_remove_to_position.insert(id, item_indexes.clone());
                    }
                }
                bucket.items[item_indexes.vec_index] = poped_item;
                let old_k =
                    std::mem::replace(&mut bucket.k_s[item_indexes.vec_index], poped_k);
                self.move_to_position.insert(
                    (bucket.items[item_indexes.vec_index].get_id()),
                    item_indexes,
                );
                old_k
            } else {
                poped_k
            };
            self.total_k -= k;
            bucket.own_k -= k;
            self.cond_update_ks(k, *bucket_index);
            return true;
        }
        false
    }

    pub fn remove_add_remove(&mut self, pos: u32) -> bool {
        if let Some(item_indexes) = self.add_remove_to_position.remove(&(pos as u64)) {
            let bucket_index = self
                .bucket_power_to_pos
                .get(&item_indexes.power_of_k)
                .unwrap();

            let bucket = &mut self.buckets_list[*bucket_index];

            let poped_item = bucket.items.pop().unwrap();
            let poped_k = bucket.k_s.pop().unwrap();
             let k = if bucket.items.len() != item_indexes.vec_index {
                let id = poped_item.get_id();
                match poped_item {
                    ItemEnum::Move(_) => {
                        self.move_to_position.insert(id, item_indexes.clone());
                    }
                    ItemEnum::AddOrRemove(_) => {
                        self.add_remove_to_position.insert(id, item_indexes.clone());
                    }
                }
                    bucket.items[item_indexes.vec_index] = poped_item;
                let old_k =
                    std::mem::replace(&mut bucket.k_s[item_indexes.vec_index], poped_k);
                // self.add_remove_to_position.insert(
                //     (bucket.items[item_indexes.vec_index].get_id()),
                //     item_indexes,
                // );
                old_k
            } else {
                poped_k
            };


            self.total_k -= k;
            bucket.own_k -= k;
            self.cond_update_ks(k, *bucket_index);
            return true;
        }
        false
    }

    pub fn update_k_if_move_exists(&mut self, move_from: u32, move_to: u32, k: f64) {
        //problem: move needs change bucket if k change makes it necessary
        if let Some(item_indexes) = self
            .move_to_position
            .get(&(move_from as u64 + ((move_to as u64) << 32)))
        {
            let bucket_index =
                if let Some(x) = self.bucket_power_to_pos.get(&item_indexes.power_of_k) {
                    *x
                } else {
                    panic!("move is supposed to be in a bucket that does not exist");
                };

            let bucket = &mut self.buckets_list[bucket_index];
            bucket.own_k -= bucket.k_s[item_indexes.vec_index];
            bucket.own_k -= k;
            bucket.k_s[item_indexes.vec_index] = k;
        }
    }

    pub fn update_k_if_item_exists(&mut self, item: ItemEnum, k: f64) {
        match item {
            ItemEnum::Move(mmove) => {
                //todo: unnecessary remove and insert could change k in place, see fn update_k_if_move_exists
                if self.remove_move(mmove.from, mmove.to) {
                    self.cond_add_item(ItemEnum::Move(mmove), k);
                }
            }
            ItemEnum::AddOrRemove(add_or_remove) => {
                if self.remove_add_remove(add_or_remove.pos) {
                    self.cond_add_item(ItemEnum::AddOrRemove(add_or_remove), k);
                }
            }
        }
    }

    fn pick_from_bucket(
        &self,
        bucket_pick: &mut SmallRng,
        coin_toss: &mut SmallRng,
        bucket: &Bucket,
    ) -> Option<ItemEnum> {
        let pick = Uniform::new_inclusive(0, bucket.items.len() - 1).sample(bucket_pick) ;
        let pot_item = &bucket.items[pick];
        if Uniform::new_inclusive(0., 2_f64.powi(bucket.bucket_power)).sample(coin_toss)
            < bucket.k_s[pick]
        {
            Some(pot_item.clone())
        } else {
            None
        }
    }

    pub fn choose_ramdom_move_kmc(
        &self,
        rng_choose: &mut SmallRng,
        bucket_pick: &mut SmallRng,
        coin_toss: &mut SmallRng,
        temp: f64,
    ) -> Option<(ItemEnum, f64)> {
        let between = Uniform::new_inclusive(0., self.total_k);
        let mut k_time_rng = between.sample(rng_choose);
        let mut cur_k = 0_f64;

        for (iter, bucket) in self.buckets_list.iter().enumerate() {
            cur_k += bucket.own_k;
            if k_time_rng <= cur_k {
                loop {
                    if let Some(x) = self.pick_from_bucket(bucket_pick, coin_toss, bucket) {
                        return Some((x, self.total_k));
                    }
                }
            }
        }
        return None;
    }

    pub fn iter(&self) -> std::slice::Iter<'_, Bucket> {
        self.buckets_list.iter()
    }
}

