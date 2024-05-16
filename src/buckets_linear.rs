use rand::distributions::{Distribution, Uniform};
use rand::rngs::SmallRng;
use std::collections::HashMap;
use std::{cell::RefCell, rc::Rc};

use super::add_remove_list;
use super::listdict;

#[derive(Clone, Debug)]
enum ItemEnum {
    Move(listdict::Move),
    AddOrRemove(add_remove_list::AtomPosChange),
}

impl ItemEnum {
    fn get_id(&self) -> u64 {
        match self {
            ItemEnum::Move(mmove) => (mmove.from as u64 + ((mmove.to as u64) << 32)),
            ItemEnum::AddOrRemove(add_or_remove) => add_or_remove.pos as u64,
        }
    }
    fn get_k(&self) -> f64 {
        match self {
            ItemEnum::Move(mmove) => mmove.k,
            ItemEnum::AddOrRemove(add_or_remove) => add_or_remove.k,
        }
    }
}

// #[derive(Clone, Debug)]
// struct Item {
//     k: f64,
//     item_type: ItemEnum,
// }

#[derive(Clone, Debug, Default)]
struct Bucket {
    is_parent_to_right: bool,
    noder_number: usize,
    bucket_power: i32,
    own_k: f64,
    // right_total_k: f64,
    items: Vec<ItemEnum>,
}

#[derive(Clone, Debug)]
struct ItemIndexes {
    power_of_k: i32,
    vec_index: usize,
}

#[derive(Clone, Debug)]
struct Buckets {
    bucket_power_to_pos: HashMap<i32, usize, ahash::RandomState>,
    item_to_position: HashMap<u64, ItemIndexes, ahash::RandomState>,
    // buckts: Vec<Node<'a>>,
    total_k: f64,
    buckets_count: u32,
    buckets_list: Vec<Bucket>,
}

impl Buckets {
    fn new() -> Buckets {
        let bucket_power_to_pos: HashMap<i32, usize, ahash::RandomState> = HashMap::default();
        let item_to_position: HashMap<u64, ItemIndexes, ahash::RandomState> = HashMap::default();
        Buckets {
            bucket_power_to_pos,
            item_to_position,
            // buckts: Vec::new(),
            total_k: 0.,
            buckets_count: 0,
            buckets_list: Vec::new(),
        }
    }

    pub fn add_item(&mut self, power_of_k: i32, k: f64, item: ItemEnum) {
        let bucket_index = if let Some(x) = self.bucket_power_to_pos.get(&power_of_k) {
            *x
        } else {
            self.cond_add_bucket(power_of_k)
        };

        self.item_to_position.insert(
            item.get_id(),
            ItemIndexes {
                power_of_k,
                vec_index: self.buckets_list[bucket_index].items.len(),
            },
        );
        self.buckets_list[bucket_index].own_k += k;
        self.buckets_list[bucket_index].items.push(item);
    }

    pub fn cond_add_bucket(&mut self, power_of_k: i32) -> usize {
        if let Some(bucket_index) = self.bucket_power_to_pos.get(&power_of_k) {
            return *bucket_index;
        }
        self.buckets_count += 1;
        println!("b count {}", self.buckets_count);
        println!("adding root_bucket");
        let bucket = Bucket {
            bucket_power: power_of_k,
            ..Default::default()
        };
        self.buckets_list.push(bucket);
        self.bucket_power_to_pos
            .insert(power_of_k, self.buckets_list.len() - 1);
        self.buckets_list.len() - 1
    }

    pub fn remove_item(&mut self, item: ItemEnum) {
        if let Some(item_indexes) = self.item_to_position.remove(&item.get_id()) {
            let bucket_index = self
                .bucket_power_to_pos
                .get(&item_indexes.power_of_k)
                .unwrap();

            let bucket = &mut self.buckets_list[*bucket_index];

            let poped_item = bucket.items.pop().unwrap();
            if bucket.items.len() != item_indexes.vec_index {
                let old_item =
                    std::mem::replace(&mut bucket.items[item_indexes.vec_index], poped_item);
                self.item_to_position.insert(
                    (bucket.items[item_indexes.vec_index].get_id()),
                    item_indexes,
                );
                self.total_k -= old_item.get_k();
                bucket.own_k -= old_item.get_k();
            } else {
                self.total_k -= poped_item.get_k();
                bucket.own_k -= poped_item.get_k();
            }
        }
    }
    pub fn update_k_if_move_exists(&mut self, item: ItemEnum) {
        if let Some(item_indexes) = self.item_to_position.get(&item.get_id()) {
            // let old_energy = std::mem::replace(&mut self.moves[*position].energy, new_energy);

            let bucket = &mut self.buckets_list[item_indexes.vec_index];

            // self.total_k += item.get_k();

            self.total_k -= bucket.items[item_indexes.vec_index].get_k();
            bucket.items[item_indexes.vec_index] = item;
            self.total_k += bucket.items[item_indexes.vec_index].get_k();
        }
    }

    fn pick_from_bucket(
        &self,
        bucket_pick: &mut SmallRng,
        coin_toss: &mut SmallRng,
        bucket: &Bucket,
    ) -> Option<ItemEnum> {
        let pick = Uniform::new_inclusive(0, 1).sample(bucket_pick) * bucket.items.len() - 1;
        let pot_item = &bucket.items[pick];
        if Uniform::new_inclusive(0., 2_f64.powi(bucket.bucket_power)).sample(coin_toss)
            < pot_item.get_k()
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
    ) -> Option<ItemEnum> {
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

        for (iter, bucket) in self.buckets_list.iter().enumerate() {
            cur_k += bucket.own_k;
            if k_time_rng <= cur_k {
                loop {
                    if let Some(x) = self.pick_from_bucket(bucket_pick, coin_toss, bucket) {
                        return Some(x);
                    }
                }
            }
        }
        return None;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // #[test]
    // fn () {
    //
    // }
}
