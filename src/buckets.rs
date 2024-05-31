use std::collections::HashMap;
use std::{cell::RefCell, rc::Rc};

use super::add_remove_list;
use super::listdict;

type BucketRef = Rc<RefCell<Bucket>>;

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
}

// #[derive(Clone, Debug)]
// struct Item {
//     k: f64,
//     item_type: ItemEnum,
// }

#[derive(Clone, Debug, Default)]
struct Bucket {
    parent: Option<BucketRef>,
    is_parent_to_right: bool,
    noder_number: usize,
    bucket_power: i64,
    own_k: f64,
    left_total_k: f64,
    left: Option<BucketRef>,
    // right_total_k: f64,
    right: Option<BucketRef>,
    items: Vec<ItemEnum>,
}

#[derive(Clone, Debug)]
struct ItemIndexes {
    power_of_k: i64,
    vec_index: usize,
}

#[derive(Clone, Debug)]
struct Buckets {
    root_bucket: Option<BucketRef>,
    bucket_power_to_pos: HashMap<i64, Rc<RefCell<Bucket>>, ahash::RandomState>,
    item_to_position: HashMap<u64, ItemIndexes, ahash::RandomState>,
    // buckts: Vec<Node<'a>>,
    total_k: f64,
    buckets_count: u32,
}

impl Buckets {
    fn new() -> Buckets {
        let bucket_power_to_pos: HashMap<i64, Rc<RefCell<Bucket>>, ahash::RandomState> =
            HashMap::default();
        let item_to_position: HashMap<u64, ItemIndexes, ahash::RandomState> = HashMap::default();
        Buckets {
            root_bucket: None,
            bucket_power_to_pos,
            item_to_position,
            // buckts: Vec::new(),
            total_k: 0.,
            buckets_count: 0,
        }
    }

    pub fn add_item(&mut self, power_of_k: i64, k: f64, item: ItemEnum) {
        let ptr_to_bucket = self
            .bucket_power_to_pos
            .get(&power_of_k)
            .unwrap_or_else(|| &self.cond_add_bucket(power_of_k));

        ptr_to_bucket.borrow_mut().own_k += k;
        ptr_to_bucket.borrow_mut().items.push(item);
        if ptr_to_bucket.borrow().is_parent_to_right {
            update_higher_k(k, k, Rc::clone(ptr_to_bucket))
        } else {
            update_higher_k(k, 0., Rc::clone(ptr_to_bucket))
        }
        self.item_to_position.insert(
            item.get_id(),
            ItemIndexes {
                power_of_k,
                vec_index: ptr_to_bucket.borrow_mut().items.len(),
            },
        );
    }

    pub fn cond_add_bucket(&mut self, power_of_k: i64) -> Rc<RefCell<Bucket>> {
        if let Some(ptr) = self.bucket_power_to_pos.get(&power_of_k) {
            return Rc::clone(ptr);
        }
        if self.root_bucket.is_some() {
            let mut current_bucket = Rc::clone(self.root_bucket.as_ref().unwrap());

            self.buckets_count += 1;
            println!("b count {}", self.buckets_count);
            println!("leading 0 {}", self.buckets_count.leading_zeros());
            println!("{:?}", (1..(31 - self.buckets_count.leading_zeros())).rev());
            for x in (1..(31 - self.buckets_count.leading_zeros())).rev() {
                println!("bin {}", (self.buckets_count >> x) & 1);
                if (self.buckets_count >> x) & 1 == 0 {
                    let new_bucket = Rc::clone(current_bucket.borrow().left.as_ref().unwrap());
                    current_bucket = new_bucket;
                } else if (self.buckets_count >> x) & 1 == 1 {
                    let new_bucket = Rc::clone(current_bucket.borrow().right.as_ref().unwrap());
                    current_bucket = new_bucket;
                } else {
                    panic!("wrong binary");
                }
            }
            if current_bucket.borrow().left.is_none() {
                println!("adding left");
                let bucket = Bucket {
                    parent: Some(Rc::clone(&current_bucket)),
                    is_parent_to_right: true,
                    bucket_power: power_of_k,
                    ..Default::default()
                };
                // let x = &mut *current_bucket.borrow_mut();
                current_bucket.borrow_mut().left = Some(Rc::new(RefCell::new(bucket)));
                self.bucket_power_to_pos.insert(
                    power_of_k,
                    Rc::clone(current_bucket.borrow_mut().left.as_ref().unwrap()),
                );
                return Rc::clone(current_bucket.borrow_mut().left.as_ref().unwrap());
            // } else if current_bucket.borrow().right.is_none() {
            } else {
                println!("adding right");
                let bucket = Bucket {
                    parent: Some(Rc::clone(&current_bucket)),
                    is_parent_to_right: false,
                    bucket_power: power_of_k,
                    ..Default::default()
                };
                current_bucket.borrow_mut().right = Some(Rc::new(RefCell::new(bucket)));
                self.bucket_power_to_pos.insert(
                    power_of_k,
                    Rc::clone(current_bucket.borrow_mut().right.as_ref().unwrap()),
                );
                return Rc::clone(current_bucket.borrow_mut().right.as_ref().unwrap());
            }
        } else {
            self.buckets_count += 1;
            println!("b count {}", self.buckets_count);
            println!("adding root_bucket");
            let bucket = Bucket {
                parent: None,
                bucket_power: power_of_k,
                ..Default::default()
            };
            self.root_bucket = Some(Rc::new(RefCell::new(bucket)));
            self.bucket_power_to_pos
                .insert(power_of_k, Rc::clone(self.root_bucket.as_ref().unwrap()));
            return Rc::clone(self.root_bucket.as_ref().unwrap());
        }
    }

    pub fn remove_item(&mut self, k: f64, power_of_k: i64, item: ItemEnum) {
        if let Some(item_indexs) = self.item_to_position.get(&item.get_id()) {
            let remove_k = -k;

            let ptr_to_bucket = self.bucket_power_to_pos.get(&item_indexs.power_of_k);

            ptr_to_bucket.borrow_mut().items.push(item);
            if ptr_to_bucket.borrow().is_parent_to_right {
                update_higher_k(remove_k, remove_k, Rc::clone(ptr_to_bucket))
            } else {
                update_higher_k(remove_k, 0., Rc::clone(ptr_to_bucket))
            }
        }
    }
}

fn update_higher_k(k: f64, left_total_k: f64, ptr: Rc<RefCell<Bucket>>) {
    println!("K {} left {}", k, left_total_k);
    ptr.borrow_mut().left_total_k += left_total_k;
    // ptr.borrow_mut().own_k += k;
    if let Some(parent) = ptr.borrow().parent.as_ref() {
        if ptr.borrow().is_parent_to_right {
            update_higher_k(k, k, Rc::clone(parent));
        } else {
            update_higher_k(k, 0., Rc::clone(parent));
        }
    }
    // else {
    //     return None;
    // }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_tree() {
        let mut bucket_tree = Buckets::new();
        bucket_tree.cond_add_bucket(19);
        bucket_tree.cond_add_bucket(21);
        bucket_tree.cond_add_bucket(100);
        bucket_tree.cond_add_bucket(110);
        bucket_tree.cond_add_bucket(120);
        bucket_tree.cond_add_bucket(121);
        bucket_tree.cond_add_bucket(221);
        // bucket_tree.insert_bucket(211);
        // bucket_tree.insert_bucket(212);

        let first_bucket = Rc::clone(
            bucket_tree
                .root_bucket
                .unwrap()
                .borrow()
                .right
                .as_ref()
                .unwrap(),
        );

        let six_bucket = first_bucket
            .borrow()
            .left
            .as_ref()
            .unwrap()
            .borrow()
            .bucket_power
            .clone();

        let seven_bucket = first_bucket
            .borrow()
            .right
            .as_ref()
            .unwrap()
            .borrow()
            .bucket_power
            .clone();
        println!("test {:?}", seven_bucket);
        // println!("test {:?}", right_bucket.borrow().left);
        assert_eq!(seven_bucket, 221);
        assert_eq!(six_bucket, 121);
        assert_eq!(bucket_tree.bucket_power_to_pos.len(), 7);
    }
    #[test]
    fn test_insert_item() {
        let mut bucket_tree = Buckets::new();
        let item = ItemEnum::Move(listdict::Move {
            from: 1,
            to: 0,
            e_diff: 13.,
            e_barr: 18.,
        });
        let item2 = ItemEnum::Move(listdict::Move {
            from: 1,
            to: 0,
            e_diff: 13.,
            e_barr: 18.,
        });
        let item3 = ItemEnum::Move(listdict::Move {
            from: 1,
            to: 0,
            e_diff: 13.,
            e_barr: 18.,
        });
        bucket_tree.add_item(10, 0.3, item);
        bucket_tree.add_item(14, 1., item2.clone());
        bucket_tree.add_item(13, 1., item2.clone());
        bucket_tree.add_item(12, 1., item2.clone());
        bucket_tree.add_item(13, 0.3, item2.clone());
        bucket_tree.add_item(11, 3., item3);

        let first_bucket_power = bucket_tree
            .root_bucket
            .as_ref()
            .unwrap()
            .borrow()
            .bucket_power
            .clone();

        let first_bucket_left_k = bucket_tree
            .root_bucket
            .as_ref()
            .unwrap()
            .borrow()
            .left_total_k
            .clone();

        assert_eq!(first_bucket_power, 10);
        assert_eq!(first_bucket_left_k, 5.0);
        assert_eq!(bucket_tree.bucket_power_to_pos.len(), 5);
    }
}
