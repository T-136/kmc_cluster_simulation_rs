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

#[derive(Clone, Debug)]
struct Item {
    k: f64,
    item_type: ItemEnum,
}

#[derive(Clone, Debug, Default)]
struct Bucket {
    parent: Option<BucketRef>,
    noder_number: usize,
    bucket_power: i64,
    own_k: f64,
    left_total_k: f64,
    left: Option<BucketRef>,
    right_total_k: f64,
    right: Option<BucketRef>,
    items: Vec<Item>,
}

#[derive(Clone, Debug)]
struct Buckets {
    root_bucket: Option<BucketRef>,
    bucket_power_to_pos: HashMap<i64, Rc<RefCell<Bucket>>, ahash::RandomState>,
    item_to_position: HashMap<u64, (usize, usize), ahash::RandomState>,
    // buckts: Vec<Node<'a>>,
    total_k: f64,
    buckets_count: u32,
}

impl Buckets {
    fn new() -> Buckets {
        let bucket_power_to_pos: HashMap<i64, Rc<RefCell<Bucket>>, ahash::RandomState> =
            HashMap::default();
        let item_to_position: HashMap<u64, (usize, usize), ahash::RandomState> = HashMap::default();
        Buckets {
            root_bucket: None,
            bucket_power_to_pos,
            item_to_position,
            // buckts: Vec::new(),
            total_k: 0.,
            buckets_count: 0,
        }
    }

    fn insert_bucket(&mut self, current_bucket: Rc<RefCell<Bucket>>, power_of_k: i64) {
        let bucket = Bucket {
            parent: Some(Rc::clone(&current_bucket)),
            bucket_power: power_of_k,
            ..Default::default()
        };
        current_bucket.borrow_mut().left = Some(Rc::new(RefCell::new(bucket)));
        self.bucket_power_to_pos.insert(
            power_of_k,
            Rc::clone(current_bucket.borrow_mut().left.as_ref().unwrap()),
        );
    }

    pub fn add_bucket(&mut self, power_of_k: i64) {
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
                    bucket_power: power_of_k,
                    ..Default::default()
                };
                // let x = &mut *current_bucket.borrow_mut();
                current_bucket.borrow_mut().left = Some(Rc::new(RefCell::new(bucket)));
            } else if current_bucket.borrow().right.is_none() {
                println!("adding right");
                let bucket = Bucket {
                    parent: Some(Rc::clone(&current_bucket)),
                    bucket_power: power_of_k,
                    ..Default::default()
                };
                current_bucket.borrow_mut().right = Some(Rc::new(RefCell::new(bucket)));
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
        }
    }

    // fn find_or_add_bucket(&mut self, power_of_k: i64) {
    //     self.bucket_power_to_pos
    //         .get(&power_of_k)
    //         .get_or_insert("fds");
    // }

    // fn add(&mut self, item: ItemEnum, k: f64) {
    //     let power_of_k: i64 = k.ln().ceil() as i64;
    //
    //     let key: u64 = match item {
    //         ItemEnum::Move(mmove) => (mmove.from as u64 + ((mmove.to as u64) << 32)),
    //         ItemEnum::AddOrRemove(add_or_remove) => add_or_remove.pos as u64,
    //     };
    //
    //     match self.item_to_position.entry(key) {
    //         std::collections::hash_map::Entry::Vacant(e) => {
    //             self.bucket_power_to_pos.get()
    //             self.moves.push();
    //             self.moves_k.push(k);
    //             e.insert(self.moves.len() - 1);
    //             self.total_k += k;
    //         }
    //         _ => return,
    //     }
    // }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_tree() {
        let mut bucket_tree = Buckets::new();
        bucket_tree.add_bucket(19);
        bucket_tree.add_bucket(21);
        bucket_tree.add_bucket(100);
        bucket_tree.add_bucket(110);
        bucket_tree.add_bucket(120);
        bucket_tree.add_bucket(121);
        bucket_tree.add_bucket(221);
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
    }
}
