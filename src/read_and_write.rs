use chemfiles::{Atom, Frame, Trajectory, UnitCell};
use fnv::FnvBuildHasher;
use fnv::FnvHashMap;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs;
use std::io::{self, BufRead};
use time_graph::instrument;

fn find_key_for_value<'a>(map: &'a HashMap<String, u8>, value: u8) -> Option<&'a str> {
    map.iter().find_map(|(key, val)| {
        if val == &value {
            Some(key.as_str())
        } else {
            None
        }
    })
}

pub fn write_occ_as_xyz(
    // trajectory: &mut Trajectory,
    save_folder: String,
    onlyocc: HashSet<u32, fnv::FnvBuildHasher>,
    xsites_positions: &Vec<[f64; 3]>,
    unit_cell: &[f64; 3],
    atom_pos: &Vec<super::AtomPosition>,
    atom_names: &HashMap<String, u8>,
) {
    let mut trajectory = Trajectory::open(save_folder.clone() + "/lowest_energy.xyz", 'w').unwrap();
    let mut xyz: Vec<[f64; 3]> = Vec::new();
    for (j, ii) in onlyocc.iter().enumerate() {
        xyz.insert(j, xsites_positions[*ii as usize]);
    }
    let mut frame = Frame::new();
    frame.set_cell(&UnitCell::new(unit_cell.clone()));

    for (i, atom) in atom_pos.iter().enumerate() {
        if atom.occ == 255 {
            continue;
        }
        frame.add_atom(
            &Atom::new(
                find_key_for_value(atom_names, atom.occ).expect(
                    format!("unknown atom number {:?}, {:?}", atom_names, atom.occ).as_str(),
                ),
            ),
            xsites_positions[i],
            None,
        );
    }

    trajectory
        .write(&frame)
        .unwrap_or_else(|x| eprintln!("{}", x));
}

#[instrument]
pub fn read_sample(input_file: &str) -> Vec<(String, [f64; 3])> {
    // if input_file.contains(".poscar") {
    //     let newatoms = Poscar::from_path(input_file).unwrap();
    //     let xyz = newatoms.scaled_cart_positions();
    //     xyz.into_owned()
    // } else if input_file.contains(".xyz") {
    let mut trajectory = Trajectory::open(input_file, 'r').unwrap();
    let mut frame = Frame::new();
    trajectory.read(&mut frame).unwrap();
    let mut atom_vec: Vec<(String, [f64; 3])> = Vec::new();
    for (i, atom) in frame.iter_atoms().enumerate() {
        atom_vec.push((atom.name(), frame.positions()[i]));
    }
    atom_vec
    // } else {
    // panic!("no .poscar or .xyz, cant read file");
    // }
}

fn fmt_scient(num: &str) -> f64 {
    let mut parts = num.split('e');

    let pre_num = parts.next().unwrap();
    let exp = parts.next().unwrap();

    let base: f64 = 10.;
    pre_num.parse::<f64>().unwrap() * base.powi(exp.parse::<i32>().unwrap())
}

pub fn read_atom_sites(input_file: &str, nsites: u32) -> Vec<[f64; 3]> {
    // println!("reading atom_sites from: {}", input_file);
    let mut xsites_positions: Vec<[f64; 3]> = Vec::with_capacity(nsites as usize);
    let pairlist = fs::File::open(input_file).expect("Should have been able to read the file");
    let lines = io::BufReader::new(pairlist);

    for line in lines.lines() {
        let r = line.unwrap();
        let list: Vec<&str> = r.split_whitespace().clone().collect();
        let temp_str_vec: [&str; 3] = [list[0], list[1], list[2]];
        let temp_vec: [f64; 3] = temp_str_vec.map(|i| fmt_scient(i));
        xsites_positions.push(temp_vec);
    }
    xsites_positions
}

pub fn read_nn(pairlist_file: &str) -> HashMap<u32, [u32; super::CN], FnvBuildHasher> {
    // println!("reading pairlists from: {}", pairlist_file);

    let pairlist = fs::File::open(pairlist_file).expect("Should have been able to read the file");

    let lines = io::BufReader::new(pairlist);
    let mut nn: HashMap<u32, [u32; super::CN], FnvBuildHasher> =
        FnvHashMap::with_capacity_and_hasher(5400, Default::default());

    for line in lines.lines() {
        let r = line.unwrap();
        let list: Vec<&str> = r.split_whitespace().clone().collect();
        let mut neighbors: [u32; super::CN] = [0; super::CN];
        let prime = list.first().clone();
        for (i, l) in list.iter().skip(1).enumerate() {
            neighbors[i] = l.parse::<u32>().unwrap()
        }
        nn.insert(prime.unwrap().parse::<u32>().unwrap(), neighbors);
    }
    nn
}

pub fn read_nnn(pairlist_file: &str) -> HashMap<u32, [u32; super::GCN], FnvBuildHasher> {
    // println!("reading pairlists from: {}", pairlist_file);

    let pairlist = fs::File::open(pairlist_file).expect("Should have been able to read the file");

    let lines = io::BufReader::new(pairlist);
    let mut nnn: HashMap<u32, [u32; super::GCN], FnvBuildHasher> =
        FnvHashMap::with_capacity_and_hasher(5400, Default::default());

    for line in lines.lines() {
        let r = line.unwrap();
        let list: Vec<&str> = r.split_whitespace().clone().collect();
        let mut neighbors: [u32; super::GCN] = [0; super::GCN];
        let prime = list.first().clone();
        for (i, l) in list.iter().skip(1).enumerate() {
            neighbors[i] = l.parse::<u32>().unwrap()
        }
        nnn.insert(prime.unwrap().parse::<u32>().unwrap(), neighbors);
    }
    nnn
}

pub fn read_surounding_moves(
    surroundin_moves_file: &str,
) -> HashMap<u64, Vec<(u32, u32)>, fnv::FnvBuildHasher> {
    let nnn_pairlist =
        fs::File::open(surroundin_moves_file).expect("Should have been able to read the file");

    let reader = io::BufReader::new(nnn_pairlist);

    let surrounding_moves_no_bit_shifting: HashMap<
        u32,
        HashMap<u32, Vec<(u32, u32)>, fnv::FnvBuildHasher>,
    > = serde_json::from_reader(reader).unwrap();

    let mut surrounding_moves: HashMap<u64, Vec<(u32, u32)>, fnv::FnvBuildHasher> =
        FnvHashMap::with_capacity_and_hasher(5400, Default::default());
    for (k1, d1) in surrounding_moves_no_bit_shifting.into_iter() {
        for (k2, d2) in d1.into_iter() {
            surrounding_moves.insert((k1 as u64 + ((k2 as u64) << 32)), d2);
        }
    }

    // let nnn_pair: HashMap<u64, [HashMap<u32, Vec<u32>, FnvBuildHasher>; 2], FnvBuildHasher> =
    //     serde_json::from_reader(reader).unwrap();

    surrounding_moves
}

pub fn read_nnn_pair_no_intersec(
    nnn_pairlist_file: &str,
) -> HashMap<
    u64,
    (
        Vec<Vec<u32>>,
        Vec<Vec<u32>>,
        Vec<(u32, Vec<u32>, Vec<u32>, Vec<u32>)>,
    ),
    fnv::FnvBuildHasher,
> {
    let nnn_pairlist =
        fs::File::open(nnn_pairlist_file).expect("Should have been able to read the file");

    let reader = io::BufReader::new(nnn_pairlist);

    let nnn_pair_no_bit_shifting: HashMap<
        u32,
        HashMap<
            u32,
            (
                Vec<Vec<u32>>,
                Vec<Vec<u32>>,
                Vec<(u32, Vec<u32>, Vec<u32>, Vec<u32>)>,
            ),
            fnv::FnvBuildHasher,
        >,
    > = serde_json::from_reader(reader).unwrap();

    let mut nnn_pair: HashMap<
        u64,
        (
            Vec<Vec<u32>>,
            Vec<Vec<u32>>,
            Vec<(u32, Vec<u32>, Vec<u32>, Vec<u32>)>,
        ),
        fnv::FnvBuildHasher,
    > = FnvHashMap::with_capacity_and_hasher(5400, Default::default());
    for (k1, d1) in nnn_pair_no_bit_shifting.into_iter() {
        for (k2, d2) in d1.into_iter() {
            nnn_pair.insert((k1 as u64 + ((k2 as u64) << 32)), d2);
        }
    }

    // let nnn_pair: HashMap<u64, [HashMap<u32, Vec<u32>, FnvBuildHasher>; 2], FnvBuildHasher> =
    //     serde_json::from_reader(reader).unwrap();

    return nnn_pair;
}

pub fn read_nn_pair_no_intersec(
    nn_pairlist_file: &str,
) -> HashMap<
    u64,
    (
        [u32; super::NN_PAIR_NO_INTERSEC_NUMBER],
        [u32; super::NN_PAIR_NO_INTERSEC_NUMBER],
        [u32; super::NN_PAIR_ONLY_INTERSEC_NUMBER],
    ),
    FnvBuildHasher,
> {
    let nn_pairlist =
        fs::File::open(nn_pairlist_file).expect("Should have been able to read the file");

    let lines = io::BufReader::new(nn_pairlist);

    let mut nn_pair: HashMap<
        u64,
        (
            [u32; super::NN_PAIR_NO_INTERSEC_NUMBER],
            [u32; super::NN_PAIR_NO_INTERSEC_NUMBER],
            [u32; super::NN_PAIR_ONLY_INTERSEC_NUMBER],
        ),
        FnvBuildHasher,
    > = FnvHashMap::with_capacity_and_hasher(32000, Default::default());

    for line in lines.lines() {
        let r = line.unwrap();
        let test: Vec<&str> = r.split_whitespace().clone().collect();
        let site: u32 = std::cmp::min(
            test[0].parse::<u32>().unwrap(),
            test[1].parse::<u32>().unwrap(),
        );
        let j: u32 = std::cmp::max(
            test[0].parse::<u32>().unwrap(),
            test[1].parse::<u32>().unwrap(),
        );
        let mut neighbors = (
            [0_u32; super::NN_PAIR_NO_INTERSEC_NUMBER],
            [0_u32; super::NN_PAIR_NO_INTERSEC_NUMBER],
            [0_u32; super::NN_PAIR_ONLY_INTERSEC_NUMBER],
        );

        for (i, l) in test.iter().skip(2).enumerate() {
            if i < 7 {
                neighbors.0[i] = l.parse::<u32>().unwrap()
            } else if i < 14 {
                neighbors.1[i - 7] = l.parse::<u32>().unwrap()
            } else {
                neighbors.2[i - 14] = l.parse::<u32>().unwrap()
            }
        }
        nn_pair
            // .entry()
            // .and_modify(|map| {
            //     map.insert(j, neighbors.clone());
            // })
            .insert(site as u64 + ((j as u64) << 32), neighbors);
        // println!("{:?}", line.unwrap());
    }

    return nn_pair;
}

pub fn read_nn_pairlists(
    nn_pairlist_file: &str,
) -> HashMap<u64, [u32; super::NN_PAIR_NUMBER], FnvBuildHasher> {
    let nn_pairlist =
        fs::File::open(nn_pairlist_file).expect("Should have been able to read the file");

    let lines = io::BufReader::new(nn_pairlist);

    let mut nn_pair: HashMap<u64, [u32; super::NN_PAIR_NUMBER], FnvBuildHasher> =
        FnvHashMap::with_capacity_and_hasher(32000, Default::default());

    for line in lines.lines() {
        let r = line.unwrap();
        let test: Vec<&str> = r.split_whitespace().clone().collect();
        let site: u32 = std::cmp::min(
            test[0].parse::<u32>().unwrap(),
            test[1].parse::<u32>().unwrap(),
        );
        let j: u32 = std::cmp::max(
            test[0].parse::<u32>().unwrap(),
            test[1].parse::<u32>().unwrap(),
        );
        let mut neighbors: [u32; super::NN_PAIR_NUMBER] = [0; super::NN_PAIR_NUMBER];

        for (i, l) in test.iter().skip(2).enumerate() {
            neighbors[i] = l.parse::<u32>().unwrap()
        }
        nn_pair
            // .entry()
            // .and_modify(|map| {
            //     map.insert(j, neighbors.clone());
            // })
            .insert(site as u64 + ((j as u64) << 32), neighbors);
        // println!("{:?}", line.unwrap());
    }

    return nn_pair;
}

pub fn read_nnn_pairlists(nn_pairlist_file: &str) -> HashMap<u64, [u32; 74], FnvBuildHasher> {
    let nn_pairlist =
        fs::File::open(nn_pairlist_file).expect("Should have been able to read the file");

    let lines = io::BufReader::new(nn_pairlist);

    let mut nn_pair: HashMap<u64, [u32; 74], FnvBuildHasher> =
        FnvHashMap::with_capacity_and_hasher(32000, Default::default());

    for line in lines.lines() {
        let r = line.unwrap();
        let test: Vec<&str> = r.split_whitespace().clone().collect();
        let site: u32 = std::cmp::min(
            test[0].parse::<u32>().unwrap(),
            test[1].parse::<u32>().unwrap(),
        );
        let j: u32 = std::cmp::max(
            test[0].parse::<u32>().unwrap(),
            test[1].parse::<u32>().unwrap(),
        );
        let mut neighbors: [u32; 74] = [0; 74];

        for (i, l) in test.iter().skip(2).enumerate() {
            neighbors[i] = l.parse::<u32>().unwrap()
        }
        nn_pair
            // .entry()
            // .and_modify(|map| {
            //     map.insert(j, neighbors.clone());
            // })
            .insert(site as u64 + ((j as u64) << 32), neighbors);
        // println!("{:?}", line.unwrap());
    }

    return nn_pair;
}
