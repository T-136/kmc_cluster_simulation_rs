use core::panic;
use fnv::FnvBuildHasher;
use std::collections::{HashMap, HashSet};

fn create_support(
    atom_pos: &mut Vec<super::AtomPosition>,
    xsites_positions: &Vec<[f64; 3]>,
    support_indices: Vec<u32>,
    nn: &HashMap<u32, [u32; 12], FnvBuildHasher>,
    iclose: u32,
) -> Vec<u8> {
    let center_of_mass: &[f64; 3] = &xsites_positions[iclose as usize];
    let mut support = Vec::new();
    let refpos = xsites_positions[0];
    for (i, xyz) in xsites_positions.iter().enumerate() {
        let mut norm = 0_f64;
        for (i2, x) in xyz.iter().enumerate() {
            norm += ((*x - center_of_mass[i2]) * support_indices[i2] as f64);
        }
        if norm.abs() < 1e-7 {
            atom_pos[i].occ = 100;
            support.push(i as u32)
        };
    }
    let mut second_layer_fixpoint = 0;
    for neighbor in nn[&iclose] {
        if atom_pos[neighbor as usize].occ == 0 {
            second_layer_fixpoint = neighbor;
        }
    }
    for (i, xyz) in xsites_positions.iter().enumerate() {
        let mut norm = 0_f64;
        for (i2, x) in xyz.iter().enumerate() {
            norm += ((*x - xsites_positions[second_layer_fixpoint as usize][i2])
                * support_indices[i2] as f64);
        }
        if norm.abs() < 1e-7 {
            atom_pos[i].occ = 100;
            support.push(i as u32)
        };
    }
    let mut nn_support = vec![0_u8; xsites_positions.len()];
    nn_support_from_supprt(atom_pos, nn);
    // panic!("fds");
    nn_support
}

pub fn nn_support_from_supprt(
    a_positions: &mut Vec<super::AtomPosition>,
    // nn_support: &mut Vec<u8>,
    nn: &HashMap<u32, [u32; 12], FnvBuildHasher>,
    // occ: &Vec<u8>,
) {
    for atom_pos in 0..a_positions.len() {
        if a_positions[atom_pos].occ != 100 {
            continue;
        }
        for neighbor in nn[&(atom_pos as u32)] {
            if a_positions[neighbor as usize].occ != 2 {
                a_positions[neighbor as usize].nn_support = 1;
            }
        }
    }
}

pub fn create_input_cluster(
    atom_pos: &mut Vec<super::AtomPosition>,
    number_of_atoms: u32,
    xsites_positions: &Vec<[f64; 3]>,
    nn: &HashMap<u32, [u32; 12], FnvBuildHasher>,
    // nsites: u32,
    support_indices: Option<Vec<u32>>,
    coating: Option<String>,
    atom_names: &HashMap<String, u8>,
) -> HashSet<u32, FnvBuildHasher> {
    let center_of_mass: [f64; 3] = {
        let mut d: [Vec<f64>; 3] = [Vec::new(), Vec::new(), Vec::new()];
        for coord in xsites_positions {
            d[0].push(coord[0]);
            d[1].push(coord[1]);
            d[2].push(coord[2]);
        }
        [
            d[0].iter().sum::<f64>() / d[0].len() as f64,
            d[1].iter().sum::<f64>() / d[1].len() as f64,
            d[2].iter().sum::<f64>() / d[2].len() as f64,
        ]
    };
    let iclose: u32 = {
        let mut minimum: f64 = f64::INFINITY;
        let mut index_of_center: u32 = 0;
        // let mut v = Vec::new();
        for (i, xsite_position) in xsites_positions.iter().enumerate() {
            let dist = (center_of_mass[0] - xsite_position[0]).powf(2.)
                + (center_of_mass[1] - xsite_position[1]).powf(2.)
                + (center_of_mass[2] - xsite_position[2]).powf(2.);
            if minimum >= dist {
                minimum = dist.clone();
                index_of_center = i as u32;
            }
        }
        index_of_center
    };
    // println!(
    //     "{:?}, {:?}, {:?}",
    //     center_of_mass, iclose, xsites_positions[iclose as usize]
    // );
    // println!("nsites: {}", atom_pos.len());
    assert_eq!(xsites_positions.len(), atom_pos.len());
    // let mut occ: Vec<u8> = Vec::with_capacity(xsites_positions.len());
    // occ = vec![0; nsites as usize];
    let mut onlyocc: HashSet<u32, FnvBuildHasher> =
        fnv::FnvHashSet::with_capacity_and_hasher(number_of_atoms as usize, Default::default());

    let nn_support =
        support_indices.map(|sup| create_support(atom_pos, xsites_positions, sup, nn, iclose));

    // occ.entry(&iclose).and_modify(1) = 1;
    if atom_pos[iclose as usize].occ == 0 {
        onlyocc.insert(iclose);
    } else {
        for neighbor in nn[&iclose] {
            if atom_pos[neighbor as usize].occ == 0 {
                onlyocc.insert(neighbor);
                break;
            }
        }
    }

    let mut ks = 1;
    assert!((number_of_atoms as usize) < atom_pos.len());

    loop {
        let mut onlyocc_temp_storag: HashSet<u32> = HashSet::new();
        for site in onlyocc.iter() {
            for j in &nn[site] {
                if !onlyocc_temp_storag.contains(&j)
                    && !onlyocc.contains(&j)
                    && atom_pos[*j as usize].occ == 0
                {
                    onlyocc_temp_storag.insert(*j);
                    // onlyocc.insert(j);
                    ks += 1;
                    if ks == number_of_atoms {
                        break;
                    }
                }
            }
            if ks == number_of_atoms {
                break;
            }
        }
        for v in onlyocc_temp_storag {
            onlyocc.insert(v);
        }
        if ks == number_of_atoms {
            break;
        }
    }
    for site in onlyocc.iter() {
        atom_pos[*site as usize].occ = atom_names[&"Pt".to_string()];
    }
    if let Some(coat_atom) = coating {
        let mut all_neigbors = Vec::new();
        for atom in onlyocc.iter() {
            all_neigbors.extend_from_slice(&nn[atom]);
        }
        for neighbor in all_neigbors {
            if atom_pos[neighbor as usize].occ == 0 {
                atom_pos[neighbor as usize].occ = atom_names[&coat_atom];
                onlyocc.insert(neighbor);
            }
        }
    }

    onlyocc
}

pub fn occ_onlyocc_from_xyz(
    atom_pos: &mut Vec<super::AtomPosition>,
    xyz: &Vec<(String, [f64; 3])>,
    xsites_positions: &[[f64; 3]],
    atom_names: &HashMap<String, u8>,
    coating: Option<String>,
    nn: &HashMap<u32, [u32; 12], FnvBuildHasher>,
) -> HashSet<u32, FnvBuildHasher> {
    let supp_metal: String = "Al".to_string();
    // let mut occ: Vec<u8> = vec![0; nsites as usize];
    let mut onlyocc: HashSet<u32, FnvBuildHasher> =
        fnv::FnvHashSet::with_capacity_and_hasher(xyz.len(), Default::default());

    for x in xyz.iter() {
        for site in 0..atom_pos.len() {
            let dist = (x.1[0] - xsites_positions[site as usize][0]).powf(2.)
                + (x.1[1] - xsites_positions[site as usize][1]).powf(2.)
                + (x.1[2] - xsites_positions[site as usize][2]).powf(2.);
            if dist < 2.15 {
                if x.0 == supp_metal {
                    atom_pos[site as usize].occ = 100;
                    continue;
                }
                atom_pos[site as usize].occ = atom_names[&x.0];
                onlyocc.insert(site as u32);
            }
        }
    }
    println!("onlyocc len: {}", onlyocc.len());
    if let Some(coat_atom) = coating {
        let mut all_neigbors = Vec::new();
        for atom in onlyocc.iter() {
            all_neigbors.extend_from_slice(&nn[atom]);
        }
        for neighbor in all_neigbors {
            if atom_pos[neighbor as usize].occ == 0 {
                atom_pos[neighbor as usize].occ = atom_names[&coat_atom];
                onlyocc.insert(neighbor);
            }
        }
    }
    println!("onlyocc len: {}", onlyocc.len());
    onlyocc
}
