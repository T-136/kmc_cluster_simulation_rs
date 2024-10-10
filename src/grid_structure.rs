use super::read_and_write;
use std::collections::HashMap;
use vasp_poscar::Poscar;

pub const NN_PAIR_NO_INTERSEC_NUMBER: usize = 7;
pub const NN_PAIR_ONLY_INTERSEC_NUMBER: usize = 4;

pub struct GridStructure {
    pub nn: HashMap<u32, [u32; super::CN], fnv::FnvBuildHasher>,
    pub nn_pair_no_intersec: HashMap<
        u64,
        (
            [u32; NN_PAIR_NO_INTERSEC_NUMBER],
            [u32; NN_PAIR_NO_INTERSEC_NUMBER],
            [u32; NN_PAIR_ONLY_INTERSEC_NUMBER],
        ),
        fnv::FnvBuildHasher,
    >,
    pub xsites_positions: Vec<[f64; 3]>,
    pub unit_cell: [f64; 3],
    pub surrounding_moves: HashMap<u64, Vec<(u32, u32)>, fnv::FnvBuildHasher>,
}

impl GridStructure {
    pub fn new(
        nn_file: String,
        nn_pair_no_int_file: String,
        atom_sites: String,
        surrounding_moves_file: String,
        grid_file: String,
    ) -> GridStructure {
        let (unit_cell, nsites) = read_and_write::unitcell_from_grid(&grid_file);
        let nn = read_and_write::read_nn(&nn_file);
        let nn_pair_no_intersec = read_and_write::read_nn_pair_no_intersec(&nn_pair_no_int_file);
        let surrounding_moves = read_and_write::read_surounding_moves(&surrounding_moves_file);

        let xsites_positions = read_and_write::read_atom_sites(&atom_sites, nsites);

        GridStructure {
            nn,
            nn_pair_no_intersec,
            surrounding_moves,
            xsites_positions,
            unit_cell,
        }
    }
}
