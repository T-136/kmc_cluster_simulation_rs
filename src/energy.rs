#[derive(Clone, Debug)]
pub enum EnergyInput {
    LinearCn(Vec<[i64; 2]>),
    Cn(Vec<[i64; 13]>),
    LinearGcn(Vec<[i64; 2]>),
    Gcn(Vec<[i64; 145]>),
}

pub fn energy_1000_calculation(
    energy: &EnergyInput,
    cn: usize,
    mut atom_typ_index: usize,
    at_support: u8,
    support_e: i64,
) -> i64 {
    atom_typ_index -= 1;
    match energy {
        EnergyInput::LinearCn(e) => {
            e[atom_typ_index][0] * cn as i64 + e[atom_typ_index][1] + at_support as i64 * support_e
        }
        EnergyInput::Cn(e) => e[atom_typ_index][cn] + at_support as i64 * support_e,
        EnergyInput::LinearGcn(e) => e[atom_typ_index][0] * cn as i64 + e[atom_typ_index][1],
        EnergyInput::Gcn(e) => e[atom_typ_index][cn],
    }
    // enrico_table(cn[*atom as usize])
    // cn[*atom as usize] as i64 * M_BETA + M_ALPHA
}

pub fn energy_diff_cn<I, O>(
    energy: &[[i64; 13]],
    cn_from_list: I,
    cn_to_list: O,
    move_from_cn: usize,
    move_to_cn: usize,
    atom_typ: usize,
    is_from_at_support: u8,
    is_to_at_support: u8,
    support_e: i64,
) -> i64
where
    I: Iterator<Item = (usize, u8)>,
    O: Iterator<Item = (usize, u8)>,
{
    // let energy = energy[atom_typ_index];
    let mut energy_diff_1000 = 0;
    for (cn_from, atom_typ) in cn_from_list {
        energy_diff_1000 -= energy[atom_typ as usize][cn_from];
        energy_diff_1000 += energy[atom_typ as usize][cn_from - 1];
    }
    for (cn_to, atom_typ) in cn_to_list {
        energy_diff_1000 -= energy[atom_typ as usize][cn_to];
        energy_diff_1000 += energy[atom_typ as usize][cn_to + 1];
    }
    energy_diff_1000 -= energy[atom_typ][move_from_cn];
    if is_from_at_support == 1 {
        energy_diff_1000 -= support_e;
    }
    energy_diff_1000 += energy[atom_typ][move_to_cn - 1];
    if is_to_at_support == 1 {
        energy_diff_1000 += support_e;
    }

    energy_diff_1000
}
pub fn energy_diff_gcn<I, O, P>(
    energy: &Vec<[i64; 145]>,
    cn_from_list: I,
    cn_to_list: O,
    intersect_change: P,
    move_from_gcn: usize,
    move_to_gcn: usize,
    atom_typ_index: usize,
) -> i64
where
    I: Iterator<Item = (usize, usize)>,
    O: Iterator<Item = (usize, usize)>,
    P: Iterator<Item = (usize, usize)>,
{
    let energy = energy[atom_typ_index];
    let mut energy_diff_1000 = 0;
    for (gcn_from_old, gcn_from_new) in cn_from_list {
        energy_diff_1000 -= energy[gcn_from_old];
        energy_diff_1000 += energy[gcn_from_new];
    }
    for (gcn_to_old, gcn_to_new) in cn_to_list {
        energy_diff_1000 -= energy[gcn_to_old];
        energy_diff_1000 += energy[gcn_to_new];
    }
    for (gcn_inter_old, gcn_inter_new) in intersect_change {
        energy_diff_1000 -= energy[gcn_inter_old];
        energy_diff_1000 += energy[gcn_inter_new];
    }
    energy_diff_1000 += energy[move_to_gcn];
    energy_diff_1000 -= energy[move_from_gcn];

    energy_diff_1000
}

pub fn energy_diff_l_cn(
    energy: &[[i64; 2]],
    cn_from: usize,
    cn_to: usize,
    neigbors_of_from: &[u8; super::NUM_ATOM_TYPES],
    neigbors_of_to: &[u8; super::NUM_ATOM_TYPES],
    atom_typ_main: usize,
    is_from_at_support: u8,
    is_to_at_support: u8,
    support_e: i64,
) -> i64 {
    let mut energy_change = 0;
    energy_change -= (cn_from as i64) * energy[atom_typ_main][0] + energy[atom_typ_main][1];
    energy_change += (cn_to as i64 - 1) * energy[atom_typ_main][0] + energy[atom_typ_main][1];
    for (atom_typ_index, neigh) in neigbors_of_from.iter().enumerate() {
        // atom_typ_index -= 1;
        energy_change -= (*neigh as i64) * energy[atom_typ_index][0];
    }
    for (atom_typ_index, neigh) in neigbors_of_to.iter().enumerate() {
        // atom_typ_index -= 1;
        if atom_typ_index == atom_typ_main {
            energy_change += (*neigh as i64 - 1) * energy[atom_typ_index][0];
        } else {
            energy_change += (*neigh as i64) * energy[atom_typ_index][0];
        }
    }
    if is_from_at_support == 1 {
        energy_change -= support_e;
    }
    if is_to_at_support == 1 {
        energy_change += support_e;
    }
    energy_change
}

pub fn energy_diff_l_gcn<I, O>(
    energy: &Vec<[i64; 2]>,
    cn_from_list: I,
    cn_to_list: O,
    move_from_cn: usize,
    move_to_cn: usize,
    atom_typ_index: usize,
) -> i64
where
    I: Iterator<Item = usize>,
    O: Iterator<Item = usize>,
{
    let energy = energy[atom_typ_index];
    let mut energy_diff_1000 = 0;
    for cn_from in cn_from_list {
        energy_diff_1000 -= energy[0] * cn_from as i64 + energy[1];
        energy_diff_1000 += energy[0] * (cn_from as i64 - 1) + energy[1];
    }
    for cn_to in cn_to_list {
        energy_diff_1000 -= energy[0] * cn_to as i64 + energy[1];
        energy_diff_1000 += energy[0] * (cn_to as i64 + 1) + energy[1];
    }
    // energy_diff_1000 -= energy[move_from_cn];
    // energy_diff_1000 += energy[move_to_cn - 1];

    2 * energy_diff_1000
}
