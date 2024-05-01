// use ahash::HashMap;
use clap::ArgGroup;
use clap::Parser;
use core::panic;
use csv::Reader;
use csv::ReaderBuilder;
use mc::alpha_energy;
use mc::energy::EnergyInput;
use mc::GridStructure;
use mc::Simulation;
use std::collections::HashMap;
use std::fs;
use std::io::BufRead;
use std::io::BufReader;
// use std::io::{self, BufRead};
use std::path::Path;
use std::sync::Arc;
use std::thread;
use std::usize;

use std::convert::TryInto;

fn vec_to_array<T, const N: usize>(v: Vec<T>) -> [T; N] {
    v.try_into()
        .unwrap_or_else(|v: Vec<T>| panic!("Expected a Vec of length {} but it was {}", N, v.len()))
}

fn fmt_scient(num: &str) -> u64 {
    let mut parts = num.split(['e', 'E']);

    let pre_num = parts.next().unwrap();
    let exp = parts.next().unwrap_or("0");
    if parts.next().is_some() {
        panic!("wrong iterations input");
    }

    let base: u64 = 10;
    pre_num.parse::<u64>().expect("wrong iterations input")
        * base.pow(exp.parse::<u32>().expect("wrong iterations input"))
}

fn collect_energy_values<const N: usize>(
    mut energy_vec: [i64; N],
    inp: String,
) -> (Vec<[i64; N]>, HashMap<String, u8>) {
    if inp.chars().next().unwrap().is_numeric() || inp.starts_with('-') {
        let mut string_iter = inp.trim().split(',');
        for x in energy_vec.iter_mut() {
            *x = string_iter
                .next()
                .unwrap()
                .trim()
                .parse::<i64>()
                .unwrap_or_else(|err| {
                    panic!(
                        "iter received from input file: {:?}, err: {}",
                        string_iter, err
                    )
                });
        }
        let mut map: HashMap<String, u8> = HashMap::new();
        map.insert("Pt".to_string(), 1_u8);
        (vec![energy_vec], map)
    } else {
        // let s = fs::read_to_string(inp.clone()).expect("can't find energy file");
        let file = fs::File::open(inp).expect("can't find energy file");
        let reader = BufReader::new(file);
        // let res: Result<Results, serde_json::Error> = serde_json::from_reader(reader);
        let res: Result<HashMap<String, Vec<i64>, fnv::FnvBuildHasher>, serde_json::Error> =
            serde_json::from_reader(reader);
        println!("{:?}", res);

        let mut energy_vec = Vec::new();

        let mut atom_names: HashMap<String, u8> = HashMap::new();
        // println!("{:?}", record);

        for (i, (atom_name, energy_line)) in res.unwrap().into_iter().enumerate() {
            println!("{:?}", energy_line);

            atom_names.insert(atom_name, i as u8 + 1);
            energy_vec.push(vec_to_array(energy_line));
        }

        (energy_vec, atom_names)
    }
}

#[derive(Parser, Debug)]
#[clap(group(
        ArgGroup::new("startstructure")
            .required(true)
            .args(&["start_cluster", "atoms_input"]),
    ))]
struct StartStructure {
    #[arg(short, long)]
    start_cluster: Option<String>,

    #[arg(short, long, value_delimiter = ',')]
    atoms_input: Option<Vec<u32>>,
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
// #[clap(group(
//         ArgGroup::new("energy")
//             .multiple(true)
//             .required(true)
//             .args(&["e_l_cn", "e_cn", "e_l_gcn", "e_gcn", "alphas"]),
//     ))]
struct Args {
    #[clap(flatten)]
    start_structure: StartStructure,

    /// folder to store results
    #[arg(short, long, default_value_t = String::from("./sim/"))]
    folder: String,

    /// iterations
    #[arg(short, long)]
    iterations: String,

    #[arg(short, long, default_value_t = 300.)]
    temperature: f64,

    // #[arg(long, allow_hyphen_values(true))]
    // e_l_cn: Option<String>,
    //
    // #[arg(long, allow_hyphen_values(true))]
    // e_cn: Option<String>,
    //
    // #[arg(long, allow_hyphen_values(true))]
    // e_l_gcn: Option<String>,
    //
    // #[arg(long, allow_hyphen_values(true))]
    // e_gcn: Option<String>,
    #[arg(long)]
    alphas: String,

    #[arg(short, long, value_delimiter = '-', default_values_t = vec!(0,1))]
    repetition: Vec<usize>,

    #[arg(long, allow_hyphen_values(true))]
    support_e: Option<i64>,

    #[arg(short, long, default_value_t = String::from("../999-pair/"))]
    grid_folder: String,

    #[arg(short, long, default_value_t = String::from("../input_cluster/bulk.poscar"))]
    core_file: String,

    #[arg(short, long, default_value_t = false)]
    write_snap_shots: bool,

    #[arg(long, default_value_t = false)]
    heat_map: bool,

    #[arg(long)]
    coating: Option<String>,

    #[arg(short, long, value_delimiter = '/', default_values_t = vec!(1,2))]
    optimization_cut_off_fraction: Vec<u64>,

    #[arg(short, long, allow_hyphen_values = true)]
    unique_levels: i32,
}

fn file_paths(
    grid_folder: String,
) -> (
    String,
    String,
    String,
    // String,
    String,
    String,
    String,
    String,
) {
    (
        format!("{}nearest_neighbor", grid_folder),
        format!("{}next_nearest_neighbor", grid_folder),
        format!("{}nn_pairlist", grid_folder),
        // format!("{}nnn_pairlist", grid_folder),
        format!("{}atom_sites", grid_folder),
        format!("{}nn_pair_no_intersec", grid_folder),
        format!("{}nnn_gcn_no_intersec.json", grid_folder),
        format!("{}surrounding_moves.json", grid_folder),
    )
}

fn unpack_atoms_input(atoms: Vec<u32>) -> (Option<u32>, Option<Vec<u32>>) {
    if atoms.len() == 1 {
        return (Some(atoms[0]), None);
    } else if atoms.len() == 4 {
        return (Some(atoms[0]), Some(vec![atoms[1], atoms[2], atoms[3]]));
    } else {
        panic!("wrong atoms input, user one of the two input options: \n number of atoms: '-a x' \n or number of atoms with miller indices: '-a x h k l' ")
    }
}
use std::env;
fn main() {
    // enable_data_collection(true);
    env::set_var("RUST_BACKTRACE", "1");
    println!("determined next-nearest neighbor list");

    let args = Args::parse();
    let save_folder: String = args.folder;
    let temperature: f64 = args.temperature;
    let unique_levels = args.unique_levels;
    if !std::path::Path::new(&save_folder).exists() {
        fs::create_dir_all(&save_folder).unwrap();
    }
    let support_e = args.support_e.unwrap_or(0);

    // let (atoms_input, sup) =
    let (atoms_input, support_indices) = if let Some(atoms_input) = args.start_structure.atoms_input
    {
        unpack_atoms_input(atoms_input)
    } else {
        (None, None)
    };

    let input_file: Option<String> = args.start_structure.start_cluster;

    #[allow(unused_variables)]
    let (
        nn_file,
        nnn_file,
        nn_pairlist_file,
        // nnn_pairlist_file,
        atom_sites,
        nn_pair_no_int_file,
        nnn_pair_no_int_file,
        surrounding_moves_file,
    ) = file_paths(args.grid_folder);

    let coating: Option<String> = args.coating;
    let niter_str = args.iterations;
    let niter = fmt_scient(&niter_str);
    let mut write_snap_shots: bool = args.write_snap_shots;
    let heat_map: bool = args.heat_map;
    if heat_map {
        write_snap_shots = true;
    }
    let bulk_file_name: String = args.core_file;
    let optimization_cut_off_fraction: Vec<u64> = args.optimization_cut_off_fraction;
    let repetition = args.repetition;

    let mut atom_names: HashMap<String, u8> = HashMap::new();
    let energy = EnergyInput::Gcn(vec![[0_i64; 145]]);
    // let (energy, mut atom_names) = if args.e_l_cn.is_some() {
    //     let (energy, atom_names) = collect_energy_values([0; 2], args.e_l_cn.unwrap());
    //     (EnergyInput::LinearCn(energy), atom_names)
    // } else if args.e_cn.is_some() {
    //     let (energy, atom_names) = collect_energy_values([0; 13], args.e_cn.unwrap());
    //     (EnergyInput::Cn(energy), atom_names)
    // } else if args.e_l_gcn.is_some() {
    //     let (energy, atom_names) = collect_energy_values([0; 2], args.e_l_gcn.unwrap());
    //     (EnergyInput::LinearGcn(energy), atom_names)
    // } else if args.e_gcn.is_some() {
    //     let (energy, atom_names) = collect_energy_values([0; 145], args.e_gcn.unwrap());
    //     (EnergyInput::Gcn(energy), atom_names)
    // } else {
    //     panic!("no energy")
    // };
    //
    // println!("energy: {:?}", energy);
    println!("{:?}", repetition);
    let mut handle_vec = Vec::new();

    let gridstructure = GridStructure::new(
        nn_file,
        // nnn_file,
        nn_pair_no_int_file,
        // nnn_pair_no_int_file,
        atom_sites,
        bulk_file_name,
        surrounding_moves_file,
    );
    let gridstructure = Arc::new(gridstructure);
    // atom_names.insert("Pt".to_string(), 1);
    // atom_names.insert("Pd".to_string(), 2);
    atom_names.insert("Al".to_string(), 100);

    // let alphas_file = "Pt_Pd.bat".to_string();
    let alphas_file = args.alphas;
    let alphas_arr = read_alphas(alphas_file, &mut atom_names);
    println!("{:?}", alphas_arr);
    let alphas = alpha_energy::Alphas::new(alphas_arr);
    // let alphas = alpha_energy::Alphas::new(alpha_energy::energy_const);
    println!(
        "alphas: {:?} \n sum alphas: {:?}",
        alphas.div_by_cn, alphas.summed_to_x_div_cn
    );
    let alphas_arc = Arc::new(alphas);

    for rep in repetition[0]..repetition[1] {
        let input_file = input_file.clone();
        let save_folder = save_folder.clone();
        let optimization_cut_off_fraction = optimization_cut_off_fraction.clone();
        let gridstructure_arc = Arc::clone(&gridstructure);
        let alphas_arc = Arc::clone(&alphas_arc);
        let support_indices = support_indices.clone();
        let atom_names = atom_names.clone();
        let coating = coating.clone();

        handle_vec.push(thread::spawn(move || {
            let mut sim = Simulation::new(
                atom_names,
                niter,
                input_file,
                atoms_input,
                temperature,
                save_folder,
                write_snap_shots,
                heat_map,
                rep,
                optimization_cut_off_fraction,
                alphas_arc,
                support_indices,
                gridstructure_arc,
                coating,
                support_e,
            );
            let exp = sim.run(unique_levels);
            sim.write_exp_file(&exp);
        }));
    }
    for handle in handle_vec {
        handle.join().unwrap();
    }
    // mc::find_simulation_with_lowest_energy(save_folder).unwrap_or_else(|err| {
    //     println!(
    //         "{:?}",
    //         format!("deleting folders with heigh energy not successful {err}")
    //     )
    // });
}

fn read_alphas(alphas_file: String, atom_names: &mut HashMap<String, u8>) -> [[[f64; 12]; 2]; 2] {
    const LINE_COUNT: usize = 14;
    // let mut x = alphas_file.split('.');
    let fiel_name = alphas_file.split('.').next().unwrap();
    for (i, metal) in fiel_name.split('_').enumerate() {
        atom_names.insert(metal.to_string(), (i + 1) as u8);
    }
    let pairlist = fs::File::open(alphas_file).expect("Should have been able to read the file");

    let lines = BufReader::new(pairlist);
    let mut alphas: [[[f64; 12]; 2]; 2] = [[[0.; 12]; 2]; 2];

    for (i, line) in lines.lines().enumerate() {
        println!("{}", i);
        let r = line.unwrap();
        let num = r.parse::<f64>().unwrap();
        println!("{}", num);
        if i < LINE_COUNT {
            if i >= 12 {
                continue;
            }
            alphas[0][0][i] = num;
        } else if i < LINE_COUNT * 2 {
            if i >= LINE_COUNT * 2 - 2 {
                continue;
            }
            alphas[1][1][i - LINE_COUNT * 1] = num;
        } else if i < LINE_COUNT * 3 {
            if i >= LINE_COUNT * 3 - 2 {
                continue;
            }
            alphas[1][0][i - LINE_COUNT * 2] = num;
        } else {
            if i >= LINE_COUNT * 4 - 2 {
                continue;
            }
            alphas[0][1][i - LINE_COUNT * 3] = num;
        }
    }
    alphas
}

// fn alphas_read_one_section(line_index: usize, section: usize,  line: String, alphas: &mut [[[0.; 12]; 2]; 2]) {
//             if line_index >= LINE_COUNT * section - 2 {
//                 continue;
//             }
//             alphas[0][1][i - LINE_COUNT * (section-1)] = line.parse::<f64>().unwrap();
//
// }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_e_file() {
        let alphas_inp = "Pt_Pd.3.bat".to_string();

        let mut atom_names: HashMap<String, u8> = HashMap::new();
        let res = read_alphas(alphas_inp, &mut atom_names);
        println!("{:?}", atom_names);

        println!("{:?}", res);
    }
}
