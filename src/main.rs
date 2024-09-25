// use ahash::HashMap;
use clap::ArgGroup;
use clap::Parser;
use core::panic;
use mc::alpha_energy;
use mc::GridStructure;
use mc::Simulation;
use std::collections::HashMap;
use std::fs;
use std::io::BufRead;
use std::io::BufReader;
use std::path::Path;
use std::sync::Arc;
use std::thread;
use std::usize;

fn prepend<T>(v: Vec<T>, s: &[T]) -> Vec<T>
where
    T: Clone,
{
    let mut tmp: Vec<_> = s.to_owned();
    tmp.extend(v);
    tmp
}

fn fmt_scient(num: &str) -> u64 {
    let mut parts = num.split(['e', 'E']);

    let pre_num = parts.next().unwrap();
    let exp = parts.next().unwrap_or("0");
    if parts.next().is_some() {
        panic!("wrong iterations input");
    }

    let mut prenum_parts = pre_num.split('.');
    let _ = prenum_parts.next().unwrap();
    let after_dot_count = if let Some(y) = prenum_parts.next() {
        y.len()
    } else {
        0
    };
    let new_pre_num = pre_num.replace('.', "");

    let base: u64 = 10;
    new_pre_num.parse::<u64>().expect("wrong iterations input")
        * base.pow(exp.parse::<u32>().expect("wrong iterations input") - after_dot_count as u32)
}

#[derive(Parser, Debug)]
#[clap(group(
        ArgGroup::new("startstructure")
            .required(true)
            .args(&["start_cluster", "atoms_input"]),
    ))]
struct StartStructure {
    #[arg(short, long)]
    /// xyz file
    start_cluster: Option<String>,

    #[arg(short, long, value_delimiter = ',')]
    /// number of atoms
    atoms_input: Option<Vec<u32>>,
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[clap(flatten)]
    start_structure: StartStructure,

    #[arg(short, long, default_value_t = String::from("./sim/"))]
    /// folder to store results
    folder: String,

    #[arg(short, long)]
    /// acceptes integer and scientific formats e.g. 1E10 or 1.5E10
    iterations: String,

    #[arg(short, long, default_value_t = 300.)]
    /// temperature at which the sumulation is run, currently the temperature has to be constatn
    /// throughout the simulation
    temperature: f64,

    #[arg(long)]
    /// path to the file with the alpha values. File name has to conatain the metals names in
    /// coorect order in the first parte of the file name before the ".", seperated by "_" e.g. Pt_Pd.bat
    alphas: String,

    #[arg(short, long, value_delimiter = '-', default_values_t = vec!(1))]
    /// How often the simulation should be repeated. This can be used to make sure the result is representative. 
    /// The input can be the number of threads e.g. 5.
    /// Each run will be started in a seperated thread.
    /// The loaded grid_folder-files will only be loaded once and used by all threads meaning the memory
    /// increase is minimal for addidtional simulation runs.
    /// Alternativley the input can be a range e.g. 5-10. This is usefull because the thread number is part of the simulation name.
    repetition: Vec<usize>,

    #[arg(short, long, value_delimiter = ',', allow_hyphen_values(true))]
    /// freez a boarder in one ore multiple directions. Directions are given by x, y, z, -x, -y or -z.
    freeze: Option<Vec<String>>,

    #[arg(long, allow_hyphen_values(true))]
    /// support energy. not supported yet.
    support_e: Option<i64>,

    #[arg(short, long, default_value_t = String::from("../303030-pair_kmc/"))]
    /// folder which contains the files for the grid previously build using the python scipt.
    grid_folder: String,

    #[arg(short, long, default_value_t = false)]
    /// writes a file with the cluster every x iterations.
    write_snap_shots: bool,

    #[arg(long)]
    /// adds one layer of atoms arround the cluster. Input muust be the atom e.g. "Pt".
    coating: Option<String>,

    #[arg(short, long, value_delimiter = '/', default_values_t = vec!(1,2))]
    /// fraction of the simulationtime where the programm starts looking for the structure with the
    /// lowest energy
    optimization_cut_off_fraction: Vec<u64>,
}

fn file_paths(
    grid_folder: String,
) -> (
    String,
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
        format!("{}bulk.poscar", grid_folder),
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

    let args = Args::parse();
    let save_folder: String = args.folder;
    let temperature: f64 = args.temperature;
    if !std::path::Path::new(&save_folder).exists() {
        fs::create_dir_all(&save_folder).unwrap();
    }
    let support_e = args.support_e.unwrap_or(0);
    let freez = args.freeze;

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
        bulk_file_name,
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
    // let bulk_file_name: String = args.core_file;
    let optimization_cut_off_fraction: Vec<u64> = args.optimization_cut_off_fraction;
    let repetition = args.repetition;

    let repetition = if repetition.len() == 1 {
        prepend(repetition, &[0])
    } else {
        repetition
    };

    let mut atom_names: HashMap<String, u8> = HashMap::new();

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

    let alphas_file = args.alphas;
    let alphas_arr = read_alphas(alphas_file, &mut atom_names);
    println!("{:?}", alphas_arr);
    let alphas = alpha_energy::Alphas::new(alphas_arr);
    println!(
        "alphas: {:?} \n sum alphas: {:?}",
        alphas.cn, alphas.summed_to_x
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
        let freez = freez.clone();

        handle_vec.push(thread::spawn(move || {
            let mut sim = Simulation::new(
                atom_names,
                niter,
                input_file,
                atoms_input,
                temperature,
                save_folder,
                write_snap_shots,
                rep,
                optimization_cut_off_fraction,
                alphas_arc,
                support_indices,
                gridstructure_arc,
                coating,
                support_e,
                freez,
            );
            let exp = sim.run();
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

///returns alphas where alhpas[x][x] is pure atom_x and alhpas[x][y] is atom_x in atom_y
fn read_alphas(alphas_file: String, atom_names: &mut HashMap<String, u8>) -> [[[f64; 12]; 2]; 2] {
    const LINE_COUNT: usize = 14;

    let path = Path::new(&alphas_file);
    let file_name = path.file_name().unwrap();
    // let mut x = alphas_file.split('.');
    let atom_names_string = file_name.to_str().unwrap().split('.').next().unwrap();
    for (i, metal) in atom_names_string.split('_').enumerate() {
        atom_names.insert(metal.to_string(), i as u8);
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
