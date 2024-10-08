// use ahash::HashMap;
use clap::ArgGroup;
use clap::Parser;
use core::panic;
use mc::alpha_energy;
use mc::GridStructure;
use mc::Simulation;
use std::collections::HashMap;
use std::env;
use std::fs;
use std::num::NonZeroU32;
use std::sync::Arc;
use std::thread;

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
    let exp_str = parts.next().unwrap_or("0");
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
    let exp = exp_str.parse::<u32>().expect("wrong iterations input");
    let new_pre_num = pre_num.replace('.', "");
    assert!(
        (after_dot_count as u32) < exp,
        "iteration has to be a whole number"
    );
    let base: u64 = 10;
    new_pre_num.parse::<u64>().expect("wrong iterations input")
        * base.pow(exp - after_dot_count as u32)
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
    /// output folder
    /// each simulation will have its own sub-folder with tempearture iteration atoms and and index
    /// number e.g. "800K_100000000I_12007A_0"
    folder: String,

    #[arg(short, long)]
    /// acceptes integer and scientific formats e.g. 1E10 or 1.5E10
    iterations: String,

    #[arg(long)]
    /// path to the json file with the alpha values.
    alphas: String,

    #[arg(short, long, value_delimiter = ',', default_values_t = vec!(300.))]
    /// temperature at which the sumulation is run
    /// Set multiple temperatures seperated by comma to run simulations with different temperature
    /// in parallel.  
    ///
    /// The repetition will be applied to each temperature. Therefore the number of paralle run
    /// simulations is the number of temperatures times the repetition.
    /// To enable all simulations to run in parallel there need to be equal or more threads then
    /// simulations running.
    temperatures: Vec<f64>,

    #[arg(short, long, value_delimiter = '-', default_values_t = vec!(1))]
    /// How often the simulation should be repeated. This can be used to make sure the result is representative.
    /// The repetition will be applied to each temperature. Therefore the number of paralle run
    /// simulations is the number of temperatures times the repetition.
    /// To enable all simulations to run in parallel there need to be equal or more threads then
    /// simulations running.
    ///
    /// The loaded grid_folder-files will only be loaded once and used by all threads meaning the memory
    /// increase is minimal for addidtional simulation runs.
    /// Alternativley the input can be a range e.g. 5-10. This is usefull because the thread number is part of the simulation name.
    repetition: Vec<usize>,

    #[arg(long, value_delimiter = ',', allow_hyphen_values(true))]
    /// freez a boarder in one ore multiple directions. Directions are given by x, y, z, -x, -y or -z.
    freeze: Option<Vec<String>>,

    #[arg(long, allow_hyphen_values(true))]
    /// support energy. not supported yet.
    support_e: Option<i64>,

    #[arg(short, long, default_value_t = String::from("../303030-pair_kmc/"))]
    /// folder which contains the files for the grid previously build using the python scipt.
    grid_folder: String,

    #[arg(long)]
    /// Set how many snapshots are saved in each simulation.
    /// Snapshots are spread out equally throughout the simulation.
    /// Saving the snapshots as binary file requires the used grid to visualize them but is more space
    /// as long as the cluster has not less then 25 atoms then there are atoms in the grid.
    write_binary_snapshots: Option<u32>,

    #[arg(short, long)]
    /// Set how many snapshots are saved in each simulation.
    /// Snapshots are spread out equally throughout the simulation.
    /// Saving the snapshots as xyz file can lead to large files. When the number of atoms is less
    /// then 25 times the number of grid positions writing to the binary format will save disk
    /// space.
    write_xyz_snapshots: Option<u32>,

    #[arg(long)]
    /// Determines how often the energy of the system and its elapsed time will be saved in the exp
    /// file. If this arg is not set the value from write_binary_snapshots/write_xyz_snapshots will
    /// be used. If none is availabel 1000 sections is the default.
    time_energy_sections: Option<u32>,

    #[arg(long)]
    /// adds one layer of atoms arround the cluster. Input muust be the coating atom e.g. "Pt".
    coating: Option<String>,
}

fn file_paths(grid_folder: String) -> (String, String, String, String, String, String) {
    (
        format!("{}nearest_neighbor", grid_folder),
        format!("{}nn_pairlist", grid_folder),
        format!("{}atom_sites", grid_folder),
        format!("{}nn_pair_no_intersec", grid_folder),
        format!("{}surrounding_moves.json", grid_folder),
        format!("{}grid_file.xyz", grid_folder),
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

fn main() {
    // enable_data_collection(true);
    env::set_var("RUST_BACKTRACE", "1");

    let args = Args::parse();
    let save_folder: String = args.folder;
    let temperatures: Vec<f64> = args.temperatures;
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
        nn_file,
        nn_pairlist_file,
        atom_sites,
        nn_pair_no_int_file,
        surrounding_moves_file,
        grid_file,
    ) = file_paths(args.grid_folder);

    let coating: Option<String> = args.coating;
    let niter_str = args.iterations;
    let niter = fmt_scient(&niter_str);
    let write_xyz_snapshots: Option<NonZeroU32> =
        NonZeroU32::new(args.write_xyz_snapshots.unwrap_or(0));
    let write_binary_snapshots: Option<NonZeroU32> =
        NonZeroU32::new(args.write_binary_snapshots.unwrap_or(0));
    let time_energy_sections: Option<NonZeroU32> =
        NonZeroU32::new(args.time_energy_sections.unwrap_or(0));
    let repetition = args.repetition;

    let repetition = if repetition.len() == 1 {
        prepend(repetition, &[0])
    } else {
        repetition
    };

    let mut atom_names: HashMap<String, u8> = HashMap::new();

    if let Ok(default_parallelism_approx) = thread::available_parallelism() {
        println!(
            "\nestimated available parallellism {}, 
                does not account for other running processes which might occupy threads,
                for accurate information read the specs of the hardware in use or do a test runs to see the performace",
            default_parallelism_approx.get()
        );
    } else {
        println!();
    };
    println!(
        "{} threads are requrired to run everything in parallel ({} repetitions * {} temperatures)\n",
        repetition.len() * temperatures.len(),
        repetition.len(),
        temperatures.len(),
    );

    let gridstructure = GridStructure::new(
        nn_file,
        nn_pair_no_int_file,
        atom_sites,
        surrounding_moves_file,
        grid_file,
    );
    let gridstructure = Arc::new(gridstructure);

    let alphas_file = args.alphas;
    let alphas = alpha_energy::Alphas::new_from_json(alphas_file, &mut atom_names);
    println!(
        "alphas: {:?} \n sum alphas: {:?}",
        alphas.cn, alphas.summed_to_x
    );
    let alphas_arc = Arc::new(alphas);

    let mut handle_vec = Vec::new();
    for temperature in temperatures {
        for rep in repetition[0]..repetition[1] {
            let input_file = input_file.clone();
            let save_folder = save_folder.clone();
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
                    write_xyz_snapshots,
                    write_binary_snapshots,
                    time_energy_sections,
                    rep,
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
    }
    for handle in handle_vec {
        handle.join().unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_e_file() {
        let alphas_inp = "Pt_Pd.3.bat".to_string();

        let mut atom_names: HashMap<String, u8> = HashMap::new();
        let res = alpha_energy::read_alphas(alphas_inp, &mut atom_names);
        println!("{:?}", atom_names);

        println!("{:?}", res);
    }
}
