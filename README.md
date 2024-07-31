# kMC cluster simulation

This kinetic Monte-Carlo simulation uses [alpha energies](https://pubs.acs.org/doi/abs/10.1021/acscatal.3c04337) to simulate the structural changes of nanoparticles alloys.  

<img src="https://github.com/user-attachments/assets/848a53df-a2df-49a5-b222-5c92cb6b5b69" width="300"  />

## Requirements

- [Rust](https://www.rust-lang.org/tools/install)
- [Python](https://www.python.org/)
- [ASE](https://wiki.fysik.dtu.dk/ase/) *used version: ase-3.22.1*


## Install from git 
```bash
git clone git@github.com:T-136/kmc_cluster_simulation_rs.git
```

## Build the binary

Build the program with:
```bash
cargo build -r
```
After the build step is complete the compiled program can be found in your project folder under "./target/release/mc".


## Run the simulation
```bash
./target/release/mc -s ./example_data/711_1pt_arround_pd.xyz -t 700  -i 1e8 -r 1  --alphas ./Pt_Pd.6.bat  -g ../303030-pair_kmc/ -w
```

```
-s, --start-cluster <START_CLUSTER>

-a, --atoms-input <ATOMS_INPUT>

-f, --folder <FOLDER>
        folder to store results [default: ./sim/]
-i, --iterations <ITERATIONS>
        iterations
-t, --temperature <TEMPERATURE>
        [default: 300]
    --alphas <ALPHAS>

-r, --repetition <REPETITION>
        [default: 0 1]
-f, --freeze <FREEZE>

    --support-e <SUPPORT_E>

-g, --grid-folder <GRID_FOLDER>
        [default: ../999-pair/]
-w, --write-snap-shots

    --heat-map

    --coating <COATING>

-o, --optimization-cut-off-fraction <OPTIMIZATION_CUT_OFF_FRACTION>
        [default: 1 2]
-h, --help
        Print help
-V, --version
        Print version
```

## Visualisation

When the "write-snap-shots" flag is activated the simulation saves 200 snapshots equally spread throughout the simulation. These snapshots can be visualized using [ASE](https://wiki.fysik.dtu.dk/ase/). This can be done using the provided python script "./python_scripts/view_snapshots.py".

```
positional arguments:
  dir                   folder of the simualtion which contains the snap_shot_sections.csv-file e.g. './sim/900K_15000000I_711A_0/'

options:
  -h, --help            show this help message and exit
  -i INDEX_SECTIONS [INDEX_SECTIONS ...], --index_sections INDEX_SECTIONS [INDEX_SECTIONS ...]
                        -1 is last, can be multiple seperated by blank
  -g GRID, --grid GRID  basis for reading data, has to be the same file used in the simualtion
  -a ALPHAS_FILE, --alphas_file ALPHAS_FILE
                        has to be the same file used in the simualtion
```

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

## License


