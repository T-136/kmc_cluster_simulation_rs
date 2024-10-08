# kMC nanoparticle simulation

This kinetic Monte-Carlo simulation uses [alpha energies](https://pubs.acs.org/doi/abs/10.1021/acscatal.3c04337) to simulate the structural changes of nanoparticles alloys.  

<img src="https://github.com/user-attachments/assets/848a53df-a2df-49a5-b222-5c92cb6b5b69" width="300"  />


## Overview 

### Accurate energy model
The [alpha energy](https://pubs.acs.org/doi/abs/10.1021/acscatal.3c04337) model gives an accurate description of atom bindings. These values have to be calculated once for the metal in question and can then be used to simulate any structure. 

### High performance
Despite the energy taking the surrounding of the neighbors into account the Simulation is able to run 10^10 Iterations in one day. Enabling a longer Simulated time.

### Multithreaded
Running multiple Simulations in different threads reduces the required memory as the largest memory usage comes from loading the precomputed grid structure. By running multiple simulations in parallel the grid structure only has to be read once. You can run multiple temperatures in parallel an run every temperature multiple times to ensure a more representative result.

### Track the Simulation using the xyz format or a compact binary
Throughout the simulation it is possible to track the energy and save snapshots. While more snapshots will give a more detailed representation of the simulation they also slow it down and increase the size of the output-file. 

### Outlook
More accurate TST representation.

## Usage

### Requirements

- [Rust](https://www.rust-lang.org/tools/install)
- [Python](https://www.python.org/) used version 3.11.10
- [ASE](https://wiki.fysik.dtu.dk/ase/) *used version: ase-3.22.1*

### Install from git 
```bash
git clone git@github.com:T-136/kmc_cluster_simulation_rs.git
```

### Build the binary

Build the program with:
```bash
cargo build -r
```
After the build step is complete the compiled program can be found in your project folder under "./target/release/mc".


### Run the simulation
```bash
./target/release/mc -s ./example_data/711_1pt_arround_pd.xyz -t 600,700  -i 1e8 -r 1  --alphas ./alpha_alt.json  -g ../303030-pair_kmc/ -w
```

Use "-h" or "--help" to see the available flags and how to use them. 
```
-s, --start-cluster <START_CLUSTER>
          xyz file
  -a, --atoms-input <ATOMS_INPUT>
          number of atoms
  -f, --folder <FOLDER>
          folder to store results [default: ./sim/]
  -i, --iterations <ITERATIONS>
          acceptes integer and scientific formats e.g. 1E10 or 1.5E10
      --alphas <ALPHAS>
          path to the file with the alpha values. File name has to conatain the metals names in coorect order in the first parte of the file name before the ".", seperated by "_" e.g. Pt_Pd.bat
  -t, --temperatures <TEMPERATURES>
          temperature at which the sumulation is run Set multiple temperatures seperated by comma to run simulations with different temperature in parallel [default: 300]
  -r, --repetition <REPETITION>
          How often the simulation should be repeated. This can be used to make sure the result is representative. The repetition will be applied to each temperature. Therefore the number of paralle run simulations is the number of temperatures times the repetition. To enable all simulations to run in parallel there need to be equal or more threads then simulations running [default: 1]
      --freeze <FREEZE>
          freez a boarder in one ore multiple directions. Directions are given by x, y, z, -x, -y or -z
      --support-e <SUPPORT_E>
          support energy. not supported yet
  -g, --grid-folder <GRID_FOLDER>
          folder which contains the files for the grid previously build using the python scipt [default: ../303030-pair_kmc/]
      --write-binary-snapshots <WRITE_BINARY_SNAPSHOTS>
          Set how many snapshots are saved in each simulation. Snapshots are spread out equally throughout the simulation. Saving the snapshots as binary file requires the used grid to visualize them but is more space as long as the cluster has not less then 25 atoms then there are atoms in the grid
  -w, --write-xyz-snapshots <WRITE_XYZ_SNAPSHOTS>
          Set how many snapshots are saved in each simulation. Snapshots are spread out equally throughout the simulation. Saving the snapshots as xyz file can lead to large files. When the number of atoms is less then 25 times the number of grid positions writing to the binary format will save disk space
      --time-energy-sections <TIME_ENERGY_SECTIONS>
          Determines how often the energy of the system and its elapsed time will be saved in the exp file. If this arg is not set the value from write_binary_snapshots/write_xyz_snapshots will be used. If none is availabel 1000 sections is the default
      --coating <COATING>
          adds one layer of atoms arround the cluster. Input muust be the coating atom e.g. "Pt"
  -h, --help
          Print help (see more with '--help')
  -V, --version
          Print version`
```

## Visualisation
When the "write-xyz-snapshots" flag is activated the simulation saves snapshots equally spread throughout the simulation. These snapshots can be visualized using [ASE](https://wiki.fysik.dtu.dk/ase/).

If the "write-binary-snapshots" flag is activated snapshots will be saved as a binary file which can be viewed  using the provided python script "./view_snapshots.py".

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
## License


