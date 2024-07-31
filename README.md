# kmc_cluster_simulation_rs

This kinetic Monte-Carlo simulation uses [alpha energies](https://pubs.acs.org/doi/abs/10.1021/acscatal.3c04337) to simulate the structural changes of nanoparticles alloys. 

![4033Aw_700K_-1_2 xyz](https://github.com/user-attachments/assets/848a53df-a2df-49a5-b222-5c92cb6b5b69)

<img src="https://github.com/user-attachments/assets/848a53df-a2df-49a5-b222-5c92cb6b5b69" width="400"  />

## Requirements

[Rust](https://www.rust-lang.org/tools/install)

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
./target/release/mc -s cube_pd_tunnel.xyz -t 900  -i 1e6 -r 0-1  -o 10/10 --alphas ./Pt_Pd.6.bat  -g ../303030-pair_kmc/ -w --support-e -310
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

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

## License


