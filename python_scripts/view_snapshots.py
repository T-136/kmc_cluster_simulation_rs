#!/opt/homebrew/bin/python3
import sys
from ase.visualize import view
import argparse
from ase import Atoms, Atom
import csv
import os


parser = argparse.ArgumentParser()
parser.add_argument(
    "dir",
    help="folder of the simualtion which contains the snap_shot_sections.csv-file e.g. './sim/900K_15000000I_711A_0/'",
    default=None,
)  # default set path later
parser.add_argument(
    "-i",
    "--index_sections",
    help="-1 is last, can be multiple seperated by blank",
    nargs="+",
    default=[-1],
    type=int,
)
parser.add_argument(
    "-g",
    "--grid",
    help="basis for reading data, has to be the same file used in the simualtion",
    default="/Users/tilman/Documents/Vertiefer/cluster_simulations/202020-pair/",
)
parser.add_argument(
    "-a",
    "--alphas_file",
    help="has to be the same file used in the simualtion",
    default=None,
)  # default set path later


exp_folder = parser.parse_args().dir
index = parser.parse_args().index_sections
grid_folder = parser.parse_args().grid
alphas_file = parser.parse_args().alphas_file

bare_alphas_file = os.path.basename(alphas_file)

atoms_list = bare_alphas_file.split(".")[0].split("_")
atoms_dict = {}
for i, atom in enumerate(atoms_list):
    atoms_dict[i] = atom
atoms_dict[100] = "Al"

print(atoms_dict)

snap_shots = []
snap_file = exp_folder + "snap_shot_sections.csv"
with open(snap_file, newline="") as csvfile:
    csv.field_size_limit(sys.maxsize)
    csv_reader = csv.reader(csvfile, delimiter=",", quotechar="|")
    for row in csv_reader:
        snap_shots.append([int(x) for x in row])

index_xyz = []
atoms_file = grid_folder + "atom_sites"
with open(atoms_file) as f:
    for line in f:
        xyz = []
        for coordinate in line.split():
            xyz.append(float(coordinate))
        index_xyz.append(xyz)

atoms = []
for i in index:
    atom_list = []
    for i, x in enumerate(snap_shots[i]):
        if x == 255:
            continue
        atom_list.append(Atom(atoms_dict[x], position=index_xyz[i]))

    atoms.append(Atoms(atom_list))
    print(len(atoms[-1]))

view(atoms)
