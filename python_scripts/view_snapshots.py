#!/opt/homebrew/bin/python3
import sys
from ase.visualize import view
import argparse
from ase import Atoms, Atom
import csv
import os
import struct


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

# for i, atom in enumerate(atoms_list):
#     atoms_dict[i] = atom
# atoms_dict[100] = "Al"

print(atoms_dict)

index_xyz = []
atoms_file = grid_folder + "atom_sites"
with open(atoms_file) as f:
    for line in f:
        xyz = []
        for coordinate in line.split():
            xyz.append(float(coordinate))
        index_xyz.append(xyz)
snap_shots = []
snap_file = exp_folder + "snap_shot_sections.csv"


class ByteBinaryFileReader:
    def __init__(self, file_path, block_size=4096):
        self.file_path = file_path
        self.block_size = block_size

    def __iter__(self):
        with open(self.file_path, "rb") as f:
            while True:
                data = f.read(self.block_size)
                if not data:
                    break
                for byte in data:
                    yield byte.to_bytes(1, "big")


binary_reader = ByteBinaryFileReader(snap_file)

line_count = 0
cur_len = 0
snap_shots.append([])
is_pre = True
pre_strig_list = []
for val in binary_reader:
    # print(val)
    if val == b"\xfe":
        is_pre = False
    if is_pre:
        pre_strig_list.append(val)
        continue

    cur_len += 1
    if cur_len >= len(index_xyz):
        len_first = cur_len - len(index_xyz) + 1
        snap_shots[-1].append(int.from_bytes(val))
        snap_shots.append([])
        cur_len = 0
        line_count += 1
        print(line_count)
    else:
        snap_shots[-1].append(int.from_bytes(val))
print(len(index_xyz))
print(f"len snapshot {len(snap_shots[1])}")
print(f"len snapshots {len(snap_shots)}")
snap_shots.pop()

print(pre_strig_list)

old_comma_index = 0
# comma_index = pre_strig_list.index(b"\x2c")
for _ in range(pre_strig_list.count(b"\x3a")):
    comma_index = pre_strig_list[old_comma_index:].index(b"\x2c") + old_comma_index
    colon_index = pre_strig_list[old_comma_index:].index(b"\x3a") + old_comma_index
    print(old_comma_index)
    print(comma_index)
    print(colon_index)
    atom_name_string = "".join(
        [x.decode() for x in pre_strig_list[(old_comma_index):(colon_index)]]
    )
    atoms_dict[int.from_bytes(pre_strig_list[colon_index + 1])] = atom_name_string
    old_comma_index = comma_index + 1
print(atoms_dict)

# with open(snap_file, "rb") as f:
#     chunk = sys.maxsize / 8
#     byte = f.readlines()
#     cur_len = 0
#     line_count = 0
#     snap_shots.append([])
#     while byte != b"":
#         # cur_len += chunk
#         for line in byte:
#             for val in line:
#                 cur_len += 1
#                 if cur_len >= len(index_xyz):
#                     len_first = cur_len - len(index_xyz) + 1
#                     # print(
#                     #     f"len_first-1 {len_first} len(byte[:len_first]) { len(byte[:len_first])}"
#                     # )
#                     # ints = struct.unpack(f"{len_first}B", byte[:len_first])
#                     snap_shots[-1].append(int(val))
#                     snap_shots.append([])
#                     cur_len = 0
#                     line_count += 1
#                     # print(
#                     #     f"chunk - len_first {chunk-(len_first)} len(byte[:len_first]) { len(byte[len_first:])}"
#                     # )
#                     # ints = struct.unpack(f"{chunk-(len_first)}B", byte[len_first:])
#                     snap_shots[-1].append(int(val))
#                 else:
#                     # ints = struct.unpack(f"{chunk}B", byte)
#                     # ints = struct.unpack(f"{chunk}B", byte)
#                     snap_shots[-1].append(int(val))

# with open(snap_file, newline="") as csvfile:
#     csv.field_size_limit(sys.maxsize)
#     csv_reader = csv.reader(csvfile, delimiter=",", quotechar="|")
#     for row in csv_reader:
#         snap_shots.append([int(x) for x in row])

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
