#!/opt/homebrew/bin/python3
from ase.visualize import view
import argparse
from ase import Atoms, Atom


def kmcAtomsView(exp_folder, grid_folder, block_size=4096):
    atoms = []
    atoms_dict = {}
    index_xyz = []
    atoms_file = grid_folder + "atom_sites"
    with open(atoms_file) as f:
        for line in f:
            xyz = []
            for coordinate in line.split():
                xyz.append(float(coordinate))
            index_xyz.append(xyz)

    snap_file = exp_folder + "snapshot_sections"
    binary_reader = BinaryFileReader(snap_file, len(index_xyz), block_size)

    old_comma_index = 0
    # comma_index = pre_strig_list.index(b"\x2c")
    for _ in range(binary_reader.pre_string_list.count(b"\x3a")):
        comma_index = (
            binary_reader.pre_string_list[old_comma_index:].index(b"\x2c")
            + old_comma_index
        )
        colon_index = (
            binary_reader.pre_string_list[old_comma_index:].index(b"\x3a")
            + old_comma_index
        )
        print(old_comma_index)
        print(comma_index)
        print(colon_index)
        atom_name_string = "".join(
            [
                x.decode()
                for x in binary_reader.pre_string_list[(old_comma_index):(colon_index)]
            ]
        )
        atoms_dict[int.from_bytes(binary_reader.pre_string_list[colon_index + 1])] = (
            atom_name_string
        )
        old_comma_index = comma_index + 1
    print(atoms_dict)

    for i in index:
        atom_list = []
        for i, x in enumerate(binary_reader.snap_shots[i]):
            if x == 255:
                continue
            atom_list.append(Atom(atoms_dict[x], position=index_xyz[i]))

        atoms.append(Atoms(atom_list))
    return atoms


class BinaryFileReader:
    snap_shots = []
    pre_string_list = []

    def __init__(self, file_path, grid_len, block_size=4096):
        self.grid_len = grid_len
        self.file_path = file_path
        self.block_size = block_size

        line_count = 0
        cur_len = 0
        self.snap_shots.append([])
        is_pre = True
        for val in self.open_file():
            if val == b"\xfe":
                is_pre = False
            if is_pre:
                self.pre_string_list.append(val)
                continue

            cur_len += 1
            if cur_len >= self.grid_len:
                len_first = cur_len - self.grid_len + 1
                self.snap_shots[-1].append(int.from_bytes(val))
                self.snap_shots.append([])
                cur_len = 0
                line_count += 1
                print(line_count)
            else:
                self.snap_shots[-1].append(int.from_bytes(val))
        print(self.grid_len)
        self.snap_shots.pop()
        print(f"len each snapshot {len(self.snap_shots[1])}")
        print(f"ammount snapshots {len(self.snap_shots)}")

    def open_file(self):
        with open(self.file_path, "rb") as f:
            while True:
                data = f.read(self.block_size)
                if not data:
                    break
                for byte in data:
                    yield byte.to_bytes(1, "big")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "dir",
        help="folder of the simualtion which contains the snap_shot_sections.csv-file e.g. './sim/900K_15000000I_711A_0/'",
        default=None,
    )
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
    )

    exp_folder = parser.parse_args().dir
    index = parser.parse_args().index_sections
    grid_folder = parser.parse_args().grid
    alphas_file = parser.parse_args().alphas_file

    atoms = kmcAtomsView(exp_folder, grid_folder)
    view(atoms)
