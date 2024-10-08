#!/usr/bin/env python
from ase import Atoms, Atom
from ase.io import Trajectory, read
import os
import numpy as np
import sys
import copy
import shutil
import json


def create_nn_dict(sites, nn_pair_folder, nsites, atoms) -> dict:
    use_pbc = True
    nn = {}
    dthr = 3.0
    # dthr= 2.52
    for site in sites:
        print(site, "nn of", nsites)
        if site not in nn:
            nn[site] = []
        # for j in range(site+1, nsites):  # start at i+1 !!!!
        # dists[ji,j]=
        range_atoms = range(site + 1, nsites)
        if len(range_atoms) == 0:
            continue
        distance_list = atoms.get_distances(site, list(range_atoms), mic=use_pbc)
        for i, j in enumerate(range_atoms):
            if j not in nn:
                nn[j] = []
            if distance_list[i] < dthr:
                nn[j] += [site]
                nn[site] += [j]
        if len(nn[site]) != 12:
            print(f"nn not == 12 nn: {len(nn[site])}")
            sys.exit()
    f = open(nn_pair_folder + "nearest_neighbor", "w")
    for site in sites:
        f.write(" ".join(str(j) for j in [site] + nn[site]) + "\n")
    f.close()
    return nn


def create_nn_pair(nn, sites, nsites, nn_pair_folder):
    nn_pair = {}
    for site in sites:
        print(site, "nn_pair of", nsites)
        nn_pair[site] = {}
        # print(nn[site])
        for nn_of_site in nn[site]:
            site_neigbors = []
            nn_of_site_neigbors = []
            for o1 in nn[site]:
                # print(o1)
                # print(f"all: {o1}")
                if o1 == nn_of_site or o1 == site:
                    continue
                if o1 not in site_neigbors:
                    site_neigbors.append(o1)

            for o2 in nn[nn_of_site]:
                if o2 == nn_of_site or o2 == site:
                    # print(o1)
                    continue
                if o2 not in nn_of_site_neigbors:
                    nn_of_site_neigbors.append(o2)

            print(len(nn_of_site_neigbors))
            min_val = min([site, nn_of_site])
            max_val = max([site, nn_of_site])
            print(site)
            print(nn_of_site_neigbors)
            if len(nn_of_site_neigbors) > 11 or len(site_neigbors) > 11:
                sys.exit()

            print([site_neigbors, nn_of_site_neigbors])
            if min_val == site:
                nn_pair[min_val][max_val] = [site_neigbors, nn_of_site_neigbors]
            if min_val == nn_of_site:
                nn_pair[min_val][max_val] = [nn_of_site_neigbors, site_neigbors]

    return nn_pair


def create_neighbor_pairs(nn, sites, nsites, nn_pair_folder):
    nnn_pair = {}
    nn_pair = {}
    neighbor = []
    for site in sites:
        print(site, "nn of", nsites)
        nnn_pair[site] = {}
        nn_pair[site] = {}
        for nn_of_site in nn[site]:
            tocheck = [site, nn_of_site]
            for ks in nn[site] + nn[nn_of_site]:
                if ks not in tocheck:
                    tocheck += [ks]

            nn_pair[min([site, nn_of_site])][max([site, nn_of_site])] = tocheck
            if [min([site, nn_of_site]), max([site, nn_of_site])] not in neighbor:
                neighbor.append([min([site, nn_of_site]), max([site, nn_of_site])])
            if [max([site, nn_of_site]), min([site, nn_of_site])] not in neighbor:
                neighbor.append([max([site, nn_of_site]), min([site, nn_of_site])])

            oldcheck = list(tocheck)
            tocheck = []
            for ks in oldcheck:
                for m in nn[ks]:
                    if m not in tocheck:
                        tocheck += [m]
            nnn_pair[min([site, nn_of_site])][max([site, nn_of_site])] = tocheck

    f = open(nn_pair_folder + "nnn_pairlist", "w")
    for site in sites:
        for j in nnn_pair[site]:
            f.write(" ".join(str(j) for j in [site, j] + nnn_pair[site][j]) + "\n")
    f.close()
    f = open(nn_pair_folder + "nn_pairlist", "w")
    for site in sites:
        for j in nn_pair[site]:
            f.write(" ".join(str(j) for j in [site, j] + nn_pair[site][j]) + "\n")
    f.close()

    neighbor.sort(key=lambda x: x[0])
    with np.printoptions(threshold=np.inf):
        print(neighbor)
    f = open(nn_pair_folder + "neighbor", "w")

    for n in neighbor:
        f.write(" ".join(str(j) for j in n) + "\n")
    f.close()

    return (nn_pair, nnn_pair, neighbor)


def sort_atoms_by_positions(atoms):
    positions = atoms.get_positions()
    new_positions = sorted(positions, key=lambda x: [x[2], x[1], x[0]])
    # my_list = np.array(positions)
    # ind = np.lexsort((my_list[:,0], my_list[:,1], my_list[:,2]))
    # sorted = my_list[ind]
    atoms.positions = new_positions
    # orderd_atoms = Atoms(cell=atoms.cell, pbc=atoms.pbc)
    # for position in sorted:
    #     atoms.append(Atom('Pt', position))

    with np.printoptions(threshold=np.inf):
        # print(orderd_atoms.get_positions())
        # print(orderd_atoms)
        print(atoms.get_positions())
        print(atoms)
    return atoms

    # positions.sort(positions, key=lambda x: x[0])
    # for number_of_slice in range(1, positions.len()/grid_size):
    #     positions.sort(key=lambda x: x[(number_of_slice-1)*positions:number_of_slice*positions][1]))
    # for number_of_slice in range(1, positions.len()/grid_size):
    #     positions.sort(key=lambda x: x[(number_of_slice-1)*positions:number_of_slice*positions][1]))


def create_nn_pair_without_intersection(nn, sites, nsites, nn_pair_folder):
    nn_pair_no_intersec = {}
    nn_pair_only_intersec = {}
    nn_pair = {}
    for site in sites:
        print(site, "nn no int-sec of", nsites)
        nn_pair[site] = {}
        nn_pair_only_intersec[site] = {}
        nn_pair_no_intersec[site] = {}
        for nn_of_site in nn[site]:
            no_intersect_list = {}
            no_intersect_list[site] = []
            no_intersect_list[nn_of_site] = []
            intersect_list = {}
            intersect_list[site] = []
            intersect_list[nn_of_site] = []
            only_intersect = {}
            only_intersect[site] = []
            only_intersect[nn_of_site] = []
            for x in nn[site]:
                if x == nn_of_site:
                    continue
                intersect_list[site].append(x)
                if x in nn[nn_of_site]:
                    only_intersect[site].append(x)
                    continue
                no_intersect_list[site].append(x)

            for x in nn[nn_of_site]:
                if x == site:
                    continue
                intersect_list[nn_of_site].append(x)
                if x in nn[site]:
                    # only_intersect[nn_of_site].append(x)
                    continue
                no_intersect_list[nn_of_site].append(x)

            for x in only_intersect[site]:
                if (
                    x
                    in no_intersect_list[
                        min([site, nn_of_site])
                        or x in no_intersect_list[max([site, nn_of_site])]
                    ]
                ):
                    print("create nnpair: intersect and no intertsect")
                    sys.exit()
            if len(only_intersect[site]) + len(no_intersect_list[nn_of_site]) != 11:
                print(len(only_intersect[site]))
                print(len(no_intersect_list[site]))
                print("only + without intertsect not 12")
                sys.exit()

            nn_pair[min([site, nn_of_site])][max([site, nn_of_site])] = [
                intersect_list[min([site, nn_of_site])],
                intersect_list[max([site, nn_of_site])],
            ]
            nn_pair_only_intersec[min([site, nn_of_site])][max([site, nn_of_site])] = (
                only_intersect[site]
            )
            nn_pair_no_intersec[min([site, nn_of_site])][max([site, nn_of_site])] = [
                no_intersect_list[min([site, nn_of_site])],
                no_intersect_list[max([site, nn_of_site])],
            ]
    f = open(nn_pair_folder + "nn_pair_no_intersec", "w")
    for site in sites:
        for j in nn_pair_no_intersec[site]:
            f.write(
                " ".join(
                    str(j)
                    for j in [site, j]
                    + nn_pair_no_intersec[site][j][0]
                    + nn_pair_no_intersec[site][j][1]
                    + nn_pair_only_intersec[site][j]
                )
                + "\n"
            )
    f.close()
    return nn_pair_no_intersec, nn_pair, nn_pair_only_intersec


def moves_k_effected_by_move(
    nn, nsites, nn_pair_folder, nn_pair_no_intersec, nn_pair_only_intersec
):
    sourrounding_moves = {}
    effected_by_both = copy.deepcopy(nn_pair_only_intersec)
    for min_site, d in nn_pair_no_intersec.items():
        # for min_site, d in nn_pair.items():
        sourrounding_moves[min_site] = {}
        for max_site, nn_of_pair in d.items():
            sourrounding_moves[min_site][max_site] = []
            print(max_site, "nnn no int-sec of", nsites)
            for neighbor_fromatom in nn_of_pair[0]:
                for target in nn[neighbor_fromatom]:
                    move = (
                        min([neighbor_fromatom, target]),
                        max([neighbor_fromatom, target]),
                    )
                    if move in sourrounding_moves[min_site][max_site]:
                        continue
                    if target == min_site or target == max_site:
                        continue
                    sourrounding_moves[min_site][max_site].append(move)

            for neighbor_toatom in nn_of_pair[1]:
                for target in nn[neighbor_toatom]:
                    move = (
                        min([neighbor_toatom, target]),
                        max([neighbor_toatom, target]),
                    )
                    if move in sourrounding_moves[min_site][max_site]:
                        continue
                    if target == min_site or target == max_site:
                        continue
                    sourrounding_moves[min_site][max_site].append(move)

    import json

    for x in sourrounding_moves[0].values():
        print(len(x))
    # with open('nnn_gcn_no_intersec.json', 'w') as fp:
    j = json.dumps(sourrounding_moves, separators=(",", ":"))
    # print(j)
    with open(nn_pair_folder + "surrounding_moves.json", "w") as f:
        print(j, file=f)
    return sourrounding_moves


def moves_k_effected_by_change(nn, nsites, nn_pair_folder, sites):
    sourrounding_moves = {}
    for site in sites:
        # for min_site, d in nn_pair.items():
        sourrounding_moves[site] = []
        print(site, "k eff by change", nsites)
        for neighbor_fromatom in nn[site]:
            for target in nn[neighbor_fromatom]:
                move = (
                    min([neighbor_fromatom, target]),
                    max([neighbor_fromatom, target]),
                )
                # print(move)
                if move in sourrounding_moves[site]:
                    continue
                if target == site:
                    continue
                print(move)
                sourrounding_moves[site].append(move)

    print(sourrounding_moves)
    # for x in sourrounding_moves[0].values():
    #     print(len(x))
    j = json.dumps(sourrounding_moves, separators=(",", ":"))
    # print(j)
    with open(nn_pair_folder + "pos_surrounding_moves.json", "w") as f:
        print(j, file=f)
    return sourrounding_moves


if __name__ == "__main__":
    # how often should the bulk file be repreated in x,y,z directions
    nlattice = [9, 9, 9]
    bulk_file = "./../input_cluster/bulk_sym.poscar"

    # bulk_file = './../input_cluster/111-bulk_fix.poscar'
    nn_pair_folder = "./202020-pair/"  # output folder

    if not os.path.exists(nn_pair_folder):
        os.mkdir(nn_pair_folder)

    shutil.copyfile(bulk_file, nn_pair_folder + "bulk.poscar")

    buld_file = read(bulk_file)
    atoms = read(bulk_file)
    atoms = atoms.repeat(nlattice)
    # atoms.translate([-atoms.cell[0][0]/2, -atoms.cell[1][1]/2, -atoms.cell[2][2]/2])
    # atoms.center()
    # new_cell = atoms.cell
    # new_cell[0][0] -= atoms.cell[0][0]/2
    # new_cell[1][1] -= atoms.cell[1][1]/2
    # new_cell[2][2] -= atoms.cell[2][2]/2
    # print(atoms.cell)
    # print(new_cell)
    # atoms.cell = new_cell
    # atoms = sorted(atoms, key=lambda x: x.get_positions()[0])
    atoms2 = sort_atoms_by_positions(atoms)
    # atoms2.write(nn_pair_folder + "total.traj")
    atoms2.write(nn_pair_folder + "grid_file.xyz")
    xsites_positions = atoms2.get_positions()

    print(xsites_positions)

    with open(nn_pair_folder + "atom_sites", "w") as f:
        np.savetxt(nn_pair_folder + "atom_sites", xsites_positions)

    nsites = len(atoms2)
    sites = range(nsites)

    nn = create_nn_dict(sites, nn_pair_folder, nsites, atoms2)
    nn_pair_no_intersec, nn_pair, nn_pair_only_intersec = (
        create_nn_pair_without_intersection(nn, sites, nsites, nn_pair_folder)
    )
    sourrounding_moves = moves_k_effected_by_move(
        nn, nsites, nn_pair_folder, nn_pair_no_intersec, nn_pair_only_intersec
    )
    pos_sourrounding_moves = moves_k_effected_by_change(
        nn, nsites, nn_pair_folder, sites
    )
