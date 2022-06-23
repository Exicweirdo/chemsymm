import os
import re
import sys
from Molecule import *
import warnings
from find_sym import find_symmetry


def read_xyz(path: str) -> Molecule:

    if os.path.splitext(path)[1] != ".xyz":
        warnings.warn("the file is not xyz file, found {} instead".format(
            os.path.splitext(path)[1]))

    with open(path, mode='r') as f:
        Natom = int(f.readline().strip(" "))
        atom_list = []
        line0 = f.readline().strip()
        mode = "normal"
        if line0:
            if line0[0] == "1":
                mode = "chemoffice"
                atom_list.append(line0.split()[1:5])

        if mode == "normal":
            for iatom in range(Natom):
                atom = f.readline().strip().split()
                for i in range(1, 4):
                    atom[i] = float(atom[i])
                atom_list.append(atom)
        else:
            for iatom in range(Natom - 1):
                atom = f.readline().strip().split()[1:5]
                for i in range(1, 4):
                    atom[i] = float(atom[i])
                atom_list.append(atom)
    return Molecule(atom_list)

# only support cartisian


def read_gjf(path: str) -> Molecule:

    if os.path.splitext(path)[1] != ".gjf":
        warnings.warn("the file is not gjf file, found {} instead".format(
            os.path.splitext(path)[1]))

    with open(path, mode='r') as f:
        atom_list = []
        while f.readline()[0] != "#":
            pass
        f.readline()  # empty
        title = f.readline()
        f.readline()  # empty
        f.readline()  # net charge and spin
        while True:
            aline = f.readline().strip()
            if aline:
                atom = aline.split()
                for i in range(1, 4):
                    atom[i] = float(atom[i])
                atom_list.append(atom)
            else:
                break
        return Molecule(atom_list)


def write_xyz(mol: Molecule, out_path=sys.stdout):
    print(len(mol), file=out_path)
    print("this file is generated by \"chemsymm\", learn more at https://github.com/Exicweirdo/chemsymm", file=out_path)
    for lis in mol.aslist():
        print("      ".join(map(str, lis)), file=out_path)
    return 0


def write_gjf(mol: Molecule, out_path=sys.stdout):
    print("# SP Test", file=out_path)
    print("", file=out_path)
    print("this file is generated by \"chemsymm\", learn more at https://github.com/Exicweirdo/chemsymm", file=out_path)
    print("", file=out_path)
    print(" 0 1", file=out_path)
    for lis in mol.aslist():
        print("      ".join(map(str, lis)), file=out_path)
    return 0


def symm_path(path: str, log, out_file=None, mode = ".xyz", print_mol=False):
    ext = os.path.splitext(path)[1]
    if ext == ".gjf":
        mol = read_gjf(path)
    else:
        mol = read_xyz(path)

    sym = find_symmetry(mol)
    print(
        f"""
path: {os.path.abspath(path)}
molecular formula: {mol.chem_formula()}
        """,
        file=log
    )
    if print_mol:
        print("coordinates:", file=log)
        print(mol, file=log)
    print("symmetry information:", file=log)
    for key, value in sym.items():
        if key == "standard_base" or key == "m_center":
            print("{}:\n {}".format(key, value), file=log)
        else:
            print("{}: {}".format(key, value), file=log)
    print("---------------------------------------------------", file=log)

    if out_file != None:
        sym_mol = mol.trans(-sym["m_center"])*sym["standard_base"].T
        sym_mol.sort()
        if os.path.exists(os.path.dirname(out_file)) == 0:
            os.makedirs(os.path.dirname(out_file))
        with open(out_file, mode='w') as out:
            if mode == ".xyz":
                write_xyz(sym_mol, out)
            elif mode == ".gjf":
                write_gjf(sym_mol, out)
    return 0


if __name__ == "__main__":
    print(read_xyz("../examples/molecules/a.xyz"))
    print(read_xyz("../examples/molecules/Untitled-1.xyz"))
    print(read_gjf("../examples/molecules/cinnamic acid.gjf"))
    print(write_xyz(read_xyz("../examples/molecules/a.xyz")))
