import os
import re
from Molecule import *
import warnings

def read_xyz(path:str) -> Molecule:

    if os.path.splitext(path)[1] != ".xyz":
        warnings.warn("the file is not xyz file, found {} instead".format(os.path.splitext(path)[1]))
    
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
                for i in range(1,4):
                    atom[i] = float(atom[i])
                atom_list.append(atom)
        else:
            for iatom in range(Natom - 1):
                atom = f.readline().strip().split()[1:5]
                for i in range(1,4):
                    atom[i] = float(atom[i])
                atom_list.append(atom)
    return Molecule(atom_list)

#only support cartisian
def read_gjf(path:str) -> Molecule:

    if os.path.splitext(path)[1] != ".gjf":
        warnings.warn("the file is not gjf file, found {} instead".format(os.path.splitext(path)[1]))
    
    with open(path, mode='r') as f:
        atom_list = []
        while f.readline()[0] != "#":
            pass
        f.readline() #empty
        title = f.readline()
        f.readline() #empty
        f.readline() #net charge and spin
        while True:
            aline = f.readline().strip()
            if aline:
                atom = aline.split()
                for i in range(1,4):
                    atom[i] = float(atom[i])
                atom_list.append(atom)
            else:
                break
        return Molecule(atom_list)

if __name__ == "__main__":
    print(read_xyz("./a.xyz"))
    print(read_xyz("./Untitled-1.xyz"))
    print(read_gjf("./cinnamic acid.gjf"))
    
