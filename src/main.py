import find_sym
from Molecule import Atomtype, Molecule, identical
from find_sym import find_symmetry
from read_file import read_gjf, read_xyz
import global_var
if __name__ == "__main__":
    global_var.set_tol(0.01)    
    sym1 = find_symmetry(read_gjf("./BH3.gjf"))
    for key, value in sym1.items():
        print("{}: {}".format(key, value))