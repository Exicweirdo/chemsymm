import os
import sys
from find_sym import find_symmetry
from read_file import read_gjf, read_xyz
import global_var
def argv_process(argv:list):
    tol = 0.01
    input_list = []
    input_flag = False
    out_file = sys.stdout
    for i, str_i in enumerate(argv):
        if input_flag == True:
            input_list.append(str_i)

        if str_i[0] == "-":
            input_flag = False
            if "-tol" == str_i:
                tol = float(argv[i+1])
            elif "-o" == str_i:
                out_file = argv[i+1]
            elif "-f" == str_i:
                input_flag = True
            elif "-h" == str_i:
                print("""
Usage: main.py -f <path(s)> -tol <tolerance number in angstrom> -o <output file>
                
path: path of .xyz file or .gjf file, default is xyz\n
tolerance: within which symmetry is checked if atoms have the same coordinate\n
output file: to write outputs of program, if not set, output will be print on command line
                """)
                return 0
            else:
                raise ValueError("Unexpected flags {}".format(str_i))

    return (tol, input_list, out_file)

if __name__ == "__main__":
    tup = argv_process(argv=sys.argv)
    if tup:
        global_var.set_tol(tup[0])
        input_list = tup[1]
        out_file = tup[2]
        if not input_list:
            print("there's no file to process")
            argv_process(["-h"])
        for path in input_list:
            if not os.path.isfile(path):
                raise IOError(f"failed to open {path}: no such file")
        for path in input_list:
            ext = os.path.splitext(path)[1]
            if ext == ".gjf":
                mol = read_gjf(path)
            else:
                mol = read_xyz(path)
            
            sym = find_symmetry(mol)
            print(
                f"""
path: {path}
molecular formula: {mol.chem_formula()}
coordinates:
                """,
                file=out_file 
            )
            print(mol, file=out_file)
            print("symmetry information:", file=out_file)
            for key, value in sym.items():
                print("{}: {}".format(key, value), file=out_file)
            print("---------------------------------------------------")