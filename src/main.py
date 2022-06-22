import time
import os
import sys

from numpy import isin
from sympy import arg
from find_sym import find_symmetry
from read_file import read_gjf, read_xyz, symm_path
import global_var
from argparse import ArgumentParser
import glob

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
def parse_arg():
    parser = ArgumentParser("find symmetry of molecules and move to standard basis")
    parser.add_argument("-f", "--fin", nargs="+", type=str, required=True,
        help="input files to be processed, support .xyz and .gjf file."
    )
    parser.add_argument("-tol", "--tolerance", action="store", 
        type=float, nargs=1, default=0.01,
        help="within which two molecule is considered the same, larger value gives higher symmetry"
    )
    parser.add_argument("-l", "--log", type=str, default=None,
     const="./sym_info_{}.log".format(time.strftime("%Y%m%d_%H_%M", time.localtime())),
        nargs="?",
        help="the file output informations will be write to, the default is sys.stdout"
    )
    parser.add_argument("-o", "--output", type=str, default=None, 
        nargs="?", const="*_sym.gjf",
        help="output symmetrized files, should correspond with input files, default is .gjf"
    )
    parser.add_argument("-p", "--printmolecule", action="store_true",
        help='whether to print molecule coordinates in log, default: False'
    )
    return parser
if __name__ == "__main__":

    parser = parse_arg()
    args = parser.parse_args() #"-l -f ../examples/molecules/*.xyz ../examples/molecules/*.gjf".split()
    global_var.set_tol(args.tolerance)
    infiles = []
    for path in args.fin:
        infiles += glob.glob(path)
    
    if args.log:
        log = open(args.log, mode="a")
        print("["+"-"*50+"]", end='\r')
    else:
        log = sys.stdout

    if args.output:
        if isinstance(args.output, str):
            outputlist = [args.output.replace("*", os.path.splitext(path)[0]) for path in infiles]
        else:
            outputlist = args.output
    else:
        outputlist = [None * len(infiles)]

    if not infiles:
        print("there's no file to process", file=log)

    else:
        for n, file in enumerate(infiles):
            if args.log:
                percent = n / len(infiles)
                filename = os.path.basename(file)
                print(" "*100, end='\r')
                print("["+"*"*int((50*percent)//1)+"-"*int(50-(50*percent)//1)+"]   processing:{}".format(filename), end='\r')
            
            symm_path(file, log, outputlist[n], print_mol=args.printmolecule, mode=".gjf")
        if args.log:
            print(f"["+"*"*50+"]")  
    

    '''tup = argv_process(argv=sys.argv)
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
path: {os.path.abspath(path)}
molecular formula: {mol.chem_formula()}
coordinates:
                """,
                file=out_file 
            )
            print(mol, file=out_file)
            print("symmetry information:", file=out_file)
            for key, value in sym.items():
                print("{}: {}".format(key, value), file=out_file)
            print("---------------------------------------------------")'''