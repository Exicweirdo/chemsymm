import numpy as np
import numpy.linalg as la
import re
import global_var


#------------------------import periodic table------------------------
elements = {}
elements_s_name = {}
with open("./elements.txt") as ele:
    for line in ele.readlines():
        tup = line.strip()
        Snum = int(re.search(r"\d+", tup).group())
        name = re.search(r"\"\D+\"", tup).group()[1:-1]
        mass = int(re.search(r"\d+ *\}", tup).group()[:-1].strip())
        elements[name] = (Snum, mass)
        elements_s_name[Snum] = name
#--------------------------------------------------------------------- 

class Atomtype:
    def __init__(self, array):
        if isinstance(array[0], str):
            self.atom_num = elements.get(array[0])[0]
            self.element = array[0]
        elif isinstance(array[0], int) or isinstance(array[0], float):
            self.atom_num = int(array[0])
            self.element = elements_s_name.get(array[0])
        self.xyz = np.asarray(array[1:4], dtype = np.float64)
    
    def mass(self):
        attr = elements.get(self.element)

        if attr == None:
            return 0 #dummy atom
        return attr[1]

    def __getitem__(self, index):
        return self.xyz.__getitem__(index)
    
    def __add__(self, array):
        return Atomtype([self.element] + (self.xyz+array).tolist())

    def __mul__(self, array):
        return Atomtype([self.element] + (array@self.xyz).tolist())

    def norm(self):
        return la.norm(self.xyz)
    
    def __eq__(self, others):
        if self.element != others.element:
            return False
        if np.any((self.xyz - others.xyz)**2>global_var.tol**2):
            return False
        else:
            return True
    
    def __lt__(self, others):
        if self.mass() < others.mass():
            return True
        elif self.mass() > others.mass():
            return False
            
        for i in range(2,-1,-1):
            if others.xyz[i] - self.xyz[i] > global_var.tol:
                return True
            elif others.xyz[i] - self.xyz[i] < -global_var.tol:
                return False
        return False
    
    def __gt__(self, others):
        if self.__lt__(others) or self.__eq__(others):
            return False
        else:
            return True

    def to_list(self):
        return [self.element, self.xyz[0], self.xyz[1], self.xyz[2]]

    def xyzm(self):
        return np.append(self.xyz, self.mass())

class Molecule(list):
    def __init__(self, arraylike):
        atoms = []
        for iarray in arraylike:
            atoms.append(Atomtype(iarray))
        super().__init__(atoms)

    def atomtypes(self):
        return [iatom.element for iatom in self]

    def __str__(self):
        output = "\tAtom\tx\t\ty\t\tz\t\t\n"
        for iatom in self:
            output += "\t{}\t{:.4e}\t{:.4e}\t{:.4e}  \n".format(
                iatom.element, iatom.xyz[0], iatom.xyz[1], iatom.xyz[2]
            )
        return(output)

    def aslist(self):
        return [iatom.to_list() for iatom in self]
    
    def asarray(self):
        return np.vstack([atom.xyzm() for atom in self]).transpose()
    def copy(self):
        return Molecule(self.aslist())

    def __add__(self, other):
        if isinstance(other, np.ndarray):
            return Molecule([(iatom + other).to_list() for iatom in self])

        elif isinstance(other, Molecule):
            return Molecule(self.aslist() + other.aslist())
    
    def __mul__(self, b):
        return Molecule([(iatom * b).to_list() for iatom in self])
    
    
    def trans(self, vec: np.ndarray):
        for iatom in self:
            iatom.xyz += vec
        return self
        
    def rotate(self, axis: np.ndarray, angle: float):
        for iatom in self:
            iatom.xyz = np.dot(rot_mat(axis, angle), iatom.xyz)
        return self

    def inverse(self):
        return self * (-1*np.eye(3))
    
    def reflect(self, norm_):
        norm = norm_/la.norm(norm_)
        return self * (np.eye(3) - 2 * np.outer(norm, norm))

#rotation matrix with an axis
def rot_mat(axis, angle):
    ax = axis/la.norm(axis)
    w = np.cos(angle/2)
    xyz = np.sin(angle/2)*ax
    return np.array([
        [1 - 2*xyz[1]**2 - 2*xyz[2]**2, 2*(xyz[0]*xyz[1] - xyz[2]*w), 2*(xyz[0]*xyz[2] + xyz[1]*w)],
        [2*(xyz[0]*xyz[1] + xyz[2]*w), 1 - 2*xyz[0]**2 - 2*xyz[2]**2, 2*(xyz[1]*xyz[2] - xyz[0]*w)],
        [2*(xyz[0]*xyz[2] - xyz[1]*w), 2*(xyz[1]*xyz[2] + xyz[0]*w), 1 - 2*xyz[0]**2 - 2*xyz[1]**2]
    ])

#whether two Molecule is the same
def identical(m1:Molecule, m2:Molecule):
    mol1 = m1.copy()
    mol2 = m2.copy()
    mol1.sort()
    mol2.sort()
    return mol1 == mol2 