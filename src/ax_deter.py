from Molecule import *
import numpy as np
import numpy.linalg as la
import global_var


# finding C2 axis and mirror in an atom set
# Parameters:
# ----------------
# mol:Molecule
#   the molecule to operate
# start:int
#   the start index of mol
# end:int
#   the end index of mol
# z_vert:bool
#   if True, only consider those C2 who is vertical to z-axis and those mirror whose normal is vertical
# 
# Returns:
# ----------------
# C2_list: list[ndarray]
#   list of C2 axis
# sigma_list: list[ndarray]
#   list of normal vectors of mirrors
def find_C2_sig(mol:Molecule, start:int, end:int, z_vert:bool = False):

    C2_list = []
    sigma_list = []
    for i in range(start, end):
        for j in range(i+1,end):
            iatom, jatom = mol[i], mol[j]
            #distance to zero point is equal
            if (la.norm(iatom.xyz)-la.norm(jatom.xyz))**2 < global_var.tol**2:
                sigma_v_norm = (iatom.xyz - jatom.xyz)/2
                C2_ax = (iatom.xyz + jatom.xyz)/2
                #C2
                #if z_vert, check if axis vertical to z axis or skip
                if (C2_ax[2]**2 < global_var.tol**2 or not z_vert):
                    if la.norm(C2_ax) <= global_var.tol:
                        C2_ax = sigma_v_norm
                    
                    if identical(mol, mol.copy().rotate(C2_ax, np.pi)):
                        C2_ax = C2_ax/la.norm(C2_ax)
                        is_new_ax = True
                        for ax in C2_list:
                            if la.norm(C2_ax-ax) < global_var.tol or la.norm(C2_ax+ax) < global_var.tol:
                                is_new_ax = False
                        if is_new_ax:
                            C2_list.append(C2_ax)
                #σ_v

                #if z_vert, check if axis vertical to z axis or skip
                if (sigma_v_norm[2]**2 < global_var.tol**2 or not z_vert) and la.norm(sigma_v_norm) > global_var.tol:
                    if identical(mol, mol.reflect(sigma_v_norm)):
                        sigma_v_norm = sigma_v_norm/la.norm(sigma_v_norm)
                        is_new_sigma = True
                        for norm in sigma_list:
                            if la.norm(sigma_v_norm - norm) < global_var.tol or la.norm(sigma_v_norm + norm) < global_var.tol:
                                is_new_sigma = False
                        if is_new_sigma:
                            sigma_list.append(sigma_v_norm)
    return C2_list, sigma_list

# find a arrays unique element wither tolerance
# return: list of counts of unique element
def unique(sorted_array):
    eqiv = []
    z0 = sorted_array[0]
    eq_count = 1
    for z in sorted_array[1:]:
        if (z - z0)**2 > global_var.tol**2:
            eqiv.append(eq_count)
            z0 = z
            eq_count = 1
        else:
            eq_count += 1
    eqiv.append(eq_count)
    return eqiv

def cubic_determine(mol:Molecule, C2_list):
    ax_list = []
    for ax in C2_list:
        if identical(mol, mol.copy().rotate(ax, np.pi/2)):
            ax_list.append(ax)
    
    if ax_list:
        return np.vstack([
            ax_list[0],
            ax_list[1],
            np.cross(ax_list[0], ax_list[1])
        ]).T
    else:
        ax1 = C2_list[0]
        for ax2 in C2_list[1:]:
            if np.dot(ax1, ax2)**2 < global_var.tol**2:
                b1 = (ax1 + ax2)/2 / (la.norm((ax1 + ax2)/2))
                b2 = (ax1 - ax2)/2 / (la.norm((ax1 - ax2)/2))
                return np.vstack([
                    b1,
                    b2,
                    np.cross(b1, b2)
                ]).T