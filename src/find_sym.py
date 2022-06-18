import numpy as np
import numpy.linalg as la
import scipy
from numba import jit
import matplotlib.pyplot as plt
import os
import re
#Tolerance <global>

import global_var
from Molecule import Atomtype, Molecule, identical
from moments import *
from ax_deter import *


def find_symmetry(mol_: Molecule)->dict:
    point_group = None
    rank = 0
    # reset center of mass to zero
    m_tot = np.sum(mol_.asarray()[3,:])
    mol = mol_ + -(m_cent(mol_.asarray()))
    # print(m_cent(mol.asarray()))
    has_i = identical(mol, mol.inverse())
    # calculate inertia tensor
    I = np.zeros([3, 3])
    
    inertia(mol.asarray(), I)
    moment1, eigv = la.eigh(I)
    # print(moment)
    #-----------------------------------------------------------------------------------#
    # decide whether the molecule belongs to spherical top / symmetric top /asymmetric top
    if sum(moment1**2) / m_tot**2 < global_var.tol**2:
        rot_type = "point"
        point_group = "SO(3)"
        rank = np.inf
        std_base = np.eye(3)
    elif moment1[0]**2 / m_tot**2 <= global_var.tol**2:
        rot_type = "linear"
        rank = np.inf
        std_base = (eigv.T[::-1, :]).T

    # for T, Td, Th, O, Oh, I and Ih
    elif (moment1[2]-moment1[0]) / m_tot < global_var.tol:
        rot_type = "spherical"
        std_base = np.eye(3)
    # symmetric top for Cn and Dn where n>=3
    # prolate top
    elif (moment1[2]-moment1[1]) / m_tot < global_var.tol:
        rot_type = "symmetric"
        std_base = (eigv.T[::-1, :]).T

    # oblate top
    elif (moment1[1]-moment1[0]) / m_tot < global_var.tol:
        rot_type = "symmetric"
        std_base = eigv

    # asymmetric top
    else:
        rot_type = "asymmetric"
        std_base = (eigv.T[::-1, :]).T

    mol = mol * std_base.T
    #----------------------------------------------------------------------------------#
    if rot_type == 'linear':
        if identical(mol, mol.copy().inverse()):
            point_group = "D_inf_h"
            
        else:
            point_group = "C_inf_v"
    
    elif rot_type == "spherical":
        moment2 = multipole(mol.asarray(), 2)
        # 5 degeneracy of I and Ih group
        if (moment2[1] - moment2[-1]) / m_tot < global_var.tol or (moment2[0] - moment2[-2]) / m_tot < global_var.tol:
            std_base = np.eye(3)
            if has_i:
                point_group = "I_h"
                rank = 120
            else:
                point_group = "I"
                rank = 60
        # T and O group
        else:
            has_octopole = max(G3moment(mol.asarray())**2) / m_tot**2 > global_var.tol**2
            #O group is the only one has octopole without i 

            #calculate if it has C2 or sigmav
            mol_no_0 = Molecule(
                        [iatom.to_list() for iatom in mol if la.norm(iatom.xyz) > global_var.tol])
            mol_no_0.sort()
            _, ele_array, ele_count = np.unique(
                mol_no_0.atomtypes(), return_index=True, return_counts=True
            )
            equiv = []
            for i, idx in enumerate(ele_array):
                equiv += unique(la.norm(mol_no_0.asarray()
                                [:3, :], axis=0)[idx:idx+ele_count[i]])
            min_i = np.argmin(ele_count)
            c2list, sigmalist = find_C2_sig(
                mol_no_0, ele_array[min_i], ele_array[min_i]+ele_count[min_i], z_vert=False)
            #T and Td group
            if has_octopole:
                if sigmalist:
                    point_group = "T_d"
                    rank = 24
                else:
                    point_group = "T"
                    rank = 12
                std_base = np.vstack([
                    c2list[0],
                    c2list[1],
                    np.cross(c2list[0], c2list[1])
                ]).T
            #Oh and Th group
            else:
                print(c2list)
                if len(c2list) == 3:
                    point_group = "T_h"
                    rank = 24
                    std_base = np.vstack([
                        c2list[0],
                        c2list[1],
                        np.cross(c2list[0], c2list[1])
                    ]).T
                else:
                    std_base = cubic_determine(mol, c2list)
                    if has_i:
                        point_group = "O_h"
                        rank = 48
                    else:
                        point_group = "O"
                        rank = 24

    elif rot_type == "symmetric":
        mol.sort()

        mol_no_z = Molecule([iatom.to_list()
                            for iatom in mol if la.norm(iatom.xyz[0:2]) > global_var.tol])
        mol_no_z.sort()
        # find equivalent set of point
        _, ele_array, ele_count = np.unique(
            mol_no_z.atomtypes(), return_index=True, return_counts=True)
        equiv = []
        for i, idx in enumerate(ele_array):
            equiv += unique(mol_no_z.asarray()[2, idx:idx+ele_count[i]])

        # decide the power of main axis
        gcd = np.gcd.reduce(equiv)
        for n in range(gcd, 1, -1):
            if gcd % n == 0:
                if identical(mol.copy().rotate(np.array([0, 0, 1]), 2*np.pi/n), mol):
                    power_m = n
                    break

        has_S2n = identical(mol, mol.reflect(
            np.array([0, 0, 1])).rotate(np.array([0, 0, 1]), np.pi / power_m))
        has_sigmah = identical(mol, mol.reflect(np.array([0, 0, 1])))
        # decide the occurence of σ_v and C_2
        min_i = np.argmin(ele_count)
        c2list, sigmalist = find_C2_sig(
            mol_no_z, ele_array[min_i], ele_array[min_i]+ele_count[min_i], z_vert=True)
        has_C2 = bool(c2list)
        has_sigmav = bool(sigmalist)
        # Dn group and analogous
        if has_C2:
            if has_sigmav:
                if power_m % 2 ^ has_i:
                    point_group = f"D_{power_m}h"
                    rank = 4*power_m
                else:
                    point_group = f"D_{power_m}d"
                    rank = 4*power_m
            else:
                point_group = f"D_{power_m}"
                rank = 2*power_m
        # Cn group and analogous
        else:
            if has_sigmav:
                point_group = f"C_{power_m}v"
                rank = 2*power_m
            elif has_sigmah:
                point_group = f"C_{power_m}h"
                rank = 2*power_m
            elif has_S2n:
                point_group = f"S{2*power_m}"
                rank = 2*power_m
            else:
                point_group = f"C_{power_m}"
                rank = power_m

    elif rot_type == "asymmetric":
        C2_xyz = [identical(mol, mol.copy().rotate(np.eye(3)[:,i], np.pi)) for i in range(3)]
        i = np.argmax(C2_xyz)
        mol = mol * np.eye(3)[[(i-2)%3, (i-1)%3, i%3], :]
        has_sigmah = identical(mol, mol.reflect(np.array([0, 0, 1])))
        if sum(C2_xyz) == 3:
            if has_sigmah:
                point_group = "D_2h"
                rank = 8
            else:
                point_group = "D_2"
                rank = 4
        elif sum(C2_xyz) == 1:
            has_sigmav = identical(mol, mol.reflect(np.array([0, 1, 0])))

            if has_sigmav:
                point_group = "C_2v"
                rank = 4
            elif has_sigmah:
                point_group = "C_2h"
                rank = 4
            else:
                point_group = "C_2"
                rank = 2

        else:
            if has_i:
                point_group = "Ci"
                rank = 2
            elif sum([identical(mol, mol.copy().reflect(np.eye(3)[:,i])) for i in range(3)]) > 0:
                point_group = "Cs"
                rank = 2
            else:
                point_group = "C1"
                rank = 1
    elif rot_type == "point":
        pass
    else:
        raise ValueError("Failed to determine the rotational top, found {} instead.".format(rot_type))

    sym = {
        "point group": point_group,
        "h": rank,
        "rotation_type": rot_type,
        "inertia moment": moment1,
        "i": has_i,
        "standard_base": std_base
    }
    return sym