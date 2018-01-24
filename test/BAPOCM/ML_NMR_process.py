from __future__ import print_function

try:
    import numpy as np
except ImportError:
    print("you need to install numpy")
    sys.exit(1)

np.set_printoptions(threshold='nan')

import math
import copy
import sys
import fileinput
import time

#print("RUNNING NMR_PROCESS.py for", sys.argv[1])

################################################################################
## REFERENCE TABLES
#periodic_table[atom number]
#covalent_radii[atom number][0,1,2]
#bond_order[atom number]
if True:
    periodic_table=[
        '',
        'H',
        'He',
        'Li',
        'Be',
        'B',
        'C',
        'N',
        'O',
        'F',
        'Ne',
        'Na',
        'Mg',
        'Al',
        'Si',
        'P',
        'S',
        'Cl',
        'Ar',
        'K',
        'Ca',
        'Sc',
        'Ti',
        'V',
        'Cr',
        'Mn',
        'Fe',
        'Co',
        'Ni',
        'Cu',
        'Zn',
        'Ga',
        'Ge',
        'As',
        'Se',
        'Br',
        'Kr',
        'Rb',
        'Sr',
        'Y',
        'Zr',
        'Nb',
        'Mo',
        'Tc',
        'Ru',
        'Rh',
        'Pd',
        'Ag',
        'Cd',
        'In',
        'Sn',
        'Sb',
        'Te',
        'I',
        'Xe',
        'Cs',
        'Ba',
        'La',
        'Ce',
        'Pr',
        'Nd',
        'Pm',
        'Sm',
        'Eu',
        'Gd',
        'Tb',
        'Dy',
        'Ho',
        'Er',
        'Tm',
        'Yb',
        'Lu',
        'Hf',
        'Ta',
        'W',
        'Re',
        'Os',
        'Ir',
        'Pt',
        'Au',
        'Hg',
        'Tl',
        'Pb',
        'Bi',
        'Po',
        'At',
        'Rn',
        'Fr',
        'Ra',
        'Ac',
        'Th',
        'Pa',
        'U',
        'Np',
        'Pu',
        'Am',
        'Cm',
        'Bk',
        'Cf',
        'Es',
        'Fm',
        'Md',
        'No',
        'Lr',
        'Rf',
        'Db',
        'Sg',
        'Bh',
        'Hs',
        'Mt',
        'Ds',
        'Rg',
        'Cn',
        'Nh',
        'Fl',
        'Mc',
        'Lv',
        'Ts',
        'Og',
        ]

    covalent_radii = [
        [   0,     0,      0],
        [0.32,	0.00,	0.00],
        [0.46,	0.00,	0.00],
        [1.33,	1.24,	0.00],
        [1.02,	0.90,	0.85],
        [0.85,	0.78,	0.73],
        [0.75,	0.67,	0.60],
        [0.71,	0.60,	0.54],
        [0.63,	0.57,	0.53],
        [0.64,	0.59,	0.53],
        [0.67,	0.96,	0.00],
        [1.55,	1.60,	0.00],
        [1.39,	1.32,	1.27],
        [1.26,	1.13,	1.11],
        [1.16,	1.07,	1.02],
        [1.11,	1.02,	0.94],
        [1.03,	0.94,	0.95],
        [0.99,	0.95,	0.93],
        [0.96,	1.07,	0.96],
        [1.96,	1.93,	0.00],
        [1.71,	1.47,	1.33],
        [1.48,	1.16,	1.14],
        [1.36,	1.17,	1.08],
        [1.34,	1.12,	1.06],
        [1.22,	1.11,	1.03],
        [1.19,	1.05,	1.03],
        [1.16,	1.09,	1.02],
        [1.11,	1.03,	0.96],
        [1.10,	1.01,	1.01],
        [1.12,	1.15,	1.20],
        [1.18,	1.20,	0.00],
        [1.24,	1.17,	1.21],
        [1.21,	1.11,	1.10],
        [1.21,	1.14,	1.06],
        [1.16,	1.07,	1.07],
        [1.14,	1.09,	1.10],
        [1.17,	1.21,	1.08],
        [2.10,	2.02,	0.00],
        [1.85,	1.57,	1.39],
        [1.63,	1.30,	1.24],
        [1.54,	1.27,	1.21],
        [1.47,	1.25,	1.16],
        [1.38,	1.21,	1.13],
        [1.28,	1.20,	1.10],
        [1.25,	1.14,	1.03],
        [1.25,	1.10,	1.06],
        [1.20,	1.17,	1.12],
        [1.28,	1.39,	1.37],
        [1.36,	1.44,	0.00],
        [1.42,	1.36,	1.46],
        [1.40,	1.30,	1.32],
        [1.40,	1.33,	1.27],
        [1.36,	1.28,	1.21],
        [1.33,	1.29,	1.25],
        [1.31,	1.35,	1.22],
        [2.32,	2.09,	0.00],
        [1.96,	1.61,	1.49],
        [1.80,	1.39,	1.39],
        [1.63,	1.37,	1.31],
        [1.76,	1.38,	1.28],
        [1.74,	1.37,	0.00],
        [1.73,	1.35,	0.00],
        [1.72,	1.34,	0.00],
        [1.68,	1.34,	0.00],
        [1.69,	1.35,	1.32],
        [1.68,	1.35,	0.00],
        [1.67,	1.33,	0.00],
        [1.66,	1.33,	0.00],
        [1.65,	1.33,	0.00],
        [1.64,	1.31,	0.00],
        [1.70,	1.29,	0.00],
        [1.62,	1.31,	1.31],
        [1.52,	1.28,	1.22],
        [1.46,	1.26,	1.19],
        [1.37,	1.20,	1.15],
        [1.31,	1.19,	1.10],
        [1.29,	1.16,	1.09],
        [1.22,	1.15,	1.07],
        [1.23,	1.12,	1.10],
        [1.24,	1.21,	1.23],
        [1.33,	1.42,	0.00],
        [1.44,	1.42,	1.50],
        [1.44,	1.35,	1.37],
        [1.51,	1.41,	1.35],
        [1.45,	1.35,	1.29],
        [1.47,	1.38,	1.38],
        [1.42,	1.45,	1.33],
        [2.23,	2.18,	0.00],
        [2.01,	1.73,	1.59],
        [1.86,	1.53,	1.40],
        [1.75,	1.43,	1.36],
        [1.69,	1.38,	1.29],
        [1.70,	1.34,	1.18],
        [1.71,	1.36,	1.16],
        [1.72,	1.35,	0.00],
        [1.66,	1.35,	0.00],
        [1.66,	1.36,	0.00],
        [1.68,	1.39,	0.00],
        [1.68,	1.40,	0.00],
        [1.65,	1.40,	0.00],
        [1.67,	0.00,	0.00],
        [1.73,	1.39,	0.00],
        [1.76,	0.00,	0.00],
        [1.61,	1.41,	0.00],
        [1.57,	1.31,	1.31],
        [1.49,	1.36,	1.26],
        [1.43,	1.21,	1.21],
        [1.41,	1.19,	1.19],
        [1.34,	1.18,	1.18],
        [1.29,	1.13,	1.12],
        [1.28,	1.12,	1.12],
        [1.21,	1.18,	1.18],
        [1.22,	1.30,	1.30],
        [1.36,	0.00,	0.00],
        [1.43,	0.00,	0.00],
        [1.62,	0.00,	0.00],
        [1.75,	0.00,	0.00],
        [1.65,	0.00,	0.00],
        [1.57,	0.00,	0.00],
        ]

    bond_order = [
        50,
        1,
        0,
        1,
        2,
        3,
        4,
        3,
        2,
        1,
        0,
        1,
        2,
        3,
        4,
        3,
        2,
        1,
        0,
        1,
        2,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        1,
        0,
        1,
        2,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        1,
        0,
        1,
        2,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        0,
        1,
        2,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        20,
        0,
        ]
################################################################################


#bond length errors: [aromatic, single, double, triple]
bond_errors =[.03, 0.18, .06, .03]

# get name of xyz file from input args
filename = sys.argv[1]

# define empty list for file input
raw_lines = []

# read file into list line by line
lines = 0
with open(filename, "r") as f:
    for line in f:
        raw_lines.append(line)
        lines += 1

# Define atoms as the number of lines in the file - 2 ( to account for header)
atoms = lines - 2

# Define arrays for output values
xyz_data  = np.zeros((atoms,3),      dtype=np.float64)
dist_mat  = np.zeros((atoms,atoms),  dtype=np.float64)
conn_mat  = np.zeros((atoms,atoms),  dtype=np.float64)
at_number = np.zeros((atoms),        dtype=np.int32)
b_angle   = np.zeros((1000,4),   dtype=np.float64)
d_angle   = np.zeros((1000,5),  dtype=np.float64)
dangle_mat= np.zeros((atoms, atoms), dtype=np.float64)

# Import xyz coordinates and atom number into arrays
for line in raw_lines:
    values = line.split()
    if len(values) == 6:

        try:
            number = int(values[0])
        except:
            continue

        xyz_data[number][0] = values[3]
        xyz_data[number][1]= values[4]
        xyz_data[number][2] = values[5]
        at_number[number] = values[1]

#define connectivity based on distance
for i in range(1,atoms):
    for j in range(i,atoms):

        if i == j:
            conn_mat[i][j] = 0
            conn_mat[j][i] = 0
            continue

        vector = xyz_data[i] - xyz_data[j]
        distance = math.sqrt(pow(vector[0],2) + pow(vector[1],2) + pow(vector[2],2))

        dist_mat[i][j] = distance
        dist_mat[j][i] = distance

        at_1 = at_number[i]
        at_2 = at_number[j]

        cutoff_single = covalent_radii[at_1][0] + covalent_radii[at_2][0] + bond_errors[1]
        cutoff_double = covalent_radii[at_1][1] + covalent_radii[at_2][1] + bond_errors[2]
        cutoff_triple = covalent_radii[at_1][2] + covalent_radii[at_2][2] + bond_errors[3]
        cutoff_ar = 1.4 + bond_errors[0]

        if distance < cutoff_triple:
            conn_mat[i][j] = 3
            conn_mat[j][i] = 3
        elif distance < cutoff_double:
            conn_mat[i][j] = 2
            conn_mat[j][i] = 2
        elif distance < cutoff_ar and at_1 == 6 and at_2 == 6:
            conn_mat[i][j] = 1.5
            conn_mat[j][i] = 1.5
        elif distance < cutoff_single:
            conn_mat[i][j] = 1
            conn_mat[j][i] = 1
        else:
            conn_mat[i][j] = 0
            conn_mat[j][i] = 0

        string = "{0:<10f} {1:<.0f} {2:<.0f} {3:<.2f} {4:<.2f} {5:<.2f} {6:<.2f} {7:<.0f}".format(distance, at_1, at_2, cutoff_single, cutoff_double, cutoff_triple, cutoff_ar,conn_mat[i][j])

        #if distance < cutoff:
            #print(string)

# get total bond order for each atom
for i in range(1, atoms):
    total_order = 0
    for bond in conn_mat[i]:
        total_order += bond

    if total_order != bond_order[at_number[i]]:
        #print("WRONG BOND ORDER, checking for aromatic rings", i, bond_order[i], total_order)
        old_bond = total_order
        ar_check = 0

        if at_number[i] == 6:
            ar_check = 1
            bonds_1 = np.transpose(np.nonzero(conn_mat[i]))
            carbon = 0
            for bond_1 in bonds_1:

                #if conn_mat[i][bond_1[0]] == 2:
                    #print("double bond", at_number[bond_1[0]])

                if at_number[bond_1[0]] == 6:
                    carbon += 1

            #print("C_counter =", carbon)

            if carbon >= 2:
                for b in range(1, atoms):
                    #print(conn_mat[i][b])
                    if conn_mat[i][b] == 2:
                        conn_mat[i][b] = 1.5
                        conn_mat[b][i] = 1.5
                        #print("C: ",conn_mat[i][b])

            total_order = 0
            j = 0
            for bond in conn_mat[i]:
                #print(bond)
                total_order += bond

            #print(total_order)

        if total_order != bond_order[at_number[i]]:
            with open("errors.log", 'a') as f_handle:
                string = "file[{4:<35s}] BOND_ORDER ERROR for atom {0:<2d} atom number {1:<2d} bond_order is {2:<2f}, should be {3:<2f}. . .".format(int(i), int(at_number[i]), total_order, bond_order[at_number[i]], filename)
                print(string, file=f_handle)
                sys.exit(3)
        #else:
            #print("FIXED")



'''
        if at_number[i] == 6:
            connections_1 = np.transpose(np.nonzero(conn_mat[i]))
            for con_1 in connections_1:
                if at_number[con_1[0]] == 6:

                    connections_2 = np.transpose(np.nonzero(conn_mat[con_1[0]]))
                    for con_2 in connections_2:
                            if at_number[con_2[0]] == 6 and con_2[0] != con_1[0]:

                                connections_3 = np.transpose(np.nonzero(conn_mat[con_2[0]]))
                                for con_3 in connections_3:
                                    if at_number[con_3[0]] == 6 and con_3[0] != con_2[0]:

                                        connections_4 = np.transpose(np.nonzero(conn_mat[con_3[0]]))
                                        for con_4 in connections_4:
                                            if at_number[con_4[0]] == 6 and con_4[0] != con_3[0]:

                                                if con_4[0] == i:
                                                    print("RING FOUND")
                                                    conn_mat[con_4[0]][con_3[0]] = 1.5
                                                    conn_mat[con_3[0]][con_2[0]] = 1.5
                                                    conn_mat[con_2[0]][con_1[0]] = 1.5
                                                    conn_mat[con_1[0]][i] = 1.5
                                                else:

                                                    connections_5 = np.transpose(np.nonzero(conn_mat[con_4[0]]))
                                                    for con_5 in connections_5:
                                                        if at_number[con_5[0]] == 6 and con_5[0] != con_4[0]:

                                                            if con_5[0] == i:
                                                                print("RING FOUND")
                                                                conn_mat[con_5[0]][con_4[0]] = 1.5
                                                                conn_mat[con_4[0]][con_3[0]] = 1.5
                                                                conn_mat[con_3[0]][con_2[0]] = 1.5
                                                                conn_mat[con_2[0]][con_1[0]] = 1.5
                                                                conn_mat[con_1[0]][i] = 1.5
                                                            else:

                                                                connections_6 = np.transpose(np.nonzero(conn_mat[con_5[0]]))
                                                                for con_6 in connections_6:
                                                                    if at_number[con_6[0]] == 6 and con_6[0] != con_5[0]:

                                                                        if con_6[0] == i:
                                                                            print("RING FOUND")
                                                                            conn_mat[con_6[0]][con_5[0]] = 1.5
                                                                            conn_mat[con_5[0]][con_4[0]] = 1.5
                                                                            conn_mat[con_4[0]][con_3[0]] = 1.5
                                                                            conn_mat[con_3[0]][con_2[0]] = 1.5
                                                                            conn_mat[con_2[0]][con_1[0]] = 1.5
                                                                            conn_mat[con_1[0]][i] = 1.5
                                                                        else:

                                                                            connections_7 = np.transpose(np.nonzero(conn_mat[con_6[0]]))
                                                                            for con_7 in connections_7:
                                                                                if at_number[con_7[0]] == 6 and con_7[0] != con_6[0]:

                                                                                    if con_7[0] == i:
                                                                                        print("RING FOUND")
                                                                                        conn_mat[con_7[0]][con_7[0]] = 1.5
                                                                                        conn_mat[con_6[0]][con_5[0]] = 1.5
                                                                                        conn_mat[con_5[0]][con_4[0]] = 1.5
                                                                                        conn_mat[con_4[0]][con_3[0]] = 1.5
                                                                                        conn_mat[con_3[0]][con_2[0]] = 1.5
                                                                                        conn_mat[con_2[0]][con_1[0]] = 1.5
                                                                                        conn_mat[con_1[0]][i] = 1.5

        total_order = 0
        for bond in conn_mat[i]:
            total_order += bond

            if total_order != bond_order[at_number[i]]:
                print("WRONG BOND ORDER, binning structure", i, bond_order[i], total_order)

                sys.exit(3)
'''

## WORK OUT BOND ANGLES
for i in range(1, atoms):
    count = 0
    bonds_1 = np.transpose(np.nonzero(conn_mat[i]))
    for bond_1 in bonds_1:
        bonds_2 = np.transpose(np.nonzero(conn_mat[bond_1[0]]))
        for bond_2 in bonds_2:
            if bond_2[0] == i:
                continue
            vector1 = xyz_data[i] - xyz_data[bond_1[0]]
            vector2 = xyz_data[bond_1[0]] - xyz_data[bond_2[0]]

            low = np.linalg.norm(vector1) * np.linalg.norm(vector2)
            high = np.inner(vector1, vector2)

            try:
                angle = math.degrees(math.acos(high / low))
            except:
                print(high, low, vector1, vector2)

            b_angle[count] = [i, bond_1[0], bond_2[0], angle]
            count + 1


count = 0

## WORK OUT dihedral ANGLES
for i in range(1, atoms):
    bonds_1 = np.transpose(np.nonzero(conn_mat[i]))
    for bond_1 in bonds_1:
        if bond_1[0] == i:
            continue
        bonds_2 = np.transpose(np.nonzero(conn_mat[bond_1[0]]))
        for bond_2 in bonds_2:
            if bond_2[0] == i or bond_2[0] == bond_1[0]:
                continue
            bonds_3 = np.transpose(np.nonzero(conn_mat[bond_2[0]]))
            for bond_3 in bonds_3:
                if bond_3[0] == i or bond_3[0] == bond_2[0] or bond_3[0] == bond_1[0]:
                    continue

                vector1 = xyz_data[i] - xyz_data[bond_1[0]]
                vector2 = xyz_data[bond_1[0]] - xyz_data[bond_2[0]]
                vector3 = xyz_data[bond_2[0]] - xyz_data[bond_3[0]]

                norma_1 = vector1/np.linalg.norm(vector1)
                norma_2 = vector2/np.linalg.norm(vector2)
                norma_3 = vector3/np.linalg.norm(vector3)

                perp_1 = np.cross(norma_1, norma_2)
                perp_2 = np.cross(norma_2, norma_3)

                m_1 = np.cross(perp_1, norma_2)

                dot_1 = np.inner(perp_1, perp_2)
                dot_2 = np.inner(m_1, perp_2)

                try:
                    angle = math.degrees(np.arctan2(dot_1,dot_2)) - 90
                except:
                    print(dot_1, dot_2)

                if angle < -180:
                    angle = angle + 360

                dangle_mat[i][bond_3[0]] = -angle
                dangle_mat[bond_3[0]][i] = -angle

                d_angle[count] = [i, bond_1[0], bond_2[0], bond_3[0], -angle]
                count += 1


##### output
name = filename.split(".")[0]

dist_file   =   name + "_dist.txt"
conn_file   =   name + "_conn.txt"
bangle_file =   name + "_bangle.txt"
dangle_file =   name + "_dangle.txt"
atno_file   =   name + "_atno.txt"

with open(dist_file, 'w') as f_handle:
    for i in range(0, atoms):
        for j in range(0, atoms):
            string = "{0:<15.4f},  ".format(dist_mat[i][j])
            print(string, file=f_handle, end='')
        print('', file=f_handle)

with open(conn_file, 'w') as f_handle:
    for i in range(0, atoms):
        for j in range(0, atoms):
            string = "{0:<5.0f},  ".format(conn_mat[i][j])
            print(string, file=f_handle, end='')
        print('', file=f_handle)

with open(bangle_file, 'w') as f_handle:
    for i in range(0, 1000):
        if not b_angle[i].any:
            continue
        for k in range(0, 4):
            string = "{0:<15.4f},  ".format(b_angle[i][k])
            print(string, file=f_handle, end='')
        print(" ", file=f_handle)
    print(i, file=f_handle)

with open(dangle_file, 'w') as f_handle:
    for j in range(0, 1000):
        if np.all(d_angle[j]==0):
            continue
        string = "{0:<5.0f}, {1:<5.0f}, {2:<5.0f}, {3:5.0f}, {4:15.4f},  ".format(d_angle[j][0], d_angle[j][1], d_angle[j][2], d_angle[j][3], d_angle[j][4])
        print(string, file=f_handle)

with open(atno_file, 'w') as f_handle:
    for i in range(0, atoms):
        print(at_number[i], file=f_handle, end=',')
        print(" ", file=f_handle)


























#


sys.exit(0)












        #
