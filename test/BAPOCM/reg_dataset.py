### reg_dataset.py
#
#   Will Gerrard (wg12385@bristol.ac.uk)
#
#   Script to convert seperate J coupling and dihedral angle files into a single column formatted CSV file
#

#import statements
from __future__ import print_function
import numpy as np
import sys


molname = sys.argv[1]
D_file = sys.argv[2]
J_file = sys.argv[3]
N_file = sys.argv[4]
Out_file = "D_J.Data"

atoms = 200


J_matrix = np.zeros((atoms,atoms), dtype=np.float64)
D_matrix = np.zeros((1000, 5), dtype=np.float64)
N_matrix = np.zeros((atoms), dtype=np.float64)


i=0
j=0
J_lines = []
with open(J_file, "r") as f:
    for line in f:
        J_lines.append(line)
        items = [x.strip() for x in line.split(',')]
        j = 0
        for val in items:
            if len(val) == 0:
                continue
            else:
                J_matrix[i][j] = float(val)
                j += 1
        i += 1

D_lines = []


i = 0
with open(D_file, "r") as f:
    for line in f:
        D_lines.append(line)
        items = [x.strip() for x in line.split(',')]
        D_matrix[i][0] = items[0]
        D_matrix[i][1] = items[1]
        D_matrix[i][2] = items[2]
        D_matrix[i][3] = items[3]
        D_matrix[i][4] = items[4]
        i += 1



i=0
j=0
N_lines = []
with open(N_file, "r") as f:
    for line in f:
        N_lines.append(line)
        items = [x.strip() for x in line.split(',')]
        j = 0
        for val in items:
            if len(val) == 0:
                continue
            else:
                N_matrix[i] = float(val)
                j += 1
            if j > 1:
                print("more than one value per line, something went wrong")
        i += 1


with open(Out_file, "a") as f:
    for i in range(0, 1000):
        D = D_matrix[i][4]
        X = D_matrix[i][0]
        Y = D_matrix[i][3]
        a = D_matrix[i][1]
        b = D_matrix[i][2]
        J = J_matrix[X-1][Y-1]
        Xn = N_matrix[X]
        Yn = N_matrix[Y]
        an = N_matrix[a]
        bn = N_matrix[b]

        if J != 0 and D != 0 and Xn == 1 and Yn == 1:
            string = "{0:<12.6f},   {1:<12.6f},   {2:<12.0f},   {3:<12.0f},   {4:<12.0f},   {5:<12.0f},   {6:<12.0f},   {7:<12.0f},   {8:<12.0f},   {9:<12.0f},   {10:<12s},".format(D, J, X, a, b, Y, Xn, an, bn, Yn, molname)
            print(string, file=f)























#
