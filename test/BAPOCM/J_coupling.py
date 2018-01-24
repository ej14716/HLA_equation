### J_coupling.py
#
#   Will Gerrard (wg12385@bristol.ac.uk)
#
#   Gets coupling constant data from NMR log files generated from gaussian NMR calculations run using NMR_setup.py
#
#   Based on awk and matlab scripts (Claire Dickson)
#

#import statements
from __future__ import print_function
import numpy as np
import sys

# Get molecule name from input argument
molname = sys.argv[1]

# Initialise input list, conformer number and atom number
raw_lines = []
conformers = 0
atoms = 0

# Define input summary file and output files
J_file = molname + "_J_Summary.txt"
out_file2 = molname + "_J_Raw.txt"

# Read through file to get number of atoms and number of conformers
with open(J_file, "r") as f:
    for line in f:
        raw_lines.append(line)
        items = line.split()
        if len(items) > 0:
            if items[0].isdigit():
                if int(items[0]) > atoms:
                    atoms = int(items[0])

# Define J constant array and array for boltzmann averaged values
j_matrix = np.zeros((atoms+1,atoms+1), dtype=np.float64)

# Read data from file into array
for line in raw_lines:
    items = line.split()
    if len(items) > 0:
        if len(items) < 2:
            continue
        if items[1].isdigit() and items[0].isdigit():
            try:
                a = int(items[0])-1
            except:
                a = "Nan"
            try:
                b = int(items[1])-1
            except:
                b = "NaN"
            try:
                c = int(items[2])-1
            except:
                c = "NaN"
            try:
                d = int(items[3])-1
            except:
                d = "NaN"
            try:
                e = int(items[4])-1
            except:
                e = "NaN"

        if items[0].isdigit() and not items[1].isdigit():
            atom = int(items[0])-1
            if a == "NaN" or len(items) < 3:
                continue
            j_matrix[atom][a] = float(items[1])
            if b == "NaN" or len(items) < 3:
                continue
            j_matrix[atom][b] = float(items[2])
            if c == "NaN" or len(items) < 4:
                continue
            j_matrix[atom][c] = float(items[3])
            if d == "NaN" or len(items) < 5:
                continue
            j_matrix[atom][d] = float(items[4])
            if e == "NaN" or len(items) < 6:
                continue
            j_matrix[atom][e] = float(items[5])

with open(out_file2, "w") as f:
        for x in range(0, atoms+1):
            for y in range(0, atoms+1):
                string = "{0:<12.6f},   ".format(j_matrix[x][y])
                print(string, end="", file=f)
            print("", file=f)

sys.exit(3)
