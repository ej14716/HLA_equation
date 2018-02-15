# ej14716
#
#
#
#


from __future__ import print_function

import math

import numpy as np

import sys

import glob

np.set_printoptions(threshold = 5000)

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass

    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass

    return False

raw_names = glob.glob("test/*_J_Raw.txt")
names = []
for items in raw_names:
	item = items.split('/')
	it = item[1].split('_')
	names.append(it[0])

for name in names:
    string1 = "test/" + name + "_dangle.txt"
    file1 = glob.glob(string1)
    file1 = file1[0]
    string2 = "test/" + name + "_atno.txt"
    file2 = glob.glob(string2)
    file2 = file2[0]
    string3 = "test/" + name + "_conn.txt"
    file3 = glob.glob(string3)
    file3 = file3[0]
    string4 = "test/" + name + "_J_Raw.txt"
    file4 = glob.glob(string4)
    file4 = file4[0]


    p1 = 13.70

    p2 = -0.73

    p3 = 0

    p4 = 0.56

    p5 = -2.47

    p6 = 16.9

    p7 = 0.14



    atom_1 = []

    atom_2 = []

    atom_3 = []

    atom_4 = []

    dihedral_angle1 = []

    atomic_number = []

    array_size = 0








    electronegativity_array = [

    0,

    0,

    50,

    50,

    50,

    50,

    0.4,

    0.85,

    1.3,

    1.7,

    50,

    50,

    50,

    50,

    -0.3,

    -0.05,

    0.4,

    0.95,

    50,

    50,

    50,

    50,

    50,

    50,

    50,

    50,

    50,

    50,

    50,

    50,

    50,

    50,

    -0.3,

    -0.1,

    0.35,

    0.75,

    50,

    50,

    50,

    50,

    50,

    50,

    50,

    50,

    50,

    50,

    50,

    50,

    50,

    50,

    -0.3,

    -0.15,

    0.1,

    0.45

    ]





    #open angles file and add them to angles list

    with open(file1, "r") as f:

    	for line in f:

    		items = line.split(",")

    		atom_1.append(int(items[0]))

    		atom_2.append(int(items[1]))

    		atom_3.append(int(items[2]))

    		atom_4.append(int(items[3]))

    		dihedral_angle1.append(float(items[4]))

    		array_size += 1



    with open(file2, "r") as f:

    	for line in f:

    		items = line.split(",")

    		atomic_number.append(int(items[0]))


    connectivityx = []




    with open(file3, "r") as f:

    	for line in f:

    		items = line.split(",")
    		connectivityx.append(items)


    connectivity_array = np.zeros((len(connectivityx), len(connectivityx[0])), dtype=np.float64)

    y=0
    x=0
    for line in connectivityx:
    	y = 0
    	for item in line:
    		if not is_number(item):
    			continue
    		connectivity_array[x][y] = item
    		y += 1
    	x += 1

    DFTJ = []

    with open(file4, "r") as f:

    	for line in f:

    		items = line.split(",")
    		DFTJ.append(items)
    y = 0

    x = 0
    
    DFTJ_array = np.zeros((len(DFTJ), len(DFTJ)), dtype=np.float64)
    
    for line in DFTJ:
    	y = 0
    	for item in line:
    		if not is_number(item):
    			continue
    		DFTJ_array[x][y] = item
    		y += 1
    	x += 1
    	

    dihedral_array = np.zeros((array_size, 8), dtype = np.float64)



    for i in range(array_size):

    	dihedral_array[i][0] = 0

    	dihedral_array[i][1] = atom_1[i]

    	dihedral_array[i][2] = atom_2[i]

    	dihedral_array[i][3] = atom_3[i]

    	dihedral_array[i][4] = atom_4[i]

    	dihedral_array[i][5] = dihedral_angle1[i]
     	dihedral_array[i][6] = i

    	if atom_1[i] > atom_4[i]:
    		low = atom_4[i]
    		high = atom_1[i]
    	else:
    		low = atom_1[i]
    		high = atom_4[i]


    	dihedral_array[i][7] = DFTJ_array[int(high)][int(low)-1]


    	if atomic_number[atom_1[i]] == 1:

    		if atomic_number[atom_2[i]] == 6:

    			if atomic_number[atom_3[i]] == 6:

    				if atomic_number[atom_4[i]] == 1:
				
					if atom_4[i] < atom_1[i]:
					
   						dihedral_array[i][0] = 1
    a = 0
    ivalues_failed = []
    bond_angle1 = np.zeros((1, 2), dtype = np.float64)
    bond_angle2 = np.zeros((1, 2), dtype = np.float64)

    sub_array = np.zeros((array_size, 5), dtype = np.int64)
    Esub1B = []
    Esub2B = []
    Esub3B = []
    Esub4B = []
    electronegativity_Barray = np.zeros((array_size, 4), dtype = np.float64)

    electronegativity_Aarray = np.zeros((array_size, 4), dtype = np.float64)

    A_array = np.zeros((array_size, 4), dtype = np.float64)

    B_array = np.zeros((array_size, 4), dtype = np.float64)

    C_array = np.zeros((array_size, 4), dtype = np.float64)

    D_array = np.zeros((array_size, 4), dtype = np.float64)

    E_array = np.zeros((array_size, 4), dtype = np.float64)

    F_array = np.zeros((array_size, 1), dtype = np.float64)

    cos_array = np.zeros((array_size, 1), dtype = np.float64)

    cossqr_array = np.zeros((array_size, 1), dtype = np.float64)

    Jvalues_array = np.zeros((array_size, 1), dtype = np.float64)
    Jvalues_array2 = np.zeros((array_size, 1), dtype = np.float64)
    HLA_array = np.zeros((array_size, 5), dtype = np.float64)
		
    #get atom labels for substituents and their relative orientations
    for i in range(array_size):
    	sub_array[i][4] = i
	SP3 = 1
    #only look at pathways  that are HCCH
    	if dihedral_array[i][0] != 0:
		bonds_8  = np.transpose(np.nonzero(connectivity_array[dihedral_array[i][4]]))
		if len(bonds_8) != 1:
			print("H2connectivity_error", name, i)
			dihedral_array[i][0] = 9
			continue
    		bonds_1 = np.transpose(np.nonzero(connectivity_array[dihedral_array[i][1]]))

    		if not len(bonds_1) == 1:
    			print("Hconnectivity_error", name, i)
    			ivalues_failed.append(i)
			dihedral_array[i][0] = 9
    			continue

    		#check H is connected to C

    		if bonds_1[0] == dihedral_array [i][2]:

    			connection1 = 1

    		else:

    			connection1 = 0

    		#if it is connected to C search for its connections

    		if connection1 == 1:

    			bonds_2 = np.transpose(np.nonzero(connectivity_array[dihedral_array[i][2]]))
    			if not len(bonds_2) == 4:
    			
				dihedral_array[i][0] = 7
    				ivalues_failed.append(i)
    				continue
    			#search each of the connections
			
    			for j in range (len(bonds_2)):
				if SP3 == 0:
					continue
    				if bonds_2[j] == dihedral_array[i][1]:

    					continue

    					#if the connection is to the third atom (C in coupling) search for what the C is connected for

    				elif bonds_2[j] == dihedral_array[i][3]:

    					bonds_3 = np.transpose(np.nonzero(connectivity_array[dihedral_array[i][3]]))

    					#loop through each of the bonds that C is connected to

    					if not len(bonds_3) == 4:
						
    						SP3 = 0
						dihedral_array[i][0] = 7
    						continue
    					for x in range (len(bonds_3)):

    						if bonds_3[x] == dihedral_array[i][2]:

    							continue

    						elif bonds_3[x] == dihedral_array[i][4]:

    							continue

    					#assign the substituents based on their orientation

    						else:

    							for k in range(array_size):
    								b = 0
    							#if the first column in dihedral array is equal to the connection we've just found and the fourth one is the 1st atom (H)in the origional coupling

    								if (dihedral_array[k][1] == bonds_3[x]) and (dihedral_array[k][4] == dihedral_array[i][1]):

    									bond_angle1[b][0] = (dihedral_array[k][5])

    							#have to take the negative of the original dihedral because we want the couplig from H4 to H1

    									bond_angle1[b][1]  = dihedral_array[i][5] - bond_angle1[b][0]





    									if (bond_angle1[b][1] > 0):

    										if (abs(bond_angle1[b][1]) < 180):

    											sub_array[i][2] = (bonds_3[x])

    										else:

    											sub_array[i][3] = (bonds_3[x])

    									if (bond_angle1[b][1]< 0):

    										if (abs(bond_angle1[b][1]) > 180) and (abs(bond_angle1[b][1]) < 360):
    											sub_array[i][2] = (bonds_3[x])
    										else:

    											sub_array[i][3] = (bonds_3[x])
									

						if SP3 == 0:
							sub_array[i][2] = 0	
							sub_array[i][3] = 0

				
    				else:

					if SP3 == 0:
						continue
						
    				#this is now for other bonds conncted to first carbon again

    					for Y in range (array_size):
    						c = 0
    					#find dihedral angle for the substituents to the fourth atom (H)

    						if (dihedral_array[Y][1] == bonds_2[j]) and (dihedral_array[Y][4] == dihedral_array[i][4]):

    							bond_angle2[c][0] = dihedral_array[Y][5]

    							bond_angle2[c][1] = dihedral_array[i][5] - bond_angle2[c][0]


    							if (bond_angle2[c][1] > 0):

    								if abs(bond_angle2[c][1]) < 180:

    									sub_array[i][0] = (bonds_2[j])

    								else:

    									sub_array[i][1] = (bonds_2[j])

    							if (bond_angle2[c][1] < 0):

    								if (abs(bond_angle2[c][1]) > 180) and (abs(bond_angle2[c][1]) < 360):

    									sub_array[i][0] = (bonds_2[j])

    								else:

    									sub_array[i][1] = (bonds_2[j])

		if SP3 == 0:
			continue
    #now need to get electronegativity values for these substituents and the beta substituents
       
    for a in range(array_size):
 	if sub_array[a][0] == 0:
		continue
       

	bonds_4 = np.transpose(np.nonzero(connectivity_array[sub_array[a][0]]))
		
   	Esub1B_1 = 0

    	for z in range(len(bonds_4)):

    		if bonds_4[z] == dihedral_array[a][2]:

    			continue

    		else:

    			#assign and sum their electronegativities

    			Esub1B_1 += electronegativity_array[atomic_number[bonds_4[z]]]


    	if sub_array[a][1] == 0:
    		continue

    	bonds_5 = np.transpose(np.nonzero(connectivity_array[sub_array[a][1]]))


    	Esub2B_1 = 0


    	for z in range(len(bonds_5)):

    		if bonds_5[z] == dihedral_array[a][2]:

    			continue

    		else:

    			#assign and sum their electronegativities

    			Esub2B_1 += electronegativity_array[atomic_number[bonds_5[z]]]

    	if sub_array[a][2] == 0:
    		continue
    	bonds_6 = np.transpose(np.nonzero(connectivity_array[sub_array[a][2]]))
	
    		
    	Esub3B_1 = 0

    	for z in range(len(bonds_6)):

    		if bonds_6[z] == dihedral_array[a][3]:

    			continue

    		else:

    			#assign and sum their electronegativities

    			Esub3B_1 += electronegativity_array[atomic_number[bonds_6[z]]]
    	if sub_array[a][3] == 0:
    		continue

    	bonds_7 = np.transpose(np.nonzero(connectivity_array[sub_array[a][3]]))
	    		
	Esub4B_1 = 0

	for z in range(len(bonds_7)):
		if bonds_7[z] == dihedral_array[a][3]:
			continue

    		else:

    			#assign and sum their electronegativities

    			Esub4B_1 += electronegativity_array[atomic_number[bonds_7[z]]]



    	#array containing the sum of the electronegativities of the Beta substituents of each alpha substituent for every value of i



    	electronegativity_Barray[a][0] = Esub1B_1

    	electronegativity_Barray[a][1] = Esub2B_1

    	electronegativity_Barray[a][2] = Esub3B_1

    	electronegativity_Barray[a][3] = Esub4B_1
    

    #electronegativity of the alpha substituents array

    for i in range(array_size):

    	electronegativity_Aarray[i][0] = electronegativity_array[atomic_number[sub_array[i][0]]]

    	electronegativity_Aarray[i][1] = electronegativity_array[atomic_number[sub_array[i][1]]]

    	electronegativity_Aarray[i][2] = electronegativity_array[atomic_number[sub_array[i][2]]]

    	electronegativity_Aarray[i][3] = electronegativity_array[atomic_number[sub_array[i][3]]]


    # x the sum of the electronegativity vales of the beta substituents by p7

    A_array	= p7 * (electronegativity_Barray)



    B_array = electronegativity_Aarray - A_array





    #do stuff with the positive and negtive substituents


    for i in range(array_size):

    	C_array[i][0] = np.radians(dihedral_array[i][5]) + (p6 * (np.absolute(B_array[i][0])))

    	C_array[i][1] = np.negative(np.radians(dihedral_array[i][5])) + (p6 * (np.absolute(B_array[i][1])))

    	C_array[i][2] = np.radians(dihedral_array[i][5]) + (p6 * (np.absolute(B_array[i][2])))

    	C_array[i][3] = np.negative(np.radians(dihedral_array[i][5])) + (p6 * (np.absolute(B_array[i][3])))


    # make aray that needs to be * by the electronegativity (B array)

    D_array = p4 + (p5 * np.square(np.cos(C_array)))


    E_array = D_array * B_array


    for k in range(array_size):
    	sum = 0
    	for l in range(4):
    		sum = sum + E_array[k][l]
    	F_array[k] = sum
    #need to chck this is right..


    for i in range(array_size):

    	cos_array[i][0]  = np.cos(np.radians(dihedral_array[i][5]))

    cossqr_array = np.square(cos_array)

    Jvalues_array = np.add((p1*cossqr_array),(p2*cos_array))

    Jvalues_array2 = np.add(Jvalues_array, F_array)



    for i in range(array_size):

    	HLA_array[i][0] = dihedral_array[i][0]

    	HLA_array[i][1] = dihedral_array[i][5]

    	HLA_array[i][2] = F_array[i]

    	HLA_array[i][3] = Jvalues_array2[i]

    	HLA_array[i][4] = dihedral_array[i][7]



    length = 0
    for i in range(array_size):
    	if dihedral_array[i][0] != 0:
    		length += 1

    HLA_final = np.zeros((length, 9), dtype=np.float64)
    name_list = []
    a = 0
    for i in range(array_size):
    	if dihedral_array[i][0] !=  0:
    		HLA_final[a][0] = dihedral_array[i][0]
    		HLA_final[a][1] = dihedral_array[i][1]
    		HLA_final[a][2] = dihedral_array[i][2]
    		HLA_final[a][3] = dihedral_array[i][3]
    		HLA_final[a][4] = dihedral_array[i][4]
    		HLA_final[a][5] = dihedral_array[i][5]
    		HLA_final[a][6] = Jvalues_array2[i]
    		HLA_final[a][7] = dihedral_array[i][7]
    		name_list.append(name)

		a +=1
  

    lengthwo7 = 0
    for i in range(array_size):
	if dihedral_array[i][0] != 0 and dihedral_array[i][0] != 7:
		lengthwo7 += 1
    a = 0 
    regression_array = np.zeros((lengthwo7, 10), dtype=np.float64)
    for i in range(array_size):
	if dihedral_array[i][0] != 0 and dihedral_array[i][0] != 7:
		regression_array[a][0] = dihedral_array[i][5]
		regression_array[a][1] = electronegativity_Aarray[i][0]
		regression_array[a][2] = electronegativity_Barray[i][0]
		regression_array[a][3] = electronegativity_Aarray[i][1]
		regression_array[a][4] = electronegativity_Barray[i][1]
		regression_array[a][5] = electronegativity_Aarray[i][2]
		regression_array[a][6] = electronegativity_Barray[i][2]
		regression_array[a][7] = electronegativity_Aarray[i][3]
		regression_array[a][8] = electronegativity_Barray[i][3]
		regression_array[a][9] = dihedral_array[i][7]
		a += 1
    outfile = "BAPOCM10.txt"

    #new_outfile = name + "_HLA.out"
   # outfile = "regression_input.txt"

    with open(outfile, "w") as f:

    	for i in range(lengthwo7):

  	#	string = "{0:<16.6f}, {1:<16.6f}, {2:<16.6f}, {3:<16.6f}, {4:<16.6f}, {5:<16.6f}, {6:<16.6f}, {7:<16.6f}, {8:<16s} ".format(HLA_final[i][0], HLA_final[i][1], HLA_final[i][2], HLA_final[i][3], HLA_final[i][4],  HLA_final[i][5],  HLA_final[i][6],  HLA_final[i][7], name_list[i])


		string = "{0:<16.6f}, {1:<16.6f}, {2:<16.6f}, {3:<16.6f}, {4:<16.6f}, {5:<16.6f}, {6:<16.6f}, {7:<16.6f}, {8:<16.6f}, {9:<16.6f}".format(regression_array[i][0], regression_array[i][1], regression_array[i][2], regression_array[i][3], regression_array[i][4],  regression_array[i][5],  regression_array[i][6], regression_array[i][7], regression_array[i][8], regression_array[i][9]) 
    		print(string, file = f)















#
