# ej14716
#
#
#
#


from __future__ import print_function

import math

import numpy as np

import sys

np.set_printoptions(threshold = 3000) 
 
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

 


file1 = "HC064_BAPOCM10_OPT_NMR_20_dangle.txt"

file2 = "HC064_BAPOCM10_OPT_NMR_20_atno.txt"

file3 = "HC064_BAPOCM10_OPT_NMR_20_conn.txt"



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

sub1 = []

sub2 = []

sub3 = []

sub4 = []




					

electronegativity_array = [

50,

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



connectivity_array = np.zeros((442, 49), dtype = np.float64)	



x = 0

y = 0

	

with open(file3, "r") as f:

	for line in f:

		items = line.split(",")
	

		for item in items:

			if not is_number(item):
				continue
			connectivity_array[x][y] = item

			y += 1

		x += 1

		y = 0



		

		

dihedral_array = np.zeros((array_size, 7), dtype = np.float64)

		

for i in range(array_size):

	dihedral_array[i][0] = 0

	dihedral_array[i][1] = atom_1[i]

	dihedral_array[i][2] = atom_2[i]

	dihedral_array[i][3] = atom_3[i]

	dihedral_array[i][4] = atom_4[i]

	dihedral_array[i][5] = dihedral_angle1[i]
 	dihedral_array[i][6] = i

	
	if atomic_number[atom_1[i]] == 1:
		
		if atomic_number[atom_2[i]] == 6:
			
			if atomic_number[atom_3[i]] == 6:
				
				if atomic_number[atom_4[i]] == 1:
					
					dihedral_array[i][0] = 1
	
a = 0
ivalues_failed = []
bond_angle1 = np.zeros((1, 2), dtype = np.float64)
bond_angle2 = np.zeros((1, 2), dtype = np.float64)

length = 0
for i in range(array_size):
	if dihedral_array[i][0] == 1:
		length = length +1
sub_array = np.zeros((length, 5), dtype = np.float64)
dihedral_coupling_array = np.zeros((length, 2), dtype = np.float64)
c = 0
for i in range(array_size):
	if dihedral_array[i][0] == 1:
		dihedral_coupling_array[c][0] = i
		dihedral_coupling_array[c][1] = dihedral_array[i][5]
		c += 1
Esub1B = []
Esub2B = []
Esub3B = []
Esub4B = []
electronegativity_Barray = np.zeros((length, 4), dtype = np.float64)

electronegativity_Aarray = np.zeros((length, 4), dtype = np.float64)	

A_array = np.zeros((length, 4), dtype = np.float64)

B_array = np.zeros((length, 4), dtype = np.float64)		

C_array = np.zeros((length, 4), dtype = np.float64)

D_array = np.zeros((length, 4), dtype = np.float64)

E_array = np.zeros((length, 4), dtype = np.float64)

F_array = np.zeros((length, 1), dtype = np.float64)

cos_array = np.zeros((length, 1), dtype = np.float64)

cossqr_array = np.zeros((length, 1), dtype = np.float64)

Jvalues_array = np.zeros((length, 1), dtype = np.float64)
Jvalues_array2 = np.zeros((length, 1), dtype = np.float64)
HLA_array = np.zeros((length, 4), dtype = np.float64)		

#get atom labels for substituents and their relative orientations

for i in range(array_size):		
	
#only look at pathways  that are HCCH			
	if dihedral_array[i][0] == 1:

		bonds_1 = np.transpose(np.nonzero(connectivity_array[dihedral_array[i][1]]))
		if not len(bonds_1) == 1:
			print("Hconnectivity_error")
			ivalues_failed.append(i)
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
				print("C1_notsp3", i )
				ivalues_failed.append(i)
				continue
			#search each of the connections

			for j in range (len(bonds_2)):

				if bonds_2[j] == dihedral_array[i][1]:
					
					continue

					#if the connection is to the third atom (C in coupling) search for what the C is connected for 

				elif bonds_2[j] == dihedral_array[i][3]:
										
					bonds_3 = np.transpose(np.nonzero(connectivity_array[dihedral_array[i][3]]))
					
					#loop through each of the bonds that C is connected to
 
					if not len(bonds_3) == 4:
						print("C2_notsp3", i)
						ivalues_failed.append(i)
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

											sub3.append(bonds_3[x])

										else:

											sub4.append(bonds_3[x])
		
									if (bond_angle1[b][1]< 0):

										if (abs(bond_angle1[b][1]) > 180) and (abs(bond_angle1[b][1]) < 360):
											sub3.append(bonds_3[x])	
										else:
	
											sub4.append(bonds_3[x])
								 
					
	
								
				else:
				

				#this is now for other bonds conncted to first carbon again
					
					for Y in range (array_size):
						c = 0
					#find dihedral angle for the substituents to the fourth atom (H)

						if (dihedral_array[Y][1] == bonds_2[j]) and (dihedral_array[Y][4] == dihedral_array[i][4]):

							bond_angle2[c][0] = dihedral_array[Y][5]

							bond_angle2[c][1] = dihedral_array[i][5] - bond_angle2[c][0]
							

							if (bond_angle2[c][1] > 0):

								if abs(bond_angle2[c][1]) < 180:

									sub1.append(bonds_2[j])
									
								else:

									sub2.append(bonds_2[j])

							if (bond_angle2[c][1] < 0):

								if (abs(bond_angle2[c][1]) > 180) and (abs(bond_angle2[c][1]) < 360):

									sub1.append(bonds_2[j])

								else:

									sub2.append(bonds_2[j])

					
						
#put the substituent labels into an array		

	
 	
		sub_array[a][0] = sub1[a]
		sub_array[a][1] = sub2[a]
		sub_array[a][2] = sub3[a]
		sub_array[a][3] = sub4[a]
		sub_array[a][4] = i
		a += 1


#now need to get electronegativity values for these substituents and the beta substituents

for a in range(length):

	bonds_4 = np.transpose(np.nonzero(connectivity_array[sub_array[a][0]]))
	
	Esub1B_1 = 0

	for z in range(len(bonds_4)):

		if bonds_4[z] == dihedral_array[sub_array[a][4]][2]:

			continue

		else:

			#assign and sum their electronegativities

			Esub1B_1 += electronegativity_array[atomic_number[bonds_4[z]]]

	bonds_5 = np.transpose(np.nonzero(connectivity_array[sub_array[a][1]]))	
	Esub2B_1 = 0

	for z in range(len(bonds_5)):

		if bonds_5[z] == dihedral_array[sub_array[a][4]][2]:

			continue
	
		else:

			#assign and sum their electronegativities

			Esub2B_1 += electronegativity_array[atomic_number[bonds_5[z]]]

	bonds_6 = np.transpose(np.nonzero(connectivity_array[sub_array[a][2]]))	

	Esub3B_1 = 0

	for z in range(len(bonds_6)):

		if bonds_6[z] == dihedral_array[sub_array[a][4]][3]:

			continue

		else:

			#assign and sum their electronegativities

			Esub3B_1 += electronegativity_array[atomic_number[bonds_6[z]]]					

	bonds_7 = np.transpose(np.nonzero(connectivity_array[sub_array[a][3]]))	

	Esub4B_1 = 0

	for z in range(len(bonds_7)):

		if bonds_7[z] == dihedral_array[sub_array[a][4]][3]:

			continue

		else:

			#assign and sum their electronegativities

			Esub4B_1 += electronegativity_array[atomic_number[bonds_7[z]]]

			

	Esub1B.append(Esub1B_1)

	Esub2B.append(Esub2B_1)

	Esub3B.append(Esub3B_1)

	Esub4B.append(Esub4B_1)
	
	

	#array containing the sum of the electronegativities of the Beta substituents of each alpha substituent for every value of i 	



	electronegativity_Barray[a][0] = Esub1B[a]

	electronegativity_Barray[a][1] = Esub2B[a]
	
	electronegativity_Barray[a][2] = Esub3B[a]

	electronegativity_Barray[a][3] = Esub4B[a]

	
#electronegativity of the alpha substituents array 	

for i in range(length):

	electronegativity_Aarray[i][0] = electronegativity_array[atomic_number[sub1[i]]]

	electronegativity_Aarray[i][1] = electronegativity_array[atomic_number[sub2[i]]]

	electronegativity_Aarray[i][2] = electronegativity_array[atomic_number[sub3[i]]]

	electronegativity_Aarray[i][3] = electronegativity_array[atomic_number[sub4[i]]]
				
# x the sum of the electronegativity vales of the beta substituents by p7

A_array	= p7 * (electronegativity_Barray)



B_array = electronegativity_Aarray - A_array





#do stuff with the positive and negtive substituents


for i in range(length):

	C_array[i][0] = np.radians(dihedral_array[sub_array[i][4]][5]) + (p6 * (np.absolute(B_array[i][0])))

	C_array[i][1] = np.negative(np.radians(dihedral_array[sub_array[i][4]][5])) + (p6 * (np.absolute(B_array[i][1])))

	C_array[i][2] = np.radians(dihedral_array[sub_array[i][4]][5]) + (p6 * (np.absolute(B_array[i][2])))

	C_array[i][3] = np.negative(np.radians(dihedral_array[sub_array[i][4]][5])) + (p6 * (np.absolute(B_array[i][3])))


# make aray that needs to be * by the electronegativity (B array)	

D_array = p4 + (p5 * np.square(np.cos(C_array)))


E_array = D_array * B_array

	
for k in range(length):
	sum = 0
	for l in range(4):
		sum = sum + E_array[k][l]
	F_array[k] = sum
#need to chck this is right..


for i in range(length):

	cos_array[i][0]  = np.cos(np.radians(dihedral_coupling_array[i][1]))

cossqr_array = np.square(cos_array)

Jvalues_array = np.add((p1*cossqr_array),(p2*cos_array))

Jvalues_array2 = np.add(Jvalues_array, F_array)



for i in range(length):

	HLA_array[i][0] = dihedral_coupling_array[i][0]

	HLA_array[i][1] = dihedral_coupling_array[i][1]

	HLA_array[i][2] = F_array[i]

	HLA_array[i][3] = Jvalues_array2[i]


print(HLA_array) 
outfile = "HLA.txt"



with open(outfile, "w") as f:

	for i in range(array_size):

		string = "{0:<16.6f}".format(HLA_array)

		print(string, file = f)
