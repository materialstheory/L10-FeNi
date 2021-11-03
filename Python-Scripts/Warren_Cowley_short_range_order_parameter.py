# Following script calculates the Warren-Cowley short range order parameters (SRO) in FCC structures

import numpy as np
import numpy.random as random


# To read the parameter file

def read_data(filename):
	with open(filename, 'r') as f:
   		data = f.readlines()[8:]						# to omit first 8 lines
	
	x_coord = []
	y_coord = []
	z_coord = []
	element = [] 

	for i in range(len(data)):
		

		x_coord.append(float(data[i].split()[0]))
		y_coord.append(float(data[i].split()[1]))	
		z_coord.append(float(data[i].split()[2]))
		element.append(str(data[i].split()[3]))

		



	return x_coord,y_coord,z_coord,element


x_coord,y_coord,z_coord,element = read_data('POSCAR')


print x_coord,y_coord,z_coord,element




L = 1.0
Total_atoms = 32

FNN = 12		# No. of first Nearest neighbours (FNN) 
SNN = 6		# No. of second Nearest neighbours (SNN)
TNN = 24 		# No. of third Nearest neighbours (TNN)

c_Fe = 0.5	# fraction of Fe 
c_Ni = 0.5	# fraction of Ni

Coordinates = {}	# To store the coordinates of each atom
Coordination_FNN = {} 	# To store coordinates of first nearest neighbours of each atom 
Coordination_SNN = {}	# To store coordinates of second nearest neighbours of each atom
Coordination_TNN = {}	# To store coordinates of third nearest neighbours of each atom


for i in range(len(x_coord)):

	Coordinates[i] = [x_coord[i],y_coord[i],z_coord[i]]




for i in range(len(x_coord)):
	
			
	# First Nearest neighbours in FCC = 12
	
	Coordination_FNN[i] =	[[np.mod(x_coord[i] + L/4.0,L),np.mod(y_coord[i] + L/4.0,L),z_coord[i]],[np.mod(x_coord[i] - L/4.0,L),np.mod(y_coord[i] - L/4.0,L),z_coord[i]],[np.mod(x_coord[i] - L/4.0,L),np.mod(y_coord[i] + L/4.0,L),z_coord[i]],[np.mod(x_coord[i] + L/4.0,L),np.mod(y_coord[i] - L/4.0,L),z_coord[i]],[np.mod(x_coord[i] + L/4.0,L),y_coord[i],np.mod(z_coord[i] + L/4.0,L)],[np.mod(x_coord[i] - L/4.0,L),y_coord[i],np.mod(z_coord[i] - L/4.0,L)],[np.mod(x_coord[i] + L/4.0,L),y_coord[i],np.mod(z_coord[i] - L/4.0,L)],[np.mod(x_coord[i] - L/4.0,L),y_coord[i],np.mod(z_coord[i] + L/4.0,L)],[x_coord[i],np.mod(y_coord[i] + L/4.0,L),np.mod(z_coord[i] + L/4.0,L)],[x_coord[i],np.mod(y_coord[i] - L/4.0,L),np.mod(z_coord[i] - L/4.0,L)],[x_coord[i],np.mod(y_coord[i] + L/4.0,L),np.mod(z_coord[i] - L/4.0,L)],[x_coord[i],np.mod(y_coord[i] - L/4.0,L),np.mod(z_coord[i] + L/4.0,L)]]
	




	# Second nearest neighbours in FCC = 6

for i in range(len(x_coord)):

	Coordination_SNN[i] = [[np.mod(x_coord[i] + L/2.0,L),y_coord[i],z_coord[i]],[np.mod(x_coord[i] - L/2.0,L),y_coord[i],z_coord[i]],[x_coord[i],np.mod(y_coord[i] + L/2.0,L),z_coord[i]],[x_coord[i],np.mod(y_coord[i] - L/2.0,L),z_coord[i]],[x_coord[i],y_coord[i],np.mod(z_coord[i] + L/2.0,L)],[x_coord[i],y_coord[i],np.mod(z_coord[i] - L/2.0,L)]]



	# Third nearest neighbours in FCC = 24

for i in range(len(x_coord)):

	Coordination_TNN[i] = [[np.mod(x_coord[i] + L/4.0,L),np.mod(y_coord[i] + L/4.0,L),np.mod(z_coord[i] + L/2.0,L)],[np.mod(x_coord[i] - L/4.0,L),np.mod(y_coord[i] - L/4.0,L),np.mod(z_coord[i] + L/2.0,L)],[np.mod(x_coord[i] - L/4.0,L),np.mod(y_coord[i] + L/4.0,L),np.mod(z_coord[i] + L/2.0,L)],[np.mod(x_coord[i] + L/4.0,L),np.mod(y_coord[i] - L/4.0,L),np.mod(z_coord[i] + L/2.0,L)],[np.mod(x_coord[i] + L/4.0,L),np.mod(y_coord[i] + L/4.0,L),np.mod(z_coord[i] - L/2.0,L)],[np.mod(x_coord[i] - L/4.0,L),np.mod(y_coord[i] - L/4.0,L),np.mod(z_coord[i] - L/2.0,L)],[np.mod(x_coord[i] - L/4.0,L),np.mod(y_coord[i] + L/4.0,L),np.mod(z_coord[i] - L/2.0,L)],[np.mod(x_coord[i] + L/4.0,L),np.mod(y_coord[i] - L/4.0,L),np.mod(z_coord[i] - L/2.0,L)],[np.mod(x_coord[i] + L/2.0,L),np.mod(y_coord[i] + L/4.0,L),np.mod(z_coord[i] + L/4.0,L)],[np.mod(x_coord[i] + L/2.0,L),np.mod(y_coord[i] - L/4.0,L),np.mod(z_coord[i] + L/4.0,L)],[np.mod(x_coord[i] - L/2.0,L),np.mod(y_coord[i] + L/4.0,L),np.mod(z_coord[i] + L/4.0,L)],[np.mod(x_coord[i] - L/2.0,L),np.mod(y_coord[i] - L/4.0,L),np.mod(z_coord[i] + L/4.0,L)],[np.mod(x_coord[i] + L/4.0,L),np.mod(y_coord[i] + L/2.0,L),np.mod(z_coord[i] + L/4.0,L)],[np.mod(x_coord[i] + L/4.0,L),np.mod(y_coord[i] - L/2.0,L),np.mod(z_coord[i] + L/4.0,L)],[np.mod(x_coord[i] - L/4.0,L),np.mod(y_coord[i] + L/2.0,L),np.mod(z_coord[i] + L/4.0,L)],[np.mod(x_coord[i] - L/4.0,L),np.mod(y_coord[i] - L/2.0,L),np.mod(z_coord[i] + L/4.0,L)],[np.mod(x_coord[i] + L/2.0,L),np.mod(y_coord[i] + L/4.0,L),np.mod(z_coord[i] - L/4.0,L)],[np.mod(x_coord[i] + L/2.0,L),np.mod(y_coord[i] - L/4.0,L),np.mod(z_coord[i] - L/4.0,L)],[np.mod(x_coord[i] - L/2.0,L),np.mod(y_coord[i] + L/4.0,L),np.mod(z_coord[i] - L/4.0,L)],[np.mod(x_coord[i] - L/2.0,L),np.mod(y_coord[i] - L/4.0,L),np.mod(z_coord[i] - L/4.0,L)],[np.mod(x_coord[i] + L/4.0,L),np.mod(y_coord[i] + L/2.0,L),np.mod(z_coord[i] - L/4.0,L)],[np.mod(x_coord[i] + L/4.0,L),np.mod(y_coord[i] - L/2.0,L),np.mod(z_coord[i] - L/4.0,L)],[np.mod(x_coord[i] - L/4.0,L),np.mod(y_coord[i] + L/2.0,L),np.mod(z_coord[i] - L/4.0,L)],[np.mod(x_coord[i] - L/4.0,L),np.mod(y_coord[i] - L/2.0,L),np.mod(z_coord[i] - L/4.0,L)]]


sum_alpha_first_coord_shell = 0
sum_alpha_second_coord_shell = 0
sum_alpha_third_coord_shell = 0

for i in range(len(x_coord)):

	sum_first_nearest_Fe_neighbours = 0
	sum_first_nearest_Ni_neighbours = 0
	sum_second_nearest_Fe_neighbours = 0
	sum_second_nearest_Ni_neighbours = 0
	sum_third_nearest_Fe_neighbours = 0
	sum_third_nearest_Ni_neighbours = 0


	for j in range(len(x_coord)):
		


		for k in range(len(Coordination_FNN[i])):
			if Coordination_FNN[i][k] == Coordinates[j]:
				 
				if element[j] == 'Fe':	
					sum_first_nearest_Fe_neighbours = sum_first_nearest_Fe_neighbours + 1
				else:
					sum_first_nearest_Ni_neighbours = sum_first_nearest_Ni_neighbours + 1




		for l in range(len(Coordination_SNN[i])):
			if Coordination_SNN[i][l] == Coordinates[j]:
		
				if element[j] == 'Fe':	
					sum_second_nearest_Fe_neighbours = sum_second_nearest_Fe_neighbours + 1
				else:
					sum_second_nearest_Ni_neighbours = sum_second_nearest_Ni_neighbours + 1



		for m in range(len(Coordination_TNN[i])):
			if Coordination_TNN[i][m] == Coordinates[j]:
		
				if element[j] == 'Fe':	
					sum_third_nearest_Fe_neighbours = sum_third_nearest_Fe_neighbours + 1
				else:
					sum_third_nearest_Ni_neighbours = sum_third_nearest_Ni_neighbours + 1



# Warren-Cowley short range order parameter

	if element[i] == 'Fe':
		alpha_first_coord_shell = 1 - sum_first_nearest_Ni_neighbours/(float(FNN)*c_Ni)
		alpha_second_coord_shell = 1- sum_second_nearest_Ni_neighbours/(float(SNN)*c_Ni)
		alpha_third_coord_shell = 1 -  sum_second_nearest_Ni_neighbours/(float(TNN)*c_Ni)
		print element[i],alpha_first_coord_shell,alpha_second_coord_shell,alpha_third_coord_shell


	else:
		alpha_first_coord_shell = 1 - sum_first_nearest_Fe_neighbours/(float(FNN)*c_Fe)
		alpha_second_coord_shell = 1- sum_second_nearest_Fe_neighbours/(float(SNN)*c_Fe)
		alpha_third_coord_shell = 1 - sum_second_nearest_Fe_neighbours/(float(TNN)*c_Fe )
		print element[i],alpha_first_coord_shell,alpha_second_coord_shell,alpha_third_coord_shell


	sum_alpha_first_coord_shell = sum_alpha_first_coord_shell + alpha_first_coord_shell

	sum_alpha_second_coord_shell = sum_alpha_second_coord_shell + alpha_second_coord_shell

	sum_alpha_third_coord_shell = sum_alpha_third_coord_shell + alpha_third_coord_shell




print "Warren-Cowley short range order parameter(First Nearest Neighbour) = ", sum_alpha_first_coord_shell/float(Total_atoms)

print "Short range order parameter(Second Nearest Neighbour) = ", sum_alpha_second_coord_shell/float(Total_atoms)

print "Short range order parameter(third Nearest Neighbour) = ", sum_alpha_third_coord_shell/float(Total_atoms)

