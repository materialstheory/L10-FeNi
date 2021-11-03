
#Generates a 32-atom supercell corresponding to a given order parameter P_z with P_x = P_y = 0 

import random
import numpy as np
#import numpy.random as random

# 4-atom based unit cell
"""

							Z
							|
							|
							|
							|
							|
							|
							|
							|
							|
							|
					Gamma		|
							|		delta
							|
							|
							|
							|
						   alpha|______________________________________________ X
						 	/
						       /
						      /
						     /
						    /
						   /			beta
						  /
					         /
						/
					       /
					      Y



"""
# To store coordinates of alpha, beta, gamma and delta sites in a 32-atom supercell


positions_alpha = {}
positions_beta = {}
positions_gamma = {}
positions_delta = {}




positions_Fe = []
positions_Ni  = []

pos_1 = [0.0,0.0,0.0]
pos_2 = [0.0,0.0,0.5]
pos_3 = [0.0,0.5,0.0]
pos_4 = [0.0,0.5,0.5]
pos_5 = [0.5,0.0,0.0]
pos_6 = [0.5,0.0,0.5]
pos_7 = [0.5,0.5,0.0]
pos_8 = [0.5,0.5,0.5]

pos_9 = [0.25,0.25,0.0]
pos_10 = [0.25,0.25,0.5]
pos_11 = [0.25,0.75,0.0]
pos_12 = [0.25,0.75,0.5]
pos_13 = [0.75,0.25,0.0]
pos_14 = [0.75,0.25,0.5]
pos_15 = [0.75,0.75,0.0]
pos_16 = [0.75,0.75,0.5]

pos_17 = [0.0,0.25,0.25]
pos_18 = [0.0,0.25,0.75]
pos_19 = [0.0,0.75,0.25]
pos_20 = [0.0,0.75,0.75]
pos_21 = [0.5,0.25,0.25]
pos_22 = [0.5,0.25,0.75]
pos_23 = [0.5,0.75,0.25]
pos_24 = [0.5,0.75,0.75]

pos_25 = [0.25,0.0,0.25]
pos_26 = [0.25,0.0,0.75]
pos_27 = [0.25,0.5,0.25]
pos_28 = [0.25,0.5,0.75]
pos_29 = [0.75,0.0,0.25]
pos_30 = [0.75,0.0,0.75]
pos_31 = [0.75,0.5,0.25]
pos_32 = [0.75,0.5,0.75]




positions_alpha = [pos_1, pos_2, pos_3, pos_4, pos_5, pos_6, pos_7, pos_8]

positions_beta = [pos_9, pos_10, pos_11, pos_12, pos_13, pos_14, pos_15, pos_16]

positions_gamma = [pos_17, pos_18, pos_19, pos_20, pos_21, pos_22, pos_23, pos_24]

positions_delta = [pos_25, pos_26, pos_27, pos_28, pos_29, pos_30, pos_31, pos_32]




# For z direction, a lattice sites --> alpha + beta

positions_alpha_beta = positions_alpha + positions_beta		# Lattice coordinates of all sites that belong to a 
positions_gamma_delta = positions_gamma + positions_delta	





# For y direction, a lattice sites --> alpha + delta


positions_alpha_delta = positions_alpha + positions_delta   	# Lattice coordinates of all sites that belong to a 






# For x direction, a lattice sites --> alpha + gamma

positions_alpha_gamma = positions_alpha + positions_gamma   # Lattice coordinates of all sites that belong to a 






all_positions = positions_alpha + positions_beta + positions_gamma + positions_delta	# Lattice coordinates of all 32-sites



# For ordered L10 structure along Z direction

########################

N_Fe_alpha = 8
N_Fe_beta = 8
N_Fe_gamma = 0
N_Fe_delta = 0 


########################



N_Fe_z_a = N_Fe_alpha + N_Fe_beta    # Total number of Fe atoms on site a (alpha + beta)

N_Fe_z_b = N_Fe_gamma + N_Fe_delta   # Total number of Fe atoms on site b (gamma + delta)




P_z = input('Enter the long range order parameter p (in multiples of 1.0/8.0) : ')  # Order parameter in Z direction
print "Long range order parameter in z direction (P_z) =", P_z


N_Fe_z_a = (P_z + 1.0)*8.0

N_Fe_z_b = 16 - N_Fe_z_a


number = 0 

while number < 10000:
	
	positions_Fe = []
	positions_Ni  = []

	positions_Fe = random.sample(positions_alpha_beta,  int(N_Fe_z_a)) + random.sample(positions_gamma_delta,  int(N_Fe_z_b))   #Distributing 16-Fe atoms on sites alpha, beta, gamma and delta


	#Checking Order parameter P_y in y direction...... 


	N_Fe_y_a = 0 # To count the number of iron atoms on lattice sites a 

	for i in positions_Fe:
	
		if i in positions_alpha_delta:
		
			N_Fe_y_a = N_Fe_y_a + 1

	P_y = N_Fe_y_a/8.0 - 1.0

	print "Long range order parameter in y direction (P_y) =", P_y


	#Checking Order parameter P_x in x direction.....


	N_Fe_x_a = 0     # To count numeber of iron atoms on lattice sites a 

	for i in positions_Fe:
	
		if i in positions_alpha_gamma:
		
			N_Fe_x_a = N_Fe_x_a + 1

	P_x = N_Fe_x_a/8.0 - 1.0

	print "Long range order parameter in x direction (P_x) =", P_x



	Fe_atoms = 16
	Ni_atoms = 16


	if P_y == 0 and P_x == 0:



		x_coord = []
		y_coord = []
		z_coord = []
		element = []



		for i in range(len(positions_Fe)):
		    	#print>>OutFile,'%.3f' % positions_Fe[i][0],'','%.3f' % positions_Fe[i][1],'','%.3f' % positions_Fe[i][2],'','Fe'

			x_coord.append(positions_Fe[i][0])
			y_coord.append(positions_Fe[i][1])	
			z_coord.append(positions_Fe[i][2])
			element.append('Fe')

		#print x_coord,y_coord,z_coord,element




		for i in all_positions:

		    if i not in positions_Fe:		# Finding sites not occupied by Fe atoms...those sites will be occupied by Ni atoms
			positions_Ni.append(i)
		
		
		for i in range(len(positions_Ni)):
		    	#print>>OutFile,'%.3f' % positions_Ni[i][0],'','%.3f' % positions_Ni[i][1],'','%.3f' % positions_Ni[i][2],'','Ni'


			x_coord.append(positions_Ni[i][0])
			y_coord.append(positions_Ni[i][1])	
			z_coord.append(positions_Ni[i][2])
			element.append('Ni')

		
		#print x_coord,y_coord,z_coord,element
		


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
		

		if (-0.2 <= sum_alpha_first_coord_shell <= 0.2 and -0.2 <= sum_alpha_second_coord_shell <= 0.2):

			# Writing to the POSCAR file 

			OutName=raw_input('Insert file name : \n')
			OutFile=open(OutName,'w')
			print>>OutFile,'NiFe L1_0 structure'
			print>>OutFile,'7.12'
			print>>OutFile,' 1.0 0.0 0.0'
			print>>OutFile,' 0.0 1.0 0.0'
			print>>OutFile,' 0.0 0.0 1.0'
			print>>OutFile,' Fe Ni'
			print>>OutFile,'',Fe_atoms,Ni_atoms
			print>>OutFile,'Direct'

			for i in range(len(positions_Fe)):
			    	print>>OutFile,'%.3f' % positions_Fe[i][0],'','%.3f' % positions_Fe[i][1],'','%.3f' % positions_Fe[i][2],'','Fe'

			for i in range(len(positions_Ni)):
			    	print>>OutFile,'%.3f' % positions_Ni[i][0],'','%.3f' % positions_Ni[i][1],'','%.3f' % positions_Ni[i][2],'','Ni'
			OutFile.close()
	else:
		print "No Luck :("
	number = number + 1


