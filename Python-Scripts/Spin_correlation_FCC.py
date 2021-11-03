"Following python script reads INCAR file and outputs the spin correlation function for first, second, and third nearest neighbour in a Face-Centered-Cubic(FCC) structure"


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



#print element


def read_data(filename):
	with open(filename, 'r') as f:


		spin = []
 		for line in f:
	
	
			if line.startswith("MAGMOM"):
				
				fields = line.split()[2:]
				for i in fields:
				
					spin.append(i)				




	return spin

spin = read_data('INCAR')  # To store the magnetic moment along with the spin directions 

spin_direction = []	# Assumption --> Collinear spins i.e. +1 or -1


# To store only the spin direction(+1 or -1) for each atom

for i in spin:

	i = int(i)

	if i == 0:
		spin_direction.append(i)
	else:
		spin_direction.append(i/np.absolute(i))
	


L = 1.0
Total_atoms = 32

Coordinates = {}	# To store the coordinates of each atom
Coordination_FNN = {} 	# To store coordinates of first nearest neighbours of each atom 
Coordination_SNN = {}	# To store coordinates of second nearest neighbours of each atom
Coordination_TNN = {}	# To store coordinates of third nearest neighbours of each atom


for i in range(len(x_coord)):

	Coordinates[i] = [x_coord[i],y_coord[i],z_coord[i]]




for i in range(len(x_coord)):
	
			
	# Nearest neighbours in FCC = 12
	
	Coordination_FNN[i] =	[[np.mod(x_coord[i] + L/4.0,L),np.mod(y_coord[i] + L/4.0,L),z_coord[i]],[np.mod(x_coord[i] - L/4.0,L),np.mod(y_coord[i] - L/4.0,L),z_coord[i]],[np.mod(x_coord[i] - L/4.0,L),np.mod(y_coord[i] + L/4.0,L),z_coord[i]],[np.mod(x_coord[i] + L/4.0,L),np.mod(y_coord[i] - L/4.0,L),z_coord[i]],[np.mod(x_coord[i] + L/4.0,L),y_coord[i],np.mod(z_coord[i] + L/4.0,L)],[np.mod(x_coord[i] - L/4.0,L),y_coord[i],np.mod(z_coord[i] - L/4.0,L)],[np.mod(x_coord[i] + L/4.0,L),y_coord[i],np.mod(z_coord[i] - L/4.0,L)],[np.mod(x_coord[i] - L/4.0,L),y_coord[i],np.mod(z_coord[i] + L/4.0,L)],[x_coord[i],np.mod(y_coord[i] + L/4.0,L),np.mod(z_coord[i] + L/4.0,L)],[x_coord[i],np.mod(y_coord[i] - L/4.0,L),np.mod(z_coord[i] - L/4.0,L)],[x_coord[i],np.mod(y_coord[i] + L/4.0,L),np.mod(z_coord[i] - L/4.0,L)],[x_coord[i],np.mod(y_coord[i] - L/4.0,L),np.mod(z_coord[i] + L/4.0,L)]]
	




	# Next nearest neighbours in FCC = 6

for i in range(len(x_coord)):

	Coordination_SNN[i] = [[np.mod(x_coord[i] + L/2.0,L),y_coord[i],z_coord[i]],[np.mod(x_coord[i] - L/2.0,L),y_coord[i],z_coord[i]],[x_coord[i],np.mod(y_coord[i] + L/2.0,L),z_coord[i]],[x_coord[i],np.mod(y_coord[i] - L/2.0,L),z_coord[i]],[x_coord[i],y_coord[i],np.mod(z_coord[i] + L/2.0,L)],[x_coord[i],y_coord[i],np.mod(z_coord[i] - L/2.0,L)]]



	# Third nearest neighbours in FCC = 24

for i in range(len(x_coord)):

	Coordination_TNN[i] = [[np.mod(x_coord[i] + L/4.0,L),np.mod(y_coord[i] + L/4.0,L),np.mod(z_coord[i] + L/2.0,L)],[np.mod(x_coord[i] - L/4.0,L),np.mod(y_coord[i] - L/4.0,L),np.mod(z_coord[i] + L/2.0,L)],[np.mod(x_coord[i] - L/4.0,L),np.mod(y_coord[i] + L/4.0,L),np.mod(z_coord[i] + L/2.0,L)],[np.mod(x_coord[i] + L/4.0,L),np.mod(y_coord[i] - L/4.0,L),np.mod(z_coord[i] + L/2.0,L)],[np.mod(x_coord[i] + L/4.0,L),np.mod(y_coord[i] + L/4.0,L),np.mod(z_coord[i] - L/2.0,L)],[np.mod(x_coord[i] - L/4.0,L),np.mod(y_coord[i] - L/4.0,L),np.mod(z_coord[i] - L/2.0,L)],[np.mod(x_coord[i] - L/4.0,L),np.mod(y_coord[i] + L/4.0,L),np.mod(z_coord[i] - L/2.0,L)],[np.mod(x_coord[i] + L/4.0,L),np.mod(y_coord[i] - L/4.0,L),np.mod(z_coord[i] - L/2.0,L)],[np.mod(x_coord[i] + L/2.0,L),np.mod(y_coord[i] + L/4.0,L),np.mod(z_coord[i] + L/4.0,L)],[np.mod(x_coord[i] + L/2.0,L),np.mod(y_coord[i] - L/4.0,L),np.mod(z_coord[i] + L/4.0,L)],[np.mod(x_coord[i] - L/2.0,L),np.mod(y_coord[i] + L/4.0,L),np.mod(z_coord[i] + L/4.0,L)],[np.mod(x_coord[i] - L/2.0,L),np.mod(y_coord[i] - L/4.0,L),np.mod(z_coord[i] + L/4.0,L)],[np.mod(x_coord[i] + L/4.0,L),np.mod(y_coord[i] + L/2.0,L),np.mod(z_coord[i] + L/4.0,L)],[np.mod(x_coord[i] + L/4.0,L),np.mod(y_coord[i] - L/2.0,L),np.mod(z_coord[i] + L/4.0,L)],[np.mod(x_coord[i] - L/4.0,L),np.mod(y_coord[i] + L/2.0,L),np.mod(z_coord[i] + L/4.0,L)],[np.mod(x_coord[i] - L/4.0,L),np.mod(y_coord[i] - L/2.0,L),np.mod(z_coord[i] + L/4.0,L)],[np.mod(x_coord[i] + L/2.0,L),np.mod(y_coord[i] + L/4.0,L),np.mod(z_coord[i] - L/4.0,L)],[np.mod(x_coord[i] + L/2.0,L),np.mod(y_coord[i] - L/4.0,L),np.mod(z_coord[i] - L/4.0,L)],[np.mod(x_coord[i] - L/2.0,L),np.mod(y_coord[i] + L/4.0,L),np.mod(z_coord[i] - L/4.0,L)],[np.mod(x_coord[i] - L/2.0,L),np.mod(y_coord[i] - L/4.0,L),np.mod(z_coord[i] - L/4.0,L)],[np.mod(x_coord[i] + L/4.0,L),np.mod(y_coord[i] + L/2.0,L),np.mod(z_coord[i] - L/4.0,L)],[np.mod(x_coord[i] + L/4.0,L),np.mod(y_coord[i] - L/2.0,L),np.mod(z_coord[i] - L/4.0,L)],[np.mod(x_coord[i] - L/4.0,L),np.mod(y_coord[i] + L/2.0,L),np.mod(z_coord[i] - L/4.0,L)],[np.mod(x_coord[i] - L/4.0,L),np.mod(y_coord[i] - L/2.0,L),np.mod(z_coord[i] - L/4.0,L)]]









# To calculate the correlation function for first nearest neighbours 

print "\n"
count_i = 0
spin_correlation_fnn = 0.0	# Spin Correlation function (nearest neighbours) 
for i in range(len(x_coord)):
	spin_sum_fnn = 0
	sum_first_nearest_Fe_neighbours = 0
	for j in range(len(x_coord)):
		
		for k in range(len(Coordination_FNN[i])):
			if Coordination_FNN[i][k] == Coordinates[j]:
				 
				spin_sum_fnn = spin_sum_fnn + spin_direction[j]	# Summing all the spins of nearest neighbours --> Sigma (j)
	
				if spin_direction[j] != 0:	# Ni atoms with 0 uB should not be counted in the neighbours !
					sum_first_nearest_Fe_neighbours = sum_first_nearest_Fe_neighbours + np.absolute(spin_direction[j])

	print sum_first_nearest_Fe_neighbours,spin_sum_fnn
	if sum_first_nearest_Fe_neighbours == 0 :

		print "Correlation function(FNN) (atom %i ) = "%(i+1),(spin_sum_fnn*spin_direction[i])	
		spin_correlation_fnn = spin_correlation_fnn + (spin_sum_fnn * spin_direction[i])	 #  Correlation(FNN) =  Sigma(i) * Sigma(j)
	else:		
		print "Correlation function(NN) (atom %i ) = "%(i+1),(spin_sum_fnn*spin_direction[i])/float(sum_first_nearest_Fe_neighbours)
		spin_correlation_fnn = spin_correlation_fnn + (spin_sum_fnn * spin_direction[i])/float(sum_first_nearest_Fe_neighbours)	 #  Correlation(DNN) =  Sigma(i) * Sigma(j)

	print spin_correlation_fnn



	if spin_sum_fnn*spin_direction[i] != 0:	# To filter atoms with zero magnetic moment (they will not contribute to the averaging in spin_correlation function)
		count_i = count_i + 1


print count_i,spin_correlation_fnn


if count_i == 0:
	print "Correlation function (FNN) = ",spin_correlation_fnn
else:
 
	print "Correlation function (FNN) = ",spin_correlation_fnn/float(count_i)






spin_correlation_snn = 0.0	# Spin Correlation function (next nearest neighbours)

print "\n"

count_ii = 0	# To count number of atoms with magnetic moments 


for ii in range(len(x_coord)):
	spin_sum_snn = 0
	sum_second_nearest_Fe_neighbours = 0.0

	for j in range(len(x_coord)):
		for k in range(len(Coordination_SNN[ii])):
			if Coordination_SNN[ii][k] == Coordinates[j]:
		
				spin_sum_snn = spin_sum_snn + spin_direction[j]


				if spin_direction[j] != 0:	# Ni atoms with 0 uB should not be counted in the neighbours !
					sum_second_nearest_Fe_neighbours = sum_second_nearest_Fe_neighbours + np.absolute(spin_direction[j])
	print "Number of Fe atoms = ", sum_second_nearest_Fe_neighbours
	print spin_sum_snn,spin_direction[ii]
	if sum_second_nearest_Fe_neighbours == 0 :

		print "Correlation function(SNN) (atom %i ) = "%(ii+1),(spin_sum_snn*spin_direction[ii])	
		spin_correlation_snn = spin_correlation_snn + (spin_sum_snn * spin_direction[ii])	 #  Correlation(NN) =  Sigma(i) * Sigma(j)
	else:		
		print "Correlation function(SNN) (atom %i ) = "%(ii+1),(spin_sum_snn*spin_direction[ii])/sum_second_nearest_Fe_neighbours
		spin_correlation_snn = spin_correlation_snn + (spin_sum_snn * spin_direction[ii])/sum_second_nearest_Fe_neighbours	 #  Correlation(NN) =  Sigma(i) * Sigma(j)


	if spin_sum_snn*spin_direction[ii] != 0:	# To filter atoms with zero magnetic moment (they will not contribute to the averaging in corrletion function)
		count_ii = count_ii + 1
	
print "Correlation function (SNN) = ",spin_correlation_snn/float(count_ii)
	









# To calculate the correlation function for third nearest neighbours 

print "\n"
count_iii = 0
spin_correlation_tnn = 0.0	# Spin Correlation function (nearest neighbours) 
for i in range(len(x_coord)):
	spin_sum_tnn = 0
	sum_third_nearest_Fe_neighbours = 0
	for j in range(len(x_coord)):
		
		for k in range(len(Coordination_TNN[i])):
			if Coordination_TNN[i][k] == Coordinates[j]:
				 
				spin_sum_tnn = spin_sum_tnn + spin_direction[j]	# Summing all the spins of nearest neighbours --> Sigma (j)
	
				if spin_direction[j] != 0:	# Ni atoms with 0 uB should not be counted in the neighbours !
					sum_third_nearest_Fe_neighbours = sum_third_nearest_Fe_neighbours + np.absolute(spin_direction[j])

	print sum_third_nearest_Fe_neighbours,spin_sum_tnn
	if sum_third_nearest_Fe_neighbours == 0 :

		print "Correlation function(FNN) (atom %i ) = "%(i+1),(spin_sum_tnn*spin_direction[i])	
		spin_correlation_tnn = spin_correlation_tnn + (spin_sum_tnn * spin_direction[i])	 #  Correlation(TNN) =  Sigma(i) * Sigma(j)
	else:		
		print "Correlation function(NN) (atom %i ) = "%(i+1),(spin_sum_fnn*spin_direction[i])/float(sum_third_nearest_Fe_neighbours)
		spin_correlation_tnn = spin_correlation_tnn + (spin_sum_tnn * spin_direction[i])/float(sum_third_nearest_Fe_neighbours)	 #  Correlation(TNN) =  Sigma(i) * Sigma(j)

	print spin_correlation_tnn


	if spin_sum_tnn*spin_direction[i] != 0:	# To filter atoms with zero magnetic moment (they will not contribute to the averaging in spin_correlation function)
		count_iii = count_iii + 1


print count_iii,spin_correlation_fnn


if count_iii == 0:
	print "Correlation function (TNN) = ",spin_correlation_tnn
else:
 
	print "Correlation function (TNN) = ",spin_correlation_tnn/float(count_iii)




