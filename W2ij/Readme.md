# W2Jij

# This repository contains the source code to calculate the magnetic exchange interactions in FeNi using the Magnetic Force Theorem. 

## Introduction

Code to calculate magnetic exchange interactions with WANNIER90 output files wannier90.up_hr.dat and wannier90.dn_hr.dat. Information about the structure of the material is read from wannier90.wout. The method is described in https://journals.aps.org/prb/abstract/10.1103/PhysRevB.101.064401 and references therein.

### Example input

##### index of atom i, start from zero!
```bash
ni=0 
```

##### With Rmax>0 specify only once per index, otherwise list all
```bash
nj=[0,1,2,3]      
```

#### Maximum interatomic distance in Angstrom 
```bash
Rmax=26.0            
```

Note: Either specify Rmax and nj uniquely corresponding to those orbitals to which it should be calculated, then all Jij interactions are calculated for which the interatomic distance Rij < Rmax. This does not consider symmetries, so all neighbours are calculated even if equivalent.
Another option is shown below. Specify each atom j explicitly, with atomic/orbital index j. The setting above will calculate all 6 NN and 12 2nd NN interactions, while that below three of the NN.


#### Indices of j'th orbitals
```bash
nj=[0]                
```

#### lattice vectors connecting origins of unit cell containing j to unit cell containing i 
```bash
Rs=np.array([[-1,-3,0]])       #
```
#### Fermi energy (Check the OUTCAR file from VASP to get E-fermi)
```bash
Ef = 8.8341             
```

#### Minimum energy for energy integration: Integral calculated from Em to Ef
```bash
Em = -2                
```

#### nr of kpoints in each direction
```bash
kmesh=[12,12,16]     
```

#### Number of energy points of energy integral
```bash
n_Ec=1000              
```

####
```bash
symdec = 1
```

#### the orbital quantum numbers of the orbitals included in wannierization, i.e. 0=s, 1=p, 2=d etc. So here the first 2s are for the Mn d, and the remaining 1's are for O p (2 Mn and 6 O per doubled perovskite unit cell). If those N atoms/orbitals included do not correspond to the first N atoms in the unit cell, then tThe second column contains the corresponding atomic indices for the atoms, basically as listed in the POSCAR if starting from VASP. In this case, the first two atoms are Sr, which are excluded in wannierization. This information is strictly speaking not needed for calculating Jij, but is used for calculating interatomic vectors/distances using information read from wannier90.wout. 
```bash
l=[[4,4,4,4,4,4,4,4], range(0,8)]       
```


Authors:
First written by Xiangzhou Zhu in 2019. 
Modified and extended by Alexander Edström in 2020.
Further modified by Ankit Izardar Sept. 2020 - getsymdec function now includes calculation 
of J_ob_mix and J_ob separately. Deltas not diagonalized
