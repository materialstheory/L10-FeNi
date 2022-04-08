#!/usr/bin/env python3
# coding: utf-8


'''
Calculate magnetic exchange interaction with WANNIER90 output files wannier90.up_hr.dat and wannier90.dn_hr.dat
Information about structure read from wannier90.wout
First written by Xiangzhou Zhu in 2019. 
Modified and extended by Alexander Edström in 2020.
Further modified by Ankit Izardar Sept. 2020 -Deltas not diagonalized
The method is described in https://journals.aps.org/prb/abstract/10.1103/PhysRevB.101.064401 and references therein. 
'''

# In[33]:

from types import *
import numpy as np
import math
from itertools import product
from scipy.integrate import quad
from scipy.integrate import simps
import matplotlib.pyplot as plt
import sys
#import DM


# In[34]:


def read_wannier90hr(hr_filename="wannier_hr.dat"):
    """
    Function for reading the seedname_hr.dat file produced by Wannier90 (http://wannier.org)

    Parameters
    ----------
    hr_filename : string
         full name of the H(R) file produced by Wannier90 (usually seedname_hr.dat)

    Returns
    -------
    nrpt : integer
        number of R vectors found in the file
    rvec_idx : np.array of integers
        Miller indices of the R vectors
    rvec_deg : np.array of floats
        weight of the R vectors
    num_wf : integer
        number of Wannier functions found
    h_of_r : list of np.array
        <w_i|H(R)|w_j> = Hamilonian matrix elements in the Wannier basis

    """

    # Read only from the "wannier_hr.dat"
    try:
        with open(hr_filename, "r") as hr_filedesc:
            hr_data = hr_filedesc.readlines()
            hr_filedesc.close()
    except IOError:
        print("The file %s could not be read!" % hr_filename)

        print("Reading %s..." % hr_filename + hr_data[0])

    try:
        # reads number of Wannier functions per spin
        num_wf = int(hr_data[1])
        nrpt = int(hr_data[2])
    except ValueError:
        print("Could not read number of WFs or R vectors")

    # allocate arrays to save the R vector indexes and degeneracies and the
    # Hamiltonian
    rvec_idx = np.zeros((nrpt, 3), dtype=int)
    rvec_deg = np.zeros(nrpt, dtype=int)
    h_of_r = [np.zeros((num_wf, num_wf), dtype=np.complex_)
              for n in range(nrpt)]
    # variable currpos points to the current line in the file
    currpos = 2
    try:
        ir = 0
        # read the degeneracy of the R vectors (needed for the Fourier
        # transform)
        while ir < nrpt:
            currpos += 1
            for x in hr_data[currpos].split():
                if ir >= nrpt:
                    raise IndexError("wrong number of R vectors??")
                rvec_deg[ir] = int(x)
                ir += 1
        # for each direct lattice vector R read the block of the
        # Hamiltonian H(R)
        for ir, jj, ii in product(range(nrpt), range(num_wf), range(num_wf)):
            # advance one line, split the line into tokens
            currpos += 1
            cline = hr_data[currpos].split()
            # check if the orbital indexes in the file make sense
            if int(cline[3]) != ii + 1 or int(cline[4]) != jj + 1:
                mpi.report(
                    "Inconsistent indices at %s%s of R n. %s" % (ii, jj, ir))
            rcurr = np.array(
                [int(cline[0]), int(cline[1]), int(cline[2])])
            if ii == 0 and jj == 0:
                rvec_idx[ir] = rcurr
                rprec = rcurr
            else:
                # check if the vector indices are consistent
                if not np.array_equal(rcurr, rprec):
                    mpi.report(
                        "Inconsistent indices for R vector n. %s" % ir)

            # fill h_of_r with the matrix elements of the Hamiltonian
            h_of_r[ir][ii, jj] = complex(float(cline[5]), float(cline[6]))

    except ValueError:
        mpi.report("Wrong data or structure in file %s" % hr_filename)

    # return the data into variables
    return nrpt, rvec_idx, rvec_deg, num_wf, h_of_r


# In[35]:


def fourier_ham(norb, n_k, rvec_idx, rvec_deg, h_of_r, kmesh):
    """
    Function for obtaining H(k) from H(R) via Fourier transform
    The R vectors and k-point mesh are read from global module variables

    Parameters
    ----------
    norb : integer
        number of orbitals
    n_k : integer
        total number of k-points in the mesh
    h_of_r : list of np.array[norb,norb]
        Hamiltonian H(R) in Wannier basis
    rvec_idx : np.array of integers
        Miller indices of the R vectors
    rvec_deg : np.array of floats
        weight of the R vectors

    Returns
    -------
    h_of_k : list of np.array[norb,norb]
        transformed Hamiltonian H(k) in Wannier basis

    """

    twopi = 2 * np.pi
    h_of_k = [np.zeros((norb, norb), dtype=np.complex_)
              for ik in range(n_k)]
    ridx = np.array(range(nrpt))

    for ik, ir in product(range(n_k), ridx):
        rdotk = twopi * np.dot(kmesh[ik], rvec_idx[ir])

        factor = (math.cos(rdotk) + 1j * math.sin(rdotk)) / float(rvec_deg[ir])
        h_of_k[ik][:, :] += factor * h_of_r[ir][:, :]
    
    return h_of_k


# In[36]:


def kmesh_build(msize=None, mmode=0):
    """
    Function for the generation of the k-point mesh.
    Right now it only supports the option for generating a full grid containing k=0,0,0.

    Parameters
    ----------
    msize : list of 3 integers
        the dimensions of the mesh
    mmode : integer
        mesh generation mode (right now, only full grid available)

    Returns
    -------
    nkpt : integer
        total number of k-points in the mesh
    kmesh : np.array[nkpt,3] of floats
        the coordinates of all k-points
    wk : np.array[nkpt] of floats
        the weight of each k-point

    """

    if mmode == 0:
        # a regular mesh including Gamma point
        # total number of k-points
        nkpt = msize[0] * msize[1] * msize[2]
        kmesh = np.zeros((nkpt, 3), dtype=float)
        ii = 0
        for ix, iy, iz in product(range(msize[0]), range(msize[1]), range(msize[2])):
            kmesh[ii, :] = [float(ix) / msize[0], float(iy) /
                            msize[1], float(iz) / msize[2]]
            ii += 1
        # weight is equal for all k-points because wannier90 uses uniform grid on whole BZ
        # (normalization is always 1 and takes into account spin degeneracy
        wk = np.ones([nkpt], dtype=float) / float(nkpt)
    elif mmode == 1:
        # a regular MP mesh including Gamma point
        # total number of k-points
        nkpt = msize[0] * msize[1] * msize[2]
        kmesh = np.zeros((nkpt, 3), dtype=float)
        ii = 0
        for ix, iy, iz in product(range(msize[0]), range(msize[1]), range(msize[2])):
            kmesh[ii, :] = [(float(ix)*2 - msize[0] - 1) / msize[0] / 2, (float(iy)*2 - msize[1] - 1) /
                            msize[1] / 2, (float(iz)*2 - msize[2] - 1) / msize[2] / 2]
            ii += 1
        # weight is equal for all k-points because wannier90 uses uniform grid on whole BZ
        # (normalization is always 1 and takes into account spin degeneracy
        wk = np.ones([nkpt], dtype=float) / float(nkpt)
    else:
        raise ValueError("Mesh generation mode not supported: %s" % mmode)

    return nkpt, kmesh, wk


# In[37]:


def delta(h_of_k_up, h_of_k_dn, ni, nj, nl, nkpt):
    """
    Function for generation of exchange splittings for atoms on site i and j

    Parameters
    ----------
    h_of_k_up : list of np.array[norb,norb]
        transformed Hamiltonian H(k) for spin up in Wannier basis
    h_of_k_dn : list of np.array[norb,norb]
        transformed Hamiltonian H(k) for spin down in Wannier basis
    norb_i : integer  NOT USED! 
        number of orbitals for the atom on site i 
    nl : array 
        contains nl[ni] is new norb_i, where ni is index of atom i
    nkpt : integer
        total number of k-points in the mesh                DIFFERENCE BETWEEN nkpt AND n_k??? I removed n_k

    Returns
    -------
    delta_i : np.array[norb_i,norb_i]
        exchange splitting for the atom on site i
    delta_j : np.array[norb_j,norb_j]
        exchange splitting for the atom on site j
    """
    Mi = sum(nl[0:ni])
    Mj = sum(nl[0:nj])

    #print("Mi = ", Mi)
    #print("Mj = ", Mj)

    #print("nl = ", nl)
    #print("nl[ni] = ", nl[ni])
    
    # atom i and j have the same number of orbitals
    delta_i = np.zeros((nl[ni], nl[ni]), dtype=np.complex_)
    delta_j = np.zeros((nl[nj], nl[nj]), dtype=np.complex_)
    for ik in range(nkpt):
        delta_i += (h_of_k_up[ik][Mi:(Mi+nl[ni]), Mi:(Mi+nl[ni])] -
                    h_of_k_dn[ik][Mi:(Mi+nl[ni]), Mi:(Mi+nl[ni])])
        delta_j += (h_of_k_up[ik][Mj:(Mj+nl[nj]), Mj:(Mj+nl[nj])] -
                    h_of_k_dn[ik][Mj:(Mj+nl[nj]), Mj:(Mj+nl[nj])])
    #print("delta_i = ", delta_i)
    #print("delta_j = ", delta_j)
    return delta_i/nkpt, delta_j/nkpt


# In[38]:

def readstruct(fwoutstr):
    """
    Read information about the crystal structure from wannier90.wout

    Parameters
    ----------
    fwoutstr : name of the file to be read.

    Returns
    -------
    A : lattice vectors
    p : atomic positions in units of the above vectors in A
    """    

    with open(fwoutstr, 'r') as fwout:  
        lines = fwout.readlines()
        lcheck = 0
        p = []
        for li in range(len(lines)):
            # for lin in fwout.readlines():
            if 'a_1' in lines[li]:
                a1 = [float(i) for i in lines[li].split()[1:]]
            if 'a_2' in lines[li]:
                a2 = [float(i) for i in lines[li].split()[1:]]
            if 'a_3' in lines[li]:
                a3 = [float(i) for i in lines[li].split()[1:]]
            if "Site       Fractional Coordinate" in lines[li]:
                N1 = li
                lcheck = 1
            if lcheck and "*---------------------------------------------------" in lines[li]:
                N2 = li
                lcheck = 0
        for li in range(N1+2, N2):
            p.append([x for x in lines[li].split()[3:6]])
    p = [[float(x[0]), float(x[1]), float(x[2])] for x in p]
    A = np.array([a1, a2, a3])
    return A, p

# In[39]


def getdist(n1, n2, p, A, R):
    """
    Calculate distance between orbital n1 and n2 in unit cells separated by -R.

    Parameters
    ----------
    n1 : orbital index of atom i        n2 : orbital index of atom j
    p : atomic positions                A : lattice vectors
    R : lattice vector connecting unit cell of j to that of i

    Returns
    -------
    Interatomic distance dij=|Rij| and the vector Rij
    """    
    
    print ("n1, n2 = ", n1, n2)

    R12 = np.matmul(p[n2][:], A) - np.matmul(p[n1][:], A) - np.matmul(R, A)
    return np.linalg.norm(R12), R12

# In[40]

def get_symdec(GF_ob_up,GF_ob_dn,delta1,delta2,t, theta): 
    """
    Orbital decomposition (also works if delta is not diagonal)
    Currently delta not diagonal

    Parameters
    ----------
    GF_ob_up : spin up GF
    GF_ob_dn : spin down GF
    t : 

    Returns
    -------
    J_ob : Contribution to Jij from each orbital
    """      

    """
    General form  --> delta1(h,k)*G_dn(k,l)*delta2(l,m)*G_up(m,k)
    Start by defining delta1 and delta2 and then G_up and G_dn based on the indices. 
    """


    GN=np.shape(GF_ob_up)
    J_ob = np.zeros((GN[1], GN[2]))
    J_int_ob = np.zeros(len(E), dtype= np.complex_)
    # print("J_int type = ", J_int.dtype)
    J_int_ob_mix = np.zeros(len(E), dtype=np.complex_)



    for x in range(GN[1]):
        for y in range(GN[2]):
            for ie in range(len(t)):
                J_int_ob[ie] = delta1[x, x]*GF_ob_dn[ie][x, y]*delta2[y, y]*GF_ob_up[ie][y, x]*1j*t[ie]

            J_ob[x, y] = -np.imag(simps(J_int_ob, theta)) / twopi 


    for ie in range(len(t)):
        
        J_int_ob_mix[ie] = (delta1[0, 4]*GF_ob_dn[ie][4,0]*delta2[0, 4]*GF_ob_up[ie][4, 0] + \
                            delta1[0, 4]*GF_ob_dn[ie][4,4]*delta2[4, 0]*GF_ob_up[ie][0, 0] + \
                            delta1[4, 0]*GF_ob_dn[ie][0,4]*delta2[4, 0]*GF_ob_up[ie][0, 4] + \
                            delta1[4, 0]*GF_ob_dn[ie][0,0]*delta2[0, 4]*GF_ob_up[ie][4, 4])*1j*t[ie]


    for k in range(GN[1]):
        for l in range(GN[2]):
            for m in range(GN[2]):
                if m !=l:
                    for ie in range(len(t)):

                        J_int_ob_mix[ie] += (delta1[k, k]*GF_ob_dn[ie][k,l]*delta2[l, m]*GF_ob_up[ie][m, k] + \
                                            delta1[l, m]*GF_ob_dn[ie][m,k]*delta2[k, k]*GF_ob_up[ie][k, l])*1j*t[ie]




    
    J_ob_mix = -np.imag(simps(J_int_ob_mix, theta)) / twopi

    
    return J_ob, J_ob_mix
    # d orbital contribution
    #print('eg', (J_ob[0, 0]+J_ob[3, 3]+J_ob[0, 3]+J_ob[3, 0]).imag)
    #print('t2g', (J_ob[1, 1]+J_ob[2, 2]+J_ob[4, 4] +
    #              (J_ob[1, 2]+J_ob[1, 4]+J_ob[2, 4])*2).imag)
    #print('eg-t2g', (J_ob[0, 1]+J_ob[0, 2]+J_ob[0, 4] +
    #                 J_ob[3, 1]+J_ob[3, 2]+J_ob[3, 4]).imag*2)


# In[41]:

def get_Jij_of_E(Em, Emax, n_E, symdec, ni, nj, Mi, Mj, R,h_of_k_up, h_of_k_dn, nkpt):
    """
    Energy resolved Jij(E), calculated along E + i eta in the energy plane.

    Parameters
    ----------
    Em : minimum energy
    Emax : maximum energy 
    n_E : nr of energy points
    symdec: print symmetry decomposed or not.
    ni : dimension of i:th orbital subspace         nj : dimension of j:th orbital subspace
    Mi : index where i orbitals start               Mj : index where j orbitals start 
    R : lattice vector from j:th unit cell to i:th
    h_of_k_up : spin up Hamiltonian                 h_of_k_dn : spin down Hamiltonian
    nkpt : nr of k-points

    Returns
    -------
    J_int : Jij(E)
    """        
    E = np.linspace(Em, Emax, n_E)

    # J of energy
    if symdec: 
        J_int = np.zeros((n_E,ni,nj), dtype=np.complex_)
    elif not symdec:
        J_int = np.zeros(n_E, dtype=np.complex_)


    # factors
    factor = np.zeros(nkpt, dtype=np.complex_)
    factor2 = np.zeros(nkpt, dtype=np.complex_)
    for ik in range(nkpt):
        rdotk = twopi*np.dot(k_mesh[ik], R)
        factor[ik] = (math.cos(rdotk) + 1j * math.sin(rdotk))
        factor2[ik] = (math.cos(rdotk) - 1j * math.sin(rdotk))

    # integrate from Em+i*eta to Emax+i*eta
    for ie in range(n_E):
        # define G(k) for diff k and diff E
        GF_k_up = [np.zeros((norb, norb), dtype=np.complex_)
                   for ik in range(nkpt)]
        GF_k_dn = [np.zeros((norb, norb), dtype=np.complex_)
                   for ik in range(nkpt)]

        for ik in range(nkpt):
            GF_k_up[ik] = np.linalg.inv((E[ie]+1j*eta)*I-h_of_k_up[ik])
            GF_k_dn[ik] = np.linalg.inv((E[ie]+1j*eta)*I-h_of_k_dn[ik])

        # integration over Brilllouin zone for E
        GF_ij_up = np.zeros((ni, ni), dtype=np.complex_)
        GF_ij_dn = np.zeros((nj, nj), dtype=np.complex_)

        for ik in range(nkpt):
            # i' and j' are different site in the unit cell
            GF_ij_up += factor2[ik] * \
                GF_k_up[ik][Mj:(Mj+nj), Mi:(Mi+ni)]/nkpt
            GF_ij_dn += factor[ik] * \
                GF_k_dn[ik][Mi:(Mi+ni),  Mj:(Mj+nj)]/nkpt 

        if symdec: 
            J_int[ie] = (
            1/(twopi)*(np.dot(np.dot(np.dot(delta1, GF_ij_dn), delta2), GF_ij_up))).imag
        elif not symdec: 
            J_int[ie] = (
            1/(twopi)*np.trace(np.dot(np.dot(np.dot(delta1, GF_ij_dn), delta2), GF_ij_up))).imag
    return real(J_int)


def get_Gij(GF_k_up, GF_k_dn, nli, nlj, Mi, Mj, R, k_mesh):
    """
    Evaluate Gij_dn and Gji_up

    Parameters
    ----------
    ni : dimension of i:th orbital subspace         nj : dimension of j:th orbital subspace
    Mi : index where i orbitals start               Mj : index where j orbitals start 
    R : lattice vector from j:th unit cell to i:th
    h_of_k_up : spin up Hamiltonian                 h_of_k_dn : spin down Hamiltonian
    nkpt : nr of k-points

    Returns
    -------
    GF_ji_up : Jij(E)
    GF_ij_dn
    """  
    nkpt=np.shape(k_mesh)[0]
    n_Ec=np.shape(GF_k_up)[0]
    GF_ji_up = np.zeros((n_Ec,nlj, nli), dtype=np.complex_)
    GF_ij_dn = np.zeros((n_Ec,nlj, nli), dtype=np.complex_)
    for ik in range(nkpt):
        rdotk = twopi*np.dot(k_mesh[ik], R)
        factor = (math.cos(rdotk) + 1j * math.sin(rdotk))
        factor2 = (math.cos(rdotk) - 1j * math.sin(rdotk))
        GF_ji_up[:] += factor2 * \
            GF_k_up[:,ik,Mj:(Mj+nlj), Mi:(Mi+nli)]/nkpt
        GF_ij_dn[:] += factor * \
            GF_k_dn[:,ik,Mi:(Mi+nli),  Mj:(Mj+nlj)]/nkpt
    return GF_ji_up, GF_ij_dn

def find_kth_atoms(Rsoc,ni, nj, A, atpos,l,R):
    """
    Find positions k closer than Rsoc to i and j,
    over which summation is performed when calculating Dij.

    Parameters
    ----------
    Rsoc : cut off distance for summation
    ni : index for atom i       nj : index for atom j 
    A : lattice vectors 
    atpos : atomic positions
    l: orbitals l[0][:] and their corresponding atomic indices l[1][:]

    Returns
    -------
    Rik : vector connecting unit cell connecting i to that containing k
    nk : index of atom k
    """     

    Rik = []
    nk = []
    lattd = [np.linalg.norm(A[0][:]), np.linalg.norm(
        A[1][:]), np.linalg.norm(A[2][:])]
    for ix, iy, iz in product(range(-1-int(Rsoc//lattd[0]), 2+int(Rsoc//lattd[0])), range(-1-int(Rsoc//lattd[1]), 2+int(Rsoc//lattd[1])), range(-1-int(Rsoc//lattd[2]), 2+int(Rsoc//lattd[2]))):
        for ik in range(len(l[0])):
            atd, Rvec = getdist(l[1][ni], l[1][ik], atpos, A, [-ix, -iy, -iz])
            if atd <= Rsoc:
                Rik.append([ix, iy, iz])
                nk.append(ik)
                #Rij.append(Rvec)
                #dij.append(atd)
    for ix, iy, iz in product(range(-1-int(Rsoc//lattd[0]), 2+int(Rsoc//lattd[0])), range(-1-int(Rsoc//lattd[1]), 2+int(Rsoc//lattd[1])), range(-1-int(Rsoc//lattd[2]), 2+int(Rsoc//lattd[2]))):
        for ik in range(len(l[0])):
            Ri=[ix - R[0], iy - R[1], iz - R[2]]
            atdj, Rvecj = getdist(l[1][nj], l[1][ik], atpos, A, [-ix, -iy, -iz])
            atdi, Rveci = getdist(l[1][ni], l[1][ik], atpos, A, [-r for r in Ri])
            if (atdj <= Rsoc) and (atdi > Rsoc):
                Rik.append(Ri)
                nk.append(ik)
                #Rij.append(Rvec)
                #dij.append(atd)               
    return Rik, nk              

# def get_Hsoc(th,ph):
#     """
#     Return the spin diagonal parts of the SOC Hamiltonian, Huu and Hdd, in basis of d-orbitals, ordered: 
#     dz2 , dxz , dyz , dx2-y2 , dxy

#     Parameters
#     ----------
#     th : spherical coordinate theta of magnetization
#     ph : spherical coordinate phi of magnetization

#     Returns
#     -------
#     Huu : Hsoc matrix coupling spin up to spin up
#     Hdd : Hsoc matrix coupling spin dn to spin dn
#     """       

#     Huu = np.array([ [      0,                   -np.sqrt(3)*1j*np.sin(th)*np.sin(ph)/2 ,  np.sqrt(3)*1j*np.sin(th)*np.cos(ph)/2,   0,                          0] , 
#                      [ np.sqrt(3)*1j*np.sin(th)*np.sin(ph)/2 ,     0 ,                    -1j*np.cos(th)/2 ,                -1j*np.sin(th)*np.sin(ph)/2,     1j*np.sin(th)*np.cos(ph)/2], 
#                      [-np.sqrt(3)*1j*np.sin(th)*np.cos(ph)/2 ,1j*np.cos(th)/2 ,                0 ,                          -1j*np.sin(th)*np.cos(ph)/2,    -1j*np.sin(th)*np.sin(ph)/2] ,
#                      [      0,                                1j*np.sin(th)*np.sin(ph)/2 ,   1j*np.sin(th)*np.cos(ph)/2 ,               0,                    -1j*np.cos(th)] ,
#                      [      0,                               -1j*np.sin(th)*np.cos(ph)/2 ,   1j*np.sin(th)*np.sin(ph)/2 ,      1j*np.cos(th) ,                 0] ])

#     Hdd = Huu
#     # Hdd = Huu.conj().T # = Huu, Hermitian.

#     return Huu, Hdd

#def eval_DM():
# def eval_DM(delta1, delta2, GF_k_up, GF_k_dn, xi, Rsoc, Em, Ef, n_Ec, ni, nj, A, atpos,l,nl,R,k_mesh):    # ENERGIES NOT NEEDED??
#     """
#     Evaluate DM interaction Dij 

#     Parameters
#     ----------
#     delta1, delta2, GF_k_up, GF_k_dn, xi, Rsoc, Em, Ef, n_Ec, ni, nj[Rn], A, atpos,l,R

#     Returns
#     -------
#     Dij : 
#     """   

#     Rik, nk = find_kth_atoms(Rsoc,ni, nj, A, atpos,l,R)  

#     Huux, Hddx = get_Hsoc(np.pi/2,0)
#     Huuy, Hddy = get_Hsoc(np.pi/2,np.pi/2)
#     Huuz, Hddz = get_Hsoc(0,0)

#     Gji_u, Gij_d=get_Gij(GF_k_up, GF_k_dn, nl[ni], nl[nj], Mi, Mj, R, k_mesh)
#     Gij_u, Gji_d=get_Gij(GF_k_up, GF_k_dn, nl[nj], nl[ni], Mj, Mi, [-x for x in R], k_mesh)
#     Dij=np.zeros((3,n_Ec), dtype=np.complex_)
#     for k in range(len(Rik)): 
#         if l[0][nk[k]] == 2:
#             Rjk=np.array(Rik[k]) + np.array(R)
#             Gki_u, Gik_d=get_Gij(GF_k_up, GF_k_dn, nl[ni], nl[nk[k]], Mi, sum(nl[0:nk[k]]), [-x for x in Rik[k]], k_mesh)
#             Gik_u, Gki_d=get_Gij(GF_k_up, GF_k_dn, nl[nk[k]], nl[ni], sum(nl[0:nk[k]]), Mi, Rik[k], k_mesh)
#             Gkj_u, Gjk_d=get_Gij(GF_k_up, GF_k_dn, nl[nj], nl[nk[k]], Mj, sum(nl[0:nk[k]]), [-x for x in Rjk], k_mesh)
#             Gjk_u, Gkj_d=get_Gij(GF_k_up, GF_k_dn, nl[nk[k]], nl[nj], sum(nl[0:nk[k]]), Mj, Rjk, k_mesh)        

            
#             Dij[0] += xi[nk[k]]*np.trace(np.matmul( delta1, np.matmul( Gik_d, np.matmul( Hddx ,  np.matmul( Gkj_d, np.matmul( delta2 , Gji_u))))), axis1=1, axis2=2 )
#             Dij[0] -= xi[nk[k]]*np.trace(np.matmul( delta1, np.matmul( Gik_u, np.matmul( Huux ,  np.matmul( Gkj_u, np.matmul( delta2 , Gji_d))))), axis1=1, axis2=2 )
#             Dij[0] += xi[nk[k]]*np.trace(np.matmul( delta1, np.matmul( Gij_d, np.matmul( delta2 , np.matmul( Gjk_u, np.matmul( Huux , Gki_u))))), axis1=1, axis2=2 )
#             Dij[0] -= xi[nk[k]]*np.trace(np.matmul( delta1, np.matmul( Gij_u, np.matmul( delta2 , np.matmul( Gjk_d, np.matmul( Hddx , Gki_d))))), axis1=1, axis2=2 )

#             Dij[1] += xi[nk[k]]*np.trace(np.matmul( delta1, np.matmul( Gik_d, np.matmul( Hddy ,  np.matmul( Gkj_d, np.matmul( delta2 , Gji_u))))), axis1=1, axis2=2 )
#             Dij[1] -= xi[nk[k]]*np.trace(np.matmul( delta1, np.matmul( Gik_u, np.matmul( Huuy ,  np.matmul( Gkj_u, np.matmul( delta2 , Gji_d))))), axis1=1, axis2=2 )
#             Dij[1] += xi[nk[k]]*np.trace(np.matmul( delta1, np.matmul( Gij_d, np.matmul( delta2 , np.matmul( Gjk_u, np.matmul( Huuy , Gki_u))))), axis1=1, axis2=2 )
#             Dij[1] -= xi[nk[k]]*np.trace(np.matmul( delta1, np.matmul( Gij_u, np.matmul( delta2 , np.matmul( Gjk_d, np.matmul( Hddy , Gki_d))))), axis1=1, axis2=2 )

#             Dij[2] += xi[nk[k]]*np.trace(np.matmul( delta1, np.matmul( Gik_d, np.matmul( Hddz ,  np.matmul( Gkj_d, np.matmul( delta2 , Gji_u))))), axis1=1, axis2=2 )
#             Dij[2] -= xi[nk[k]]*np.trace(np.matmul( delta1, np.matmul( Gik_u, np.matmul( Huuz ,  np.matmul( Gkj_u, np.matmul( delta2 , Gji_d))))), axis1=1, axis2=2 )
#             Dij[2] += xi[nk[k]]*np.trace(np.matmul( delta1, np.matmul( Gij_d, np.matmul( delta2 , np.matmul( Gjk_u, np.matmul( Huuz , Gki_u))))), axis1=1, axis2=2 )
#             Dij[2] -= xi[nk[k]]*np.trace(np.matmul( delta1, np.matmul( Gij_u, np.matmul( delta2 , np.matmul( Gjk_d, np.matmul( Hddz , Gki_d))))), axis1=1, axis2=2 )

#     return Dij




# START OF MAIN PROGRAM

"""
Input prameters
--------------------------

norb: total number of orbitals in wannier 
    For SrMnO3 with lattice vectors (0,a,a), (a,0,a), (a,a,0)
    Consider Mn d and O p orbitals, total orbitals are 2*(5d+3*3p)=28 (for each spin)

norb_i: number of orbitals for atom on site i 
        for Mn, there are 5 d orbitals
        NOT USED! INSTEAD SPECIFY l ARRAY! 
    
Rs: Interatomic lattice vectors used in Fourier transformation of Green's function
    R = (R_i-R_i')-(R_j-R_j'), where R_i, R_j are position of Mn, R_i', R_j' are positions of 
    i:th and j:th atoms in the unit cell
    For first nearest neighbor in SrMnO3, R = [1,0,0], second R = [1,-1,-1], ...

ni: index of orbital i 

nj: index/indices of orbital j 
    
Ef: Fermi energy

Em: Lower limit of energy integration 

kmesh:  of k-point mesh

kmode: Specifies how k-mesh is generated. 0 or 1. 

n_Ec: size of energy mesh in semicircle path

es: distance between the path and real axis
    can be chosen as small as possible
Rmax: Maximum interatomic distance in Angstrom. Either specify this or set it to zero and list nj, Rs. 
"""

"""
INPUT
"""
# Set default values of input parameters
#norb = 5                    # Total basis size, calculated below
ni = 0                      # Index of orbitals of atom i   (We calculate interaction Jij between i and j)
nj = [0]                    # Index/indices of orbitals on atom j 
Rs = np.array([1, 0, 0]) 
Rmax = -5.0
Ef = 0.0000
Em = -8
kmesh = [1, 1, 1]
n_Ec = 3000
# es = 0.0000
l = [2]
kmode = 0
symdec = 1
outname = "Jij.out"

# Input for calculating DM interactions
evalDM=0        # Calculate DM interactions or not? 
Rsoc=5.0        # Cut-off distance [Å], for SOC summation 
xi=[0.0]        # SOC constants in eV, 

# Rectangular path input, for plotting Jij(E) with imaginary energy E + i eta
E_res = 0
if E_res:
    n_E = 1000
    Emax = Ef
    eta = 0.1

inputf = open(sys.argv[1], 'r')
for lin in inputf.readlines():
    exec(lin)
inputf.close()

# If not specified, add atomic indices of the orbitals
if len(np.shape(l)) == 1:
    l=[l,range(len(l))]


nl = [2*x+1 for x in l[0]]
print("nl = " ,nl)
Mi = sum(nl[0:ni])
#print("Mi = ",Mi)
norb = sum(nl)
#print("norb =", norb)
if (len(nj) != np.shape(Rs)[0]) and Rmax <= 0.0:    # Check that the lengths of nj and Rs are same.
    raise Exception('Number of j and R do not coincide!')

"""
END INPUT
"""

I = np.eye(norb, norb)
twopi = 2 * np.pi

# In[42]:

# Generate k-mesh
nkpt, k_mesh, wk = kmesh_build(kmesh, kmode)
# read Hamiltonian spin up
nrpt, rvec_idx, rvec_deg, num_wf, h_of_r_up = read_wannier90hr("wannier90.up_hr.dat")
#print("num_wf = ", num_wf)
# if num_wf != norb: 
#     raise Exception('Number of orbitals considered do not match the number of Wannier functions ' +
#                     str(norb) + ' vs ' + str(num_wf))
# Fourier transform Hamiltonian
h_of_k_up = fourier_ham(norb, nkpt, rvec_idx, rvec_deg, h_of_r_up, k_mesh)




# Hamiltonian spin down
nrpt, rvec_idx, rvec_deg, num_wf, h_of_r_dn = read_wannier90hr("wannier90.dn_hr.dat")
h_of_k_dn = fourier_ham(norb, nkpt, rvec_idx, rvec_deg, h_of_r_dn, k_mesh)


'''
Setup semi-circular path over which integration in the complex plane is done
'''

# radius and center of the circle
r = (Ef-Em)/2
C0 = (Ef+Em)/2

# Energy mesh along the integral path
theta = np.linspace(0, np.pi, n_Ec)
E = np.zeros((n_Ec), dtype=np.complex_)
t = np.zeros((n_Ec), dtype=np.complex_)
for it in range(n_Ec):
    E[it] = math.cos(theta[it])*r+C0+1j*math.sin(theta[it])*r
    t[it] = E[it]-C0

# define J as function of E/theta
J_int = np.zeros(n_Ec, dtype=np.complex_)
if evalDM:
    Dij_E = np.zeros((n_Ec, 3), dtype=np.complex_)


# If the Rij, nj considered are chosen by Rmax, they are determined here.
# The interatomic vectors and distances are stored in Rij and dij=|Rij|.
A, atpos = readstruct('wannier90.wout')
Rij = []
dij = []
if Rmax > 0.0:
    Rs = []
    njs=nj
    nj = []
    #nj, Rs, dij, Rij = make_R(A, atpos, Rmax, ni, njs)
    lattd = [np.linalg.norm(A[0][:]), np.linalg.norm(
        A[1][:]), np.linalg.norm(A[2][:])]
    for ix, iy, iz in product(range(-1-int(Rmax//lattd[0]), 2+int(Rmax//lattd[0])), range(-1-int(Rmax//lattd[1]), 2+int(Rmax//lattd[1])), range(-1-int(Rmax//lattd[2]), 2+int(Rmax//lattd[2]))):
        for nnj in njs:
            atd, Rvec = getdist(l[1][ni], l[1][nnj], atpos, A, [ix, iy, iz])
            if atd <= Rmax and atd != 0:
                Rs.append([ix, iy, iz])
                nj.append(nnj)
                Rij.append(Rvec)
                dij.append(atd)
else:
    for Rn in range(np.shape(Rs)[0]):
        atd, Rvec = getdist(l[1][ni], l[1][nj[Rn]], atpos, A, Rs[Rn])
        Rij.append(Rvec)
        dij.append(atd)

# define G(k) for chosen k and E
GF_k_up = np.zeros((n_Ec, nkpt, norb, norb), dtype=np.complex_)
GF_k_dn = np.zeros((n_Ec, nkpt, norb, norb), dtype=np.complex_)
for ie in range(n_Ec):
    # G(k)=(E+i*eta-H(k))^-1
    for ik in range(nkpt):
        GF_k_up[ie,ik] = np.linalg.inv((E[ie])*I-h_of_k_up[ik])
        GF_k_dn[ie,ik] = np.linalg.inv((E[ie])*I-h_of_k_dn[ik])

outf = open(outname, 'w')
# Loop over different R vectors to calculate the different Jij
for Rn in range(np.shape(Rs)[0]):
    R = Rs[Rn]
    Mj = sum(nl[0:nj[Rn]])

    # exchange splitting
    delta1, delta2 = delta(h_of_k_up, h_of_k_dn, ni, nj[Rn], nl, nkpt)
    GF_ji_up, GF_ij_dn = get_Gij(GF_k_up, GF_k_dn, nl[ni], nl[nj[Rn]], Mi, Mj, R, k_mesh)    

    # To get d-bands contribution to J's (Total bands --> s(1), p(3), d(5))
    delta1_sp = delta1[:3,:3] 
    delta2_sp = delta2[:3,:3]

    delta1_d = delta1[4:,4:] 
    delta2_d = delta2[4:,4:]


    GF_ji_up_sp = GF_ji_up[:, :3, :3]
    GF_ij_dn_sp = GF_ij_dn[:, :3, :3]
 

    GF_ji_up_d = GF_ji_up[:, 4:, 4:]
    GF_ij_dn_d = GF_ij_dn[:, 4:, 4:]
 
    # Evaluate Jij(E) ( and Dij(E) )
   
   
    J_int_sp = np.trace(np.matmul(np.matmul(np.matmul(delta1_sp, GF_ij_dn_sp), delta2_sp), GF_ji_up_sp), axis1=1, axis2=2 )*1j*t
    J_int_d = np.trace(np.matmul(np.matmul(np.matmul(delta1_d, GF_ij_dn_d), delta2_d), GF_ji_up_d), axis1=1, axis2=2 )*1j*t
    J_int = np.trace(np.matmul(np.matmul(np.matmul(delta1, GF_ij_dn), delta2), GF_ji_up), axis1=1, axis2=2 )*1j*t
    # print ("J_int = ", J_int)
    # if evalDM: 
    #     Dij_E=eval_DM(delta1, delta2, GF_k_up, GF_k_dn, xi, Rsoc, Em, Ef, n_Ec, ni, nj[Rn], A, atpos,l, nl,R, k_mesh)*1j*t     

    # Integration over energy
    # add minus sign with AFM reference state
    J_sp = -np.imag(simps(J_int_sp, theta)) / (twopi)
    J_d = -np.imag(simps(J_int_d, theta)) / (twopi)
    J = -np.imag(simps(J_int, theta)) / (twopi)
    # print ("J = ", J)
    # if evalDM: 
    #     Dij = -np.real(simps(np.real(Dij_E), theta) ) / twopi


    # Below, results are printed. A summary to stdout and then most results are printed in Jij.out
    if not evalDM:
        print("Rij = " + '{:08.6f}'.format(Rij[Rn][0]) + "  " + '{:08.6f}'.format(Rij[Rn][1]) + "  " '{:08.6f}'.format(Rij[Rn][2]) +
          ",  |Rij| = " + '{:08.6f}'.format(dij[Rn]) + ",  Jij(sp) = " + '{:10.8f}'.format(J_sp*1000) + "meV" +",  Jij(d) = " + '{:10.8f}'.format(J_d*1000) + " meV" + ",  Jij(total) = " + '{:10.8f}'.format(J*1000) + " meV")
    elif evalDM:
        print("Rij = " + '{:08.6f}'.format(Rij[Rn][0]) + "  " + '{:08.6f}'.format(Rij[Rn][1]) + "  " '{:08.6f}'.format(Rij[Rn][2]) +
          ",  |Rij| = " + '{:08.6f}'.format(dij[Rn]) + ",  Jij = " + '{:10.8f}'.format(J*1000) + " meV" + ",  Dij = (" + '{:10.8f}'.format(Dij[0]*1000) + ", " + '{:10.8f}'.format(Dij[1]*1000) + ", " + '{:10.8f}'.format(Dij[2]*1000) + ") meV")

    # Print to Jij.out
    outstr = str(ni) + "  " + str(nj[Rn]) + "    " + str(R[0]) + "  " + str(R[1]) + "  " + str(R[2]) + "    " + '{:08.6f}'.format( 
        Rij[Rn][0]) + "  " + '{:08.6f}'.format(Rij[Rn][1]) + "  " + '{:08.6f}'.format(Rij[Rn][2]) + "   " + '{:08.6f}'.format(dij[Rn]) + "  " + '{:10.8f}'.format(J*1000) 
    if evalDM:
        outstr = outstr + "    " + '{:10.8f}'.format(Dij[0]*1000) + " " + '{:10.8f}'.format(Dij[1]*1000) + " " + '{:10.8f}'.format(Dij[2]*1000) 
    if symdec:
        # if np.count_nonzero([round(np.absolute(x),5) for x in np.reshape(delta1 - np.diag(np.diagonal(delta1)),np.size(delta1))]) != 0 or np.count_nonzero([round(np.absolute(x),5) for x in np.reshape(delta2 - np.diag(np.diagonal(delta2)),np.size(delta2))])!=0:
        #   raise Exception('Require diagonal exchange splitting for symdec.')
        J_ob, J_ob_mix = get_symdec(GF_ji_up, GF_ij_dn,delta1,delta2,t,theta)
        J_ob_sp = J_ob[:3,:3]
        J_ob_d = J_ob[4:,4:]
        print("J_ob = ", J_ob)
        print("J_ob sum (meV) = ", np.sum(J_ob) * 1000)
        print("J_ob sum_d_only (meV) = ", np.sum(J_ob_d) * 1000)
        print("J_ob sum_sp_only (meV) = ", np.sum(J_ob_sp) * 1000)
        # print("J_ob sum of diagonal elements (meV) = ", np.trace(J_ob) * 1000)
        print("J_ob_mix (meV) = ", J_ob_mix * 1000)
        print("Total (J_ob + J_ob_mix) ", np.sum(J_ob) * 1000 + J_ob_mix * 1000)
    #     Jorbstr=" "
    #     for ii, ij in product(range(np.shape(J_ob)[0]),range(np.shape(J_ob)[1])):
    #       Jorbstr = Jorbstr + "   " + '{:12.10f}'.format(J_ob[ii,ij]*1000) 
    #     outstr = outstr + Jorbstr 
    # outstr = outstr + "\n"
    # outf.write(outstr)        

    # if E_res: # Energy resolved Jij printed to E_Jij.dat 
    #     Print("WARNING: The function printing Jij as function of energy has not been carefully tested!")
    #     Earr=np.linspace(Em, Emax, n_E)
    #     Eresdat=Earr
    #     for Rn in range(np.shape(Rs)[0]):
    #         J_E = get_Jij_of_E(Em, Emax, n_E, symdec, nl[ni], nl[nj[Rn]],Mi, Mj, R,h_of_k_up, h_of_k_dn, nkpt)
    #         if symdec: 
    #             Eresdat=np.vstack((Eresdat, J_E.reshape(np.shape(J_E)[0],np.shape(J_E)[1]*np.shape(J_E)[2]).T))
    #         elif not symdec:
    #             Eresdat=np.vstack((Eresdat, J_E))
    #     np.savetxt('E_Jij.dat', Eresdat.T)


    

