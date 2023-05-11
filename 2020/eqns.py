"""
2020 Modelling and Visualisation Past Paper
"""

import matplotlib
matplotlib.use('TKAgg')

import sys
import math
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import pandas as pd

cmap = plt.colormaps['gnuplot']

nstep=1000
lx=int(sys.argv[1])
ly=lx
p=float(sys.argv[2]) #probability S -> I
r=1000 #number of resamplings

def nearest_neighbours(lattice, i, j):

    """
    Returns array of nearest neighbour statuses (-1=I, 0=S), given an indexed position 
    in the lattice. Periodic boundary conditions are applied.
    
    :param lattice: NxN array of sites
    :param i: index i
    :param j: index j
    :return: 1x4 array containing the values of the nearest neighbours of
    lattice[i,j]
    """
    
    lx=len(lattice[0])

    NN_0=lattice[np.mod(i-1, lx), j]

    NN_1=lattice[i, np.mod(j-1, lx)]
        
    NN_2=lattice[np.mod(i+1, lx), j]
        
    NN_3=lattice[i, np.mod(j+1, lx)]
        
    NN_array=np.array([NN_0, NN_1, NN_2, NN_3])
        
    return NN_array

def NN_update(site, p):
    """
    Infects a given site with probability p
    """
    r=random.random()
    if r<p:
        site=-1

    return site

def I_test(lattice, i, j, p):
    """
    Returns updated status of an infected site and its nearest neighbours
    
    :param site: selected site from NxN lattice
    :param p: probability that an infected site will recover
    :return: updated status of chosen site and its updated nearest neighbour
    """

    r1=random.random()
    if r1>p:
        lattice[i,j]=0
       
    #select random nearest neighbour
    r2=random.randint(0,3)

    if r2==0:
        lattice[np.mod(i-1, lx), j]=NN_update(lattice[np.mod(i-1, lx), j], p)

    elif r2==1:
        lattice[i, np.mod(j-1, lx)]=NN_update(lattice[i, np.mod(j-1, lx)], p)

    elif r2==2:
        lattice[np.mod(i+1, lx), j]=NN_update(lattice[np.mod(i+1, lx), j], p)

    elif r2==3:
        lattice[i, np.mod(j+1, lx)]=NN_update(lattice[i, np.mod(j+1, lx)], p)

    return lattice

def init_random(lx):
    """
    Returns lx*lx lattice of randomly generated infected and 
    susceptible sites. 
    
    :param lx: axis size of desired array
    :return: randomly generated lattice
    """        
    ly=lx

    lattice=np.zeros((lx,ly), dtype=float) #NxN array to hold sites

    #initialise lattice randomly
    for i in range(lx):
        for j in range(ly):
            x=random.random()
            if(x<=(0.5)): lattice[i,j]=-1 #Infected
            if(x>(0.5)): lattice[i,j]=0 #Susceptible

    return lattice

def update_lattice(lattice, p):
    """
    Returns updated lattice based on infection/recovery rules
    """

    #select lattice randomly
    itrial=np.random.randint(0,len(lattice[0]))
    jtrial=np.random.randint(0,len(lattice[0]))

    site=lattice[itrial, jtrial]

    if site==-1:
        lattice=I_test(lattice, itrial, jtrial, p)

    return lattice

def number_I(lattice):
    """
    Returns the fraction of infected sites for a given lattice
    """
    lx=len(lattice[0])
    N=lx*lx

    n_infected=0
    for x in range(lx):
        for y in range(lx):
            n_infected+=abs(lattice[x,y])

    return n_infected
