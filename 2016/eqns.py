import sys
import math
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import pandas as pd

J=-1.0 #set constant J to -1
k_B=1.0 #Boltzmann constant
nstep=10000
lx=50

def nearest_neighbours(spin, i, j):

    """
    Returns array of nearest neighbour spins, given an indexed position in the 
    spin lattice. Periodic boundary conditions are applied.
    
    :param spin: NxN array of spins
    :param i: index i
    :param j: index j
    :return: 1x4 array containing the spin values of the nearest neighbours of
    spin[i,j]
    """
    
    N=len(spin[0])-1
    
    lx=len(spin[0])

    NN_0=spin[np.mod(i-1, lx), j]

    NN_1=spin[i, np.mod(j-1, lx)]
        
    NN_2=spin[np.mod(i+1, lx), j]
        
    NN_3=spin[i, np.mod(j+1, lx)]
        
    NN_array=np.array([NN_0, NN_1, NN_2, NN_3])
            
    return NN_array

def delta_E(spin, i, j, J, h):
    
    """
    Computes change in energy which arises from flipping the spin at position
    [i, j] in the NxN array of spins.
    
    :param spin: NxN array of spins
    :param i: index i
    :param j: index j
    :param J: magnetism constant (set to -1)
    :return: float giving the change in energy due to flipping of spin[i,j]
    """

    #generate array of nearest neighbour spins
    NN_spins=nearest_neighbours(spin, i, j)
    
    #change in energy due to flipping spin [i,j]
    E_change=2*J*spin[i,j]*np.sum(NN_spins)+2*h*spin[i,j]
    
    return E_change

def Glauber(spin, T, h):

    """
    Applies Glauber dynamics to a randomly selected spin and performs the Metropolis 
    algorithm.
    
    :param spin: NxN array of spins
    :param T: temperature
    :return: updated spin array with Glauber dynamics and Metropolis algorithm applied
    """

    #select spin randomly
    itrial=np.random.randint(0,len(spin[0]))
    jtrial=np.random.randint(0,len(spin[0]))
    spin_new=-spin[itrial,jtrial] #flipped trial spin at selected site

    #compute delta E eg via function (PBC applied)
    energy_change=delta_E(spin, itrial, jtrial, J, h) 

    #perform metropolis test

    if energy_change<=0:
        spin[itrial,jtrial]=spin_new 
    
    else:
        test=random.random() 
        boltz_factor=math.exp(-energy_change/(k_B*T))  
        
        if test<boltz_factor:
            spin[itrial,jtrial]=spin_new

    return spin

def init_random(lx):
    """
    Initialises lx*lx lattice with randomly generated spins (up=+1, down=-1)
    """
    spin=np.zeros((lx,lx),dtype=float) #NxN array to hold spins

    #initialise spins randomly
    for i in range(lx):
        for j in range(lx):
            x=random.random()
            if(x<0.5): spin[i,j]=-1 #down
            if(x>=0.5): spin[i,j]=1 #up
    
    return spin

def magnetisation(spin):

    """
    Calculates the absolute value of the total magnetisation of the spin array.
    
    :param spin: NxN array of spins
    :return: float describing absolute total magnetisation, M
    """

    return np.sum(spin)

def mag_stag(spin):
    sng=np.ones((lx,lx))

    for x in range(lx):
        for y in range(lx):
            sng[x,y]=(-1)**(x+y)

    return np.sum(np.multiply(spin, sng))
    
