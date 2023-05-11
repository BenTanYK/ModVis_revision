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
nstep=2000

lx=50
T=1

P=25
t=10000
h_o=10

from eqns import nearest_neighbours, init_random, magnetisation, mag_stag

def h_field(P,t,n,h_o):
    """
    Returns spatially and temporally dependent h-field
    """
    h=h_o*np.ones((lx,lx))

    for x in range(lx):
        for y in range(lx):
            h[x,y]=np.cos(2*math.pi*x/P)*np.cos(2*math.pi*y/P)*np.cos(2*math.pi*n/t)

    return h

def delta_E(spin, i, j, J, h):
    
    """
    Computes change in energy which arises from flipping the spin at position
    [i, j] in the NxN array of spins.

    """

    #generate array of nearest neighbour spins
    NN_spins=nearest_neighbours(spin, i, j)
    
    #change in energy due to flipping spin [i,j]
    E_change=2*J*spin[i,j]*np.sum(NN_spins)+2*np.multiply(h[i,j], spin[i,j])
    
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

def main():

    lattice=init_random(lx)

    for n in range(nstep): #nstep sweeps

        h=h_field(P, t, n, h_o)

        for m in range(lx*lx): #2500 attempted flips per sweep
            
            lattice=Glauber(lattice, T, h)

        if n%10==0:
            plt.cla()
            im=plt.imshow(lattice, cmap='gnuplot', animated=True, vmin=-1, vmax=1 )
            plt.draw()
            plt.pause(0.0001)
            
            print('Number of sweeps: ' + str(n))

main()
