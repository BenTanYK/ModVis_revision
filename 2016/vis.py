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

J=-1.0 #set constant J to -1
k_B=1.0 #Boltzmann constant
nstep=1000

lx=50
T=1
h=float(sys.argv[1])

from eqns import nearest_neighbours, delta_E, Glauber, init_random

def main():

    lattice=init_random(lx)

    for n in range(nstep): #nstep sweeps

        for m in range(lx*lx): #2500 attempted flips per sweep
            
            lattice=Glauber(lattice, T, h)

        if n%100==0:
            plt.cla()
            im=plt.imshow(lattice, cmap='gnuplot', animated=True, vmin=-1, vmax=1 )
            plt.draw()
            plt.pause(0.0001)
            
            print('Number of sweeps: ' + str(n))
main()




