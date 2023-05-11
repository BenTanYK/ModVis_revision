"""
2020 Modelling and Visualisation Past Paper parts e) and f)
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

nstep=300
lx=int(sys.argv[1])
ly=lx
p=float(sys.argv[2]) #probability S -> I

from eqns import NN_update, I_test, update_lattice, number_I

def init_condition(lx):
    """
    Initialises lattice of inactive sites with one active site placed
    at some random position
    """
    lattice=np.zeros((lx,lx))
    x=random.randint(0,lx-1)
    y=random.randint(0,lx-1)
    lattice[x,y]=-1

    return lattice

def main():

    df=pd.DataFrame()
    df['Number of sweeps']=10*np.arange(nstep/10)

    for m in range(100):
        lattice=init_condition(lx)
        number_infected=[]

        for n in range(nstep): #nstep sweeps

            for k in range(lx*ly): #2500 attempted flips per sweep
                
                lattice=update_lattice(lattice, p)

            if n%10==0:

                if number_I(lattice)==0:
                    number_infected.append(0)
                else:
                    number_infected.append(1)


        column_name='Run' + str(m)
        df[str(column_name)]=np.array(number_infected)

        print('Number of simulations performed: ' + str(m))

    df.to_excel('survival.xlsx')

main()    


        
