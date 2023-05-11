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

from eqns import NN_update, I_test, update_lattice, number_I

def var(lx, p):

    """
    Performs nstep sweeps of the simulation and calculates the number of infected cells
    at every sweep. The function returns the variance of the number of infected sites.

    :param lx: size of initially random array of sites
    :param p: probability for a nearest neighbour S -> I
    :return: average number of infected sites, variance in fraction of infected sites
    """
    ly=lx

    lattice=init_random(lx)
     
    n_infected=[]

    for n in range(nstep): #nstep sweeps

        for m in range(lx*ly): #2500 attempted flips per sweep

            lattice=update_lattice(lattice, p)
                      
        n_I=number_I(lattice) #calculate fraction of infected sites every sweep
        
        if n_I==0: #detect absorbing phase and break loop
            n_infected.append(0.0)
            break

        else:
            n_infected.append(n_I)

    if n_infected[-1]==0: #if absorbing state has been reached, set variance to zero
        var=0
        n_infected=np.zeros(nstep)

    else:
        n_I_equil=np.array(n_infected)[100:] #remove first 100 results due to equilibration time
        av=np.average(n_I_equil)
        av_sq=np.average(np.square(n_I_equil))

        var=(av_sq-av**2)/(lx*ly)

    return n_infected, var

def boostrap_var(n_infected, lx, r):

    """
    Calculates the error in the calculated variance for a given value of p using the
    bootstrap method.
    
    :param n_infected: array containing number of infected sites over time
    :param N: total number of sites in the simulation
    :param r: number of resamplings
    :return: float giving the error in the calculated scaled heat capacity, based
    on the bootstrap resampling method
    """
    var_list=[]
    ly=lx    

    for x in range(r):

        n_resampled=[]

        for n in range(len(n_infected)):

            rand_int=random.randint(0, len(n_infected)-1)

            n_resampled.append(n_infected[rand_int])

        av=np.average(np.array(n_resampled))

        av_sq=np.average(np.square(np.array(n_resampled)))

        var=(av_sq-av**2)/(lx*ly)

        var_list.append(var)
    
    var_av=np.average(np.array(var_list))
    var_sq_av=np.average(np.square(np.array(var_list)))

    return np.sqrt(var_sq_av-var_av**2)

#Define two different main functions for the various exam questions

def main_b():
    
    lattice=init_random(lx)
     
    # fig = plt.figure()
    # im=plt.imshow(lattice, cmap=cmap, animated=True)
    # plt.colorbar()

    n_sweeps=[]
    infected_frac=[]

    for n in range(nstep): #nstep sweeps

        for m in range(lx*ly): #2500 attempted flips per sweep
            
            lattice=update_lattice(lattice, p)

    #plot measurements every 10 sweeps
        if n%10==0:
    #       show animation
            # plt.cla()
            # im=plt.imshow(lattice, cmap=cmap, animated=True, vmin=-1, vmax=0 )
            # plt.draw()
            # plt.pause(0.0001)

            print(n)
            n_sweeps.append(n)
            infected_frac.append(number_I(lattice)/(lx**2))

    plt.title('Infected fraction over time')
    plt.ylabel('Fraction of infected sites')
    plt.xlabel('Number of sweeps')
    plt.scatter(np.array(n_sweeps), np.array(infected_frac))
    plt.show()

def main_c():

    variances=[]
    av_frac=[]
    error_var=[]

    for n in range(31):

        p=0.55+0.005*n

        n_infected, v=var(lx, p) #compute variance for given p
        error_v=boostrap_var(n_infected, lx, r)
        av=np.average(n_infected)

        variances.append(v)
        error_var.append(error_v)
        av_frac.append(av)

        print('run ' + str(n) + "/31")

    #initialise dataframe and export to excel
    df=pd.DataFrame()
    df['p']=np.linspace(0.55, 0.70, 31)
    df['variance']=np.array(variances)
    df['variance error']=np.array(error_var)
    df['Average infected fraction']=np.array(av_frac)
    df.to_excel("2020_variance.xlsx")
    print(df)

main_c()
