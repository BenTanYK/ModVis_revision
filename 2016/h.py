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
nstep=200

lx=50
T=1

from eqns import nearest_neighbours, delta_E, Glauber, init_random, magnetisation, mag_stag

def main():

    h_list=[]
    M_av_list=[]
    var_M_list=[]
    M_stag_av_list=[]
    var_M_s_list=[]

    for k in range(21):
        print('Number of runs:' + str(k))

        h=0+0.5*k

        lattice=init_random(lx)
        mag_list=[]
        mag_stag_list=[]

        for n in range(nstep): #nstep sweeps

            for m in range(lx*lx): #2500 attempted flips per sweep
                
                lattice=Glauber(lattice, T, h)

            if n%10==0 and n>100: #calculate M and M_s every 10 sweeps, allow for equilibration time of 100 sweeps
                M=magnetisation(lattice)
                M_s=mag_stag(lattice)

                mag_list.append(M)
                mag_stag_list.append(M_s)

        M_av=np.average(np.array(mag_list))
        M_stag_av=np.average(np.array(mag_stag_list))

        M_sq_av=np.average(np.square(np.array(mag_list)))
        M_s_sq_av=np.average(np.square(np.array(mag_stag_list)))

        var_M=M_sq_av-M_av**2
        var_M_s=M_s_sq_av-M_stag_av**2

        h_list.append(h)
        M_av_list.append(M_av)
        var_M_list.append(var_M)
        M_stag_av_list.append(M_stag_av)
        var_M_s_list.append(var_M_s)

    #export data to pd DataFrame
    df=pd.DataFrame()
    df['h']=np.array(h_list)
    df['M_av']=np.array(M_av_list)
    df['var M_av']=np.array(var_M_list)
    df['M_s_av']=np.array(M_stag_av_list)
    df['var M_s_av']=np.array(var_M_s_list)

    df.to_excel('magnetisation.xlsx')
   
main()

