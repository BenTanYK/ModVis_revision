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

lx=int(sys.argv[1])
R=float(sys.argv[2])

dx=2
dt=1
nstep=10000

F=0.035
k=0.06
D1=0.2
D2=0.1

def init_U(lx):
    """
    Initialises lx*lx array containing initial values of U, subject to
    the given initial conditions.
    """

    lattice=np.zeros((lx,lx))

    for i in range(lx):
        for j in range(lx):
            noise=0.02*random.random()-0.01

            c=int(lx/2)

            x=c-i
            y=c-j

            r=x**2+y**2

            if r>R:
                lattice[i,j]=1.0+noise

            else:
                lattice[i,j]=0.5+noise
    
    return lattice

def init_V(lx):
    """
    Initialises lx*lx array containing initial values of V, subject to
    the given initial conditions.
    """

    lattice=np.zeros((lx,lx))

    for i in range(lx):
        for j in range(lx):
            noise=0.02*random.random()-0.01

            c=int(lx/2)

            x=c-i
            y=c-j

            r=x**2+y**2

            if r>R:
                lattice[i,j]=0.01+noise

            else:
                lattice[i,j]=0.25+noise
    
    return lattice

def U_new(U, V, D1, F):
    """
    Returns updated U lattice using basic Euler time-integration scheme
    """
    U_new=np.copy(U)+dt*(D1*laplacian(U, dx)-np.multiply(U, np.square(V))+F*(np.ones((lx, lx))-U))

    return U_new

def V_new(U, V, D2, F, k):
    """
    Returns updated V lattice using basic Euler time-integration scheme
    """
    V_new=np.copy(V)+dt*(D2*laplacian(V, dx)+np.multiply(U, np.square(V))-(F+k)*V)

    return V_new

def laplacian(lattice, dx):
    """
    Returns lattice of grad squared values
    """
    grad_sq = (1/dx**2)*np.roll(lattice,1,axis=1) + np.roll(lattice,-1,axis=1) + np.roll(lattice,1,axis=0) + np.roll(lattice,-1,axis=0) - 4*lattice
    return grad_sq

def var_U(U):
    """
    Output spatial variance for a given U lattice. 
    """

    lx=int(len(U[0])) #lattice size length
    N=lx*lx #total number of lattice sites

    U_av=np.sum(U)/N

    U_sq=np.square(U)
    U_sq_av=np.sum(U_sq)/N

    var=U_sq_av-U_av**2

    return var

def main_b():

    U=init_U(lx)
    V=init_V(lx)

    for n in range(nstep):

        U_updated=U_new(U, V, D1, F)
        V_updated=V_new(U, V, D2, F, k)

        U=np.copy(U_updated)
        V=np.copy(V_updated)

        if n%10==0:

            plt.cla()
            im=plt.imshow(U_updated, cmap='viridis', animated=True)
            plt.colorbar()
            plt.draw()
            plt.pause(0.0001)
            #plt.savefig('U.png')
            plt.clf()

def main_c():

    F_list=[]
    var_list=[]

    for m in range(8):
        
        F=0.020+0.005*m
        F_list.append(F)

        U=init_U(lx)
        V=init_V(lx)

        for n in range(nstep):

            U_updated=U_new(U, V, D1, F)
            V_updated=V_new(U, V, D2, F, k)

            U=np.copy(U_updated)
            V=np.copy(V_updated)

            if n%100==0:
                print(n)

                plt.cla()
                im=plt.imshow(U_updated, cmap='viridis', animated=True)
                plt.colorbar()
                plt.draw()
                plt.pause(0.0001)
                #plt.savefig('U.png')
                plt.clf()

        var=var_U(U)
        var_list.append(var)

        df=pd.DataFrame()
        df['F']=np.array(F_list)
        df['Variance']=np.array(var_list)
        df.to_excel('c.xlsx')



main_b()






    
