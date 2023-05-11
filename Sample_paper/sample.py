"""
ModVis sample paper
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

lx=int(sys.argv[1])
D=1
kappa=float(sys.argv[2])
sigma=10
v_o=float(sys.argv[3])

dx=1
dt=0.1
nstep=8000

def init_phi(lx, phi_o):
    """
    Returns randomly generated lx*ly array based on the specified initial condition, phi_o=0.5.
    """    
    #initialise empty NxN array 
    phi=np.zeros((lx,lx),dtype=float) 

    #initialise values of phi randomly
    for i in range(lx):
        for j in range(lx):
            r=random.random()*0.2-0.1
            phi[i,j]=phi_o+r    

    return phi

def grad(lattice, dx): 
    """
    Uses dell operator on NxN lattice to estimate the gradient field 
    
    :param lattice: NxN lattice
    :return: gradient vector field of NxN lattice 
    """    

    grad_x=(1/(2*dx))*(np.roll(lattice,1,axis=0)-np.roll(lattice,-1,axis=0))
    grad_y=(1/(2*dx))*(np.roll(lattice,1,axis=1)-np.roll(lattice,-1,axis=1))

    return np.array([grad_x, grad_y])

def grad_sq(lattice, dx):
    """
    Estimates laplacian of a NxN lattice 
    
    :param lattice: NxN lattice
    :param dx: spatial discretisation step 
    :return: lattice of laplacian values
    """  
    laplacian = (1/dx**2)*(np.roll(lattice,1,axis=0) + np.roll(lattice,-1,axis=0) + np.roll(lattice,1,axis=1) + np.roll(lattice,-1,axis=1) - 4*lattice)
    return laplacian

def init_rho(lx, sigma):
    """
    Calcuates distance to lattice centre at all points and returns rho lattice and 
    lattice of radial distances.
    """
    R=np.zeros((lx,lx))

    for i in range(lx):
        for j in range(lx):

            c=int(lx/2)

            x=c-i
            y=c-j

            r_sq=x**2+y**2

            R[i,j]=r_sq

    rho=np.exp(R/sigma**2)    

    return rho

def v_x(lx, v_o):
    """
    Returns v_x(y) field for an NxN lattice
    """
    N=lx*lx
    v_x=np.zeros((lx,lx))

    for x in range(lx):
        for y in range(lx):
            v_x[x,y]=-v_o*np.sin(2*math.pi*y/N)

    return v_x

def update_phi(phi, rho, kappa, dx, dt):
    """
    Returns updated phi lattice based on the first order Euler time integration
    """
    phi_new=dt*(D*grad_sq(phi, dx)+rho-kappa*phi)+phi

    return phi_new

def update_phi_vel(phi, rho, kappa, v_o, dx, dt):
    """
    Returns updated phi lattice based on the first order Euler time integration.
    A y-dependent v_x field is also introduced.
    """
    phi_new=dt*(D*grad_sq(phi, dx)+rho-kappa*phi-np.multiply(v_x(lx, v_o), grad(phi, dx)[0]))+phi

    return phi_new

def main():

    rho=init_rho(lx, sigma)
    phi=init_phi(lx, 0.5)
    
    averages=[]
    times=[]

    for n in range(nstep):

        phi_new=update_phi_vel(np.copy(phi), rho, kappa, v_o, dx, dt)

        phi=np.copy(phi_new)

        if n%7950==0 and n!=0:
            plt.cla()
            im=plt.imshow(phi, cmap='gnuplot', animated=True)
            plt.colorbar()
            plt.draw()
            plt.pause(5.0)
            #plt.savefig('v_o=0.5.png')
            plt.clf()

            print(n)

            averages.append(np.average(phi))
            times.append(n)

    # plt.title('Average value of phi')
    # plt.xlabel('Number of updates, n')
    # plt.ylabel('Average value of phi')
    # plt.plot(np.array(times), np.array(np.array(averages)))
    # plt.savefig('Average_phi')

    # df=pd.DataFrame(phi)
    # df.to_excel('phi.xlsx')

   
    
    
    

main()