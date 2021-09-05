import numpy as np
import matplotlib.pyplot as plt

def FCC(x1,x2,x3,a0):
    #generate the position of particles in a FCC lattice
    vector1=np.array([0.5,0.5,0],dtype='float64')
    vector2=np.array([0,0.5,0.5],dtype='float64')
    vector3=np.array([0.5,0,0.5],dtype='float64')
    position=(x1*vector1+x2*vector2+x3*vector3)*a0
    return position

def distance(position1,position2,size):
    #calculate the distance between two particles
    #size: size of the simulation box
    position1=position1.astype('float64')
    position2=position2.astype('float64')
    dist=np.abs(position1-position2)
    for i in range(len(dist)):
        if dist[i]>size/2.0:
            dist[i]=dist[i]-size/2.0
    distance=np.sqrt(dist[0]**2+dist[1]**2+dist[2]**2)
    return distance

def LJpotential(epsilon,sigma,distance,cutoff,mode):
    #cutoff: cutoff radius, how many times it is of sigma.
    #mode: 
    #mode1 is the potential demonstrated in the requirement of the assignment
    #mode2 is the continuous version of the potential
    #mode3 is the continuous and smooth version
    x=sigma/distance
    Xc=sigma/cutoff
    Urc=4*epsilon*(np.power(Xc,12)-np.power(Xc,6))
    Vx=4*epsilon*(12*np.power(cutoff,11)-6*np.power(cutoff,5))* \
            (-sigma/(np.power(cutoff,2)))
    if mode==1:
        if distance<=cutoff:
            Ur=4*epsilon*(np.power(x,12)-np.power(x,6))
        else:
            Ur=0
    if mode==2:
        if distance<=cutoff:
            Ur=4*epsilon*(np.power(x,12)-np.power(x,6))-Urc
        else:
            Ur=0
    if mode==3:
        if distance<=cutoff:
            Ur=Ur=4*epsilon*(np.power(x,12)-np.power(x,6))-Urc-Vx
        else:
            Ur=0
    return Ur
    
def SumPotential(a0,size,epsilon,sigma,cutoff,mode):
    #a0:lattice parameter  size:size of the simulation box
    #mode: the same variable in LJpotential function
    n=1
    sumpotential=0.0
    while a0*n<2*cutoff:
        n=n+1
    #one particle in the center
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if i==0 and j==0 and k==0:
                    continue
                position=FCC(i,j,k,a0)
                dist=distance(position,np.array([0.0,0.0,0.0]),size)
                sumpotential=sumpotential+LJpotential(epsilon,sigma, \
                                                      dist,cutoff,mode)
                #print(sumpotential)
    return sumpotential

def find(start,stop,step,size,epsilon,sigma,cutoff,mode):
    X=np.arange(start*sigma,stop*sigma,step*sigma,dtype='float64')
    X1=X/sigma
    Y=np.zeros(len(X),dtype='float64')
    minimum,minix,minidex = 0.0, 0.0, 0.0
    choice=np.array([0.0,0.0],dtype='float64')
    for i in range(len(X)):
        Y[i]=SumPotential(X[i],size*sigma,epsilon,sigma,cutoff*sigma,mode)
    Y1=Y/epsilon
    for i in range(len(Y)):
        if Y[i]<minimum:
            minimum=Y[i]
            minix=X[i]
            minidex=X1[i]
    #unit: angstrom
    choice[0],choice[1]=minix/(1e-10),minidex
    plt.plot(X1,Y1,'b-')
    plt.show()
    return choice

start,stop,step,size,epsilon,sigma,cutoff,mode= \
    1.55, 1.6, 0.00001, 100, 1.67*1e-14, 3.4*1e-10, 3, 3
print(find(start,stop,step,size,epsilon,sigma,cutoff,mode))

'''
eps,sig,cut,step=1.67e-14,3.4e-8,3*3.4e-8,1e-10
X=np.arange(1.2*sig,3*sig,step)
X1=np.arange(1.2,3,step/sig)
Y=np.zeros(len(X))
for i in range(len(X)):
    Y[i]=SumPotential(X[i],20*sig,eps,sig,cut,3)
#print(np.array([X,Y]))
plt.plot(X1,Y,'b.')
plt.show() 
'''               
'''     
a=np.array([3,0,1])
b=np.array([0,1,5])
size=3
print(distance(a,b,size))'''

'''
X=np.arange(0.95*3.4*1e-4,3.4*3*1e-4,1e-5)
Y=np.zeros(len(X))
for i in range(len(X)):
    Y[i]=LJpotential(1.67e-8,3.4*1e-4,X[i],3.4*2*1e-4,3)

plt.plot(X,Y,'r-')
plt.show()
'''






