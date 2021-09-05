# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 15:00:52 2020

@author: kongfh
"""
import numpy as np
import matplotlib.pyplot as plt
import math

def position(x,y,a0):
    #2D lattice
    position=np.array([x*1.0*a0,y*1.0*a0])
    return position

def criteria(E1,E2,temperature):
    energydiff=0.0
    move,probability=int(0),0.0
    #E1:energy of the current state. E2:energy of the prepared next state
    if temperature>1e-10:
        beta=1.0/((1.38e-23)*temperature)
        randnum=np.random.rand()
        if E2-E1>0:
            probability=np.exp(-beta*(E2-E1))
        else:
            probability=1.0
        if randnum<=probability:
            move=int(1)
            energydiff=E2-E1
    else:
        if E2<E1:
            move=int(1)
            energydiff=E2-E1
    #print(energydiff)
    return move,energydiff

def distance(position1,position2,size):
    #size: real value
    vector=np.abs(position1-position2)
    for i in range(len(vector)):
        while vector[i]>size*0.5:
            vector[i]=np.abs(vector[i]-size)
    distance=np.sqrt(vector[0]**2+vector[1]**2)
    return distance

def LJpotential(position1,position2,epsilon,sigma,size,cutoff):
    #position1 and position2 are actual positions
    #cutoff:scaled in terms of sigma
    vector=np.abs(position1-position2)
    for i in range(len(vector)):
        while vector[i]>size*0.5:
            vector[i]=np.abs(vector[i]-size)
    distance=np.sqrt(vector[0]**2+vector[1]**2)
    Urc=4*epsilon*(np.power(1.0/cutoff,12)-np.power(1.0/cutoff,6))
    Ur=4*epsilon*(np.power(sigma/distance,12)-np.power(sigma/distance,6))
    if distance<sigma*cutoff:
        result=Ur-Urc
    else:
        result=0.0
    return result

def simulationbox(N,percent,a0,epsilonaa,epsilonbb,epsilonab,sigma,cutoff):
    #2D lattice, 2D box, N*N particles
    #A and B, 2 components. percent: the percentage of A. A--(-1) B--1
    #a0:distance between neighboring particles, scaled in terms of sigma
    totalene=0.0
    info=np.zeros((N,N,3),dtype='float64')
    neighbor=np.zeros((N,N,N*N,2),dtype='int32')
    countneighbor=np.zeros((N,N),dtype='int32')
    dist,size=a0*sigma,a0*sigma*N
    numA=round(N*N*percent)
    numB=N*N-numA
    for i in range(N):
        for j in range(N):
            pos=position(i,j,dist)
            info[i][j][0],info[i][j][1]=pos[0],pos[1]
            #print("{} {}".format(numA,numB))
            if numA>0:
                info[i][j][2]=-1
                numA=numA-1
            else:
                info[i][j][2]=1
                numB=numB-1
    for i in range(N):
        for j in range(N):
            for k in range(N):
                for l in range(N):
                    if i==k and j==l:
                        continue
                    far=distance(position(i,j,dist),position(k,l,dist),size)
                    if far-cutoff*sigma<0:
                        neighbor[i][j][countneighbor[i][j]]=np.array([k,l])
                        countneighbor[i][j]=countneighbor[i][j]+1
                        if info[i][j][2]<0 and info[k][l][2]<0:
                            totalene=totalene+LJpotential(position(i,j,dist),position(k,l,dist), \
                                 epsilonaa,sigma,size,cutoff)
                        if info[i][j][2]>0 and info[k][l][2]>0:
                            totalene=totalene+LJpotential(position(i,j,dist),position(k,l,dist), \
                                 epsilonbb,sigma,size,cutoff)
                        if info[i][j][2]*info[k][l][2]<0.0:
                            totalene=totalene+LJpotential(position(i,j,dist),position(k,l,dist), \
                                 epsilonab,sigma,size,cutoff)
    return info,neighbor,countneighbor,totalene

def onestep(info,neighbor,countneighbor,N,percent,a0, \
            epsilonaa,epsilonbb,epsilonab,sigma,cutoff,temperature,totalene):
    #2D lattice, 2D box, N*N particles
    #A and B, 2 components. percent: the percentage of A. A--(-1) B--1
    #a0:distance between neighboring particles, scaled in terms of sigma
    info2=np.zeros((N,N,3),dtype='float64')
    #3 entries: realpositionx,realpositiony, A or B
    for i in range(N):
        info2[i]=info[i]
    dist,size=a0*sigma,a0*sigma*N
    choose,select=np.random.randint(0,N,(2)),np.random.randint(0,N,(2))
    while (choose[0]==select[0] and choose[1]==select[1])or( \
             info[choose[0]][choose[1]][2]==info[select[0]][select[1]][2]):
        choose,select=np.random.randint(0,N,(2)),np.random.randint(0,N,(2))
    info2[choose[0]][choose[1]][2]=info2[choose[0]][choose[1]][2]*(-1.0)
    info2[select[0]][select[1]][2]=info2[select[0]][select[1]][2]*(-1.0)
    energy1,energy2=0.0,0.0
    #calculate energy1
    for i in range(countneighbor[choose[0]][choose[1]]):
        pos=position(neighbor[choose[0]][choose[1]][i][0], \
                    neighbor[choose[0]][choose[1]][i][1],dist)
        value1=info[choose[0]][choose[1]][2]
        value2=info[neighbor[choose[0]][choose[1]][i][0]][neighbor[choose[0]][choose[1]][i][1]][2]
        if value1*value2<0:
            energy1=energy1+LJpotential(position(choose[0],choose[1],dist), \
                                        pos,epsilonab,sigma,size,cutoff)
        if value1<0 and value2<0:
            energy1=energy1+LJpotential(position(choose[0],choose[1],dist), \
                                        pos,epsilonaa,sigma,size,cutoff)
        if value1>0 and value2>0:
            energy1=energy1+LJpotential(position(choose[0],choose[1],dist), \
                                        pos,epsilonbb,sigma,size,cutoff)
    for i in range(countneighbor[select[0]][select[1]]):
        pos=position(neighbor[select[0]][select[1]][i][0], \
                     neighbor[select[0]][select[1]][i][1],dist)
        value1=info[select[0]][select[1]][2]
        value2=info[neighbor[select[0]][select[1]][i][0]][neighbor[select[0]][select[1]][i][1]][2]
        if value1*value2<0:
            energy1=energy1+LJpotential(position(select[0],select[1],dist), \
                                    pos,epsilonab,sigma,size,cutoff)
        if value1<0 and value2<0:
            energy1=energy1+LJpotential(position(select[0],select[1],dist), \
                                    pos,epsilonaa,sigma,size,cutoff)
        if value1>0 and value2>0:
            energy1=energy1+LJpotential(position(select[0],select[1],dist), \
                                    pos,epsilonbb,sigma,size,cutoff)
    #calculate energy2
    for i in range(countneighbor[choose[0]][choose[1]]):
        pos=position(neighbor[choose[0]][choose[1]][i][0], \
                     neighbor[choose[0]][choose[1]][i][1],dist)
        value1=info2[choose[0]][choose[1]][2]
        value2=info2[neighbor[choose[0]][choose[1]][i][0]][neighbor[choose[0]][choose[1]][i][1]][2]
        if value1*value2<0:
            energy2=energy2+LJpotential(position(choose[0],choose[1],dist), \
                                        pos,epsilonab,sigma,size,cutoff)
        if value1<0 and value2<0:
            energy2=energy2+LJpotential(position(choose[0],choose[1],dist), \
                                        pos,epsilonaa,sigma,size,cutoff)
        if value1>0 and value2>0:
            energy2=energy2+LJpotential(position(choose[0],choose[1],dist), \
                                        pos,epsilonbb,sigma,size,cutoff)
    for i in range(countneighbor[select[0]][select[1]]):
        pos=position(neighbor[select[0]][select[1]][i][0], \
                     neighbor[select[0]][select[1]][i][1],dist)
        value1=info2[select[0]][select[1]][2]
        value2=info2[neighbor[select[0]][select[1]][i][0]][neighbor[select[0]][select[1]][i][1]][2]
        if value1*value2<0:
            energy2=energy2+LJpotential(position(select[0],select[1],dist), \
                                        pos,epsilonab,sigma,size,cutoff)
        if value1<0 and value2<0:
            energy2=energy2+LJpotential(position(select[0],select[1],dist), \
                                        pos,epsilonaa,sigma,size,cutoff)
        if value1>0 and value2>0:
            energy2=energy2+LJpotential(position(select[0],select[1],dist), \
                                        pos,epsilonbb,sigma,size,cutoff)
    res=criteria(energy1,energy2,temperature)
    move,energydiff=res[0],res[1]
    if move==1:
        info=info2
    totalene=totalene+energydiff
    return info,neighbor,countneighbor,totalene

def MonteCarlo(N,percent,a0,epsilonaa,epsilonbb,epsilonab, \
               sigma,cutoff,temperature,rounds):
    initial=simulationbox(N,percent,a0,epsilonaa,epsilonbb,epsilonab,sigma,cutoff)
    info,neighbor,countneighbor,totalene=initial[0],initial[1],initial[2],initial[3]
    #display(info,epsilonaa,epsilonbb,epsilonab,sigma)
    energy=np.zeros(rounds,dtype='float64')
    for q in range(rounds):
        if q==0:
            display(info,epsilonaa,epsilonbb,epsilonab,sigma,q+1,percent,temperature)
        energy[q]=totalene
        print("{}  ".format(q),end='')
        state=onestep(info,neighbor,countneighbor,N,percent,a0, \
            epsilonaa,epsilonbb,epsilonab,sigma,cutoff,temperature,totalene)
        info,neighbor,countneighbor,totalene=state[0],state[1],state[2],state[3]
        if (q+1)%1000==0:
            display(info,epsilonaa,epsilonbb,epsilonab,sigma,q+1,percent,temperature)
    #display(info,epsilonaa,epsilonbb,epsilonab,sigma)
    displayene(energy,percent,temperature)
    return info,energy

def display(info,epsilonaa,epsilonbb,epsilonab,sigma,rounds,percent,temperature):
    N=len(info)
    listax,listay,listbx,listby=list(),list(),list(),list()
    #info2=np.zeros((N,N,3),dtype='float64')
    for i in range(N):
        for j in range(N):
            if info[i][j][2]<0.0:
                listax.append(info[i][j][0]/sigma)
                listay.append(info[i][j][1]/sigma)
            else:
                listbx.append(info[i][j][0]/sigma)
                listby.append(info[i][j][1]/sigma)
    plt.figure(figsize=(10,8))
    plt.gca().set_aspect('equal')
    sizea=19/(N/20)
    plt.plot(listax,listay,'o',color="red",markersize=sizea)
    plt.plot(listbx,listby,'o',color="#FFB90F",markersize=sizea)
    plt.xlabel("aa:{:.3f}     bb:{:.3f}     ab:{:.3f}\nA: red      B: yellow".format( \
                  epsilonaa/1e-21,epsilonbb/1e-21,epsilonab/1e-21),fontsize=15)
    plt.ylabel("(sigma)",fontsize=15)
    plt.xticks(size=15)
    plt.yticks(size=15)
    plt.title("Monte Carlo —— A-B alloy (step {})\nT={:.1f}K   A:{:.1f}%   B:{:.1f}%".format( \
                            rounds,temperature,percent*100,100-percent*100),fontsize=20)
    plt.savefig("F:/MCHW/HWc30{}.jpg".format(rounds),dpi=600)
    plt.show()
    
def displayene(energy,percent,temperature):
    X=np.arange(len(energy))
    plt.figure(figsize=(10,8))
    plt.plot(X,energy,'b-')
    plt.xlabel("the number of steps",fontsize=15)
    plt.ylabel("total energy (Joule)",fontsize=15)
    plt.xticks(size=12)
    plt.yticks(size=12)
    plt.title("Total Energy\nT={:.1f}K   A:{:.1f}%   B:{:.1f}%".format( \
                           temperature,percent*100,100-percent*100),fontsize=20)
    #plt.ylim(-2.16e-17,-2e-17)
    plt.savefig("F:/MCHW/totaleneHWc30.jpg",dpi=600)
    plt.show()

#def MonteCarlo(N,percent,a0,epsilonaa,epsilonbb,epsilonab, \
#               sigma,cutoff,temperature,rounds):
MonteCarlo(40,0.3,1.12,1.67e-21,2.227e-21,1.928e-21,3.4e-10,3.0,20,10000)
#display(result,3.4e-10)




#result=simulationbox(4,0.5,1.5,3.4e-10,2.0)
#print(result)
#info,neighbor,countneighbor=result[0],result[1],result[2]
#onestep(info,neighbor,countneighbor,4,0.5,1.5, \
#              1.67e-21,4*1.67e-21/3,1.93e-21,3.4e-10,2.0,100)                    
            

'''
X=np.linspace(3.4e-10*0.95,3.4e-10*2.5,200)
Y=np.zeros(200)
for i in range(200):
    Y[i]=LJpotential(np.array([X[i],0]),np.array([0.0,0.0]),1.67e-21,3.4e-10,2.0)
plt.plot(X,Y,'r-')
plt.show()
'''

















