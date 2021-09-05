import numpy as np
import matplotlib.pyplot as plt
import math

def info():
    #collect the mass, elastic constants
    N=eval(input("Please input the number of mass points:\n"))
    INFO=np.zeros(3*N-2,dtype='float64')
    lista=list(map(eval, input("Please input the mass of each \
mass points:\n(from left to right, separated by comma)\n").split(",")))
    listb=list(map(eval, input("Please input the elastic constant of each spring:\n \
(from left to right, separated by comma)\n").split(",")))
    listc=list(map(eval,(input("Please input the length of each spring when \
potential energy is zero:\n(from left to right, separated by comma)\n").split(","))))
    INFO=np.array(lista+listb+listc).astype('float64')
    #print("INFO:{}".format(INFO))
    return INFO

def initial(INFO):
    #initialization, position, velocity
    N=(len(INFO)+2)//3
    STATE=np.zeros((2,N),dtype="float64")
    lista=list(map(eval, input("Please input the position of each mass point:\n \
(from left to right, separated by comma)\n").split(",")))
    listb=list(map(eval, input("Please input the velocity of each mass point:\n \
(from left to right, separated by comma)\n").split(",")))
    STATE=np.array([lista,listb]).astype("float64")
    return STATE
#print(initial(info()))

def energy(INFO,STATE):
    N=(len(INFO)+2)//3
    ENERGY=np.zeros((2,N),dtype='float64')
    ENERGY[1][N-1]=0.0
    #first N entries:kinetic energy. second N entries: potential energy. The rest: sum of energy
    listinfo=INFO.tolist()
    mass=np.array((listinfo)[0:N]).astype('float64')
    K=np.array((listinfo)[N:2*N-1]).astype('float64')
    x0=(listinfo)[2*N-1:3*N-2]
    #ENERGY[0]=0.5*mass*STATE[1]**2
    for i in range(N):
        #kinetic energy
        ENERGY[0][i]=0.5*mass[i]*STATE[1][i]**2
        if i==N-1:
            continue
        #potential energy
        #print((STATE[0][i+1]-STATE[0][i])-x0[i])
        ENERGY[1][i]=0.5*K[i]*((STATE[0][i+1]-STATE[0][i])-x0[i])**2
    TOTAL=np.sum(ENERGY)
    return ENERGY,TOTAL

def acc(INFO,STATE):
    N=(len(INFO)+2)//3
    FORCE=np.zeros(N-1,dtype='float64')
    ACC=np.zeros(N,dtype='float64')
    listinfo=INFO.tolist()
    mass=np.array((listinfo)[0:N]).astype('float64')
    K=np.array((listinfo)[N:2*N-1]).astype('float64')
    x0=(listinfo)[2*N-1:3*N-2]
    #print("K: {}".format(K))
    for i in range(N-1):
        dev=(STATE[0][i+1]-STATE[0][i])-x0[i]
        FORCE[i]=K[i]*dev
        #when compressed, force is negative. when streched, force is positive
        #only the interaction between neighboring mass points is considered
    #print("force {}".format(FORCE))
    for i in range(N):
        if i==0:
            ACC[i]=FORCE[i]/mass[i]
        elif i==N-1:
            ACC[i]=-FORCE[i-1]/mass[i]
        else:
            ACC[i]=(FORCE[i]-FORCE[i-1])/mass[i]
    return ACC

def masscenter(INFO,STATE):
    N=(len(INFO)+2)//3
    listinfo=INFO.tolist()
    mass=np.array((listinfo)[0:N]).astype('float64')
    summass=np.sum(mass)
    masscenter=np.zeros(2,dtype='float64')
    center=np.zeros((2,N),dtype='float64')
    for i in range(N):
        masscenter[0]=masscenter[0]+mass[i]*STATE[0][i]/summass
        masscenter[1]=masscenter[1]+mass[i]*STATE[1][i]/summass
    for i in range(N):
        center[0][i]=masscenter[0]
        center[1][i]=masscenter[1]
    return center
     
def verlet(INFO,STATE,timestep):
    N=(len(INFO)+2)//3
    #STATE: (2,N),STATE[0]:positions STATE[1]:velocities
    STATE2=np.zeros((2,N),dtype="float64")
    ACC=acc(INFO,STATE)
    ACC2=np.zeros(N,dtype='float64')
    STATE2[0]=STATE[0]+timestep*STATE[1]+0.5*(timestep**2)*ACC
    ACC2=acc(INFO,STATE2)
    STATE2[1]=STATE[1]+0.5*timestep*(ACC+ACC2)
    return STATE2

def euler(INFO,STATE,timestep):
    N=(len(INFO)+2)//3
    #STATE: (2,N),STATE[0]:positions STATE[1]:velocities
    STATE2=np.zeros((2,N),dtype="float64")
    ACC=acc(INFO,STATE)
    STATE2[0]=STATE[0]+timestep*STATE[1]
    STATE2[1]=STATE[1]+ACC*timestep
    return STATE2

def main(timestep,rounds):
    #Velocity Verlet
    INFO=info()
    STATE=initial(INFO)
    N=(len(INFO)+2)//3
    storagepos=np.zeros((rounds,N),dtype='float64')
    storagevel=np.zeros((rounds,N),dtype='float64')
    centerpos=np.zeros((rounds,N),dtype='float64')
    centervel=np.zeros((rounds,N),dtype='float64')
    storagekin=np.zeros((rounds,N),dtype='float64')
    storagepot=np.zeros((rounds,N-1),dtype='float64')
    storagetotalene=np.zeros(rounds,dtype='float64')
    for p in range(rounds):
        storagepos[p]=STATE[0]
        storagevel[p]=STATE[1]
        centerpos[p]=masscenter(INFO,STATE)[0]
        centervel[p]=masscenter(INFO,STATE)[1]
        ener=energy(INFO,STATE)
        storagekin[p]=ener[0][0]
        storagepot[p]=ener[0][1][0:-1]
        storagetotalene[p]=ener[1]
        STATE=verlet(INFO,STATE,timestep)
    return storagepos,storagevel,storagekin,storagepot,storagetotalene,centerpos,centervel

def main2(timestep,rounds):
    #Euler
    INFO=info()
    STATE=initial(INFO)
    N=(len(INFO)+2)//3
    storagepos=np.zeros((rounds,N),dtype='float64')
    storagevel=np.zeros((rounds,N),dtype='float64')
    centerpos=np.zeros((rounds,N),dtype='float64')
    centervel=np.zeros((rounds,N),dtype='float64')
    storagekin=np.zeros((rounds,N),dtype='float64')
    storagepot=np.zeros((rounds,N-1),dtype='float64')
    storagetotalene=np.zeros(rounds,dtype='float64')
    for p in range(rounds):
        storagepos[p]=STATE[0]
        storagevel[p]=STATE[1]
        centerpos[p]=masscenter(INFO,STATE)[0]
        centervel[p]=masscenter(INFO,STATE)[1]
        ener=energy(INFO,STATE)
        storagekin[p]=ener[0][0]
        storagepot[p]=ener[0][1][0:-1]
        storagetotalene[p]=ener[1]
        STATE=euler(INFO,STATE,timestep)
    return storagepos,storagevel,storagekin,storagepot,storagetotalene,centerpos,centervel

def display(timestep,rounds,result):
    position,velocity,kinetic,potential,totalene,centerpos,centervel= \
        result[0],result[1],result[2],result[3],result[4],result[5],result[6]
    N=len(result[0][0])
    X=np.arange(rounds)
    position=position.swapaxes(0,1)
    velocity=velocity.swapaxes(0,1)
    centerpos=centerpos.swapaxes(0,1)
    centervel=centervel.swapaxes(0,1)
    kinetic=kinetic.swapaxes(0,1)
    potential=potential.swapaxes(0,1)
    
    figpos,axpos=plt.subplots(figsize=(9,6))
    #figpos.figure(figsize=(6,4))
    for i in range(N):
        string='{}'.format(i+1)
        axpos.plot(X,position[i]-centerpos[i],label=string)
    axpos.set_xlabel("the number of steps\ntimestep={}".format(timestep),fontsize=10)
    axpos.set_ylabel("position",fontsize=10)
    axpos.set_title("Position (relative to mass center)",fontsize=12)
    plt.legend(loc="upper right")
    plt.savefig("F:/Oscillator2/caterpillar.jpg",dpi=500)
    plt.show()
    
    figvel,axvel=plt.subplots(figsize=(9,6))
    for i in range(N):
        string='{}'.format(i+1)
        axvel.plot(X,velocity[i]-centervel[i],label=string)
    axvel.set_xlabel("the number of steps\ntimestep={}".format(timestep),fontsize=10)
    axvel.set_ylabel("velocity",fontsize=10)
    axvel.set_title("Velocity (relative to mass center)",fontsize=12)
    plt.legend(loc="upper right")
    #plt.savefig("F:/Oscillator2/velocitye.jpg",dpi=500)
    plt.show()
    
    figkin,axkin=plt.subplots(figsize=(9,6))
    for i in range(N):
        string='{}'.format(i+1)
        axkin.plot(X,kinetic[i],label=string)
    axkin.set_xlabel("the number of steps\ntimestep={}".format(timestep),fontsize=10)
    axkin.set_ylabel("kinetic energy",fontsize=10)
    axkin.set_title("Kinetic Energy",fontsize=12)
    plt.legend(loc="upper right")
    #plt.savefig("F:/Oscillator2/kinetice.jpg",dpi=500)
    plt.show()
    
    figpot,axpot=plt.subplots(figsize=(9,6))
    for i in range(N-1):
        string='spring {}'.format(i+1)
        axpot.plot(X,potential[i],label=string)
    axpot.set_xlabel("the number of steps\ntimestep={}".format(timestep),fontsize=10)
    axpot.set_ylabel("potential energy",fontsize=10)
    axpot.set_title("Potential Energy",fontsize=12)
    plt.legend(loc="upper right")
    #plt.savefig("F:/Oscillator2/potentiale.jpg",dpi=500)
    plt.show()
    
    figtot,axtot=plt.subplots(figsize=(9,6))
    axtot.plot(X,totalene,'r-')
    axtot.set_xlabel("the number of steps\ntimestep={}".format(timestep),fontsize=10)
    axtot.set_ylabel("total energy",fontsize=10)
    axtot.set_title("Total Energy",fontsize=12)
    #plt.savefig("F:/Oscillator2/totalenee.jpg",dpi=500)
    plt.show()
    
    figlen,axlen=plt.subplots(figsize=(9,6))
    length=position[len(position)-1]-position[0]
    axlen.plot(X,length,'b-')
    axlen.set_xlabel("the number of steps\ntimestep={}".format(timestep),fontsize=10)
    axlen.set_ylabel("the total length",fontsize=10)
    axlen.set_title("Total Length",fontsize=12)
    plt.savefig("F:/Oscillator2/totallength.jpg",dpi=500)
    plt.show()
    '''
    pos=(position-centerpos).swapaxes(0,1)
    y=np.zeros(N)
    for i in range(rounds):
        plt.plot(pos[i],y,'ro')
        plt.xlim(-6,6)
        plt.show()
    '''
    
result=main(0.01,10000)
#result=main2(0.1,200)
display(0.01,10000,result)
    
#INFO=info()
#STATE=initial(INFO)
#print(masscenter(INFO,STATE))
#print(energy(INFO,STATE))
#print(verlet(INFO,STATE,0.1))












