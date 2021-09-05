import math
import numpy as np
import matplotlib.pyplot as plt

#epsilon, sigma : International System of Units

def position(x,y,z,sigma):
    #denote the position of particles
    #sigma is the important variable shown in the expression of L-J potential
    base=np.array([[1,0,0],[0,1,0],[0,0,1]],dtype='float64')
    position=(x*base[0]+y*base[1]+z*base[2])*sigma
    return position

def displacement(position1,position2):
    #calculate the displacement of particles
    #absolute value of displacement, no PBC involved
    vector=position2-position1
    displace=np.sqrt(vector[0]**2+vector[1]**2+vector[2]**2)
    return displace

def distance(position1,position2,size,sigma):
    #calculate the distance between two particles
    #size: size of the simulation box, which is scaled in terms of sigma
    #sigma: the important variable in the expression of L-J potential
    dist=np.abs(position1-position2)
    for i in range(len(dist)):
        while dist[i]>size*sigma*0.5:
            dist[i]=np.abs(dist[i]-size*sigma)
    distance=np.sqrt(dist[0]**2+dist[1]**2+dist[2]**2)
    return distance    

def force(position1,position2,size,epsilon,sigma):
    #calculate the interactive force between two particles
    #position1 and position2 are vectors
    #the interactive force is a vector, the output is a vector
    #position1:studied particle
    #position2: the particle generating the force field
    vector=position1-position2
    for i in range(len(vector)):
        if vector[i]>=0:
            while vector[i]>size*sigma/2.0:
                vector[i]=vector[i]-size*sigma
        else:
            while vector[i]<-size*sigma/2.0:
                vector[i]=vector[i]+size*sigma
    dist=np.sqrt(vector[0]**2+vector[1]**2+vector[2]**2)
    x=sigma/dist
    absforce=4*epsilon*(12*np.power(x,11)-6*np.power(x,5))*(sigma/(dist**2))
    xforce=absforce*(vector[0]/dist)
    yforce=absforce*(vector[1]/dist)
    zforce=absforce*(vector[2]/dist)
    force=np.array([xforce,yforce,zforce])
    return force

def LJpotential(epsilon,sigma,distance,cutoff):
    #distance: scaled in terms of sigma
    #cutoff: cutoff radius, how many times it is of sigma.
    x=sigma/distance
    if distance<=cutoff*sigma:
        Ur=4*epsilon*(np.power(x,12)-np.power(x,6))
    else:
        Ur=0
    return Ur

def kinetic(velocity,mass):
    #calculate the kinetic energy of a particle
    #mass: molar mass
    mreal=mass*(1e-3)/(6.02e23)
    kinetic=0.5*mreal*(velocity[0]**2+velocity[1]**2+velocity[2]**2)
    return kinetic
    
def dT(epsilon,sigma,mass,scale):
    #generate timestep, namely, delta t
    #mass: molar mass, g/mol
    #recommended scale: 0.001
    numerator=sigma
    denominator=np.sqrt(2*epsilon/(mass*1e-3/(6.02e23)))
    timestep=(numerator/denominator)*scale
    return timestep

def initial(epsilon,sigma,mass,size,N,temperature):
    #simulation setup： initialization
    #kb: Boltzmann's Constant
    #N: the number of particles
    #mass: molar mass, g/mol
    #vmean: not the average of the absolute value of velocity
    #vmean is the mean of Vx, Vy or Vz
    #the simulation box is placed in the first octant of 3D space
    kb=1.380649e-23
    dev=np.sqrt(kb*temperature/(mass*1e-3/(6.02e23)))
    startposition=np.zeros((N,3),dtype='float64')
    startvelocity=np.random.normal(0,dev,(N,3))
    #the simulation box will be diveded. part is the number of units of each axis
    part=int(math.floor(np.power(N,1.0/3))+1.0)
    unit=size*sigma/part
    count=int(N)
    i=0
    for m in range(part):
        if count==0:
            break
        for n in range(part):
            if count==0:
                break
            for q in range(part):
                startposition[i]=np.array([m*unit,n*unit,q*unit])+ \
                        0.5*unit*np.array([1.0,1.0,1.0])
                #startposition[i]=np.array([m*unit,n*unit,q*unit])+ \
                #        0.2*unit*np.random.rand(3)
                #startposition[i]=np.array([m*unit,n*unit,q*unit])+ \
                #    0.1*unit*np.array([1.0,1.0,1.0])
                #print('count={}'.format(count))
                i=i+1
                count=count-1
                #print("{} {} {}".format(m,n,q))
                if count==0:
                    break
    initial=np.array([startposition,startvelocity])
    return initial

def main(epsilon,sigma,mass,cutoff,skin,size,scale,N,temperature,rounds):
    #microstate: current position, next position, current velocity, next velocity
    #cutoff: scaled in terms of sigma
    #scale: related to time step. recommended:0.001
    #mass: molar mass  mreal: real mass of a particle
    #acc: acceleration, currentacc and nextacc'
    #rounds: the number of time steps
    jud=0 #whether the particle goes out of the box or not. 0:no 1:yes
    l,m,n=0,0,0
    kb=1.380649e-23
    mreal=mass*(1e-3)/(6.02e23)
    density=N*1.0/(np.power(size*sigma,3))
    #set up storage space
    storage=np.zeros((rounds,N,3))
    currentposition=np.zeros((N,3),dtype='float64')
    nextposition=np.zeros((N,3),dtype='float64')
    currentvelocity=np.zeros((N,3),dtype='float64')
    nextvelocity=np.zeros((N,3),dtype='float64')
    currentacc=np.zeros((N,3),dtype='float64')
    nextacc=np.zeros((N,3),dtype='float64')
    displace=np.zeros(N,dtype='float64')
    timestep=dT(epsilon,sigma,mass,scale)
    #sumkinetic: the sum of kinetic energy of all particles
    #TEMP: current temperature of the system
    sumkinetic,sumpotential,sumenergy,eneperparticle,TEMP,PRES=0.0,0.0,0.0,0.0,0.0,0.0
    #figure:temperature
    X=np.arange(rounds)
    Ytemp=np.zeros(rounds)
    Ypres=np.zeros(rounds)
    Yene=np.zeros(rounds)
    #initialization
    initialstate=initial(epsilon,sigma,mass,size,N,temperature)
    currentposition,currentvelocity=initialstate[0],initialstate[1]
    #create linked cells
    #at the same time, calculate the first currentacc
    unitlength=size
    divisor=1
    while unitlength/divisor>=cutoff:
        divisor=divisor+1
    unitlength=size/(divisor*1.0)
    print('unitlength={}'.format(unitlength))
    print('divisor={}'.format(divisor))
    linkedcell=np.full((divisor,divisor,divisor,N),-1,dtype='int32')
    #antilink: see which cell the particle is in
    antilink=np.zeros((N,3),dtype='int32')
    #the size of a subcell is obtained
    pointer=np.zeros((divisor,divisor,divisor),dtype='int32')
    #print(currentposition/sigma)
    for i in range(N):
        x,y,z=currentposition[i][0]/sigma,currentposition[i][1]/sigma, \
            currentposition[i][2]/sigma
        l,m,n=int(math.floor(x/unitlength)),int(math.floor(y/unitlength)), \
            int(math.floor(z/unitlength))
        antilink[i]=np.array([l,m,n])
        #print("{} {} {}".format(l,m,n))
        #print(currentposition[i]/sigma)
        #print(unitlength)
        linkedcell[l][m][n][(pointer[l][m][n])]=i
        #print(i)
        pointer[l][m][n]=pointer[l][m][n]+1
    #print(linkedcell)
    #print(antilink)
    #print(pointer)
    #count=0
    #print(currentacc)
    for i in range(N):
        #decide which cell the particle is located in
        l,m,n=antilink[i][0],antilink[i][1],antilink[i][2]
        #if (l<divisor-1 and m<divisor-1 and n<divisor-1 and l>0 and m>0 and n>0):
        for a in range(l-1,l+2):
            for b in range(m-1,m+2):
                for c in range(n-1,n+2):
                    #count=count+1
                    #print(count)
                    z1,z2,z3=a,b,c
                    #periodic boundary condition
                    if z1<0:
                        z1=divisor-1
                    if z1>divisor-1:
                        z1=0
                    if z2<0:
                        z2=divisor-1
                    if z2>divisor-1:
                        z2=0
                    if z3<0:
                        z3=divisor-1
                    if z3>divisor-1:
                        z3=0
                    #print("{} {} {}\n".format(z1,z2,z3))
                    for j in range (pointer[z1][z2][z3]):
                        position1=currentposition[i]
                        position2=currentposition[linkedcell[z1][z2][z3][j]]
                        if i==linkedcell[z1][z2][z3][j]:
                            continue
                        currentacc[i]=currentacc[i]+(1.0/mreal)*force(position1, \
                                                        position2,size,epsilon,sigma)
    #print(currentacc)
    #initialization finished        
    #start running (for p in range: biggest loop)
    jud=int(0)
    for p in range(rounds):
        print(p,end=' ')
        #print('round {}!'.format(p))
        #print('timestep={}'.format(timestep))
        #print(currentposition,end='\n')
        #print(currentvelocity,end='\n')
        #calculate nextposition
        storage[p]=currentposition
        nextposition=currentposition+currentvelocity*timestep+0.5* \
            np.power(timestep,2)*currentacc
        #calculate the displacement of each particle
        for i in range(len(nextposition)):
            displace[i]=displacement(nextposition[i],currentposition[i])
        lista=displace.tolist()
        lista.sort(reverse=True)
        displace=np.array(lista)
        #print(displace/sigma)
        for i in range(len(nextposition)):
            for j in range(3):
                if nextposition[i][j]>0.0:
                    while nextposition[i][j]>size*sigma:
                        nextposition[i][j]=nextposition[i][j]-size*sigma
                        jud=1
                else:
                    while nextposition[i][j]<0.0:
                        nextposition[i][j]=nextposition[i][j]+size*sigma    
        if (displace[0]+displace[1]>skin*sigma or (p%100==0 and p!=0) or jud==1):
            print("\n")
            #recreate linked cell
            linkedcell=np.full((divisor,divisor,divisor,N),-1,dtype='int32')
            #antilink: see which cell the particle is in
            antilink=np.zeros((N,3),dtype='int32')
            pointer=np.zeros((divisor,divisor,divisor),dtype='int32')
            for i in range(N):
                x,y,z=nextposition[i][0]/sigma,nextposition[i][1]/sigma, \
                        nextposition[i][2]/sigma
                l,m,n=int(math.floor(x/unitlength)),int(math.floor(y/unitlength)), \
                        int(math.floor(z/unitlength))
                antilink[i]=np.array([l,m,n])
                #print("{} {} {}".format(l,m,n))
                #print(currentposition[i]/sigma)
                #print(unitlength)
                #print('{} {} {}'.format(l,m,n))
                linkedcell[l][m][n][(pointer[l][m][n])]=i
                #print(i)
                pointer[l][m][n]=pointer[l][m][n]+1
            #print(linkedcell)
            #print(antilink)
            #print(pointer)
            #storage=np.zeros(N,dtype='int32')
            #count=0
            #print(currentacc)
            for i in range(N):
                #decide which cell the particle is located in
                l,m,n=antilink[i][0],antilink[i][1],antilink[i][2]
                #if (l<divisor-1 and m<divisor-1 and n<divisor-1 and l>0 and m>0 and n>0):
                for a in range(l-1,l+2):
                    for b in range(m-1,m+2):
                        for c in range(n-1,n+2):
                            #count=count+1
                            #print(count)
                            z1,z2,z3=a,b,c
                            #periodic boundary condition
                            if z1<0:
                                z1=divisor-1
                            if z1>divisor-1:
                                z1=0
                            if z2<0:
                                z2=divisor-1
                            if z2>divisor-1:
                                z2=0
                            if z3<0:
                                z3=divisor-1
                            if z3>divisor-1:
                                z3=0
                            #print("{} {} {}\n".format(z1,z2,z3))
                            for j in range (pointer[z1][z2][z3]):
                                position1=nextposition[i]
                                position2=nextposition[linkedcell[z1][z2][z3][j]]
                                if i==linkedcell[z1][z2][z3][j]:
                                    continue
                                nextacc[i]=nextacc[i]+(1.0/mreal)*force(position1, \
                                                                position2,size,epsilon,sigma)
                                sumpotential=sumpotential+LJpotential(epsilon,sigma, \
                                        distance(position1,position2,size,sigma),cutoff)
            sumpotential=sumpotential/2.0
                        
        else:
            for i in range(N):
                l,m,n=antilink[i][0],antilink[i][1],antilink[i][2]
                for a in range(l-1,l+2):
                    for b in range(m-1,m+2):
                        for c in range(n-1,n+2):
                            #count=count+1
                            #print(count)
                            z1,z2,z3=a,b,c
                            #periodic boundary condition
                            if z1<0:
                                z1=divisor-1
                            if z1>divisor-1:
                                z1=0
                            if z2<0:
                                z2=divisor-1
                            if z2>divisor-1:
                                z2=0
                            if z3<0:
                                z3=divisor-1
                            if z3>divisor-1:
                                z3=0
                            #print("{} {} {}\n".format(z1,z2,z3))
                            for j in range (pointer[z1][z2][z3]):
                                position1=nextposition[i]
                                position2=nextposition[linkedcell[z1][z2][z3][j]]
                                if i==linkedcell[z1][z2][z3][j]:
                                    continue
                                nextacc[i]=nextacc[i]+(1.0/mreal)*force(position1, \
                                                                position2,size,epsilon,sigma)
                                sumpotential=sumpotential+LJpotential(epsilon,sigma, \
                                        distance(position1,position2,size,sigma),cutoff)
            sumpotential=sumpotential/2.0
        #second,calculate nextvelocity, namely, velocity(t+delta(t))
        nextvelocity=currentvelocity+0.5*timestep*(currentacc+nextacc)
        #calculate sum of kinetic energy and temperature
        for i in range(len(nextvelocity)):
            sumkinetic=sumkinetic+kinetic(nextvelocity[i],mass)
        TEMP=2*sumkinetic/(3.0*N*kb)
        PRES=(2.0/3)*density*(sumkinetic/N)/1000.0
        sumenergy=sumkinetic+sumpotential
        eneperparticle=sumenergy/N
        #print(TEMP,end=' ')
        #figure preparation:
        Ytemp[p]=TEMP
        Ypres[p]=PRES
        Yene[p]=eneperparticle
        #renew the microstate:
        currentposition=nextposition
        currentvelocity=nextvelocity
        currentacc=nextacc
        nextacc=np.zeros((N,3),dtype='float64')
        displace=np.zeros(N,dtype='float64')
        sumkinetic,sumpotential,sumenergy=0.0,0.0,0.0
        jud=0
    plt.title("Temperature (N={}) initial T={}K\nε={}J σ={}m box size={}σ\n \
    molar mass={}g/mol cutoff={}σ\n".format(N,temperature,epsilon,sigma, \
                  size,mass,cutoff))
    plt.xlabel("the number of steps\ntime step={}".format(timestep))
    plt.ylabel("temperature (Kelvin)")
    plt.plot(X,Ytemp,'r-')
    plt.savefig("F:/MD V2/temp.jpg",dpi=500)
    plt.show()
    plt.title("Pressure (N={}) initial T={}K\nε={}J σ={}m box size={}σ\n \
    molar mass={}g/mol cutoff={}σ\n".format(N,temperature,epsilon,sigma, \
                  size,mass,cutoff))
    plt.xlabel("the number of steps\ntime step={}".format(timestep))
    plt.ylabel("Pressure (KPa)")
    plt.plot(X,Ypres,'b-')
    plt.savefig("F:/MD V2/pres.jpg",dpi=500)
    plt.show()
    plt.title("Energy per Particle (N={}) initial T={}K\nε={}J σ={}m box size={}σ\n \
    molar mass={}g/mol cutoff={}σ\n".format(N,temperature,epsilon,sigma, \
                  size,mass,cutoff))
    plt.xlabel("the number of steps\ntime step={}".format(timestep))
    plt.ylabel("Energy (Joule)")
    plt.plot(X,Yene,'g-')
    plt.savefig("F:/MD V2/ener.jpg",dpi=500)
    plt.show()
    
    for i in range(rounds):
        #3D display!
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        microstate=storage[i].swapaxes(0,1)
        ax.scatter(microstate[0]/sigma,microstate[1]/sigma,microstate[2]/sigma, \
                   color='r',marker='o')
        ax.set_title("Current Position, step {}".format(i+1))
        ax.set_xlabel('(sigma)')
        ax.set_ylabel('(sigma)')
        ax.set_zlabel('(sigma)')
        plt.savefig("F:/MD V2/box 3D/{}.jpg".format(i+1),dpi=500)
        plt.show()
    
#epsilon,sigma,mass,cutoff,skin,size,scale,N,temperature,rounds
main(1.67e-21,3.4e-10,39.9,2.0,0.6,6.525,0.001,100,83.95,10000)









