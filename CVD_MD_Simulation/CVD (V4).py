import numpy as np
import matplotlib.pyplot as plt
from mayavi import mlab
import math

system={'mode':1,'epsilon':0.0168,'sigma':2.2,'latticeCu':3.6149,'radiusCu':1.2781, \
        'radiusC':0.71,'number':5,'layerCu':2,'numx':6, 'numy':7,'numz':10,'cutoffCu':6.6, \
        'cutoffC':2.0,'De':6.325,'S':1.29,'beta':1.5,'Re':1.315,'R1':1.7, \
        'R2':2.0,'delta':0.80469,'a0':0.011304,'c0':19.0,'d0':2.5, \
        'initialtemp':300, 'molarmassC':12.011,'molarmassCu':63.546,'rounds':10000,'skin':1.0, \
        'timescale':1,'roughness':0.0,'bottom':2.0,'height':4.0,'initialloc':0.45, \
        'kb':1.3806e-23,'zone':3.0, \
        'angle':np.pi/3.0,'epsilonCC':0.0037,'sigmaCC':3.4,'defaultangle':np.pi/3.0, \
        'address':'F:/CVDresult/300K h=4.0 rand0.1 10000'}
#lattice parameter: Cu 361.49pm C(graphene):142pm
#temperature: kelvin  pressure:Pa  time: fs(1e-15s)
#velocity: angstrom/fs,  acc: angstrom/(fs^2)

def drawatom(center,radius,color):
    #create a point set in spherical coordinates
    dtheta,dphi=np.pi/40,np.pi/40
    theta,phi=np.mgrid[-np.pi-dtheta:np.pi+dtheta:dtheta, \
                     -np.pi-dphi:np.pi+dphi:dphi]
    #spherical coordinates-->cartesian coordinates
    x=radius*np.sin(theta)*np.cos(phi)+center[0]
    y=radius*np.sin(theta)*np.sin(phi)+center[1]
    z=radius*np.cos(theta)+center[2]
    atom=mlab.mesh(x,y,z)
    if color==0:
        #copper atom
        atom.module_manager.scalar_lut_manager.lut_mode = 'autumn'
    if color==1:
        #carbon atom
        atom.module_manager.scalar_lut_manager.lut_mode = 'winter'

def boxedge(boxsize):
    sizex,sizey,sizez=boxsize[0],boxsize[1],system['zone']*system['sigma']
    X=np.linspace(0,sizex,10)
    Y=np.linspace(0,sizey,10)
    Z=np.linspace(0,sizez,10)
    X1=np.full(10,sizex)
    Y1=np.full(10,sizey)
    Z1=np.full(10,sizez)
    Zero=np.zeros(10)
    
    mlab.plot3d(X,Zero,Zero)
    mlab.plot3d(X,Y1,Zero)
    mlab.plot3d(X,Y1,Z1)
    mlab.plot3d(X,Zero,Z1)
    
    mlab.plot3d(Zero,Y,Zero)
    mlab.plot3d(X1,Y,Zero)
    mlab.plot3d(X1,Y,Z1)
    mlab.plot3d(Zero,Y,Z1)
    
    mlab.plot3d(Zero,Zero,Z)
    mlab.plot3d(X1,Zero,Z)
    mlab.plot3d(X1,Y1,Z)
    mlab.plot3d(Zero,Y1,Z)
    
    numx,numy,numz=int(sizex),int(sizey),int(sizez)
    for i in range(numz):
        Z2=np.full(10,(i+1)*1.0)
        mlab.plot3d(X,Zero,Z2)
        mlab.plot3d(Zero,Y,Z2)
    for i in range(numx):
        X2=np.full(10,(i+1)*1.0)
        mlab.plot3d(X2,Zero,Z)
    for i in range(numy):
        Y2=np.full(10,(i+1)*1.0)
        mlab.plot3d(Zero,Y2,Z)

def display(boxsize,copperloc,carbonloc,name):
    print("\nLOADING...\n")
    try:
        engine = mayavi.engine
    except NameError:
        from mayavi.api import Engine
        engine = Engine()
        engine.start()
    if len(engine.scenes) == 0:
        engine.new_scene()
    j=0
    for i in range(len(copperloc)):
        drawatom(copperloc[i],system['radiusCu'],0)
    for i in range(len(carbonloc)):
        if carbonloc[i][2]<system['zone']*system['sigma']:
            drawatom(carbonloc[i],system['radiusC'],1)
            module_manager=engine.scenes[0].children[len(copperloc)+j].children[0].children[0]
            module_manager.scalar_lut_manager.reverse_lut=True
            j=j+1
    boxedge(boxsize)
    scene = engine.scenes[0]
    scene.scene.background = (0.0, 0.0, 0.0)
    scene.scene.show_axes = True
    scene.scene.magnification = 5
    scene.scene.camera.position = [10.441775588983637, 12.780601516556253, 62.211631014469155]
    scene.scene.camera.focal_point = [10.441775588983637, 12.780601516556253, 1.6299182559532115]
    scene.scene.camera.view_angle = 30.0
    scene.scene.camera.view_up = [0.0, 1.0, 0.0]
    scene.scene.camera.clipping_range = [50.03568296027767, 74.00311821873089]
    scene.scene.camera.compute_view_plane_normal()
    scene.scene.render()
    scene.scene.save('{}/{}1.jpg'.format(system['address'],name))
    scene.scene.isometric_view()
    scene.scene.camera.position = [34.331447175007646, 36.67027310258027, 25.519589841977226]
    scene.scene.camera.focal_point = [10.441775588983637, 12.780601516556253, 1.6299182559532115]
    scene.scene.camera.view_angle = 30.0
    scene.scene.camera.view_up = [0.0, 0.0, 1.0]
    scene.scene.camera.clipping_range = [5.607655074187303, 86.50558329057532]
    scene.scene.camera.compute_view_plane_normal()
    scene.scene.render()
    scene.scene.save('{}/{}2.jpg'.format(system['address'],name))
    scene.scene.camera.position = [49.6033084103893, 13.220646972407467, 14.983941340364805]
    scene.scene.camera.focal_point = [10.441775588983637, 12.780601516556253, 1.6299182559532115]
    scene.scene.camera.view_angle = 30.0
    scene.scene.camera.view_up = [-0.3227273826678203, -0.002145387394995258, 0.9464895317906605]
    scene.scene.camera.clipping_range = [15.385615742497164, 74.19714536388366]
    scene.scene.camera.compute_view_plane_normal()
    scene.scene.render()
    scene.scene.save('{}/{}3.jpg'.format(system['address'],name))

def distance(position1,position2,boxsize):
    sizex,sizey=boxsize[0],boxsize[1]
    dist=np.abs(position1-position2)
    while dist[0]>sizex*0.5:
        dist[0]=np.abs(dist[0]-sizex)
    while dist[1]>sizey*0.5:
        dist[1]=np.abs(dist[1]-sizey)
    distance=np.sqrt(dist[0]**2+dist[1]**2+dist[2]**2)
    return distance

def displacement(position1,position2):
    #calculate the displacement of particles
    #absolute value of displacement, no PBC involved
    vector=position2-position1
    displace=np.sqrt(vector[0]**2+vector[1]**2+vector[2]**2)
    return displace

def Gtheta(i,j,k,carbonloc,boxsize):
    a0,c0,d0=system['a0'],system['c0'],system['d0']
    sizex,sizey=boxsize[0],boxsize[1]
    positioni,positionj,positionk=np.zeros(3),np.zeros(3),np.zeros(3)
    for m in range(3):
        positioni[m]=carbonloc[i][m]
        positionj[m]=carbonloc[j][m]
        positionk[m]=carbonloc[k][m]
    vectorij=positionj-positioni
    vectorik=positionk-positioni
    if vectorij[0]>0:
        while vectorij[0]>sizex*0.5:
            positionj[0]=positionj[0]-sizex
            vectorij=positionj-positioni
    else:
        while vectorij[0]<-sizex*0.5:
            positionj[0]=positionj[0]+sizex
            vectorij=positionj-positioni
    if vectorij[1]>0:
        while vectorij[1]>sizey*0.5:
            positionj[1]=positionj[1]-sizey
            vectorij=positionj-positioni
    else:
        while vectorij[1]<-sizey*0.5:
            positionj[1]=positionj[1]+sizey
            vectorij=positionj-positioni
    if vectorik[0]>0:
        while vectorik[0]>sizex*0.5:
            positionk[0]=positionk[0]-sizex
            vectorik=positionk-positioni
    else:
        while vectorik[0]<-sizex*0.5:
            positionk[0]=positionk[0]+sizex
            vectorik=positionk-positioni
    if vectorik[1]>0:
        while vectorik[1]>sizey*0.5:
            positionk[1]=positionk[1]-sizey
            vectorik=positionk-positioni
    else:
        while vectorik[1]<-sizey*0.5:
            positionk[1]=positionk[1]+sizey
            vectorik=positionk-positioni
    modij=np.sqrt(vectorij[0]**2+vectorij[1]**2+vectorij[2]**2)
    modik=np.sqrt(vectorik[0]**2+vectorik[1]**2+vectorik[2]**2)
    costheta=np.sum(vectorij*vectorik)/(modij*modik)
    theta=np.arccos(costheta)
    costheta=np.cos(theta+system['angle'])
    Gtheta=a0*(1.0+(c0/d0)**2-c0**2/(d0**2+(1+costheta)**2))
    return Gtheta

def cutfunction(rxx):
    #cutoff function in the multibody potential
    fcut=int(0) 
    if rxx<=system['R1']:
        fcut=1.0
    elif rxx>system['R1'] and rxx<system['R2']:
        fcut=0.5*(1.0+np.cos((np.pi*(rxx-system['R1']))/(system['R2']-system['R1'])))
    else:
        fcut=int(0)
    return fcut

def dercut(rxx):
    #derivative of cutfuction
    dcut=0.0
    R1,R2=system['R1'],system['R2']
    if rxx<=R1:
        dcut=0.0
    elif rxx>R1 and rxx<R2:
        dcut=-0.5*np.sin((np.pi/(R2-R1))*(rxx-R1))*(np.pi/(R2-R1))
    else:
        dcut=0.0
    return dcut

def substrate():
    A0,number,layer=system['latticeCu'],system['number'],system['layerCu']
    #Cu FCC substsate (distance unit: angstroms)
    #lattice parameter: a0  (sqrt(2)/2)*a0: distance between closest neighbors
    sizex,sizey=number*(0.5*np.sqrt(6.0)*A0),number*np.sqrt(2.0)*A0
    sizez=(sizex+sizey)*system['height']
    #draw the edges of the simulation box
    rawsize=np.array([sizex,sizey,sizez])
    vector1=A0*np.array([np.sqrt(6)/2,0.0,0.0])
    vector2=A0*np.array([0.0,np.sqrt(2)/2,0.0])
    vector3=A0*np.array([0.0,0.0,-np.sqrt(3.0)])
    pointa=A0*np.array([np.sqrt(6.0)/6,0.0,np.sqrt(3.0)/3])
    point0=A0*np.array([0.0,0.0,0.0])
    point1=A0*np.array([np.sqrt(6.0)/12,np.sqrt(2.0)/4,-np.sqrt(3.0)/3])
    point2=A0*np.array([np.sqrt(6.0)/6,0.0,-2*np.sqrt(3.0)/3])
    startpoint=np.array([point0,point1,point2])
    startpoint1=np.array([point0+vector1/2.0+vector2/2.0,point1+vector1/2.0+vector2/2.0, \
                          point2+vector1/2.0+vector2/2.0])
    copperloc=list()
    for i in range(number):
        for j in range(int(2*number)):
            for k in range(layer):
                temp0,temp1=k//3,k%3
                start=temp0*vector3+startpoint[temp1]
                start1=temp0*vector3+startpoint1[temp1]
                point=start+i*vector1+j*vector2
                point1=start1+i*vector1+j*vector2
                if point[0]<=sizex+startpoint[temp1][0]-0.9999*np.sqrt(3.0)*A0/2 and \
                        point[1]<=sizey+startpoint[temp1][1]-0.9999*A0/2 and \
                        point[0]<sizex and point[1]<sizey:
                    copperloc.append(point)
                if point1[0]<=sizex+startpoint1[temp1][0]-0.99999*np.sqrt(3.0)*A0/2 and \
                        point1[1]<=sizey+startpoint1[temp1][1]-0.9999*A0/2 and \
                        point1[0]<sizex and point[1]<sizey:
                    copperloc.append(point1)
    start=pointa
    start1=pointa+vector1/2.0+vector2/2.0
    for i in range(number):
        for j in range(int(2*number)):
            random=np.random.rand()
            point=start+i*vector1+j*vector2
            point1=start1+i*vector1+j*vector2
            if point[0]<=sizex+start[0]-0.9999*np.sqrt(3.0)*A0/2 and \
                        point[1]<=sizey+start[1]-0.9999*A0/2 and \
                        point[0]<sizex and point[1]<sizey and random>(1.0-system['roughness']):
                copperloc.append(point)
            if point1[0]<=sizex+start1[0]-0.99999*np.sqrt(3.0)*A0/2 and \
                    point1[1]<=sizey+start1[1]-0.9999*A0/2 and \
                    point1[0]<sizex and point[1]<sizey and random>(1.0-system['roughness']):
                copperloc.append(point1)
    copperloc=np.array(copperloc)
    print("The initialization of Cu substrate is finished.\n")
    return copperloc,rawsize

def carbon(rawsize):
    #mass: molar mass, g/mol
    #time: fs(1e-15s), distance: angstrom
    zbase,cutoffC,initialtemp,mass=system['bottom']*system['sigma'], \
        system['cutoffC'],system['initialtemp'],system['molarmassC']
    N=int(system['numx']*system['numy']*system['numz'])
    carbonloc=list()
    kb=1.380649e-23
    dev=np.sqrt(kb*initialtemp/(mass*1e-3/(6.02e23)))
    carbonvel=1e-5*np.random.normal(0,dev,(N,3))
    for i in range(N):
        if carbonvel[i][2]>0:
            carbonvel[i]=-carbonvel[i]
    sizex,sizey,sizez=rawsize[0],rawsize[1],rawsize[2]
    neighborC=np.full((N,N),-1,dtype='int32')
    numneighborC=np.zeros(N,dtype='int32')
    #initialize the positions of carbon atoms
    unitx,unity,unitz=sizex/system['numx'],sizey/system['numy'], \
        (sizez-zbase)/(system['numz'])
    for i in range(system['numx']):
        for j in range(system['numy']):
            for k in range(system['numz']):
                array1=np.array([unitx*(i+system['initialloc']),unity*(j+system['initialloc']), \
                                 unitz*k+zbase])
                randarray=np.array([unitx*(1.0-2*system['initialloc']), \
                                    unity*(1.0-2*system['initialloc']),0.0])*np.random.rand(3)
                position=array1+randarray
                carbonloc.append(position)
    carbonloc=np.array(carbonloc)
    boxsize=np.array([sizex,sizey,sizez])
    print('sizex={:.2f}  sizey={:.2f}  sizez={:.2f}'.format(sizex,sizey,sizez));
    print('Ncarbon={}'.format(system['numx']*system['numy']*system['numz']));
    #create neighborlist
    for i in range(N):
        for j in range(i+1,N):
            position1,position2=carbonloc[i],carbonloc[j]
            dist=distance(position1,position2,boxsize)
            if system['mode']==1:
                if dist<=cutoffC:
                    neighborC[i][numneighborC[i]]=j
                    neighborC[j][numneighborC[j]]=i
                    numneighborC[i]=numneighborC[i]+1
                    numneighborC[j]=numneighborC[j]+1
            if system['mode']==2:
                if dist<=2.0*system['sigmaCC']:
                    neighborC[i][numneighborC[i]]=j
                    neighborC[j][numneighborC[j]]=i
                    numneighborC[i]=numneighborC[i]+1
                    numneighborC[j]=numneighborC[j]+1
    print("The distribution of carbon atoms is finished.\n")
    return boxsize,carbonloc,carbonvel,neighborC,numneighborC

def LJpotential(atom,distance):
    if atom==1:
        #Cu-C interaction
        epsilon,sigma,cutoff=system['epsilon'],system['sigma'],system['cutoffCu']
    if atom==2:
        #C-C interaction
        epsilon,sigma,cutoff=system['epsilonCC'],system['sigmaCC'],2.0*system['sigmaCC']
    #distance: Angstroms
    #cutoff: cutoff radius, Angstroms
    x=sigma/distance
    xc=sigma/cutoff
    Urc=4*epsilon*(np.power(xc,12)-np.power(xc,6))
    if distance<=cutoff:
        Ur=4*epsilon*(np.power(x,12)-np.power(x,6))-Urc
    else:
        Ur=0.0
    return Ur

def LJacc(atom,boxsize,position1,position2):
    #calculate the interactive force between two particles
    #the interactive force is a vector, the output is a vector
    #position1:studied particle (carbon atoms)
    #position2: the particle generating the force field (copper atoms)
    sizex,sizey=boxsize[0],boxsize[1]
    if atom==1:
        #Cu-C interaction
        epsilon,sigma=system['epsilon'],system['sigma']
    if atom==2:
        #C-C interaction
        epsilon,sigma=system['epsilonCC'],system['sigmaCC']
    vector=position1-position2
    if vector[0]>=0:
        while vector[0]>sizex/2.0:
            vector[0]=vector[0]-sizex
    else:
        while vector[0]<-sizex/2.0:
            vector[0]=vector[0]+sizex
    if vector[1]>=0:
        while vector[1]>sizey/2.0:
            vector[1]=vector[1]-sizey
    else:
        while vector[1]<-sizey/2.0:
            vector[1]=vector[1]+sizey
            
    dist=np.sqrt(vector[0]**2+vector[1]**2+vector[2]**2)
    x=sigma/dist
    absforce=4*epsilon*(12*np.power(x,11)-6*np.power(x,5))*(sigma/(dist**2))
    xforce=absforce*(vector[0]/dist)
    yforce=absforce*(vector[1]/dist)
    zforce=absforce*(vector[2]/dist)
    force=np.array([xforce,yforce,zforce])
    acc=(force/(system['molarmassC']*0.001/(6.02e23)))*(1.602e-9/1e20)
    return acc

def multipotential(boxsize,i,j,carbonloc,neighborC,numneighborC):
    positioni,positionj=carbonloc[i],carbonloc[j]
    multipot,Bij,Bji=0.0,0.0,0.0
    #energy: eV distance: angstrom
    rij=distance(positioni,positionj,boxsize)
    fcut=cutfunction(rij)
    #VR:repulsive term  VA:attractive term
    #neighborC=np.full((N,N),-1,dtype='int32')
    #numneighborC=np.zeros(N,dtype='int32')
    if fcut!=0:
        VR=(system['De']/(system['S']-1.0))*fcut* \
            np.exp((-np.sqrt(2.0*system['S']))*system['beta']*(rij-system['Re']))
        VA=((system['De']*system['S'])/(system['S']-1.0))*fcut* \
            np.exp((-np.sqrt(2.0/system['S']))*system['beta']*(rij-system['Re']))
        for p in range(numneighborC[i]):
            k=neighborC[i][p]
            if k==j:
                continue
            positionk=carbonloc[k]
            rik=distance(positioni,positionk,boxsize)
            Bij=Bij+Gtheta(i,j,k,carbonloc,boxsize)*cutfunction(rik)
        Bij=pow(Bij+1.0,-system['delta'])
        for p in range(numneighborC[j]):
            k=neighborC[j][p]
            if k==i:
                continue
            positionk=carbonloc[k]
            rjk=distance(positionj,positionk,boxsize)
            Bji=Bji+Gtheta(j,i,k,carbonloc,boxsize)*cutfunction(rjk)
        Bji=pow(Bji+1.0,-system['delta'])
        B=(Bij+Bji)/2.0
        multipot=VR-B*VA
    return multipot

def multiacc(boxsize,i,j,carbonloc,neighborC,numneighborC):
    #i: studied atom  j: other atom
    sizex,sizey=boxsize[0],boxsize[1]
    De,S,beta,Re=system['De'],system['S'],system['beta'],system['Re']
    multiforce=np.array([0.0,0.0,0.0])
    Bij,Bji=0.0,0.0
    positioni,positionj=carbonloc[i],carbonloc[j]
    vector=positioni-positionj
    if vector[0]>=0:
        while vector[0]>sizex/2.0:
            vector[0]=vector[0]-sizex
    else:
        while vector[0]<-sizex/2.0:
            vector[0]=vector[0]+sizex
    if vector[1]>=0:
        while vector[1]>sizey/2.0:
            vector[1]=vector[1]-sizey
    else:
        while vector[1]<-sizey/2.0:
            vector[1]=vector[1]+sizey
    rij=np.sqrt(vector[0]**2+vector[1]**2+vector[2]**2)
    if cutfunction(rij)!=0:
        factor0=De/(S-1.0)
        factor1=np.exp(-np.sqrt(2.0*S)*beta*(rij-Re))
        factor2=np.exp(-np.sqrt(2.0/S)*beta*(rij-Re))
        DVR=factor0*((-np.sqrt(2.0*S)*beta*factor1)*cutfunction(rij)+factor1*dercut(rij))
        DVA=factor0*((-np.sqrt(2.0/S)*beta*factor2)*cutfunction(rij)+factor2*dercut(rij))
        for p in range(numneighborC[i]):
            k=neighborC[i][p]
            if k==j:
                continue
            positionk=carbonloc[k]
            rik=distance(positioni,positionk,boxsize)
            Bij=Bij+Gtheta(i,j,k,carbonloc,boxsize)*cutfunction(rik)
        Bij=pow(Bij+1.0,-system['delta'])
        for p in range(numneighborC[j]):
            k=neighborC[j][p]
            if k==i:
                continue
            positionk=carbonloc[k]
            rjk=distance(positionj,positionk,boxsize)
            Bji=Bji+Gtheta(j,i,k,carbonloc,boxsize)*cutfunction(rjk)
        Bji=pow(Bji+1.0,-system['delta'])
        B=(Bij+Bji)/2.0
        absmultiforce=-(DVR-B*DVA)
        multiforce[0]=absmultiforce*(vector[0]/rij)
        multiforce[1]=absmultiforce*(vector[1]/rij)
        multiforce[2]=absmultiforce*(vector[2]/rij)
    acc=(multiforce/(system['molarmassC']*0.001/(6.02e23)))*(1.602e-9/1e20)
    return acc

def kinetic(velocity,mass):
    velocity=velocity*1e5
    #calculate the kinetic energy of a particle
    #mass: molar mass
    mreal=mass*(1e-3)/(6.02e23)
    kinetic=0.5*mreal*(velocity[0]**2+velocity[1]**2+velocity[2]**2)
    kinetic=kinetic/(1.602e-19)
    return kinetic

def initialization():
    #lattice parameter: Cu 361.49pm C(graphene):142pm      
    copperloc,rawsize=substrate()
    boxsize,carbonloc,carbonvel,neighborC,numneighborC=carbon(rawsize)
    Ncarbon,Ncopper=len(carbonloc),len(copperloc)
    #neighborCu: the neighboring copper atoms of each carbon atom
    #numneighbor: the number of neighboring copper atoms of each carbon atom
    neighborCu=np.full((Ncarbon,Ncopper),-1,dtype='int32')
    numneighborCu=np.zeros(Ncarbon,dtype='int32')
    #carbonacc: acceleration of all carbon atoms
    carbonacc=np.zeros((Ncarbon,3),dtype='float64')
    print("Ncopper={}  Ncarbon={}".format(Ncopper,Ncarbon))
    #how many neighboring copper atom each carbon atom has
    for i in range(Ncarbon):
        for j in range(Ncopper):
            positionC,positionCu=carbonloc[i],copperloc[j]
            dist=distance(positionC,positionCu,boxsize)
            if dist<system['cutoffCu']:
                neighborCu[i][numneighborCu[i]]=j
                #print(neighborCu[numneighbor[i]])
                numneighborCu[i]=numneighborCu[i]+1
    for i in range(Ncarbon):
        for p in range(numneighborC[i]):
            #C-C interaction
            j=neighborC[i][p]
            if system['mode']==1:
                carbonacc[i]=carbonacc[i]+multiacc(boxsize,i,j,carbonloc,neighborC, \
                                                 numneighborC)
            if system['mode']==2:
                carbonacc[i]=carbonacc[i]+LJacc(2,boxsize,carbonloc[i],carbonloc[j])
        for q in range(numneighborCu[i]):
            #Cu-C interaction
            h=neighborCu[i][q]
            dist=distance(carbonloc[i],copperloc[h],boxsize)
            carbonacc[i]=carbonacc[i]+LJacc(1,boxsize,carbonloc[i],copperloc[h])
    display(boxsize,copperloc,carbonloc,'initial ')
    return boxsize,copperloc,carbonloc,carbonvel,carbonacc, \
           neighborC,numneighborC,neighborCu,numneighborCu

def main():
    boxsize,copperloc,carbonloc,carbonvel,carbonacc, \
        neighborC,numneighborC,neighborCu,numneighborCu=initialization()
    rounds,dT,skin=system['rounds'],system['timescale'],system['skin']
    sizex,sizey,sizez=boxsize[0],boxsize[1],boxsize[2]
    N,NCu=len(carbonloc),len(copperloc)
    potentialCC,potentialCuC=0.0,0.0
    sumkinetic,TEMP=0.0,0.0
    jud=0
    
    #see which particle is on the surface
    surface=list()
    #see which particle is in the reacting zone
    countreact,supplement=int(0),int(0)
    Nconstant=countreact
    numinside=np.zeros(rounds,dtype='int32')
    
    storagepos=np.zeros((rounds,N,3))
    storagepot=np.zeros((rounds),dtype='float64')
    storagetemp=np.zeros((rounds),dtype='float64')
    currentposition=np.zeros((N,3),dtype='float64')
    nextposition=np.zeros((N,3),dtype='float64')
    currentvelocity=np.zeros((N,3),dtype='float64')
    nextvelocity=np.zeros((N,3),dtype='float64')
    currentacc=np.zeros((N,3),dtype='float64')
    nextacc=np.zeros((N,3),dtype='float64')
    currentstatus=np.zeros(N,dtype='int32')
    displace=np.zeros(N,dtype='float64')
    for i in range(N):
        for j in range(len(carbonloc[0])):
            currentposition[i][j]=carbonloc[i][j]
            currentvelocity[i][j]=carbonvel[i][j]
            currentacc[i][j]=carbonacc[i][j]
        if currentposition[i][2]<=1.1*system['sigma']:
            currentstatus[i]=0
            surface.append(i)
        elif currentposition[i][2]<system['zone']*system['sigma']:
            currentstatus[i]=1
            countreact=countreact+1
        else:
            currentstatus[i]=2
            currentacc[i]=np.array([0.0,0.0,0.0])
    Nconstant=countreact
    print('Nconstant={}'.format(Nconstant))
    contreact2=countreact;supplement2=supplement;
    #print(Nconstant)
    #display(boxsize,copperloc,currentposition)
    print("Start Running!")
    for m2 in range(rounds):
        print(m2,end='   ')
        storagepos[m2]=currentposition
        nextposition=currentposition+currentvelocity*dT+0.5*(dT**2)*currentacc
        for i in range(len(nextposition)):
            if currentstatus[i]==2:
                continue
            displace[i]=displace[i]+displacement(nextposition[i],currentposition[i])
        for i in range(len(nextposition)):
            if nextposition[i][0]>0.0:
                while nextposition[i][0]>sizex:
                    nextposition[i][0]=nextposition[i][0]-sizex
            else:
                while nextposition[i][0]<0.0:
                    nextposition[i][0]=nextposition[i][0]+sizex
            if nextposition[i][1]>0.0:
                while nextposition[i][1]>sizey:
                    nextposition[i][1]=nextposition[i][1]-sizey
            else:
                while nextposition[i][1]<0.0:
                    nextposition[i][1]=nextposition[i][1]+sizey
            if nextposition[i][2]>sizez:
                nextposition[i][2]=sizez
                currentvelocity[i][2]=-currentvelocity[i][2]
        if (1.0*np.max(displace)>skin or (m2%1000==0 and m2!=0)):
            print(np.max(displace))
            if 2.0*np.max(displace)>skin:
                print("A")
            if m2%1000==0 and m2!=0:
                print("B")
            displace=np.zeros(N,dtype='float64')
            print("\n")
            neighborC=np.full((N,N),-1,dtype='int32')
            numneighborC=np.zeros(N,dtype='int32')
            neighborCu=np.full((N,NCu),-1,dtype='int32')
            numneighborCu=np.zeros(N,dtype='int32')
            for i in range(N):
                if currentstatus[i]==2:
                    continue
                for j in range(i+1,N):
                    if currentstatus[j]==2:
                        continue
                    position1,position2=nextposition[i],nextposition[j]
                    dist=distance(position1,position2,boxsize)
                    if system['mode']==1:
                        if dist<=system['cutoffC']:
                            neighborC[i][numneighborC[i]]=j
                            neighborC[j][numneighborC[j]]=i
                            numneighborC[i]=numneighborC[i]+1
                            numneighborC[j]=numneighborC[j]+1
                    if system['mode']==2:
                        if dist<=2.0*system['sigmaCC']:
                            neighborC[i][numneighborC[i]]=j
                            neighborC[j][numneighborC[j]]=i
                            numneighborC[i]=numneighborC[i]+1
                            numneighborC[j]=numneighborC[j]+1
            for i in range(N):
                if currentstatus[i]==2:
                    continue
                for j in range(NCu):
                    positionC,positionCu=nextposition[i],copperloc[j]
                    dist=distance(positionC,positionCu,boxsize)
                    if dist<system['cutoffCu']:
                        neighborCu[i][numneighborCu[i]]=j
                        numneighborCu[i]=numneighborCu[i]+1
        for i in range(N):
            if currentstatus[i]==2:
                continue
            if countreact>0:
                if currentstatus[i]==1:
                    sumkinetic=sumkinetic+kinetic(currentvelocity[i],system['molarmassC'])
                    countreact=countreact-1
            if countreact==0:
                if currentstatus[i]==0 and supplement>0:
                    sumkinetic=sumkinetic+kinetic(currentvelocity[i],system['molarmassC'])
                    supplement=supplement-1
            for p in range(numneighborC[i]):
                j=neighborC[i][p]
                if system['mode']==1:
                    potentialCC=potentialCC+multipotential(boxsize,i,j,nextposition,neighborC, \
                                                              numneighborC)
                    nextacc[i]=nextacc[i]+multiacc(boxsize,i,j,nextposition,neighborC, \
                                                             numneighborC)
                if system['mode']==2:
                    potentialCC=potentialCC+LJpotential(2,distance(nextposition[i],nextposition[j], \
                                                             boxsize))
                    nextacc[i]=nextacc[i]+LJacc(2,boxsize,nextposition[i],nextposition[j])
            for q in range(numneighborCu[i]):
                h=neighborCu[i][q]
                dist=distance(carbonloc[i],copperloc[h],boxsize)
                potentialCuC=potentialCuC+LJpotential(1,dist)
                nextacc[i]=nextacc[i]+LJacc(1,boxsize,nextposition[i],copperloc[h])
        TEMP=(sumkinetic*1.602e-19*2.0)/(3*Nconstant*system['kb'])
        Lambda=np.sqrt(system['initialtemp']/TEMP)
        totalpot=potentialCC+potentialCuC
        storagepot[m2]=totalpot
        storagetemp[m2]=TEMP
        nextvelocity=currentvelocity+0.5*dT*(currentacc+nextacc)
        #clear everything up
        if len(surface)!=0:
            randchoice=np.random.randint(0,len(surface),supplement2)
            jud=1
        for i in range(len(nextvelocity)):
            for j in range(len(nextvelocity[0])):
                currentposition[i][j]=nextposition[i][j]
                currentvelocity[i][j]=nextvelocity[i][j]
                currentacc[i][j]=nextacc[i][j]
                if j==2:
                    if currentposition[i][2]<system['sigma']*1.1 and (i not in surface):
                        surface.append(i)
                        currentstatus[i]=0
                        absvel1=np.sqrt(currentvelocity[i][0]**2+currentvelocity[i][1]**2+currentvelocity[i][2]**2)
                        absvel2=np.sqrt(currentvelocity[i][0]**2+currentvelocity[i][1]**2)
                        currentvelocity[i][2]=0.0
                        currentvelocity[i]=currentvelocity[i]*(absvel1/absvel2)
                        if jud==1:
                            if i in randchoice:
                                currentvelocity=currentvelocity*Lambda
                    elif currentposition[i][2]<system['zone']*system['sigma']:
                        currentstatus[i]=1
                        countreact=countreact+1
                        currentvelocity[i]=currentvelocity[i]*Lambda
                    else:
                        currentstatus[i]=2
        #numinside[m2]=len(surface)+countreact
        numinside[m2]=len(surface)
        if countreact>Nconstant:
            countreact=Nconstant;supplement=0;
            countreact2=countreact;supplement2=supplement;
        else:
            supplement=Nconstant-countreact
            countreact2=countreact;supplement2=supplement;
        nextacc=np.zeros((N,3),dtype='float64')
        potentialCC,potentialCuC,sumkinetic=0.0,0.0,0.0
        jud=0
        if m2%500==0 and m2>0:
            display(boxsize,copperloc,currentposition,'round {} '.format(m2))
    display(boxsize,copperloc,currentposition,'final ')
    return storagepos,storagepot,storagetemp,numinside

def figuretemp(storagetemp):
    X=np.arange(len(storagetemp))
    fig,ax=plt.subplots(figsize=(9,8))
    ax.plot(X,storagetemp)
    ax.set_xlabel("the number of steps\none step = {:.2f} fs".format(system['timescale']), \
                  fontsize=15)
    ax.set_ylabel("temperature (kelvin)",fontsize=15)
    ax.set_title("Temperature\nT={:.1f}K\n".format(system['initialtemp']),fontsize=20)
    plt.savefig("{}/CVDtemp.jpg".format(system['address']),dpi=500)
    plt.show()
    
def figurenum(numinside):
    X=np.arange(len(numinside))
    num=np.zeros(len(numinside),dtype='int32')
    num[0]=0
    for i in range(len(numinside)-1):
        num[i+1]=numinside[i+1]-numinside[i]
    fig,ax=plt.subplots(figsize=(9,8))
    ax.plot(X,num,'r-',markersize=1)
    ax.set_xlabel("the number of steps\none step = {:.2f} fs  T={:.1f}K".format(system['timescale'],system['initialtemp']), \
                  fontsize=15)
    ax.set_ylabel("deposition rate (per step)",fontsize=15)
    ax.set_title("Deposition Rate (per step)\n average={:.3f}\n".format(np.mean(num)),fontsize=20)
    plt.savefig("{}/CVDnumber.jpg".format(system['address']),dpi=500)
    plt.show()
    
#initialization()

storagepos,storagepot,storagetemp,numinside=main()
figuretemp(storagetemp)
figurenum(numinside)











