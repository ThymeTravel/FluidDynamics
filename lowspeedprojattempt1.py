import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import math
import random
import time
import numpy as pi
import sympy
from sympy.abc import K, x, y


#basil aranda is a G. wish me luck
#So we beat on, boats against the current, 
#borne back ceaselessly into the past. - F. SCOTT FITZGERALD
#Set variables

N = 4;
U = 1;
a = 1;
alpha = 0;
Vinfinity = 1;
R = 1;
iteration = 0;
iterationZetaEta = 0;
iterationThetaI = 0;
iterationVi = 0;
iterationVj = 0;
iterationV_n = 0;
iterationB = 0;
iterationCp = 0;
i = 0;
j = 0;

K = 1;

upperGamma = K*(2*math.pi);



x = [0]*N;
y = [0]*N;
zeta = [0]*N;
eta = [0]*N;
theta = [0]*N;
V_n = [0]*N;
B = [0]*N;
lambdasolve = [0]*N;
Cp = [0]*N;
Vtan = [0]*N;
CbarAsolve = [0]*N;
CpMaybe = [0]*N;


phi = np.zeros((N,N));
r = np.zeros((N,N));
C = np.zeros((N,N));
Cbar = np.zeros((N,N));
lambdaLow = np.zeros((N,N));
CbarA = np.zeros((N,N));


assignI = 0;
beta = 2*math.pi/N;
alpha = 0;

#Code Start

print("Basil Aranda Fluids Project");

while(assignI < N):

    x[assignI] = -1*(R * math.cos((2*assignI+1)*beta/2))
    
    y[assignI] = (R * math.sin((2*assignI+1)*beta/2))
    
    assignI += 1;

x = np.around(x, decimals = 8);
y = np.around(y, decimals = 8);

#Panel Length
ds = 2*R*math.sin(beta/2);
#print(ds)
#print(x)
#print(y)

#x1, y1 = [-1, 12], [1, 4]
#x2, y2 = [1, 10], [3, 2]

#while iteration < N-1:
    #plt.plot(x[iteration], y[iteration], x[iteration + 1], y[iteration + 1], linewidth = 2.0, marker = 'o')
    #iteration += 1;
'''
plt.plot(x, y, 'b')
plt.plot([x[3], y[0]],[x[0], y[3]], 'b')
angle = np.linspace(0 , 2 * math.pi) 
 
radius = 1
 
x = radius * np.cos( angle ) 
y = radius * np.sin( angle ) 
 
plt.plot(x, y)
plt.title( 'Somehow Working' ) 
plt.grid()
plt.show() 
'''
omega = np.arange(2*math.pi/N, 2*math.pi+math.pi/(2*N),  2*math.pi/N);
#print(omega)

while(iterationZetaEta < N-1):
    zeta[iterationZetaEta] = 1/2*(x[iterationZetaEta] + x[iterationZetaEta +1]);
    eta[iterationZetaEta] = 1/2*(y[iterationZetaEta] + y[iterationZetaEta +1]);
    iterationZetaEta +=1;
zeta[N-1] = 1/2*(x[0] + x[N-1]);
eta[N-1] = 1/2*(y[0] + y[N-1]);

#print(zeta)
#print(eta)

while (iterationThetaI < N-1):
    theta[iterationThetaI] = math.atan2((y[iterationThetaI+1]-y[iterationThetaI]),
    (x[iterationThetaI+1]-x[iterationThetaI]));
    iterationThetaI +=1;
theta[N-1] = math.atan2((y[0]-y[N-1]),(x[0]-x[N-1]));

while(i < N):
    j = 0;

    while(j < N):
        phi[i,j] = math.atan2((eta[j]-eta[i]),(zeta[j]-zeta[i]));
        r[i,j] = math.sqrt((zeta[j]-zeta[i])**2+(eta[j]-eta[i])**2);
        if(i!=j):
            C[i,j] = ((math.sin(theta[i]-phi[i,j])/(2*math.pi*r[i,j]))*ds);
            Cbar[i,j] = ((math.cos(theta[i]-phi[i,j])/(2*math.pi*r[i,j]))*ds);
        else:
            C[i,j] = 0.5;
            Cbar[i,j] = 0;
        j+=1;
    i+=1;

#print(theta)



while(iterationVi < N):

    iterationVj = 0;
    while(iterationVj < N):
        if(iterationVi == iterationVj):
            lambdaLow[iterationVi, iterationVj] = -0.5;
            #CbarA[iterationVi,iterationVj] = 0;
        else:
            lambdaLow[iterationVi, iterationVj] = ds/(2*math.pi)*math.sin(phi[iterationVi][iterationVj] - theta[iterationVi])/r[iterationVi][iterationVj];
            #CbarA[iterationVi,iterationVj] = ds/(2*math.pi)*math.cos(phi[iterationVi][iterationVj] - theta[iterationVi])/r[iterationVi][iterationVj];
        iterationVj +=1;


    iterationVi+=1;



while(iterationV_n < N):
    V_n[iterationV_n] = -1 * Vinfinity * math.sin(theta[iterationV_n] - alpha);
    B[iterationV_n] = -1 * Vinfinity * math.cos(theta[iterationV_n] - alpha);
    iterationV_n +=1;



V_n = np.around(V_n, decimals = 3);
B = np.around(B, decimals = 3);


lambdasolve = np.linalg.solve(lambdaLow, V_n);
lambdasolve= np.around(lambdasolve, decimals = 8);


#print(lambdasolve)

while(iterationB < N):
    Vtan[iterationB] =  -2 * math.sin(theta[iterationB]) + upperGamma/(2*math.pi);
    iterationB += 1;


print(Vtan)



iterationVi = 0; iterationVj = 0;





while(iterationVi < N):

    iterationVj = 0;
    while(iterationVj < N):
        if(iterationVi == iterationVj):
           
            CbarA[iterationVi,iterationVj] = 0;
        else:
           
            CbarA[iterationVi,iterationVj] = lambdasolve[iterationVj]*ds/(2*math.pi)*math.cos(phi[iterationVi][iterationVj] - theta[iterationVi])/r[iterationVi][iterationVj];
        iterationVj +=1;


    iterationVi+=1;




CbarAsolve = np.sum(CbarA, axis = 1);
CbarAsolve= np.around(CbarAsolve, decimals = 8);


VMaybe = np.add(CbarAsolve, B)



print(B)
print(VMaybe)
#CbarAsolve = np.linalg.solve(CbarA, B);
#CbarAsolve= np.around(CbarAsolve, decimals = 8);

#print(Vtan)

while(iterationCp < N):
    
    Cp[iterationCp] = 1 - (Vtan[iterationCp]/Vinfinity)**2;
    CpMaybe[iterationCp] = 1 - VMaybe[iterationCp]**2;
    iterationCp += 1;


#print(r)
print(Cp)
print(CpMaybe)

omegaDeg = np.arange(0,360,1);
iterationOmegaDeg = 0;

CpExact = [0]*360
CpVTAN = [0]*N;
#print(omegaDeg)

while(iterationOmegaDeg < 360):
    CpExact[iterationOmegaDeg] = 1- 4 * (math.sin(math.radians(omegaDeg[iterationOmegaDeg])))**2 + 4*K/(R*Vinfinity)*math.sin(math.radians(omegaDeg[iterationOmegaDeg])) - (K/(R*Vinfinity)**2);
    #CpVTAN[iterationOmegaDeg] = 1 - (Vtan[iterationOmegaDeg]/Vinfinity)**2
    iterationOmegaDeg += 1;
    
#omegaDeg = linspace(0, )


plt.plot(omegaDeg, CpExact)
plt.plot(omega*180/math.pi, Cp,"o")
plt.show()