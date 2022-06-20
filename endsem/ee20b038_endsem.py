#Gudivada Harshad Kumar
from matplotlib.pylab import *
l = 0.5; c = 2.9979e8; mu0 = 4e-7*pi;Im = 1.0; a = 0.01
N=100; lam = l*4; f = c/lam; k = 2*pi/lam; dz = l/N
#Now we define array of points from 0 to 2*N which actually indicate 0 to 2*N -1
i = 0
z_1 = zeros(2*N +1)
while i<2*N +1:
    z_1[i] = (i-N)*(dz)
    i+=1
z = z_1
z_1 = delete(z_1,N)
z_1 = z_1[1:-1]
u = z_1 #Locations of unknown curretns
#Defining a function that generates identity matrix M
def Mat(N,a):
    M = (1/(2*pi*a))*identity(2*N -2)
    return M
#Creating vectors Rz and Ru and RiN
zi, zj = meshgrid(z, z)
Rz = sqrt((zi - zj) ** 2 + ones([2 * N + 1, 2 * N + 1], dtype=complex) * a ** 2)
ui, uj = meshgrid(u, u)
Ru = sqrt((ui - uj) ** 2 + ones([2 * N - 2, 2 * N - 2], dtype=complex) * a ** 2)
RiN = delete(Rz[N], [0, N, 2 * N])
j = 1j
#Defining vectors P_B and and P 

P = (mu0 / (4 * pi)) * (exp(-k * Ru * j) / Ru) * dz
P_B = (mu0 / (4 * pi)) * (exp(-k * RiN * j) / RiN) * dz


#Defining the matrices Qij, Qb

Q = -(a / mu0) * P * ((-k * j / Ru) + (-1 / Ru ** 2))
Qb = -(a / mu0) * P_B * ((-k * j / RiN) + (-1 / RiN ** 2))

#Finding the matrix J and I eventually
J = matmul(inv(Mat(N,a)-Q),Qb)*(Im)
I = J
I = insert(I,0,0)
I = append(I,0)
I = insert(I,N,1)
#Plotting the calculated values of I 
figure(num=0,figsize=(7,7))
plot(z,abs(I),label='Iapprox')
xlabel(r'$z$')
ylabel(r'$I$')
title("Approximate value of current")
legend()
grid(True)
#Calculating and Plotting the exact values of I 
I_act = zeros(2*N + 1)
i=0
while i<2*N +1:
    if z[i]<=0:
        I_act[i] = Im*sin(k*(l+z[i]))
    elif z[i]>0:
        I_act[i] = Im*sin(k*(l-z[i]))
    i+=1
figure(num=1,figsize=(7,7))
plot(z,I_act,color='r',label="Iactual")
xlabel(r'$z$')
ylabel(r'$I$')
title("Actual current calulated from formula ")
legend()
grid(True)
show()