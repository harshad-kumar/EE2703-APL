#Gudivada Harshad Kumar EE20B038
from matplotlib.pylab import *
import scipy.special as sp
import numpy as np
#defining the function
def g(t,A,B):
	return A*sp.jn(2,t) + B*t
   
#claculating matrix M to check equality for question6
    
#loading fitting.dat and extracting data and time    
try:
    data = np.loadtxt("fitting.dat",dtype=float)
    t = np.array(data[:,0])
    d = np.array(data[:,1:])
except:
    print("File not found in the directory")
    exit(0)
sigma = logspace(-1,-3,9)
#plotting out figure0 as asked in q4
figure(num=0,figsize=(7,7))
for i in range(9):
	plot(t,d[:,i],label = r'$\sigma$=%.3f'%sigma[i])    
y0 = g(t,1.05,-0.105)
plot(t,y0,color = 'brown',label ="truevalue")
title(r'Q4:Data to be fitted the theory',size=16)
xlabel(r'$t \rightarrow$',fontsize =12)
ylabel(r'$f(t)+noise \rightarrow$',fontsize = 12)
grid(True)
legend()
show()
#error plot in first set of data along with actual plot
figure(num=1,figsize=(7,7))
plot(t,y0,color = 'black',label="f(t)")
errorbar(t[::5],d[::5,0],sigma[0],fmt='ro',label ="Errorbar")
title(r'Q5:Data points for $\sigma=0.1$ along with exact function',size=16)
xlabel(r'$t \rightarrow$',fontsize=12)
grid(True)
legend()
show()
#As mentioned in Q6 create matrix M
J2 = sp.jv(2,t) 
M = c_[J2,t]
P = np.array([1.05,-0.105])
y = np.dot(M,P)
equi = y-y0
if all(equi<=(1e-15)):                  #due to some precision errors it may happen that few values go very close to zero but not zero
    print("Both the solutions are equal")
else:
    print("Both the solutions are not equal")
#Computing the mean squared error    
A = np.linspace(0,2,21)
B = np.linspace(-0.2,0,21)
epsilon = np.zeros((21,21),float)
for i in range(21):
    for j in range(21):
        epsilon[i][j] = np.mean(np.square(y0 - g(t,A[i],B[j])))
# to plot a contour of ε[i][j] 
figure(num=2,figsize=(7,7))
Contour = contour(A,B,epsilon[:,:],levels=np.linspace(0,0.3,13))
clabel(Contour, inline = True, fontsize = 6)       
title(r'Q8:contour plot of $\epsilon_{ij}$',size=16)
xlabel(r'A$\rightarrow$',fontsize=12)
ylabel(r'B$\rightarrow$',fontsize=12)
plot(1.05,-0.105,'ro')
annotate('Exact Location',xy=(1.05,-0.105))
show()
#Calculating the error in estimate of A and B and plotting the variation of error with noice
Ex = np.empty((9,1))
Ey = np.empty((9,1))
for i in range(9):
	AB = linalg.lstsq(M,data[:,i+1],rcond=None)
	Ex[i] = abs(AB[0][0]-P[0])
	Ey[i] = abs(AB[0][1]-P[1])      
figure(num=3,figsize=(7,7))
grid(True)
title(r'Q10:Variation of error with noise',size =16)
xlabel(r'Noise Standard Deviation $\rightarrow$',fontsize=12)
ylabel(r'MS Error$\rightarrow$',fontsize=12)
plot(sigma, Ex, linestyle = 'dashed', marker = 'o', color = 'r', label = 'Aerr')
plot(sigma, Ey, linestyle = 'dashed', marker = 'o', color = 'g', label = 'Berr')
legend()
show()
#Reploting the above curve using loglog
figure(num=4,figsize=(7,7))
grid(True)
title(r'Q11:Variation of error with noise',size=16)
xlabel(r'$\sigma_{n}\rightarrow$')
ylabel(r'$MS Error$\rightarrow$')
loglog(sigma, Ex,'ro', label = 'Aerr')
loglog(sigma, Ey,'go', label = 'Berr')
stem(sigma,Ex,'-ro')
stem(sigma,Ey,'-go')
legend()
show()




      
      

