#Gudivada Harshad Kumar EE20B038
import scipy.integrate as integrate
import scipy
from matplotlib.pylab import *
import numpy as np
#Defining all the required functions
def e(x):
    return exp(x)
def cos_cos(x):
    return cos(cos(x))
#Define a function that helps in converting a given function to periodic 
def periodic(fun):
    def newfun(x):
        return fun(np.remainder(x,2*np.pi))
    return newfun
#Defning the functions needed to find fourier transform
def a_f(x,k,f):
    return f(x)*np.cos(k*x)/np.pi
def b_f(x,k,f):
    return f(x)*np.sin(k*x)/np.pi    
#plotting both the graphs exp(x)(on semilog) and cos(cos(x))
n = np.arange(1,52)
figure(num=1,figsize=(7,7))
x = np.linspace(-2*np.pi,4*np.pi,1000)
title(r'Figure 1')
xlabel(r'$x\rightarrow$',size=12)
ylabel(r'$e^x\rightarrow$',size=12)
grid(True)
semilogy(x,periodic(e)(x),label = 'semilog($e^x$)')
legend()
show()
figure(num=2,figsize=(7,7))
title(r'Figure 2')
xlabel(r'$x\rightarrow$',size=12)
ylabel(r'$cos(cos(x))\rightarrow$',size=12)
grid(True)
plot(x,cos_cos(x),label='cos(cos(x))')
legend()
show()
#Finding the fourier coefficients of the given function cos(cos(x))
a_cos_cos = np.zeros(26)
b_cos_cos = np.zeros(25)
a_cos_cos[0] = integrate.quad(cos_cos,0,2*np.pi)[0]/(2*np.pi)
for i in range(1,26):
    a_cos_cos[i] = integrate.quad(a_f,0,2*np.pi,args=(i,cos_cos))[0]
for j in range(25):
    b_cos_cos[j] = integrate.quad(b_f,0,2*np.pi,args=(j+1,cos_cos))[0]
Coeff_cos_cos = np.zeros(51)
Coeff_cos_cos[0] = a_cos_cos[0]
for  i in range(1,26):
    Coeff_cos_cos[2*i-1] += a_cos_cos[i]
    Coeff_cos_cos[2*i] += b_cos_cos[i-1]
#Finding the fourier coefficients of the given function e^x
a_exp = np.zeros(26)
b_exp = np.zeros(25)
a_exp[0] = integrate.quad(e,0,2*np.pi)[0]/(2*np.pi)
for i in range(1,26):
    a_exp[i] = integrate.quad(a_f,0,2*np.pi,args=(i,e))[0]
for j in range(25):
    b_exp[j] = integrate.quad(b_f,0,2*np.pi,args=(j+1,e))[0]
Coeff_exp = np.zeros(51)
Coeff_exp[0] = a_exp[0]
for  i in range(1,26):
    Coeff_exp[2*i-1] += a_exp[i]
    Coeff_exp[2*i] += b_exp[i-1]
figure(3,figsize=(7,7))
semilogy(n,abs(Coeff_exp),'or',label='Coefficients')
title(r'Coefficients of fourier series of $e^x$ on a semilog scale')
xlabel(r'$n\rightarrow$',size=12)
ylabel(r'log(coeff)$\rightarrow$',size=12)
grid(True)
legend()
show()
figure(4,figsize=(7,7))
loglog(n,abs(Coeff_exp),'or',label='Coefficients')
title(r'Coefficients of fourier series of $e^x$ on a loglog scale')
xlabel(r'$n\rightarrow$',size=12)
ylabel(r'log(coeff)$\rightarrow$',size=12)
grid(True)
legend()
show()
figure(5,figsize=(7,7))
semilogy(n,abs(Coeff_cos_cos),'or',label='Coefficients')
title(r'Coefficients of fourier series of $cos(cos(x))$ on a semilog scale')
xlabel(r'$n\rightarrow$',size=12)
ylabel(r'log(coeff)$\rightarrow$',size=12)
grid(True)
legend()
show()
figure(6,figsize=(7,7))
loglog(n,abs(Coeff_cos_cos),'or',label='Coefficients')
title(r'Coefficients of fourier series of $cos(cos(x))$ on a loglog scale')
xlabel(r'$n\rightarrow$',size=12)
ylabel(r'log(coeff)$\rightarrow$',size=12)
grid(True)
legend()
show()
#Constructing the matrices by defining a function to use least squares approach
def matrixAb(x,f):
    b = f(x)
    A = np.zeros((x.shape[0],51))
    A[:,0] = 1
    for i in range(1,26):
        A[:,2*i-1]=np.cos(i*x)
        A[:,2*i]=np.sin(i*x)
    return A,b
x = np.linspace(0,2*np.pi,num=400,endpoint = True) 
A_coscos,b_coscos = matrixAb(x,cos_cos)
A_ex,b_ex = matrixAb(x,e) 
#Solving the problem using least squares approach 
c_coscos = scipy.linalg.lstsq(A_coscos,b_coscos)[0]
c_exp = scipy.linalg.lstsq(A_ex,b_ex)[0]
#These solutions are to be plotted with corresponding plots
figure(7,figsize=(7,7))
semilogy(n,abs(Coeff_exp),'or',label='TrueCoeff')
semilogy(n,abs(c_exp),'og',label = 'BestfitCoeff')
title(r'Coefficients of fourier series of $e^x$ on a semilog scale')
xlabel(r'$n\rightarrow$',size=12)
ylabel(r'log(coeff)$\rightarrow$',size=12)
grid(True)
legend()
show()
figure(8,figsize=(7,7))
loglog(n,abs(Coeff_exp),'or',label='TrueCoeff')
loglog(n,abs(c_exp),'og',label='BestfitCoeff')
title(r'Coefficients of fourier series of $e^x$ on a loglog scale')
xlabel(r'$n\rightarrow$',size=12)
ylabel(r'log(coeff)$\rightarrow$',size=12)
grid(True)
legend()
show()
figure(9,figsize=(7,7))
semilogy(n,abs(Coeff_cos_cos),'or',label='TrueCoeff')
semilogy(n,abs(c_coscos),'og',label='BestfitCoeff')
title(r'Coefficients of fourier series of $cos(cos(x))$ on a semilog scale')
xlabel(r'$n\rightarrow$',size=12)
ylabel(r'log(coeff)$\rightarrow$',size=12)
grid(True)
legend()
show()
figure(10,figsize=(7,7))
loglog(n,abs(Coeff_cos_cos),'or',label='TrueCoeff')
loglog(n,abs(c_coscos),'og',label='BestfitCoeff')
title(r'Coefficients of fourier series of $cos(cos(x))$ on a loglog scale')
xlabel(r'$n\rightarrow$',size=12)
ylabel(r'log(coeff)$\rightarrow$',size=12)
grid(True)
legend()
show()
#Calculating absolute error.
def error (u,v):
    p = abs(u-v)
    p1 = abs(p.max())
    return p1
print("The error in Coefficients of e^x =",error(Coeff_exp,c_exp))
print("The error in Coefficients of cos(cos(x)) =",error(Coeff_cos_cos,c_coscos))
#Calculating Ac for both functions
y_cos_cos = np.dot(A_coscos,c_coscos)
y_exp = np.dot(A_ex,c_exp)
x = np.linspace(-2*np.pi,4*np.pi,1000)
x_1 = np.linspace(0,2*np.pi,num=400)
figure(11,figsize=(7,7))
semilogy(x,periodic(e)(x),label = 'Actual')
semilogy(x_1,y_exp,'go',label = 'Bestfit')
title(r'Figure 1 modified')
xlabel(r'$x\rightarrow$',size=12)
ylabel(r'$e^x\rightarrow$',size=15)
grid(True)
legend()
show()
figure(12,figsize=(7,7))
plot(x,cos_cos(x),label='Actual')
plot(x_1,y_cos_cos,'go','bestfit')
title(r'Figure 2')
xlabel(r'$x\rightarrow$',size=12)
ylabel(r'$cos(cos(x))\rightarrow$',size=12)
grid(True)
legend()
show()







