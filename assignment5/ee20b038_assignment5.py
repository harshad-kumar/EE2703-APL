#EE20B038 GUDIVADA HARSHAD KUMAR
from matplotlib.pylab import *
import sys
import mpl_toolkits.mplot3d.axes3d as p3
#Taking inputs from users if provided
if len(sys.argv) == 5:
    Nx = int(sys.argv[1])
    Ny = int(sys.argv[2])
    radius = int(sys.argv[3])
    Niter  = int(sys.argv[4])
#Setting up the required parameters when no inputs are given.
else:
    Nx = 25
    Ny = 25
    radius = 8
    Niter = 1500
#Creating matrice phi and initiaising it 
phi = zeros((Ny,Nx))
x = linspace(-0.5,0.5,Nx)
y = linspace(-0.5,0.5,Ny)
Y,X = meshgrid(y,x)
R = (X**2) + (Y**2)
ii = where(R <= (0.35*0.35)) #potential in this area is 1 hence
phi[ii] = 1.0
n = arange(Niter)
#Contour plot of potential with label
figure(num =0,figsize = (7,7))
plot(ii[0]/Nx-0.48,ii[1]/Ny-0.48,'ro',label="V = 1")
contourf(X,Y,phi)
colorbar()
title("Initial Potential Contour")
xlim(-0.5,0.5)
ylim(-0.5,0.5)
xlabel(r'$X\rightarrow$')
ylabel(r'$Y\rightarrow$')
grid(True)
legend()
errors = zeros(Niter)
#Updating the potential 
for i in range(Niter):
    phiold = phi.copy()
    phi[1:-1,1:-1] = 0.25*(phiold[1:-1,0:-2] + phiold[1:-1,2:] + phiold[0:-2,1:-1] + phiold[2:,1:-1])
#Set appropriate boundary conditions
    phi[:,0]=phi[:,1] 
    phi[:,Nx-1]=phi[:,Nx-2] 
    phi[0,:]=phi[1,:] 
    phi[Ny-1,:]=0
    phi[ii] = 1
    errors[i]=(abs(phi-phiold)).max()
#Current densities
Jx,Jy = (1/2*(phi[1:-1,0:-2]-phi[1:-1,2:]),1/2*(phi[:-2,1:-1]-phi[2:,1:-1]))
log_y_mat = log(transpose(errors))       #modelling the error in a exponential function
x_mat = c_[arange(0,Niter),ones(shape = [Niter,1])]
E1 = lstsq(x_mat,log_y_mat,rcond=None)[0]     #lstsq fit for whole error array
E2 = lstsq(x_mat[500:-1],log_y_mat[500:-1],rcond=None)[0]     #lstsq fit for error values post 500 iterations
#Plotting the error values vs iteration on semilog plot
figure(num=1,figsize=(7,7))
semilogy(n,errors,label = "Actual")
semilogy(n[::50],errors[::50],'ro',label = "Every 50th point")
title("Error versus iteration in a semilog plot")
xlabel(r'$Iteration\rightarrow$',size=12)
ylabel(r'$Error\rightarrow$',size=12)
legend()
grid(True)
#Plotting the error values vs iteration on loglog plot
figure(num=2,figsize=(7,7))
loglog(n,errors,label="Actual")
loglog(n[::50],errors[::50],'ro',label = "Every 50th point")
title("Error versus iteration in a loglog plot")
xlabel(r'$Iteration\rightarrow$',size=12)
ylabel(r'$Error\rightarrow$',size=12)
legend()
grid(True)
#Plot of errors calculated by least squares method
figure(num=3,figsize=(7,7))
semilogy(arange(0,Niter,50),exp(dot(x_mat,E1))[::50],'bo',markersize=8,label='Fit 1')
semilogy(n,exp(dot(x_mat,E1)),label = "Best fit all errors")
title("Best fit plot of errors")
xlabel(r'$Iteration\rightarrow$',size=12)
ylabel(r'$Error\rightarrow$',size=12)
legend()
grid(True)
#Plot of errors post 500 iterations found by least squares method
figure(num=4,figsize=(7,7))
semilogy(arange(500,Niter,50),exp(dot(x_mat[500:-1],E2))[::50],'ro',markersize=8,label='Fit 2')
semilogy(n[500:-1],exp(dot(x_mat[500:-1],E2)),label = "Best fit post 500 errors")
xlabel(r'$Iteration\rightarrow$',size=12)
ylabel(r'$Error\rightarrow$',size=12)
legend()
grid(True)
#Plot of fit1 fit2 and errors in same graph
figure(num=5,figsize=(7,7))
loglog(arange(0,Niter,50),exp(dot(x_mat,E1))[::50],'bo',markersize=8,label='Fit 1')
loglog(arange(500,Niter,50),exp(dot(x_mat[500:-1],E2))[::50],'ro',markersize=8,label='Fit 2')
loglog(n,errors,label =errors)
xlabel(r'$Iteration\rightarrow$',size=12)
ylabel(r'$Error\rightarrow$',size=12)
legend()
grid(True)
#Contour plot of potential in 2-D
figure(num=6,figsize=(7,7))
contourf(Y,X[::-1],phi)
plot(ii[0]/Nx-0.48,ii[1]/Ny-0.48,'ro',label="V = 1")
title("Contour plot of potential in 2-D")
xlabel(r'$X\rightarrow$')
ylabel(r'$Y\rightarrow$')
colorbar()
grid(True)
legend()
#3-D plot of the potential
fig1=figure(num=7,figsize=(7,7))
ax=p3.Axes3D(fig1) 
title("The 3-D surface plot of the potential")
xlabel(r'$X\rightarrow$')
ylabel(r'$Y\rightarrow$')
surf = ax.plot_surface(Y, X[::-1], phi, rstride=1, cstride=1, cmap=cm.jet)
fig1.colorbar(surf)
#plotting current densities 
figure(num=8,figsize=(7,7))
title("Vector plot of current flow")
quiver(Y[1:-1,1:-1],-X[1:-1,1:-1],-Jx[:,::-1],-Jy)
plot(ii[0]/Nx-0.48,ii[1]/Ny-0.48,'ro')
title("Vector plot of the current flow")
xlabel(r'$X\rightarrow$')
ylabel(r'$Y\rightarrow$')

show()



