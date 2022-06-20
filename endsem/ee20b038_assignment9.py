from matplotlib.pylab import *
from mpl_toolkits.mplot3d import Axes3D
#Question 1 plotting spectrum of sin(sqrt(2)t) without windowing 
t = linspace(-pi,pi,65); t = t[:-1]
dt = t[1]-t[0]; fmax = 1/dt
y = sin(sqrt(2)*t)
y[0] = 0 
y = fftshift(y) 
Y = fftshift(fft(y))/64.0
w = linspace(-pi*fmax,pi*fmax,65); w = w[:-1]
figure(num=0,figsize=(7,7))
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-10,10])
ylabel(r"$|Y|\rightarrow$")
title(r"Spectrum of $\sin\left(\sqrt{2}t\right)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-10,10])
ylabel(r"Phase of $Y\rightarrow$",size=16)
xlabel(r"$\omega\rightarrow$",size=16)
grid(True)
#Plotting spectrum of sin(sqrt(2)t) after windowing 
t = linspace(-4*pi,4*pi,257); t = t[:-1]
dt = t[1]-t[0]; fmax = 1/dt
n = arange(256)
wnd = fftshift(0.54+0.46*cos(2*pi*n/256))
y = sin(sqrt(2)*t)*wnd
y[0] = 0 
y = fftshift(y) 
Y = fftshift(fft(y))/256.0
w = linspace(-pi*fmax,pi*fmax,257); w = w[:-1]
figure(num=1,figsize = (7,7))
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-4,4])
ylabel(r"$|Y|\rightarrow$")
title(r"Improved Spectrum of $\sin\left(\sqrt{2}t\right)$")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-4,4])
ylabel(r"Phase of $Y\rightarrow$")
xlabel(r"$\omega\rightarrow$")
grid(True)
#Question 2 plotting spectrum of cos^3(wt) with and without windowing 
y = cos(0.86*t)**3
y1 = y*wnd
y[0]=0
y1[0]=0
y = fftshift(y)
y1 = fftshift(y1)
Y = fftshift(fft(y))/256.0
Y1 = fftshift(fft(y1))/256.0
figure(num=2,figsize=(7,7))
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-4,4])
ylabel(r"$|Y|\rightarrow$")
title(r"Spectrum of $\cos^{3}(0.86t)$ without Hamming window")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-4,4])
ylabel(r"Phase of $Y\rightarrow$")
xlabel(r"$\omega\rightarrow$")
grid(True)

figure(num=3,figsize=(7,7))
subplot(2,1,1)
plot(w,abs(Y1),lw=2)
xlim([-4,4])
ylabel(r"$|Y|\rightarrow$")
title(r"Spectrum of $\cos^{3}(0.86t)$ with Hamming window")
grid(True)
subplot(2,1,2)
plot(w,angle(Y1),'ro',lw=2)
xlim([-4,4])
ylabel(r"Phase of $Y\rightarrow$")
xlabel(r"$\omega\rightarrow$")
grid(True)
#Question 3, considering delta = 0.5 and w0 = 1.5
w0 = 1.5
delta = 0.5
t = linspace(-pi,pi,129)[:-1]
dt = t[1]-t[0]; fmax = 1/dt
n = arange(128)
wnd = fftshift(0.54+0.46*cos(2*pi*n/128))
y = cos(w0*t + delta)*wnd
y[0]=0
y = fftshift(y)
Y = fftshift(fft(y))/128.0
w = linspace(-pi*fmax,pi*fmax,129); w = w[:-1]
figure(num=4,figsize=(7,7))
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-4,4])
ylabel(r"$|Y|\rightarrow$")
title(r"Spectrum of $\cos(w_0t+\delta)$ with Hamming window")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-4,4])
ylabel(r"Phase of $Y\rightarrow$")
xlabel(r"$\omega\rightarrow$")
grid(True)
# w0 is calculated by finding the weighted average where w>0
def est_omega(w,Y):
    ii1 = where(w>0)
    omega = (sum(abs(Y[ii1])**2*w[ii1])/sum(abs(Y[ii1])**2))#weighted average
    return omega

def est_delta(w,Y,sup = 1e-4,window = 1):
    ii2=np.where(np.logical_and(np.abs(Y)>sup, w>0))[0]
    np.sort(ii2)
    points=ii2[1:window+1]
    return (np.sum(np.angle(Y[points]))/len(points)) #weighted average for first 2 points
print ("calculated omega =", est_omega(w,Y))
print ("calculated delta = ",est_delta(w,Y))
#Question4 adding “white gaussian noise” to the signal in previous question 
y = (cos(w0*t + delta) + 0.1*randn(128))*wnd
y[0]=0
y = fftshift(y)
Y = fftshift(fft(y))/128.0
figure(num=5,figsize=(7,7))
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-4,4])
ylabel(r"$|Y|\rightarrow$")
title(r"Spectrum of a noisy $\cos(w_0t+\delta)$ with Hamming window")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-4,4])
ylabel(r"Phase of $Y\rightarrow$")
xlabel(r"$\omega\rightarrow$")
grid(True)
print ("calculated omega with noise =", est_omega(w,Y))
print ("calculated delta with noise = ",est_delta(w,Y))
#Question5 plotting spectrum of a chirped signal 
t = linspace(-pi,pi,1025); t = t[:-1]
dt = t[1]-t[0]; fmax = 1/dt
n = arange(1024)
wnd = fftshift(0.54+0.46*cos(2*pi*n/1024))
y = cos(16*t*(1.5 + t/(2*pi)))*wnd
y[0]=0
y = fftshift(y)
Y = fftshift(fft(y))/1024.0
w = linspace(-pi*fmax,pi*fmax,1025); w = w[:-1]
figure(num=6,figsize=(7,7))
subplot(2,1,1)
plot(w,abs(Y),lw=2)
xlim([-100,100])
ylabel(r"$|Y|\rightarrow$")
title(r"Spectrum of chirped function")
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro',lw=2)
xlim([-100,100])
ylabel(r"Phase of $Y\rightarrow$")
xlabel(r"$\omega\rightarrow$")
grid(True)
#Question6 plotting "time frequency" plot for the chirped signal
t=linspace(-pi,pi,1025);t=t[:-1]
t_arrays=split(t,16)

Y_mags=zeros((16,64))
Y_angles=zeros((16,64))
for i in range(len(t_arrays)):
	n = arange(64)
	wnd = fftshift(0.54+0.46*cos(2*pi*n/64))
	y = cos(16*t_arrays[i]*(1.5 + t_arrays[i]/(2*pi)))*wnd
	y[0]=0
	y = fftshift(y)
	Y = fftshift(fft(y))/64.0
	Y_mags[i] = abs(Y)
	Y_angles[i] = angle(Y)
t = t[::64]	
w = linspace(-fmax*pi,fmax*pi,64+1); w = w[:-1]
t,w = meshgrid(t,w)

fig = figure(7)
ax = fig.add_subplot(111, projection='3d')
surf=ax.plot_surface(w,t,Y_mags.T,cmap='viridis',linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5, aspect=5)
ax.set_title('surface plot');
ylabel(r"$\omega\rightarrow$")
xlabel(r"$t\rightarrow$")

show()