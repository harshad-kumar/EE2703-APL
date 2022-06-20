from matplotlib.pylab import *
#Q1: Working the given examples 
#DFT of sin(5t)
t_1 = linspace(0,2*pi,129); t_1 = t_1[:-1]
y_1 = sin(5*t_1)
Y_1 = fftshift(fft(y_1))/128
w_1 = linspace(-64,63,128)
figure(num=0,figsize=(7,7))
subplot(2,1,1)
plot(w_1,abs(Y_1))
title(r"Spectrum of $\sin(5t)$",loc='left')
ylabel(r"$|Y(\omega)|\rightarrow$")
xlabel(r"$\omega\rightarrow$")
xlim([-10,10])
grid(True)
subplot(2,1,2)
plot(w_1,angle(Y_1),'ro')
ii = where(abs(Y_1)>1e-3)
plot(w_1[ii],angle(Y_1[ii]),'go')
title(r"Phase of $\sin(5t)$",loc = 'left')
ylabel(r"$\angle Y(\omega)\rightarrow$")
xlabel(r"$\omega\rightarrow$")
xlim([-10,10])
grid(True)
show()

#DFT of (1 + 0.1cos(t))cos(10t)
t_2 = linspace(-4*pi,4*pi,513); t_2 = t_2[:-1]
y_2 = (1 + 0.1*cos(t_2))*cos(10*t_2)
Y_2 = fftshift(fft(y_2))/512
w_2 = linspace(-64,64,513); w_2 = w_2[:-1]
figure(num=1,figsize=(7,7))
subplot(2,1,1)
plot(w_2,abs(Y_2))
title(r"Spectrum of $(1 + 0.1*cos(t))*cos(10*t)$",loc='left')
ylabel(r"$|Y(\omega)|\rightarrow$")
xlabel(r"$\omega\rightarrow$")
xlim([-15,15])
grid(True)
subplot(2,1,2)
plot(w_2,angle(Y_2),'ro')
title(r"Phase of $(1 + 0.1*cos(t))*cos(10*t)$",loc='left')
ylabel(r"$\angle Y(\omega)\rightarrow$")
xlabel(r"$\omega\rightarrow$")
xlim([-15,15])
grid(True)
show()
#Q2:Generating spectrum for sin^3t and cos^3t
#DFT of sin^3t
t_3 = linspace(-4*pi,4*pi,513); t_3 = t_3[:-1]
y_3 = (3*sin(t_3) - sin(3*t_3))/4
Y_3 = fftshift(fft(y_3))/512
w_3 = linspace(-64,64,513); w_3 = w_3[:-1] 
figure(num=2,figsize=(7,7))
subplot(2,1,1)
plot(w_3,abs(Y_3))
title(r"Spectrum of $sin^{3}(t)$",loc='left')
ylabel(r"$|Y(\omega)|\rightarrow$")
xlabel(r"$\omega\rightarrow$")
xlim([-15,15])
grid(True)
subplot(2,1,2)
plot(w_3,angle(Y_3),'ro')
title(r"Phase of $sin^{3}(t)$",loc='left')
ylabel(r"$\angle Y(\omega)\rightarrow$")
xlabel(r"$\omega\rightarrow$")
xlim([-15,15])
grid(True)
show()
#DFT of cos^3t
y_4 = (3*cos(t_3) + cos(3*t_3))/4
Y_4 = fftshift(fft(y_4))/512; w_4 = w_3
figure(num=3,figsize=(7,7))
subplot(2,1,1)
plot(w_4,abs(Y_4))
title(r"Spectrum of $cos^{3}(t)$",loc='left')
ylabel(r"$|Y(\omega)|\rightarrow$")
xlabel(r"$\omega\rightarrow$")
xlim([-15,15])
grid(True)
subplot(2,1,2)
plot(w_4,angle(Y_4),'ro')
title(r"Phase of $cos^{3}(t)$",loc='left')
ylabel(r"$\angle Y(\omega)\rightarrow$")
xlabel(r"$\omega\rightarrow$")
xlim([-15,15])
grid(True)
show()
#Q3:Generate the spectrum of cos(20t +5cos(t)). Plot phase where magnitude is significant
y_5 = cos(20*t_3 + 5*cos(t_3))
Y_5 = fftshift(fft(y_5))/512; w_5 = w_4
figure(num=4,figsize=(7,7))
subplot(2,1,1)
plot(w_5,abs(Y_5))
title(r"Spectrum of $cos(20t + 5cos(t))$",loc='left')
ylabel(r"$|Y(\omega)|\rightarrow$")
xlabel(r"$\omega\rightarrow$")
xlim([-40,40])
grid(True)
subplot(2,1,2)
ii = where(abs(Y_5)>1e-3)
plot(w_5[ii],angle(Y_5[ii]),'go')
title(r"Phase of $cos(20t + 5cos(t))$",loc = 'left')
ylabel(r"$\angle Y(\omega)\rightarrow$")
xlabel(r"$\omega\rightarrow$")
xlim([-40,40])
grid(True)
show()
#Q4:find out the most accurate DFT of a Gaussian exp(-0.5t^2)
T = 2*pi
N = 128
iter = 0
tolerance = 1e-15
error = tolerance + 1
while error > tolerance:
    t = linspace(-T/2,T/2,N+1)[:-1]
    w = N/T * linspace(-pi,pi,N+1)[:-1]
    y = exp(-0.5*t**2) 
    iter = iter + 1

    Y = fftshift(fft(y))*T/(2*pi*N)
    Y_actual = (1/sqrt(2*pi))*exp(-0.5*w**2)
    error = mean(abs(abs(Y)-Y_actual))

    T = T*2
    N = N*2
print(" Error: %g \n Iteration: %g" % (error,iter))
print(" Best value of T: %g*pi \n Best value of N: %g"%(T/pi,N))
figure(num=5,figsize=(7,7))
subplot(2,1,1)
plot(w,abs(Y))
title(r"Spectrum of Guassian $exp(-0.5t^2)$",loc='left')
ylabel(r"$|Y(\omega)|\rightarrow$")
xlabel(r"$\omega\rightarrow$")
xlim([-10,10])
grid(True)
subplot(2,1,2)
plot(w,angle(Y),'ro')
title(r"Phase of $exp(-0.5t^2)$",loc='left')
ylabel(r"$\angle Y(\omega)\rightarrow$")
xlabel(r"$\omega\rightarrow$")
xlim([-10,10])
grid(True)
show()

