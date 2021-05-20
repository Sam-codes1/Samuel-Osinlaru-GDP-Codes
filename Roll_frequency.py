# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 14:06:29 2021

@author: Sam_O
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
import scipy.integrate as spi
from scipy.integrate import odeint
from scipy.signal import find_peaks

# def model(y,t):
#     theta1=theta[0]
#     theta2=theta[1]
#     dtheta1_dt= theta2
B= 93034274.1651228
w=0.36
A44=261247.614   
# I= 21024100.587*1000
# I = 902667534.498*1000
# I=594972946.799*1000
# I=4.65e10
L=300
B=50
# T=11.25
T=5.2
# GM= 14
GM=53
# vol= 106920.376
vol=36816.0
p=1025
g= 9.81
I=(20.6)**2 *p*vol + (0.5*(20.6)**2 *p*vol)
bo = p*g*vol
# bo=108000*1000
#vol = bo/(p*g)
C44=GM*vol*p*g
B4=0.00897*2*np.sqrt(C44*I)
w0= np.sqrt((bo*GM*g)/I)
# phi=(B4)/I
# # phi=0
# δ =phi/(2*w0)
# D= np.sqrt(1- δ**2)
# t= np.linspace(0,20,20)

# o= np.exp(-(phi/2) *t) *(20*np.cos(w0*D*t) +
#                             (δ*20/D)*np.sin(w0*D*t))
#plt.plot(t,o)

def f(u,t):
    return(u[1],-B4/I*u[1]-B4/I*u[0] - C44/I*u[0])
o0=[-np.radians(25), 0]
ts=np.linspace(0,100,20000)
us = odeint(f, o0,ts)
os = np.degrees(us[:,0])
# plt.plot(ts,os)
plt.xlabel("time")
plt.ylabel("Amplitude Radians")
hi=find_peaks(os)[0]
time = [ts[i] for i in hi]
T_P = [ (time[i]-time[i - 1]) for i in range(1, len(hi)) ]
t_p = np.mean(T_P)

def area(lbk):
    x = 16345 + lbk*3*0.887
    return x

def k(u,t):
    return(u[1],-0.070647*2*np.sqrt(C44*I)/I*u[1]-0.070647*2*np.sqrt(C44*I)/I*u[0] - C44/I*u[0])
o0=[-np.radians(25), 0]
ts=np.linspace(0,100,20000)
us = odeint(k, o0,ts)
os1 = np.degrees(us[:,0])

def k2(u,t):
    return(u[1],-0.09127*2*np.sqrt(C44*I)/I*u[1]-0.09127*2*np.sqrt(C44*I)/I*u[0] - C44/I*u[0])
o0=[-np.radians(25), 0]
ts=np.linspace(0,100,20000)
us = odeint(k2, o0,ts)
os2 = np.degrees(us[:,0])

def k3(u,t):
    return(u[1],-0.11191 *2*np.sqrt(C44*I)/I*u[1]-0.11191*2*np.sqrt(C44*I)/I*u[0] - C44/I*u[0])
o0=[-np.radians(25), 0]
ts=np.linspace(0,100,20000)
us = odeint(k3, o0,ts)
os3 = np.degrees(us[:,0])

def k4(u,t):
    return(u[1],-0.13256 *2*np.sqrt(C44*I)/I*u[1]-0.13256*2*np.sqrt(C44*I)/I*u[0] - C44/I*u[0])
o0=[-np.radians(25), 0]
ts=np.linspace(0,100,20000)
us = odeint(k4, o0,ts)
os4 = np.degrees(us[:,0])

hi=find_peaks(os4)[0]
time = [ts[i] for i in hi]
T_P = [ (time[i]-time[i - 1]) for i in range(1, len(hi)) ]
t_p = np.mean(T_P)
# r=2*np.pi/t_p

plt.plot(ts,os,label="No Bilge Keel")
plt.plot(ts,os1,label="60m Keel")
plt.plot(ts,os2,label="80m Keel")
plt.plot(ts,os3,label="100m Keel")
plt.plot(ts,os4,label="120m Keel")
plt.legend()
plt.xlabel("time(s)")
plt.ylabel("Amplitude (degs)")
plt.title("Roll Decay- Lightship")