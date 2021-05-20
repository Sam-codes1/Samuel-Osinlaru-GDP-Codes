# -*- coding: utf-8 -*-
"""
Created on Sat Apr 17 18:31:33 2021

@author: Sam_O
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
import scipy.integrate as spi
from scipy.integrate import odeint
from scipy.signal import find_peaks
#inputs
m=20000
# # b=0.05
g=9.81
L=34.6
h=45
a= 0.5
p = 450
k= np.pi/L
B=45
w=0.5
v=4.36e-6
o_c= np.radians(12)
E = np.pi/4 *p*g*B* a**2/k
E1 = (g**2 * a**2)/w**2 * ((np.pi**2)/4) * np.sqrt((v*p*np.sinh(2*k*h))/2*w*np.cosh(k*h)**2)
E2 = (g**2 * a**2)/w**2 * ((np.pi)/2) * np.sqrt(v*p/2*w) * (B*k/np.cosh(k*h)**2) * (np.sinh(2*k*h)/2 - k*h)
E3 = (g**2 * a**2)/w**2 * ((np.pi**2)/2) * np.sqrt(v*p/2*w) * (B*k/np.cosh(k*h)**2) 
row = (E1 +E2 +E3)/(2*E)

def pen(theta,t):
    a= 1
    b= 15
    c= 25e-4
    d=12
    L=0.8
    # row=0.000273
    # m2= m*(8*L*np.tanh(np.pi*L/L))/np.pi**3 *h
    m1=1
    theta1=theta[0]
    theta2=theta[1]
    dtheta1_dt= theta2
    R0 = math.sin(theta1)
    At =-0.5* w**2 * np.sin(w*t)
    dtheta2_dt=  -row*theta1 + ((1/L)*(-g*R0 + At)) - ((a*(theta1/o_c)**b) + (c*(theta1/o_c)**(2*d)) * theta2**m1) 
    dtheta_dt = [dtheta1_dt,dtheta2_dt]    
    return dtheta_dt

# initial conditions
theta_0 = [0,o_c]

#time points
t= np.linspace(0,4000,4000)

# soliving ode
theta_results = odeint(pen, theta_0, t)
Angular_Displacement = []
Angular_velocity = []
for i in range (0,len(theta_results)):
    angular_displacement = theta_results[i,0]
    Angular_Displacement.append(angular_displacement)
    
    angular_velocity = theta_results[i,1]
    Angular_velocity.append(angular_velocity)

# " Impact Force"
# theta = np.linspace(0,45,46)
# F= 0.5* (theta/45)**19
# plt.plot(theta,F)
T=2*np.pi/w
plt.xlim(0,100)
plt.plot(t,Angular_velocity)
plt.xlabel("time(s)")
plt.ylabel("Angular Velocity (rads/s)")
# i=np.interp(25, t, Angular_velocity)

hi=find_peaks(Angular_velocity)[0]
time = [t[i] for i in hi]
T_P = [ (time[i]-time[i - 1]) for i in range(1, len(hi)) ]
t_p = np.mean(T_P)
