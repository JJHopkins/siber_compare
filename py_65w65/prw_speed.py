#!/usr/bin/env python
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as pl

# Input dielectric response data
eiz_x = np.loadtxt('data/eiz_x_output_eV.txt') # len(n), LDS, ~material response in perpendicular direction
eiz_z = np.loadtxt('data/eiz_z_output_eV.txt') # len(n), LDS, ~material response in parallel direction
eiz_w = np.loadtxt('data/eiz_w_output_eV.txt') # len(n), LDS, ~response of water, intervening medium

# Constants
c = 2.99e8   # light speed in m/s
coeff = 2.411e14 # in rad/s
Temp = 297 # K 
kbT = Temp * 1.3807e-23 # in J

# Matsubara frequencies
ns = np.arange(0.,500.)  # index for Matsubara sum
zs = ns * coeff           # thermal frequencies, they are multiples of n
Ls = np.arange(1e-9,1e-6,1e-9) #np.linspace(1.0e-9, 1.0e-6, 20) # separation distance between 2 cyclinders

def Aiz(perp, par,med):
	return (2.0*(perp-med)*med)/((perp+med)*(par-med))

def Delta(par,med):
	return (par - med)/med

def Pn(e,zn,l,v):
	return np.sqrt(e)*zn*l*(1./v)

a =  Aiz(eiz_x,eiz_z,eiz_w)
delta = Delta(eiz_z,eiz_w)

p = np.zeros(shape = (len(Ls),len(ns)))

for i,L in enumerate(Ls):
	for j,n in enumerate(ns):
		p[i,j] = Pn(eiz_w[j],zs[j],L,c)

t = np.linspace(0.,100000.,100)
y = np.linspace(1.000001,100001.000001,100)
T,Y = np.meshgrid(t,y)
A = np.zeros(shape = (len(Ls),len(ns)))

for i,L in enumerate(Ls):
	for j,n in enumerate(ns):
		p[i,j] = Pn(eiz_w[j],zs[j],L,c)
		f = (1./np.sqrt(Y**2 -1.0))*np.exp(-2.*Y*p[i,j]*np.sqrt(T*T+1))*(T/np.sqrt(T*T+1))*(2.*(1.-3.*a[j])*(1.-3.*a[j])*T**4 + 4.*(1.+2.*a[j]+2.*a[j]+3*a[j]*a[j])*T**2 + 4.*(1.+a[j])*(1.+a[j]) +  (T**4 + 4.*T**2+4.)*(1.-a[j])*(1.-a[j]) )
        print i,j,f	
        Ft = np.sum(f,axis = 1)
        Fty = np.sum(Ft)
        #	print j,Fty
        A[i,j] = delta[j]*delta[j]*(p[i,j]**5)*Fty
        #print A
        sum_A = (kbT/(12.*np.pi)) * np.sum(A, axis = 1)
np.savetxt('140316_65w65_prw_sum_A.txt',sum_A)
np.savetxt('140316_65w65_prw_Ls.txt',Ls)

pl.figure()
pl.loglog(Ls,sum_A)
pl.show()
