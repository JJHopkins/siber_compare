#$ {\bf Free energy between two skewed cylinders (CG-10 in water). Full retarded result, function of separation $\ell$ and angle $\theta$} \\
#$ Equation 12: $G(\ell,\theta) = - \frac{ (\pi R_1^{2})(\pi R_2^{2}) }{2 \pi~\ell^{4} \sin{\theta}} \left( {\cal A}^{(0)}(\ell) + {\cal A}^{(2)}(\ell) \cos 2\theta \right)$ \\
#$ $G(\ell,\theta) = - \frac{k_BT}{64 \pi} \frac{ \pi^2 R_1^{2} R_2^{2} }{\ell^{4} \sin{\theta}} {\sum_{n=0}^{\infty}}' \Delta_{1,\parallel} \Delta_{2,\parallel} ~p_n^{4} ~\int_0^{\infty} t dt ~\frac{e^{- 2 p_n \sqrt{t^{2} + 1}}}{(t^{2} + 1)} \tilde g(t, a_1(i \omega_n), a_2(i \omega_n), \theta),$ \\
#$ with $\tilde g(t, a_1, a_2, \theta) &=& 2 \left[ (1+3a_1)(1+3a_2) t^{4} + 2 (1+2a_1+2a_2+3a_1a_2) t^{2}  + 2(1+a_1)(1+a_2)\right] + \nonumber \\
#$ & & ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ + (1-a_1)(1-a_2)(t^{2} + 2)^2 \cos 2\theta.$ \\
#!/usr/bin/python
import numpy as np
import scipy.optimize as opt
from scipy.integrate import trapz
import matplotlib.pyplot as pl
from matplotlib import axis as ax
# use pyreport -l file.py
from pylab import show
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.mplot3d import Axes3D
from pylab import pause
from matplotlib.backends.backend_pdf import PdfPages

# Input dielectric response data
eiz_x = np.loadtxt('data/eiz_x_output_eV.txt') # len(n), LDS, ~material response in perpendicular direction
eiz_z = np.loadtxt('data/eiz_z_output_eV.txt') # len(n), LDS, ~material response in parallel direction
eiz_w = np.loadtxt('data/eiz_w_output_eV.txt') # len(n), LDS, ~response of water, intervening medium

x90_A0, y90_A0 = np.loadtxt('data/65-W-65-PERPA0.PRN',unpack=True, usecols = [0,1]) # water as intervening medium
x90_A2, y90_A2 = np.loadtxt('data/65-W-65-PERPA2.PRN',unpack=True, usecols = [0,1]) # water as intervening medium

# Constants
r_1 = 1.0e-9
r_2 = 1.0e-9
c = 2.99e8   # light speed in m/s
coeff = 2.411e14 # in rad/s
Temp = 297 # K 
kbT = Temp * 1.3807e-23 # in J
# Matsubara frequencies
ns = np.arange(0.,500.)  # index for Matsubara sum
zs = ns * coeff           # thermal frequencies, they are multiples of n
Ls = np.arange(1e-9,1e-6,1e-9) #np.linspace(1.0e-9, 1.0e-6, 20) # separation distance between 2 cyclinders
thetas = [np.pi/8,np.pi/4,np.pi/3,np.pi/2]
#ts = np.linspace(0.,100000.,1000)

dt = 10.0
ts = np.arange(0.0,100000.,dt)

def Aiz(perp, par,med):
	return (2.0*(perp-med)*med)/((perp+med)*(par-med))
def ys(a,time,eizw,L, N):
	term0 = ( time    / (time*time+1.0)               )
	term1 =	( time**4 * 2.0*(1. + 3.*a)*(1.+3.*a)     )
	term2 = ( time**2 * 4.0*(1. + 2.0*a+2.0*a+3.0*a*a))
	term3 = (          4.0*(1. + a)*(1.0 + a)         )
	term4 = (-2.0 * np.sqrt(eizw)* L * coeff * N / c * np.sqrt(time*time + 1.0))
	#print 'ys term0', term0
	#print 'ys term1', term1
	#print 'ys term2', term2
	#print 'ys term3', term3
	#print 'ys term4', term4
	#print '----'
	return (term0) * np.exp(term4)*( (term1) + (term2) + (term3))#* term5
def y_2s(a,time,eizw, L, N): 
	term0 = (time    / (time*time+1.0)                               ) 
	term1 = ((1.- a)*(1.- a)*(time * time  + 2.0)*(time * time + 2.0))                     
	term2 = (-2.0 * np.sqrt(eizw)* L * coeff * N / c * np.sqrt(time*time + 1.0))
	#print 'y_2s term0', term0
	#print 'y_2s term1', term1
	#print 'y_2s term2', term2
	#print '----'
	return term0 * term1* np.exp(term2) #* term3
def As(eizz,eizw,L,N,Y): 
	term1 = (((eizz-eizw)/eizw)*((eizz-eizw)/eizw))
	term2 = (Y * eizw *eizw * (coeff*N)**4 * L**4 / (c**4))
	#term3 = Y
	#print 'As term1 = ', term1
	#print 'As term2 = ', term2
	##print 'As term3 = ', term3
	#print '----'
	return  term1 *  term2# * term3
def A_2s(eizz,eizw, L , N ,Y):		
	term1 = (((eizz-eizw)/eizw)*((eizz-eizw)/eizw))
	term2 = (Y * eizw *eizw * (coeff*N)**4 * L**4 / (c**4))
	#term3 = Y
	#print 'A_2s term1 = ', term1
	#print 'A_2s term2 = ', term2
	##print 'A_2s term3 = ', term3
	#print '----'
	return (term1 * term2)# * term3


y = np.zeros(shape=(len(ns),len(Ls)))
y_2 = np.zeros(shape=(len(ns),len(Ls)))
A   = np.zeros(shape=(len(ns),len(Ls)))
A_2 = np.zeros(shape=(len(ns),len(Ls)))
EL = np.zeros(len(Ls))
G_l_t_dt = np.zeros(shape=(len(Ls),len(thetas)))
A2_theta = np.zeros(shape=(len(Ls),len(thetas)))

aiz = []
aiz = Aiz(eiz_x,eiz_z, eiz_w) # of length = len(ns)

for k,length in enumerate(Ls):
	sum_A = np.empty(len(Ls))
	sum_A_2 = np.empty(len(Ls))
	for j,n in enumerate(ns):
		# Integral:
		y[j,k]   = trapz(ys(aiz[j],ts,eiz_w[j],length,n),ts,dt)
		y_2[j,k] = trapz(y_2s(aiz[j],ts,eiz_w[j],length,n),ts,dt)
		#print 'dt Integral   y = ',i,k,j, y
		#print 'dt Integral y_2 = ',i,k,j, y_2
		#print '----'
		#print 'N terms for A0 = '  , As(eiz_z[j],eiz_w[j],length,n,y)
		#print 'N terms for A2 = ', A_2s(eiz_z[j],eiz_w[j],length,n,y_2)
		#print '----'
		A[j,k]   = As(eiz_z[j],eiz_w[j],length,n,y[j,k])
		A_2[j,k] = A_2s(eiz_z[j],eiz_w[j],length,n,y_2[j,k])# * np.cos(2.0*theta)  		
		A[0] = (1./2)*A[0]  		
		A_2[0] = (1./2)*A_2[0]  		
	sum_A = np.sum(A,axis=0)
	#print 'sum of A0 = ', k,j,sum_A
	sum_A_2 = np.sum(A_2,axis=0)
	#print 'sum of A2 = ', k,j,sum_A_2
	#print '----'
	#print 'shape sum_A_2 = ', np.shape(sum_A_2)
	#sys.exit()
for k,length in enumerate(Ls):
	for i, theta in enumerate(thetas):
		EL[k] = 1./(length*length*length*length)
		A2_theta[k,i] = sum_A_2[k]* np.cos(2.0*theta)
		#G_l_t_dt[k,i] = (1.602e-19 / 4.11e-21) * (1./32) * EL[k]*np.pi*r_1*r_1*r_2*r_2*(sum_A[k] + sum_A_2[k]* np.cos(2.0*theta) )/(2.0*np.sin(theta))# (1e21)*
		G_l_t_dt[k,i] = (1./32) * EL[k]*np.pi*r_1*r_1*r_2*r_2*(sum_A[k] + sum_A_2[k]* np.cos(2.0*theta) )/(2.0*np.sin(theta))# (1e21)*
np.savetxt('065_A0.txt',sum_A)
np.savetxt('065_A2.txt',sum_A_2)
#np.savetxt('G_srw.txt',G_l_t_dt)
np.savetxt('Ls.txt',Ls)

pl.figure()
pl.loglog(1e9*Ls,(1e21*kbT/32)*sum_A  ,'b-', label = r'$\mathcal{A^{(0)}}(\ell)$')
pl.loglog(1e9*Ls,(1e21*kbT/32)*sum_A_2,'g-', label = r'$\mathcal{A^{(2)}}(\ell)$')

#pl.loglog(x90_A0, y90_A0,'k-' , label = r'GH $\mathcal{A^{(0)}}(\ell)$')
#pl.loglog(x90_A2, y90_A2,'k--', label = r'GH $\mathcal{A^{(2)}}(\ell)$')

pl.xlabel(r'$\mathrm{separation}\,\ell\,\,\,\rm{[nm]}$', size = 20)
pl.ylabel(r'$\mathrm{\mathcal{A^{(0)},\,\,A^{(2)}}}\,\,\,\rm{[zJ]}$', size = 20)
pl.title(r'$\mathrm{Hamaker \, coeff.s \,:\,skewed,\,retarded,\,water}$', size = 20)
#pl.legend(loc = 'best')
#pl.axis([1e-9,1e-6,1e-24,1e-19])
pl.minorticks_on()
pl.ticklabel_format(axis = 'both')
pl.grid(which = 'both')
pl.tick_params(which = 'both',labelright = True)
pl.savefig('plots/skew_ret_water/140309_65w65_GH_skew_ret_A0_A2.pdf')
show()

#ls4 = 1e9*Ls[ 2]#2]
#ls5 = 1e9*Ls[12]#4]
#ls6 = 1e9*Ls[22]#6]
#ls1 = 1e9*Ls[32]#8]
#ls2 = 1e9*Ls[42]#12]
#ls3 = 1e9*Ls[52]#16]

#fig = pl.figure()
#ax = fig.add_axes([0.1,0.1,0.8,0.8])
##pl.semilogy(thetas, G_l_t_dt)
#ax.semilogy(thetas, G_l_t_dt[ 2,:], label = r'$\ell$ = %1.2f nm' %ls4)
#ax.semilogy(thetas, G_l_t_dt[12,:], label = r'$\ell$ = %1.2f nm' %ls5)
##ax.semilogy(thetas, G_l_t_dt[22,:], label = r'$\ell$ = %1.2f nm' %ls6)
#ax.semilogy(thetas, G_l_t_dt[32,:], label = r'$\ell$ = %1.2f nm' %ls1)
##ax.semilogy(thetas, G_l_t_dt[42,:], label = r'$\ell$ = %1.2f nm' %ls2)
#ax.semilogy(thetas, G_l_t_dt[52,:], label = r'$\ell$ = %1.2f nm' %ls3)
##ax.semilogy(0,0,'', label = r'$G_\theta = cos(2\theta)/2sin(\theta)$')
#pl.xlabel(r'$Angle\,\,\mathrm{[radians]}$', size = 20)
#pl.ylabel(r'$-G(\ell,\theta)\,\,\mathrm{[k_{B}T]}$', size = 20)
##pl.axis([0,1.7,1e-10,1.0])
##pl.title(r'$\mathrm{-G(\ell,\theta)\,vs.\,angle:\,skewed,\,retarded,\,water}$', size = 20)
#pl.legend(loc = 'lower left')
##pl.savefig('plots/skew_ret_water/skew_ret_water_G_vs_theta.pdf')
##show()
#pl.minorticks_on()
#pl.ticklabel_format(axis = 'both')
#pl.grid(which = 'both')
#pl.tick_params(which = 'both',labelright = True)
#pl.savefig('plots/skew_ret_water/140306_290w290_skew_ret_G_vs_theta_fixed_l.pdf')
#show()

#pl.figure()
##pl.loglog(Ls, G_l_t_dt)#, label = labels[i])
#pl.loglog(Ls, G_l_t_dt[:,1], label = r'$\theta = \pi/4$')
#pl.loglog(Ls, G_l_t_dt[:,2], label = r'$\theta = \pi/3$')
#pl.loglog(Ls, G_l_t_dt[:,3], label = r'$\theta = \pi/2$')
#pl.xlabel(r'$\ell\,\,\mathrm{[m]}$', size = 24)
#pl.ylabel(r'$-G(\ell,\theta)\,\,\mathrm{[k_{B}T]}$', size = 20)
##pl.axis([1.0e-9, 1.0e-6,1e-16,1e3])
##pl.title(r'$\mathrm{-G(\ell,\theta)\,vs.\,separation:\,skewed,\,retarded,\,water}$', size = 20)
#pl.legend(loc = 'best')
#pl.minorticks_on()
#pl.ticklabel_format(axis = 'both')
#pl.grid(which = 'both')
#pl.tick_params(which = 'both',labelright = True)
#pl.savefig('plots/skew_ret_water/140306_65w65_skew_ret_G_vs_l.pdf')
#show()


