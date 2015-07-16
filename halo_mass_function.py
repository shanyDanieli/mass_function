#!/usr/bin/python
# coding: latin-1


# Shany Danieli
# July 15, 2015
# Yale University

""" 
The program calculates the halo mass function for a simulated dataset
and/or an analytical fit over a given mass range for a set of specified cosmological parameters.
"""


####import modules########################################################################
import numpy as np
from math import *
import matplotlib.pyplot as plt
from yt.mods import *
from yt.analysis_modules.halo_mass_function.api import *
import os, sys
##########################################################################################


"""
Main program
"""

# General useful quantities
Lbox = 250.0
Vbox = Lbox**3
Ndivs = 5


"""
Data processing 
"""

data = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/noab_mock.dat",float) # no assembly bias
#data = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/sham_mock.dat",float) # with assembly bias
halo_ID = data[:,0] # Halo ID
Mhost = data[:,2] # Host mass in units of (Msun/h)

halo_ID = np.asarray(halo_ID)
Mhost = np.asarray(Mhost)	

print "len(Mhost)"
print len(Mhost)

mass_mid = [11.875 + 0.25*x for x in range(10)]
mass_bins = np.asarray([11.75 + 0.25*x for x in range(11)]) # list of theleft boundaries defining the host halo mass bins)
counter = np.zeros(10)

for i in range(10):
	for j in range(len(Mhost)):
		if (mass_bins[i] < log(Mhost[j],10) < mass_bins[i+1]):
			counter[i] = counter[i] + 1
			
			
dndlogM = [log(x,10) for x in counter/Vbox]

"""
# plotting the halo mass function
#plt.errorbar(x_cen,y_cen, yerr=yerr, capsize=4, ls='none', color='red', elinewidth=2,marker='o',markerfacecolor='red')
plt.loglog([10**x for x in mass_mid],dndlogM,color='k',marker='x')
#plt.xlim(7.5, 12)
plt.xlabel(r'$M(h^{-1}M_{\bigodot}))$',fontsize=15)
plt.ylabel(r'$dn/dlog(M)(h^3Mpc^{-3})$',fontsize=15)
plt.show()
#plt.title(label)
#plt.savefig('plots/clf_cen_fit_'+str(M_host_min)+'-'+str(M_host_max)+'.png')	
#plt.close() 
"""


hmf = HaloMassFcn(omega_baryon0=0.042, omega_matter0=0.27,
                  omega_lambda0=0.73, hubble0=0.7, sigma8 = 0.82, this_redshift=0,
                  log_mass_min=11.25, log_mass_max=14.25, fitting_function=1)
hmf.write_out(prefix='hmf', analytic=True, simulated=True)
plt.loglog(hmf.masses_analytic, hmf.n_cumulative_analytic)
plt.show()

	
	


