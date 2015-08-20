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
from __future__ import division
import numpy as np
from math import *
import matplotlib.pyplot as plt
from pylab import *
from yt.mods import *
from yt.analysis_modules.halo_mass_function.api import *
import h5py
import os, sys
##########################################################################################


"""
Main program
"""

# General useful quantities
Lbox = 250.0
Vbox = Lbox**3


"""
Data processing 
"""

#data = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/noab_mock_fixed.dat",float) # no assembly bias
data = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/noab_mock.dat",float) # no assembly bias
#data = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/noab_mock_downsized.dat",float) # no assembly bias - downsized for a halo mass range []
#data = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/sham_mock.dat",float) # with assembly bias


"""
Analytical hmf
"""

num_bins_model = 15 # number of mass bins
# importing data from the analytical calculation done with Tinker's code (c code)
masses_model = loadtxt("/Users/si224/Documents/2014/vanDenBosch/Code/mass_function/MF_Code_Tinker/tinker_"+str(num_bins_model)+"_bins.dndM",float)[:,0]
hmf_model = loadtxt("/Users/si224/Documents/2014/vanDenBosch/Code/mass_function/MF_Code_Tinker/tinker_"+str(num_bins_model)+"_bins.dndM",float)[:,1]
factor = log(10)*masses_model  # Converting from dndM to dndlogM
hmf_model = factor*hmf_model

print "masses_model:"
print np.log10(masses_model)
print "\n"
print "hmf_model"
print hmf_model
print "\n\n"




"""
Computing hmf from data - NEW
"""

V = 250.0**3 # box volume

# Host mass bins
num_bins = 15
Mmin = 11.75 # left edge
Mmax = 14.25 # right edge
dlogM = (Mmax-Mmin)/num_bins
Mass_bins = np.arange(Mmin,Mmax+dlogM/2,dlogM)
Mass_centers = (Mass_bins+dlogM/2).tolist()
del Mass_centers[-1]

print "\n\n"
print "Mass_centers:"
print Mass_centers
print "\n\n"

"""
# Host mass bins
num_bins = 20
#Mmin = 10.4
#Mmax = 15.0
Mmin = 11.85
Mmax = 14.15
dlogM = (Mmax-Mmin)/(num_bins-1)
Mass_centers = np.arange(Mmin,Mmax+dlogM/2,dlogM)
Mass_bins = (Mass_centers-dlogM/2).tolist()
Mass_bins.append(Mass_bins[-1]+dlogM)
"""


# randomize - without assembly bias
#data1_cen = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/splitted/cen_mag.txt",float)
#data2_cen = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/splitted/cen.txt",float)
data1_cen = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/downsized/cen_mag_ds.txt",float)
data2_cen = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/downsized/cen_ds.txt",float)

host_mass_cen = data2_cen[:,2]
counts = np.histogram(np.log10(host_mass_cen), Mass_bins)[0]
plt.close()
dndlogM = counts/(dlogM*V)



print "diff:"
print 100*(dndlogM-hmf_model)/dndlogM



plt.plot(Mass_centers,np.log10(dndlogM),marker='o', color='c',linestyle='None',markeredgecolor='k',markeredgewidth=1,markersize=7,label='mock - noab')
plt.plot(np.log10(masses_model),np.log10(hmf_model),marker='o', color='b',linestyle='None',markeredgecolor='k',markeredgewidth=1,markersize=7,label='model - tinker')
#plt.plot(np.log10(masses_model),np.log10(hmf_model),label='model',c='k', linewidth=1.0)
plt.legend()
plt.xlabel(r'M(h$^{-1}$M$_{\bigodot}$))',fontsize=15)
plt.ylabel(r'dn/dlogM(h$^3$Mpc$^{-3}$)',fontsize=15)
#plt.axvline(x=11.5,linewidth=1.0, color='m')
plt.savefig('hmf_mock_vs_model.png')
plt.show()
np.savetxt('hmf_data.dat',dndlogM)

print "************************* DONE *************************"




"""
Computing hmf from data
"""
"""
num_bins_data = 3
dm = (14.25-11.75)/(num_bins_data-1)

#print "dm:"
#print dm

# Downsizing the sample
Mhalo_tmp = data[:,2] # Host mass in units of (Msun/h)
Mhalo = np.asarray([Mhalo_tmp[i] for i in range(len(Mhalo_tmp)) if (11.75-dm/2)<log10(Mhalo_tmp[i])<=(14.25+dm/2)])
upid_tmp = data[:,1] # Host mass in units of (Msun/h)
upid = np.asarray([upid_tmp[i] for i in range(len(upid_tmp)) if (11.75-dm/2)<log10(Mhalo_tmp[i])<=(14.25+dm/2)])
N = len(upid)
print "len:"
print N

Mhost=[]
for i in range(len(Mhalo)):
    if upid[i]<0:
        Mhost.append(Mhalo[i])
Mhost = np.log10(Mhost)
Mmin = 11.75
Mmax = 14.25
dlogM = (Mmax-Mmin)/(num_bins_data-1)
#dM = 10**dlogM
bin_centers_data = np.arange(Mmin,Mmax+dlogM/2,dlogM)
bins = (bin_centers_data-dlogM/2).tolist()
bins.append(bins[-1]+dlogM)
n = np.histogram(Mhost, bins)[0]
dndlogM = n/(dlogM*Vbox)
plt.close()
np.savetxt('hmf_data.dat',dndlogM)

print "dndlogM:"
print dndlogM
"""














"""
Computing hmf from data for the SMF
"""
"""
num_bins_data = 15


f = h5py.File('/Users/si224/Documents/2014/vanDenBosch/Data/Duncan/sm_9.5_s0.0_sfr_c-0.25_250.hdf5', 'r') #open mock in 'read' mode
catalogue = 'sm_9.5_s0.0_sfr_c-0.25_250'
mock = f.get(catalogue)

halo_ID = mock['id'] 
upid = mock['pid']
Mhalo = mock['mvir'] 
Mhost=[]
N = len(halo_ID)

for i in range(len(Mhalo)):
    if upid[i]<0:
        Mhost.append(Mhalo[i])
Mhost = np.log10(Mhost)
Mmin = 9.0
Mmax = 14.0
dlogM = (Mmax-Mmin)/(num_bins_data-1)
#dM = 10**dlogM
bin_centers_data = np.arange(Mmin,Mmax+dlogM/2,dlogM)
bins = (bin_centers_data-dlogM/2).tolist()
bins.append(bins[-1]+dlogM)
n = np.histogram(Mhost, bins)[0]
dndlogM = n/(dlogM*Vbox)
plt.close()
np.savetxt('hmf_data.dat',dndlogM)

plt.loglog([10**x for x in bin_centers_data] ,dndlogM,linestyle='None',marker='.',label='Simulation')
plt.show()
"""




	
"""
diff = (dndlogM - hmf_model)/dndlogM
print "diff:"
print diff


# Plotting 
plt.figure()
plt.subplot(2,1,1)
plt.loglog([10**x for x in bin_centers_data] ,dndlogM,linestyle='None',marker='.',label='Simulation')
plt.loglog(masses_model,hmf_model,label='Tinker')
plt.legend()
#plt.xlim(10**11.5,10**14.5)
#plt.ylim(10**-5,10**0)
plt.xlabel(r'$M(h^{-1}M_{\bigodot}))$',fontsize=15)
plt.ylabel(r'$dn/dlog(M)(h^3Mpc^{-3})$',fontsize=15)
np.savetxt('hmf_data.dat',dndlogM)
#plt.savefig('hmf_11.75-14.25.png')

plt.subplot(2,1,2)
plt.semilogx(masses_model,diff,linestyle='None',marker='o',markersize=5,markeredgecolor='k', markerfacecolor='m')
plt.ylim(-0.4,0.4)
plt.axhline(linewidth=0.5, color='k')
plt.axhspan(-0.05, 0.05, facecolor='0.5',linewidth=0)
plt.axhspan(-0.05, -0.1, facecolor='0.7',linewidth=0)
plt.axhspan(0.05, 0.1, facecolor='0.7',linewidth=0)
plt.axhspan(0.1, 0.2, facecolor='0.9',linewidth=0)
plt.axhspan(-0.1, -0.2, facecolor='0.9',linewidth=0)
plt.plot(masses_model,diff,linestyle='None',marker='o',markersize=5,markeredgecolor='k', markerfacecolor='m')
plt.xlabel(r'$M(h^{-1}M_{\bigodot}))$',fontsize=15)
plt.ylabel(r'$(mock-tinker)/mock$',fontsize=15)
#plt.savefig('mockVsTinker_'+str(num_bins_model)+'.png')
plt.show()
"""

"""
plt.savefig('mockVsTinker_'+str(num_bins_model)+'.png')
plt.show()
#np.savetxt('mass_mid_bins.dat',bin_centers)
#plt.close()
"""

"""
#bins = np.arange(11.75,14.25,dlogM)
#bin_centers = (bins[:-1]+bins[1:])/2.0
#dx = bins[1:]-bins[:-1]

#diff = (hmf_analytic-dndlogM)/hmf_analytic
#diff = (hmf_analytic-dndlogM)
"""









"""
# Analytical calculation using the code from the yt-project
hmf = HaloMassFcn(omega_baryon0=0.042, omega_matter0=0.27,
                  omega_lambda0=0.73, hubble0=0.7, sigma8 = 0.82, this_redshift=0,
                  log_mass_min=11.75, log_mass_max=14.25, fitting_function=5,num_sigma_bins = num_bins)
hmf.write_out(prefix='hmf', analytic=True, simulated=True)
#plt.loglog(hmf.masses_analytic, hmf.n_cumulative_analytic,label = 'Tinker')
#plt.title('hmf_Tinker')
#plt.xlabel(r'$M(h^{-1}M_{\bigodot}))$',fontsize=15)
#plt.ylabel(r'$dn/dlog(M)(h^3Mpc^{-3})$',fontsize=15)
#plt.show()
#plt.savefig('hmf_Press-Schechter.png')
#plt.close()
"""




























