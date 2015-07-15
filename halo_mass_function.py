#!/usr/bin/python

# Shany Danieli
# July 15, 2015
# Yale University

""" 
The program calculates the halo mass function of a set of data.
"""


####import modules########################################################################
import numpy as np
from math import *
import matplotlib.pyplot as plt
##########################################################################################


def CLF(M,L_bins,dL,Mhost_min,Mhost_max,Mhost):
    """
    Calculate the galaxy conditional luminosity function (CLF).
    
    Parameters
    ----------
    M: array_like
        absolute magnitude of galaxies
    
    L_bins: array_like, optional
        bin edges to use for the luminosity function

    dL: int
        bin width

    Mhalo_min: int, log (Mhalo)
        the lower limit of the host halo mass bin 

    Mhalo_max: int, log (Mhalo)
        the upper limit of the host halo mass bin

    Mhost: array_like
        the mass of the host halos

    
    Returns
    -------
    CLF, L_bins: np.array, np.array
    """

    M_bin = np.asarray([M[i] for i in range(len(M)) if  Mhost_min<log(Mhost[i],10)<=Mhost_max]) # First luminosity bin
    Nhosts = len(M_bin)
    Msun = 4.76 # The Sun's absolute magnitude
    L = ((Msun-M_bin)/2.5) # Calculated the luminosity 
    counts = np.histogram(L,L_bins)[0]
    CLF = counts/(dL*Nhosts)
    return CLF,L_bins



def jackknife_errors_CLF(pos,Phi,Ndivs,Lbox,M,L_bins,dL,Mhost_min,Mhost_max,Mhost):
    """
    Calculate the errors to the CLF using Jackknife resampling method
    
    Parameters
    ----------
    pos: array-like
        Npts x 3 numpy array containing 3d positions of Npts.
        
    Phi: array-like
        the CLF for the full box.

    Ndivs: int
        the number by which we divide every dimension in the Jackknife resampling.

    Lbox: int
        the length of one dimension of the full box.
        
    M: array-like
        the absolute magnitude of the galaxies.
        
    L_bins: array_like
    	numpy array of the left boundaries defining the luminosity bins in which pairs are counted. 
    
	dL: float
		the luminosity bin size.
		
	Mhost_min: float
		the lower bound for the halo mass bin.
		
	Mhost_max: float		
		the upper bound for the halo mass bin.
		
	Mhost: array-like
		Npts numpy array containing the halo masses of Npts.	

			
    Returns
    -------
    Jackknife errors for the different luminosity bins: np.array
    """

    n_subBox = Ndivs*Ndivs*Ndivs # The number of sub volumes for the Jackknife resampling
    V_subBox = Vbox - Vbox/n_subBox # The volume of a Jackknife sample
    N = len(pos) 
    delta = Lbox/Ndivs
    
    # Indices for the galaxies positions
    index = np.asarray([floor(pos[i,0]/delta) + (floor(pos[i,1]/delta)*Ndivs) + (floor(pos[i,2]/delta)*Ndivs*Ndivs) + 1 for i in range(N)]) # index for the position of particle2
    M_sub_sample = [] # keeps the absolute magnitude for the sub-samples
    Mhost_sub_sample = [] # keeps the halo mass for the sub-samples
    CLF_all = []  # keeps the values of the CLF for the full sample and for each of the sub-samples
    CLF_all.append(Phi)
    for k in range(1,n_subBox+1): # run over the sub-samples
        for i in range(0,N): # runs over all the points (galaxies)
                if (index[i] != k): # the point is inside the sub-box
                    M_sub_sample.append(M[i]) # then add to sub-box list
                    Mhost_sub_sample.append(Mhost[i])
        CLF_sub,L_bins = CLF(M_sub_sample,L_bins,dL,Mhost_min,Mhost_max,Mhost_sub_sample)
#        print "CLF for sub-samples:"
#        print CLF_sub
        CLF_all.append(CLF_sub)
        M_sub_sample = []
        Mhost_sub_sample = []

	n_subBox = float(n_subBox)
    full = np.asarray(CLF_all[0]) # the CLF for the full sample
    sub_samples = np.asarray(CLF_all[1:]) # the CLF for the Jackknife sub-samples
    after_subtraction =  sub_samples - np.mean(sub_samples,axis=0)
    squared = after_subtraction**2
    error2 = ((n_subBox-1)/n_subBox)*squared.sum(axis=0)
    errors = error2**0.5
    return errors



def compute_CLF_jackknife(L_bins,dL,M_cen,Mhost_min,Mhost_max,Mhost_cen):
	M = 10**((Mhost_max+Mhost_min)/2)
	Phi_cen,L_bins = CLF(M_cen,L_bins,dL,Mhost_min,Mhost_max,Mhost_cen)
	errors = jackknife_errors_CLF(pos_cen,Phi_cen,Ndivs,Lbox,M_cen,L_bins,dL,Mhost_min,Mhost_max,Mhost_cen)
	y_cen = [log(i,10) for i in Phi_cen if i>0]
	x_cen = [L_bins[i]+dL/2 for i in range(len(L_bins)-1) if Phi_cen[i]>0]
	delta_y_cen = [errors[i] for i in range(len(L_bins)-1) if Phi_cen[i]>0]
	yerr_cen = [i/(10**j) for i,j in zip(delta_y_cen,y_cen)]
	
	return x_cen,y_cen,yerr_cen 
	
	
def compute_CLF_jackknife_sat(L_bins,dL,M_sat,Mhost_min,Mhost_max,Mhost_sat):
	M = 10**((Mhost_max+Mhost_min)/2)
	Phi_sat,L_bins = CLF(M_sat,L_bins,dL,Mhost_min,Mhost_max,Mhost_sat)
	errors = jackknife_errors_CLF(pos_sat,Phi_sat,Ndivs,Lbox,M_sat,L_bins,dL,Mhost_min,Mhost_max,Mhost_sat)
	y_sat = [log(i,10) for i in Phi_sat if i>0]
	x_sat = [L_bins[i]+dL/2 for i in range(len(L_bins)-1) if Phi_sat[i]>0]
	delta_y_sat = [errors[i] for i in range(len(L_bins)-1) if Phi_sat[i]>0]
	yerr_sat = [i/(10**j) for i,j in zip(delta_y_sat,y_sat)]	
	return x_sat,y_sat,yerr_sat



def plot_clf_fit_cen(x_cen,y_cen, yerr,M_host_min,M_host_max):
	label = str(M_host_min) +r'$ < M_{h}  \leq $' + str(M_host_max)
	mid_halo_mass = (M_host_min+M_host_max)/2
	print "mid_halo_mass"
	print mid_halo_mass
	plt.errorbar(x_cen,y_cen, yerr=yerr, capsize=4, ls='none', color='red', elinewidth=2,marker='o',markerfacecolor='red')
	x_fit_cen = np.arange(x_cen[-len(x_cen)]-0.1,x_cen[-1]+0.1,0.05)
	y_fit_plot_cen = [log(n,10) for n in peval(x_fit_cen, fitted_para_cen,10**mid_halo_mass)]
	plt.plot(x_fit_cen,y_fit_plot_cen,color='black')
	plt.xlim(7.5, 12)
	plt.xlabel(r'$\log[L/(h^{-2}L_{\bigodot})])$',fontsize=15)
	plt.ylabel(r'$\log(\Phi(L) d\log L / group)$',fontsize=15)
	plt.title(label)
	plt.savefig('plots/clf_cen_fit_'+str(M_host_min)+'-'+str(M_host_max)+'.png')
	plt.close() 


def plot_clf_fit_sat(x_sat,y_sat, yerr,M_host_min,M_host_max):
	label = str(M_host_min) +r'$ < M_{h}  \leq $' + str(M_host_max)
	mid_halo_mass = (M_host_min+M_host_max)/2
	plt.errorbar(x_sat,y_sat, yerr=yerr, capsize=4, ls='none', color='red', elinewidth=2,marker='o',markerfacecolor='red')
	x_fit_sat = np.arange(x_sat[-len(x_sat)]-0.1,x_sat[-1]+0.1,0.05)
	y_fit_plot_sat = [log(n,10) for n in peval_sat(x_fit_sat, fitted_para_sat,fitted_para_cen,10**mid_halo_mass)]
	plt.plot(x_fit_sat,y_fit_plot_sat,color='black')
	plt.xlim(7.5, 12)
	plt.xlabel(r'$\log[L/(h^{-2}L_{\bigodot})])$',fontsize=15)
	plt.ylabel(r'$\log(\Phi(L) d\log L / group)$',fontsize=15)
	plt.title(label)
	plt.savefig('plots/clf_sat_fit_'+str(M_host_min)+'-'+str(M_host_max)+'.png')
	plt.close() 




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

""" 7/15/2015 - maybe I don't need it 
# Reads the data from the mock 
data1 = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/erased_assembly_bias_Mr_gr_model.dat",float) # no assembly bias
#data1 = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/conditional_abundance_matching_Mr_gr_model.dat",float) # with assembly bias
M = np.asarray(data1[:,7]) # The absoulute magnitude
pos = data1[:,range(1,4)] # the galaxies positions
N = len(M) # the number of data points (galaxies)
"""

data = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/noab_mock.dat",float) # no assembly bias
#data = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/sham_mock.dat",float) # with assembly bias
halo_ID = data[:,0] # Halo ID
Mhost = data[:,2] # Host mass in units of (Msun/h)

print "min(Mhost):"
print min(Mhost)
print "max(Mhost):"
print max(Mhost)

halo_ID = np.asarray(halo_ID)
Mhost = np.asarray(Mhost)

Mass_bins = [10**11.875 for i in range(len(x_cen_0))] + [10**12.125 for i in range(len(x_cen_1))] + [10**12.375 for i in range(len(x_cen_2))] + [10**12.625 for i in range(len(x_cen_3))] +  [10**12.875 for i in range(len(x_cen_4))] + [10**13.125 for i in range(len(x_cen_5))]  + [10**13.375 for i in range(len(x_cen_6))]+[10**13.625 for i in range(len(x_cen_7))]+[10**13.875 for i in range(len(x_cen_8))]+[10**14.125 for i in range(len(x_cen_9))] 


