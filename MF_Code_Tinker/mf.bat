%  Cosmological Parameters
%----------------------------------------------------------------------------

GAMMA		0.2
OMEGA_M		0.27
SIGMA_8		0.82
RHO_CRIT	2.775E11		% h^2 M_sol/Mpc^3
SPECTRAL_INDX	0.95
HUBBLE		0.7
OMEGA_B		0.042	
DELTA_CRIT	1.686
ITRANS		5
TF_file		transfunc.WMAP3
DELTA_HALO	360
REDSHIFT	0

root_filename	tinker_15_bins     % [root_filename].dndM will be output file

%-------------------------
% TRANSFER FUNCTION OPTIONS
%-------------------------
% ITRANS = 4 is for the Efstathiou, Bond, White transfer function. Uses GAMMA ONLY.
% ITRANS = 5 is Eisenstein & Hu, uses OMEGA_M, OMEGA_B, HUBBLE
% ITRANS = 11 is an exterior file, must be in CMBFAST format, uses TF_file for filename
%
