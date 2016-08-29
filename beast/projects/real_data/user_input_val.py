"""
BEAST USER INPUT TEMPLATE

MF: this is a python template based on Heddy's ascii file
"""

#-------------------------
#Define your project
#-------------------------
project_name = 'val_b1'
subproject_name = 'test5_1000'

#-------------------------
#Define your observation catalog
#-------------------------
obsfile = '/home/arab/beast_data/brick1_dec11_b1_best_fit_param.fits'
#all PHAT filters
filters = 'hst_wfc3_f275w hst_wfc3_f336w hst_acs_wfc_f475w hst_acs_wfc_f814w hst_wfc3_f110w hst_wfc3_f160w'.upper().split()
distanceModulus = 24.4   # mag

#-------------------------
#Define the model SED grid space
#-------------------------
#---Dust parameters
extLaw_name = 'RvFbumpLaw'

#-Av
#Av_min = 0.0
#Av_max = 4.0
#Av_step = 0.5
Av_min = 0.0
Av_max = 10.0
Av_step = 0.1


#-Rv
#Rv_min = 2.0
#Rv_max = 6.0
#Rv_step = 0.5
Rv_min = 1.5
Rv_max = 6.0
Rv_step = 0.5

#-fb
#f_bump_min = 0.0
#f_bump_max = 1.0
#f_bump_step = 0.1
f_bump = 1.0

#---Stellar parameters
#-Ages
#log_age_min = 6.0
#log_age_max = 10.0
#log_age_step = 0.1
log_age_min = 6.0
log_age_max = 10.0
log_age_step = 0.1

#-Masses
#log_mass_min = -0.8
#log_mass_max = 2.0
#log_mass_step = 0.01
log_mass_min = -0.8
log_mass_max = 2.0
log_mass_step = 0.01#0.01

#-Metallicity
# this could / should be iterable I think
#Z_min = 0.02
#Z_max = 0.03
#Z_step = 0.01
Z = 0.02

#-------------------------
# Intermediate product storage
#-------------------------
dir = "/astro/dust_kg/harab/beast/projects/real_data/"
isofile = 'beast/libs/iso.proposal.fits'
lnp_outname = dir+'lnp.{0}.{1}.hd5'.format(project_name, subproject_name)
stat_outname = dir+'stat.{0}.{1}.hd5'.format(project_name, subproject_name)
res_outname = dir+'res.{0}.{1}.fits'.format(project_name, subproject_name)
spectral_grid_fname = dir+'test5.spectral.grid.fits'
sed_grid_fname = dir+'test5.sed.grid.fits'





