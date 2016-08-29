"""
BEAST USER INPUT TEMPLATE

MF: this is a python template based on Heddy's ascii file
"""

#-------------------------
#Define your project
#-------------------------
project_name = 'phat_b15'
subproject_name = 'coarse_1000'

#-------------------------
#Define your observation catalog
#-------------------------
datafile = '/home/arab/beast_data/b15f08_sixband.data.v1.fits'
#all PHAT filters
filters = 'hst_wfc3_f275w hst_wfc3_f336w hst_acs_wfc_f475w hst_acs_wfc_f814w hst_wfc3_f110w hst_wfc3_f160w'.upper().split()
distance_modulus = 24.4   # mag

#-------------------------
#Define the model SED grid space
#-------------------------
#---Dust parameters
extLaw_name = 'RvFbumpLaw'

#-Av
Av_min = 0.0
Av_max = 2.0
Av_step = 0.5

#-Rv
Rv_min = 2.0
Rv_max = 6.0
Rv_step = 0.5

#-fb
f_bump_min = 0.0
f_bump_max = 1.0
f_bump_sep = 0.1

#---Stellar parameters
#-Ages
log_age_min = 6.0
log_age_max = 10.0
log_age_step = 0.1

#-Masses
log_mass_min = -0.8
log_mass_max = 2.0
log_mass_step = 0.01

#-Metallicity
# this could / should be iterable I think
Z = 0.004

#-------------------------
# Intermediate product storage
#-------------------------
isofile = 'beast/libs/iso.proposal.fits'
lnp_outname = 'lnp.{0}.{1}.hd5'.format(project_name, subproject_name)
stat_outname = 'stat.{0}.{1}.hd5'.format(project_name, subproject_name)
res_outname = 'res.{0}.{1}.fits'.format(project_name, subproject_name)
spectral_grid_fname = 'coarse.spectral.grid.fits'
sed_grid_fname = 'medium.coarse.sed.grid.fits'





