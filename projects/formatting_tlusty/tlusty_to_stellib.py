import numpy as np
import glob
import re
from sedfitter.tools import progressbar
from sedfitter.core import stellib
from sedfitter.core import grid 
from sedfitter.external.eztables import table
from scipy.interpolate import InterpolatedUnivariateSpline

#The TLUSTY models come as a bunch of separate text files. I'm assuming 
# you have these text files, organized into the folder structure
# path_to_tlusty_models/(O or B)STAR*/(metallicity letter)models/(individual models).flux
# The metallicity letters are listed later on. 
#
#I haven't optimized anything, since this will hopefully be something you 
# only need to run once. It'll take about an hour or two to convert all of the 
# TLUSTY atmospheres.
#
path_to_tlusty_models = '/Users/K/Documents/PHAT_SED/Stellar_Atmospheres/TLUSTY' 
out_path = 'sedfitter/libs/tlusty.grid.fits'
filepaths = glob.glob(path_to_tlusty_models+'/*STAR*/*models/*.flux')

def edges_from_centers(centers):
  length = centers.size
  edges = np.zeros(length+1)
  edges[1:-1] = 0.5*(centers[1:] + centers[:-1])
  edges[0] = centers[1] - (centers[2] - centers[1])
  edges[-1] = centers[-2] + (centers[-2] - centers[-3])
  return edges

def sanitize_tlusty_input(frequency, flux):
  #This method deals with two problems in the TLUSTY input files-
  # rows that are not in order and consecutive frequencies that
  # are so close together that the difference is smaller than the 
  # cutoff imposed by the fixed width of the frequency column.
  #I deal with duplicates by shifting the frequency a tiny bit higher. 
  # This isn't in any way correct, but is a fairly small effect, especially
  # considering the coarseness of the Kurucz wavelength grid relative to the
  # TLUSTY one. 
  order = np.argsort(frequency)
  frequency = frequency[order]
  flux = flux[order]
  duplicates = np.where(frequency[1:]==frequency[:-1])[0]
  duplicate_frequencies = np.unique(frequency[duplicates])
  duplicate_sequences = [np.where(frequency==freq)[0] for freq in duplicate_frequencies]
  for duplicate_sequence in duplicate_sequences:
    frequency[duplicate_sequence] = np.linspace(frequency[duplicate_sequence[0]], 
                                                frequency[duplicate_sequence[-1]+1],
                                                len(duplicate_sequence), 
                                                endpoint=False)
  return frequency, flux

kurucz = stellib.Kurucz()
c = 2.998e18 #in Angstroms per second
ck_waves = kurucz.wavelength
ck_reversed_frequencies = c/ck_waves
#ck_frequencies = ck_reversed_frequencies[::-1]
#ck_frequency_bounds = edges_from_centers(ck_frequencies)
ck_wave_bounds = edges_from_centers(ck_waves)
ck_wave_spans = ck_wave_bounds[1:]-ck_wave_bounds[:-1]

Z_sol = 0.02
Z_dict = {'C':2., 'G':1., 'L':1./2., 'S':1./5., 'T':1./10., 'V':1./30., 'W':1./50., 'X':1./100., 'Y':1./1000., 'Z':0.}

spectrum = np.zeros([len(filepaths),len(ck_waves)], dtype='float')
log_g_vals = np.zeros(len(filepaths), dtype='float')
T_eff_vals = np.zeros(len(filepaths), dtype='float')
Z_vals = np.zeros(len(filepaths), dtype='float')

p = re.compile('([CGLSTVWXYZ])([\d]+)g([\d]+)')

with progressbar.PBar(len(filepaths), txt='Reading %d files' %(len(filepaths))) as Pbar:
  for e, filepath in enumerate(filepaths):
    f = open(filepath)
    text = (" ".join(f.readlines())).replace('D','').replace('+','e+').replace('-','e-')
    f.close()
    tl_fluxes_and_frequencies = np.double(np.array(text.split()))
    tl_frequencies = tl_fluxes_and_frequencies[::2][::-1]
    tl_fluxes = tl_fluxes_and_frequencies[1::2][::-1]
    tl_frequencies, tl_fluxes = sanitize_tlusty_input(tl_frequencies, tl_fluxes)
    tl_wavelengths = c/tl_frequencies[::-1]
    tl_fluxes = tl_fluxes[::-1] * (c/(tl_wavelengths*tl_wavelengths))
    #finite_sel = np.where(np.isfinite(tl_fluxes))
    #tl_fluxes = tl_fluxes[finite_sel]
    #tl_wavelengths = tl_wavelengths[finite_sel]
    spl = InterpolatedUnivariateSpline(tl_wavelengths, tl_fluxes, k=1)
    vect_int = np.vectorize(spl.integral)
    integrated_flux = vect_int(ck_wave_bounds[:-1], ck_wave_bounds[1:])
    flux_per_wavelength = integrated_flux/ck_wave_spans
    spectrum[e] = flux_per_wavelength*4*np.pi
    m = p.findall(filepath)[0]
    Z_vals[e] = Z_dict[m[0]]
    T_eff_vals[e] = np.float(m[1])
    log_g_vals[e] = np.float(m[2])/100.

    Pbar.update(e)

g = table.Table(name='Tlusty')
g.addCol('Z', Z_vals*Z_sol, dtype='float')
g.addCol('Teff', T_eff_vals, unit='K', dtype='float')
g.addCol('logG', log_g_vals, dtype='float')


tlusty = grid.MemoryGrid(ck_waves, seds=spectrum, grid=g)
tlusty._backend.writeFITS(out_path, clobber=True)
