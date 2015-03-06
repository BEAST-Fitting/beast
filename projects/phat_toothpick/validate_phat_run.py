""" Script to validate a BEAST run against an existing PHAT run
  Exactly the same data needs to be used.

.. history::
    Written 10Nov14 by KDG.

"""
import numpy as np
import matplotlib.pyplot as pyplot
from beast.external.eztables import Table

if __name__ == '__main__':

    #oldresults = Table("/astro/dust_kg/kgordon/new_grid4/new_grid4_stats.fits")
    #oldresults = Table("b15_nov15_test/frozen_version_SD_b15_27_stats_new_grid2.fits")
    #oldresults = Table("b15_nov15_test/b15_nov15_test_stats.fits")
    #oldresults = Table("/astro/dust_kg/beast_share/b15_out/output_catalog_bSD_b15_27.fits")
    #oldresults = Table("b15_late_jan15_test_cmd23/b15_late_jan15_test_cmd23_stats_w_trimgrid.fits")
    #oldresults = Table("b15_late_jan15_test_cmd27/b15_late_jan15_test_cmd27_stats.fits")
    oldresults = Table("b15_late_jan15_test_small/b15_late_jan15_test_small_stats.fits")

    #newresults = Table("b15_late_jan15_test_small/b15_late_jan15_test_small_stats.fits")
    #newresults = Table("b15_late_jan15_test/b15_late_jan15_test_stats.fits")
    #newresults = Table("production_b15/production_b15_sd27_A_stats.fits")
    #newresults = Table("BEAST_production/old_production_b15/production_b15_sd27_A_stats.fits")
    #newresults = Table("b15_late_jan15_test_cmd23/b15_late_jan15_test_cmd23_stats.fits")
    #newresults = Table("b15_late_jan15_test_cmd27_cut/b15_late_jan15_test_cmd27_cut_stats.fits")
    newresults = Table("b15_late_jan15_test_small/b15_late_jan15_test_small_memory_stats.fits")

    print oldresults.keys()
    print newresults.keys()

    fig, ax = pyplot.subplots(nrows=2, ncols=4, figsize=(20,10))

    ax[0,0].plot(oldresults['logT_Exp'],oldresults['Av_Exp'],'bo',label='old')
    ax[0,0].set_xlabel('logT_Exp')
    ax[0,0].set_ylabel('A(V)_Exp')
    ax[0,0].set_ylim(0.0,4.5)
    ax[0,0].legend()

    ax[0,1].plot(newresults['logT_Exp'],newresults['Av_Exp'],'ro',label='new')
    ax[0,1].set_xlabel('logT_Exp')
    ax[0,1].set_ylabel('A(V)_Exp')
    ax[0,1].set_ylim(0.0,4.5)
    ax[0,1].legend()

    ax[1,0].plot(oldresults['logT_p50'],0.5*(oldresults['logT_p84']-oldresults['logT_p16']),'bo')
    ax[1,0].set_xlabel('logT_p50')
    ax[1,0].set_ylabel('logT_unc')
    ax[1,0].set_ylim(0.0,0.45)

    ax[1,1].plot(newresults['logT_p50'],0.5*(newresults['logT_p84']-newresults['logT_p16']),'ro')
    ax[1,1].set_xlabel('logT_p50')
    ax[1,1].set_ylabel('logT_unc')
    ax[1,1].set_ylim(0.0,0.45)

    ax[0,2].plot(oldresults['logT_Exp'],oldresults['logT_p50'],'bo')
    ax[0,2].plot(newresults['logT_Exp'],newresults['logT_p50'],'ro')
    ax[0,2].set_xlabel('logT_Exp')
    ax[0,2].set_ylabel('logT_p50')

    ax[1,2].plot(oldresults['Av_Exp'],oldresults['Av_p50'],'bo')
    ax[1,2].plot(newresults['Av_Exp'],newresults['Av_p50'],'ro')
    ax[1,2].set_xlabel('A(V)_Exp')
    ax[1,2].set_ylabel('A(V)_p50')

    ax[0,3].plot(oldresults['logT_Best'],newresults['logT_Best'],'go',label='Best')
    ax[0,3].plot(oldresults['logT_p50'],newresults['logT_p50'],'ro',label='p50')
    ax[0,3].plot(oldresults['logT_Exp'],newresults['logT_Exp'],'bo',label='Exp')
    ax[0,3].set_xlabel('logT(old)')
    ax[0,3].set_ylabel('logT(new)')
    ax[0,3].legend(loc=2)

    ax[1,3].plot(oldresults['Av_Best'],newresults['Av_Best'],'go',label='Best')
    ax[1,3].plot(oldresults['Av_p50'],newresults['Av_p50'],'ro',label='p50')
    ax[1,3].plot(oldresults['Av_Exp'],newresults['Av_Exp'],'bo',label='Exp')
    ax[1,3].set_xlabel('Av(old)')
    ax[1,3].set_ylabel('Av(new)')
    ax[1,3].legend(loc=2)

    pyplot.show()
