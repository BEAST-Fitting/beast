""" Script to merge the single camera PHAT ASTs into a single file
      basically, mocking up the 6 band ASTs files by merging and renaming columns

.. history::
    Written 5Nov14 by KDG.

"""
import numpy as np
#import tables
from beast.external.eztables import Table

def merge_phat_asts(uvfile, optfile, irfile, outfile):

    out_dict = {}
    filters = ['F275W','F336W','F475W','F814W','F110W','F160W']

    uv_asts = Table(uvfile)
    opt_asts = Table(optfile)
    ir_asts = Table(irfile)
    asts = [uv_asts,opt_asts,ir_asts]

    # get the maximum number of ASTs
    #  the different cameras can have different #s of ASTs
    lens = np.array([len(uv_asts),len(opt_asts),len(ir_asts)])
    max_len = max(lens)
    
    bad_val = 99.9999999
    magin = np.full(max_len,bad_val)
    magout = np.full(max_len,bad_val)
    for i in range(6):
        magin[:] = bad_val
        magout[:] = bad_val
        magstr = str(i%2 + 1)
        magin[0:lens[i/2]] = asts[i/2]['MAG'+magstr+'IN']
        magout[0:lens[i/2]] = asts[i/2]['MAG'+magstr+'OUT']

        out_dict[filters[i]+'_IN'] = np.array(magin)
        out_dict[filters[i]+'_VEGA'] = np.array(magout)

    out_table = Table(out_dict)

    out_table.write(outfile)

if __name__ == '__main__':

    Merge_PHAT_ASTs('PHAT_camera_AST/fake_stars_b15_27_uv.fits','PHAT_camera_AST/fake_stars_b15_27_opt.fits',
                    'PHAT_camera_AST/fake_stars_b15_27_ir.fits','PHAT_camera_AST/fake_stars_b15_27_all.hd5')
