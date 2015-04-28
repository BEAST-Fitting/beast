from __future__ import print_function
import argparse
import tables

import numpy as np

if __name__ == '__main__':
    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("lnp_file", help="filename for lnp file")
    args = parser.parse_args()

    f = tables.openFile(args.lnp_file, 'r')

    nobs = f.root._v_nchildren - 2

    print('# stars = ', nobs)
    for obj in range(nobs):
        lnps = f.getNode('/star_{0:d}/lnp'.format(obj)).read().astype(np.float64)
        indx = f.getNode('/star_{0:d}/idx'.format(obj)).read().astype(int)
        print(obj, np.amin(lnps), np.amax(lnps), np.amin(indx), np.amax(indx))
    
    f.close()
