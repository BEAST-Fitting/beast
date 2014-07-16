__doc__ = """
Running fits in batch mode on chunks of input data

    > python run_parallel <chunk#> [-prep] [-?]

    -prep   prepare chunks from observation file given in the datamodel
    -?      display this message
"""
import sys
import datamodel
import noisemodel
from pipeline import prepare_individual_inputs, run_chunk_fit

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Wrong number of arguments. See usage -?')
        sys.exit(1)

    if '-?' in sys.argv[1]:
        print __doc__
    elif '-prep' in sys.argv[1]:
        # 14000 stars per file is HDF5 performance limited
        prepare_individual_inputs(datamodel.obsfile, 14000)
    else:
        chunk = int(sys.argv[1])
        # ==================================================================
        # MODULE TO CHANGE
        # note: import at the final stage only to avoid reading files if not
        # running
        # ==================================================================
        g = '{project:s}/{project:s}_seds.grid.hd5'.format(project=datamodel.project)
        ast = noisemodel.get_noisemodelcat(datamodel.astfile)
        run_chunk_fit(datamodel.project, g, chunk, ast=ast,
                      obsfile=datamodel.obsfile)
