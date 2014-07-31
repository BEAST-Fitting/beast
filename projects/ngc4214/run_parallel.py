__doc__ = """
Running fits in batch mode on chunks of input data

    > python run_parallel <chunk#> [--models] [--prep] [-?]

    --models        generates the models only (needed once)
    --prep          prepare chunks from observation file given in the datamodel
    -?, --help      display this message
"""
import sys
import datamodel
import noisemodel
from pipeline import prepare_individual_inputs, run_chunk_fit, make_models

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Wrong number of arguments. See usage -?')
        sys.exit(1)

    if ('-?' in sys.argv[1]) or ('--help' in sys.argv[1]):
        print(__doc__)
    elif '-models' in sys.argv[1]:
        make_models()
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
        noisefile = '{project:s}/{project:s}_noisemodel.hd5'.format(project=datamodel.project)
        noise = noisemodel.get_noisemodelcat(noisefile)

        run_chunk_fit(datamodel.project, g, chunk, noise=noise, obsfile=datamodel.obsfile)
