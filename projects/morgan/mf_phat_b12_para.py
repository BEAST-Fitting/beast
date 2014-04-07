"""
Run the module function that makes chunks first

"""
import sys


def prepare_input_chunks():
    import mf_phat_b12 as datamod
    return datamod.prepare_individual_inputs(datamod.obsfile, 14000)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('must have a chunk number or -prep option')
        sys.exit(1)
    if '-prep' in sys.argv[1]:
        prepare_input_chunks()
    else:
        chunk = int(sys.argv[1])
        # ==================================================================
        # MODULE TO CHANGE
        # note: import at the final stage only to avoid reading files if not
        # running
        import mf_phat_b12 as datamod
        # ==================================================================
        g = '{project:s}/{project:s}_seds.grid.hd5'.format(project=datamod.project)
        datamod.run_chunk_fit(datamod.project, g, chunk)
