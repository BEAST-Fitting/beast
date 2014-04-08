"""
Run the module function that makes chunks first

"""
import sys


def prepare_input_chunks():
    import mf_phat_b12 as datamod
    return datamod.prepare_individual_inputs(datamod.obsfile, 14000)


def merge_inputs_outputs(project=None):
    import mf_phat_b12 as datamod
    import glob
    from beast.external.eztables import Table
    from beast.tools.pbar import Pbar

    obs_base = datamod.obsfile.split('.')
    obs_ext = obs_base[-1]
    obs_base = obs_base[:-1]

    lst = glob.glob('.'.join(obs_base) + '.part*.' + obs_ext)
    if len(lst) == 0:
        raise ValueError('cannot find any chunk. Did you run prepare_individual_inputs?')

    if project is None:
        outname = datamod.project[:]
    else:
        outname = project[:]

    t_final = None
    for chunk, input_fname in Pbar(len(lst)).iterover(enumerate(lst)):
        ti = Table(input_fname)
        s_file = '{0:s}/{0:s}.part{1:d}_stats.fits'.format(outname, chunk)
        to = Table(s_file)

        for k in to.keys():
            ti.addCol(k, to[k])

        if t_final is None:
            t_final = ti
        else:
            t_final.stack(ti.data)

    t_final.write('{0:s}/{0:s}.final_stats.fits'.format(outname))
    return t_final


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
