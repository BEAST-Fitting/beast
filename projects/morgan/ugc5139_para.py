import sys


def merge_inputs_outputs(project=None):
    import ugc5139 as datamod
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
        ti = Table(input_fname, silent=True)
        s_file = '{0:s}/{0:s}.part{1:d}_stats.fits'.format(outname, chunk)
        to = Table(s_file, silent=True)

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
        print('must have a chunk number')
        sys.exit(1)
    chunk = int(sys.argv[1])
    import ugc5139
    g = '{project:s}/{project:s}_seds.grid.hd5'.format(project=ugc5139.project)
    ugc5139.run_chunk_fit(ugc5139.project, g, chunk)
