from graphviz import Digraph
import numpy as np
from collections import defaultdict

__all__ = ["plot_graphic_file_flow"]


def plot_graphic_file_flow(n_sd=1, n_sub=3, savefig="png"):
    """
    Plot a model showing what files are created at each stage of a BEAST
    production run.

    Parameters
    ----------
    n_sd : int (default=1)
        number of source density bins to show

    n_sub : int (default=3)
        number of source density sub-bins to show

    savefig : str (default='png')
        set to the file extension of desired file to save image of model
    """

    # initialize graph
    # (note: can't do dictionary for nodes like plot_graphic_model.py, because
    # we need to exert more control over which things are lined up)
    graph = Digraph(node_attr={'shape': 'box'})

    # first layer
    with graph.subgraph() as sg:
        sg.attr(rank='same')
        sg.node('sed','SEDgrid')
        sg.node('phot','phot')
        sg.node('fake','phot_fake')

    # additional layers
    sg1 = Digraph(node_attr={'shape': 'box'})
    sg2 = Digraph(node_attr={'shape': 'box'})
    sg3 = Digraph(node_attr={'shape': 'box'})
    sg4 = Digraph(node_attr={'shape': 'box'})
    sg5 = Digraph(node_attr={'shape': 'box'})
    _ = [x.attr(rank='same') for x in [sg1, sg2, sg3, sg4, sg5]]

    # initialize dict of edges
    edges = defaultdict(list)

    # iterate through source density bins
    for s in range(n_sd):

        curr_sd = f'SD{s}-{s+1}'

        # files for this source density bin
        sg1.node(f'phot{s}', f'phot_{curr_sd}')
        sg1.node(f'fake{s}', f'phot_fake_{curr_sd}')
        sg2.node(f'obs{s}', f'obsmodel_{curr_sd}')

        # make edges
        edges['phot'].append(f'phot{s}')
        edges['fake'].append(f'fake{s}')
        edges[f'fake{s}'].append(f'obs{s}')

        # files for combo of source density and sub-bin
        for b in range(n_sub):
            curr_sub = f'sub{b}'

            # -- nodes --
            sg3.node(f'phot{s}s{b}', f'phot_{curr_sd}_{curr_sub}')
            sg3.node(f'sed{s}s{b}', f'SEDgrid_{curr_sd}_{curr_sub}_trim')
            sg3.node(f'obs{s}s{b}', f'obsmodel_{curr_sd}_{curr_sub}_trim')
            sg4.node(f'lnps_{s}s{b}', f'likelihoods for\n{curr_sd}_{curr_sub}')
            sg5.node(f'stat{s}s{b}', f'stats_{curr_sd}_{curr_sub}')
            sg5.node(f'pdf1d{s}s{b}', f'pdf1d_{curr_sd}_{curr_sub}')
            sg5.node(f'pdf2d{s}s{b}', f'pdf2d_{curr_sd}_{curr_sub}')
            sg5.node(f'lnp{s}s{b}', f'lnp_{curr_sd}_{curr_sub}')


            # -- edges --
            # phot to sub-phot
            edges[f'phot{s}'].append(f'phot{s}s{b}')
            # SED to trimmed SED
            edges['sed'].append(f'sed{s}s{b}')
            # obsmodel to trimmed obsmodel
            edges[f'obs{s}'].append(f'obs{s}s{b}')
            # photometry + trimmed SED + trimmed obsmodel to likelihoods
            edges[f'phot{s}s{b}'].append(f'lnps_{s}s{b}')
            edges[f'sed{s}s{b}'].append(f'lnps_{s}s{b}')
            edges[f'obs{s}s{b}'].append(f'lnps_{s}s{b}')
            # likelihoods to output files
            edges[f'lnps_{s}s{b}'] += [
                f'stat{s}s{b}', f'pdf1d{s}s{b}', f'pdf2d{s}s{b}', f'lnp{s}s{b}'
            ]


    # add the rest of the edges
    for ckey in edges.keys():
        for cval in edges[ckey]:
            graph.edge(ckey, cval)

    # append subgraphs
    graph.subgraph(sg1)
    graph.subgraph(sg2)
    graph.subgraph(sg3)
    graph.subgraph(sg4)
    graph.subgraph(sg5)

    # make the graph flow horizontally
    graph.graph_attr['rankdir'] = 'LR'

    # save it
    graph.render("beast-file-flow", format=savefig)

if __name__ == "__main__":
    plot_graphic_file_flow(n_sd=2, n_sub=3)
