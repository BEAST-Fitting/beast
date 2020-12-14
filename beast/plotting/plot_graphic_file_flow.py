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

    # initialize dict of nodes
    nodes = {}

    # initialize dict of edges
    edges = defaultdict(list)

    # iterate through source density bins
    for s in range(n_sd):

        curr_sd = f'SD{s}-{s+1}'

        # files for this source density bin
        nodes[f'phot{s}'] = f'phot_{curr_sd}'
        nodes[f'fake{s}'] = f'phot_fake_{curr_sd}'
        nodes[f'obs{s}'] = f'obsmodel_{curr_sd}'

        # make edges
        edges['phot'].append(f'phot{s}')
        edges['fake'].append(f'fake{s}')
        edges[f'fake{s}'].append(f'obs{s}')

        # files for combo of source density and sub-bin
        for b in range(n_sub):
            curr_sub = f'sub{b}'

            # -- nodes --
            nodes[f'phot{s}s{b}'] = f'phot_{curr_sd}_{curr_sub}'
            nodes[f'sed{s}s{b}'] = f'SEDgrid_{curr_sd}_{curr_sub}_trim'
            nodes[f'obs{s}s{b}'] = f'obsmodel_{curr_sd}_{curr_sub}_trim'
            nodes[f'lnps_{s}s{b}'] = f'likelihoods for\n{curr_sd}_{curr_sub}'
            nodes[f'stat{s}s{b}'] = f'stats_{curr_sd}_{curr_sub}'
            nodes[f'pdf1d{s}s{b}'] = f'pdf1d_{curr_sd}_{curr_sub}'
            nodes[f'pdf2d{s}s{b}'] = f'pdf2d_{curr_sd}_{curr_sub}'
            nodes[f'lnp{s}s{b}'] = f'lnp_{curr_sd}_{curr_sub}'

            # -- edges --
            # phot to sub-phot
            edges[f'phot{s}'].append(f'phot{s}s{b}')
            # sub-phot to trimmed SED/obsmodel
            #edges[f'phot{s}s{b}'] += [f'sed{s}s{b}', f'obs{s}s{b}']
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

    #import pdb; pdb.set_trace()
    # initialize graph
    graph = Digraph()

    # first layer
    with graph.subgraph() as sg:
        sg.attr(rank='same')
        sg.node('sed','SEDgrid')
        sg.node('phot','phot')
        sg.node('fake','phot_fake')

    # add the rest of the nodes
    for ckey in nodes.keys():
        graph.node(ckey, nodes[ckey])

    # add the rest of the edges
    for ckey in edges.keys():
        for cval in edges[ckey]:
            graph.edge(ckey, cval)

    # make the graph flow horizontally
    graph.graph_attr['rankdir'] = 'LR'

    # save it
    graph.render("beast-file-flow", format=savefig)

if __name__ == "__main__":
    plot_graphic_file_flow(n_sd=2, n_sub=3)
