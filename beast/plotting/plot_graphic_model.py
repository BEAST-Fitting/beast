from graphviz import Digraph
import numpy as np

__all__ = ["create_graphic_model", "plot_graphic_model"]


def create_graphic_model(nodes, edges, gtype):
    """
    Create a graphic model given nodes and edges

    Parameters
    ----------
    nodes : dict
        for each node {key, text, math}
    edges : dict
        for each edge {key, text, math}
    gtype : str [default="text"]
        "text" for a verbose version, "math" for a compact version
    """
    mod = Digraph()
    if gtype == "math":
        tindx = 1
    else:
        tindx = 0
    for ckey in nodes.keys():
        if ckey == "Like":
            cstyle = "filled"
        else:
            cstyle = None
        mod.node(ckey, nodes[ckey][tindx], style=cstyle)

    for ckey in edges.keys():
        for cval in np.atleast_1d(edges[ckey]):
            mod.edge(ckey, cval)

    return mod


def plot_graphic_model(gtype="text", savefig="png"):
    """
    Plot the graphical model of the BEAST.

    Parameters
    ----------
    gtype : str [default="text"]
        "text" for a verbose version, "math" for a compact version
    savefig : str
        set to the file extension of desired file to save image of model
    """

    nodes = {
        "BF": ("Band Fluxes", "<F<sup>Mod</sup><sub>i</sub>(&theta;)>"),
        "L": ("Instrisinc Spectrum", "<F<sub>&lambda;</sub>(&theta;<sub>star</sub>)>"),
        "D": ("Dust Extinction", "<D<sub>&lambda;</sub>(&theta;<sub>dust</sub>)>"),
        "LD": (
            "Dust Extinguished Spectrum",
            "<F<sup>Mod</sup><sub>&lambda;</sub>(&theta;)>",
        ),
        "F": ("Band Responses", "<B<sub>i</sub>(&lambda;)>"),
        "M": ("mass\nM", "M"),
        "t": ("age\nT", "t"),
        "Z": ("metallicity\nZ", "Z"),
        "d": ("distance\nd", "d"),
        "Av": ("dust column\nA(V)", "A(V)"),
        "Rv": ("grain size\nR(V)", "R(V)"),
        "fA": ("<f<SUB>A</SUB>>", "<f<SUB>A</SUB>>"),
        "IMF": ("Initial\nMass\nFunction", "IMF"),
        "SFH": ("Star\nFormation\nHistory", "SFH"),
        "AMR": ("Age\nMetallicity\nRelation", "AMR"),
        "DP": ("measured\ndistance", "P(d)"),
        "FC": ("Foreground\nDust", "FD"),
        "GC": ("Instrinsic\nDust", "ID"),
        "Obs": (
            "Observation Model",
            "<&mu;<sub>i</sub>(F<sup>Mod</sup><sub>i</sub>), &sigma;<sub>i</sub>(F<sup>Mod</sup><sub>i</sub>)>",
        ),
        "Like": ("Observed Band Fluxes", "Observed Band Fluxes"),
        "AST": ("Artifical Star Tests", "ASTs"),
    }

    edges = {
        "F": "BF",
        "L": "LD",
        "D": "LD",
        "LD": "BF",
        "M": "L",
        "t": "L",
        "Z": "L",
        "d": "L",
        "Av": "D",
        "Rv": "D",
        "fA": "D",
        "IMF": "M",
        "SFH": ("M", "t"),
        "AMR": ("Z", "t"),
        "DP": "d",
        "FC": ("Av", "Rv", "fA"),
        "GC": ("Av", "Rv", "fA"),
        "Obs": "Like",
        "BF": ("Like", "AST"),
        "AST": "Obs",
    }

    beast = create_graphic_model(nodes, edges, gtype)

    beast.render(f"beast-graphic-{type}", format=savefig)


def plot_file_flow(savefig="png"):
    """
    Plot a model showing what files are created at each stage of a BEAST
    production run.

    Parameters
    ----------
    savefig : str
        set to the file extension of desired file to save image of model
    """

    nodes = {
        "sed": ["SEDgrid"],
        "phot1": ["phot_SD0-1"],
        "fake1": ["phot_fake_SD0-1"],
        "obs1": ["obsmodel_SD0-1"],
        "phot1s0": ["phot_SD0-1_sub0"],
        "phot1s1": ["phot_SD0-1_sub1"],
        "phot1s2": ["phot_SD0-1_sub2"],
        "sed1s0": ["SEDgrid_SD0-1_sub0_trim"],
        "sed1s1": ["SEDgrid_SD0-1_sub1_trim"],
        "sed1s2": ["SEDgrid_SD0-1_sub2_trim"],
        "obs1s0": ["obsmodel_SD0-1_sub0_trim"],
        "obs1s1": ["obsmodel_SD0-1_sub1_trim"],
        "obs1s2": ["obsmodel_SD0-1_sub2_trim"],
        "stat1s0": ["stats_SD0-1_sub0"],
        "stat1s1": ["stats_SD0-1_sub1"],
        "stat1s2": ["stats_SD0-1_sub2"],
        "pdf1d1s0": ["pdf1d_SD0-1_sub0"],
        "pdf1d1s1": ["pdf1d_SD0-1_sub1"],
        "pdf1d1s2": ["pdf1d_SD0-1_sub2"],
        "pdf2d1s0": ["pdf2d_SD0-1_sub0"],
        "pdf2d1s1": ["pdf2d_SD0-1_sub1"],
        "pdf2d1s2": ["pdf2d_SD0-1_sub2"],
        "lnp1s0": ["lnp_SD0-1_sub0"],
        "lnp1s1": ["lnp_SD0-1_sub1"],
        "lnp1s2": ["lnp_SD0-1_sub2"],
    }

    edges = {
        "phot1": ("phot1s0", "phot1s1", "phot1s2"),
        "fake1": "obs1",
        "phot1s0": ("sed1s0", "obs1s0"),
        "phot1s1": ("sed1s1", "obs1s1"),
        "phot1s2": ("sed1s2", "obs1s2"),
        "sed": ("sed1s0", "sed1s1", "sed1s2"),
        "obs1": ("obs1s0", "obs1s1", "obs1s2"),
        "sed1s0": ("stat1s0", "pdf1d1s0", "pdf2d1s0", "lnp1s0"),
        "obs1s0": ("stat1s0", "pdf1d1s0", "pdf2d1s0", "lnp1s0"),
        "sed1s1": ("stat1s1", "pdf1d1s1", "pdf2d1s1", "lnp1s1"),
        "obs1s1": ("stat1s1", "pdf1d1s1", "pdf2d1s1", "lnp1s1"),
        "sed1s2": ("stat1s2", "pdf1d1s2", "pdf2d1s2", "lnp1s2"),
        "obs1s2": ("stat1s2", "pdf1d1s2", "pdf2d1s2", "lnp1s2"),
    }

    beast = create_graphic_model(nodes, edges, "text")
    beast.graph_attr['rankdir'] = 'LR'
    beast.render("beast-file-flow", format=savefig)

if __name__ == "__main__":
    #plot_graphic_model("text")
    #plot_graphic_model("math")
    plot_file_flow()
