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


if __name__ == "__main__":
    plot_graphic_model("text")
    plot_graphic_model("math")
