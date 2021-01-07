"""
Add whichever filter that is not currently available in the BEAST.
"""
from beast.config import __ROOT__
from beast.observationmodel import phot

# 3rd party import
from astropy.io as ascii
import tables
import argparse

__default_filtlist__ = __ROOT__ + "/filters.hd5"
__default_vega__ = __ROOT__ + "vega.hd5"

__all__ = [
    "add_filter",
]


class __newFilterTable__(tables.IsDescription):
    """ define table to store filter dataset """

    WAVELENGTH = tables.FloatCol(pos=0)
    THROUGHPUT = tables.FloatCol(pos=1)


def append_filter(
    throughphutfile,
    tablename,
    filtername,
    observatory=None,
    instrument=None,
    comment=None,
    filterLib=__default_filtlist__,
    updateVegaLib=True,
):
    """
    Add a filter that is not currently available in the BEAST to
    the ./beast/filters.hd5.
    The input ascii file has to have two columns: wavelength (in 
    Angstrom) and corresponding throughput. 

    Parameters
    ----------
    throughputfile: string
        ascii file name of a filter throughput curve
        
    tablename: str
        table name in the library

    filtername: str
        name of the filter (F225W, FUV, ...)

    observatory: str (default=None)
        observatory of the filter (Ground, HST, Spitzer, ...)

    instrument: str (default=None)
        instrument associated with the filter (WFC3, ACS_WFC, ...)

    comment: str, optional (default=None)
        optinal comment to keep with the filter

    filterLib: str, optional (default='~/.beast/filters.hd5')
        filter library file to use

    updateVegaLib: bool (default=True)
        if set calls the update function to the vega library
    """

    # read in the filter curve
    curve = ascii.read(throughputfile)
    wave, thrp = curve[tmp.colnames[0]], curve[tmp.colnames[1]]
 
    # generate filter instance
    newfilt = phot.Filter(wave, thrp, name=filtername)

    # open the BEAST filter library in 'append' mode
    filtlist = tables.open_file(filterLib, mode='r+')

    # add a content for the filter
    contentTab = filtlist.get_node("/content")
    if contentTab.read_where('TABLENAME == "{0}"'.format(tablename)).size > 0:
        print(
            "% {0}: Filter {1} already exists. Returning".format(sys.argv[0], tablename)
        )
        return

    newRow["TABLENAME"] = tablename
    newRow["OBSERVATORY"] = observatory
    newRow["INSTRUMENT"] = instrument
    newRow["NAME"] = newfilt.name
    newRow["NORM"] = newfilt.norm
    newRow["CWAVE"] = newfilt.cl
    newRow["PWAVE"] = newfilt.lpivot
    newRow["COMMENT"] = commnet
    newRow.append()
    contentTab.flush()

    # add a table for the filter
    newTab = filtlist.create_table(
        "/filters",
        tablename,
        __newFilterTable__,
        title=newfilt.name,
        expectedrows=newfilt.wavelength.size,
    )

    newRow = newTab.row
    for i in range(newfilt.wavelength.size):
        newRow["WAVELENGTH"] = newfilt.wavelength[i]
        newRow["THROUGHPUT"] = newfilt.transmit[i]
        newRow.append()
    newTab.flush()
    filtlist.flush()
    filtlist.close()

    print("% {0}: Filter {1} added to {2}".format(sys.argv[0], name, filterLib))

    if updateVegaLib:
        appendVegaFilter(newfilt)


def appendVegaFilter(filtInst, VegaLib=__default_vega__):
    """
    Add filter properties to the Vega library

    Parameters
    ----------
    filtInst: filter instance
        filter instance to get properties from and store information with Vega

    VegaLib: str
        Vega Library
    """

    # read the Vega flux library in 'append' mode.
    vtab = tables.open_file(VegaLib, "r+")
    vl = vtab.root.spectrum[:]["WAVELENGTH"]
    vf = vtab.root.spectrum[:]["FLUX"]

    # check if Vega flux information already exists for the filter 
    sedTab = vtab.get_node("/sed")
    if sedTab.read_where('FNAME == "{0}"'.format(filtInst.name)).size > 0:
        print(
            "% {0}: Filter {1} already exists. Returning".format(
                sys.argv[0], filtInst.name
            )
        )
        return

    # get Vega flux and mag
    flux = filtInst.getFlux(vl, vf)
    mag = -2.5 * np.log10(flux)

    # append a new Tab
    sedTab.append([(
        filtInst.name,
        filtInst.cl,
        flux,
        mag)
        ]
    )
    
    sedTab.flush()
    vtab.close()

    print("% {0}: Filter {1} added to {2}".format(sys.argv[0], filtInst.name, VegaLib))

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("thrpfile", help="filename of filter throughput")
    parser.add_argument(
        "--observatory",
        type=str,
        default=None,
        help="If set, save the Observatory",
    )
