""" FITS Backend
    read/write handles: units, column comments, aliases, header keywords
"""

import os
import inspect
localpath = '/'.join(os.path.abspath(inspect.getfile(inspect.currentframe())).split('/')[:-1])
import astropy.io.fits as pyfits
import numpy as np
from .basebackend import BaseBackend
from ..core.tableheader import TableHeader
from ..table import Table

# mapping from TFORM data type to numpy data type (code)
# L: Logical (Boolean)
# B: Unsigned Byte
# I: 16-bit Integer
# J: 32-bit Integer
# K: 64-bit Integer
# E: Single-precision Floating Point
# D: Double-precision Floating Point
# C: Single-precision Complex
# M: Double-precision Complex
# A: Character
fitstypeDict = {}
fitstypeDict[np.bool_]    = 'B'
fitstypeDict[np.int8]     = 'I'
fitstypeDict[np.uint8]    = 'I'
fitstypeDict[np.int16]    = 'I'
fitstypeDict[np.uint16]   = 'J'
fitstypeDict[np.int32]    = 'J'
fitstypeDict[np.uint32]   = 'K'
fitstypeDict[np.int64]    = 'K'
fitstypeDict[np.uint64]   = 'E'
fitstypeDict[np.float16]  = 'E'
fitstypeDict[np.float32]  = 'F'
fitstypeDict[np.float64]  = 'D'
if hasattr(np, 'float128'):
    fitstypeDict[np.float128] = 'D'
fitstypeDict[np.complex]  = 'M'
fitstypeDict[np.str]      = 'A'
fitstypeDict[np.string_]  = 'A'
fitstypeDict[str]         = 'A'
fitstypeDict[np.unicode_] = 'A'


def __getFitsFmt__(col):
    try:
        ft = fitstypeDict[col.dtype.type]
    except KeyError:
        raise Exception("Type conversion not found for %s " % col.dtype.type )
    if ft == 'A':
        _len = np.max([len(k) for k in col])
        return(str(_len) + 'A')
    return ft


#==============================================================================
class fitsBackend(BaseBackend):
#==============================================================================
    def __init__(self):
        """ constructor """
        BaseBackend.__init__(self, tableType='fits')

    def readHeader(self, hdu):
        header = TableHeader()
        alias = []
        genTerms = [ 'XTENSION', 'BITPIX', 'NAXIS', 'NAXIS1',
                 'NAXIS2', 'PCOUNT', 'GCOUNT', 'TFIELDS'  ]
        fieldTerms = ['TTYPE', 'TFORM', 'TUNIT', 'ALIAS']
        for k in hdu.header:
            if (not k in genTerms) & (not k[:5] in fieldTerms):
                header[k] = hdu.header[k]
            elif (k[:5] == 'ALIAS'):
                alias.append( hdu.header[k].split('=') )
        if 'EXTNAME' in hdu.header:
            header['NAME'] = hdu.header['EXTNAME']
        return header, alias

    def _readColComments(self, card):
        ttype = str(card).split('/')
        if len(ttype) > 1:
            return ' '.join(ttype[-1].split())
        else:
            return ''

    def readData(self, hdu):
        colDef = hdu.columns
        names = [ k.name for k in colDef ]
        units = [ k.unit for k in colDef ]
        #comms = [ k.comment for k in colDef ]
        comms = [ self._readColComments(k) for k in hdu.header.cards['TTYPE*'] ]
        data  = hdu.data.view(np.rec.recarray)
        return data, names, units, comms

    def read(self, filename, extension=1, **kwargs):
        #TODO: read from buffer
        hdu = pyfits.open(filename)
        header, alias = self.readHeader(hdu[extension])
        data, names, units, comm = self.readData(hdu[extension])
        hdu.close()
        tab = Table()
        for k, v in header.items():
            tab.header[k] = v
        for i, k in enumerate(names):
            tab.add_column(k, data[k], description=comm[i], dtype=data.dtype[i])
        for k in alias:
            tab.set_alias(k[0], k[1])
        return tab

    def writeColComment(self, header, colName, comment):
        cards = header.ascard['TTYPE*']
        refs = {}
        for k in cards:
            refs[k.value] = k.key
        header.update(refs[colName], colName, comment=comment)

    def write(self, tab, output='exportedData.fits', fmt=None,
              name=None, comments=None, units=None, clobber=False, append=True,
              global_attrs=None, silent=False, hdr_only=False, keep_open=False, hdr0=None):
        """
        export data to a FITS file

        inputs:
            data -- data dictionnary to export

        outputs:
            output -- output file (def: exportedData.dat)

        keywords:
            fmt     -- force a given column format
                     (alphabetic data name order)
            name    -- extention name of the fits file (def: DATA)
            header  -- dictonnary of keywords and corresponding
                   values
            comments-- list of column comments
            units   -- list of column units
            global_attrs   -- dictonnary of keywords to add to the
                      main container. (or update; No need
                      to include other keywords to do so.)
            clobber -- overwrite if set (def: False)
            append  -- add data to existing file (def: True)
        """
        #if hasattr(output, 'write'):
        #    unit = output
        #else:
        #    unit = pyfits.open(output, 'w')
        #    if (not os.path.isfile(output)) & append:
        #        if not silent:
        #            print "Warning: %s does not seem to exist" % output
        #            print "         Creating a new file."
        #        append = False

        keys = list(tab.keys())

        #get formats, units, descriptions
        fmt = [ __getFitsFmt__(tab[k]) for k in keys ]
        comments = [None] * tab.ncols
        units = [None] * tab.ncols
        for i, k in enumerate(keys):
            comm = tab.columns[k].description
            _unit = tab.columns[k].unit
            if comm not in ['', None, 'None']:
                comments[i] = comm
            if _unit not in ['', None, 'None']:
                units[i] = _unit

        #generate FITS table from a set of columns
        cols = [ pyfits.Column( name=k, array=tab[k],   format=fmt[i], unit=units[i] ) for i, k in enumerate(keys)  ]

        hdr = pyfits.new_table(cols)

        #Update table header
        if tab.header['NAME'] not in ['', 'None', None]:
            hdr.header.update( 'EXTNAME', tab.header['NAME'] )
        if not hasattr(output, 'write'):
            hdr.header.update( 'FILENAME', output.split('/')[-1] )

        for k, v in tab.header.items():
            if v not in ['', 'None', None]:
                if (k != 'COMMENT') & (k != 'HISTORY'):
                    hdr.header.update(k, v)
                else:
                    txt = v.split('\n')
                    for j in txt:
                        if k == 'COMMENT':
                            hdr.header.add_comment(j)
                        elif k == 'HISTORY':
                            hdr.header.add_history(j)
        for i, k in enumerate(keys):
            if comments[i] not in  ['', 'None', None]:
                self.writeColComment(hdr.header, k, comments[i])

        for k, v in enumerate(tab._aliases.items()):
            hdr.header.update( 'ALIAS%d' % k, '='.join(v) )

        if hdr_only:
            return hdr

        if hasattr(output, 'write'):
            hdr.writeto(output)
        elif not append:
            hdr.writeto(output, clobber=clobber)
            if not silent:
                print("Data exported into %s" % output)
        else:
            if not hasattr(output, 'write'):
                if (hdr0 is None):
                    retHdr = True
                    hdr0 = pyfits.open(output, mode='append')
                else:
                    retHdr = False

                hdr0.append(hdr)
                if not keep_open:
                    hdr0.flush()
                    hdr0.close()
                else:
                    if retHdr:
                        return hdr0
                if not silent:
                    print("Data added into %s" % output)
        if not hasattr(output, 'write'):
            if global_attrs is not None:
                if hdr0 is None:
                    hdr0 = pyfits.open(output, mode='update')
                for k in global_attrs:
                    hdr0[0].header.update(k, global_attrs[k])
                hdr0.flush()
                hdr0.close()
                if not silent:
                    print("Keywords added to main table into %s" % output)
