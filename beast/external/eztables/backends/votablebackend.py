"""
Simplistic VOTABLE XML Parser.

This package aims to parse VOTABLE independently from any non-standard package.
It uses the python xml package and numpy to generate a named array with meta
descriptions.

Currently very little verifications are applied to the xml structure.

VO Format descriptions
.. _BINARY: http://www.ivoa.net/Documents/PR/VOTable/VOTable-20040322.html#ToC27
.. _COOSYS: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC19
.. _DESCRIPTION: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC19
.. _FIELD: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC24
.. _FIELDref: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC31
.. _FITS: http://fits.gsfc.nasa.gov/fits_documentation.html
.. _GROUP: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC31
.. _ID: http://www.w3.org/TR/REC-xml/#id
.. _INFO: http://www.ivoa.net/Documents/VOTable/20040811/REC-VOTable-1.1-20040811.html#ToC19
.. _LINK: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC22
.. _multidimensional arrays: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC12
.. _numerical accuracy: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC26
.. _PARAM: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC24
.. _PARAMref: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC31
.. _RESOURCE: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC21
.. _TABLE: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC23
.. _TABLEDATA: http://www.ivoa.net/Documents/PR/VOTable/VOTable-20040322.html#ToC25
.. _unified content descriptor: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC28
.. _unique type: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC29
.. _units: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC27
.. _VALUES: http://www.ivoa.net/Documents/REC/VOTable/VOTable-20040811.html#ToC30
.. _VOTABLE: http://www.ivoa.net/Documents/PR/VOTable/VOTable-20040322.html#ToC9
"""

from xml.dom import minidom
import numpy as np


class Node(object):
    def __init__(self, node):
        self.value = node.nodeValue
        if node.attributes is not None:
            self.attributes = {}
            for k, v in list(node.attributes.items()):
                self.attributes[k] = v
        else:
            self.attributes = None
        self.__getChilds__(node)
        self.source = node

    def __getChilds__(self, node):
        if not node.hasChildNodes():
            return

        for nk in node.childNodes:
            if nk.nodeName not in ['#text']:
                suf = 1
                key = '%s' % (nk.nodeName)
                while key in list(self.__dict__.keys()):
                    key = '%s_%d' % (nk.nodeName, suf)
                    suf += 1
                setattr( self, key, Node(nk) )
            else:
                txt = nk.nodeValue
                if txt not in ['', '\n', '\n\n']:
                    setattr( self, 'text', nk.nodeValue )

    @property
    def childrenNames(self):
        return [ k for k, v in list(self.__dict__.items()) if isinstance(v, Node) ]

    @property
    def children(self):
        return [ v for k, v in list(self.__dict__.items()) if isinstance(v, Node) ]

    def items(self):
        return [ (k, v) for k, v in list(self.__dict__.items()) if isinstance(v, Node) ]

    def __repr__(self):
        txt = 'XML Node: %s \n' % object.__repr__(self)
        txt += 'children:   %s \n' % ', '.join(self.childrenNames)
        if (self.attributes is not None):
            if (len(self.attributes) > 0):
                txt += 'attributes:\n'
                for k, v in list(self.attributes.items()):
                    txt += '   %s: %s\n' % (k, v)
        if self.value is not None:
            txt += 'value:\n  %s' % str(self.value)
        return txt

    def tree(self, level=0):
        """ return ascii tree representation of the nodes
        KEYWORDS:
            level   int     defines the initial indentation
        """
        txt = ''
        if level > 0:
            indent = ' ' * 3 * ( level - 1 ) + ' +--'
        else:
            indent = ' ' * 2 * (level)
        for name, node in list(self.items()):
            txt += '{}{}\n'.format(indent, name)
            txt += '{}'.format(node.tree(level + 1))
        return txt


def genericMissing(_type, val=0, repr='nan'):
    """
    This function create a class of missing data respecting the requested types
    Instances will act as traditional np.nan (in particular when using int)

    INPUT:
        _type   type    type for the value to generate
    OUTPUT:
        v      NA       missing value instance of type _type
    KEYWORD:
        val     type()  a value by default to use for numeric conversion
        repr    str     the representation of the value
    """
    try:
        return _type('nan')
    except TypeError:
        return type(_type)('nan')
    except ValueError:

        class NA(_type):
            """Subclass the expected type to generate a value acting as a NaN
                Any operation will return the same object without any change.
            """

            def __repr__(self):
                return repr

            def __cmp__(self, other):
                return False

            def __add__(self, other):
                return self

            def __mul__(self, other):
                return self

            def __div__(self, other):
                return self

            def __sub__(self, other):
                return self

            def __pow__(self, other):
                return self

            def __float__(self):
                return np.nan

            def __str__(self):
                return str(np.nan)

        return NA(val)


def parse(str):
    """ Return a Node-tree from a xml str """
    return Node( minidom.parseString(str) )


def votable(obj):
    """ Return a named numpy array (equivalent to a recarray) as well as header
    information and meta-descriptions of the columns
    INPUT:
        obj     str     if a string is provided, parse will be called
                Node    if obj is a Node object, it will assume a result from parse
    OUPUTS:
        rtab    ndarray a numpy named array
        header  dict    dictionary containing votable descriptions
        meta    tuple   (ucds, units, description) of individual columns
    """
    if type(obj) in [str, np.str, np.unicode0, np.unicode_, str, np.string_, np.str_, np.string0]:
        node = parse(obj)
    else:
        node = obj
    assert('VOTABLE' in node.childrenNames), 'Expected XML tree root named VOTABLE'

    header = {}
    for k, v in list(node.VOTABLE.attributes.items()):
        header[k] = v

    for k, v in list(node.VOTABLE.RESOURCE.attributes.items()):
        header[k] = v

    if 'DESCRIPTION' in node.VOTABLE.RESOURCE.TABLE.childrenNames:
        header['DESCRIPTION'] = node.VOTABLE.RESOURCE.TABLE.DESCRIPTION.text

    cols = [ k for k in node.VOTABLE.RESOURCE.TABLE.childrenNames if k[:5] == 'FIELD' ]
    idx  = list(range(len(cols)))
    for k, v in list(node.VOTABLE.RESOURCE.TABLE.items()):
        if len(k.split('_')) > 1:
            idx = int(k.split('_')[-1])
        else:
            idx = 0
        cols[idx] = v

    colnames  = [ str(v.attributes.get('ID')) for v in cols ]
    coltypes  = [ str(v.attributes.get('datatype')) for v in cols ]
    colwidths = [ int(v.attributes.get('width')) for v in cols ]
    colucds   = [ str(v.attributes.get('ucd')) for v in cols ]
    colunits  = [ str(v.attributes.get('unit')) for v in cols ]
    coldescs  = [ v.DESCRIPTION.text for v in cols if 'DESCRIPTION' in v.childrenNames ]
    for e, ktype in enumerate(coltypes):
            #the width on decimals needs to be odd, and does not matter on strings
            if ktype in ['char', 'str', 'string']:
                coltypes[e] = '|S%d' % colwidths[e]
            else:
                coltypes[e] = np.typeDict[ktype]

    data  = node.VOTABLE.RESOURCE.TABLE.DATA.TABLEDATA
    dtype = np.dtype( [ (colnames[e], ktype) for e, ktype in enumerate(coltypes) ])
    rtab = np.zeros( len(data.childrenNames), dtype=dtype )

    for nk, line in list(data.items()):
        if nk[:2] == 'TR':
            lidx = 0
            if len(nk.split('_')) > 1:
                lidx = int(nk.split('_')[-1])

            for lk, v in list(line.items()):
                if lk[:2] == 'TD':
                    cidx = 0
                    if len(lk.split('_')) > 1:
                        cidx = int(lk.split('_')[-1])
                    if hasattr(v, 'text'):
                        rtab[lidx][cidx] = v.text
                    else:
                        rtab[lidx][cidx] = genericMissing(coltypes[cidx])

    if 'DEFINITIONS' in node.VOTABLE.childrenNames:
        if 'COOSYS' in node.VOTABLE.DEFINITIONS.childrenNames:
            header['EPOCH'] = node.VOTABLE.DEFINITIONS.COOSYS.attributes['epoch']
            header['COORDS'] = node.VOTABLE.DEFINITIONS.COOSYS.attributes['system']

    return rtab, header, (colucds, colunits, coldescs)
