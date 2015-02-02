""" This file implements a Table class
    that is designed to be the basis of any format

Requirements
------------

* FIT format:
    * astropy:
        provides a replacement to pyfits
        pyfits can still be used instead but astropy is now the default

* HDF5 format:
    * pytables

RuntimeError will be raised when writing to a format associated with missing
package.


.. code-block::python

    >>> t = SimpleTable('path/mytable.csv')
    # get a subset of columns only
    >>> s = t.get('M_* logTe logLo U B V I J K')
    # set some aliases
    >>> t.set_alias('logT', 'logTe')
    >>> t.set_alias('logL', 'logLLo')
    # make a query on one or multiple column
    >>> q = s.selectWhere('logT logL', '(J > 2) & (10 ** logT > 5000)')
    # q is also a table object
    >>> q.plot('logT', 'logL', ',')
    # makes a simple plot
    >>> s.write('newtable.fits')
    # export the initial subtable to a new file
"""
from __future__ import (absolute_import, division, print_function)

__version__ = '3.0'

import sys
import numpy as np
from numpy.lib import recfunctions
from copy import deepcopy
import re
import itertools

try:
    from astropy.io import fits as pyfits
except ImportError:
    import pyfits
except:
    pyfits = None

try:
    import tables
except ImportError:
    tables = None

# ==============================================================================
# Python 3 compatibility behavior
# ==============================================================================
# remap some python 2 built-ins on to py3k behavior or equivalent
# Most of them become generators
import operator

PY3 = sys.version_info[0] > 2

if PY3:
    iteritems = operator.methodcaller('items')
    itervalues = operator.methodcaller('values')
    basestring = (str, bytes)
else:
    range = xrange
    from itertools import izip as zip
    iteritems = operator.methodcaller('iteritems')
    itervalues = operator.methodcaller('itervalues')
    basestring = (str, unicode)


# ==============================================================================
# Specials -- special functions
# ==============================================================================

def pretty_size_print(num_bytes):
    """
    Output number of bytes in a human readable format

    keywords
    --------
    num_bytes: int
        number of bytes to convert

    returns
    -------
    output: str
        string representation of the size with appropriate unit scale
    """
    if num_bytes is None:
        return

    KiB = 1024
    MiB = KiB * KiB
    GiB = KiB * MiB
    TiB = KiB * GiB
    PiB = KiB * TiB
    EiB = KiB * PiB
    ZiB = KiB * EiB
    YiB = KiB * ZiB

    if num_bytes > YiB:
        output = '%.3g YB' % (num_bytes / YiB)
    elif num_bytes > ZiB:
        output = '%.3g ZB' % (num_bytes / ZiB)
    elif num_bytes > EiB:
        output = '%.3g EB' % (num_bytes / EiB)
    elif num_bytes > PiB:
        output = '%.3g PB' % (num_bytes / PiB)
    elif num_bytes > TiB:
        output = '%.3g TB' % (num_bytes / TiB)
    elif num_bytes > GiB:
        output = '%.3g GB' % (num_bytes / GiB)
    elif num_bytes > MiB:
        output = '%.3g MB' % (num_bytes / MiB)
    elif num_bytes > KiB:
        output = '%.3g KB' % (num_bytes / KiB)
    else:
        output = '%.3g Bytes' % (num_bytes)

    return output


def _fits_read_header(hdr):
    """
    Convert pyfits header into dictionary with relevant values

    Parameters
    ----------

    hdr: pyftis.Header
        fits unit

    Returns
    -------
    header: dict
        header dictionary

    alias: dict
        aliases

    units: dict
        units

    comments: dict
        comments/description of keywords
    """
    header = {}
    alias = {}
    units = {}
    comments = {}

    # generic cards
    genTerms = ['XTENSION', 'BITPIX', 'NAXIS', 'NAXIS1',
                'NAXIS2', 'PCOUNT', 'GCOUNT', 'TFIELDS',
                'EXTNAME']
    fieldTerms = ['TTYPE', 'TFORM', 'TUNIT', 'ALIAS']

    # read col comments
    for k, name, comment in hdr.ascard['TTYPE*']:
        comments[name] = comment
        u = hdr.get(k.replace('TYPE', 'UNIT'), None)
        if u is not None:
            units[name] = u

    for k, val, _ in hdr.ascard['ALIAS*']:
        al, orig = val.split('=')
        alias[al] = orig

    # other specific keywords: COMMENT, HISTORY
    header_comments = []
    header_history = []
    for k, v in hdr.items():
        if (k not in genTerms) and (k[:5] not in fieldTerms):
            if (k == 'COMMENT'):
                header_comments.append(v)
            elif (k == 'HISTORY'):
                header_history.append(v)
            else:
                header[k] = v

    # COMMENT, HISTORY polish
    if len(header_comments) > 0:
        header['COMMENT'] = '\n'.join(header_comments)
    if len(header_history) > 0:
        header['HISTORY'] = '\n'.join(header_history)

    if 'EXTNAME' in hdr:
        header['NAME'] = hdr['EXTNAME']

    return header, alias, units, comments


def _fits_generate_header(tab):
    """ Generate the corresponding fits Header that contains all necessary info

    Parameters
    ----------

    tab: SimpleTable instance
        table

    Returns
    -------
    hdr: pyfits.Header
        header instance
    """
    # get column cards

    cards = []

    # names units and comments
    for e, k in enumerate(tab.keys()):
        cards.append(('TTYPE{0:d}'.format(e + 1), k, tab._desc.get(k, '')))
        u = tab._units.get(k, '')
        if u not in ['', 'None', None]:
            cards.append(('TUNIT{0:d}'.format(e + 1), tab._units.get(k, ''),
                          'unit of {0:s}'.format(k)))

    # add aliases
    for e, v in enumerate(tab._aliases.items()):
        cards.append( ('ALIAS{0:d}'.format(e + 1), '='.join(v), '') )

    if tab.header['NAME'] not in ['', 'None', None, 'No Name']:
        cards.append(('EXTNAME', tab.header['NAME'], ''))

    hdr = pyfits.Header(cards)

    for k, v in tab.header.items():
        if (v not in ['', 'None', None]) & (k != 'NAME'):
            if (k != 'COMMENT') & (k != 'HISTORY'):
                hdr.update(k, v)
            else:
                txt = v.split('\n')
                for j in txt:
                    if k == 'COMMENT':
                        hdr.add_comment(j)
                    elif k == 'HISTORY':
                        hdr.add_history(j)
    return hdr


def _fits_writeto(filename, data, header=None, output_verify='exception',
                  clobber=False, checksum=False):
    """
    Create a new FITS file using the supplied data/header.
    Patched version of pyfits to correctly include provided header

    Parameters
    ----------
    filename : file path, file object, or file like object
        File to write to.  If opened, must be opened in a writeable binary
        mode such as 'wb' or 'ab+'.

    data : array, record array, or groups data object
        data to write to the new file

    header : `Header` object, optional
        the header associated with ``data``. If `None`, a header
        of the appropriate type is created for the supplied data. This
        argument is optional.

    output_verify : str
        Output verification option.  Must be one of ``"fix"``, ``"silentfix"``,
        ``"ignore"``, ``"warn"``, or ``"exception"``.  May also be any
        combination of ``"fix"`` or ``"silentfix"`` with ``"+ignore"``,
        ``+warn``, or ``+exception" (e.g. ``"fix+warn"``).  See :ref:`verify`
        for more info.

    clobber : bool, optional
        If `True`, and if filename already exists, it will overwrite
        the file.  Default is `False`.

    checksum : bool, optional
        If `True`, adds both ``DATASUM`` and ``CHECKSUM`` cards to the
        headers of all HDU's written to the file
    """

    hdu = pyfits.convenience._makehdu(data, header)
    hdu.header.update(header.cards)
    if hdu.is_image and not isinstance(hdu, pyfits.PrimaryHDU):
        hdu = pyfits.PrimaryHDU(data, header=header)
    hdu.writeto(filename, clobber=clobber, output_verify=output_verify,
                checksum=checksum)


def _fits_append(filename, data, header=None, checksum=False, verify=True,
                 **kwargs):
    """
    Append the header/data to FITS file if filename exists, create if not.

    If only ``data`` is supplied, a minimal header is created.
    Patched version of pyfits to correctly include provided header

    Parameters
    ----------
    filename : file path, file object, or file like object
        File to write to.  If opened, must be opened for update (rb+) unless it
        is a new file, then it must be opened for append (ab+).  A file or
        `~gzip.GzipFile` object opened for update will be closed after return.

    data : array, table, or group data object
        the new data used for appending

    header : `Header` object, optional
        The header associated with ``data``.  If `None`, an appropriate header
        will be created for the data object supplied.

    checksum : bool, optional
        When `True` adds both ``DATASUM`` and ``CHECKSUM`` cards to the header
        of the HDU when written to the file.

    verify : bool, optional
        When `True`, the existing FITS file will be read in to verify it for
        correctness before appending.  When `False`, content is simply appended
        to the end of the file.  Setting ``verify`` to `False` can be much
        faster.

    kwargs
        Any additional keyword arguments to be passed to
        `astropy.io.fits.open`.
    """

    name, closed, noexist_or_empty = pyfits.convenience._stat_filename_or_fileobj(filename)

    if noexist_or_empty:
        #
        # The input file or file like object either doesn't exits or is
        # empty.  Use the writeto convenience function to write the
        # output to the empty object.
        #
        _fits_writeto(filename, data, header, checksum=checksum, **kwargs)
    else:
        hdu = pyfits.convenience._makehdu(data, header)
        hdu.header.update(header.cards)

        if isinstance(hdu, pyfits.PrimaryHDU):
            hdu = pyfits.ImageHDU(data, header)

        if verify or not closed:
            f = pyfits.convenience.fitsopen(filename, mode='append')
            f.append(hdu)

            # Set a flag in the HDU so that only this HDU gets a checksum when
            # writing the file.
            hdu._output_checksum = checksum
            f.close(closed=closed)
        else:
            f = pyfits.convenience._File(filename, mode='append')
            hdu._output_checksum = checksum
            hdu._writeto(f)
            f.close()


def _ascii_read_header(fname, comments='#', delimiter=None, commentedHeader=True,
                       *args, **kwargs):
    """
    Read ASCII/CSV header

    Parameters
    ----------
    fname: str or stream
        File, filename, or generator to read.
        Note that generators should return byte strings for Python 3k.

    comments: str, optional
        The character used to indicate the start of a comment;
        default: '#'.

    delimiter: str, optional
        The string used to separate values.  By default, this is any
        whitespace.

    commentedHeader: bool, optional
        if set, the last line of the header is expected to be the column titles

    Returns
    -------
    nlines: int
        number of lines from the header

    header: dict
        header dictionary

    alias: dict
        aliases

    units: dict
        units

    comments: dict
        comments/description of keywords

    names: sequence
        sequence or str, first data line after header, expected to be the column
        names.
    """
    if hasattr(fname, 'read'):
        stream = fname
    else:
        stream = open(fname, 'r')

    header = {}
    alias = {}
    units = {}
    desc = {}

    def parseStrNone(v):
        """ robust parse """
        _v = v.split()
        if (len(_v) == 0):
            return None
        else:
            _v = ' '.join(_v)
            if (_v.lower()) == 'none' or (_v.lower() == 'null'):
                return None
            else:
                return _v

    done = False
    oldline = None
    lasthdr = None
    nlines = 0
    header.setdefault('COMMENT', '')
    header.setdefault('HISTORY', '')
    while done is False:
        line = stream.readline()[:-1]  # getting rid of '\n'
        nlines += 1
        if (line[0] == comments):  # header part
            if (len(line) > 2):
                if line[1] == comments:  # column meta data
                    # column meta is expected to start with ##
                    k = line[2:].split('\t')
                    colname = k[0].strip()
                    colunit = None
                    colcomm = None
                    if len(k) > 1:
                        colunit = parseStrNone(k[1])
                    if len(k) > 2:
                        colcomm = parseStrNone(k[2])

                    if colunit is not None:
                        units[colname] = colunit
                    if colcomm is not None:
                        desc[colname] = colcomm
                else:
                    # header is expected as "# key \t value"
                    k = line[1:].split('\t')
                    if len(k) > 1:
                        key = k[0].strip()  # remove trainling spaces
                        val = ' '.join(k[1:]).strip()

                        if key in ('', None, 'None', 'NONE', 'COMMENT'):
                            header['COMMENT'] = header['COMMENT'] + '\n' + val
                        if key in ('HISTORY', ):
                            header['HISTORY'] = header['HISTORY'] + '\n' + val
                        elif 'alias' in key.lower():
                            # take care of aliases
                            al, orig = val.split('=')
                            alias[al] = orig
                        else:
                            header[key] = val
                        lasthdr = key
                    else:
                        header['COMMENT'] = header['COMMENT'] + '\n' + line[1:]
        else:
            done = True
            if commentedHeader and (oldline is not None):
                names = oldline.split(delimiter)
                nlines -= 1
                if lasthdr == names[0]:
                    header.pop(lasthdr)
            else:
                names = line.split(delimiter)
        oldline = line[1:]

    if not hasattr(fname, 'read'):
        stream.close()
    else:
        stream.seek(stream.tell() - len(line))
        nlines = 0  # make sure the value is set to the current position

    return nlines, header, units, desc, alias, names


def _hdf5_write_data(filename, data, tablename=None, mode='w', append=False,
                     header={}, units={}, comments={}, aliases={}, **kwargs):
    """ Write table into HDF format

    Parameters
    ----------
    filename : file path, or tables.File instance
        File to write to.  If opened, must be opened and writable (mode='w' or 'a')

    data: recarray
        data to write to the new file

    tablename: str
        path of the node including table's name

    mode: str
        in ('w', 'a') mode to open the file

    append: bool
        if set, tends to append data to an existing table

    header: dict
        table header

    units: dict
        dictionary of units

    alias: dict
        aliases

    comments: dict
        comments/description of keywords

    .. note::
        other keywords are forwarded to :func:`tables.openFile`
    """

    if hasattr(filename, 'read'):
        raise Exception("HDF backend does not implement stream")

    if append is True:
        mode = 'a'

    if isinstance(filename, tables.File):
        if (filename.mode != mode) & (mode != 'r'):
            raise tables.FileModeError('The file is already opened in a different mode')
        hd5 = filename
    else:
        hd5 = tables.openFile(filename, mode=mode)

    # check table name and path
    tablename = tablename or header.get('NAME', None)
    if tablename in ('', None, 'Noname', 'None'):
        tablename = '/data'

    w = tablename.split('/')
    where = '/'.join(w[:-1])
    name = w[-1]
    if where in ('', None):
        where = '/'
    if where[0] != '/':
        where = '/' + where

    if append:
        try:
            t = hd5.getNode(where + name)
            t.append(data.astype(t.description._v_dtype))
            t.flush()
        except tables.NoSuchNodeError:
            print(("Warning: Table {0} does not exists.  \n A new table will be created").format(where + name))
            append = False

    if not append:
        t = hd5.createTable(where, name, data, **kwargs)

        # update header
        for k, v in header.items():
            if (k == 'FILTERS') & (float(t.attrs['VERSION']) >= 2.0):
                t.attrs[k.lower()] = v
            else:
                t.attrs[k] = v
        if 'TITLE' not in header:
            t.attrs['TITLE'] = name

        # add column descriptions and units
        for e, colname in enumerate(data.dtype.names):
            _u = units.get(colname, None)
            _d = comments.get(colname, None)
            if _u is not None:
                t.attrs['FIELD_{0:d}_UNIT'] = _u
            if _d is not None:
                t.attrs['FIELD_{0:d}_DESC'] = _d

        # add aliases
        for i, (k, v) in enumerate(aliases.items()):
            t.attrs['ALIAS{0:d}'.format(i)] = '{0:s}={1:s}'.format(k, v)

        t.flush()

    if not isinstance(filename, tables.File):
        hd5.flush()
        hd5.close()


def _hdf5_read_data(filename, tablename=None, silent=False, *args, **kwargs):
    """ Generate the corresponding ascii Header that contains all necessary info

    Parameters
    ----------
    filename: str
        file to read from

    tablename: str
        node containing the table

    silent: bool
        skip verbose messages

    Returns
    -------
    hdr: str
        string that will be be written at the beginning of the file
    """
    source = tables.openFile(filename, *args, **kwargs)

    if tablename is None:
        node = source.listNodes('/')[0]
        tablename = node.name
    else:
        if tablename[0] != '/':
            node = source.getNode('/' + tablename)
        else:
            node = source.getNode(tablename)
    if not silent:
        print("\tLoading table: {0}".format(tablename))

    hdr = {}
    aliases = {}

    # read header
    exclude = ['NROWS', 'VERSION', 'CLASS', 'EXTNAME', 'TITLE']
    for k in node.attrs._v_attrnames:
        if (k not in exclude):
            if (k[:5] != 'FIELD') & (k[:5] != 'ALIAS'):
                hdr[k] = node.attrs[k]
            elif k[:5] == 'ALIAS':
                c0, c1 = node.attrs[k].split('=')
                aliases[c0] = c1

    empty_name = ['', 'None', 'Noname', None]
    if node.attrs['TITLE'] not in empty_name:
        hdr['NAME'] = node.attrs['TITLE']
    else:
        hdr['NAME'] = '{0:s}/{1:s}'.format(filename, node.name)

    # read column meta
    units = {}
    desc = {}

    for (k, colname) in enumerate(node.colnames):
        _u = getattr(node.attrs, 'FIELD_{0:d}_UNIT'.format(k), None)
        _d = getattr(node.attrs, 'FIELD_{0:d}_DESC'.format(k), None)
        if _u is not None:
            units[colname] = _u
        if _d is not None:
            desc[colname] = _d

    data = node[:]

    source.close()

    return hdr, aliases, units, desc, data


def _ascii_generate_header(tab, comments='#', delimiter=' ',
                           commentedHeader=True):
    """ Generate the corresponding ascii Header that contains all necessary info

    Parameters
    ----------

    tab: SimpleTable instance
        table

    comments: str
        string to prepend header lines

    delimiter: str, optional
        The string used to separate values.  By default, this is any
        whitespace.

    commentedHeader: bool, optional
        if set, the last line of the header is expected to be the column titles

    Returns
    -------
    hdr: str
        string that will be be written at the beginning of the file
    """
    hdr = []

    if comments is None:
        comments = ''

    # table header
    length = max(map(len, tab.header.keys()))
    fmt = '{{0:s}} {{1:{0:d}s}}\t{{2:s}}'.format(length)
    for k, v in tab.header.items():
        for vk in v.split('\n'):
            if len(vk) > 0:
                hdr.append(fmt.format(comments, k.upper(), vk.strip()))

    # column metadata
    hdr.append(comments)  # add empty line
    length = max(map(len, tab.keys()))
    fmt = '{{0:s}}{{0:s}} {{1:{0:d}s}}\t{{2:s}}\t{{3:s}}'.format(length)
    for colname in tab.keys():
        unit = tab._units.get(colname, 'None')
        desc = tab._desc.get(colname, 'None')
        hdr.append(fmt.format(comments, colname, unit, desc))

    # aliases
    if len(tab._aliases) > 0:
        hdr.append(comments)  # add empty line
        for k, v in tab._aliases.items():
            hdr.append('{0:s} alias\t{1:s}={2:s}'.format(comments, k, v))

    # column names
    hdr.append(comments)
    if commentedHeader:
        hdr.append('{0:s} {1:s}'.format(comments, delimiter.join(tab.keys())))
    else:
        hdr.append('{0:s}'.format(delimiter.join(tab.keys())))

    return '\n'.join(hdr)


def _latex_writeto(filename, tab, comments='%'):
    """ Write the data into a latex table format

    Parameters
    ----------
    filename: str
        file or unit to write into

    tab: SimpleTable instance
        table

    comments: str
        string to prepend header lines

    delimiter: str, optional
        The string used to separate values.  By default, this is any
        whitespace.

    commentedHeader: bool, optional
        if set, the last line of the header is expected to be the column titles
    """
    txt = "\\begin{table}\n\\begin{center}\n"

    # add caption
    tabname = tab.header.get('NAME', None)
    if tabname not in ['', None, 'None']:
        txt += "\\caption{{{0:s}}}\n".format(tabname)

    # tabular
    txt += '\\begin{{tabular}}{{{0:s}}}\n'.format('c' * tab.ncols)
    txt += tab.pprint(delim=' & ', fields='MAG*', headerChar='', endline='\\\\\n', all=True, ret=True)
    txt += '\\end{tabular}\n'

    # end table
    txt += "\\end{center}\n"

    # add notes if any
    if len(tab._desc) > 0:
        txt += '\% notes \n\\begin{scriptsize}\n'
        for e, (k, v) in enumerate(tab._desc.items()):
            if v not in (None, 'None', 'none', ''):
                txt += '{0:d} {1:s}: {2:s} \\\\\n'.format(e, k, v)
        txt += '\\end{scriptsize}\n'
    txt += "\\end{table}\n"
    if hasattr(filename, 'write'):
        filename.write(txt)
    else:
        with open(filename, 'w') as unit:
            unit.write(txt)


def _convert_dict_to_structured_ndarray(data):
    """convert_dict_to_structured_ndarray

    Parameters
    ----------

    data: dictionary like object
        data structure which provides iteritems and itervalues

    returns
    -------
    tab: structured ndarray
        structured numpy array
    """
    newdtype = []
    for key, dk in iteritems(data):
        _dk = np.asarray(dk)
        dtype = _dk.dtype
        # unknown type is converted to text
        if dtype.type == np.object_:
            if len(data) == 0:
                longest = 0
            else:
                longest = len(max(_dk, key=len))
                _dk = _dk.astype('|%iS' % longest)
        if _dk.ndim > 1:
            newdtype.append((str(key), _dk.dtype, (_dk.shape[1],)))
        else:
            newdtype.append((str(key), _dk.dtype))
    tab = np.rec.fromarrays(itervalues(data), dtype=newdtype)
    return tab


def __indent__(rows, header=None, units=None, headerChar='-',
               delim=' | ', endline='\n', **kwargs):
    """Indents a table by column.

    Parameters
    ----------
    rows: sequences of rows
        one sequence per row.

    header: sequence of str
        row consists of the columns' names

    units: sequence of str
        Sequence of units

    headerChar: char
        Character to be used for the row separator line

    delim: char
        The column delimiter.

    returns
    -------
    txt: str
        string represation of rows
    """
    length_data = list(map(max, zip(*[list(map(len, k)) for k in rows])))
    length = length_data[:]

    if (header is not None):
        length_header = list(map(len, header))
        length = list(map(max, zip(length_data, length_header)))

    if (units is not None):
        length_units = list(map(len, units))
        length = list(map(max, zip(length_data, length_units)))

    if headerChar not in (None, '', ' '):
        rowSeparator = headerChar * (sum(length) + len(delim) * (len(length) - 1)) + endline
    else:
        rowSeparator = ''

    # make the format
    fmt = ['{{{0:d}:{1:d}s}}'.format(k, l) for (k, l) in enumerate(length)]
    fmt = delim.join(fmt) + endline
    # write the string
    txt = rowSeparator
    if header is not None:
        txt += fmt.format(*header)  # + endline
        txt += rowSeparator
    if units is not None:
        txt += fmt.format(*units)  # + endline
        txt += rowSeparator
    for r in rows:
        txt += fmt.format(*r)  # + endline
    txt += rowSeparator
    return txt


def pprint_rec_entry(data, num=0, keys=None):
        """ print one line with key and values properly to be readable

        Parameters
        ----------
        data: recarray
            data to extract entry from

        num: int, slice
            indice selection

        keys: sequence or str
            if str, can be a regular expression
            if sequence, the sequence of keys to print
        """
        if (keys is None) or (keys == '*'):
            _keys = data.dtype.names
        elif type(keys) in basestring:
            _keys = [k for k in data.dtype.names if (re.match(keys, k) is not None)]
        else:
            _keys = keys

        length = max(map(len, _keys))
        fmt = '{{0:{0:d}s}}: {{1}}'.format(length)
        data = data[num]

        for k in _keys:
            print(fmt.format(k, data[k]))


def pprint_rec_array(data, idx=None, fields=None, ret=False, all=False,
                     headerChar='-', delim=' | ', endline='\n' ):
        """ Pretty print the table content
            you can select the table parts to display using idx to
            select the rows and fields to only display some columns
            (ret is only for insternal use)

        Parameters
        ----------
        data: array
            array to show

        idx: sequence, slide
            sub selection to print

        fields: str, sequence
            if str can be a regular expression, and/or list of fields separated
            by spaces or commas

        ret: bool
            if set return the string representation instead of printing the result

        all: bool
            if set, force to show all rows

        headerChar: char
            Character to be used for the row separator line

        delim: char
            The column delimiter.
        """
        if (fields is None) or (fields == '*'):
            _keys = data.dtype.names
        elif type(fields) in basestring:
            if ',' in fields:
                _fields = fields.split(',')
            elif ' ' in fields:
                _fields = fields.split()
            else:
                _fields = [fields]
            lbls = data.dtype.names
            _keys = []
            for _fk in _fields:
                _keys += [k for k in lbls if (re.match(_fk, k) is not None)]
        else:
            lbls = data.dtype.names
            _keys = []
            for _fk in _fields:
                _keys += [k for k in lbls if (re.match(_fk, k) is not None)]

        nfields = len(_keys)
        nrows = len(data)
        fields = list(_keys)

        if idx is None:
            if (nrows < 10) or (all is True):
                rows = [ [ str(data[k][rk]) for k in _keys ] for rk in range(nrows)]
            else:
                _idx = range(6)
                rows = [ [ str(data[k][rk]) for k in _keys ] for rk in range(5) ]
                if nfields > 1:
                    rows += [ ['...' for k in range(nfields) ] ]
                else:
                    rows += [ ['...' for k in range(nfields) ] ]
                rows += [ [ str(data[k][rk]) for k in fields ] for rk in range(-5, 0)]
        elif isinstance(idx, slice):
            _idx = range(idx.start, idx.stop, idx.step or 1)
            rows = [ [ str(data[k][rk]) for k in fields ] for rk in _idx]
        else:
            rows = [ [ str(data[k][rk]) for k in fields ] for rk in idx]

        out = __indent__(rows, header=_keys, units=None, delim=delim,
                         headerChar=headerChar, endline=endline)
        if ret is True:
            return out
        else:
            print(out)


# ==============================================================================
# SimpleTable -- provides table manipulations with limited storage formats
# ==============================================================================
class SimpleTable(object):
    """ Table class that is designed to be the basis of any format wrapping
    around numpy recarrays

    Attributes
    ----------

    fname: str or object
        if str, the file to read from. This may be limited to the format
        currently handled automatically. If the format is not correctly handled,
        you can try by providing an object.__

        if object with a structure like dict, ndarray, or recarray-like
            the data will be encapsulated into a Table

    caseless: bool
        if set, column names will be caseless during operations

    aliases: dict
        set of column aliases (can be defined later :func:`set_alias`)

    units: dict
        set of column units (can be defined later :func:`set_unit`)

    desc: dict
        set of column description or comments (can be defined later :func:`set_comment`)

    header: dict
        key, value pair corresponding to the attributes of the table
    """

    def __init__(self, fname, *args, **kwargs):

        dtype = kwargs.pop('dtype', None)
        self.caseless = kwargs.get('caseless', False)
        self._aliases = kwargs.get('aliases', {})
        self._units = kwargs.get('units', {})
        self._desc = kwargs.get('desc', {})

        if (type(fname) == dict) or (dtype in [dict, 'dict']):
            self.header = fname.pop('header', {})
            self.data = _convert_dict_to_structured_ndarray(fname)
        elif (type(fname) in (str,)) or (dtype is not None):
            if (type(fname) in (str,)):
                extension = fname.split('.')[-1]
            else:
                extension = None
            if (extension == 'csv') or dtype == 'csv':
                kwargs.setdefault('delimiter', ',')
                kwargs.setdefault('comments', '#')
                commentedHeader = kwargs.pop('commentedHeader', False)
                n, header, units, comments, aliases, names = _ascii_read_header(fname, commentedHeader=commentedHeader, **kwargs)
                kwargs.setdefault('names', names)
                kwargs.setdefault('skip_header', n)
                self.data = np.recfromcsv(fname, *args, **kwargs)
                self.header = header
                self._units.update(**units)
                self._desc.update(**comments)
                self._aliases.update(**aliases)
                kwargs.setdefault('names', True)
            elif (extension in ('tsv', 'dat', 'txt')) or dtype in ('tsv', 'dat', 'txt'):
                kwargs.setdefault('delimiter', None)
                kwargs.setdefault('comments', '#')
                commentedHeader = kwargs.pop('commentedHeader', True)
                n, header, units, comments, aliases, names = _ascii_read_header(fname, commentedHeader=commentedHeader, **kwargs)
                kwargs.setdefault('names', names)
                kwargs.setdefault('skip_header', n)
                self.data = np.recfromtxt(fname, *args, **kwargs)
                self.header = header
                self._units.update(**units)
                self._desc.update(**comments)
                self._aliases.update(**aliases)
            elif (extension == 'fits') or dtype == 'fits':
                if pyfits is None:
                    raise RuntimeError('Cannot read this format, Astropy or pyfits not found')
                if ('extname' not in kwargs) and ('ext' not in kwargs) and (len(args) == 0):
                    args = (1, )
                self.data = np.array(pyfits.getdata(fname, *args, **kwargs))
                header, aliases, units, comments = _fits_read_header(pyfits.getheader(fname, *args, **kwargs))
                self.header = header
                self._desc.update(**comments)
                self._units.update(**units)
                self._aliases.update(**aliases)
            elif (extension in ('hdf5', 'hd5', 'hdf')) or dtype in (extension in ('hdf5', 'hd5', 'hdf')):
                if tables is None:
                    raise RuntimeError('Cannot read this format, pytables not found')
                hdr, aliases, units, desc, data = _hdf5_read_data(fname, *args, **kwargs)
                self.data = data
                self.header = hdr
                self._units.update(**units)
                self._desc.update(**desc)
                self._aliases.update(**aliases)
            else:
                raise Exception('Format {0:s} not handled'.format(extension))
        elif type(fname) == np.ndarray:
            self.data = fname
            self.header = {}
        elif type(fname) == pyfits.FITS_rec:
            self.data = np.array(fname)
            self.header = {}
        elif type(fname) == SimpleTable:
            cp = kwargs.pop('copy', True)
            if cp:
                self.data = deepcopy(fname.data)
                self.header = deepcopy(fname.header)
                self._aliases = deepcopy(fname._aliases)
                self._units = deepcopy(fname._units)
                self._desc = deepcopy(fname._desc)
            else:
                self.data = fname.data
                self.header = fname.header
                self._aliases = fname._aliases
                self._units = fname._units
                self._desc = fname._desc
        else:
            raise Exception('Type {0!s:s} not handled'.format(type(fname)))
        if 'NAME' not in self.header:
            if type(fname) not in basestring:
                self.header['NAME'] = 'No Name'
            else:
                self.header['NAME'] = fname

    def pprint_entry(self, num, keys=None):
        """ print one line with key and values properly to be readable

        Parameters
        ----------
        num: int, slice
            indice selection

        keys: sequence or str
            if str, can be a regular expression
            if sequence, the sequence of keys to print
        """
        if (keys is None) or (keys == '*'):
            _keys = self.keys()
        elif type(keys) in basestring:
            _keys = [k for k in (self.keys() + tuple(self._aliases.keys()))
                     if (re.match(keys, k) is not None)]
        else:
            _keys = keys

        length = max(map(len, _keys))
        fmt = '{{0:{0:d}s}}: {{1}}'.format(length)
        data = self[num]

        for k in _keys:
            print(fmt.format(k, data[self.resolve_alias(k)]))

    def pprint(self, idx=None, fields=None, ret=False, all=False,
               full_match=False, headerChar='-', delim=' | ', endline='\n',
               **kwargs):
        """ Pretty print the table content
            you can select the table parts to display using idx to
            select the rows and fields to only display some columns
            (ret is only for insternal use)

        Parameters
        ----------

        idx: sequence, slide
            sub selection to print

        fields: str, sequence
            if str can be a regular expression, and/or list of fields separated
            by spaces or commas

        ret: bool
            if set return the string representation instead of printing the result

        all: bool
            if set, force to show all rows

        headerChar: char
            Character to be used for the row separator line

        delim: char
            The column delimiter.
        """
        if full_match is True:
            fn = re.fullmatch
        else:
            fn = re.match

        if (fields is None) or (fields == '*'):
            _keys = self.keys()
        elif type(fields) in basestring:
            if ',' in fields:
                _fields = fields.split(',')
            elif ' ' in fields:
                _fields = fields.split()
            else:
                _fields = [fields]
            lbls = self.keys() + tuple(self._aliases.keys())
            _keys = []
            for _fk in _fields:
                _keys += [k for k in lbls if (fn(_fk, k) is not None)]
        else:
            lbls = self.keys() + tuple(self._aliases.keys())
            _keys = []
            for _fk in _fields:
                _keys += [k for k in lbls if (fn(_fk, k) is not None)]

        nfields = len(_keys)

        fields = list(map( self.resolve_alias, _keys ))

        if idx is None:
            if (self.nrows < 10) or all:
                rows = [ [ str(self[k][rk]) for k in _keys ] for rk in range(self.nrows)]
            else:
                _idx = range(6)
                rows = [ [ str(self[k][rk]) for k in _keys ] for rk in range(5) ]
                if nfields > 1:
                    rows += [ ['...' for k in range(nfields) ] ]
                else:
                    rows += [ ['...' for k in range(nfields) ] ]
                rows += [ [ str(self[k][rk]) for k in fields ] for rk in range(-5, 0)]
        elif isinstance(idx, slice):
            _idx = range(idx.start, idx.stop, idx.step or 1)
            rows = [ [ str(self[k][rk]) for k in fields ] for rk in _idx]
        else:
            rows = [ [ str(self[k][rk]) for k in fields ] for rk in idx]

        if len(self._units) == 0:
            units = None
        else:
            units = [ '(' + str( self._units.get(k, None) or '') + ')' for k in fields ]

        out = __indent__(rows, header=_keys, units=units, delim=delim,
                         headerChar=headerChar, endline=endline)
        if ret is True:
            return out
        else:
            print(out)

    def write(self, fname, **kwargs):
        """ write table into file

        Parameters
        ----------
        fname: str
            filename to export the table into

        .. note::
            additional keywords are forwarded to the corresponding libraries
            :func:`pyfits.writeto` or :func:`pyfits.append`
            :func:`np.savetxt`
        """
        extension = fname.split('.')[-1]
        if (extension == 'csv'):
            comments = kwargs.pop('comments', '#')
            delimiter = kwargs.pop('delimiter', ',')
            commentedHeader = kwargs.pop('commentedHeader', False)
            hdr = _ascii_generate_header(self, comments=comments, delimiter=delimiter,
                                         commentedHeader=commentedHeader)
            np.savetxt(fname, self.data, delimiter=delimiter, header=hdr,
                       comments='', **kwargs)
        elif (extension in ['txt', 'dat']):
            comments = kwargs.pop('comments', '#')
            delimiter = kwargs.pop('delimiter', ' ')
            commentedHeader = kwargs.pop('commentedHeader', True)
            hdr = _ascii_generate_header(self, comments=comments, delimiter=delimiter,
                                         commentedHeader=commentedHeader)
            np.savetxt(fname, self.data, delimiter=delimiter, header=hdr,
                       comments='', **kwargs)
        elif (extension == 'fits'):
            hdr0 = kwargs.pop('header', None)
            append = kwargs.pop('append', False)
            hdr = _fits_generate_header(self)
            if hdr0 is not None:
                hdr.update(**hdr0)
            if append:
                _fits_append(fname, self.data, hdr, **kwargs)
            else:
                # patched version to correctly include the header
                _fits_writeto(fname, self.data, hdr, **kwargs)
        elif (extension in ('hdf', 'hdf5', 'hd5')):
            _hdf5_write_data(fname, self.data, header=self.header,
                             units=self._units, comments=self._desc,
                             aliases=self._aliases, **kwargs)
        else:
            raise Exception('Format {0:s} not handled'.format(extension))

    def set_alias(self, alias, colname):
        """
        Define an alias to a column

        Parameters
        ----------
        alias: str
            The new alias of the column

        colname: str
            The column being aliased
        """
        if (colname not in self.keys()):
            raise KeyError("Column {0:s} does not exist".format(colname))
        self._aliases[alias] = colname

    def reverse_alias(self, colname):
        """
        Return aliases of a given column.

        Given a colname, return a sequence of aliases associated to this column
        Aliases are defined by using .define_alias()
        """
        _colname = self.resolve_alias(colname)
        if (_colname not in self.keys()):
            raise KeyError("Column {0:s} does not exist".format(colname))

        return tuple([ k for (k, v) in self._aliases.iteritems() if (v == _colname) ])

    def resolve_alias(self, colname):
        """
        Return the name of an aliased column.

        Given an alias, return the column name it aliases. This
        function is a no-op if the alias is a column name itself.

        Aliases are defined by using .define_alias()
        """
        # User aliases
        if hasattr(colname, '__iter__') & (type(colname) not in basestring):
            return [ self.resolve_alias(k) for k in colname ]
        else:
            if self.caseless is True:
                maps = dict( [ (k.lower(), v) for k, v in self._aliases.items() ] )
                maps.update( (k.lower(), k) for k in self.keys() )
                return maps.get(colname.lower(), colname)
            else:
                return self._aliases.get(colname, colname)

    def set_unit(self, colname, unit):
        """ Set the unit of a column referenced by its name

        Parameters
        ----------
        colname: str
            column name or registered alias

        unit: str
            unit description
        """
        if isinstance(unit, basestring) and isinstance(colname, basestring):
            self._units[self.resolve_alias(colname)] = str(unit)
        else:
            for k, v in zip(colname, unit):
                self._units[self.resolve_alias(k)] = str(v)

    def set_comment(self, colname, comment):
        """ Set the comment of a column referenced by its name

        Parameters
        ----------
        colname: str
            column name or registered alias

        comment: str
            column description
        """
        if isinstance(comment, basestring) and isinstance(colname, basestring):
            self._desc[self.resolve_alias(colname)] = str(comment)
        else:
            for k, v in zip(colname, comment):
                self._desc[self.resolve_alias(k)] = str(v)

    def keys(self, regexp=None, full_match=False):
        """
        Return the data column names or a subset of it

        Parameters
        ----------
        regexp: str
            pattern to filter the keys with

        full_match: bool
            if set, use :func:`re.fullmatch` instead of :func:`re.match`

        Try to apply the pattern at the start of the string, returning
        a match object, or None if no match was found.

        returns
        -------
        seq: sequence
            sequence of keys
        """
        if (regexp is None) or (regexp == '*'):
            return self.colnames
        elif type(regexp) in basestring:
            if full_match is True:
                fn = re.fullmatch
            else:
                fn = re.match

            if regexp.count(',') > 0:
                _re = regexp.split(',')
            elif regexp.count(' ') > 0:
                _re = regexp.split()
            else:
                _re = [regexp]

            lbls = self.colnames + tuple(self._aliases.keys())
            _keys = []
            for _rk in _re:
                _keys += [k for k in lbls if (fn(_rk, k) is not None)]

            return _keys
        elif hasattr(regexp, '__iter__'):
            _keys = []
            for k in regexp:
                _keys += self.keys(k)
            return _keys
        else:
            raise ValueError('Unexpected type {0} for regexp'.format(type(regexp)))

    @property
    def name(self):
        """ name of the table given by the Header['NAME'] attribute """
        return self.header.get('NAME', None)

    @property
    def colnames(self):
        """ Sequence of column names """
        return self.data.dtype.names

    @property
    def ncols(self):
        """ number of columns """
        return len(self.colnames)

    @property
    def nrows(self):
        """ number of lines """
        return len(self.data)

    @property
    def nbytes(self):
        """ number of bytes of the object """
        n = sum(k.nbytes if hasattr(k, 'nbytes') else sys.getsizeof(k)
                for k in self.__dict__.values())
        return n

    def __len__(self):
        """ number of lines """
        return self.nrows

    @property
    def shape(self):
        """ shape of the data """
        return self.data.shape

    @property
    def dtype(self):
        """ dtype of the data """
        return self.data.dtype

    def __getitem__(self, v):
        return np.asarray(self.data.__getitem__(self.resolve_alias(v)))

    def get(self, v, full_match=False):
        """ returns a table from columns given as v

        this function is equivalent to :func:`__getitem__` but preserve the
        Table format and associated properties (units, description, header)

        Parameters
        ----------
        v: str
            pattern to filter the keys with

        full_match: bool
            if set, use :func:`re.fullmatch` instead of :func:`re.match`

        """
        new_keys = self.keys(v)
        t = self.__class__(self[new_keys])
        t.header.update(**self.header)
        t._aliases.update((k, v) for (k, v) in self._aliases.items() if v in new_keys)
        t._units.update((k, v) for (k, v) in self._units.items() if v in new_keys)
        t._desc.update((k, v) for (k, v) in self._desc.items() if v in new_keys)
        return t

    def __setitem__(self, k, v):
        if k in self:
            return self.data.__setitem__(self.resolve_alias(k), v)
        else:
            object.__setitem__(self, k, v)

    def __getattr__(self, k):
        try:
            return self.data.__getitem__(self.resolve_alias(k))
        except:
            return object.__getattribute__(self, k)

    def __iter__(self):
        return self.data.__iter__()

    def iterkeys(self):
        """ Iterator over the columns of the table """
        for k in self.colnames:
            yield k

    def itervalues(self):
        """ Iterator over the lines of the table """
        for l in self.data:
            yield l

    def info(self):
        s = "\nTable: {name:s}\n       nrows={s.nrows:d}, ncols={s.ncols:d}, mem={size:s}"
        s = s.format(name=self.header.get('NAME', 'Noname'), s=self,
                     size=pretty_size_print(self.nbytes))

        s += '\n\nHeader:\n'
        vals = list(self.header.items())
        length = max(map(len, self.header.keys()))
        fmt = '\t{{0:{0:d}s}} {{1}}\n'.format(length)
        for k, v in vals:
            s += fmt.format(k, v)

        vals = [(k, self._units.get(k, ''), self._desc.get(k, ''))
                for k in self.colnames]
        lengths = [(len(k), len(self._units.get(k, '')), len(self._desc.get(k, '')))
                   for k in self.colnames]
        lengths = list(map(max, (zip(*lengths))))

        s += '\nColumns:\n'

        fmt = '\t{{0:{0:d}s}} {{1:{1:d}s}} {{2:{2:d}s}}\n'.format(*(k + 1 for k in lengths))
        for k, u, c in vals:
            s += fmt.format(k, u, c)

        print(s)

        if len(self._aliases) > 0:
            print("\nTable contains alias(es):")
            for k, v in self._aliases.items():
                print('\t{0:s} --> {1:s}'.format(k, v))

    def __repr__(self):
        s = object.__repr__(self)
        s += "\nTable: {name:s}\n       nrows={s.nrows:d}, ncols={s.ncols:d}, mem={size:s}"
        return s.format(name=self.header.get('NAME', 'Noname'), s=self,
                        size=pretty_size_print(self.nbytes))

    def __getslice__(self, i, j):
        return self.data.__getslice__(i, j)

    def __contains__(self, k):
        return (k in self.colnames) or (k in self._aliases)

    def __array__(self):
        return self.data

    def __call__(self, *args, **kwargs):
        if (len(args) > 0) or (len(kwargs) > 0):
            return self.evalexpr(*args, **kwargs)
        else:
            return self.info()

    def sort(self, keys, copy=False):
        """
        Sort the table inplace according to one or more keys. This operates on
        the existing table (and does not return a new table).

        Parameters
        ----------

        keys: str or seq(str)
            The key(s) to order by

        copy: bool
            if set returns a sorted copy instead of working inplace
        """
        if not hasattr(keys, '__iter__'):
            keys = [keys]

        if copy is False:
            self.data.sort(order=keys)
        else:
            t = self.__class__(self, copy=True)
            t.sort(keys, copy=False)
            return t

    def match(self, r2, key):
        """ Returns the indices at which the tables match
        matching uses 2 columns that are compared in values

        Parameters
        ----------
        r2:  Table
            second table to use

        key: str
            fields used for comparison.

        Returns
        -------
        indexes: tuple
            tuple of both indices list where the two columns match.
        """
        return np.where( np.equal.outer( self[key], r2[key] ) )

    def stack(self, r, defaults=None):
        """
        Superposes arrays fields by fields inplace

        Parameters
        ----------
        r: Table
        """
        if not hasattr(r, 'data'):
            raise AttributeError('r should be a Table object')
        self.data = recfunctions.stack_arrays([self.data, r.data], defaults,
                                              usemask=False, asrecarray=True)

    def join_by(self, r2, key, jointype='inner', r1postfix='1', r2postfix='2',
                defaults=None, asrecarray=False, asTable=True):
        """
        Join arrays `r1` and `r2` on key `key`.

        The key should be either a string or a sequence of string corresponding
        to the fields used to join the array.
        An exception is raised if the `key` field cannot be found in the two input
        arrays.
        Neither `r1` nor `r2` should have any duplicates along `key`: the presence
        of duplicates will make the output quite unreliable. Note that duplicates
        are not looked for by the algorithm.

        Parameters
        ----------
        key: str or seq(str)
            corresponding to the fields used for comparison.

        r2: Table
            Table to join with

        jointype: str in {'inner', 'outer', 'leftouter'}
            * 'inner'     : returns the elements common to both r1 and r2.
            * 'outer'     : returns the common elements as well as the elements of r1 not in r2 and the elements of not in r2.
            * 'leftouter' : returns the common elements and the elements of r1 not in r2.

        r1postfix: str
            String appended to the names of the fields of r1 that are present in r2

        r2postfix:  str
            String appended to the names of the fields of r2 that are present in r1

        defaults:   dict
            Dictionary mapping field names to the corresponding default values.

        Returns
        -------
        tab: Table
            joined table

        .. note::

            * The output is sorted along the key.

            * A temporary array is formed by dropping the fields not in the key
              for the two arrays and concatenating the result. This array is
              then sorted, and the common entries selected. The output is
              constructed by filling the fields with the selected entries.
              Matching is not preserved if there are some duplicates...
        """
        arr = recfunctions.join_by(key, self, r2, jointype=jointype,
                                   r1postfix=r1postfix, r2postfix=r2postfix,
                                   defaults=defaults, usemask=False,
                                   asrecarray=True)

        return SimpleTable(arr)

    @property
    def empty_row(self):
        """ Return an empty row array respecting the table format """
        return np.rec.recarray(shape=(1,), dtype=self.data.dtype)

    def add_column(self, name, data, dtype=None, unit=None, description=None):
        """
        Add one or multiple columns to the table

        Parameters
        ----------
        name: str or sequence(str)
           The name(s) of the column(s) to add

        data: ndarray, or sequence of ndarray
            The column data, or sequence of columns

        dtype: dtype
            numpy dtype for the data to add

        unit: str
            The unit of the values in the column

        description: str
            A description of the content of the column
        """

        _data = np.array(data, dtype=dtype)
        dtype = _data.dtype

        # unknown type is converted to text
        if dtype.type == np.object_:
            if len(data) == 0:
                longest = 0
            else:
                longest = len(max(data, key=len))
                _data = np.asarray(data, dtype='|%iS' % longest)

        dtype = _data.dtype

        if len(self.data.dtype) > 0:
            # existing data in the table
            self.data = recfunctions.append_fields(self.data, name, _data,
                                                   dtypes=dtype, usemask=False,
                                                   asrecarray=True)
        else:
            if _data.ndim > 1:
                newdtype = (str(name), _data.dtype, (_data.shape[1],))
            else:
                newdtype = (str(name), _data.dtype)
            self.data = np.array(_data, dtype=[newdtype])

        if unit is not None:
            self.set_unit(name, unit)

        if description is not None:
            self.set_unit(name, description)

    def append_row(self, iterable):
        """
        Append one row in this table.

        see also: :func:`stack`

        Parameters
        ----------
        iterable: iterable
            line to add
        """
        if (len(iterable) != self.ncols):
            raise AttributeError('Expecting as many items as columns')
        r = self.empty_row
        for k, v in enumerate(iterable):
            r[0][k] = v
        self.stack(r)

    def remove_columns(self, names):
        """
        Remove several columns from the table

        Parameters
        ----------
        names: sequence
            A list containing the names of the columns to remove
        """
        self.pop_columns(names)

    def pop_columns(self, names):
        """
        Pop several columns from the table

        Parameters
        ----------

        names: sequence
            A list containing the names of the columns to remove

        Returns
        -------

        values: tuple
            list of columns
        """

        if not hasattr(names, '__iter__') or type(names) in basestring:
            names = [names]

        p = [self[k] for k in names]

        _names = set([ self.resolve_alias(k) for k in names ])
        self.data = recfunctions.drop_fields(self.data, _names)
        for k in names:
            self._aliases.pop(k, None)
            self._units.pop(k, None)
            self._desc.pop(k, None)

        return p

    def find_duplicate(self, index_only=False, values_only=False):
        """Find duplication in the table entries, return a list of duplicated
        elements Only works at this time is 2 lines are *the same entry* not if
        2 lines have *the same values*
        """
        dup = []
        idd = []
        for i in range(len(self.data)):
            if (self.data[i] in self.data[i + 1:]):
                if (self.data[i] not in dup):
                    dup.append(self.data[i])
                    idd.append(i)
        if index_only:
            return idd
        elif values_only:
            return dup
        else:
            return zip(idd, dup)

    def evalexpr(self, expr, exprvars=None, dtype=float):
        """ evaluate expression based on the data and external variables
            all np function can be used (log, exp, pi...)

        Parameters
        ----------
        expr: str
            expression to evaluate on the table
            includes mathematical operations and attribute names

        exprvars: dictionary, optional
            A dictionary that replaces the local operands in current frame.

        dtype: dtype definition
            dtype of the output array

        Returns
        -------
        out : NumPy array
            array of the result
        """
        _globals = {}
        for k in ( list(self.colnames) + list(self._aliases.keys()) ):
            if k in expr:
                _globals[k] = self[k]

        if exprvars is not None:
            if (not (hasattr(exprvars, 'keys') & hasattr(exprvars, '__getitem__' ))):
                raise AttributeError("Expecting a dictionary-like as condvars")
            for k, v in ( exprvars.items() ):
                _globals[k] = v

        # evaluate expression, to obtain the final filter
        r    = np.empty( self.nrows, dtype=dtype)
        r[:] = eval(expr, _globals, np.__dict__)

        return r

    def where(self, condition, condvars=None, *args, **kwargs):
        """ Read table data fulfilling the given `condition`.
        Only the rows fulfilling the `condition` are included in the result.

        Parameters
        ----------
        condition: str
            expression to evaluate on the table
            includes mathematical operations and attribute names

        condvars: dictionary, optional
            A dictionary that replaces the local operands in current frame.

        Returns
        -------
        out: ndarray/ tuple of ndarrays
        result equivalent to :func:`np.where`

        """
        ind = np.where(self.evalexpr(condition, condvars, dtype=bool ), *args, **kwargs)
        return ind

    def select(self, fields, indices=None, **kwargs):
        """
        Select only a few fields in the table

        Parameters
        ----------
        fields: str or sequence
            fields to keep in the resulting table

        indices: sequence or slice
            extract only on these indices

        returns
        -------
        tab: SimpleTable instance
            resulting table
        """
        _fields = self.keys(fields)

        if fields == '*':
            if indices is None:
                return self
            else:
                tab = self.__class__(self[indices])
                for k in self.__dict__.keys():
                    if k not in ('data', ):
                        setattr(tab, k, deepcopy(self.__dict__[k]))
                return tab
        else:
            d = {}
            for k in _fields:
                _k = self.resolve_alias(k)
                if indices is not None:
                    d[k] = self[_k][indices]
                else:
                    d[k] = self[_k]
            d['header'] = deepcopy(self.header)
            tab = self.__class__(d)
            for k in self.__dict__.keys():
                if k not in ('data', ):
                    setattr(tab, k, deepcopy(self.__dict__[k]))
            return tab

    def selectWhere(self, fields, condition, condvars=None, **kwargs):
        """ Read table data fulfilling the given `condition`.
            Only the rows fulfilling the `condition` are included in the result.

        Parameters
        ----------
        fields: str or sequence
            fields to keep in the resulting table

        condition: str
            expression to evaluate on the table
            includes mathematical operations and attribute names

        condvars: dictionary, optional
            A dictionary that replaces the local operands in current frame.

        Returns
        -------
        tab: SimpleTable instance
            resulting table
        """
        if condition in [True, 'True', None]:
            ind = None
        else:
            ind = self.where(condition, condvars, **kwargs)

        tab = self.select(fields, indices=ind)

        return tab

    def groupby(self, key):
        """
        Create an iterator which returns (key, sub-table) grouped by each value
        of key(value)

        Parameters
        ----------
        key: str
            expression or pattern to filter the keys with

        Returns
        -------
        key: str or sequence
            group key

        tab: SimpleTable instance
           sub-table of the group
           header, aliases and column metadata are preserved (linked to the
           master table).
        """
        _key = self.keys(key)
        getter = operator.itemgetter(*_key)

        for k, grp in itertools.groupby(self.data, getter):
            t = self.__class__(np.dstack(grp))
            t.header = self.header
            t._aliases = self._aliases
            t._units = self._units
            t._desc = self._desc
            yield (k, t)

    def stats(self, fn=None, fields=None, fill=None):
        """ Make statistics on columns of a table

        Paramters
        ---------
        fn: callable or sequence of callables
            functions to apply to each column
            default: (np.mean, np.std, np.nanmin, np.nanmax)

        fields: str or sequence
            any key or key expression to subselect columns
            default is all columns

        fill: value
            value when not applicable
            default np.nan

        returns
        -------
        tab: Table instance
            collection of statistics, one column per function in fn and 1 ligne
            per column in the table
        """
        from collections import OrderedDict

        fn = (stats.mean, stats.std,
              stats.min, stats.max,
              stats.has_nan)

        d = OrderedDict()
        d.setdefault('FIELD', [])
        for k in fn:
            d.setdefault(k.__name__, [])

        if fields is None:
            fields = self.colnames
        else:
            fields = self.keys(fields)

        if fill is None:
            fill = np.nan

        for k in fields:
            d['FIELD'].append(k)
            for fnk in fn:
                try:
                    val = fnk(self[k])
                except:
                    val = fill
                d[fnk.__name__].append(val)

        return self.__class__(d, dtype=dict)

    # method aliases
    remove_column = remove_columns

    # deprecated methods
    addCol = add_column
    addLine = append_row
    setComment = set_comment
    setUnit = set_unit
    delCol = remove_columns


class stats(object):
    @classmethod
    def has_nan(s, v):
        return (True in np.isnan(v))

    @classmethod
    def mean(s, v):
        return np.nanmean(v)

    @classmethod
    def max(s, v):
        return np.nanmax(v)

    @classmethod
    def min(s, v):
        return np.nanmin(v)

    @classmethod
    def std(s, v):
        return np.nanstd(v)

    @classmethod
    def var(s, v):
        return np.var(v)

    @classmethod
    def p16(s, v):
        try:
            return np.nanpercentile(v, 16)
        except AttributeError:
            return np.percentile(v, 16)

    @classmethod
    def p84(s, v):
        try:
            return np.nanpercentile(v, 84)
        except AttributeError:
            return np.percentile(v, 84)

    @classmethod
    def p50(s, v):
        try:
            return np.nanmedian(v)
        except AttributeError:
            return np.percentile(v, 50)


# =============================================================================
# Adding some plotting functions
# =============================================================================

try:
    import pylab as plt

    def plot_function(tab, fn, *args, **kwargs):
        """ Generate a plotting method of tab from a given function

        Parameters
        ----------
        tab: SimpleTable instance
            table instance

        fn: str or callable
            if str, will try a function in matplotlib
            if callable, calls the function directly

        xname: str
            expecting a column name from the table

        yname: str, optional
            if provided, another column to use for the plot

        onlywhere: sequence or str, optional
            if provided, selects only data with this condition
            the condition can be a ndarray slice or a string.
            When a string is given, the evaluation calls :func:`SimpleTable.where`

        ax: matplotlib.Axes instance
            if provided make sure it uses the axis to do the plots if a mpl
            function is used.

        Returns
        -------
        r: object
            anything returned by the called function
        """
        if not hasattr(fn, '__call__'):
            ax = kwargs.pop('ax', None)
            if ax is None:
                ax = plt.gca()
            _fn = getattr(ax, fn, None)
            if _fn is None:
                raise AttributeError('function neither callable or found in matplotlib')
        else:
            _fn = fn

        onlywhere = kwargs.pop('onlywhere', None)
        if type(onlywhere) in basestring:
            select = tab.where(onlywhere)
        else:
            select = onlywhere

        _args = ()
        for a in args:
            if (hasattr(a, '__iter__')):
                try:
                    b = tab[a]
                    if select is not None:
                        b = b.compress(select)
                    if (len(b.dtype) > 1):
                        b = list((b[k] for k in b.dtype.names))
                    _args += (b, )
                except Exception as e:
                    print(e)
                    _args += (a, )
            else:
                _args += (a, )

        return _fn(*_args, **kwargs)

    def attached_function(fn, doc=None, errorlevel=0):
        """ eclare a function as a method to the class table"""

        def _fn(self, *args, **kwargs):
            try:
                return plot_function(self, fn, *args, **kwargs)
            except Exception as e:
                if errorlevel < 1:
                    pass
                else:
                    raise e

        if doc is not None:
            _fn.__doc__ = doc

        return _fn

    SimpleTable.plot_function = plot_function
    SimpleTable.plot = attached_function('plot', plt.plot.__doc__)
    SimpleTable.hist = attached_function('hist', plt.hist.__doc__)
    SimpleTable.hist2d = attached_function('hist2d', plt.hist2d.__doc__)
    SimpleTable.hexbin = attached_function('hexbin', plt.hexbin.__doc__)
    SimpleTable.scatter = attached_function('scatter', plt.scatter.__doc__)

    # newer version of matplotlib
    if hasattr(plt, 'violinplot'):
        SimpleTable.violinplot = attached_function('violinplot', plt.violinplot.__doc__)
    if hasattr(plt, 'boxplot'):
        SimpleTable.boxplot = attached_function('boxplot', plt.boxplot.__doc__)

except Exception as e:
    print(e)
