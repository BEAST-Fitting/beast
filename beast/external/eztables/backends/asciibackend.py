
import os
import inspect
import re
localpath = '/'.join(os.path.abspath(inspect.getfile(inspect.currentframe())).split('/')[:-1])
import numpy as np
from .basebackend import BaseBackend
from ..core.tableheader import TableHeader
from ..table import Table

#==============================================================================
# ASCII/CSV TABLE MANAGERS
#==============================================================================


#==============================================================================
class csvBackend(BaseBackend):
#==============================================================================
    def __init__(self):
        """ constructor """
        BaseBackend.__init__(self, tableType='csv')

    def readHeader(self, filename, comment='#', noheader=False, delimiter=',', *args, **kwargs):
        """ read the CSV header """
        header = None

        def parseStrNone(v):
            _v = v.split()
            if (len(_v) == 0):
                return None
            else:
                _v = ' '.join(_v)
                if (_v.lower()) == 'none' or (_v.lower() == 'null'):
                    return None
                else:
                    return _v

        if hasattr(filename, 'read'):
            stream = filename
        else:
            stream = open(filename, 'r')

        #get emtpy header
        description = TableHeader()

        #ColumnHeader(dtype, unit=None, description=None, null=None, format=None, aliases=None)
        colInfo = {}
        aliases = []
        nHeadLines = 0
        header = None
        while header is None:
            line = stream.readline()[:-1]
            nHeadLines += 1
            if line[0] == comment:
                if line[1] != comment:
                    #get table header
                    k = line[1:].split('\t')
                    key = k[0].split()[0]  # remove trailing spaces
                    #check aliases
                    if key[:5] != 'alias':
                        for cv in k[1:]:
                            description[key] = cv
                    else:
                        aliases.append(k[1].split('='))
                else:
                    #get columns meta
                    k = line[2:].split('\t')
                    colName = k[0].split()[0]
                    if len(k) > 1:
                        colUnit = parseStrNone(k[1])
                    else:
                        colUnit = None
                    if len(k) > 2:
                        colComm = parseStrNone(k[2])
                    else:
                        colComm = None
                    if len(k) > 3:
                        colNull = parseStrNone(k[3])
                    else:
                        colNull = None
                    if len(k) > 4:
                        colfmt  = parseStrNone(k[4])
                    else:
                        colfmt = None

                    colInfo[colName] = (colUnit, colComm, colNull, colfmt)
            else:
                header = line.split(delimiter)
                if noheader:
                    header = ['Col%s' % k for k in range(len(header))]

        if not hasattr(filename, 'read'):
            stream.close()

        if not 'NAME' in list(description.keys()):
            description['NAME'] = filename.split('/')[-1]
        return nHeadLines, description, colInfo, header, aliases

    def readData(self, filename, skiprows, *args, **kwargs):
        """ returns the recarray of the data """
        return np.recfromcsv(filename, skip_header=skiprows - 1, *args, **kwargs)

    def read(self, filename, delimiter=',', noheader=False, skiprows=0, comment='#', *args, **kwargs):
        """
        Read Csv file with header or not. Especially useful in association with
        exportdata module.
        So far it uses also the np.recfromcsv method
        """

        #Get header
        nHeadLines, description, colInfo, header, aliases = self.readHeader(filename, comment, noheader, delimiter)

        #get data
        if hasattr(filename, 'read'):
            if not filename.closed:
                filename.seek(0)
        skip = skiprows + (nHeadLines - int(noheader))
        d = self.readData(filename, skiprows=skip, *args, **kwargs)

        #generate an empty table and fill it
        tab = Table()
        tab.header = description

        #add each column
        for k in range( len(header) ):
            colName = header[k]
            if colName in colInfo:
                colUnit, colComm, colNull, colfmt = colInfo[colName]
            else:
                colUnit, colComm, colNull, colfmt = (None, None, None, None)

            if colName in list(tab.keys()):
                i = 1
                while '%s_%d' % (colName, i) in d.dtype.names:
                    i += 1
                colName = '%s_%d' % (colName, i)

            tab.add_column(colName, d[ re.sub(r"[/().']", '', colName.lower())],
                           unit=colUnit or '',
                           null=colNull or '',
                           description=colComm or '',
                           format=colfmt,
                           dtype=d.dtype[k] )
        #set aliases
        for k in aliases:
            tab.set_alias(k[0], k[1])

        return tab

    def writeHeader(self, unit, header, comment='#'):
        """ Write File Header definition into the opened unit
        e.g. >>> self.writeHeader(unit, header, comment='#')
            # NAME    tablename
            # KEY     Value
        """
        keys = np.sort(list(header.keys()))
        for key in keys:
            val = header[key]
            for kval in str(val).split('\n'):
                unit.write('%s %s\t%s\n' % (comment, key.upper(), kval) )

    def writeColHeader(self, unit, cols, comment='#'):
        """ Write column description into opened buffer unit
        e.g. >>> self.writeColHeader(unit, cols, comment='#')
            ## Column     Unit    Comment   Null    Fmt
        """
        keys = list(cols.keys())
        unit.write('%s%s %20s\t%10s\t%s\t%s\t%5s\n' % (comment, comment, 'Column', 'Unit', 'Comment', 'Null', 'Fmt') )
        for k in keys:
            hdr = cols[k]
            txt = "%20s\t%10s\t%s\t%s\t%5s" % (k, hdr.unit, hdr.description, hdr.null, hdr.format)
            unit.write('%s%s %s\n' % (comment, comment, txt) )

    def writeAliasesDef(self, unit, aliases, comment='#'):
        """ Write aliases into the header
        e.g. >>> self.writeAliasesDef(unit, aliases, comment='#')
            # alias   col_alias=col1
        """
        for k, v in aliases.items():
            unit.write('%s alias\t%s=%s\n' % (comment, k, v) )

    def writeColDef(self, unit, cols, delimiter=',', comment=''):
        """ Write column definition into the opened unit
            This corresponds to the first line of the data for
            standards compatibility
        e.g. >>> self.writeColDef(unit, cols, comment='#')
            #col1,col2,col3,col4,...,coln
        """
        unit.write( comment + delimiter.join( list(cols.keys()) ) + '\n' )

    def writeData(self, unit, data, fmt, delimiter=','):
        """ Write data part into the opened unit """
        size = data.shape[0]
        for ik in range(size):
            unit.write( fmt % data[ik].tolist() )
            unit.write("\n")

    def write(self, tab, output='exportedData.csv', header=True,
              delimiter=',', comment='#', keep=False, verbose=False, **kwargs):
        """
        export data to a comma separated value file

        inputs:
            data -- data dictionnary to export

        outputs:
            output -- output file (def: exportedData.dat)

        keywords:
            header    -- header object describing the file
            delimiter -- delimiter to use (def: ',')
            comment   -- comment character for header (def: '#')
            keep      -- keeps unit opened
            unit      -- uses opened stream, if provided
        """

        if hasattr(output, 'write'):
            unit = output
        else:
            unit = open(output, 'w')
        if header:
            self.writeHeader(unit, tab.header, comment=comment)
            self.writeColHeader(unit, tab.columns, comment=comment)
            self.writeAliasesDef(unit, tab._aliases, comment=comment)
            self.writeColDef(unit, tab.columns, comment='', delimiter=delimiter)

        fmt  = delimiter.join(['%' + tab.columns[k].format for k in tab.columns])
        self.writeData( unit, tab.data[list(tab.keys())], fmt, delimiter=delimiter)

        if hasattr(output, 'write') or keep:
            return unit
        else:
            unit.close()


#==============================================================================
class asciiBackend(BaseBackend):
#==============================================================================
    def __init__(self):
        """ constructor """
        BaseBackend.__init__(self, tableType='csv')

    def readHeader(self, filename, comment='#', noheader=False, delimiter=None, *args, **kwargs):
        """ read the CSV header """
        header = None

        def parseStrNone(v):
            _v = v.split()
            if (len(_v) == 0):
                return None
            else:
                _v = ' '.join(_v)
                if (_v.lower()) == 'none' or (_v.lower() == 'null'):
                    return None
                else:
                    return _v

        if hasattr(filename, 'read'):
            stream = filename
        else:
            stream = open(filename, 'r')

        #get emtpy header
        description = TableHeader()

        #ColumnHeader(dtype, unit=None, description=None, null=None, format=None, aliases=None)
        colInfo = {}
        aliases = []
        nHeadLines = 0
        header = None
        oldline = ''
        while header is None:
            line = stream.readline()[:-1]
            nHeadLines += 1
            if line[0] == comment:
                if line[1] != comment:
                    #get table header
                    k = line[1:].split('\t')
                    if k[0] != '':
                        key = k[0].split()[0]  # remove trailing spaces
                    else:
                        key = k[1].split()[0]
                    #check aliases
                    if key[:5] != 'alias':
                        for cv in k[1:]:
                            description[key] = cv
                    else:
                        aliases.append(k[1].split('='))
                else:
                    #get columns meta
                    k = line[2:].split('\t')
                    colName = k[0].split()[0]
                    if len(k) > 1:
                        colUnit = parseStrNone(k[1])
                    else:
                        colUnit = None
                    if len(k) > 2:
                        colComm = parseStrNone(k[2])
                    else:
                        colComm = None
                    if len(k) > 3:
                        colNull = parseStrNone(k[3])
                    else:
                        colNull = None
                    if len(k) > 4:
                        colfmt  = parseStrNone(k[4])
                    else:
                        colfmt = None

                    colInfo[colName] = (colUnit, colComm, colNull, colfmt)
                oldline = line
            else:
                header = oldline[1:].split(delimiter)
                if noheader:
                    header = ['Col%s' % k for k in range(len(header))]

        if not hasattr(filename, 'read'):
            stream.close()

        if not 'NAME' in list(description.keys()):
            description['NAME'] = filename.split('/')[-1]
        return nHeadLines, description, colInfo, header, aliases

    def readData(self, filename, skiprows, *args, **kwargs):
        """ returns the recarray of the data """
        return np.recfromtxt(filename, skip_header=skiprows - 1, deletechars="[/().']", *args, **kwargs)

    def read(self, filename, delimiter=None, noheader=False, skiprows=0, comment='#', *args, **kwargs):
        """
        Read Csv file with header or not. Especially useful in association with
        exportdata module.
        So far it uses also the np.recfromtxt method
        """

        #Get header
        nHeadLines, description, colInfo, header, aliases = self.readHeader(filename, comment, noheader, delimiter)

        #get data
        if hasattr(filename, 'read'):
            if not filename.closed:
                filename.seek(0)
        skip = skiprows + (nHeadLines - int(noheader))
        d = self.readData(filename, skiprows=skip, names=header, *args, **kwargs)

        #generate an empty table and fill it
        tab = Table()
        tab.header = description

        #add each column
        for k in range( len(header) ):
            colName = header[k]
            if colName in colInfo:
                colUnit, colComm, colNull, colfmt = colInfo[colName]
            else:
                colUnit, colComm, colNull, colfmt = (None, None, None, None)

            if colName in list(tab.keys()):
                i = 1
                while '%s_%d' % (colName, i) in d.dtype.names:
                    i += 1
                colName = '%s_%d' % (colName, i)
            _key = re.sub(r"[/().']", '', colName)
            if _key in d.dtype.names:
                tab.add_column(colName, d[_key],
                               unit=colUnit or '',
                               null=colNull or '',
                               description=colComm or '',
                               format=colfmt,
                               dtype=d.dtype[k] )
            else:
                tab.add_column(colName, d[ _key ],
                               unit=colUnit or '',
                               null=colNull or '',
                               description=colComm or '',
                               format=colfmt,
                               dtype=d.dtype[k] )
        #set aliases
        for k in aliases:
            tab.set_alias(k[0], k[1])

        return tab

    def writeHeader(self, unit, header, comment='#'):
        """ Write File Header definition into the opened unit
        e.g. >>> self.writeHeader(unit, header, comment='#')
            # NAME    tablename
            # KEY     Value
        """
        keys = np.sort(list(header.keys()))
        for key in keys:
            val = header[key]
            for kval in str(val).split('\n'):
                unit.write('%s %s\t%s\n' % (comment, key.upper(), kval) )

    def writeColHeader(self, unit, cols, comment='#'):
        """ Write column description into opened buffer unit
        e.g. >>> self.writeColHeader(unit, cols, comment='#')
            ## Column     Unit    Comment   Null    Fmt
        """
        keys = list(cols.keys())
        unit.write('%s%s %20s\t%10s\t%s\t%s\t%5s\n' % (comment, comment, 'Column', 'Unit', 'Comment', 'Null', 'Fmt') )
        for k in keys:
            hdr = cols[k]
            txt = "%20s\t%10s\t%s\t%s\t%5s" % (k, hdr.unit, hdr.description, hdr.null, hdr.format)
            unit.write('%s%s %s\n' % (comment, comment, txt) )

    def writeAliasesDef(self, unit, aliases, comment='#'):
        """ Write aliases into the header
        e.g. >>> self.writeAliasesDef(unit, aliases, comment='#')
            # alias   col_alias=col1
        """
        for k, v in aliases.items():
            unit.write('%s alias\t%s=%s\n' % (comment, k, v) )

    def writeColDef(self, unit, cols, delimiter=',', comment=''):
        """ Write column definition into the opened unit
            This corresponds to the first line of the data for
            standards compatibility
        e.g. >>> self.writeColDef(unit, cols, comment='#')
            #col1,col2,col3,col4,...,coln
        """
        unit.write( comment + delimiter.join( list(cols.keys()) ) + '\n' )

    def writeData(self, unit, data, fmt, delimiter=None):
        """ Write data part into the opened unit """
        size = data.shape[0]
        for ik in range(size):
            unit.write( fmt % data[ik].tolist() )
            unit.write("\n")

    def write(self, tab, output='exportedData.txt', header=True,
              delimiter=None, comment='#', verbose=False, keep=False, **kwargs):
        """
        export data to a comma separated value file

        inputs:
            data -- data dictionnary to export

        outputs:
            output -- output file (def: exportedData.dat)

        keywords:
            header    -- header object describing the file
            delimiter -- delimiter to use (def: ',')
            comment   -- comment character for header (def: '#')
            keep      -- keeps unit opened
            unit      -- uses opened stream, if provided
        """
        if delimiter is None:
            delimiter = ' '

        if hasattr(output, 'write'):
            unit = output
        else:
            unit = open(output, 'w')
        if header:
            self.writeHeader( unit, tab.header, comment=comment)
            self.writeColHeader( unit, tab.columns, comment=comment)
            self.writeAliasesDef( unit, tab._aliases, comment=comment)
            self.writeColDef( unit, tab.columns, comment='#', delimiter=delimiter)

        fmt  = delimiter.join(['%' + tab.columns[k].format for k in tab.columns])
        self.writeData( unit, tab.data[list(tab.keys())], fmt, delimiter=delimiter)

        if hasattr(output, 'write') or keep:
            return unit
        else:
            unit.close()
