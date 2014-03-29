#!/usr/bin/env python
from __future__ import print_function
import sys
import os
if sys.version_info.major <= 3:
    import urllib2 as request
else:
    from urllib import request

from .tools import progressbar as pb
from .config import libs_server, libs, libsdir


class Reporter:
    """ Report Status during the download of files

    Attributes
    ----------
    pbar: ProgressBar
        instance of progress bar configured to indicate the downloading
        progression and status.
    """
    def __init__(self, fname):
        widgets = ['Downloading %s: ' % fname,
                   pb.Percentage(), ' ',
                   pb.Bar(marker='#', left='[', right=']'),
                   ' ', pb.ETA(), ' ', pb.FileTransferSpeed()]
        self.pbar = pb.ProgressBar(widgets=widgets, maxval=100)
        self.pbar.start()

    def __call__(self, bytes_so_far, chunk_size, total_size):
        if self.pbar.maxval != total_size:
            self.pbar.maxval = total_size
        percent = float(bytes_so_far)
        if bytes_so_far >= total_size:
            self.pbar.finish()
        else:
            self.pbar.update(percent)


def chunk_read(response, buf, chunk_size=8192, report_hook=None):
    total_size = response.info().getheader('Content-Length').strip()
    total_size = int(total_size)
    bytes_so_far = 0

    while 1:
        chunk = response.read(chunk_size)
        bytes_so_far += len(chunk)

        if not chunk:
            break
        else:
            buf.write(chunk)
        if report_hook:
            report_hook(bytes_so_far, chunk_size, total_size)

    return bytes_so_far


def urlretrieve(url, reporthook, f):
    """ Retrieve data from an url """
    response = request.urlopen(url)
    return chunk_read(response, f, report_hook=reporthook)


def _dl(urltxt, dest, fname=None):
    if fname is None:
        fname = dest
    with open(dest, 'wb') as f:
        s = urlretrieve(urltxt, Reporter(fname), f)
    fmt = '%6.2f %s'
    units = ['B', 'K', 'M', 'G', 'T', 'P']
    for u in units:
        if s < 1000:
            break
        s /= 1000
    size = fmt % (s, u)
    print('\nDownloaded %s bytes' % size)


def dl_install_all_libs(server, libs):
    libdir = os.path.abspath(libsdir)
    if os.path.exists(libdir):
        if not os.path.isdir(libdir):
            print(libdir, 'exists and is not a directory.')
            return 1
    else:
        os.makedirs(libdir)

    if server[-1] != '/':
        server += '/'
    if libdir[-1] != '/':
        libdir += '/'
    for k in libs:
        url = server + libs[k]
        _dl(url, libdir + libs[k], libs[k])

if __name__ == '__main__':
    print("""
            Installing libraries
          """)
    dl_install_all_libs(libs_server, libs)
