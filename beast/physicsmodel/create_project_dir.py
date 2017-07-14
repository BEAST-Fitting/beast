from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import os

__all__ = ['create_project_dir']

def create_project_dir(project):
    """ 
    Create the project directory if necessary

    Parameters
    ----------
    project: str
        project name

    Returns
    -------
    dirname: str
        <project>/<project>

    Raises
    ------
    Exception
        if the name already exists as a file instead of a directory
    """
    outdir = project
    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            raise Exception('Output directory %s already exists but is not a directory' % outdir)
    else:
        os.mkdir(outdir)
    return '%s/%s' % (outdir, project)
