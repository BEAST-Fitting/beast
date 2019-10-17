import os
import inspect
import sys

localpath = "/".join(
    os.path.abspath(inspect.getfile(inspect.currentframe())).split("/")[:-1]
)
