""" sqllite, mysql, psql """

import os
import inspect

import numpy as np

localpath = "/".join(
    os.path.abspath(inspect.getfile(inspect.currentframe())).split("/")[:-1]
)

# from basebackend import BaseBackend

# SQLite
try:
    import sqlite3
    sqlite3_installed = True
except Exception:
    sqlite3_installed = False

# Type conversion dictionary
type_dict = {}
type_dict[np.bool_] = "BOOL"
type_dict[np.uint8] = "TINYINT"
type_dict[np.uint16] = "SMALLINT"
type_dict[np.uint32] = "INT"
type_dict[np.uint64] = "BIGINT"
type_dict[np.int8] = "TINYINT"
type_dict[np.int16] = "SMALLINT"
type_dict[np.int32] = "INT"
type_dict[np.int64] = "BIGINT"
type_dict[np.float16] = "FLOAT"
type_dict[np.float32] = "FLOAT"
type_dict[np.float64] = "DOUBLE PRECISION"
type_dict[np.str] = "TEXT"
type_dict[np.string_] = "TEXT"
type_dict[str] = "TEXT"
