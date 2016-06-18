from . import table
import numpy as np
from StringIO import StringIO

t = table.Table()
t.add_column('a', np.arange(10, dtype=float))
t.add_empty_column('t2', np.int8, shape=(10,))
t.add_empty_column('t3', np.int8, shape=(10,))
t.add_empty_column('t4', np.float16, shape=(10,), format='0.3f')
t.set_alias('t1','a')


f = StringIO()

t.write(f, 'txt')

f.seek(0)

print f.read()
del f
