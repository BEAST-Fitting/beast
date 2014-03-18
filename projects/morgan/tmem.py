""" Testing memmap process for parallel processing of big data """
import multiprocessing
import numpy as np
import os
import sys

dtype = np.int32
shape = 1000 * 1000 * 1000
sz = 4 * shape
fn = 'tmp.memmap'

print 'PID', os.getpid()


def gen_data():
    if not os.path.exists(fn) or os.stat(fn).st_size != sz:
        print 'Writing file...'
        mf = np.memmap(fn, mode='w+', dtype=dtype, shape=shape)
        mf[:] = np.random.uniform(size=shape, high=2 ** 31).astype(dtype)
        # mf.close()
        del mf
    print 'done'


def task1(val):
    _mf = np.memmap(fn, mode='r', dtype=dtype, shape=shape)
    for i in range(1000 * 100):
        idx = np.random.randint(0, len(_mf), 1)
        val = int(_mf[idx])
    print '.',
    sys.stderr.flush()


gen_data()

print 'Group id:'
print os.getgid()
p = multiprocessing.Pool()
print 'running...'
p.map(task1, xrange(1000))
print 'done'
