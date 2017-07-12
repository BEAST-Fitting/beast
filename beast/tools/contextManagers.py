from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import time
import math


class timeit(object):
    """ Time a block of your code.
    Usage:
        with timeit(text):
            <code>
    print how long took the execution of this code.
    """
    def __init__(self, text=None):
        self.text = text or ''

    def __enter__(self):
        print("Timing %s" % (self.text))
        self.start = time.time()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__pretty_print(time.time() - self.start)

    def __pretty_print(self, t):
        units = ["s", "ms", 'us', "ns"]
        scaling = [1, 1e3, 1e6, 1e9]
        if t > 0.0 and t < 1000.0:
            order = min(-int(math.floor(math.log10(t)) // 3), 3)
        elif t >= 1000.0:
            order = 0
        else:
            order = 3

        print("%s Execution time: %.3g %s" % (self.text, t * scaling[order], units[order]))
