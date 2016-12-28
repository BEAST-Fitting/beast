""" Common helpers in the grid package """

__all__ = ['isNestedInstance', 'pretty_size_print']

def isNestedInstance(obj, cl):
    """ Test for sub-classes types
        I could not find a universal test

        keywords
        --------
        obj: object instance
            object to test

        cl: Class
            top level class to test

        returns
        -------
        r: bool
            True if obj is indeed an instance or subclass instance of cl
    """
    tree = []
    for k in cl.__subclasses__():
        tree += k.__subclasses__()
    tree += cl.__subclasses__() + [ cl ]
    return issubclass(obj.__class__, tuple(tree))


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
