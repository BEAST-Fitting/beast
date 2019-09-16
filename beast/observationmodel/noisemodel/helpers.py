import numpy as np

__all__ = ["convert_dict_to_structured_ndarray"]


def convert_dict_to_structured_ndarray(data):
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
    for key, dk in data.items():
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
    tab = np.rec.fromarrays(iter(data.values()), dtype=newdtype)
    return tab
