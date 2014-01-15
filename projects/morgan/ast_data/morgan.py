import numpy as np
import pylab as plt

from beast.external.eztables import Table

import notebookrc
#make good figures
notebookrc.ezrc(16)

from faststats.faststats.plot import corrplot

ir_file = '12055_M31-B21-F11-IR_F110W_F160W_gst.fake.fits'
uv_file = '12055_M31-B21-F11-UVIS_F275W_F336W_gst.fake.fits'
vis_file = '12055_M31-B21-F11-WFC_F475W_F814W_gst.fake.fits'


def get_input_table(f):
    """get_input_table

    keywords
    --------

    f: str
        file to open

    returns
    -------

    tab: eztable.Table
        table instance
    """
    return Table(f)


def get_info(tab):
    """get_info - retrieve useful columns from the AST files

    Output information contains:
        MAG1IN      Input mag in filter 1
        MAG1OUT     Output mag in filter 1
        MAG2IN      Input mag in filter 1
        MAG2OUT     Output mag in filter 2
        MAG1_STD    Output uncertainty in filter 1
        MAG2_STD    Output uncertainty in filter 2
        BIAS1       Output - Input in filter 1
        BIAS2       Output - Input in filter 2

    keywords
    --------

    tab: eztable.Table instance
        table to extract from
        (can be replaced by a direct access at some point)

    returns
    -------
    res: ndarray
        array that contains the output array

    fields: sequence of strings
        string corresponding to the ordered content of the array

    names: sequence of strings
        strings used to make nice plot labels
    """

    #extract useful info
    fields = 'RA DEC MAG1IN MAG1OUT MAG2IN MAG2OUT MAG1_STD MAG2_STD'.split()
    res = np.asarray([ tab[k] for k in fields ])

    #add sub products
    sub = np.asarray([ res[1] - res[0], res[3] - res[2]])
    res = np.vstack([res, sub])

    fields += 'BIAS1 BIAS2'.split()
    names = r'RA DEC MAG1$_{IN}$ MAG1$_{OUT}$ MAG2$_{IN}$ MAG2$_{OUT}$ $\sigma_{MAG1}$ $\sigma_{MAG2}$ $\mu_{MAG1}$ $\mu_{MAG_2}$'.split()
    return res, fields, names


#tab = get_input_table(vis_file)
#res, fields, names = get_info(tab)
#del tab

## plot magnitude correlations
#plt.figure()
#plt.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.9)
#corrplot(res, names=names)
#plt.title('Magnitude correlations')

"""
Interestingly, mag1 et mag 2 IN are 100% correlated. This translate into a
high correlation in the output magnitudes and of course biases and
dispersions.
"""

## plotting triangle plot
#plt.figure()
#d = { nk: res[ek] for ek, nk in enumerate(fields) }
#notebookrc.plotCorr(d, d.keys(), lbls=names)
#plt.show()


