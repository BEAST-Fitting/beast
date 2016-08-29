from knn import KDTreeInterpolator
import numpy as np
import pylab as plt

from beast.external.eztables import Table

import notebookrc
#make good figures
notebookrc.ezrc(16)

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


def NDinterp(x, y):
    """ generate an interpolator but taking care of non-finite values """
    ind = np.isfinite(y)
    inter = KDTreeInterpolator(x[ind], y[ind])
    return inter


#lazy functions

def toFlux(mag):
    return 10 ** (-0.4 * mag)


def toMag(flux):
    return -2.5 * np.log10(flux)


if __name__ == '__main__':
    tab = get_input_table(vis_file)
    res, fields, names = get_info(tab)
    d = { nk: res[ek] for ek, nk in enumerate(fields) }
    del tab

    # get biases in flux instead of mags
    F1IN = toFlux(d['MAG1IN'])
    F1OUT = toFlux(d['MAG1OUT'])
    B1F = F1IN - F1OUT
    B1mag = toMag(B1F)

    F2IN = toFlux(d['MAG2IN'])
    F2OUT = toFlux(d['MAG2OUT'])
    B2F = F2IN - F2OUT
    B2mag = toMag(B2F)

    fit_npts = 400

    # in mag first

    #make an array of points (Npoints, Ndims)
    data = np.asarray([d['MAG1IN'], d['MAG2IN']]).T
    fn1 = NDinterp(data, d['BIAS1'])
    fn2 = NDinterp(data, d['BIAS2'])

    fdata = np.asarray([F1IN, F2IN]).T
    ffn1 = NDinterp(fdata, B1F)
    ffn2 = NDinterp(fdata, B2F)

    #x1 = np.linspace(d['MAG1IN'][ind1].min(), d['MAG1IN'][ind1].max(), fit_npts)
    #x2 = np.linspace(d['MAG2IN'][ind2].min(), d['MAG2IN'][ind2].max(), fit_npts)
    #fx1 = np.exp(np.linspace(np.log(F1IN[ind1].min()), np.log(F1IN[ind1].max()), fit_npts))
    #fx2 = np.exp(np.linspace(np.log(F2IN[ind2].min()), np.log(F2IN[ind2].max()), fit_npts))

    #mdata = np.asarray( [ k for k in np.nditer(np.ix_(x1, x2))] )

    # take input data in mags and add noise
    # noise is set to 1 mag dispersion
    mdata = data + np.random.normal(0, 1., data.shape)

    plt.figure(figsize=(20, 16))

    ax = plt.subplot(221)
    ax.plot(d['MAG1IN'], d['BIAS1'], ',', alpha=0.2)
    ax.plot(mdata[:, 0], fn1(mdata), 'k,', lw=2)
    ax.set_ylim(-3, 2)
    notebookrc.hide_axis(['top', 'right'], ax=ax)
    ax.set_xlabel('MAG1$_{IN}$')
    ax.set_ylabel('$\mu_{MAG1}$')

    ax = plt.subplot(222)
    ax.plot(d['MAG2IN'], d['BIAS2'], ',', alpha=0.2)
    ax.plot(mdata[:, 1], fn2(mdata), 'k,', lw=2)
    ax.set_ylim(-3, 2)
    notebookrc.hide_axis(['top', 'right'], ax=ax)
    ax.set_xlabel('MAG2$_{IN}$')
    ax.set_ylabel('$\mu_{MAG2}$')

    # in flux
    #mdata = np.asarray( [ k for k in np.nditer(np.ix_(x1, x2))] )
    mdata = fdata * np.random.normal(1., 2., fdata.shape)

    ax = plt.subplot(223)
    ax.loglog(F1IN, B1F, ',', alpha=0.2)
    ax.loglog(mdata[:, 0], ffn1(mdata), 'k,')
    notebookrc.hide_axis(['top', 'right'], ax=ax)
    ax.set_xlabel('F1$_{IN}$')
    ax.set_ylabel("$\mu'_{F1}$")

    ax = plt.subplot(224)
    ax.loglog(F2IN, B2F, ',', alpha=0.2)
    ax.loglog(mdata[:, 1], ffn2(mdata), 'k,')
    notebookrc.hide_axis(['top', 'right'], ax=ax)
    ax.set_xlabel('F2$_{IN}$')
    ax.set_ylabel("$\mu'_{F2}$")
