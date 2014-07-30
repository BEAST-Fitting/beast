import numpy as np
from scipy.spatial import cKDTree
import tables
# from functools import wraps

from beast.external.eztables import Table
from beast.core.grid import FileSEDGrid
#from beast.external.ezunits import unit
from beast.core.vega import Vega


# def generator(func):
#     """ Allow to clearly read when a function is a generator in the code
#     Do nothing more.
#     """
#     @wraps(func)
#     def wrap(*args, **kwargs):
#         return func(*args, **kwargs)
#     return wrap


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
    for key, dk in data.iteritems():
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
    tab = np.rec.fromarrays(data.itervalues(), dtype=newdtype)
    return tab


def _prepare_x(x):
        """ make sure the data is correctly presented """
        xi = np.asarray(x)
        shape = xi.shape

        if len(shape) < 1:
            raise ValueError('single value passed')

        if len(shape) == 1:
            return xi[..., np.newaxis]
        else:
            return xi


def nearest_neighbors(x, k=10,eps=0.):
    """ kd-tree for quick nearest-neighbor lookup

    ..note::
        using the default leafsize=10

    Parameters
    ----------
    x : array_like, last dimension self.m
        An array of points to query.

    k : integer
        The number of nearest neighbors to return.

    eps : non-negative float
        Return approximate nearest neighbors; the kth returned value
        is guaranteed to be no further than (1+eps) times the
        distance to the real k-th nearest neighbor.

    p : float, 1<=p<=infinity, default=2
        Which Minkowski p-norm to use.
        1 is the sum-of-absolute-values "Manhattan" distance
        2 is the usual Euclidean distance
        infinity is the maximum-coordinate-difference distance

    distance_upper_bound : nonnegative float, default=np.inf
        Return only neighbors within this distance.  This is used to prune
        tree searches, so if you are doing a series of nearest-neighbor
        queries, it may help to supply the distance to the nearest neighbor
        of the most recent point.

    Returns
    -------
    i : ndarray of ints
        The locations of the neighbors in self.data.
        If `x` has shape tuple+(self.m,), then `i` has shape tuple+(k,).
        Missing neighbors are indicated with self.n.
    """
    tree = cKDTree(x)
    d, ind = tree.query(x, k=k, eps=eps)
    return ind


def get_info(tab):
    """retrieve useful columns from the AST files

    Output information contains:
        MAG1IN      Input mag in filter 1
        MAG1OUT     Output mag in filter 1
        MAG2IN      Input mag in filter 1
        MAG2OUT     Output mag in filter 2
        MAG1_STD    Output uncertainty in filter 1
        MAG2_STD    Output uncertainty in filter 2
        BIAS1       Output - Input in filter 1
        BIAS2       Output - Input in filter 2

    Parameters
    ----------

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
    sub = np.asarray([ res[2] - res[3], res[4] - res[5]])
    res = np.vstack([res, sub])

    fields += 'BIAS1 BIAS2'.split()
    names = 'RA DEC MAG1$_{IN}$ MAG1$_{OUT}$ MAG2$_{IN}$ MAG2$_{OUT}$ $\sigma_{MAG1}$ $\sigma_{MAG2}$ $\mu_{MAG1}$ $\mu_{MAG_2}$'.split()
    return res, fields, names


# def get_AST_dict(files):
#     """
#     Returns a separate dictionnary containing AST data for each camera
#
#     Parameters
#     ----------
#     files: sequence
#         expecting one file per camera
#
#     Returns
#     -------
#     lst: sequence(dict)
#         sequence of dictionaries (one per input file)
#     """
#
#     res = []
#     for fk in files:
#         tab = Table(fk)
#         res, fields, names = get_info(tab)
#         res.append(dict([(nk, res[ek]) for ek, nk in enumerate(fields)]))
#
#     return res


# @generator
# def gen_AST_dict(files):
#     """
#     Returns a separate dictionnary containing AST data for each input file
#
#     Parameters
#     ----------
#     files: sequence
#         expecting one file per camera
#
#     Returns
#     -------
#     lst: sequence(dict)
#         sequence of dictionaries (one per input file)
#     """
#     res = []
#     for fk in files:
#         tab = Table(fk)
#         res, fields, names = get_info(tab)
#         yield dict([(nk, res[ek]) for ek, nk in enumerate(fields)])


#lazy functions
def toFlux(mag):
    return 10 ** (-0.4 * mag)


# def toMag(flux):
#     return -2.5 * np.log10(flux)


def toFlux2(mag, band, distanceModulus=24.47):
    """ convert vega magnitude to apparent flux at a given distance """

    filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W',
               'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']

    match = False
    for e, k in enumerate(filters):
        if band.lower() == k.split('_')[-1].lower():
            match = True
            break

    if not match:
        raise('Band {0:s} not found in filters'.format(band))

    with Vega() as v:
        _, bb, _ = v.getFlux([filters[e]])
    return bb[0] * 10 ** (0.4 * (distanceModulus - mag))


def compute_stddev(mag_in, mag_out, k=10, eps=0,
                   completeness_mag_cut=80, name_prefix=None):
    """
    Computes standard deviation and store the result in a dictionary

    Parameters
    ----------
    mag_in: ndarray
        input magnitudes

    mag_out: ndarray
        output magnitudes

    completeness_mag_cut: float
        magnitude at which consider a star not recovered

    k: Integer
        Number of nearest neighbors taken in the standard deviation computation

    eps: non-negative float
        precision on the NN search

    name_prefix: str
        if set, all output names in the final structure will start with this
        prefix.

    Returns
    -------
    d: np.recarray
        named array containing the statistics


    Method
    ------
    Statistics are computed for each input artificial star.
    For each input star, we find the k-NN in input magnitude space and compute
    their mean and variance also in magnitude and completeness fractions.

    Finally, we return the flux conversion of the statistics: bias and stddev
    dispersion per input star averaged over k-NN.

    """
    n_asts = len(mag_in)

    if name_prefix is None:
        name_prefix = ''
    else:
        if name_prefix[-1] != '_':
            name_prefix += '_'

    mag_in_prep = _prepare_x(mag_in)
    NN_mag = nearest_neighbors(mag_in_prep, k=k, eps=eps)

    std_mag = np.zeros(n_asts, dtype=float)
    completeness = np.zeros(n_asts, dtype=float)

    for i in range(n_asts):
        nn_i = mag_out[NN_mag[i]]
        std_mag[i] = np.sqrt( (1. / k) * ( (nn_i - mag_in[i]) ** 2).sum(axis=0))
        completeness[i] = np.sum(nn_i < completeness_mag_cut, axis=0)

    completeness /= float(k)

    mag_inup = mag_in + std_mag
    mag_indown = mag_in - std_mag

    f_in = toFlux(mag_in)
    f_out = toFlux(mag_out)
    f_up = toFlux(mag_indown)
    f_down = toFlux(mag_inup)
    std_flux = 0.5 * (f_up - f_down)

    d = {name_prefix + 'MAG_STD': std_mag,
         name_prefix + 'MAG_BIAS': mag_out - mag_in,
         name_prefix + 'MAG_IN': mag_in,
         name_prefix + 'MAG_OUT': mag_out,
         name_prefix + 'FLUX_STD': std_flux,
         name_prefix + 'FLUX_BIAS': f_out - f_in,
         name_prefix + 'FLUX_IN': f_in,
         name_prefix + 'FLUX_OUT': f_out,
         name_prefix + 'COMPLETENESS': completeness}

    return convert_dict_to_structured_ndarray(d)


def toothpick_ast_interp(D, sedgrid):
    """
    Interpolate the results of the ASTs on a sed model grid

    Parameters
    ----------
    D: dictionary
        output of extract_bias_std

    sedgrid: beast.core.grid type
        model grid to interpolate AST results on

    Returns
    -------
    (bias,sigma): tuple
        Two 1d numpy arrays
    """

    N = sedgrid.grid.nrows
    filters = sedgrid.filters
    #filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W',
    #           'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
    M = len(filters)

    flux = sedgrid.seds
    with Vega() as v:
        names,vega_flux,lamb = v.getFlux(filters)

    bias = np.empty((N,M), dtype=float)
    sigma = np.empty((N,M), dtype=float)

    keys = ['uv', 'opt', 'ir']

    for i in range(M):
        ii = (i % 2) + 1
        key = keys[(i // 2)]

        var1 = D[key]['F{0:d}IN'.format(ii)]
        var2 = D[key]['BIAS{0:d}_FLUX'.format(ii)]
        var3 = D[key]['STD{0:d}_FLUX'.format(ii)]

        arg_sort = np.argsort(var1)

        bias[:, i] = np.interp(flux[:, i], var1[arg_sort], var2[arg_sort])
        sigma[:, i] = np.interp(flux[:, i], var1[arg_sort], var3[arg_sort])

    return (bias, sigma)


# absflux calibration covariance matrix for NGC4214 specific filters
ngc4214_calibration_covariance = [1.80e-4,  1.37e-4, 6.02e-5, 2.44e-5, 1.23e-6, -4.21e-6,
                                  1.37e-4,  1.09e-4, 5.72e-5, 3.23e-5, 1.65e-5, 1.32e-5,
                                  6.02e-5,  5.72e-5, 5.07e-5, 4.66e-5, 4.40e-5, 4.34e-5,
                                  2.44e-5,  3.23e-5, 4.66e-5, 5.42e-5, 5.87e-5, 5.99e-5,
                                  1.23e-6,  1.65e-5, 4.40e-5, 5.87e-5, 6.98e-5, 7.33e-5,
                                  -4.21e-6, 1.32e-5, 4.34e-5, 5.99e-5, 7.33e-5, 7.81e-5 ]

# PHAT values:
# abs_calib = np.array([1.19,1.04,0.71,0.74,0.84,0.88]) * 0.01


def make_toothpick_noise_model(outname, astfile, sedgrid, covariance=None, **kwargs):
    """ toothpick noise model assumes that every filter is independent with
    any other.
    """
    #outname = dir_project + 'noise_model_b%s_%s.hd5' % (brick,subdivision)

    # read mag_in, mag_out
    model = ToothPick_MultiFilterASTs(astfile, sedgrid.filters)
    # compute k-NN statistics: bias, stddev, completeness
    model.compute_knn_statistics(k=10, eps=0, completeness_mag_cut=80)
    # evaluate the noise model for all the models in sedgrid
    bias, sigma, compl = model(sedgrid)
    # save to disk/mem

    # absolute flux calibration uncertainties
    if covariance is not None:
        if covariance.ndim == 1:
            abs_calib_2 = covariance[:] ** 2
        else:   # assumes a cov matrix
            abs_calib_2 = np.diag(covariance)

        noise = np.sqrt(abs_calib_2 * sedgrid.seds ** 2 + sigma ** 2)
    else:
        noise = sigma

    with tables.openFile(outname, 'w') as outfile:
        outfile.createArray(outfile.root,'bias', bias)
        outfile.createArray(outfile.root,'error', noise)
        outfile.createArray(outfile.root,'completeness', compl)


class NoiseModel(object):
    def __init__(self, astfile, *args, **kwargs):
        self.astfile = astfile
        self.load_data()

    def load_data(self):
        self.data = Table(self.astfile)


class ToothPick_MultiFilterASTs(NoiseModel):
    def __init__(self, astfile, filters, *args, **kwargs):
        NoiseModel.__init__(self, astfile, *args, **kwargs)
        self.setFilters(filters)
        self.set_data_mappings()

        self._fluxes = None
        self._biases = None
        self._sigmas = None
        self._compls = None

    def set_data_mappings(self):
        """ hard code mapping directly with the interface to PHAT-like ASTs

        .. note::

            it makes it trivial to update this function for other input formats
        """
        try:
            for k in self.filters:
                #self.data.set_alias(k, k.split('_')[-1].lower() + '_rate')
                self.data.set_alias(k + '_out', k.split('_')[-1].lower() + '_vega')
                self.data.set_alias(k + '_in', k.split('_')[-1].lower() + '_in')
        except Exception as e:
            print(e)
            print('Warning: Mapping failed. This could lead to wrong results')

    def compute_knn_statistics(self, k=10, eps=0, completeness_mag_cut=80):
        """

        .. see also: :func:`compute_stddev`
        """

        shape = len(self.data), len(self.filters)

        self._biases = np.empty( shape, dtype=float)
        self._sigmas = np.empty( shape, dtype=float)
        self._compls = np.empty( shape, dtype=float)

        for e, k in enumerate(self.filters):

            mag_in = self.data[k + '_in']
            mag_out = self.data[k + '_out']

            d = compute_stddev(mag_in, mag_out, k=k, eps=eps,
                               completeness_mag_cut=completeness_mag_cut)

            self._fluxes[:, e] = d['FLUX_IN'] / self.vega_flux[e]
            self._sigmas[:, e] = d['FLUX_STD'] / self.vega_flux[e]
            self._biases[:, e] = d['FLUX_BIAS'] / self.vega_flux[e]
            self._compls[:, e] = d['COMPLETENESS']
            del d

    def setFilters(self, filters):
        """ set the filters and update the vega reference for the conversions

        Parameters
        ----------
        filters: sequence
            list of filters using the internally normalized namings
        """
        self.filters = filters

        #Data "rates" are normalized to Vega already, fits are not using vega

        # for optimization purpose: pre-compute
        #   getting vega mags, require to open and read the content of one file.
        #   since getObs, calls getFlux, for each star you need to do this expensive
        #   op.
        with Vega() as v:
            # name, vega_flux, lamb = v.getFlux(filters)
            _, vega_flux, _ = v.getFlux(filters)

        self.vega_flux = vega_flux

    def interpolate(self, sedgrid):
        """
        Interpolate the results of the ASTs on a model grid

        Parameters
        ----------
        sedgrid: beast.core.grid type
            model grid to interpolate AST results on

        Returns
        -------
        bias: ndarray
            bias table of the models

        sigma: ndarray
            dispersion table of the models

        comp: ndarray
            completeness table per model
        """
        flux = sedgrid.seds
        N, M = flux.shape

        if M != len(self.filter):
            raise AttributeError('the grid of models does not seem to be defined with the same number of filters')

        bias = np.empty((N, M), dtype=float)
        sigma = np.empty((N, M), dtype=float)
        compl = np.empty((N, M), dtype=float)

        for i in range(M):

            _fluxes = self._fluxes[:, i]
            arg_sort = np.argsort(self._fluxes[:, i])
            _fluxes = _fluxes[arg_sort]

            bias[:, i] = np.interp(flux[:, i], _fluxes, self._biases[arg_sort, i] )
            sigma[:, i] = np.interp(flux[:, i], _fluxes, self._sigmas[arg_sort, i])
            compl[:, i] = np.interp(flux[:, i], _fluxes, self._compls[arg_sort, i])

        return (bias, sigma, compl)


class ToothPick_perCamereASTs(NoiseModel):
    def __init__(self, astfiles, filters, *args, **kwargs):

        NoiseModel.__init__(self, astfiles, *args, **kwargs)
        self.models = [ ToothPick_MultiFilterASTs(astfile) for astfile in astfiles ]
        self.setFilters(filters)

    def set_data_mappings(self):
        """ hard code mapping directly with the interface to PHAT-like ASTs

        .. note::

            it makes it trivial to update this function for other input formats
        """
        #TODO: update the mapping to stick to the initial PHAT version
        pass

    def compute_knn_statistics(self, k=10, eps=0, completeness_mag_cut=80):
        """

        .. see also: :func:`compute_stddev`
        """
        for model in self.models:
            model.compute_knn_statistics(k=k, eps=eps,
                                         completeness_mag_cut=completeness_mag_cut)

    def setFilters(self, filters):
        """ set the filters and update the vega reference for the conversions

        Parameters
        ----------
        filters: sequence
            list of filters using the internally normalized namings
        """
        self.filters = filters
        for model in self.models:
            model.setFilters(filters)

    @property
    def vega_flux(self):
        return self.models[0].vega_flux


    def interpolate(self, sedgrid):
        """
        Interpolate the results of the ASTs on a model grid

        Parameters
        ----------
        sedgrid: beast.core.grid type
            model grid to interpolate AST results on

        Returns
        -------
        bias: ndarray
            bias table of the models

        sigma: ndarray
            dispersion table of the models

        comp: ndarray
            completeness table per model
        """
        flux = sedgrid.seds
        N, M = flux.shape

        if M != len(self.filter):
            raise AttributeError('the grid of models does not seem to be defined with the same number of filters')

        bias = np.empty((N, M), dtype=float)
        sigma = np.empty((N, M), dtype=float)
        compl = np.empty((N, M), dtype=float)

        for i in range(M):

            _fluxes = self._fluxes[:, i]
            arg_sort = np.argsort(self._fluxes[:, i])
            _fluxes = _fluxes[arg_sort]

            bias[:, i] = np.interp(flux[:, i], _fluxes, self._biases[arg_sort, i] )
            sigma[:, i] = np.interp(flux[:, i], _fluxes, self._sigmas[arg_sort, i])
            compl[:, i] = np.interp(flux[:, i], _fluxes, self._compls[arg_sort, i])

        return (bias, sigma, compl)


def ast_interp(D, sedgrid):
    """
    Interpolate the results of the ASTs on a model grid

    Parameters
    ----------
    D: dictionary
        output of extract_bias_std

    sedgrid: beast.core.grid type
        model grid to interpolate AST results on

    Returns
    -------
    (bias,sigma): tuple
        Two 1d numpy arrays
    """

    N = sedgrid.grid.nrows
    filters = sedgrid.filters
    #filters = ['HST_WFC3_F275W', 'HST_WFC3_F336W', 'HST_ACS_WFC_F475W',
    #           'HST_ACS_WFC_F814W', 'HST_WFC3_F110W', 'HST_WFC3_F160W']
    M = len(filters)

    flux = sedgrid.seds
    with Vega() as v:
        names,vega_flux,lamb = v.getFlux(filters)

    bias = np.empty((N,M), dtype=float)
    sigma = np.empty((N,M), dtype=float)

    keys = ['uv', 'opt', 'ir']

    for i in range(M):
        ii = (i % 2) + 1
        key = keys[(i // 2)]

        var1 = D[key]['F{0:d}IN'.format(ii)]
        var2 = D[key]['BIAS{0:d}_FLUX'.format(ii)]
        var3 = D[key]['STD{0:d}_FLUX'.format(ii)]

        arg_sort = np.argsort(var1)

        bias[:, i] = np.interp(flux[:, i], var1[arg_sort], var2[arg_sort])
        sigma[:, i] = np.interp(flux[:, i], var1[arg_sort], var3[arg_sort])

    return (bias, sigma)

if __name__ == '__main__':
    brick = '15'       # Brick index
    subdivision = '15'  # Index of the subregion in brick

    project = 'b%s_%s' % (brick, subdivision)

    # Directory where to write the noise model file =  Project directory
    dir_project = '/astro/dust_kg/harab/beast/projects/prod/%s/' % project
    # Model Grid
    sedgrid = FileSEDGrid( dir_project + '/' + project + '_seds.grid.hd5')

    # absolute flux calibration uncertainties
    abs_calib = [1.19,1.04,0.71,0.74,0.84,0.88]
    abs_calib_flux = sedgrid.seds * abs_calib * 0.01

    D = extract_bias_std(brick, subdivision)
    (bias, sig) = ast_interp(D, sedgrid)

    noise = np.sqrt(abs_calib_flux ** 2 + sig ** 2)

    outname = dir_project + 'noise_model_b%s_%s.hd5' % (brick,subdivision)

    with tables.openFile(outname, 'w') as outfile:
        outfile.createArray(outfile.root,'bias',bias)
        outfile.createArray(outfile.root,'error',noise)
