import numpy as np
from astropy.table import Table, Column

input_column = 'value'
bin_colname = 'bin'


class DensityMap:
    """
    Class which helps with using a density map consistently, and allows
    for reliable writing and reading. The file format in question is the
    one produced by create_background_density_map. It should readable as
    an astropy Table, where every entry represents a tile of the
    background map, with the following tile properties as columns:
    'i_ra': right ascencion (RA) index
    'i_dec': declination (DEC) index
    'value': the density value
    'min_ra': edge of tile with lowest RA
    'max_ra': edge of tile with highest RA
    'min_dec': edge of tile with lowest DEC
    'max_dec': edge of tile with highest DEC
    and the following metadata
    table.meta['ra_grid']
    table.meta['dec_grid']
    containing the edges of the RA and DEC bins.
    For the time being, just look at the source code of
    create_background_density_map to see how this is constructed exactly.
    """

    def __init__(self, tile_data):
        """
        density_map: Table
            of the format described above (extra columns are allowed)
        --or--
        density_map: string
            path to file containing this table
        """
        if isinstance(tile_data, str):
            self.tile_data = Table.read(tile_data)
        else:
            self.tile_data = tile_data

        self.ra_grid = self.tile_data.meta['ra_grid']
        self.dec_grid = self.tile_data.meta['dec_grid']

        self.min_i_ra = min(self.tile_data['i_ra'])
        self.max_i_ra = max(self.tile_data['i_ra'])
        self.min_i_dec = min(self.tile_data['i_dec'])
        self.max_i_dec = max(self.tile_data['i_dec'])

        # map index pairs to table rows
        self.tile_for_ij = {}
        for r in range(len(self.tile_data)):
            ij_ra_dec = (self.tile_data[r]['i_ra'], self.tile_data[r]['i_dec'])
            self.tile_for_ij[ij_ra_dec] = r

    def write(self, fname):
        """
        Write this map to file fname (in fits format)
        """
        self.tile_data.write(fname, format='hdf5', path='tile_data',
                             overwrite=True)

    def tile_for_position(self, ra, dec):
        """
        Finds which tile a certain ra,dec fits into
        """
        # Get index pair
        i_ra = np.searchsorted(self.ra_grid[:-1], ra, side='right') - 1
        i_ra = max(i_ra, self.min_i_ra)
        i_ra = min(i_ra, self.max_i_ra)

        i_dec = np.searchsorted(self.dec_grid[:-1], dec, side='right') - 1
        i_dec = max(i_dec, self.min_i_ra)
        i_dec = min(i_dec, self.max_i_ra)

        # Use index pair to row index map
        return self.tile_for_ij[(i_ra, i_dec)]

    def min_ras_decs(self):
        """
        Return a tuple, containing the coordinates of the bottom left
        corner for all the tiles (RAs, DECs)
        """
        return self.tile_data['min_ra'], self.tile_data['min_dec']

    def delta_ras_decs(self):
        """
        Return a tuple, containing the widths (RA) and heights (DEC) of
        each tile
        """
        return (self.tile_data['max_ra'] - self.tile_data['min_ra'],
                self.tile_data['max_dec'] - self.tile_data['min_dec'])

    def value(self, tile_index):
        """
        Return the map value at the given tile index
        """
        return self.tile_data[input_column][tile_index]

    def tile_vals(self):
        """
        Return all the values of the tiles
        """
        return self.tile_data[input_column]


class BinnedDensityMap(DensityMap):
    """
    Subclass which adds an extra column, which groups tiles by density
    bin. It is recommended to not use the constructor directly. In the
    'create' function, you can choose how many density bins you want,
    while the 'read' function reads an existing binned density map from
    file.
    """

    def __init__(self, tile_data, bins=None):
        DensityMap.__init__(self, tile_data)

        if bins is None:
            # Check if the bins are already there, and return
            if bin_colname not in self.tile_data.colnames:
                raise Exception('{} column not yet calculated. '
                                'Please use \'create\' function instead'.format(bin_colname))
            return

        if bin_colname in self.tile_data.colnames:
            print('{} column already there, overwriting it.'.format(bin_colname))
            self.tile_data[bin_colname] = bins
        else:
            c = Column(name=bin_colname, data=bins)
            self.tile_data.add_column(c)

        self.bin_indices_used = np.sort(np.unique(bins))

    def create(density_map, N_bins=None):
        """
        Creates a binned density map from a DensityMap file, or from an
        astropy table loaded from it. The tiles are grouped into
        N_bins density bins.
        If N_bins is none, each tile is treated as a separate bin.
        """
        # Use the base class to decide what to do with density_map (can
        # be file or table object)
        binned_density_map = DensityMap(density_map)

        # Create the extra column here
        if N_bins is None:
            bins = np.array(range(len(binned_density_map.tile_data)))

        else:
            # Create the density bins
            # [min, ., ., ., max]
            tile_densities = binned_density_map.tile_data[input_column]
            min_density = np.amin(tile_densities)
            max_density = np.amax(tile_densities)
            bin_edges = np.linspace(min_density - 0.01 * abs(min_density),
                                    max_density + 0.01 * abs(max_density),
                                    N_bins + 1)

        # Find which bin each tile belongs to
        # e.g. one of these numbers: 0 [1, 2, 3, 4, 5] 6
        # We have purposely chosen our bin boundaries so that no points fall
        # outside (or on the edge) of the [1,5] range
        bins = np.digitize(
            binned_density_map.tile_data[input_column], bin_edges)

        # Upgrade to this subclass, and return
        return BinnedDensityMap(binned_density_map.tile_data, bins)

    def read(density_map_fname):
        return BinnedDensityMap(density_map_fname)

    def bin_indices_used(self):
        return self.tile_data[bin_colname]

    def bin_for_position(self, ra, dec):
        """
        Finds which density bin a certain ra,dec fits into, and
        returns its index.
        """
        t = self.tile_for_position(ra, dec)
        return self.tile_data[bin_colname][t]

    def value_foreach_tile(self):
        return self.tile_data[input_column]

    def bin_foreach_tile(self):
        return self.tile_data[bin_colname]

    def tiles_foreach_bin(self):
        """
        Invert the above function. Result is a list of lists. For each
        bin, the tiles that where grouped into it are listed.
        """
        b_per_tile = self.bin_foreach_tile()
        return [np.nonzero(b_per_tile == b)[0] for b in
                self.bin_indices_used]
