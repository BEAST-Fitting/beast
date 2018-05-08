##############################
Guide to working with subgrids
##############################

Concept
=======

The idea of this approach is that we split up the model grid into
non-overlapping subgrids, and are then able to run each step of a full BEAST run
for each grid individually. The calculations for the individual subgrids can
then be run in parallel, for a speed boost without memory overhead, or
sequentially, to use less memory. Of course, a combination of the two is also
possible, splitting a grid into many subgrids, and running a couple of subgrid
calculations at the same time.

At the end of the calculation, we will possess partial PDFs and statistics of
each subgrid, which can be merged into a single 1dpdf file and a single stats
file. By taking into account the weights of each subgrid correctly, the
resulting file should be equivalent to the result of a BEAST run on the full
grid.

Workflow
========

To make use of this functionality, no extra data or changes to the datamodel
file are needed, but you will need a custom run script that makes use of a set
of newly implemented functions, and new options for existing functions. We will
now give a summary of what such a run script has to pay attention to, for each
of the steps in a BEAST run. (An example script might be provided later).

Most the new functions can be found in ``beast.tools.subgridding_tools``.

Please refer to the regular example code for the single grid implementation of
these steps (``beast/examples/phat_small``).

Physics model
-------------

First the spectral (stellar) grid is created, using ``make_iso_table``,
``make_spectral_grid`` and ``add_stellar_priors``. Then, the extinction
parameters are applied to this grid, and an extinguished SED grid is obtained,
using ``make_extinguished_sed_grid``.

The splitting of the grids has to happen somewhere in this function.
Technically, ``split_grid`` can be either after obtaining the spectral grid with
prior weights, or after obtaining the complete SED grid. The former makes more
sense however, because then ``make_extinguished_sed_grid`` can be run for
individual spectral subgrids, which avoids the memory impact of creating the
complete SED grid. This choice also allows the user to run the construction of
the grid in parallel.

.. tip:: The ``split_grid`` function returns the file names of the newly created
   subgrids. It is very useful to save these to a text file, so that they can be
   used in the other steps.

AST input list
--------------

This is the only step where the complete SED grid is needed. The subgrids can be
merged into a single file using ``merge_grids``. Just provide an output name,
and a list of file names pointing to all the subgrids. The rest of the AST input
list generation needs no changes once the full grid file is available.

Observations/Noise models
-------------------------

Here we will create separate noise model files, one for each subgrid. Nothing
special happens here, e.g. just call ``make_toothpick_noise_model`` for each
subgrid using the same AST results file, providing adequate output names for the
resulting noise models. It is safe to run this in parallel.

Trimming of the physics and noise models
----------------------------------------

The same as the above applies here. Just make sure that the
subgrid/subnoisemodel files are paired correctly.

Fitting & merging the results
-----------------------------

Compatibility
~~~~~~~~~~~~~

To make sure that the results of the fitting routine for the individual grids
are compatible, there are several subtleties which come into play here. Firstly,
it needs to be made sure that the 1dpdfs are compatible: their number of bins
and the values for the bin centers need to be exactly the same. To ensure this,
we need to fix three values for each quantity:

1) the minimum value
2) the maximum value
3) the number of unique values

This is why a new optional argument is provided in the main fitting function,
``summary_table_memory``, which allows the user to override the min, max and
number of unique values for all of the quantities.

The option is called ``grid_info_dict``, and needs to be a nested dictionary of
a certain format. ``subgridding_tools`` contains a function called
``reduce_grid_info`` which will generate this dictionary for you. Just provide
the filenames to all the (trimmed) subgrids and their (trimmed) noisemodels.

This dictionary has an entry like this for each quantity (``Rv`` in this example):

.. code:: python

   grid_info_dict['Rv'] = {'min': 0, 'max': 10, 'num_unique': 20}

Fit
~~~

When the info described above has been collected, you can start calling
``summary_table_memory`` for each of the subgrids, each time providing a trimmed
subgrid/trimmed subnoisemodel pair, and adequate filenames for the output. The
rest of the arguments can be identical the fit on each subgrid. However, be sure
to set ``do_not_normalize`` to ``True``, see note below.

Merge
~~~~~

When all the subgrid fits have been successfully completed, the merge step can be
started. To do this, just gather all the filenames for the pdf1d and stats
files, and pass them to ``merge_pdf1d_stats``.

.. note::

   The main fitting function needed to be modified so that the `Pmax` values
   that it stores (which are the maximum log likelihood, needed to calculate the
   `Best` values) are compatible between subgrids. This meant getting rid of
   some forms of normalization (specifically, the prior weight normalization
   needed to be disabled). Setting ``do_not_normalize`` should have no effect on
   the result actually, so we might remove this option altogether and make it
   the default behavior.

.. note::

   To calculate the expectation values, another modification to the same function
   has been done. It now stores a measure for the total weight of the subgrid,
   `total_log_norm`. This value is equal to ``log(sum(exp(lnp)))``, and is
   calculated by taking the log of the normalization factor used in the code
   (because ``sum(exp(lnp)) / normalization = 1``). By comparing this value
   between subgrids, we are able to calculate a weighted average for each
   expectation value, which should be close to the one that would be obtained
   by fitting over the whole grid at once.
