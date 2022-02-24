##############
Append filters
##############

Users can add any filters into their local ``beast/libs/filters.hd5``. 
The list of default filters in the BEAST filter library can be found `<here https://beast.readthedocs.io/en/latest/beast_setup.html#beast-filters>`_.

If you need to include fluxes in particular bands in your SED fitting or BEAST
predictions but those particular bands do not exist in the BEAST filter library,
you should append those filters to the BEAST filter library before generating 
your physics model grid. To do this, you need an ascii file that has two columns:
wavelegnth in angstrom and thoughput. The BEAST filter library requires basic information
about the filter. This includes the tablename for the filter information (e.g., HST_WFC3_F275W), 
observatorty associated with the filter (e.g., HST), instrument of the filter (e.g., WFC3), 
and name of the filter (e.g., F275W). 

Command to add a filter to the BEAST filter library.

  .. code-block:: console

     $ python -m beast.observationmodel.phot filter_curve_file --tablename tablename
       --observatory observatory --instrument instrument --filtername filtername
