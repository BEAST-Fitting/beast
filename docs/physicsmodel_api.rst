
Physics Model
*************

Stars
=====

.. automodapi:: beast.physicsmodel.stars.stellib

.. automodapi:: beast.physicsmodel.stars.isochrone

Dust
====

.. automodapi:: beast.physicsmodel.dust.extinction

.. automodapi:: beast.physicsmodel.dust.extinction_extension

Weights
=======

Priors are implemented as weights to allow for fast integration.
Handling the grid spacing is done with grid weights allowing the priors values
to be independent of the grid spacing.

.. automodapi:: beast.physicsmodel.priormodel

.. automodapi:: beast.physicsmodel.grid_weights_stars

.. automodapi:: beast.physicsmodel.grid_and_prior_weights

Grid
====

.. automodapi:: beast.physicsmodel.grid

.. automodapi:: beast.physicsmodel.helpers.gridbackends

.. automodapi:: beast.physicsmodel.creategrid
