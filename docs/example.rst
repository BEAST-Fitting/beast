###############
Running Example
###############

You can find examples of BEAST runs in the
`beast-examples repository <https://github.com/BEAST-Fitting/beast-examples>`_

Inside each example, there is a run_beast*.py script.

phat_small Example
------------------

This example is based on a *very* small amount of old PHAT data.

If you installed Python through AstroConda, first activate the correct
AstroConda environment

.. code-block:: console

    $ source activate astroconda

Verify that the current default Python is version 3

.. code-block:: console

    $ python --version

Next, bring up the BEAST help message, which describes the available switch
options, with

.. code-block:: console

    $ ./run_beast.py -h

You should be presented with the following options::

  -h, --help              show this help message and exit
  -p, --physicsmodel      Generate the model grid
  -o, --observationmodel  Calculate the noise model
  -t, --trim              Trim the model and noise grids
  -f, --fit               Fit the observed data
  -r, --resume            Resume a run

Now launch a sample BEAST run (with flags set to run through the full
sequence of generation of physics model, observation model generation, trimming
of the grid, and fitting to the observed data) using

.. code-block:: console

  $ ./run_beast.py -potf

If the BEAST is running correctly, this command should run without errors
and should have written the output files into 'beast_example_phat/'. The result
can be plotted using

.. code-block:: console

    $ python beast/plotting/plot_indiv_fit.py beast_example_phat/beast_example_phat

The argument for this script is the prefix of the output files. The output
should look like this

.. image:: beast_example_phat_ifit_starnum_0.png
