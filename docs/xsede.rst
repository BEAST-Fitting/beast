.. _beast_xsede:

##########################
Running the BEAST on XSEDE
##########################

(Before reading this section, be sure you're familiar with the BEAST
:ref:`production run workflow<beast_standard_workflow>`)

Running the BEAST with a finely-spaced grid requires considerable computational
resources, so you may choose to use `XSEDE <https://www.xsede.org/>`_.  This
page gives an overview of running the BEAST on XSEDE based on the team's
experience with METAL.  It includes applying for an allocation, using the
`slurm` queue system, and documentation for the `XSEDE BEAST wrapper
<https://github.com/BEAST-Fitting/beast-examples/tree/master/metal_xsede>`_
in `beast-examples <https://github.com/BEAST-Fitting/beast-examples>`_.

The XSEDE online `documentation <https://portal.xsede.org/documentation-overview>`_
is quite extensive, and their help desk is very helpful and responsive.  Note
that XSEDE also periodically runs free online workshops for different topics,
several of which BEAST team members have attended.


*****************
XSEDE Allocations
*****************

Very broadly, these are the steps you follow to use XSEDE resources:

* Get a `startup allocation <https://portal.xsede.org/allocations/startup>`_.
  There's a convenient request form that only requires a short justification
  (both for the science and to explain why you don't currently have access to
  sufficient resources).
* Run the BEAST on enough of your data to get a good estimate of the resources
  you'll need for a full production run.  Though if you're only doing a few
  fields, the startup allocation may be enough for your needs!
* Submit a proposal for a `research allocation <https://portal.xsede.org/allocations/research>`_.
  Proposals are accepted every 3 months.  Be sure to carefully read the
  proposal requirements and/or watch the webinar, because it's not always clear
  what documents are required for what proposal types (if in doubt, ask the
  helpdesk!).  You're welcome to reference the `METAL proposal
  <https://www.overleaf.com/read/ysmvjxbbrtvf>`_.

For METAL, we used a combination of Bridges Regular and Bridges Large.

* Bridges Regular: Charges by CPU usage (e.g., using 5 CPUs for 3 hours charges
  15 CPU-hours).  Each CPU comes with 4.5GB of memory.
* Bridges Large: Charges by memory usage (e.g., using 2 TB for 4 hours charges
  8 TB-hours).  Each 45GB comes with 1 CPU.  The minimum memory you can request
  for a given job is 128GB.

You can use your time either by submitting scripts to the `slurm` queue (see
below) or by doing an interactive session.  In either case, your usage is charged
based on how long you're using the requested resources: if you request to use
Bridges Large for 4 hours with 500GB of memory, but your code only uses 250GB,
you'll still get charged 2 TB-hour.  However, if your code finishes after only 2
hours (regardless of memory usage), you'd get charged 1 TB-hour.  So you'll
need to be strategic in requesting enough memory to accomplish your task, but
not so much that you waste your allocation.  Overestimating the time isn't a
problem, as long as it's not so large that you get stuck waiting in the queue
(e.g., >30 hours).


*************
Using `slurm`
*************

If you're running on XSEDE or another system that uses the slurm queue, you may
wish to use `write_sbatch_file.py`.  This will create a job file that can be
submitted with ``sbatch``. More information about how this file is constructed
can be found in the TACC user guide
`here <https://portal.tacc.utexas.edu/archives/stampede#slurm-job-control>`_.

Here is an example call to `write_sbatch_file.py` that shows some of its
functionality.  Note that for XSEDE, `many different programs
<https://portal.xsede.org/software>`_ are already installed, so you just need to
load the relevant module (e.g., anaconda).

 .. code-block:: console

    $ # create submission script
    $ python -m beast.tools.write_sbatch_file \
      'sbatch_file.script' './path/to/job/beast_batch_fit_X.joblist' \
      '/path/to/files/projectname/' \
      --modules 'module load anaconda3' 'source activate beast_v1.4' \
      --queue LM --run_time 2:30:00 --mem 250GB


This creates a file ``sbatch_file.script`` with these contents:

 .. code-block:: console

    #!/bin/bash

    #SBATCH -J beast      # Job name
    #SBATCH -p LM            # Queue name
    #SBATCH -t 2:30:00      # Run time (hh:mm:ss)
    #SBATCH --mem 250GB      # Requested memory

    # move to appropriate directory
    cd /path/to/files/projectname/

    # Load any necessary modules
    # Loading modules in the script ensures a consistent environment.
    module load anaconda3
    source activate beast_v1.4

    # Launch a job
    ./path/to/job/beast_batch_fit_X.joblist


Then the file can be submitted:

 .. code-block:: console

    $ sbatch sbatch_file.script


*******************
BEAST XSEDE wrapper
*******************

This section will go through the `METAL XSEDE example
<https://github.com/BEAST-Fitting/beast-examples/tree/master/metal_xsede>`_.
The wrapper `run_beast_xsede.py` follows the
:ref:`production run workflow<beast_standard_workflow>`,
but at relevant steps, writes out `sbatch` files that the user can then submit
to the `slurm` queue.  The example has addition supplementary files that will
be described here, too.

The XSEDE workflow generally goes as follows:

 1. Type ``sbatch submit_beast_wrapper.script`` to submit the wrapper.
 2. This will run the wrapper.  Once it reaches a step that writes `sbatch`
    file(s), it will stop and write out a text file with the commands to run.
    The wrapper is set up to loop through fields, so once it gets to that point
    for one field, it'll continue on to the next field, until it's looped through
    all fields and written all necessary `sbatch` files.
 3. Submit the `sbatch` commands.
 4. Once those have finished running, do ``sbatch submit_beast_wrapper.script``
    to submit the wrapper again.  It'll see that new files exist, and progress
    along the workflow until it reaches the next set of sbatch files.
 5. Repeat steps 3 and 4 until everything is done!

For the wrapper `run_beast_xsede.py` itself, here is what happens when it runs:

 1. Make source density and background maps.  Determine which one has the most
    dynamic range, and choose that one to split observations.
 2. Write out a `beast_settings` file for the field.
 3. Make SED grid
    * If all subgrids don't exist: Write an `sbatch` script to make any missing
      SED subgrids, then go to the next field.  For METAL, different fields have
      different combinations of filters, so this step is really copying out the
      necessary columns from the master grid file (details below).
    * If all subgrids do exist: Continue onto step 4.
 4. Make quality cuts to photometry and fake stars
 5. 
