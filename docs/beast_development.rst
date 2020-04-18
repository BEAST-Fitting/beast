.. _beast_development:

#################
BEAST Development
#################

Basic Workflow
==============

You are encouraged to help maintain and improve the ``beast``. Before doing so,
please familiarize yourself with basic version control and Git workflow
concepts. A good description of the standard fork-branch-pull workflow can be
found here:

- https://gist.github.com/Chaser324/ce0505fbed06b947d962

Some other useful guides are:

- https://guides.github.com/introduction/flow/
- https://lifehacker.com/5983680/how-the-heck-do-i-use-github
- https://homes.cs.washington.edu/~mernst/advice/version-control.html
- https://www.youtube.com/watch?v=y_YKHXuJ-ak

In short, the fork-branch-pull workflow for contributing to the ``beast``
project is:

- Create your own 'fork' of the official ``beast`` repository

- Create purpose-specific 'branches' off your 'fork'

- Make changes or additions within the branches

- Contribute your modified codes to the ``beast`` project or share them with
  your collaborators via 'pull requests'

- Keep your fork updated to benefit from continued development of the
  official version and to minimize version conflicts

- Resolve version conflicts as much as possible before sending pull requests


Development Install
===================

If you plan on modifying the ``beast`` in addition to running the code, it may
be useful to create a development installation. First, create a fork of the
official ``beast`` repository and clone it:

.. code-block:: console

   $ git clone https://github.com/YourName/beast.git

Optionally, you can rename this cloned copy:

.. code-block:: console

   $ git clone https://github.com/YourName/beast.git beast-YourName

Set the value of the fork's 'upstream' to the official distribution so you
can incorporate changes made by others to your development fork. In the clone
of your fork, run the following:

.. code-block:: console

   $ git remote add upstream https://github.com/BEAST-Fitting/beast.git

In order to run a development installation, navigate to the directory in your
``beast`` repository that contains `setup.py`, and run:

.. code-block:: console

   $ pip install -e .

Alternatively, you can perform a development install directly through Python
with:

.. code-block:: console

   $ python setup.py develop


Adding Branches
===============

- Make sure you are in the directory for your fork of the beast. You will be on
  branch `master` by default.

- Create and switch to a branch (here named `beast-dev1`):

  .. code-block:: console

     $ git checkout -b beast-dev1

- To see a list of all branches of the fork, with '*' indicating which branch
  you are currently working on:

  .. code-block:: console

     $ git branch

- To 'upload' this branch to your fork:

  .. code-block:: console

     $ git push origin beast-dev1

- To revert back to your fork's `master` branch:

  .. code-block:: console

     $ git checkout master


Making Changes
==============

It is recommended that branches have a single purpose; for example, if you are working
on adding a test suite, improving the fitting algorithm, and speeding up some task,
those should be in separate branches (e.g. `add-test-suite`, `improve-fitting-algorithm`
and `beast-dev1`).

- Switch to the branch you wish to work off of:

  .. code-block:: console

     $ git checkout beast-dev1

- Make changes to the existing files as you wish and/or create new files.

- To see what changes have been made at any time:

  .. code-block:: console

     $ git status

- To stage any new or edited file (e.g., 'newfile.py') in preparation for committing:

  .. code-block:: console

     $ git add newfile.py

- To add all edited files (*not recommended* unless you are sure of all your changes):

  .. code-block:: console

     $ git add -A

- To 'commit' all changes after adding desired files:

  .. code-block:: console

     $ git commit -m 'brief comments describing changes'

- Commit messages should be short but descriptive.

- To see the status of your changed files:

  .. code-block:: console

     $ git status

- To view any differences between a file and the last committed version:

  .. code-block:: console

     $ git diff PathToFile/filename

- To undo all changes made to a specific file since the last commit:

  .. code-block:: console

     $ git checkout PathToFile/filename

- To sync changes made to the branch locally with your GitHub repository:

  .. code-block:: console

     $ git push origin beast-dev1


Test Changes
============

It is a good idea to test that your changes have not caused problems.  In the
base ``beast`` directory the following commands may be run to do this.

Run existing tests, including a regression test against a full ``beast`` model
run.  Once the command below has finished, the coverage of the tests can
be viewed in a web browser by pointing to files in the `htmlconv` subdirectory
(which gets produced when the tests are run).

  .. code-block:: console

     $ python setup.py test --remote-data --coverage

Make sure the documentation can be created.

  .. code-block:: console

     $ python setup.py build_docs

The resulting HTML files are placed in the `docs/_build/html` subdirectory, and
can be viewed in a web browser.


Submitting a Pull Request
=========================

Once you have changes that you'd like to contribute back to the upstream branch,
you can open a pull request for review. Pull requests can be submitted at
https://github.com/BEAST-Fitting/beast/pulls. If you push any commits to your
origin repository in a development branch (`beast-dev1`), then a "Compare &
pull request" button should appear at the top of this site. Briefly describe the
changes/additions you made in the comments section and submit the pull request.


Staying up-to-date
==================

The ``beast`` project's official repository will be updated from time to time
to accommodate bug fixes, improvements and new features. You can keep your
fork's `master` repository up-to-date with the following steps:

- Switch to your fork's `master` branch:

  .. code-block:: console

     $ git checkout master

- Fetch the project's up-to-date distribution:

  .. code-block:: console

     $ git fetch upstream

- Merge the official (upstream) `master` branch with your fork's `master` branch:

  .. code-block:: console

     $ git merge upstream/master

- Sync this change with your origin repository:

  .. code-block:: console

     $ git push origin master


BEAST on Slack
==============

There is a ``beast`` space on Slack.  Email kgordon@stsci.edu for an invite.


Visualizing Repository Commits
==============================

The commits to the ``beast`` repository can be visualized using `gource`.  This
creates a movie showing the time evolution of the code and who made the
changes.

Version created 22 Jan 2018:  <http://stsci.edu/~kgordon/beast/beast_repo.mp4>

Command to create it:

    .. code-block:: console

        $ gource -s .06 -1280x720 --auto-skip-seconds .1 --multi-sampling  --stop-at-end --key --highlight-users --hide mouse,progress --file-idle-time 0 --max-files 0  --background-colour 000000 --font-size 22 --title "This is beast" --output-ppm-stream - --output-framerate 30 | avconv -y -r 30 -f image2pipe -vcodec ppm -i - -b 65536K beast_repo.mp4
