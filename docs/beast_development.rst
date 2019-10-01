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

It is much easier to perform development if any changes you make to the code are
immediately reflected in how the ``beast`` runs (in contrast to needed to perform a
new install each team). This can be achieved by using a development install.

This is most easily achieved via a pip development install. Navigate to the
directory in your your ``beast`` repository that contains `setup.py`, and run:

  .. code-block:: console

     $ mv beast beast-official
     $ git clone https://github.com/YourName/beast.git
     $ mv beast beast-YourName
     $ mv beast-official beast

- Set the value of the fork's 'upstream' to the official distribution so you
  can incorporate changes made by others to your development fork. In the clone
  of your fork, run the following:

  .. code-block:: console

     $ git remote add upstream https://github.com/BEAST-Fitting/beast.git

Adding Branches
===============

- Make sure you are in the directory for your fork of the beast. You will be on
  branch 'master' by default.

- Create and switch to a branch (here named 'beast-dev1'; generally it's good
  practice to give branches names related to their purpose)

  .. code-block:: console

     $ git checkout -b beast-dev1

- Instead, if you want to create first a branch and then switch to it:

  .. code-block:: console

     $ git branch beast-dev1
     $ git checkout beast-dev1

- To see a list of all branches of the fork, with '*' indicating which branch you are
  currently working on:

  .. code-block:: console

     $ git branch

- To 'upload' this branch to your fork:

  .. code-block:: console

     $ git push origin beast-dev1

- To revert back to your fork's master branch:

  .. code-block:: console

     $ git checkout master


Development Install
==============================

To perform development, and see your changes reflected immediately in your
installed copy of the BEAST, you can perform a development install. This can
either be performed via a pip development install, by navigating to the
directory that contains `setup.py` and running:

  .. code-block:: console

     $ pip install -e .

Alternatively, you can perform a development install directly though Python
with:

  .. code-block:: console

     $ python setup.py develop


Making Changes
==============

It is recommended that branches have a single purpose; for example, if you are working
on adding a test suite, on improving the fitting algorithm and on speeding up some task,
those should be in separate branches (e.g.) `add-test-suite`, `improve-fitting-algorithm`
and `beast-dev1`.

- Anywhere below 'beast-YourName', switch to the branch you wish to work off of:

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

- To see the status of or commit changes of a single file:

  .. code-block:: console

     $ git status PathToFile/filename
     $ git commit PathToFile/filename

- To undo all changes made to a file since last commit:

  .. code-block:: console

     $ git checkout PathToFile/filename

- To sync changes made to the branch locally with your GitHub repo:

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

The resulting HTML files are placed in `docs/_build/html` subdirectory, and
can be viewed in a web browser.


Collaborating and Contributing
==============================

BEAST on Slack
==============

There is a ``beast`` space on Slack.  Email kgordon@stsci.edu for an invite.


Visualizing Repository Commits
==============================

The commits to the ``beast`` repository can be visualized using `gource`.  This
creates a movie showing the time evolution of the code and who make the
changes.

Version created 22 Jan 2018:  <http://stsci.edu/~kgordon/beast/beast_repo.mp4>

Command to create it:

    .. code-block:: console

        $ gource -s .06 -1280x720 --auto-skip-seconds .1 --multi-sampling  --stop-at-end --key --highlight-users --hide mouse,progress --file-idle-time 0 --max-files 0  --background-colour 000000 --font-size 22 --title "This is ``beast``" --output-ppm-stream - --output-framerate 30 | avconv -y -r 30 -f image2pipe -vcodec ppm -i - -b 65536K ``beast``_repo.mp4
