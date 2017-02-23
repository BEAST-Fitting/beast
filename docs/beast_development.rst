BEAST Development
=================

You are encouraged to help maintain and improve the BEAST. Before doing so,
please familiarize yourself with basic version control and Git workflow
concepts using one or more of these guides:

- < https://guides.github.com/introduction/flow/ >
- < https://lifehacker.com/5983680/how-the-heck-do-i-use-github >
- < https://homes.cs.washington.edu/~mernst/advice/version-control.html >
- < https://www.youtube.com/watch?v=y_YKHXuJ-ak >


Fork the BEAST distro
=====================

- The main BEAST repository lives at < https://github.com/karllark/beast >.
  The master branch of this repository is the version that is distributed.

- Log in to your github account, and on the top right corner of the BEAST
  repository page click on the ``Fork`` button. This will create a copy of the
  repository in your github accout.

- Clone a copy of your fork to your local computer. If you have a copy of
  the official BEAST distro, you may need to rename it; cloning will
  automatically name the folder ``beast``.

- Example of cloning your fork into ``beast-YourName`` while keeping the
  official distribution in ``beast``:

  .. code:: shell
  $ mv beast beast-official
  
  $ git clone https://github.com/YourName/beast.git

  $ mv beast beast-YourName

  $ mv beast-official beast

- Set the value of the fork's ``upstream`` to the official distribution so you
  can incorporate changes made by others to your development fork. In the clone
  of your fork, run the following:

  .. code:: shell
  $ git remote set-url --add upstream https://github.com/karllark/beast.git
 
   
Adding Branches
===============

- Make sure you are in the directory for your fork of the beast. You will be on
  branch ``master`` by default.

- Create and switch to a branch (here named ``beast-dev1``; generally it's good
  practice to give branches names related to their purpose)

  .. code:: shell
  $ git checkout -b beast-dev1
	  
- Instead, if you want to create first a branch and then switch to it:

  .. code:: shell
  $ git branch beast-dev1

  $ git checkout beast-dev1

- To see a list of all branches of the fork, with ``*`` indicating which branch you are
  currently working on:

  .. code:: shell
  $ git branch

- To ``upload`` this branch to your fork:

  .. code:: shell
  $ git push origin beast-dev1

- To revert back to your fork's master branch:

  .. code:: shell
  $ git checkout master

    
Making Changes
==============

It is recommended that branches have a single purpose; for example, if you are working
on adding a test suite and on improving the fitting algorithm, those should be in
branches (e.g.) ``add-test-suite`` or ``improve-fitting-algorithm`` or ``beast-dev1``

- Anywhere below ``beast-YourName``, switch to the branch you wish to work off of:

  .. code:: shell
  $ git checkout beast-dev1

- Make changes to the existing files as you wish and/or create new files.

- To see what changes have been made at any time:

  .. code:: shell
  $ git status

- To stage any new or edited file (e.g., ``newfile.py``) in preparation for committing:

  .. code:: shell
  $ git add newfile.py

- To add all edited files (*not recommended* unless you are sure of all your changes):

  .. code:: shell
  $ git add -A

- To ``commit`` all changes after adding desired files:

  .. code:: shell
  $ git commit -m ``brief comments describing changes``

- Commit messages should be short but descriptive.
    
- To see the status of or commit changes of a single file:

  .. code:: shell
  $ git status PathToFile/filename

  $ git commit PathToFile/filename
	  
- To undo all changes made to a file since last commit:

  .. code:: shell
  $ git checkout PathToFile/filename

- To sync changes made to the branch locally with your GitHub repo:

  .. code:: shell
  $ git push origin beast-dev1


Collaborating and Contributing
==============================

Once you have changes that you'd like to contribute back to the project or share
with collaborators, you can open a pull request. It is a good idea to check with
the projects or your collaborators which branch of their BEAST repo you should
send the pull requests. 

Note: Generally in git-lingo, ``Pull`` is to ``download`` what ``Push`` is
to ``upload``. When you are making a ``pull request``, you are requesting
that your contributions are ``pulled`` from the other side. So you are not
pushing it, but the other party is pulling it :-)

- Use ``git add``, ``git commit`` and ``git push`` as summarized earlier to
  sync your local edits with your github repo

- From the github page of your fork of BEAST, e.g.,
  <https://github.com/rubab1/beast/branches>
  click on ``Branches``. Next to the name of the branch on which you
  commited/pushed the changes, click on ``New pull request``. Verify that
  names of the target repo (``base fork``) and branch (``master``) *to* which
  you want to send the pull request, and those of your repo (``head fork``)
  and your branch (``compare``) *from* which you are sending the pull request
  match what you intend to do.

- In the comments section briefly describe the changes/additions you made
  and submit the pull request.

- It is at the other party's (project, collaborator etc.) discretion to
  accept the changes and merge them with their repo.

    
Staying up-to-date
==================

The BEAST project's official repository will be updated from time to time
to accommodate bug fixes, improvements and new features. You may keep your
fork's master repo up to date with the following steps.

It is highly recommended that you do this if you intend to contribute
changes back to the project. Creating new branches off of an up-to-date
fork-master minimizes the chances of conflicting contributions, duplicative
efforts and other complications.

- Switch to your fork's master branch:

  .. code:: shell
  $ git checkout master

- Fetch the project's up-to-date distribution:

  .. code:: shell
  $ git fetch upstream

- Merge the project-master (upstream) with your fork's master (master):

  .. code:: shell
  $ git merge upstream/master

- Sync this change with your GitHub repo:

  .. code:: shell
  $ git push origin master


- Any branch created off of the fork's master now will start from the
  correct BEAST distro and *not* contain any changes made to any prior
  branch, unless those changes have been incorporated into the official
  distro via an accepted pull request and merge


Managing Conflicts via Re-basing
================================

Let's consider a situation where a fork's master has been updated. A local
branch (e.g., beast-dev1) was created before the update and it has changes
that hadn't been contributed back to the project. As a results, there may
be conflicting versions of some files. The following steps can resolve this.


- Follow the instructions under ``staying up-to-date`` to update your fork's
  master. *Do not* skip the ``push``.

- Switch to the branch you wish to re-base:

  .. code:: shell
  $ git checkout beast-dev1

- *DO NOT SKIP THIS* Make a backup and push it to your gitHub repo:

  .. code:: shell
  $ git branch beast-dev1-backup beast-dev1

  $ git push origin beast-dev1-backup

- Fetch the project's up-to-date distribution:

  .. code:: shell
  $ git fetch upstream
    
- ``Re-base`` the branch:

  .. code:: shell
  $ git rebase upstream/master

  - This step may continue to fail until you resolve all conflicts

  - Once all conflicts have been resolved and the re-base goes through
    without any error message, push the changes to your gitHub repo:

  .. code:: shell
  $ git push origin beast-dev1
    
  - If something goes wrong during re-base, you can start over:

    .. code:: shell
    $ git rebase --abort

  - If the re-base goes fine but later you wish to restore the backup:

    .. code:: shell
    $ git reset --hard beast-dev1-backup
    
- Once all conflicts have been resolved and the re-base goes through,
  you can delete the backup branch:

  .. code:: shell
  $ git branch -D beast-dev1-backup


Managing Conflicts without Re-basing
====================================

If re-basing a branch on an upstream master keeps failing, an alternative  
is that instead of re-basing a branch, you can resolve the conflicts
manually. This is less elegant but simpler / easier for beginners.
Here are the general steps to follow.

- Merge your fork's master with upstream/master, and push the master

- Create a new branch from updated fork-master, and push the new branch
  
- Switch to and backup the older branch with conflicts, push the backup
  
- Check the differences between the two branches and merge the two branches,
  or edit files on the newer branch to resolve differences
  
- Commit and push the newer branch
  
- Example:

  - Do the preparatory steps

    .. code:: shell
    $ git checkout master

    $ git fetch upstream

    $ git merge upstream/master

    $ git push origin master

    $ git checkout -b beast-dev2

    $ git push origin beast-dev2

    $ git branch beast-dev1-backup beast-dev1

    $ git push origin beast-dev1-backup

    $ git diff beast-dev1 beast-dev2
     
  - Now you can either try to merge the branches:

    .. code:: shell
    $ git checkout beast-dev2

    $ git merge beast-dev1

  - Or manually edit files under beast-dev2 to resolve differences

  - Finally, push the uodated new branch into your gitHub repo:
    (Note: an error free push confirms that all conflicts have been
    resolved both locally and on the gitHub repo)

    .. code:: shell
    $ git push origin beast-dev2
