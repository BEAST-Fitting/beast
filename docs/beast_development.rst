ABEAST Development
=================

You are encouraged to help maintain and improve the BEAST.
It is recommened that you work on a 'branch' of your own
'fork' of the official BEAST distro, and then send 'pull
requests' to the project to contribute your changes.


Fork the BEAST distro
=====================

- Log in to your gitHub account, and on this page's
  < https://github.com/karllark/beast >
  top right corner click on 'fork' to create your own
  'master' copy of the official BEAST distro.

- Clone a copy of your fork.

- CAUTION: do not do it at the same location where you
    cloned the official BEAST distro, or at least move the
    official distro temporarily.

- Example of cloning the 'fork' into beast-YourName at location
    where the official distro lives in 'beast'
    - $ mv beast beast-official
    - $ git clone https://github.com/YourName/beast.git
    - $ mv beast beast-YourName
    - $ mv beast-official beast


Branching Off A Fork
====================

- Change directory into you fork's clone
  - $ cd beast-YourName

- Create a branch (e.g., named beast-dev1)
  - $ git checkout -b beast-dev1

- Now 'git branch' anywhere below 'beast-YourName' will show  all
  branches of the fork, with '*' indicating which branch you are
  currently working on. Here, 'master' refers to the fork's master.

- To 'upload'  this branch on your GitHub fork's repo -
  - $ git push origin beast-dev1

- To revert back to your fork's master -
  - $ git checkout master

    
Making Changes
==============

It is recommended that you create a new branch for every file
you edit or every new file you create. You may very well work
on multiple files within the same branch, but then you no longer
have the option of 'contributing' changes to just one file.

However, if you are make interdependant changes to multiple file,
then they should be part of the same branch and contributed back
to the project at the same time. 

- Anywhere below 'beast-YourName', switch to the branch you wish
  to work off of
  - $ git checkout beast-dev1

- Make changes to the existing files as you wish and/or create
  new files.

- To put any new file (e.g., newfile.py) under git version control
  - $ git add newfile.py

- To 'commit' all changes
  - $ git commit
    
- To see what changes have been made since last commit
  - $ git status

- To see the status of or commit changes of a single file
  - $ git status PathToFile/filename
  - $ git commit PathToFile/filename

- To undo all changes made to a file since last commit
  - $ git checkout PathToFile/filename

- To sync changes made to the branch locally with your GitHub repo
  - $ git push origin beast-dev1


Collaborating and Contributing
==============================

Once you have changes that you'd like to contribute back to the
project or share with collaborators, it recommneded that you use
gitHub 'pull requests'. It is a good idea to check with the project
or your collaborators which branch of their BEAST repo you should
send the pull requests. 

Note: Generally in git-lingo, 'Pull' is to 'download' what 'Push' is
to 'upload'. When you are making a 'pull request', you are requesting
that your contributions are 'pulled' from the other side. So you are not
pushing it, but the other party is pulling it :-)

- Use 'git add', 'git commit' and 'git push' as summarized earlier to
  sync your local edits with your gitHub repo

- From the gitHub page of your fork of BEAST, e.g.,
  < https://github.com/rubab1/beast/branches >
  click on 'Branches'. Next to the name of the branch on which you
  commited/pushed the changes, click on 'New pull request'. Verify that
  names of the target repo ('base fork') and branch ('master') *to* which
  you want to send the pull request, and those of your repo ('head fork')
  and your branch ('compare') *from* which you are sending the pull request
  match what you intend to do.

- In the comments section briefly describe the changes/additions you made
  and submit the pull request

- It is at the other party's (project, collaborator etc.) discretion if to
  accept the changes and merge them with their repo.

    
Staying up-to-date
==================

The BEAST project's official repository will be updated from time to time
to accomodate bug fixes, improvements and new features. You may keep your
fork's master repo up to date with the following steps.

It is highly recommeneded that you do this if you intend to contribute
changes back to the project. Creating new branches off of an up-to-date
fork-master minimizes the chances of conflicting contributions, duplicative
efforts and other complications.

- Switch to your fork's master branch
  - $ git checkout master

- Fetch the project's up-to-date distribution
  - $ git fetch upstream

- Merge the project-master (upstream) with your fork's master (master)
  - $ git merge upstream/master

- Any branch created off of the fork's master now will start from the
  currect BEAST distro and *not* contain any changes made to any prior
  branch, unless those changes have been incorporated into the official
  distro via an accepted pull request and merge


Managing Conflicts / Re-basing
==============================

Let's consider a situation where a fork's master has been updated. A local
branch (e.g., beast-dev1) was created before the update and it has changes
that hadn't been contributed back to the project. As a results, there may
be conflicting versions of some files. The following steps can resolve this.

- Make sure that you are on the correct branch
  - $ git checkout beast-dev1

- *DO NOT SKIP THIS* Make a backup.
  - $ git branch beast-dev1-backup beast-dev1

- Fetch the project's up-to-date distribution
  - $ git fetch upstream
    
- 'Re-base' the branch
  - $ git rebase upstream/master
  - This step may continue to fail until you resolve all conflicts
  - If something goes wrong during re-base, you can start over
    - $ git rebase --abort
  - If the re-base goes fine but later you wish to restore the backup
    - git reset --hard beast-dev1-backup
    
- Once all conflicts have been resolved and the re-base goes through,
  you can delete the backup branch (not recommended)
  - $ git branch -D beast-dev1-backup

- Instead of re-basing a branch, you can do this instead. This is less
  elegant but simpler / easier for beginners
  - Backup your current branch
  - Update and push your fork's master (see 'Staying up to date')
  - Create a new branch from updated fork-master
  - Resolve conflicts between the two branches
  - Commit and push the newer branch
 
