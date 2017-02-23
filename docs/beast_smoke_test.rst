Checking if the BEAST installation worked fine
==============================================

There is a small example script named 'run_beast.py' located at
'beast/beast/examples/phat_small' running which is a quick check
of if the BEAST installation is working. 

- In 'beast/beast/examples/phat_small', place a soft link
  named 'beast' pointing two levels up
  - $ cd beast/beast/examples/phat_small
  - $ ln -s ../../ beast

- Verify that the default python installation is version 2.7
  - $ python --version
  - If you installed through AstroConda, first activate
    activate the correct AstroConda environment, e.g.,
    - $ source activate astroconda2

- Check out the basic gelp message of 'run_beast.py'
  - $ ./run_beast.py -h

- Make BEAST say 'Hello World'!
  - $ ./run_beast.py -potf
