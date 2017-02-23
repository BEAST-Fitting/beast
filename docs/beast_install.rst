Obtaining the BEAST
===================

- There are two alternative ways to obtain the BEAST
  
  - The easiest way is to click on the green 'Clone or download'
    button on this page < https://github.com/karllark/beast > and
    then select 'Download Zip'.  Now simply unzip that file.

  - It may be better to checkout the BEAST repo from GitHub,
    as it will make it easier to stay up to date with code
    improvements.
    - $ git clone https://github.com/karllark/beast.git

  - In either case, you will ahve a directory named 'beast' and a
    sub-directory 'beast/beast'
	
- Now, obtain the BEAST 'libs' file and copy it into 'beast/beast',
  or place a soft link pointing to it, e.g.,
  - $ cd beast/beast
  - $ wget -r <www.TBD.edu/libs>
  - OR
  - ln -s /usr/local/bin/beast-libs libs
