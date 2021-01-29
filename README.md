# UoPfluid

## Known issues

Sometimes the latest version of gfortran compiler (the default in Ubuntu 20.04 is gfortran-9) can not compile the solvers. To get around this issue the gfortran version must be downgraded to gfortran 5. 

To do that first delete all installed version of gfortran and install gfortran-5, with the following commands.

> $ sudo apt-get remove --auto-remove gfortran
> $ sudo apt-get install gfortran-5

Make sure in the Makefile to call the correct compiler:

CC = gfortran-5
