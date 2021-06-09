# UoPfluid

A Computational Fluid Dynamic multi-solver code, that is targeted in the solution of fluid dynamics, heat transfer and magnetohydrodynamics equations of motion. The code is written in object-oriented Fortran, while the flow is visualized with Python scripts, using the matplotlib tool. CPU parallelization is enabled with the OpenMP protocol.

skataaaaaa

## Solvers and Capabilities


## Installation

To install the source code in your local system, browse to the desired directory and run the command

```bash
$ git clone https://github.com/GeorgeVafakos/UoPfluid.git
```
The `UoPfluid` code needs the `GNU gfortran` compiler and the `matplotlib` python library. If the `matplotlib` library isn't already isntalled in the system you can install it by the commands

```bash
$ sudo apt-get install python3-matplotlib
```

## Troubleshooting

Sometimes the latest version of `gfortran` compiler fails to compile the solvers. To get around this issue install version 5 of `gfortran`, as this problem started from version 6. The current version can be found with the `gfortran --version` command.

For Ubuntu 20.04 LTS and later versions, the following lines must be added to the `/etc/apt/sources.list` file:

```bash
deb http://archive.ubuntu.com/ubuntu bionic main universe multiverse restricted
deb http://security.ubuntu.com/ubuntu/ bionic-security main multiverse universe restricted
deb http://archive.ubuntu.com/ubuntu bionic-updates main multiverse universe restricted
```

And update the system

```bash
$ sudo apt-get update
```

To install version 5 `gfortran` of do that first delete all currently installed versions and install `gfortran-5`, with the following commands.

```bash
$ sudo apt-get remove --auto-remove gfortran
$ sudo apt-get install gfortran-5
```

After installing `gfortran-5` make it the default version with the commands:

```bash
$ sudo update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-5 60
$ sudo update-alternatives --config gfortran
```
And choose version 5 as the default version. 

Alternatively, if you don't want version 5 to be the default vertsion of `gfortran` in the system, change the `Makefile` first line to call the correct compiler:

```Makefile
CC = gfortran-5
```

## Published work using UoPfluid

1. G.P. Vafakos & P.K. Papadopoulos (2020). ‘A grid stretching technique for efficient capturing of the MHD boundary layers’, *Fusion Engineering and Design*, vol. 154, https://doi.org/10.1016/j.fusengdes.2020.111477

## License

The `UoPfluid` code is issued under the [MIT](https://choosealicense.com/licenses/mit/) license. 
