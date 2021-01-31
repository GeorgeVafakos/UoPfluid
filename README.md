# UoPfluid

A Computational Fluid Dynamic multi-solver code, that is targeted in the solution of fluid dynamics, heat transfer and magnetohydrodynamics equations of motion. The code is written in object-oriented Fortran, while the flow is visualized with Python scripts, using the matplotlib tool. CPU parallelization is enabled with the OpenMP protocol.


## Solvers and Capabilities


## Installation

To install the source code in your local system, browse to the desired directory and run the command

```bash
$ git clone https://github.com/GeorgeVafakos/UoPfluid.git
```
The `UoPfluid` code needs the `GNU gfortran` compiler and the `matplotlib` python library. The `matplotlib` is included in the `anaconda` package.


## Troubleshooting

Sometimes in Ubuntu 20.04 LTS the latest version of `gfortran` compiler fails to compile the solvers. To get around this issue install an older version of `gfortran`, such as version 5. The current version can be found with the `gfortran --version` command.

To do that first delete all installed versions of `gfortran` and install `gfortran-5`, with the following commands.

```bash
$ sudo apt-get remove --auto-remove gfortran
$ sudo apt-get install gfortran-5
```

Make sure that the `Makefile` calls the correct compiler:

```Makefile
CC = gfortran-5
```

## Published work using UoPfluid

1. G.P. Vafakos & P.K. Papadopoulos (2020). ‘A grid stretching technique for efficient capturing of the MHD boundary layers’, *Fusion Engineering and Design*, vol. 154, https://doi.org/10.1016/j.fusengdes.2020.111477

## License

The `UoPfluid` code is issued under the [MIT](https://choosealicense.com/licenses/mit/) license. 
