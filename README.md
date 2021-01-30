# UoPfluid

A Computational Fluid Dynamic multi-solver code, that is targeted in the solution of fluid dynamics, heat transfer and magnetohydrodynamics equations of motion.

The code is written in object-oriented Fortran, while the flow is visualized with Python scripts, using the matplotlib tool. CPU parallilazation is enabled with the OpenMP protocol.


## Solvers and Capabilities

### Fluid Dynamics
#### Fluid Dynamics
##### Fluid Dynamics
###### Fluid Dynamics




## Troubleshooting

It has been reported that sometimes in Ubuntu 20.04 LTS the latest version of `gfortran` compiler fails to compile the solvers. To get around this issue install an older version of `gfortran`, such as version 5.h The current version can be found with the `gfortran --version` command.

To do that first delete all installed versions of `gfortran` and install `gfortran-5`, with the following commands.

```bash
$ sudo apt-get remove --auto-remove gfortran
$ sudo apt-get install gfortran-5
```

Make sure that the `Makefile` calls the correct compiler:

```Makefile
CC = gfortran-5
```

## License

The `UoPfluid` code is issued under the [MIT](https://choosealicense.com/licenses/mit/) license. 