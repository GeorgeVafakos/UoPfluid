# UoPfluid

## Platform


## Known issues

Sometimes the latest version of `gfortran` compiler (the default in Ubuntu 20.04 LTS is *Version 9.3*) fails to compile the solvers. To get around this issue install an older version of `gfortran`, such as *Version 5*. You can check the current version with the `gfortran --version` command.

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