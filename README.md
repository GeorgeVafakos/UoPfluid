# UoPfluid



## Known issues

Sometimes the latest version of `gfortran` compiler (the default in Ubuntu 20.04 LTS is version 9.3) fails to compile the solvers. To get around this issue the gfortran version must be downgraded to version 5. 

To do that first delete all installed versions of `gfortran` and install `gfortran-5`, with the following commands.

```bash
$ sudo apt-get remove --auto-remove gfortran
$ sudo apt-get install gfortran-5
```

Make sure in the `Makefile` to call the correct compiler:

```Makefile
CC = gfortran-5
```

## License
The `UoPfluid` code is issued under the [MIT](https://choosealicense.com/licenses/mit/) license. 