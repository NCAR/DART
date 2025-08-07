# pyfortran

Python package for Fortran-based utilities, including a program to determine subroutine usage in a Fortran module.

## Installation

Standard installation:

   ```sh
   pip install pyfortran
   ```

To install `pyfortran` in editable/development mode using pip:

   ```sh
   pip install -e DART/pytools/pyfortran
   ```

where DART is the path to DART.

> **Note:** If you update DART with `git pull`, using the editable install (`pip install -e DART/pytools/pyfortran`) ensures your installed package always matches the latest source code.


## Usage

After installation, you can use the `check_fortran_module_usage.py` script to display unused subroutines in the specified module:

```sh
check_fortran_module_usage <path_to_fortran_file>
```

To view the help

```sh
check_fortran_module_usage --help
```