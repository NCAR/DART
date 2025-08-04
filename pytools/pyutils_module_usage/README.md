# pyutils_module_usage

Python package for determining subroutine usage in a Fortran module.

## Installation

To install `pyutils_module_usage` in editable/development mode using pip:

   ```sh
   pip install -e DART/pytools/pyutils_module_usage
   ```

where DART is the path to DART.

> **Note:** If you update DART with `git pull`, using the editable install (`pip install -e DART/pytools/pyutils_module_usage`) ensures your installed package always matches the latest source code.


## Usage

After installation, you can use the `check_fortran_module_usage.py` script to display unused subroutines in the specified module:

```sh
check_fortran_module_usage <path_to_fortran_file>
```

To view the help

```sh
check_fortran_module_usage --help
```
