# pyqceff

Python package for reading and displaying QCEFF table information.

## Installation

To install `pyqceff` in editable/development mode using pip:

   ```sh
   pip install -e DART/pytools/pyqceff
   ```

where DART is the path to DART.

> **Note:** If you update DART with `git pull`, using the editable install (`pip install -e DART/pytools/pyqceff`) ensures your installed package always matches the latest source code.


## Usage

After installation, you can use the `display_qceff_table.py` script to read and display QCEFF tables:

```sh
display_qceff_table <path_to_csv_file> [--details] [QTY]
```

To view the help

```sh
display-qceff-table --help
```


