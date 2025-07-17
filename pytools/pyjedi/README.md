# pyjedi

Python package for JEDI utilities.

## Installation

To install `pyjedi` in editable/development mode using pip:

   ```sh
   pip install -e DART/pytools/pyjedi
   ```

where DART is the path to DART.

> **Note:** If you update DART with `git pull`, using the editable install (`pip install -e DART/pytools/pyjedi`) ensures your installed package always matches the latest source code.


## Usage

After installation, you can use the `ioda2df.py` function to read a IODA observation file into a pandas dataframe:

```sh
python3
>> import pyjedi as pj
df = pj.ioda2df('path_to_ioda_file')
```


