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

After installation, you can use the `ioda2obsq` script to convert a IODA observation file into a DART obs_seq.out file:

```sh
ioda2obsq yamlConfigFile inputIodaFile outputObsqFile
```

To view the help
```sh
ioda2obsq -h
```

The yamConfigFile is used to tell ioda2obsq which variables, and their type, from the input IODA file to convert, as well as define the vertical coordinate.
Here is a sample from the radiosonde test case:
```yaml
---

ioda to obsq converter:
  observation variables:
    - name: airTemperature
      type: RADIOSONDE_TEMPERATURE
    - name: specificHumidity
      type: RADIOSONDE_SPECIFIC_HUMIDITY
    - name: windEastward
      type: RADIOSONDE_U_WIND_COMPONENT
    - name: windNorthward
      type: RADIOSONDE_V_WIND_COMPONENT

  observation category:
    name: conventional
    vertical coordinate:
      name: MetaData/pressure
      units: "pressure (Pa)"
```

Note "observation variables" is the list of variables you intend to assimilate in the DA job.
