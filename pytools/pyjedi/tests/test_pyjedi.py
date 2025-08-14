import subprocess
import sys
import os
import pyjedi

def test_ioda2obsq_radiosonde():
    # Set up paths to all the necessary files
    testDir = os.path.dirname(__file__)
    command = 'ioda2obsq'
    config_path = os.path.join(testDir, "obsq.radiosonde.yaml")
    ioda_path = os.path.join(testDir, "sondes_obs_2018041500.200.nc4")
    obsq_path = os.path.join(testDir, "obs_seq.radiosonde.out")

    # Run the converter, ioda2obsq
    result = subprocess.run(
        [command, config_path, ioda_path, obsq_path],
        capture_output=True,
        text=True
    )

    # Check results
    print(result.stdout)
    if (result.returncode != 0):
        print("Error messages from ioda2obsq: \n", result.stderr)

    assert result.returncode == 0, "Non-zero return code detected"

    assert "RADIOSONDE_TEMPERATURE" in result.stdout, "Did not find 'RADIOSONDE_TEMPERATURE' in ioda2obsq output"
    assert "RADIOSONDE_SPECIFIC_HUMIDITY" in result.stdout, "Did not find 'RADIOSONDE_SPECIFIC_HUMIDITY' in ioda2obsq output"
    assert "RADIOSONDE_U_WIND_COMPONENT" in result.stdout, "Did not find 'RADIOSONDE_U_WIND_COMPONENT' in ioda2obsq output"
    assert "RADIOSONDE_V_WIND_COMPONENT" in result.stdout, "Did not find 'RADIOSONDE_V_WIND_COMPONENT' in ioda2obsq output"
    assert "Converted 169 observations" in result.stdout, "Did not find 'Converted 169 observations' in ioda2obsq output"

def test_ioda2obsq_amsua():
    # Set up paths to all the necessary files
    testDir = os.path.dirname(__file__)
    command = 'ioda2obsq'
    config_path = os.path.join(testDir, "obsq.amsua.yaml")
    ioda_path = os.path.join(testDir, "amsua_n19_obs_2018041500.10.nc4")
    obsq_path = os.path.join(testDir, "obs_seq.amsua.out")

    # Run the converter, ioda2obsq
    result = subprocess.run(
        [command, config_path, ioda_path, obsq_path],
        capture_output=True,
        text=True
    )

    # Check results
    print(result.stdout)
    if (result.returncode != 0):
        print("Error messages from ioda2obsq: \n", result.stderr)

    assert result.returncode == 0, "Non-zero return code detected"

    assert "NOAA_19_AMSUA_TB" in result.stdout, "Did not find 'NOAA_19_AMSUA_TB' in ioda2obsq output"
    assert "Converted 70 observations" in result.stdout, "Did not find 'Converted 169 observations' in ioda2obsq output"

if __name__ == "__main__":
    for func in [test_ioda2obsq_radiosonde, test_ioda2obsq_amsua]:
        try:
            print(f"****** Testing {func.__name__} ******")
            func()
        except AssertionError as ae:
            print(f"AssertionError in {func.__name__}: {ae}")
        except Exception as e:
            print(f"Exception in {func.__name__}: {e}")

