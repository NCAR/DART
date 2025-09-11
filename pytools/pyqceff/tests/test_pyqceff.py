import subprocess
import sys
import os

def test_display_qceff_table_all_distributions():
    # Create a temporary CSV file with all possible distributions and bounds

    script_path = os.path.join(os.path.dirname(__file__), "../src/pyqceff/display_qceff_table.py")
    csv_path = os.path.join(os.path.dirname(__file__), "test.csv")

    # Run the script and capture output
    result = subprocess.run(
        [sys.executable, script_path, csv_path, "--details"],
        capture_output=True,
        text=True
    )


    # Check that output contains expected distribution names and bounds
    assert "QTY_ONE" in result.stdout
    assert "QTY_TWO" in result.stdout
    assert "QTY_THREE" in result.stdout
    assert "QTY_FOUR" in result.stdout
    assert "QTY_FIVE" in result.stdout
    assert "QTY_SIX" in result.stdout
    assert "QTY_SEVEN" in result.stdout
    assert "QTY_EIGHT" in result.stdout
    assert "QTY_NINE" in result.stdout

    # Optionally, check for errors
    assert result.returncode == 0