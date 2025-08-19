#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""
Script to read and pretty print QCEFF table information from a CSV file.
"""

import csv
import sys
from typing import Dict, List, Any
import argparse


def read_qceff_table(filename: str) -> Dict[str, Any]:
    """
    Read the QCEFF table CSV file and parse its contents.
    
    Args:
        filename: Path to the CSV file

    Returns:
        Dictionary containing parsed table data
    """
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        header1 = next(reader)  # Read the first line (version info)
        headers = next(reader)  # Read the second line (column headers)
        rows = list(reader)

    # Header
    version_info = header1[0]  # "QCEFF table version: X"

    return {
        'version': version_info,
        'headers': headers,
        'data': rows
    }

   
def display_quantity_info(qty_name: str, row: List[str]) -> None:
    """
    Display information for a single quantity in a nicely formatted way.
    
    Args:
        qty_name: Name of the quantity
        row: Row data from CSV
    """
    print(f"\n{'='*80}")
    print(f"QUANTITY: {qty_name}")
    print(f"{'='*80}")
    
    # obs_error_info
    print(f"\nOBS ERROR INFO:")
    print(f"  Bounded Below: {row[1]}")
    print(f"  Bounded Above: {row[2]}")
    print(f"  Lower Bound:   {row[3]}")
    print(f"  Upper Bound:   {row[4]}")

    # probit_inflation
    print(f"\nPROBIT INFLATION:")
    print(f"  Dist Type:     {row[5]}")
    print(f"  Bounded Below: {row[6]}")
    print(f"  Bounded Above: {row[7]}")
    print(f"  Lower Bound:   {row[8]}")
    print(f"  Upper Bound:   {row[9]}")

    # probit_state
    print(f"\nPROBIT STATE:")
    print(f"  Dist Type:     {row[10]}")
    print(f"  Bounded Below: {row[11]}")
    print(f"  Bounded Above: {row[12]}")
    print(f"  Lower Bound:   {row[13]}")
    print(f"  Upper Bound:   {row[14]}")

    # probit_extended_state
    print(f"\nPROBIT EXTENDED STATE:")
    print(f"  Dist Type:     {row[15]}")
    print(f"  Bounded Below: {row[16]}")
    print(f"  Bounded Above: {row[17]}")
    print(f"  Lower Bound:   {row[18]}")
    print(f"  Upper Bound:   {row[19]}")

    # obs_inc_info
    print(f"\nOBS INC INFO:")
    print(f"  Filter Kind:   {row[20]}")
    print(f"  Bounded Below: {row[21]}")
    print(f"  Bounded Above: {row[22]}")
    print(f"  Lower Bound:   {row[23]}")
    print(f"  Upper Bound:   {row[24]}")


def format_bounds(below, above, lower, upper):
    """
    Format bounds for display
    Returns (left_bracket, right_bracket, lower_str, upper_str)
    """
    if below.strip().lower() == '.true.':
        left_bracket = '['
        bound_below = True
    else:
        left_bracket = '('
        bound_below = False
    if above.strip().lower() == '.true.':
        right_bracket = ']'
        bound_above = True
    else:
        right_bracket = ')'
        bound_above = False
    
    lower_str = lower.strip()
    upper_str = upper.strip()
    # Handle lower bound
    if bound_below:
        lower_str = lower_str
    else:
        lower_str = '-inf'
    # Handle upper bound
    if bound_above:
        try:
            if float(upper_str) == -888888:
                upper_str = 'inf'
        except ValueError:
            if upper_str == '' or upper_str.lower() == 'none':
                upper_str = 'inf'
    else:
        upper_str = 'inf'
    
    return left_bracket, right_bracket, lower_str, upper_str

def format_dist(dist, bounded_below, bounded_above, lower, upper):
    """
    Format the distribution type and bounds for display.
    """   
    dist_upper = dist.strip().upper()
    dist_upper = dist_upper.replace('_DISTRIBUTION', '')
    # Shorten BOUNDED_NORMAL_RH_DISTRIBUTION
    dist_upper = dist_upper.replace('BOUNDED_NORMAL_RH', 'BNRH')

    # Truncate long names
    dist_short = dist_upper[:18] + '..' if len(dist_upper) > 20 else dist_upper
    # If normal_distribution, do not show bounds
    if dist_upper == 'NORMAL':
        return dist_short
    # GAMMA_DISTRIBUTION: lower bound at 0
    if dist_upper == 'GAMMA':
        return f"{dist_short} [0,inf)"
    # BETA_DISTRIBUTION: [0,1]
    if dist_upper == 'BETA':
        return f"{dist_short} [0,1]"
    # LOG_NORMAL_DISTRIBUTION: lower bound at 0
    if dist_upper == 'LOG_NORMAL':
        return f"{dist_short} [0,inf)"
    # UNIFORM_DISTRIBUTION: no bounds
    if dist_upper == 'UNIFORM':
        return dist_short
    # Format bounds
    left_bracket, right_bracket, lower_str, upper_str = format_bounds(bounded_below, bounded_above, lower, upper)
    return f"{dist_short} {left_bracket}{lower_str},{upper_str}{right_bracket}"

def format_kind(kind, bounded_below, bounded_above, lower, upper):
    """
    Format the filter kind, with bounds if appropriate.
    """
    kind_upper = kind.strip().upper()
    kind_short = kind_upper[:18] + '..' if len(kind_upper) > 20 else kind_upper
    if kind_upper == 'EAKF': # no bounds
        return kind_short
    if kind_upper == 'ENKF': # no bounds
        return kind_short
    if kind_upper == 'UNBOUNDED_RHF': # no bounds
        return kind_short
    if kind_upper == 'GAMMA_FILTER':
        return f"{kind_short} [0,inf)"
    if kind_upper == 'BOUNDED_NORMAL_RHF' or kind_upper == 'KDE_FILTER':
        left_bracket, right_bracket, lower_str, upper_str = format_bounds(bounded_below, bounded_above, lower, upper)
        return f"{kind_short} {left_bracket}{lower_str},{upper_str}{right_bracket}"
    # Default: just the kind name
    return kind_short

def format_obs_error(obs_error_below, obs_error_above, obs_error_low, obs_error_up):
    """
    Format the observation error information for display.
    """
    obs_error_below = obs_error_below.strip()
    obs_error_above = obs_error_above.strip()
    obs_error_low = obs_error_low.strip()
    obs_error_up = obs_error_up.strip()

    left_bracket, right_bracket, lower_str, upper_str = format_bounds(obs_error_below, obs_error_above, obs_error_low, obs_error_up)
    return f"{left_bracket}{lower_str},{upper_str}{right_bracket}"

def display_summary_table(data: List[List[str]]) -> None:

    """
    Display a summary table of all quantities.
    
    Args:
        data: List of data rows from CSV
    """
    print(f"\n{'='*154}")
    
    # Header
    qty_width = 31
    col_width = 24
    short_col_width = 15
    print(f"{'QUANTITY':<{qty_width}} {'OBS_ERROR':<{short_col_width}} {'PROBIT_INFL':<{col_width}} {'PROBIT_STATE':<{col_width}} {'PROBIT_EXT':<{col_width}} {'FILTER_KIND':<{col_width}}")
    print(f"{'-'*qty_width} {'-'*short_col_width} {'-'*col_width} {'-'*col_width} {'-'*col_width} {'-'*col_width}")

    for row in data:
        qty_name = row[0].upper()
        # obs_error_info
        obs_error_below = row[1] if len(row) > 1 else ''
        obs_error_above = row[2] if len(row) > 2 else ''
        obs_error_low = row[3] if len(row) > 3 else ''
        obs_error_up = row[4] if len(row) > 4 else ''
        obs_error_info = format_obs_error(obs_error_below, obs_error_above, obs_error_low, obs_error_up)

        # Probit inflation
        probit_infl_dist = row[5] if len(row) > 5 else 'N/A'
        probit_infl_below = row[6] if len(row) > 6 else ''
        probit_infl_above = row[7] if len(row) > 7 else ''
        probit_infl_low = row[8] if len(row) > 8 else ''
        probit_infl_up = row[9] if len(row) > 9 else ''
        probit_infl = format_dist(probit_infl_dist, probit_infl_below, probit_infl_above, probit_infl_low, probit_infl_up)

        # Probit state
        probit_state_dist = row[10] if len(row) > 10 else 'N/A'
        probit_state_below = row[11] if len(row) > 11 else ''
        probit_state_above = row[12] if len(row) > 12 else ''
        probit_state_low = row[13] if len(row) > 13 else ''
        probit_state_up = row[14] if len(row) > 14 else ''
        probit_state = format_dist(probit_state_dist, probit_state_below, probit_state_above, probit_state_low, probit_state_up)

        # Probit extended state
        probit_ext_dist = row[15] if len(row) > 15 else 'N/A'
        probit_ext_below = row[16] if len(row) > 16 else ''
        probit_ext_above = row[17] if len(row) > 17 else ''
        probit_ext_low = row[18] if len(row) > 18 else ''
        probit_ext_up = row[19] if len(row) > 19 else ''
        probit_ext = format_dist(probit_ext_dist, probit_ext_below, probit_ext_above, probit_ext_low, probit_ext_up)

        filter_kind = row[20] + '..' if len(row[20]) > 20 else row[20]
        obs_inc_below = row[21] if len(row) > 21 else ''
        obs_inc_above = row[22] if len(row) > 22 else ''
        obs_inc_low = row[23] if len(row) > 23 else ''
        obs_inc_up = row[24] if len(row) > 24 else ''
        obs_inc_info = format_kind(filter_kind, obs_inc_below, obs_inc_above, obs_inc_low, obs_inc_up)

        print(f"{qty_name:<{qty_width}} {obs_error_info:<{short_col_width}} {probit_infl:<{col_width}} {probit_state:<{col_width}} {probit_ext:<{col_width}} {obs_inc_info:<{col_width}}")

def get_parser():
    parser = argparse.ArgumentParser(description="Display QCEFF table summary and details.")
    parser.add_argument("csv_file", help="Path to QCEFF table CSV file")
    parser.add_argument(
        "--details",
        nargs="?",
        const=True,
        metavar="QTY",
        help="Print detailed information for all quantities, or for a specific quantity QTY"
    )
    parser.usage = "%(prog)s csv_file [--details [QTY]]"
    return parser


def main():
    """Main function to run the script."""
    
    parser = get_parser()
    # Enforce csv_file is always first argument, but allow --help/-h anywhere
    if any(arg in ('--help', '-h') for arg in sys.argv[1:]):
        parser.parse_args()  # argparse will handle help and exit
        return
    if len(sys.argv) > 1 and not sys.argv[1].startswith('-'):
        args = parser.parse_args()
    else:
        parser.print_usage()
        print("\nError: csv_file must be the first argument.")
        sys.exit(2)

    filename = args.csv_file

    try:
        # Read the CSV file
        table_data = read_qceff_table(filename)

        print(f"File: {filename}")
        print(f"Version: {table_data['version']}")

        display_summary_table(table_data['data'])

        if args.details:
            print(f"\n{'='*80}")
            print("DETAILED INFORMATION")
            print(f"{'='*80}")

            if args.details is True:
                # Print all quantities
                for row in table_data['data']:
                    qty_name = row[0].upper()
                    display_quantity_info(qty_name, row)
            else:
                # Print only the specified quantity
                found = False
                for row in table_data['data']:
                    qty_name = row[0].upper()
                    # Accept QTY_NAME, NAME, or just the name (case-insensitive)
                    arg = args.details.upper()
                    # Remove QTY_ prefix for comparison if present
                    qty_name_noprefix = qty_name
                    if qty_name.startswith('QTY_'):
                        qty_name_noprefix = qty_name[4:]
                    if (
                        qty_name == arg or
                        qty_name_noprefix == arg or
                        ('QTY_' + arg) == qty_name or
                        qty_name_noprefix == arg.replace('QTY_', '')
                    ):
                        display_quantity_info(qty_name, row)
                        found = True
                        break
                if not found:
                    print(f"Quantity '{args.details}' not found.")

            print(f"\n{'='*80}")
            print("END OF REPORT")
            
        print(f"{'='*154}")

    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error processing file: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
