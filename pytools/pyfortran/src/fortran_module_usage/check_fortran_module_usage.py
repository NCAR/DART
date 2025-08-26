#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
import re
from collections import defaultdict

def get_parser():
    import argparse
    
    parser = argparse.ArgumentParser(description="Check for unused subroutines in a given Fortran module.")
    parser.add_argument('fortran_file', help="Path to the Fortran file to check")
    parser.usage = "%(prog)s fortran_file"
    return parser

# Join continuation lines ending with '&'
def join_continued_lines(lines):
    joined_lines = []
    buffer = ''
    for line in lines:
        stripped = line.rstrip()
        if '&' in stripped:
            position = stripped.find('&')
            buffer += stripped[:position] + ' '
        else:
            buffer += stripped
            joined_lines.append(buffer)
            buffer = ''
    if buffer:
        joined_lines.append(buffer)
    return joined_lines

# Find unused subroutines written in the module
def find_unused_subroutines(lines):
    # Find all subroutine definitions and extract routines
    subroutine_pattern = re.compile(r'^\s*subroutine\s+(\w+)', re.IGNORECASE)
    subroutines = []
    for idx, line in enumerate(lines):
        m = subroutine_pattern.match(line)
        if m:
            subroutines.append(m.group(1))

    # Check if each subroutine is called in the file
    usage = defaultdict(bool)
    for subroutine in subroutines:
        pattern = re.compile(r'\bcall {}\b'.format(re.escape(subroutine)))
        for line in lines:
            # Skip lines that are commented out
            stripped_line = line.strip()
            if stripped_line.startswith('!'):
                continue
            # Assuming subroutines marked as public are used in other modules
            if subroutine in line and 'public' in line:
                usage[subroutine] = True
                break
            if subroutine in line and 'module procedure' in line:
                usage[subroutine] = True
                break
            if pattern.search(line):
                usage[subroutine] = True
                break
    unused = [s for s in subroutines if not usage[s]]
    return unused

# Find unused routines from other modules
def find_unused_routines_from_other_modules(lines):
    # Find all 'use ... only :' statements and extract routines
    use_pattern = re.compile(r'^\s*use\s+\w+.*only\s*:\s*(.*)', re.IGNORECASE)
    routines = []
    use_lines = []
    for idx, line in enumerate(lines):
        m = use_pattern.match(line)
        if m:
            items = [item.strip().replace('&','') for item in m.group(1).split(',')]
            routines.extend([item for item in items if item])
            use_lines.append(idx)

    # Remove duplicates and ignore "operator( )" routines and routines with "=>"
    routines = sorted(set([r for r in routines if r and not r.strip().lower().startswith('operator(') and '=>' not in r]))

    # Check if each routine is used in the file
    usage = defaultdict(bool)
    for routine in routines:
        pattern = re.compile(r'\b{}\b'.format(re.escape(routine)))
        for line in lines:
            # Skip lines that are commented out
            stripped_line = line.strip()
            if stripped_line.startswith('!'):
                continue
            if routine in line and 'use' in line:
                continue
            if pattern.search(line):
                usage[routine] = True
                break
    # Find unused routines
    unused = [r for r in routines if not usage[r]]
    return unused

def main():
    import sys
    
    parser = get_parser()
    # Allow --help/-h anywhere
    if any(arg in ('--help', '-h') for arg in sys.argv[1:]):
        parser.parse_args()  # argparse will handle help and exit
        return
    args = parser.parse_args()

    with open(args.fortran_file, 'r') as f:
        lines = f.readlines()
        lines = join_continued_lines(lines)
    
    unused = find_unused_subroutines(lines)
    print("Routines written in '{}' NOT USED:".format(args.fortran_file))
    for subroutine in unused:
        print(subroutine)

    unused = find_unused_routines_from_other_modules(lines)
    print("'use ... only :' routines from '{}' NOT USED:".format(args.fortran_file))
    for subroutine in unused:
        print(subroutine)
    
if __name__ == "__main__":
    main()