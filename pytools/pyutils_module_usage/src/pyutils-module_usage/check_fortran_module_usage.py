#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
import re
from collections import defaultdict

# Join continuation lines ending with '&'
def join_continued_lines(lines):
    joined_lines = []
    buffer = ''
    for line in lines:
        stripped = line.rstrip()
        if stripped.endswith('&'):
            buffer += stripped[:-1] + ' '
        else:
            buffer += stripped
            joined_lines.append(buffer)
            buffer = ''
    if buffer:
        joined_lines.append(buffer)
    return joined_lines

# Find unused subroutines written in the module
def find_unused_subroutines(fortran_file):
    with open(fortran_file, 'r') as f:
        lines = f.readlines()
    lines = join_continued_lines(lines)

    # Find all subroutine definitions
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
                print("line", line)
                usage[subroutine] = True
                break
            if pattern.search(line):
                usage[subroutine] = True
                break

    unused = [s for s in subroutines if not usage[s]]

    print("Routines written in '{}' NOT USED:".format(fortran_file))
    for subroutine in unused:
        print(subroutine)

    return unused, subroutines, lines

# Find unused routines from other modules
def find_unused_routines_from_other_modules(fortran_file):
    with open(fortran_file, 'r') as f:
        lines = f.readlines()
    lines = join_continued_lines(lines)

    # Find all routines listed in 'use ... only : ...'
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

    print("usage[]: ", usage)
    unused = [r for r in routines if not usage[r]]

    print("'use ... only :' routines from '{}' NOT USED:".format(fortran_file))
    for routine in unused:
        print(routine)

    return unused, use_lines, lines

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Check for unused subroutines written in Fortran file.")
    parser.add_argument('fortran_file', help="Path to the Fortran file to check")
    args = parser.parse_args()

    find_unused_subroutines(args.fortran_file)
    
    parser = argparse.ArgumentParser(description="Check for unused 'use ... only :' routines in Fortran file.")
    parser.add_argument('fortran_file', help="Path to the Fortran file to check")
    args = parser.parse_args()

    find_unused_routines_from_other_modules(args.fortran_file)
    
if __name__ == "__main__":
    main()