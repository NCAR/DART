import re

def extract_public_routines(fortran_file):
    """
    Extracts all routines made public in a given Fortran module.
    """
    with open(fortran_file, 'r') as f:
        lines = f.readlines()

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

    lines = join_continued_lines(lines)

    # Regex to find the public statement
    public_pattern = re.compile(r'^\s*public\s*::\s*(.*)', re.IGNORECASE)
    public_routines = []

    for line in lines:
        match = public_pattern.match(line)
        if match:
            print('match: ', match)
            print('line: ', line)
            if '!' in line.strip():
                line = line + line
            # Split the routines listed after 'public ::' and clean them
            items = []
            for item in match.group(1).split(','):
                print("item: ", item)
                print("item.strip(): ", item.strip())
                if '!' in item.strip():
                    print("! FOUND !!!!!!!!!!!!!!!!!!!: ", item.strip())
                    continue
                items.append(item.strip())
            public_routines.extend(items)

    print("Public routines in '{}':".format(fortran_file))
    for routine in public_routines:
        print(routine)

    return public_routines

# Example usage
if __name__ == "__main__":
    fortran_file = '/Users/masmith/Desktop/pytools/DART/assimilation_code/modules/utilities/utilities_mod.f90'
    public_routines = extract_public_routines(fortran_file)
