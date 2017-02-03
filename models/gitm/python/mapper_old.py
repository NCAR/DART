# #!/Library/Frameworks/Python.framework/Versions/2.7/bin/python - specific to me, so don't use it
""" Find arg[1].f90 file in the arg[2] directory (or its subdirs), generate list of "use"d .f90 files, and write them as a JSON string to w.json.
Example: python mapper.py main /Users/alex/GITM2
"""
#    DART $Id$
#    Copyright (C) 2013 by
#    Alexey Morozov <alexeymor gmail>
__author__ = """Alexey Morozov <alexeymor gmail>"""

#the goal string for hierarchical edge d3.js plot, to be put into w.json is:
# [{"name":"analytics","size":3,"imports":["vis"]},{"name":"vis","size":3,"imports":[]}]

#import networkx as nx #didn't work because want "name" instead of "nodes"
#g=dict([('sape', 4139), ('guido', 4127), ('jack', 4098)]) #didn't work when printing - need " instead of '

#so I'll construct the graph string by hand - not the most elegant solution, but at least I'll learn how to work with strings and stacks in python

import sys #for getting command line args

fn = sys.argv[1] #where do I start from (find all use statements in this file and then go into those files and so on)?, appends .f90
hd = sys.argv[2] #where can I search for *.f90 files?

import os # modified http://code.activestate.com/recipes/577027/ 

def find_file(filename, path):
    for root, dirs, names in os.walk(path):
        if filename in names:
            return os.path.join(root, filename)
#    raise 'File not found'



def read_file(hd, stack, visited, files, depends, size):
    while len(stack)>0:
        fn = stack.pop() #pop the last file off the stack and analyze it
        visited.append(fn) #put this file on the visited list
        path = find_file(fn + '.f90', hd) #find the file
        loc_dep = [] #local dependencies - dependencies for the current file
        try:
            f = open(path, 'r')
            files.append(fn) #put this file on the files list
            size.append(os.path.getsize(path))
            for line in f:
                if 'use' in line:
                    line_trunc = line[line.find('use')+3 : line.find(',')].strip() #remove "use", anything after comma, any spaces
                    loc_dep.append(line_trunc)
                    if (line_trunc not in visited) and (line_trunc not in stack): 
                        stack.append(line_trunc)
            f.close()
        except TypeError:
            pass
        depends.append(loc_dep)
        print ' stack right now is ', stack
        (files, depends, size)=read_file(hd, stack, visited, files, depends, size) #call r_f recursively until the stack is empty
    return(files, depends, size)



(files, depends, size) = read_file(hd, [fn], [], [], [], [])
#files is a list containing file names used in the particular code
#depends is a list containing lists of depencies for each file in "files"
#size is a list of file sizes for each file in "files"

#merge the lists into a JSON-like string:
# [{"name":"analytics","size":3,"imports":["vis"]},{"name":"vis","size":3,"imports":[]}]
string = '['
for f,d,s in zip(files, depends, size):
    df = ''
    for d1 in d: #flatten the dependencies first into "vis","hi" form
        df = df + '"' + d1 + '",'
    df = df[:-1] #remove the last comma
    string = string + '{"name":"' + f + '","size":' + str(s) + ',"imports":[' + df + ']},'
string = string[:-1] + ']' #remove the last comma and add the last bracket

print string

f = open('w.json', 'w')
f.write(string)
f.close

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
