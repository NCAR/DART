# #!/Library/Frameworks/Python.framework/Versions/2.7/bin/python - specific to me, so don't use it
""" Find arg[1].f90 file in the arg[2] directory (or its subdirs), generate list of "use"d .f90 files, and write them as a JSON string to he.json. If the resulting pdf doesn't look like what you expected (too few nodes or a node you expected is missing, search this script's standard output for "ignored" and see if any of the "right" files was replaced by "wrong" ones (typical confusion comes from main.f90 in GITM2/src and main.f90 in GITM2/util/DATAREAD/srcMagnetogram/, also timing.f90 and ModKind.f90). If there are such mishandled files, you can move the wrong files out of the path, run the script, and put them back.

Examples: 
python mapper.py main.f90 ../GITM2/ #method 1 - check uses and runs 
python mapper.py filter ../../../ #method 2 - check only uses (need to change read_file2( after #contains to read_file( )
"""
#    DART $Id$
#    Copyright (C) 2013 by
#    Alexey Morozov <alexeymor gmail>
__author__ = """Alexey Morozov <alexeymor gmail>"""

#the goal string for hierarchical edge d3.js plot, to be put into he.json is:
# [{"name":"analytics","size":3,"imports":["vis"]},{"name":"vis","size":3,"imports":[]}]

#import networkx as nx #didn't work because want "name" instead of "nodes"
#g=dict([('sape', 4139), ('guido', 4127), ('jack', 4098)]) #didn't work when printing - need " instead of '

#so I'll construct the graph string by hand - not the most elegant solution, but at least I'll learn how to work with strings and stacks in python

import subprocess #for using grep
import random
import sys #for getting command line args
import os # modified http://code.activestate.com/recipes/577027/ 

def build_index(path): #build index first, so that I don't have to build it for every file
    global root_st, dirs_st, names_st
    for root, dirs, names in os.walk(path):
        root_st.append(root)
        dirs_st.append(dirs)
        names_st.append(names)

def file_search(st,p):
    '''find all instances of string s in a file that are not a comment and return whatever word follows after that string before a '(' or a ',' '''
    t = subprocess.check_output("grep " + st + " " + p, shell=True).split('\n') #split by lines
    t1 = []
    for s in t:
        if ('!' not in s[:s.find(st)]): #if this line is not commented out
            t1.append(s) #... add it to the list of lines
    t1 = [s.split() for s in t1] #split each line by spaces
    t1 = filter(None, t1) #remove empty lines
    t1 = [s[1] for s in t1] #keep only the second entry, as we don't care for the word use or call or module or subroutine itself, but only the name that follows it
    for s in t1:
        if ('(' in s):
            t1[t1.index(s)] = s[:s.find('(')] #remove anything after bracket (including)
    for s in t1: #separate loop is needed, because s changes in the loop
        if (',' in s):
            t1[t1.index(s)] = s[:s.find(',')] #remove anything after comma (including)
    t1 = list(set(t1)) #remove duplicates
    t1 = filter(None, t1) #remove empty entries
    return t1


def files_flat():
    ''' extract from each file names of defined and called subroutines and modules '''
    global paths, n_lines, n_bytes, subroutines, calls, modules, uses
    for root, dirs, names in zip(root_st, dirs_st, names_st):
        for n in names:
            if '.f90' in n:
                if (n in paths): #if a file with this name was already seen, ignore subsequent files with the same name (even if contents are different)
                    print os.path.join(root, n), 'is ignored, instead using', paths[n]
                elif ('Jupiter' not in n) and ('LV-426' not in n) and ('Mars' not in n) and ('Titan' not in n) and ('.svn-base' not in n): #if this file was not considered before and jupiter, mars, titan and lv-426 are not in the name (GITM-specific - just want earth-related files) and this is not an SVN copy of some file (DART-specific)
                    paths[n] = os.path.join(root, n)
                    print n #display what file we are about to read
                    n_lines[n] = int(subprocess.check_output("wc -l < " + paths[n], shell=True)) 
                    n_bytes[n] = os.path.getsize(paths[n])
                    try:                                        
                        t = file_search('subroutine',paths[n])
                        del(t[t.index('subroutine')]) #remove "subroutine" (it comes from 'end sub')
                        subroutines[n] = t
                    except:
                        subroutines[n] = []
                    try:
                        calls[n] = file_search('call',paths[n])
                    except:
                        calls[n] = []
                    try:
                        t = file_search('module',paths[n])
                        del(t[t.index('module')]) #remove "module" (it comes from 'end module')
                        modules[n] = t
                    except:
                        modules[n] = []
                    try:
                        uses[n] = file_search('use',paths[n])
                    except:
                        uses[n] = []
#    raise 'File not found'    

def read_file2(hd, stack, visited, files, depends, size, color):
    while len(stack)>0:
        fn = stack.pop() #pop the last file off the stack and analyze it
        visited.append(fn) #put this file on the visited list
        path = paths[fn] #find the file
        loc_dep = [] #local dependencies - dependencies for the current file
#        try:
        files.append(fn) #put this file on the files list
        color.append(str(hex(random.randint(0,16**6-1)))[2:] ) #put this file on the files list
            #print color[-1]
        size.append(n_bytes[fn]) #or can do n_lines[fn]
        print fn
        for call in calls[fn]:
            for child in paths:
                if call in subroutines[child]:
                    if (child not in loc_dep): loc_dep.append(child) #if the child is not listed as a dependency already
                    if (child not in visited) and (child not in stack): stack.append(child) #if child hasn't been visited yet and not scheduled yet
                    break #take only the first file listed in paths that has the needed subroutine
        for use in uses[fn]:
            for child in paths:
                if use in modules[child]:
                    if (child not in loc_dep): loc_dep.append(child) #if the child is not listed as a dependency already
                    if (child not in visited) and (child not in stack): stack.append(child) #if child hasn't been visited yet and not scheduled yet
                    break #take only the first file listed in paths that has the needed subroutine
        depends.append(loc_dep)
        #        except TypeError: #if file does not exist
        #            pass #don't add anything to files or depends or size
        print ' stack right now is ', stack
        (files, depends, size, color)=read_file2(hd, stack, visited, files, depends, size, color) #call r_f recursively until the stack is empty
    return(files, depends, size, color)


def dictinvert(d): #invert dictionary http://code.activestate.com/recipes/252143-invert-a-dictionary-one-liner/
    inv = {}
    for k, v in d.iteritems():
        keys = inv.setdefault(v, [])
        keys.append(k)
    return inv




def find_file(filename, path):
    for root, dirs, names in zip(root_st, dirs_st, names_st):
        if filename in names:
            return os.path.join(root, filename)
#    raise 'File not found'    


def read_file(hd, stack, visited, files, depends, size, color):
    while len(stack)>0:
        fn = stack.pop() #pop the last file off the stack and analyze it
        visited.append(fn) #put this file on the visited list
        path = find_file(fn + '.f90', hd) #find the file
        loc_dep = [] #local dependencies - dependencies for the current file
        try:
            f = open(path, 'r')
            files.append(fn) #put this file on the files list
            color.append(str(hex(random.randint(0,16**6-1)))[2:] ) #put this file on the files list
            #print color[-1]
            size.append(os.path.getsize(path))
            for line in f:
                if ('use' in line): #if a use statement 
                    if ('!' not in line[:line.find('use')+3]): #if not a comment
                        line_trunc = line[line.find('use')+3 : line.find(',')].strip() #remove "use", anything after comma, any spaces
                        path = find_file(line_trunc + '.f90', hd) #find the file
                        try: #only if file exists, add it to list of dependencies 
                            f2 = open(path, 'r')
                            if (line_trunc not in loc_dep): loc_dep.append(line_trunc)
                            if (line_trunc not in visited) and (line_trunc not in stack): 
                                stack.append(line_trunc)
                            f2.close()
                        except:
                            pass
#                if ('call' in line): #if a call statement 
#                    if ('!' not in line[:line.find('call')+4]): #if not a comment
#                        line_trunc = line[line.find('call')+4 : line.find('(')].strip() #remove "call", anything after parenthesis, any spaces
#                        path = find_file(line_trunc + '.f90', hd) #find the file
#                        try: #only if file exists, add it to list of dependencies 
#                            f2 = open(path, 'r')
#                            if (line_trunc not in loc_dep): loc_dep.append(line_trunc)
#                            if (line_trunc not in visited) and (line_trunc not in stack): 
#                                stack.append(line_trunc)
#                            f2.close()
#                        except:
#                            pass
            depends.append(loc_dep)
            f.close()
        except TypeError: #if file does not exist
            pass #don't add anything to files or depends or size
        print ' stack right now is ', stack
        (files, depends, size, color)=read_file(hd, stack, visited, files, depends, size, color) #call r_f recursively until the stack is empty
    return(files, depends, size, color)


#contains (the imports and function definitions are over, so start the actual execution)

fn = sys.argv[1] #where do I start from (find all use statements in this file and then go into those files and so on)?, appends .f90
hd = sys.argv[2] #where can I search for *.f90 files?

root_st = [] #st for store, global variable
dirs_st = [] 
names_st = [] 

paths = {} #path dictionary - given name return path 
n_lines = {} #number of lines dict
n_bytes = {} #number of bytes

subroutines = {} #subroutines defined in this file
calls = {} #calls to subroutines made in this file

modules = {} #modules defined in this file
uses = {} #modules used in this file

build_index(hd) #build index first, so that I don't have to build it for every file

print len(root_st), len(dirs_st), len(names_st)

files_flat()

(files, depends, size, color) = read_file2(hd, [fn], [], [], [], [], [])
#files is a list containing file names used in the particular code
#depends is a list containing lists of depencies for each file in "files"
#size is a list of file sizes for each file in "files"


#merge the lists into a string understood by Graphviz's "dot"
#digraph G { main -> parse; main -> init; }
string = 'digraph G { size="7.75,10.25"; rankdir=LR; weight=1.2; nodesep=0.1; \n'
ranks = {fn:0} #dictionary of ranks for all files in files, first file gets rank of 0, every file used directly in the first file gets rank 1, every file used in file which is used in first file gets rank 2 and so on recursively
i=1 #number of edges
for f,d,s,c in zip(files, depends, size, color):
    for d1 in d: #flatten the dependencies first into "vis","hi" form
        string = string + '"'+f+'"' + '->' + '"'+d1+'"' + '[color="#' + c + '"];' 
        i=i+1
        if (d1 not in ranks): ranks[d1]=ranks[f]+1
    string = string + '\n' 

inv_ranks = dictinvert(ranks)
for r in inv_ranks:
    s=''
    for f in inv_ranks[r]:
        s = s + '"' + f + '" [color="#' + color[files.index(f)] + '"] '
    string = string + '{rank=same ' + s + '};\n'

string = string[:-2] + '}' #remove the last carriage return and semicolon and add the last bracket
f = open('dot.gv', 'w')
f.write(string)
f.close()
print 'number of nodes is', len(files) #print number of nodes
print 'number of edges is', i #print number of edges



#merge the lists into a JSON-like string: (for hierarchical edge layout like in http://mbostock.github.io/d3/talk/20111116/bundle.html 
# [{"name":"analytics","size":3,"imports":["vis"]},{"name":"vis","size":3,"imports":[]}]
string = '['
for f,d,s in zip(files, depends, size):
    df = ''
    for d1 in d: #flatten the dependencies first into "vis","hi" form
        df = df + '"' + d1 + '",'
    df = df[:-1] #remove the last comma
    string = string + '{"name":"' + f + '","size":' + str(s) + ',"imports":[' + df + ']},\n' 

string = string[:-2] + ']' #remove the last carriage return and comma and add the last bracket
#print string
f = open('d.json', 'w')
f.write(string)
f.close()


#merge the lists into a JSON-like string: (for hierarchical edge layout like in http://mbostock.github.io/d3/talk/20111116/bundle.html 
# [{"name":"analytics","size":3,"imports":["vis"]},{"name":"vis","size":3,"imports":[]}]
string = '['
for f,d,s in zip(files, depends, size):
    df = ''
    for d1 in d: #flatten the dependencies first into "vis","hi" form
        df = df + '"' + d1[:d1.find('.')] + '",'
    df = df[:-1] #remove the last comma
    string = string + '{"name":"' + f[:f.find('.')] + '","size":' + str(s) + ',"imports":[' + df + ']},\n' 

string = string[:-2] + ']' #remove the last carriage return and comma and add the last bracket
#print string
f = open('g.json', 'w')
f.write(string)
f.close()


#merge the lists into a more complex JSON-like string: (for hierarchical edge layout like in http://mbostock.github.io/d3/talk/20111116/bundle.html
# [{"name":"flare","size":3,"imports":["flare.vis"]},{"name":"flare.vis","size":3,"imports":[]}]
long_names={fn:fn[:fn.find('.')]} #store filenames in flare.vis.a.b.c format (flare depends on vis which depends on a, which depends on b ... on c and throw away all the ".f90" extensions.
string = '['
for f,d,s in zip(files, depends, size):
    df = ''
    for d1 in d: #flatten the dependencies first into "vis","hi" form
        if (d1 not in long_names): long_names[d1]=long_names[f]+'.'+d1[:d1.find('.')]
        df = df + '"' + long_names[d1] + '",'
    df = df[:-1] #remove the last comma
    string = string + '{"name":"' + long_names[f] + '","size":' + str(s) + ',"imports":[' + df + ']},\n' 

string = string[:-2] + ']' #remove the last carriage return and comma and add the last bracket
#print string
f = open('gl.json', 'w')
f.write(string)
f.close()


#merge the lists into another JSON-like string: (for forced layout like here http://bl.ocks.org/mbostock/1153292
# [{source: "main", target: "ModInputs"},{source: "main", target: "ModTime"}]
string = '['
for f,d,s,c in zip(files, depends, size, color):
    df = ''
    for d1 in d: #flatten the dependencies
        df = df + '{source: "' + f + '", target: "' + d1 + '", color: "#' + c + '"},\n'
    string = string + df 
string = string[:-2] + ']' #remove the last carriage return and comma and add the last bracket

#print string

f = open('forced.json', 'w')
f.write(string)
f.close()

# now create a css file for colors of nodes 
# .main {fill: #fdd700;} .ModInputs {fill: #fdd700;}
string = ''
for f, s, c in zip(files, size, color):
    string = string + '.' + f + ' {fill: #' + c + '; r: ' + str(s) + ';}\n' 

#print string

f = open('forced.css', 'w')
f.write(string)
f.close()


#mv he.json he1.json 
#mv dot.gv dot1.gv 
#mv forced.css forced1.css 
#mv forced.json forced1.json 
#dot -Teps dot1.gv -o gitm.eps

#mv he.json he2.json 
#mv dot.gv dot2.gv 
#mv forced.css forced2.css 
#mv forced.json forced2.json 
#dot -Teps dot2.gv -o dart.eps

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
