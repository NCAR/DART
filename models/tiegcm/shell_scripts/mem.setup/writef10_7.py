#!/usr/bin/python3
import sys

original_stdout = sys.stdout

def replace_f10_7(string):
    f10_7 = string.split()[2]
    return string.replace(f10_7, str(netcdf)) 

netcdf=sys.argv[1]

with open('tiegcm_res5.0.inp') as f:
    lines = f.readlines()
f.close

with open('tiegcm_res5.0.inp', 'w') as g:
    sys.stdout = g 
    for line in lines:
      if 'F107 ' in line:
         print(replace_f10_7(line), end='')
      else:
         print(line, end='')
   
    sys.stdout = original_stdout

g.close
