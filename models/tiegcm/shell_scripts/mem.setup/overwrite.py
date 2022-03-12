#!/usr/bin/python3
import sys

original_stdout = sys.stdout

def incr_day_hour_min(string):
    (day, hour, second) = string.split()[2:5]
    return string.replace(day, str(int(day)+1))

def incr_day(string):
    day = string.split()[2]
    return string.replace(day, str(int(day)+1))

with open('tiegcm_res5.0.inp') as f:
    lines = f.readlines()
f.close

with open('tiegcm_res5.0.inp', 'w') as g:
    sys.stdout = g
    for line in lines:

      if 'SOURCE =' in line:
         continue 
      elif 'SOURCE_START =' in line:
         continue
      elif ' START ' in line:
         print(incr_day_hour_min(line), end='')
      elif ' STOP ' in line:
         print(incr_day_hour_min(line), end='')
      elif 'SECSTART ' in line:
         print(incr_day_hour_min(line), end='')
      elif 'SECSTOP ' in line:
         print(incr_day_hour_min(line), end='')
      elif 'START_DAY ' in line:
         print(incr_day(line), end='')
      else:
         print(line, end='')
   
    sys.stdout = original_stdout

g.close
