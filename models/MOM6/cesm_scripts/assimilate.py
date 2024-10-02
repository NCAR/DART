#!/usr/bin/env python3
import sys
from CIME.case import Case

def assimilate(caseroot):
    with Case(caseroot) as case:
         print ("rundir is ",case.get_value("RUNDIR"))
         print ("hello Helen")

if __name__ == "__main__":
    caseroot = sys.argv[1]
    assimilate(caseroot)
