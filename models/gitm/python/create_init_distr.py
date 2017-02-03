'''
write a sequence of argv[1] normally distributed random numbers with mean argv[2] and std.dev argv[3] into argv[4] (ASCII text file)
Example:
python create_init_distr.py 20 -16.44 0.3 fens.txt
'''
from math import *
import random
import sys #for getting command line args

noe = int(sys.argv[1]) #number of ensemble members
mu  = float(sys.argv[2]) #mean of the normal distribution
sig = float(sys.argv[3]) #standard deviation of the normal distribution
f = open(sys.argv[4], 'w') #open the desired text file to write output in there

for i in range(noe):
    #r = random.gauss(log10(4e-17), log10(8e-17/4e-17)) #-16.4+-0.3 seems reasonable, 4e-17 is mu, 8e-17 is mu+1sigma
    r = random.gauss(mu, sig) 
    f.write(str(r)+'\n')

f.close()
    
# DART $Id$
# from Alexey Morozov
#
# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$
