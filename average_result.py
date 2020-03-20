#!/usr/bin/python
import sys
import os
import os.path
import commands
import time
import math
import copy
import pprint
import collections
import numpy  as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

try:
        infilename = sys.argv[1]

except:
     print "Usage:",sys.argv[0], "infile"; sys.exit(1)

##### Variable Initializations ##########
ifile = open(infilename,'r') # open file for reading



# duplicate has been removed

generalstates      = ['bound', 'nonspecific', 'dissociated', 'specificparallel']



transitiontimes  = []
transitiontimesall = []
# initialize the region

#if initialstate in allowablepair:
#   i = 3    #  
#else:
#   i = 0       
#if finalstate in allowablepair:
#   j = 4   #  
#else:
#   j = 1  
counter = float(0)
totoaltime = float(0)

for line in ifile:
    columns = line.split() 
    if columns[0][0] != '#':
                currenttime   = float(columns[3])
                totoaltime    = totoaltime+currenttime
                counter       = counter+1
                transitiontimesall.append(currenttime)   #  average and stderr  not necessart within range


print    'average  (ns)                             occur', 
print str(totoaltime/(counter*1000)),               counter
tt = np.array(transitiontimesall)
#print 'av and stderr (ns)', np.mean(tt)/1000, np.std(tt)/1000















