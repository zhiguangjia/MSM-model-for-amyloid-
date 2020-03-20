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

try:
        infilename = sys.argv[1]
        infilename2 = sys.argv[2]
        outfilename = sys.argv[3]
        simulation_number = int(sys.argv[4])
except:
     print "Usage:",sys.argv[0], "infile (average time)   infile2 (total time)   outputfile"; sys.exit(1)

##### Variable Initializations ##########
ifile = open(infilename,'r') # open file for reading
ifile2 = open(infilename2,'r') # open file for reading
ofile = open(outfilename,'w') # open file for writing

ofile.write( "event key       time of entered          average lifetime          total  residence time  (ns)            total  residence time (% of total time)    "  +"\n")


average_time_key = {}

for line in ifile:
    columns = line.split()  

    if columns[0][0] != '#':
                key   = str(columns[0] + ' ' + columns[1] + ' ' + columns[2] + ' ' + columns[3])
                # average_time_key[key] = []
                # build the key dictenory  1-d is enough (currently)
                counter = float(0)
                totaltime = float(0) 
                t        = 0
                counter = float(0)

                while t < simulation_number -1:  
 
                  currenttime    = float(columns[t+4])
                  if   currenttime != 0 :
                      # skip if the MSM simulation not experence this state 
                      totaltime       = currenttime  + totaltime   
                      counter = counter+1 
                  t = t+1
   
                if  counter  != 0 :   
                     #averagetime   = totaltime /float(counter)
                     averagetime   = totaltime /float(simulation_number -1)
                    
                     average_time_key[key]  = averagetime
                else : 
                     averagetime   = 0
                     average_time_key[key]  = float(0)
                     # skip if the MSM simulation not experence this state 

averaged_totaltime = {}
crosstime          = {}


for line in ifile2:
    columns = line.split()  

    if columns[0][0] != '#':
                key   = str(columns[0] + ' ' + columns[1] + ' ' + columns[2] + ' ' + columns[3])

                counter = float(0)
                totaltime = float(0) 
                t        = 0
                counter = float(0)

                while t < simulation_number -1:  

                  currenttime    = float(columns[t+4])
                  #print currenttime ,  key
                  if   currenttime != 0 :
                      # skip if the MSM simulation not experence this state 
                      totaltime       = currenttime  + totaltime   
                      counter = counter+1 
                  t = t+1
                  
                if  counter  != 0 :   
                     #averagetime   = totaltime /float(counter)
                     averaged_totaltime[key]   = totaltime /float(simulation_number-1)

                     crosstime[key]     = (averaged_totaltime[key]/(average_time_key[key]))
                       # math.ceil take upper intact    math.floor take lower one   round add 1 when > 0.5
                     print key,  (average_time_key[key])      
                   
                     #crosstime[key]      = ceil(counter/float(samulation_number-1))    #  this is rong, enter more than one still count as one

                else : 
                     averaged_totaltime[key]  = 0
                     crosstime[key]           = 0
                #print key, ' cout time  ',  counter

simulationm_time = float(0) 
for key in  averaged_totaltime:
    simulationm_time +=  averaged_totaltime[key]                            
 
print '  total time ', simulationm_time
               
#for key in  average_time_key:
# for clear order
ifile.close()
ifile = open(infilename,'r') # open file for reading

for line in ifile:
    columns = line.split()

    if columns[0][0] != '#':
                key   = str(columns[0] + ' ' + columns[1] + ' ' + columns[2] + ' ' + columns[3])

                

                ofile.write(key + "           " + str( format( average_time_key[key], '.4f')  )  + "              " +  str( format(crosstime[key], '.4f') ) + "                     " +  str( format( averaged_totaltime[key]/1000, '.2f') ) + "                 " +  str( format( averaged_totaltime[key]/simulationm_time, '.4f') ) + "   \t  "  +"\n")

                #                        # math.ceil take upper intact    math.floor take lower one   round add 1 when > 0.5










