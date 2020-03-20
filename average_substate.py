#!/usr/bin/python
import sys,string
import numpy as np
import math
import collections
from collections import OrderedDict

#try:
#     infilenumber = int(sys.argv[1])
#except:
#     print "Usage:",sys.argv[0], "infile minv1 maxv1 minv2 maxv2 devisions1 devisions2 frame outfile"; sys.exit(1)

###  input reading

infilename = str(sys.argv[1])
outfilename_ratio = str(sys.argv[2])
outfilename_total = str(sys.argv[3])
#simulation_number = int(sys.argv[3])

try:
   cutoff  = float(sys.argv[3])
except:
   cutoff  = float(0) 

value_array          = OrderedDict()
average_value_array  = OrderedDict()

value_array['antiparallel even even 0']  = OrderedDict()
value_array['antiparallel odd odd 0']    = OrderedDict()
value_array['antiparallel even even 2']  = OrderedDict()
value_array['antiparallel even even -2'] = OrderedDict()
value_array['antiparallel even odd 2']   = OrderedDict()
value_array['antiparallel even odd -2']  = OrderedDict()
value_array['antiparallel odd even 2']   = OrderedDict()
value_array['antiparallel odd even -2']  = OrderedDict()
value_array['parallel even even 0']      = OrderedDict()
value_array['parallel odd odd 0']        = OrderedDict()
value_array['parallel odd odd 2']        = OrderedDict()
value_array['parallel odd odd -2']       = OrderedDict()
value_array['parallel even odd 2']       = OrderedDict()
value_array['parallel even odd -2']      = OrderedDict()
value_array['parallel odd even 2']       = OrderedDict()
value_array['parallel odd even -2']      = OrderedDict()

average_value_array['antiparallel even even 0']  = OrderedDict()
average_value_array['antiparallel odd odd 0']    = OrderedDict()
average_value_array['antiparallel even even 2']  = OrderedDict()
average_value_array['antiparallel even even -2'] = OrderedDict()
average_value_array['antiparallel even odd 2']   = OrderedDict()
average_value_array['antiparallel even odd -2']  = OrderedDict()
average_value_array['antiparallel odd even 2']   = OrderedDict()
average_value_array['antiparallel odd even -2']  = OrderedDict()
average_value_array['parallel even even 0']      = OrderedDict()
average_value_array['parallel odd odd 0']        = OrderedDict()
average_value_array['parallel odd odd 2']        = OrderedDict()
average_value_array['parallel odd odd -2']       = OrderedDict()
average_value_array['parallel even odd 2']       = OrderedDict()
average_value_array['parallel even odd -2']      = OrderedDict()
average_value_array['parallel odd even 2']       = OrderedDict()
average_value_array['parallel odd even -2']      = OrderedDict()

for resid in value_array :

    value_array[resid][6] = []
    value_array[resid][5] = []
    value_array[resid][4] = []
    value_array[resid][3] = []
    value_array[resid][2] = []
    value_array[resid][1] = []
    value_array[resid][0] = []

    average_value_array[resid][6] = 0
    average_value_array[resid][5] = 0
    average_value_array[resid][4] = 0
    average_value_array[resid][3] = 0
    average_value_array[resid][2] = 0
    average_value_array[resid][1] = 0
    average_value_array[resid][0] = 0

formarted_name = OrderedDict()

formarted_name['antiparallel even even 0']  = 'antiparallel even even 0 '
formarted_name['antiparallel odd odd 0']    = 'antiparallel odd odd 0   '
formarted_name['antiparallel even even 2']  = 'antiparallel even even 2 '
formarted_name['antiparallel even even -2'] = 'antiparallel even even -2'
formarted_name['antiparallel even odd 2']   = 'antiparallel even odd 2  ' 
formarted_name['antiparallel even odd -2']  = 'antiparallel even odd -2 ' 
formarted_name['antiparallel odd even 2']   = 'antiparallel odd even 2  ' 
formarted_name['antiparallel odd even -2']  = 'antiparallel odd even -2 ' 
formarted_name['parallel even even 0']      = 'parallel even even 0     '    
formarted_name['parallel odd odd 0']        = 'parallel odd odd 0       '    
formarted_name['parallel odd odd 2']        = 'parallel odd odd 2       '    
formarted_name['parallel odd odd -2']       = 'parallel odd odd -2      '    
formarted_name['parallel even odd 2']       = 'parallel even odd 2      '    
formarted_name['parallel even odd -2']      = 'parallel even odd -2     '    
formarted_name['parallel odd even 2']       = 'parallel odd even 2      '    
formarted_name['parallel odd even -2']      = 'parallel odd even -2     '   

ifile = open(infilename,'r') # open file for reading
for line in ifile:
    columnes = line.split()
    if len(columnes) >= 11  and columnes[0][0] != '#' and columnes[0][0] != '!' and columnes[0][0] != '@' :
        resid     = str(columnes[0] + ' ' + columnes[1] + ' ' + columnes[2] + ' ' + columnes[3] )
                     #  antiparallel           even                even                 0 
                     #      orientation        incoming            core               shift   
        value_FCL_6     = float(columnes[4]) 
        value_FCL_5     = float(columnes[5])
        value_FCL_4     = float(columnes[6])
        value_FCL_3     = float(columnes[7])
        value_FCL_2     = float(columnes[8])
        value_FCL_1     = float(columnes[9])
        value_FCL_0     = float(columnes[10])
        value_array[resid][6].append(value_FCL_6)     
        value_array[resid][5].append(value_FCL_5)
        value_array[resid][4].append(value_FCL_4)
        value_array[resid][3].append(value_FCL_3)
        value_array[resid][2].append(value_FCL_2)
        value_array[resid][1].append(value_FCL_1)
        value_array[resid][0].append(value_FCL_0)

ifile.close
   

##### Variable Initializations ##########

#outfilename = 'temp_average.dat'
ofile_ratio = open(outfilename_ratio,'w') # open file for writing
ofile_total = open(outfilename_total,'w') # open file for writing

kB = 3.2976268E-24     # cal/K
#kB = 1.3806488E-23    # j/K
An = 6.02214179E23
T = float(300)
R = 1.987  # (cal/mol.degree)

total_average_residuece_time = float(0)

for resid in value_array:
    for FCL in value_array[resid] :
         mean  = np.mean(value_array[resid][FCL])
         #error = np.std(value_array[resid])
         total_average_residuece_time    += mean
         average_value_array[resid][FCL]  = mean
            
for resid in average_value_array:
    ofile_ratio.write(str(formarted_name[resid]) + "      0    " )
    ofile_total.write(str(formarted_name[resid]) + "      0    " )

    for FCL in average_value_array[resid] :
         temp_average = average_value_array[resid][FCL]/total_average_residuece_time
         temp_average = format(temp_average, '.4f')

         temp_total   = format(average_value_array[resid][FCL]/1000, '.2f')

         ofile_ratio.write( str( temp_average ) + "    " )
                            #  we use ratio, for finaly this should be energy
         ofile_total.write( str( temp_total ) + "    " )
                            # output total time in ns, which will be connected to nonspecific state
    ofile_ratio.write( "# "  +"\n") 
    ofile_total.write( "# "  +"\n")



ofile_ratio.close()
ofile_total.close()









