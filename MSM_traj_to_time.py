#!/usr/bin/python
import sys,string
import numpy as np
import math
import collections
import re

from collections import OrderedDict

try:
        infilename = sys.argv[1]

except:
     print "Usage:",sys.argv[0], "infile outfile"; sys.exit(1)
##### Variable Initializations ##########


ifile = open(infilename,'r') # open file for readind

#   1  for x, 2 for  y

#kB = 3.2976268E-24
#An = 6.02214179E23
#T = float(_temp)

##  manually define x coordination , which is  in fact the states

state_time       = OrderedDict()

#   *****  note in later, this file will be modified
# allowableregisters          = ['2-8', '4-6', '6-4', '8-2', '3-7', '5-5', '7-3', '6-2', '4-4', '2-6', '8-4', '6-6', '4-8', '7-2', '5-4', '3-6', '7-4', '5-6', '3-8', '8-3', '6-5', '4-7', '6-3', '4-5', '2-7']
# allowableregisters_parallel =  ['2-2', '4-4', '6-6', '8-8', '3-3', '5-5', '7-7', '2-4', '4-6', '6-8', '3-4', '5-6', '7-8', '4-2', '6-4', '8-6', '3-2', '5-4', '7-6', '2-3', '4-5', '6-7', '4-3', '6-5', '8-7'  ]



inregisters   = ['2-8', '4-6', '6-4', '8-2', '3-7', '5-5', '7-3' ] 
inregisters_o = [ '3-7', '5-5', '7-3']
inregisters_e = ['2-8', '4-6', '6-4', '8-2']
anti_ee_2     = ['6-2', '4-4', '2-6']
anti_ee_n2    = ['8-4', '6-6', '4-8']
anti_oe_2     = ['7-2', '5-4', '3-6']
anti_oe_n2    = ['7-4', '5-6', '3-8']
anti_eo_2     = ['6-3', '4-5', '2-7']
anti_eo_n2    = ['8-3', '6-5', '4-7']
para_eo_0     = ['2-2', '4-4', '6-6', '8-8']
para_eo_2     = ['2-4', '4-6', '6-8']
para_eo_n2    = ['4-2', '6-4', '8-6']
para_oe_0     = ['3-3', '5-5', '7-7']
para_oo_2     = ['3-4', '5-6', '7-8']
para_oo_n2    = ['3-2', '5-4', '7-6']
para_ee_2     = ['2-3', '4-5', '6-7']
para_ee_n2    = ['4-3', '6-5', '8-7']

state_time['nonspecific_b']              =  0.0
state_time['nonspecific_13']             =  0.0
state_time['nonspecific_12']             =  0.0   
state_time['nonspecific_10']             =  0.0   
state_time['nonspecific_8']              =  0.0   
state_time['nonspecific_7']              =  0.0   
state_time['nonspecific_6']              =  0.0   
state_time['nonspecific_5']              =  0.0   
state_time['nonspecific_4']              =  0.0   
state_time['nonspecific_3']              =  0.0   
state_time['nonspecific_2']              =  0.0   
state_time['nonspecific_1']              =  0.0   
state_time['nonspecific_0']              =  0.0   
state_time['nonspecific_14']             =  0.0   
state_time['nonspecific_11']             =  0.0   
state_time['nonspecific_9']              =  0.0   


state_time['anti_ee_2']    = 0.0 
state_time['anti_ee_n2']   = 0.0 
state_time['anti_oe_2']    = 0.0 
state_time['anti_oe_n2']   = 0.0 
state_time['anti_eo_2']    = 0.0 
state_time['anti_eo_n2']   = 0.0 
state_time['para_eo_0']    = 0.0 
state_time['para_eo_2']    = 0.0 
state_time['para_eo_n2']   = 0.0 
state_time['para_oe_0']    = 0.0 
state_time['para_oo_2']    = 0.0 
state_time['para_oo_n2']   = 0.0 
state_time['para_ee_2']    = 0.0 
state_time['para_ee_n2']   = 0.0 

state_time['inregisters_e']  = 0.0
state_time['inregisters_o']  = 0.0


begin_flat    = 'off' 
first_cycle   = 'yes'
current_time  = 0.0
last_time     = 0.0
registered    = 'no'
current_orientation = 'none'
current_state       = 'none'
   
for line in ifile:
    columnes = line.split()
    if len(columnes) >= 2  and columnes[0][0] != '#' and columnes[0][0] != '!' and columnes[0][0] != '@' :
      #  print line   

        if str(columnes[0] + ' ' + columnes[1] ) == 'current state' :

      #      print('enter 2nd elif')
            current_state = str(columnes[2])

            #* this is marker for  _b, nonspecifc ...
            registered    = 'no'   # may modify in next elif


        elif len(columnes) >= 3 and  str(columnes[0] + ' ' + columnes[1]  + ' ' + columnes[2]) == 'new cycle, previous:' :

       #     print('enter 3rd elif')

            #* this is in fact the 1st line ...
            #  not read as 1st because it need modify previous elif
 
            new_columnes  = re.split(r'[;,\[\]\s]\s*', line)   # split with delimiters comma, semicolon and space

            current_strand = [int(new_columnes[4]), int(new_columnes[5]), int(new_columnes[6]), int(new_columnes[7]), int(new_columnes[8]), int(new_columnes[9]), int(new_columnes[10])]
            
         #   print(current_strand,  '     de bug,  can it find here ?' )
            if current_strand.count(1) > 0:
                registered = 'yes'
                
                current_state = str(current_strand.count(0))   #  for Hbond state, state in dex is FCL
              #  current_Hbond = float(current_strand.count(0))
        #        print('tttt')

        elif len(columnes) >= 3 and ( registered == 'yes' and ( str(columnes[-2]) == 'specificstate' ) ) :

        #    print('enter 4th elif')

            #* this is marker for Hbond state and register

            current_orientation      = columnes[-3]
            current_registered_state = columnes[-4]  # e.g. 2-8 ...

        elif len(columnes) == 3 and ( registered == 'yes' and str(columnes[0] + ' ' + columnes[1]) ==  'cureent orientation' ) :
        # another place may define orientation
            if columnes[2] != 'none':
                # this happens when last H break
                current_orientation      = columnes[2]

        elif len(columnes) >= 3 and  str(columnes[0]) == 'strandpointleft'  :    

         #   print('enter 5th elif')

            current_left_FCL         = float(columnes[1])  - 2
            if str(columnes[2]) == 'strandpointright' :
                current_right_FCL        = 8  - float(columnes[3])

        elif len(columnes) >= 2 and str(columnes[0]) == 'strandpointright' :   

          #  print('enter 6th elif')

            current_right_FCL        = 8  - float(columnes[1])    


            #  when single bond:
            #      strandpointleft 6 strandpointright 6
            #  
            #  when more bond:
            #      strandpointleft  4 strandleftend   2 
            #      strandpointright  6 strandrightend 8
            #  


        elif len(columnes) >= 3 and ( registered == 'yes' and ( str(columnes[0]) == 'SS' ) ) :   
            current_left_SS          = float(columnes[2]) 
            current_right_SS         = float(columnes[4])  

        elif len(columnes) >= 3 and ( registered == 'yes' and ( str(columnes[0]) == 'SS_left' ) ) :
            current_left_SS          = float(columnes[1])
            current_right_SS         = float(columnes[3])

            #*  a typo in origin file, SS_left became SS left

        elif columnes[0] == 'Time(ps):' :

          #  last_time    = current_time
            current_time = float(columnes[1])/1000   # ps to ns

        elif len(columnes) == 4 and ( columnes[0] ==  'Time' and columnes[1] ==  'elapsed(ps)' ):
            timeelapsed = float(columnes[3])/1000   # ps to ns


        elif str(columnes[0] + ' ' + columnes[1] ) == 'end cycle'   :
            #* cycle end, begin write
           # print (current_state, 'debugdebugdebugdebugdebugdebug' ) 

            if   current_state == 'dissociate' :
                 pass

            elif current_state == 'nonspecific_b' :

                 state_time['nonspecific_b']  += timeelapsed  

            elif registered  == 'no' and current_state in ['nonspecific_14', 'nonspecific_11', 'nonspecific_9' ] :

                 state_time[current_state]  += timeelapsed  

            elif registered  == 'no' and current_state not in ['nonspecific_14', 'nonspecific_11', 'nonspecific_9' ] :

                 state_time[current_state]  += timeelapsed

            elif current_orientation  == 'antiparallel' and ( current_registered_state in inregisters_e ):
                 state_time['inregisters_e']  += timeelapsed

            elif current_orientation  == 'antiparallel' and ( current_registered_state in inregisters_o ):
                 state_time['inregisters_o']  += timeelapsed

            elif current_orientation  == 'antiparallel' :

                 if current_registered_state in anti_ee_2 :

                     state_time['anti_ee_2']  += timeelapsed

                 elif current_registered_state in anti_ee_n2 :  

                     state_time['anti_ee_n2']  += timeelapsed

                 elif current_registered_state in anti_oe_2 : 

                     state_time['anti_oe_2']  += timeelapsed
                 elif current_registered_state in anti_oe_n2 : 

                     state_time['anti_oe_n2']  += timeelapsed

                 elif current_registered_state in anti_eo_2 : 

                     state_time['anti_eo_2']  += timeelapsed

                 elif current_registered_state in anti_eo_n2 : 

                     state_time['anti_eo_n2']  += timeelapsed

                 else:
                     print('antiparallel', current_registered_state, 'not found') 


            elif current_orientation  == 'parallel' :

                 if current_registered_state in para_eo_0 :
                     state_time['para_eo_0']  += timeelapsed

                 elif current_registered_state in para_eo_2 : 
                     state_time['para_eo_2']  += timeelapsed

                 elif current_registered_state in para_eo_n2 : 
                     state_time['para_eo_n2']  += timeelapsed

                 elif current_registered_state in para_oe_0 : 
                     state_time['para_oe_0']  += timeelapsed

                 elif current_registered_state in para_oo_2 : 
                     state_time['para_oo_2']  += timeelapsed

                 elif current_registered_state in para_oo_n2 : 
                     state_time['para_oo_n2']  += timeelapsed

                 elif current_registered_state in para_ee_2 : 
                     state_time['para_ee_2']  += timeelapsed

                 elif current_registered_state in para_ee_n2 : 
                     state_time['para_ee_n2']  += timeelapsed

                 else:
                     print('parallel', current_registered_state, 'not found')

            else:

                 print('miss:', current_state, current_orientation, current_registered_state )



            #  now write general file: 

    else :
        pass

print('# unit ns')
for state in state_time:
    print state, '                           '  , state_time[state] 













