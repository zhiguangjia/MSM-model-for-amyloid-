#!/usr/bin/python

import sys
import os
import os.path
import commands
import time as TIME
import math
import copy
import pprint
import random
####################################################################################
#
#  Global variables 
#  These are put here as they need manually modified for each protein 
#
####################################################################################


#  abbreviations:
##
#  1  'Ts' the residue in (incoming) strand 
#     'Tc' the residue in core
#
#  2  'shift'  the relative movement between strand and core, an in-register state should have 0  shift
#
#  3  task/type of MSM:  onec2b   incoming strand has a single hbond pair with core, MSM exaime how long take it to reach full bound state (or dissociate state)
#                        ub2b     incoming strand start from b circle, and can reac full bound state or dissociate state
#                        b2ub     incoming strand start from full bound state, examine how long take it to dissociated (= Toff)            
#
#  4  'SS' secondary structure, currently only beta non-beta transition is enabled
#
#  5  'FCL' free chain length (how many 'free/unbound' residues beteen last bbound one and terminus of incoming peptide )

###########  overall flow:
#
#
# 1 from input, find target (ub2b ... b2ub ...)   
#               find inital state and registry (bound or unbound ?  which hydrogen bond formed? )
#               find final state 'quit condition' (e.g. fully bound/dessociate)
#               
# 2 read the transition rate
#    2.1 )  the raw input contains inital sub-state (e.g. single bond 2-8), final state , average transtion time, occures (how many time recorded)
#    2.2 )  in the raw file, we also has fine classfication the of the transitions, such as inregistered/misregistered
#                                                                                        or NB to different signle H bond state
#           These has been averaged before next step, as shown in SI, this classfication is overkilled 
#    2.3 )  also, during this process, we now how many non-spcific state we have
#           this is from cluster analysis and could not be pre-defined 
#    2.4 )  the final rate data is a dictionary. The key contains the inial and final sub state (e.g., for a bound )
#    2.5 )  SS transition rate will also be read in this step, transitions with different SS state could not be comnbined (manuscript SI show SS state have large impact on the transition rate)
#
# 3 build transition matrix
#    3.1 )  As first step, all sub-states is defined
#    3.2 )  fill the rate from step 2 to here, note different transistions may have same ratte (as long as the key are same ). 
#           Note, for a symmetric peptide, this will be a very common situation
#    
# 4 MSM
#    from the intial state, start simulation. 
#


####################################################################################


### pdb information, when write trajectory option activated, this information will be used as templete for writing pdb  (see func 'writepdb')

model = 1        #frame number
pdbstrand = []   #pdb templates for trajectory output
pdbstrand_parallel  = []
pdbcore = []
pdbstrand.append("ATOM    158  CA  LYS     2       0.000   2.000   0.000  1.00  0.00      A8")
pdbstrand.append("ATOM    171  CA  LEU     3       1.200   2.000   0.000  1.00  0.00      A3")
pdbstrand.append("ATOM    180  CA  VAL     4       2.400   2.000   0.000  1.00  0.00      A3")
pdbstrand.append("ATOM    188  CA  PHE     5       3.600   2.000   0.000  1.00  0.00      A3")
pdbstrand.append("ATOM    200  CA  PHE     6       4.800   2.000   0.000  1.00  0.00      A3")
pdbstrand.append("ATOM    212  CA  ALA     7       6.000   2.000   0.000  1.00  0.00      A3")
pdbstrand.append("ATOM    218  CA  GLU     8       7.200   2.000   0.000  1.00  0.00      A2")
pdbstrand_parallel.append("ATOM    158  CA  LYS     2       7.200   2.000   0.000  1.00  0.00      A8")
pdbstrand_parallel.append("ATOM    171  CA  LEU     3       6.000   2.000   0.000  1.00  0.00      A3")
pdbstrand_parallel.append("ATOM    180  CA  VAL     4       4.800   2.000   0.000  1.00  0.00      A3")
pdbstrand_parallel.append("ATOM    188  CA  PHE     5       3.600   2.000   0.000  1.00  0.00      A3")
pdbstrand_parallel.append("ATOM    200  CA  PHE     6       2.400   2.000   0.000  1.00  0.00      A3")
pdbstrand_parallel.append("ATOM    212  CA  ALA     7       1.200   2.000   0.000  1.00  0.00      A3")
pdbstrand_parallel.append("ATOM    218  CA  GLU     8       0.000   2.000   0.000  1.00  0.00      A2")
pdbcore.append("ATOM    234  CA  LYS     2       7.200   0.000   0.000  1.00  99.0      A2")
pdbcore.append("ATOM    247  CA  LEU     3       6.000   0.000   0.000  1.00  0.00      A4")
pdbcore.append("ATOM    256  CA  VAL     4       4.800   0.000   0.000  1.00  0.00      A4")
pdbcore.append("ATOM    264  CA  PHE     5       3.600   0.000   0.000  1.00  0.00      A4")
pdbcore.append("ATOM    276  CA  PHE     6       2.400   0.000   0.000  1.00  0.00      A4")
pdbcore.append("ATOM    288  CA  ALA     7       1.200   0.000   0.000  1.00  0.00      A4")
pdbcore.append("ATOM    294  CA  GLU     8       0.000   0.000   0.000  1.00  0.00      A8")

###  predefined registry   

allowableregisters          = ['2-8', '4-6', '6-4', '8-2', '3-7', '5-5', '7-3', '6-2', '4-4', '2-6', '8-4', '6-6', '4-8', '7-2', '5-4', '3-6', '7-4', '5-6', '3-8', '8-3', '6-5', '4-7', '6-3', '4-5', '2-7']
# formart Ts:Tc 
allowableregisters_register = ['0',   '0',    '0',  '0',   '0',    '0',  '0',    '2',  '2',   '2',   '-2',  '-2',   '-2', '2',    '2',   '2',   '-2',  '-2',  '-2',  '-2',  '-2', '-2',   '2',   '2',   '2']
# above added can be call as  allowableregisters_register[allowableregisters.index('8-3')]

inregisters                 = ['2-8', '4-6', '6-4', '8-2', '3-7', '5-5', '7-3']
allowableregisters_parallel          = ['2-2', '4-4', '6-6', '8-8', '3-3', '5-5', '7-7', '2-4', '4-6', '6-8', '3-4', '5-6', '7-8', '4-2', '6-4', '8-6', '3-2', '5-4', '7-6', '2-3', '4-5', '6-7', '4-3', '6-5', '8-7'  ]
allowableregisters_parallel_register = ['0',   '0',   '0',   '0',   '0',   '0',   '0',   '2',   '2',   '2',   '2',   '2',   '2',   '-2',  '-2',  '-2',  '-2',  '-2',  '-2',  '2',    '2',  '2',   '-2',  '-2',  '-2'  ] 
# above added can be call as   allowableregisters_parallel_register[allowableregisters_parallel.index('2-2')]

inregisters_parallel        = ['2-2', '4-4', '6-6', '8-8', '3-3', '5-5', '7-7']

allowablestate = [ 'nonspecific_b', 'dissociate' ]    #  of couse there is a "specific" state, but it is default and this list is used only after state leave "specific"  
                                    #  the  sub-nonspecific states will be filled into array later with their matrix
                                    #  the  nonspecific_b means not contact but in b circle
                                    #  the  dissociate , equlivent with  nonspecific_q  in McCammon model
allowableregister = ['antiparallel even even 0',  'antiparallel odd odd 0',  'antiparallel even even 2',  'antiparallel even even -2',  'antiparallel even odd 2',  'antiparallel even odd -2',  'antiparallel odd even 2',  'antiparallel odd even -2',  'parallel even even 0',  'parallel odd odd 0',  'parallel odd odd 2',  'parallel odd odd -2',  'parallel even odd 2',  'parallel even odd -2',  'parallel odd even 2','parallel odd even -2' ]
#################################################################################################
#function for writing pdb using pdbstrand and pdbcore as templates

def writepdb(strand,  model, time, timeelapsed, specificstate, orientation, strandpointleft, corepointleft):
    j = 0
    while j < 8 :  # repeat a same image for n times so final movie are more smooth
        model += 1
        j    += 1
        pdbfile.write('REMARK ' + str(time) + ' ' + str(timeelapsed) + ' ' + str(specificstate)  + ' ' + orientation + '\n')
        pdbfile.write('MODEL ' + str(model) + '\n')

        if float(strandpointleft)/2 - int(float(strandpointleft/2)) == 0:
                strandfcl =  strandpointleft-2
                currentstrend  =  'even'
        else :
                strandfcl =  strandpointleft-2
                currentstrend  =  'odd'

        if float(corepointleft)/2 - int(float(corepointleft/2)) == 0 and orientation == 'antiparallel':
                corefcl   =  8 - corepointleft 
                currentcore   =  'even'
        elif float(corepointleft)/2 - int(float(corepointleft/2)) != 0 and orientation == 'antiparallel':
                corefcl   =  7 - corepointleft   
                currentcore   =  'odd'
        elif float(corepointleft)/2 - int(float(corepointleft/2)) == 0 and orientation == 'parallel':
                corefcl   =  corepointleft -3                  #   note in paralle, the even correspond odd in antiparallel
                currentcore   =  'even'
        elif float(corepointleft)/2 - int(float(corepointleft/2)) != 0 and orientation == 'parallel':
                corefcl   =  corepointleft -2
                currentcore   =  'odd'



        for i in range(0,len(strand)):
                pdbcolumns = pdbstrand[i].split()
                pdbcolumns_parallel = pdbstrand_parallel[i].split()

                #  x axix  control shift

                if  specificstate in allowablestate :
                      strand_x = float(pdbcolumns[5])
                elif orientation ==  'antiparallel'  :
                      strand_x = float(pdbcolumns[5])
                      shift  =  strandfcl - corefcl
                      strand_x = strand_x - 1.2*shift
                elif orientation ==  'parallel'  :
                      strand_x = float(pdbcolumns_parallel[5])
                      shift  =  strandfcl - corefcl
                      strand_x = strand_x + 1.2*shift
                else :
                      strand_x = float(pdbcolumns[5])


                if strand[i] == 1 and specificstate ==  'specific':
                      pdbcolumns[6] = '1.700'
                elif specificstate in allowablestate and specificstate != 'dissociate':
                      pdbcolumns[6] = '3.700'
                elif specificstate == 'dissociate':
                      pdbcolumns[6] = '5.700'

                # when write pdb, we do not need change x y position unless the orientation changed (then need read a different templete)
 

                if (orientation == 'antiparallel' or  orientation == 'none') and i == 0:
                      pdbcolumns[9] = str(99.0)
                elif orientation == 'parallel' and i == len(strand):
                      pdbcolumns[9] = str(99.0)

                
                 
                pdbfile.write(str(pdbcolumns[0]) + '    ')
                pdbfile.write(str(pdbcolumns[1]) + '  ')
                pdbfile.write(str(pdbcolumns[2]) + '  ')
                pdbfile.write(str(pdbcolumns[3]) + '     ')
                pdbfile.write(str(pdbcolumns[4]) + '      ')
                #  x axix  control shift
                if str(strand_x)[0] != '-':
                      pdbfile.write(' ' + str(strand_x) + '00   ')
                else :
                      pdbfile.write(str(strand_x) + '00   ')
                #  y axix
                pdbfile.write(str(pdbcolumns[6]) + '   ')

                pdbfile.write(str(pdbcolumns[7]) + '  ')
                pdbfile.write(str(pdbcolumns[8]) + '  ')
                pdbfile.write(str(pdbcolumns[9]) + '      ')
                pdbfile.write(str(pdbcolumns[10]))
                pdbfile.write('\n')


        for line in pdbcore:
            pdbfile.write(line + '\n')
        pdbfile.write('ENDMODEL\n')
    return model

######################################################################################
#function that diplays cartoon of peptide when run in verbose mode

def printstrand(strand, restypes, initialcontactstrand):
    for resid in reversed(restypes):
        if resid == 'NOCONTACT':
            print ' ',
        else:
            break
    for i in range(0,len(strand)):
        if strand[i] == 0:
            print '_',
        else:
            print ' ',
    print
    firstbondon = 'true'
    leftoff = 'false'
    rightoff = 'false'
    printstring = ''
    for resid in reversed(restypes):
        if resid == 'NOCONTACT':
            printstring = printstring + '  '
        else:
            break
    for i in range(0,len(strand)):
        if strand[i] == 0:
            printstring = printstring + '  '
        if strand[i] == 1:
            printstring = printstring + '_ '
    if strand[0] == 1 and strand[1] == 0:
        printstring = '/ ' + printstring[2:len(printstring)]
    if strand[len(strand)-1] == 1 and strand[len(strand)-2] == 0:
        printstring = printstring[0:len(printstring)-2] + '\\ '
    for i in range(2,len(printstring)-2):
        if printstring[i-2] == ' ' and printstring[i] == '_' and printstring[i+2] == ' ':
            printstring = printstring[0:i] + '| ' + printstring [i+2:len(printstring)]
        if printstring[i-2] == ' ' and printstring[i] == '_':# and restypes[(i/2)-1] == 'NOCONTACT':
            printstring = printstring[0:i] + '\\ ' + printstring [i+2:len(printstring)]
        if printstring[i+2] == ' ' and (printstring[i] == '_' or printstring[i] == '\\'):
            printstring = printstring[0:i] + '/ ' + printstring [i+2:len(printstring)]
    print printstring
    for resid in restypes:
        if resid == 'NOCONTACT': #and initialcontactstrand % 2 != 1:
            print ' ',
        else:
            break
    for i in range(0, len(strand)):
        print '_',
    print
    print
    return


##########################################################################################
# given the initial contact postiion, fill the "restype" array with all possible contacts

def findrestypes(initialcontactstrand, initialcontactcore, resnames):
    restypes = [resnames[str(initialcontactstrand)] + ':' + resnames[str(initialcontactcore)]]

    residstrand = initialcontactstrand - 1
    residcore = initialcontactcore + 1
    while residstrand >= 2:
        if residcore <= 8:
            restypes.insert(0, resnames[str(residstrand)] + ':' + resnames[str(residcore)])
        else:
            restypes.insert(0, 'NOCONTACT')
        residstrand = residstrand - 1
        residcore = residcore + 1

    residstrand = initialcontactstrand + 1
    residcore = initialcontactcore - 1
    while residstrand <= 8:
        if residcore >= 2:
            restypes.append(resnames[str(residstrand)] + ':' + resnames[str(residcore)])
        else:
            restypes.append('NOCONTACT')
        residstrand = residstrand + 1
        residcore = residcore - 1

    #if strand is bonding using odd resids, ends of peptide are unbound
    #if initialcontactstrand % 2 == 1:
    #    restypes[6] = 'NOCONTACT'
    #    restypes[0] = 'NOCONTACT'

    return restypes

##########################################################################################
#  find the shift of current state

def find_state_shift(key):
    ##  note the key formart is   key = state + ' ' + str(singcontact)  + ' antiparallel'
    columns          =   key.split()                                      #
    orientation      =   columns[2]                                       # 
    tempregister     =   columns[1]      
    if orientation == 'antiparallel' :    
        if  float(tempregister[2:3])/2 - int(float(tempregister[2:3])/2) == 0:
               tempcoreface = "even"               
        else:
               tempcoreface = "odd"
        if  float(tempregister[0:1])/2 - int(float(tempregister[0:1])/2) == 0:
               tempstrandface = "even"
        else:
               tempstrandface = "odd"
        #  now define shift
        tempshift = allowableregisters_register[allowableregisters.index(tempregister)]
  
          

    if orientation == 'parallel' :
        if  float(tempregister[2:3])/2 - int(float(tempregister[2:3])/2) == 0:
               tempcoreface = "odd"
        else:
               tempcoreface = "even"
        if  float(tempregister[0:1])/2 - int(float(tempregister[0:1])/2) == 0:
               tempstrandface = "odd"
        else:
               tempstrandface = "even"
        tempshift = allowableregisters_parallel_register[allowableregisters_parallel.index(tempregister)]
    
    fullyregister =  str( orientation + ' ' + tempstrandface + ' ' + tempcoreface + ' ' + str(tempshift))    
    return fullyregister

def main(model, register):


    f = open(dir + '/timedata', 'w')


#   peptide sequences
    wtresnames = {}
    wtresnames['2'] = 'LYS'
    wtresnames['3'] = 'LEU'
    wtresnames['4'] = 'VAL'
    wtresnames['5'] = 'PHE'
    wtresnames['6'] = 'PHE'
    wtresnames['7'] = 'ALA'
    wtresnames['8'] = 'GLU'
    
    cha19resnames = copy.deepcopy(wtresnames)
    cha19resnames['5'] = 'CHA'

    cha20resnames = copy.deepcopy(wtresnames)
    cha20resnames['6'] = 'CHA'

    cha1920resnames = copy.deepcopy(wtresnames)
    cha1920resnames['5'] = 'CHA'
    cha1920resnames['6'] = 'CHA'

#"dir"
    dirwt = ['wt']
    dircha19 = ['cha19']
    dircha1920 = ['cha1920']
    dircha20 = ['cha20']
    

    if dir in dirwt:
        resnames = wtresnames
    elif dir in dircha19:
        resnames = cha19resnames
    elif dir in dircha1920:
        resnames = cha1920resnames
    elif dir in dircha20:
        resnames = cha20resnames
    else:
        print 'special case unrecognized'
        sys.exit()

#initial assignment of states
    #strand: states of contacts; default is unbound
    strand = [0]*len(wtresnames)   #   replace   strand = [0, 0, 0, 0, 0, 0, 0] --- a more general way

    restypes = []  # array of possible strand-core contact types (eg. "LYS:GLU")
    #add single contact if desired
    if onec2b == 'true':
        initialcontactstrand = int(register[0])
        initialcontactcore = int(register[2])
        strand[initialcontactstrand-2] = 1
        timeofbinding = 0        
        initialcontactstrandleft  = initialcontactstrand
        initialcontactstrandright = initialcontactstrand 
        initialcontactcoreleft    = initialcontactcore
        initialcontactcoreright   = initialcontactcore
        if float(register[2:3])/2 - int(float(register[2:3])/2) == 0:
            tag1 = "e"
        else:
            tag1 = "o"
        specificstate    =  'specific'   #  this  move to initial part, as it's depend on b2ub/onec2b or  ub2b
        orientation      =  'antiparallel'

    if ub2b == 'true':
        #from random import choice
        #register = choice(allowableregisters)
        initialcontactstrand = int(register[0])
        initialcontactcore = int(register[2])
        timeofbinding = 0

        initialcontactstrandleft  = initialcontactstrand
        initialcontactstrandright = initialcontactstrand
        initialcontactcoreleft    = initialcontactcore
        initialcontactcoreright   = initialcontactcore
        specificstate    =  'dissociate'   #  this  move to initial part, as it's depend on b2ub/onec2b or  ub2b
                                            #  change in initial state, use 
        orientation      =  'none' 
        if float(register[2:3])/2 - int(float(register[2:3])/2) == 0:
            tag1 = "e"
        else:
            tag1 = "o"

        if verbose == 'true':
            print 'register', register
           
    if b2ub == 'true':
        # 
        #        when b2ub is true, we need from input get which state we are now (e.g. odd, in-registed fully bound,.... )
        # 
        #        The following figure include all antiparallel state 
  
        #  shift 0                                   shift  +2                             shift  -2 ( not sure this is the same direction in the peper defined )
        # left  2   4   6   8  right     left      2   4   6   8  right                 left    2   4   6   8          right
        # left  8   6   4   2  right     left  8   6   4   2      right                 left        8   6   4   2      right
        # 
        #  shift 0                                   shift  +2                             shift  -2
        # left    3   5   7    right     left         3   5   7    right                left      3   5   7        right             ##  notice current this +/2 is not allowed
        # left    7   5   3    right     left     7   5   3        right                left          7   5   3    right
        #
        #           shift +2                         shift  -2                        
        # left  2   4   6   8  right     left      2   4   6   8   right         even
        # left  7   5   3      right     left          7   5   3   right         odd
        #
        #           shift +2                         shift  -2 
        # left        3   5   7    right     left  3   5   7          right     odd
        # left    8   6   4   2    right     left  8   6   4   2      right     even
        
        initialcontactstrand = int(register[0])
        initialcontactcore = int(register[2])
        specificstate    =  'specific'   #  this  move to initial part, as it's depend on b2ub/onec2b or  ub2b
        orientation      =  'antiparallel'

        if  float(register[2:3])/2 - int(float(register[2:3])/2) == 0:
            tag1 = "e"
            coreleftend  = 8
            corerightend = 2
        else:
            tag1 = "o"
            coreleftend  = 7
            corerightend = 3
        if  register in  inregisters :
            b2ubregisterstate = 'inregister'
        else:
            b2ubregisterstate = 'misregister'
        if  float(register[0])/2 - int(float(register[0])/2) == 0: 
            tag0 = "e"
            strand = [1, 1, 1, 1, 1, 1, 1] 
            strandleftend  = 2
            strandrightend = 8
        else:
            tag0 = "o"
            strandleftend  = 3
            strandrightend = 7
            strand = [0, 1, 1, 1, 1, 1, 0]
        if  tag0 == tag1  : #  common situation
             strandlength      = initialcontactstrand - strandleftend
             corelength        = coreleftend - initialcontactcore
             shift             = corelength - strandlength
             print 'shift', shift
             b2ubstate         = str(tag0 + '-' + tag1 + '-s' + str(shift))
             initialcontactstrandleft  =  max(strandleftend, strandleftend-shift)
             initialcontactstrandright =  min(strandrightend, strandrightend-shift)
             initialcontactcoreleft    =  min(coreleftend, coreleftend-shift)
             initialcontactcoreright   =  max(corerightend, corerightend-shift)
             print 'initialcontactstrandleft', initialcontactstrandleft, 'initialcontactstrandright', initialcontactstrandright, 'initialcontactcoreleft', initialcontactcoreleft, 'initialcontactcoreright', initialcontactcoreright
             tt = 0             
             while tt < (abs(shift)) :
                   if shift > 0:
                       stateleft  = 0
                       stateright = shift
                       strand[strandrightend-tt-2] = 0  # element -2 is resid
                   else:
                       strand[strandleftend+tt-2] = 0
                       stateleft  = -shift
                       stateright = 0
                   tt=tt+1
             print strand          
        if  tag0 != tag1  :  #  note, as the even -odd or odd-even has no inregister state, it will be more triky to set
             strand = [0, 0, 0, 0, 0, 0, 0]
             strandlength      = initialcontactstrand - strandleftend
             corelength        = coreleftend - initialcontactcore
             if tag0  == 'e' :
               shift             = (corelength - strandlength +1 ) *2 
             else :
               shift             = (corelength - strandlength -1 ) *2     
             b2ubstate         = str(tag0 + '-' + tag1 + '-s' + str(shift))
             initialcontactstrandleft  =  initialcontactstrand
             initialcontactstrandright =  initialcontactstrand
             initialcontactcoreleft    =  initialcontactcore
             initialcontactcoreright   =  initialcontactcore
             strand[initialcontactstrandleft-2] = 1            

             # left side
             while  initialcontactstrandleft >= strandleftend +2 and initialcontactcoreleft <= coreleftend - 2 :
                   initialcontactstrandleft = initialcontactstrandleft -2 
                   initialcontactcoreleft   = initialcontactcoreleft +2                      
                   strand[initialcontactstrandleft-1] = 1  #  middle part need be filled
                   strand[initialcontactstrandleft-2] = 1  # element -2 is resid
             while  initialcontactstrandright <= strandrightend-2 and initialcontactcoreright >= corerightend + 2 :
                   initialcontactstrandright = initialcontactstrandright + 2
                   initialcontactcoreright   = initialcontactcoreright   -2 
                   strand[initialcontactstrandright-3] = 1   #  middle part need be filled
                   strand[initialcontactstrandright-2] = 1

             print 'initialcontactstrandleft', initialcontactstrandleft, 'initialcontactstrandright', initialcontactstrandright, 'initialcontactcoreleft', initialcontactcoreleft, 'initialcontactcoreright', initialcontactcoreright
             print strand
             
        timeofbinding = 0

#############################################################################################################
#                        ******************    read transition time matrix   ***************************
#############################################################################################################

    transitionnames  = []
    transitiontimes  = {}
    transitionevents = {}
    SS_beta_proba    = {}
      

    for line in open(transitiondatafile, 'rU'):
        columns = line.split()
        if len(columns) <= 1 :
            pass 
        elif len(columns) >=2 and  columns[0][0] != '#' and columns[0] == 'antiparalleldata':            
            # add SS state for bond formation
            if int(columns[1]) > int(columns[2]) :
                key = str(columns[1] + ' ' + columns[2] + ' ' + columns[4] ) + ' antiparallel ' +  str(columns[5])  
                       #         2                   0            GLU:LYS                           SS state      form from 0?1  
            else :
                key = str(columns[1] + ' ' + columns[2] + ' ' + columns[4])  + ' antiparallel ' +  str(columns[5])
                #         0                   2                   GLU:LYS                            SS state      break into 0?1
                # e.g.  a input line  "0 2 inregister GLU:LYS 1308.15726966 15998"  
                #For bond breaks    , Ts:Tc refers to the residue pair of the breaking bond
                #For bond formations, Ts:Tc refers to the residue pair of the bond that is forming

                # note the  in / mis  is skipped !     # same as in Alex old code  L317


            if key not in transitionnames:
                transitionnames.append(key)
                #print 'Found new transition:', columns
                transitiontimes[key] = float(columns[6])
                transitionevents[key] = int(columns[7])
            else:
                #average times if there is the same transition for an inregister and misregistered strand
                transitiontimes[key] = ((transitiontimes[key] * transitionevents[key]) + (float(columns[6]) * int(columns[7])))/float(transitionevents[key] + int(columns[7]))
                transitionevents[key] = transitionevents[key] + int(columns[7])         

        if len(columns) >=2 and  columns[0][0] != '#' and columns[0] == 'paralleldata':
            # add SS state for bond formation
            if int(columns[1]) > int(columns[2]) :
                key = str(columns[1] + ' ' + columns[2] + ' ' + columns[4]) + ' parallel ' +  str(columns[5])
                       #         2                   0            GLU:LYS                      SS state  form from 0?1      
            else :
                key = str(columns[1] + ' ' + columns[2] + ' ' + columns[4]) + ' parallel ' +  str(columns[5])
                #         0                   2                   GLU:LYS                      SS state      break into 0?1
                # e.g. of full line  "0 2 inregister GLU:LYS 1308.15726966 15998"  
                #For bond breaks, Ts:Tc refers to the residue pair of the breaking bond
                #For bond formations, Ts:Tc refers to the residue pair of the bond that is forming

            if key not in transitionnames:
                transitionnames.append(key)
                #print 'Found new transition:', columns
                transitiontimes[key] = float(columns[6])
                transitionevents[key] = int(columns[7])
            else:
                #average times if there is the same transition for an inregister and misregistered strand
                transitiontimes[key] = ((transitiontimes[key] * transitionevents[key]) + (float(columns[6]) * int(columns[7])))/float(transitionevents[key] + int(columns[7]))
                transitionevents[key] = transitionevents[key] + int(columns[7])         # this line added may 4 2016, will not change result for 16-22, but important for 1-40


   #     elif columns[0] == 'nonrestraindata' : 
   #         key = str(columns[1] + ' ' + columns[2] + ' ' + columns[4] ) 
   #                # initial-state     final-state     parallel/antiparallel/none
   #                # a residue pair?
   #                                              ##  note at Aug 20  10:57  , new column has not  add to average yet.
   #         if key not in transitionnames:
   #             transitionnames.append(key)
   #             #print 'Found new transition:', columns
   #             transitiontimes[key] = float(columns[5])
   #             transitionevents[key] = int(columns[6])
   #         else:
   #             #average times if there is the same transition for an inregister and misregistered strand
   #             transitiontimes[key] = ((transitiontimes[key] * transitionevents[key]) + (float(columns[5]) * int(columns[6])))/float(transitionevents[key] + int(columns[6]))
   #             transitionevents[key] = transitionevents[key] + int(columns[6])         # this line added may 4 2016, will not change result for 16-22, but important for 1-40
                 

        elif len(columns) >=2 and columns[0] == 'SS' :         ## 0 (nonbeta) <=> 1(beta)  transition rate at free state
            key = str(columns[1] + ' ' + columns[3] + ' ' + columns[4] + ' ' + columns[5]  +  ' SS')
                  #    resid          previous SS state     final SS state     left/right
            transitiontimes[key] = float(columns[6])
        elif len(columns) >=2 and columns[0] == 'SS_free' :    ## when neighbor residue form H bond, the state/probability of being "beta" of this residue (which became the free residue can form next hbond)
            res = int(columns[1])
                  #    resid     
            SS_beta_proba[res] = float(columns[2])


        elif len(columns) >=2 and  columns[0] == 'MSM_statelist' : 
            allowablestate.append(columns[1])
 
        elif len(columns) >=2 and columns[0] == 'MSM_transitions' :
            key = str(columns[1] + ' ' + columns[2] + ' ' + columns[3] ) 
                     # initial-state     final-state     parallel/antiparallel/none
                     # if it is a residue pair, it will use residue name 

            transitionnames.append(key)
            transitiontimes[key] = float(columns[4])
            transitionevents[key] = float(columns[5])
 

    #####   print all transition event in the input file, check the state

    for key in transitionnames:

            print  key,  '  ', transitiontimes[key],  '  ', transitionevents[key] ,  '\n '



#############################################################################################################
#                        ******************    read transition time matrix  finished    ***************************
#############################################################################################################


##################################################   end of predefine  allowable events ########################################################################################################


##############################################################################################################################################################################################
   ###############################################         Matrix construction begin     ########################################################################################################

#
# the term 'neighbor' used here means neighbor state of current state 
# 
#  We keep tracking the terminal state of the peptide, 
#  Forr example, in antiparallel orientation
#  
#  For a strand, strandpointleft contact with the corepointright  
#

######################
#
#  We build five mattrix for five different situation 
#
#  1)  nonspcific:  in peptide still contact but no antiparallel  hbond;   
#  2)   singcontact, only 1 pair   antiparallel ;    
#  3)   > 1 pair  antiparallel 
#  4)  singcontact, only 1 pair   parallel     ;    
#  5)   > 1 pair  parallel


###################    note    we only define   the bond formation / dissociate transition here, the SS change were appended in each MSM step


    ktotnonboundstate      = {}
    possiblenonboundevent = {}
    possiblenonboundtime  = {}
   
    for state in  allowablestate:               ##  block 1)     
       possiblenonboundevent[state] =  {}                         #  build a  diconary for each state        
       possiblenonboundtime[state]  =  {}                         #  build an array for time
       ktotnonboundstate[state]     = 0                           # inital value                   
       
       #  now  built the matrix start fram each of the added nonspecifc-staes

       if state == 'dissociate'  :
          #  Q => B enter event, new in v2018.2
          for nextstate in  allowablestate:
            if state != nextstate :                   # identical as nextstate!= nonspecific                      
               key = state + ' ' + nextstate  + ' none'
               if key in transitionnames :                                                             #  as dissociate --> parallel are removed; they are rare but happend fast, will ruin data    
                      ktotnonboundstate[state]  = 1/transitiontimes[key] + ktotnonboundstate[state]    #  add transition rate
                      possiblenonboundevent[state][key] = nextstate                                    #  add event
                      possiblenonboundtime[state][key]  = 1/transitiontimes[key]                       #  add time
               else :
                      pass

       if state != 'dissociate'  :
          for nextstate in  allowablestate:           # if it is nonspecific, it can exit to dissociate, or any nonbound satate 
            if state != nextstate :                   # identical as nextstate!= nonspecific                      
               key = state + ' ' + nextstate  + ' none'
               if key in transitionnames :                                                             #  as dissociate --> parallel are removed; they are rare but happend fast, will ruin data    
                      ktotnonboundstate[state]  = 1/transitiontimes[key] + ktotnonboundstate[state]    #  add transition rate
                      possiblenonboundevent[state][key] = nextstate                                    #  add event
                      possiblenonboundtime[state][key]  = 1/transitiontimes[key]                       #  add time
               else :
                   #   print 'key in use', key, 'key in input'
                   #   wait = str(raw_input('This key not been sampled or error in dealing transition data, script will exit ! '))
                   #   sys.exit()
                      pass
                      #  change in 2018_v2, not average transition is observed

          ##  next two "for" loop is single bond form event

          for singcontact in allowableregisters:      # if it is nonspecific, it can go  any single bond stage   (if it not in list then sample not enough!)
               if  float(singcontact[2:3])/2 - int(float(singcontact[2:3])/2) == 0 :
                   coretag = "e"
               else :
                   coretag = "o"        
               if  tag1 == coretag :      #  we  need confirm the even/odd type are same for core in allowableregisters  and  core in current simulation 
                                          #  note tag1 is decided from  -b2ub or -onec2b input                
                      tempkey = state + ' ' + resnames[singcontact[0:1]] + ':' + resnames[singcontact[2:3]] + ' antiparallel'
                      # tempkey is defined by resdue type, more than one transition key can ref to same tempkey  (very common if a symatric peptide used)
                      key = state + ' ' + str(singcontact)  + ' antiparallel'
                      # real key is resid, has unique marker
                      if tempkey not in transitionnames:
                        #  print 'key in use', key, 'key in input', tempkey  
                        #  wait = str(raw_input('This key not been sampled or error in dealing transition data, script will exit ! '))
                        #  sys.exit()
                          pass
                         # not all transition is observed , e.g. some transition will not happen between a nb state and a bond state/certein nb state
                      else:
                          ktotnonboundstate[state]  = 1/transitiontimes[tempkey] + ktotnonboundstate[state]      #  add transition rate
                          possiblenonboundevent[state][key] = singcontact                               #  add event ; note resanem may duplicate, here we need use resid, which is unique   
                                                                                                       ##  also note different resid has same resname are treated same (have same rate) !
                          possiblenonboundtime[state][key]  = 1/transitiontimes[tempkey]                #  add time

          for singcontact in allowableregisters_parallel:      # if it is nonspecific, it can go  any single bond stage -- now parallel one   (if it not in list then sample not enough!)
               if  float(singcontact[2:3])/2 - int(float(singcontact[2:3])/2) == 0 :
                   coretag = "e"
               else :
                   coretag = "o"
               if  tag1 != coretag :      # we  need confirm the even/odd type are same for core in allowableregisters  and  core in current simulation  
                                          # this reflected in the SI 
                                          # important,  a odd core in parallel is an even core in anti parallel  (in code we treat it like this, so it more easy to find the state) 
                                          #             in manuscript, we finaaly use same core-type
                      tempkey = state + ' ' + resnames[singcontact[0:1]] + ':' + resnames[singcontact[2:3]] + ' parallel'
                      # tempkey is defined by resdue type, more than one key can ref to same tempkey  (very common if a symatric peptide used)
                      key = state + ' ' + str(singcontact)  + ' parallel'
                      # real key is resid, has unique marker
                      if tempkey not in transitionnames:
                    #      print 'key in use', key, 'key in input', tempkey
                    #      wait = str(raw_input('This key not been sampled or error in dealing transition data, script will exit ! '))
                    #      sys.exit()
                          pass
                         # not all transition is observed , e.g. some transition will not happen between a nb state and a bond state/certein nb state
                      else:
                          ktotnonboundstate[state]  = 1/transitiontimes[tempkey] + ktotnonboundstate[state]      #  add transition rate
                          possiblenonboundevent[state][key] = singcontact                               #  add event ; note resanem has duplicate, here we need use resid   
                                                                                                        #  also note different resid has same resname are treated as same !
                          possiblenonboundtime[state][key]  = 1/transitiontimes[tempkey]                #  add time

       print  'from state', state,  'pssible', possiblenonboundevent[state]

               
# only 1 hbond ;  for each etate , have 1 or 2 neighbor for hbond, >1 neighbor for break  , this is seperated for further extended a non-specific interaction part

# different from above,  in above  possiblenonboundtime[state][key]   state has no oritation and key sefin destini orientation
#                           herer  possiblesingleevent[singcontact][key]   singcontact has duplicate one in parallel and anti parallel, so we add n orentation label

    ktotsinglestate     = {}  
    possiblesingleevent = {}
    for singcontact in allowableregisters:      ##  block 2)
       #  print 'singcontact', singcontact, 'singcontact[0:1]', singcontact[0:1], 'resnames[singcontact[0:1]]', resnames[singcontact[0:1]], 'resnames[singcontact[2:3]]', resnames[singcontact[2:3]]        
       #  first define core is an even or an odd, and the end residue of core part depends on it 
       if float(singcontact[2:3])/2 - int(float(singcontact[2:3])/2) == 0:
            tag1 = "e"
            coreleftend  = 8
            corerightend = 2
       else:
            tag1 = "o"
            coreleftend  = 7
            corerightend = 3

       FCLleft  = int(singcontact[0:1])-2
       FCLright = 8-int(singcontact[0:1])
       detailedcontact = singcontact +  ' antiparallel'
             
       # initialize the value
       ktotsinglestate[detailedcontact]  = 0
       possiblesingleevent[detailedcontact]      =  {}
                                                 #  build a  diconary for each singcontact

       #  dissociate        #  not condistion limit, 1 hbond can always dissociate
       for state in  allowablestate:
           #  there are multiple NB states 

                     # it contains dissociate, but don't worry, it won't in transition list read (if it is , some thingt wrong)

           key = resnames[singcontact[0:1]] + ':' + resnames[singcontact[2:3]] + ' ' + state + ' ' + 'antiparallel'  
                                                                         #  note, here antiparallel not the finial destiniy, but type marker, other wise will duplicate with paallel  
          
 
           if key in transitionnames:
               #  change (add "if") in 2018_v2, not all transition is observed
               ktotsinglestate[detailedcontact]     = 1/transitiontimes[key]    #  add transition rate
               possiblesingleevent[detailedcontact][key] = state                #  add event



       # left side forma event (no break event, as for single H-bond, break means back to NB (see above)) 
                    
       ##  add SS tag 
       ##      for bond formation. could form from coil and beta, the rate determine by SS state of the last free residue (and of couse also dertermine by residue type and FCL)
       ##  also delete the key not used in MSM part
       #
       #   if there is free AA in the left,add situation left bond form
       if FCLleft >= 2 and int(singcontact[2:3])  < coreleftend :       
          nextFCL      = FCLleft-2
          nextstrandAA = int(singcontact[0:1])-2
          nextcoreAA   = int(singcontact[2:3])+2
          ##  form from coil
          key = str(str(FCLleft) + ' ' + str(nextFCL) + ' ' + resnames[str(nextstrandAA)] + ':' + resnames[str(nextcoreAA)]) + ' antiparallel 0' 
          #print 'singcontact', singcontact, 'key', key, 'transitiontimes[key]', transitiontimes[key]  
          if key in transitionnames  :
              ktotsinglestate[detailedcontact] = 1/transitiontimes[key] + ktotsinglestate[detailedcontact]
              possiblesingleevent[detailedcontact][key] = 'left'             #  add event  
          ##  form from beta
          key = str(str(FCLleft) + ' ' + str(nextFCL) + ' ' + resnames[str(nextstrandAA)] + ':' + resnames[str(nextcoreAA)]) + ' antiparallel 1'
          #print 'singcontact', singcontact, 'key', key, 'transitiontimes[key]', transitiontimes[key]  
          if key in transitionnames  :
              ktotsinglestate[detailedcontact] = 1/transitiontimes[key] + ktotsinglestate[detailedcontact]
              possiblesingleevent[detailedcontact][key] = 'left'             #  add event

       #  right bond form
       if FCLright >= 2 and int(singcontact[2:3]) >corerightend :  
          nextFCL      = FCLright-2
          nextstrandAA = int(singcontact[0:1])+2
          nextcoreAA   = int(singcontact[2:3])-2
          ##  form from coil
          key = str(str(FCLright) + ' ' + str(nextFCL) + ' ' + resnames[str(nextstrandAA)] + ':' + resnames[str(nextcoreAA)]) + ' antiparallel 0'
          if key in transitionnames  :  
              ktotsinglestate[detailedcontact] = 1/transitiontimes[key]  + ktotsinglestate[detailedcontact]        
              possiblesingleevent[detailedcontact][key] = 'right'
          ##  form from beta
          key = str(str(FCLright) + ' ' + str(nextFCL) + ' ' + resnames[str(nextstrandAA)] + ':' + resnames[str(nextcoreAA)]) + ' antiparallel 1'
          if key in transitionnames  :
              ktotsinglestate[detailedcontact] = 1/transitiontimes[key]  + ktotsinglestate[detailedcontact]
              possiblesingleevent[detailedcontact][key] = 'right'

          #print 'singcontact', singcontact, 'key', key, 'transitiontimes[key]', transitiontimes[key]
       #print  'possibleeventsunbind[detailedcontact]', possiblesingleevent[detailedcontact]

## > 1 hbond; state difined by "strend number"-"core number" +'side (left/right)'
    ktotsidestate = {}
    possiblesideevent = {}

    # the  possiblesideevent[detailedcontact] is set as a dictionary (key contain SS information). Later in delete necessary the forming event in coil terminal (real MSM part) we will call this

    for contact in allowableregisters:      ##  block 3)
       ##  first define core is an even or an odd, and the end residue of core part depends on it
       if float(contact[2:3])/2 - int(float(contact[2:3])/2) == 0:
            tag1 = "e"
            coreleftend  = 8
            corerightend = 2
       else:
            tag1 = "o"
            coreleftend  = 7
            corerightend = 3

       temp_strandleftend = 3    #  these two only used to prevent break of last 1 bond, which should be classed to above 1 hbond condition
       temp_strandrightend = 7   #
       FCLleft  = int(contact[0:1])-2
       FCLright = 8-int(contact[0:1])
       #print 'contact ', contact  
       #   for each contact, it generate two key, one for left and one for right
       ##  left part, start from break, notice not all allowableregisters can take here,, e. g. 8-2 will not have left break
       detailedcontact = str(contact + ' left antiparallel')
       ktotsidestate[detailedcontact] = 0
       possiblesideevent[detailedcontact]= {}   # change in SS

       #  left break
       if int(contact[0:1]) < temp_strandrightend and int(contact[2:3]) > corerightend:  #  will not consider the condition that peptide only have a single contact in the right, 
          nextFCL      = FCLleft+2                

          ##  after breake, the new free AA have beta or non-beta SS state (two type event)
          #   however, the transistion rate are same (see discussion in SI) 
         
          temp_break_time_beta = float(0)
          temp_break_time_coil = float(0) 
          temp_break_number_beta = float(0) # how many break into beta, use to calculate probability
          temp_break_number_coil = float(0) # how many break into beta, use to calculate probability
          # break into coil
          temp_key = str(str(FCLleft)) + ' ' + str(nextFCL) + ' ' +  resnames[contact[0:1]] + ':' + resnames[contact[2:3]] + ' antiparallel 0'
          if temp_key in  transitionnames  :   
              temp_break_time_coil   = transitiontimes[temp_key]           
              temp_break_number_coil = transitionevents[temp_key]
          # break into beta
          temp_key = str(str(FCLleft)) + ' ' + str(nextFCL) + ' ' +  resnames[contact[0:1]] + ':' + resnames[contact[2:3]] + ' antiparallel 1'
          if temp_key in  transitionnames  :
              temp_break_time_beta   = transitiontimes[temp_key]
              temp_break_number_beta = transitionevents[temp_key]
          # merge beta/coil event
          if temp_break_number_coil + temp_break_number_beta != 0:
   
              temp_beta_pro =  temp_break_number_beta/float(temp_break_number_coil + temp_break_number_beta)
              key = str(str(FCLleft)) + ' ' + str(nextFCL) + ' ' +  resnames[contact[0:1]] + ':' + resnames[contact[2:3]] + ' antiparalle ' + str(temp_beta_pro)
              temp_break_time = (temp_break_time_coil*temp_break_number_coil + temp_break_time_beta*temp_break_number_beta )/float(temp_break_number_coil + temp_break_number_beta)
              ktotsidestate[detailedcontact] = 1/temp_break_time
              possiblesideevent[detailedcontact][key]  =  'left_break'      # change in SS               
 
              # also add this to overall key library
              transitiontimes[key] = temp_break_time
              transitionevents[key] = temp_break_number_coil + temp_break_number_beta 
            
          else:
              print 'break event both 0, 1 for ', str(str(str(FCLleft)) + ' ' + str(nextFCL) + ' ' +  resnames[contact[0:1]] + ':' + resnames[contact[2:3]]), ' are lost' 
              quit()


       #  left  form #  if there is free AA in the left, add left bond form event
       if FCLleft >= 2 and int(contact[2:3]) < coreleftend :
          nextFCL      = FCLleft-2
          nextstrandAA = int(contact[0:1])-2
          nextcoreAA   = int(contact[2:3])+2
          ##  form from coil
          key = str(str(FCLleft) + ' ' + str(nextFCL) + ' ' + resnames[str(nextstrandAA)] + ':' + resnames[str(nextcoreAA)])  + ' antiparallel 0'
          #print 'singcontact', singcontact, 'key', key, 'transitiontimes[key]', transitiontimes[key]  
          if key in transitionnames  :
              ktotsidestate[detailedcontact] = 1/transitiontimes[key] + ktotsidestate[detailedcontact]
              possiblesideevent[detailedcontact][key]  =  'left_form'      # change in SS
          ##  form from beta
          key = str(str(FCLleft) + ' ' + str(nextFCL) + ' ' + resnames[str(nextstrandAA)] + ':' + resnames[str(nextcoreAA)])  + ' antiparallel 1'
          #print 'singcontact', singcontact, 'key', key, 'transitiontimes[key]', transitiontimes[key]  
          if key in transitionnames  :
              ktotsidestate[detailedcontact] = 1/transitiontimes[key] + ktotsidestate[detailedcontact]
              possiblesideevent[detailedcontact][key]  =  'left_form'      # change in SS


       #print 'detailedcontact', detailedcontact, 'FCLleft', FCLleft, 'possiblesideevent[detailedcontact]', possiblesideevent[detailedcontact]
       #print '                                 ktotsidestate[detailedcontact]', ktotsidestate[detailedcontact]

       ##  right part, start from break, notice not all allowableregisters can take here,, e.g. 2-8 (on right terminus) will not have right break 
       detailedcontact = contact + ' ' + 'right' + ' antiparallel'
       ktotsidestate[detailedcontact] = 0
       possiblesideevent[detailedcontact]  = {}   # change in SS
      
       # right break
       if int(contact[0:1]) > temp_strandleftend and int(contact[2:3]) < coreleftend :  #  remove condition that peptide only have a single contact in the left
          nextFCL      = FCLright+2

          ##  after breake, the new free AA have beta or non-beta SS state (two type event)
          #   however, the transistion rate are same (see discussion in SI) 

          temp_break_time_beta = float(0)
          temp_break_time_coil = float(0)
          temp_break_number_beta = float(0) # how many break into beta, use to calculate probability
          temp_break_number_coil = float(0) # how many break into beta, use to calculate probability
          # break into coil
          temp_key = str(str(FCLright)) + ' ' + str(nextFCL) + ' ' +  resnames[contact[0:1]] + ':' + resnames[contact[2:3]] + ' antiparallel 0'
          if temp_key in  transitionnames  :
              temp_break_time_coil   = transitiontimes[temp_key]
              temp_break_number_coil = transitionevents[temp_key]
          # break into beta
          temp_key = str(str(FCLright)) + ' ' + str(nextFCL) + ' ' +  resnames[contact[0:1]] + ':' + resnames[contact[2:3]] + ' antiparallel 1'
          if temp_key in  transitionnames  :
              temp_break_time_beta   = transitiontimes[temp_key]
              temp_break_number_beta = transitionevents[temp_key]
          # merge beta/coil event (take average)
          if temp_break_number_coil + temp_break_number_beta != 0:
              temp_beta_pro =  temp_break_number_beta/float(temp_break_number_coil + temp_break_number_beta)
              key = str(str(FCLright)) + ' ' + str(nextFCL) + ' ' +  resnames[contact[0:1]] + ':' + resnames[contact[2:3]] + ' antiparalle ' + str(temp_beta_pro)
              temp_break_time = (temp_break_time_coil*temp_break_number_coil + temp_break_time_beta*temp_break_number_beta )/float(temp_break_number_coil + temp_break_number_beta)
              ktotsidestate[detailedcontact] = 1/temp_break_time
              possiblesideevent[detailedcontact][key]  =  'right_break'      # change in SS               

              # also add this to overall key library
              transitiontimes[key] = temp_break_time
              transitionevents[key] = temp_break_number_coil + temp_break_number_beta
          else:
              print 'break event both 0, 1 for ', str(str(str(FCLright)) + ' ' + str(nextFCL) + ' ' +  resnames[contact[0:1]] + ':' + resnames[contact[2:3]]), ' are lost' 
              quit()

       ##  add SS tag  simulation_use_fullybind_to_unbind_v2.py and   simulation_v5.2.py , bond could form from coil and beta, the rate determine by SS state
       ##  will delete the key not used in MSM part

       #  right bond form, no mater in- or mis- register
       if FCLright >= 2 and int(contact[2:3]) > corerightend:
          nextFCL      = FCLright-2
          nextstrandAA = int(contact[0:1])+2
          nextcoreAA   = int(contact[2:3])-2
          ##  form from coil
          key = str(str(FCLright) + ' ' + str(nextFCL) + ' ' + resnames[str(nextstrandAA)] + ':' + resnames[str(nextcoreAA)]) + ' antiparallel 0'
          if key in transitionnames  :
              ktotsidestate[detailedcontact] = 1/transitiontimes[key] + ktotsidestate[detailedcontact]
              possiblesideevent[detailedcontact][key]  =  'right_form'      # change in SS
          ##  form from beta
          key = str(str(FCLright) + ' ' + str(nextFCL) + ' ' + resnames[str(nextstrandAA)] + ':' + resnames[str(nextcoreAA)]) + ' antiparallel 1'
          if key in transitionnames  :
              ktotsidestate[detailedcontact] = 1/transitiontimes[key] + ktotsidestate[detailedcontact]
              possiblesideevent[detailedcontact][key]  =  'right_form'      # change in SS
       #print 'detailedcontact', detailedcontact, 'FCLleft', FCLleft, 'possiblesideevent[detailedcontact]', possiblesideevent[detailedcontact]



######   parallel part 
# only 1 hbond ;  for each etate , have three neighbor, this is seperated for further extended a non-specific interaction part

#                           herer  possiblesingleevent[singcontact][key]   singcontact has duplicate one in parallel and anti parallel ,  and key is also important --it also has duplicate

    #ktotsinglestate_parallel     = {}
    #possiblesingleevent_parallel = {}
    for singcontact in allowableregisters_parallel:      ##  block 2)
       #print 'singcontact', singcontact, 'singcontact[0:1]', singcontact[0:1], 'resnames[singcontact[0:1]]', resnames[singcontact[0:1]], 'resnames[singcontact[2:3]]', resnames[singcontact[2:3]]        
       ##  first define core is an even or an odd, and the end residue of core part depends on it
       if float(singcontact[2:3])/2 - int(float(singcontact[2:3])/2) == 0:
            tag1 = "e"
            coreleftend  = 2
            corerightend = 8
       else:
            tag1 = "o"
            coreleftend  = 3
            corerightend = 7     #   reverse to antiparallel


       #    important , although we have seen define in MS script :  set misreg_even_strand = ( 3 5 2 4 6 4 6 8 )
       #                                                            set misreg_even_core = ( 5 7 3 5 7 3 5 7 )
       #    although core 3 5 7 are set to even in simulation, in here, we relay on the resid not the hbond form, so the  3 5 7 are odd state here
  

       FCLleft  = int(singcontact[0:1])-2
       FCLright = 8-int(singcontact[0:1])
       detailedcontact = singcontact + ' parallel'
       #  dissociate        #  not condistion limit, 1 hbond can always dissociate

       # initialize the value
       ktotsinglestate[detailedcontact]  = 0
       possiblesingleevent[detailedcontact]      =  {}
                                                 #  build a  diconary for each singcontact

       for state in  allowablestate:
           #  change in 2018_v2, now there are multiple states
           key = resnames[singcontact[0:1]] + ':' + resnames[singcontact[2:3]] + ' ' + state +  ' ' +  'parallel'
           if key in transitionnames:
               ktotsinglestate[detailedcontact]     = 1/transitiontimes[key]    #  add transition rate
               possiblesingleevent[detailedcontact][key] = state               #  add event


       # if there is free AA in the left,add situation left bond form, no mater in- or mis- register, if not, leave blank 
       if FCLleft >= 2 and int(singcontact[2:3])  >  coreleftend :   #   reverse to antiparallel
          nextFCL      = FCLleft-2
          nextstrandAA = int(singcontact[0:1])-2
          nextcoreAA   = int(singcontact[2:3])-2                     #   reverse to antiparallel
          ##  form from coil          
          key = str(str(FCLleft) + ' ' + str(nextFCL) + ' ' + resnames[str(nextstrandAA)] + ':' + resnames[str(nextcoreAA)])  +  ' parallel 0'
          #print 'singcontact', singcontact, 'key', key, 'transitiontimes[key]', transitiontimes[key]  
          if key in transitionnames  :
              ktotsinglestate[detailedcontact] = 1/transitiontimes[key] + ktotsinglestate[detailedcontact]
              possiblesingleevent[detailedcontact][key] = 'left'             #  add event  
          ##  form from beta
          key = str(str(FCLleft) + ' ' + str(nextFCL) + ' ' + resnames[str(nextstrandAA)] + ':' + resnames[str(nextcoreAA)])  +  ' parallel 1'
          #print 'singcontact', singcontact, 'key', key, 'transitiontimes[key]', transitiontimes[key]  
          if key in transitionnames  :
              ktotsinglestate[detailedcontact] = 1/transitiontimes[key] + ktotsinglestate[detailedcontact]
              possiblesingleevent[detailedcontact][key] = 'left'             #  add event  

       #  right bond form, average the mis and in-  state data
       if FCLright >= 2 and int(singcontact[2:3]) < corerightend :   #   reverse to antiparallel
          nextFCL      = FCLright-2
          nextstrandAA = int(singcontact[0:1])+2
          nextcoreAA   = int(singcontact[2:3])+2                     #   reverse to antiparallel
          ##  form from coil
          key = str(str(FCLright) + ' ' + str(nextFCL) + ' ' + resnames[str(nextstrandAA)] + ':' + resnames[str(nextcoreAA)]) +  ' parallel 0'
          if key in transitionnames  :
              ktotsinglestate[detailedcontact] = 1/transitiontimes[key]  + ktotsinglestate[detailedcontact]
              possiblesingleevent[detailedcontact][key] = 'right'
          ##  form from beta
          key = str(str(FCLright) + ' ' + str(nextFCL) + ' ' + resnames[str(nextstrandAA)] + ':' + resnames[str(nextcoreAA)]) +  ' parallel 1'
          if key in transitionnames  :
              ktotsinglestate[detailedcontact] = 1/transitiontimes[key]  + ktotsinglestate[detailedcontact]
              possiblesingleevent[detailedcontact][key] = 'right'

       #print 'detailedcontact', detailedcontact, 'key', key, 'transitiontimes[key]', transitiontimes[key], 'which side:', possiblesingleevent[detailedcontact][key]
       #print  'possibleeventsunbind[detailedcontact]', possiblesingleevent[detailedcontact]

## > 1 hbond; state difined by "strend number"-"core number" +'side (left/right)'
    # ktotsidestate = {}             ##  already defined
    # possiblesideevent = {}

    # after add SS, the  possiblesideevent[detailedcontact]  were change from [] to {}. so we can easy delete necessary the forming event in coil terminal (real MSM part)

    for contact in allowableregisters_parallel:      ##  block 3)
       ##  first find current core is an even or an odd, and the end residue of core part depends on it
       if float(contact[2:3])/2 - int(float(contact[2:3])/2) == 0:
            tag1 = "e"
            coreleftend  = 2
            corerightend = 8
       else:
            tag1 = "o"
            coreleftend  = 3
            corerightend = 7  #   reverse to antiparallel

       temp_strandleftend = 3    #  these two only used to prevent break of last 1 bond, which should be classed to above 1 hbond condition
       temp_strandrightend = 7   #
       FCLleft  = int(contact[0:1])-2
       FCLright = 8-int(contact[0:1])
       #print 'contact ', contact  
       # for each contact, it generate two key, one for left and one for right
       ##  left part, start from break, notice not all allowableregisters can take here,, e. g. 8-2 will not have left break
       detailedcontact = contact + ' '  + 'left' + ' parallel'
       ktotsidestate[detailedcontact] = 0
       possiblesideevent[detailedcontact] = {}   # change in SS  

       #  left break
       if int(contact[0:1]) < temp_strandrightend and int(contact[2:3]) < corerightend :  #  remove condition that peptide only have a single contact in the right, 
          nextFCL      = FCLleft+2                #  notice in fact only peptide have > 2 contacts will call this part, 
                                                  #  however keep this setting as we could not cll empty event in the tranisition event list (mnot matrix)

          ##  need load breake into different SS state

          temp_break_time_beta = float(0)
          temp_break_time_coil = float(0)
          temp_break_number_beta = float(0) # how many break into beta, use to calculate probability
          temp_break_number_coil = float(0) # how many break into beta, use to calculate probability
          # break into coil
          temp_key = str(str(FCLleft)) + ' ' + str(nextFCL) + ' ' +  resnames[contact[0:1]] + ':' + resnames[contact[2:3]] +  ' parallel 0'
          if temp_key in  transitionnames  :           
              temp_break_time_coil   = transitiontimes[temp_key]
              temp_break_number_coil = transitionevents[temp_key]
          # break into beta
          temp_key = str(str(FCLleft)) + ' ' + str(nextFCL) + ' ' +  resnames[contact[0:1]] + ':' + resnames[contact[2:3]] +  ' parallel 1'
          if temp_key in  transitionnames  :
              temp_break_time_beta   = transitiontimes[temp_key]
              temp_break_number_beta = transitionevents[temp_key]
          # merge beta/coil event
          if temp_break_number_coil + temp_break_number_beta != 0:
              temp_beta_pro =  temp_break_number_beta/float(temp_break_number_coil + temp_break_number_beta)
              key = str(str(FCLleft)) + ' ' + str(nextFCL) + ' ' +  resnames[contact[0:1]] + ':' + resnames[contact[2:3]] +  ' parallel '+ str(temp_beta_pro)
              temp_break_time = (temp_break_time_coil*temp_break_number_coil + temp_break_time_beta*temp_break_number_beta )/float(temp_break_number_coil + temp_break_number_beta)
              ktotsidestate[detailedcontact] = 1/temp_break_time
              possiblesideevent[detailedcontact][key]  =  'left_break'      # change in SS         
              # also add this to overall key library
              transitiontimes[key] = temp_break_time
              transitionevents[key] = temp_break_number_coil + temp_break_number_beta
          else:
              print 'break event both 0, 1 for ', str(str(str(FCLleft)) + ' ' + str(nextFCL) + ' ' +  resnames[contact[0:1]] + ':' + resnames[contact[2:3]]), ' are lost'
              quit()

       #  left formation
       ##  add SS tag  simulation_use_fullybind_to_unbind_v2.py and   simulation_v5.2.py , bond could form from coil and beta, the rate determine by SS state
       ##  will delete the key not used in MSM part

       # if there is free AA in the left,add situation left bond form, no mater in- or mis- register, if not, leave blank 
       if FCLleft >= 2 and int(contact[2:3]) > coreleftend : #   reverse to antiparallel
          nextFCL      = FCLleft-2
          nextstrandAA = int(contact[0:1])-2
          nextcoreAA   = int(contact[2:3])-2                 #   reverse to antiparallel
          ##  form from coil
          key = str(str(FCLleft) + ' ' + str(nextFCL) + ' ' + resnames[str(nextstrandAA)] + ':' + resnames[str(nextcoreAA)]) +  ' parallel 0'
          #print 'singcontact', singcontact, 'key', key, 'transitiontimes[key]', transitiontimes[key]  
          if key in transitionnames  : 
              ktotsidestate[detailedcontact] = 1/transitiontimes[key] + ktotsidestate[detailedcontact]
              possiblesideevent[detailedcontact][key] =  'left_form'      # change in SS
          ##  form from beta
          key = str(str(FCLleft) + ' ' + str(nextFCL) + ' ' + resnames[str(nextstrandAA)] + ':' + resnames[str(nextcoreAA)]) +  ' parallel 1'
          #print 'singcontact', singcontact, 'key', key, 'transitiontimes[key]', transitiontimes[key]  
          if key in transitionnames  :
              ktotsidestate[detailedcontact] = 1/transitiontimes[key] + ktotsidestate[detailedcontact]
              possiblesideevent[detailedcontact][key] =  'left_form'      # change in SS

       #print 'detailedcontact', detailedcontact, 'FCLleft', FCLleft, 'possiblesideevent[detailedcontact]', possiblesideevent[detailedcontact]

       ##  right part, start from break, notice not all allowableregisters can take here,, e.g. 2-8 will not have right break
       detailedcontact = contact + ' ' + 'right' + ' parallel'
       ktotsidestate[detailedcontact] = 0
       possiblesideevent[detailedcontact] = {}   # change in SS
 
       # right break
       if int(contact[0:1]) > temp_strandleftend and int(contact[2:3]) > coreleftend :  #  remove condition that peptide only have a single contact in the left
          nextFCL      = FCLright+2

          ##  need load breake into different SS state, however if there is no such key , it is fine as some pair won't break into beta or coil
          # in v2.5, merge break event
          temp_break_time_beta = float(0)
          temp_break_time_coil = float(0)
          temp_break_number_beta = float(0) # how many break into beta, use to calculate probability
          temp_break_number_coil = float(0) # how many break into beta, use to calculate probability
          # break into coil
          temp_key = str(str(FCLright)) + ' ' + str(nextFCL) + ' ' +  resnames[contact[0:1]] + ':' + resnames[contact[2:3]] +  ' parallel 0'
          if temp_key in  transitionnames  :
              temp_break_time_coil   = transitiontimes[temp_key]
              temp_break_number_coil = transitionevents[temp_key]
          # break into beta
          temp_key = str(str(FCLright)) + ' ' + str(nextFCL) + ' ' +  resnames[contact[0:1]] + ':' + resnames[contact[2:3]] +  ' parallel 1'
          if temp_key in  transitionnames  :
              temp_break_time_beta   = transitiontimes[temp_key]
              temp_break_number_beta = transitionevents[temp_key]
          # merge beta/coil event
          if temp_break_number_coil + temp_break_number_beta != 0:
              temp_beta_pro =  temp_break_number_beta/float(temp_break_number_coil + temp_break_number_beta)
              key = str(str(FCLright)) + ' ' + str(nextFCL) + ' ' +  resnames[contact[0:1]] + ':' + resnames[contact[2:3]] +  ' parallel ' + str(temp_beta_pro)
              temp_break_time = (temp_break_time_coil*temp_break_number_coil + temp_break_time_beta*temp_break_number_beta )/float(temp_break_number_coil + temp_break_number_beta)
              ktotsidestate[detailedcontact] = 1/temp_break_time
              possiblesideevent[detailedcontact][key]  =  'right_break'   # change in SS

              # also add this to overall key library
              transitiontimes[key] = temp_break_time
              transitionevents[key] = temp_break_number_coil + temp_break_number_beta
              
          else:
              print 'break event both 0, 1 for ', str(str(str(FCLleft)) + ' ' + str(nextFCL) + ' ' +  resnames[contact[0:1]] + ':' + resnames[contact[2:3]]), ' are lost'
              quit()

       #  right bond formation

       if FCLright >= 2 and int(contact[2:3]) < corerightend:   #   reverse to antiparallel
          nextFCL      = FCLright-2
          nextstrandAA = int(contact[0:1])+2
          nextcoreAA   = int(contact[2:3])+2                    #   reverse to antiparallel
          ##  form from coil
          key = str(str(FCLright) + ' ' + str(nextFCL) + ' ' + resnames[str(nextstrandAA)] + ':' + resnames[str(nextcoreAA)]) +  ' parallel 0'
          if key in transitionnames  :
              ktotsidestate[detailedcontact] = 1/transitiontimes[key] + ktotsidestate[detailedcontact]
              possiblesideevent[detailedcontact][key]  =  'right_form'      # change in SS
          ##  form from beta
          key = str(str(FCLright) + ' ' + str(nextFCL) + ' ' + resnames[str(nextstrandAA)] + ':' + resnames[str(nextcoreAA)]) +  ' parallel 1'
          if key in transitionnames  :
              ktotsidestate[detailedcontact] = 1/transitiontimes[key] + ktotsidestate[detailedcontact]
              possiblesideevent[detailedcontact][key]  =  'right_form'      # change in SS

       #print 'detailedcontact', detailedcontact, 'FCLleft', FCLleft, 'possiblesideevent[detailedcontact]', possiblesideevent[detailedcontact]
  
##################################################   end of predefine  allowablevents ########################################################################################################

                   
####################################################################################
# Actual MSM simulation

    # the following setup aim at tracking substates transitions, the total residence time (of sub states) and average life time
    state_residence_time = {}
    state_occur          = {}
    sub_state_residence_time = {} ##  add in  simulation_2018_v2.6.py 
    sub_state_occur          = {}
    for key in allowableregister :
        state_residence_time[key]  = float(0)
        state_occur[key]           = float(0)
        sub_state_residence_time[key]  = {}
        sub_state_occur[key]           = {}
       # for Hbond_unmber in [1,2,3,4]  :
       #      sub_state_residence_time[key][Hbond_unmber]  = float(0)   
       #      sub_state_occur[key][Hbond_unmber]           = float(0)
        for FCL_unmber in [0,1,2,3,4,5,6]  :
             sub_state_residence_time[key][FCL_unmber]  = float(0) 
             sub_state_occur[key][FCL_unmber]           = float(0) 

    ##  the non specific states recorder (name like nonspecific_b nonspecific_1 .... )
    for key in allowablestate:
        state_residence_time[key]  = float(0)
        state_occur[key]           = float(0)

    state_tracking_flag      = 'NA'
    # this flag in  begin  means it is in tracking an event  ,  and then it reset to NA after put life time to state_tracking_timestep
    # it should be set to end if start from nb state,           and then it reset to NA after reset state_residence_time   
    # 'NA' means tracking on a state/sub-state ; 'begin' means first enter a state/sub-state, need record initial time
                                   
  
    state_tracking_timestep = float(0)
    # this will be updated every time it leaves an non specific state, use to output lifge time of  nonspecific states

    if b2ub == 'true'  or  onec2b  == 'true' :
      # in b2ub, a paire will still be given and this can be used to initialize tempfullregister,  in case the peptide direct dissociatre
      tempfullregister = find_state_shift('NA ' +str(register)+ ' ' + str(orientation) )
      #  find_state_shift not read first state
      state_tracking_flag       = 'begin' 


    ##############  initilize all parameters
    sub_state_lasttime = float(0) 
    if strand == [0, 0, 0, 0, 0, 0, 0]:
        stateleft = 8
        stateright = 0

    time        = 0
    timeelapsed = 0
    lasttime    = 0             #   this is a record when last time leave nonspecific and form specific again
    lasttime    =  time         #   this use to trace how much time spend on last leave a stage, only used by ub2b and b2ub
    b2ubmarker  =  'true'       #   we also want to trace the time b2ub spend when first enter non specifc state
    b2ub_fbcounter =  0         #   we also  count how many times from fully bind to unbind will back to a fb state
    NCcounter      =  0         #   we also  count how many times  will cross nonspecific state
    NCtime         =  float(0)  #            .... and hownamy time spend on it
    NC_last_time   =  float(0)  #            .... count from everay time it break a single hygrogen bond
    eventkey = ''
    
    if verbose == 'true':
        print
        print 'Initial Configuration'
        print 'Time(ps): 0'
        printstrand(strand, restypes, initialcontactstrand)
        print 'restypes', restypes
        print 'strand', strand
        wait = str(raw_input('Press Enter to continue'))
 
    ##  origin there is a read transition matrix file part here, move above     
    loopcounter = 0
    #simulation will continue until desired event occurs (ie full unbinding) at which point the loop will break

    ##  initialize the strandpoint and corepoint  
    strandpointleft  =  initialcontactstrandleft  
    strandpointright =  initialcontactstrandright
    corepointleft    =  initialcontactcoreleft
    corepointright   =  initialcontactcoreright         
    # specificstate    =  'specific'   #  this  move to initial part, as it's depend on b2ub/onec2b or  ub2b
    lastregister     =  register + ' antiparallel'


    #  for SS
    if ub2b != 'true':
        if strandpointleft  ==  strandleftend  :
            SS_left = 999
            #  a null value for those reaching edges 
        else :
            SS_left = 0
        if strandpointright  ==  strandrightend  :
            SS_right = 999
            #  a null value for those reaching edges 
        else :
            SS_right = 0
    else:   
            SS_right = 0    
            SS_left = 0


    # add in  simulation_2018_v2.6.py 
    prev_specificstate = 'nonspecific_b'

    if pdboutput == 'true' :
        model = writepdb(strand, loopcounter, time, timeelapsed,  specificstate, orientation, strandpointleft, corepointleft)

    registerstate = 'misregister'
    if register in inregisters:
        registerstate = 'inregister'

    while True:
        #generate 2 random numbers required for gillespie algorithm
        randomeventnumber = random.random()
        randomtimenumber = random.random()
        #STEPS FOR WHEN UNBOUND PEPTIDE BINDS TO AGGREGATE

        if loopcounter != 0 and pdboutput == 'true'  :
           model = writepdb(strand, model, time, timeelapsed,  specificstate, orientation, strandpointleft, corepointleft)

        loopcounter = loopcounter+1   
        #the following steps will find all possible events given current state
        #the ktot* below will be calculated and summed as ktot
        ktotleft = 0
        ktotright = 0
        ktotunbind = 0

        possibleeventsleft = []
        possibleeventsright = []
        possibleeventsunbind = []

        if verbose == 'true':
            print 'POSSIBLE EVENTS:'

        keepchecking = 'true'
        leftevent = 'false'
        rightevent = 'false'
        unbindevent = 'false'
        sumtransitiontimes = 0

        prev_strand  = strand
        print 'new cycle, previous:', prev_strand

        if strand.count(1) == 0 and specificstate in allowablestate  :        #    and specificstate != 'dissociate':                       
           ktot        =  ktotnonboundstate[specificstate]                                                        
                          ## similar to next block of single hbond break, here state change code are build-in
           print  'nonbound block'                 
                          ## note, the for even/odd core, only contact on same core were put in L534 building region

           print  'current state', specificstate
           print  'passible events'
           print   possiblenonboundevent[specificstate]
           print  'ktot', ktot

           # unlike Hbond, we record the nonspecific passing time/counts  when leave it (in Hbond we count when enter it)
           # we add it here because current state will exsit any way, no mater it dissociate or form H bond, we just add 1 
           state_occur[specificstate] = float(state_occur[specificstate]) + 1
           prev_specificstate         = specificstate

           for key in possiblenonboundevent[specificstate]:
               #print 'key', key, 'specificstate', specificstate 
               #print 1/possiblenonboundtime[specificstate][key], '(ps), 1/possiblenonboundtime[specificstate][key] or lifetime', 1/possiblenonboundtime[specificstate][key]               
               sumtransitiontimes = sumtransitiontimes + possiblenonboundtime[specificstate][key]    # for this one, time already storded when build matrix
               #print 'sumtransitiontimes', sumtransitiontimes, 'ST/ktot', sumtransitiontimes/ktot, 'randomeventnumber', randomeventnumber
               #print 'key', key, '1/transitiontimes[key]', 1/transitiontimes[key], 'sumtransitiontimes', sumtransitiontimes,
               #print 'ktot and ST/ktot', ktot, sumtransitiontimes/ktot, 'randomeventnumber', randomeventnumber
               if   randomeventnumber <= (sumtransitiontimes/ktot) and keepchecking == 'true':
                    eventkey        = str(key) 
                    keepchecking    = 'false'
                    #  check which event happens here , note possiblenonboundevent[state][key]  is the next state, could be any nonbond state or a single hbond
                    if   possiblenonboundevent[specificstate][key] in allowablestate:          #  not form hbond,  still in nonbond state      
                         specificstate   =   possiblenonboundevent[specificstate][key]         #  the contant for possiblenonboundevent[state][key] is next state                          
                         if  specificstate  == 'dissociate' :
                               unbindevent    = 'true'
                    else :                                                 
                         #  form single hbond, the next state is specific possiblenonboundevent[state][key] is the forming single bond
                         #  reminder : key here formart  key = state + ' ' + str(singcontact)  + ' antiparallel'
                         #           and contains is   possiblenonboundevent[state][key] = singcontact    , content has no orint information    
                    
                         strandpointleft  =  int(possiblenonboundevent[specificstate][key][0:1]) 
                         corepointleft    =  int(possiblenonboundevent[specificstate][key][2:3])
                         strandpointright =  int(possiblenonboundevent[specificstate][key][0:1])
                         corepointright   =  int(possiblenonboundevent[specificstate][key][2:3])
                         lastregister     =  possiblenonboundevent[specificstate][key]
                         specificstate    =   'specific'                                        # note, this need to be changed after strandpointleft has been defined, 
                         columns          =   key.split()                                      #
                         orientation      =   columns[2]                                       # 
                         tempfullregister =   find_state_shift(key)
                         lastregister     =  str(lastregister + ' ' +  orientation)
                                             #  fix bug in oct 2 2019, previous forgot orientation


                         print  'state tracer ',  tempfullregister
                         state_occur[tempfullregister] = float(state_occur[tempfullregister]) + 1
                         #  now this register state has been recorded, when it back to non specfic, the residence time will also be recorded

                         state_tracking_flag                   = 'begin'  # this flag became true when it is in tracking an event  
                         strand[strandpointleft - 2] = 1                  # as define strandpointleft/corepointleft relay on current state not next state  
                         # reminder : the index in strand[0*7], so the index = real resid -2
                         if ub2b == 'true':
                                  lasttime         =  time                                              # this is a label of last time leave a nonspecific state, used at last part
                                                                                                        # this keep updating until the last time which end at correct fb state
                                                                                                        # however only the success binding need this
                         NCcounter     += 1                                                                #  counter for how many time cross nonspecifc
                         NCtime        += time - NC_last_time                                              #  add the time spend on this non-specific journey to the total non-specific time 
                         ##   here, a  first hbond pair formed,  assign its left/right neighbor  SS state
                         ##   note this does not take time here, the SS transition time comes as a possible event in single state
                         
                         print 'a single hbond forms !! ' , strandpointleft, corepointleft
                         print 'require reassign the strand end as orientation/odd-even may change'
                         if  float(strandpointleft)/2 - int(float(strandpointleft)/2) == 0:
                             strandleftend  = 2
                             strandrightend = 8
                         else:
                             strandleftend  = 3
                             strandrightend = 7

                         print ' new  strend end: ', strandleftend, strandrightend   
                         if strandpointleft  ==  strandleftend  :
                             SS_left = 999
                             #  a null value for those reaching edges 
                         #elif strandpointleft  ==  strandleftend + 2  :
                         #    SS_left = 1
                         #  previous we set  terminal always coil  , which not true
                         else:
                             # determine the state, use random number (0-1) and probability from simulation   
                             SSrandom      = random.random()
                             if SSrandom   > SS_beta_proba[strandpointleft-2] :
                                 SS_left = 0
                             else:
                                 SS_left = 1
                         if strandpointright  ==  strandrightend  :
                             SS_right = 999
                             #  a null value for those reaching edges 
                         #elif strandpointright ==  strandrightend - 2  :
                         #    SS_right = 1
                         #  previous we set  terminal always coil  , which not true
                         else:
                             # determine the state, use random number (0-1) and probability from simulation   
                             SSrandom      = random.random()
                             if SSrandom   > SS_beta_proba[strandpointright+2] :
                                 SS_right = 0
                             else:
                                 SS_right = 1

                    break  

        #  single hbond state are put seperately, is is junction between variers nonspecific state (a lot deaden event havn't been put here) and potential parallel event
        elif strand.count(1) == 1:
           stateindex  =  str(strandpointleft) + '-' + str(corepointleft) + ' '  + str(orientation)

           print  ' single hbond 3333333333333333333333333333333333333333333333333333333333333333333' 
           print  'strandpointleft', strandpointleft , 'strandpointright', strandpointright
           print  'SS_left', SS_left, 'SS_right ', SS_right
           lastregister     =   stateindex

           #  remider: key formart  key = resnames[singcontact[0:1]] + ':' + resnames[singcontact[2:3]] + ' ' +  'nonspecific' + ' antiparallel'
           #                                                                                                     'nonspecific' now refers multiple different states/clusters after v2018.2
           #                        possiblesingleevent[detailedcontact][key] = 'left'

           #   now, append SS time, event to the pre-defined bound event
                # key for SS transition:
                      #    resid          previous SS state     final SS state   + ' SS'

           #  first , copy origi transition here
           ktot                   =   ktotsinglestate[stateindex]                    ## note: definded ~ L377
           possible_event_current =   possiblesingleevent[stateindex].copy()           ##   important !  dictionary copy use "=" is deep copy, modify the copy, the origin one will change !

           print  'origin events'
           print   possible_event_current
           print  'origin ktot', ktot
           print  'ktotsinglestate[stateindex]', ktotsinglestate[stateindex]

           # now, add event SS transiton LEFT, remove event if could not form Hbond
           if strandpointleft  <=  strandleftend :
           #if strandpointleft  ==  strandleftend or strandpointleft  ==  strandleftend + 2  :
           #  previous we set  terminal always coil  , which not true
               pass
           elif SS_left == 1  :
               # need to add the beta to coil convert  
               key                                 =  str(str(strandpointleft-2) + ' 1 0 left SS')
               ktot                               +=  1/transitiontimes[key]
               possible_event_current[key]         =  'left_beta_to_coil' 
               
               ## delete form from coil
               if strandpointleft  >= strandleftend + 2  :
                     key_counter = 0
                     for key  in  possible_event_current :
                         columns = key.split()
                         if possible_event_current[key] == 'left' and int(columns[4]) == 0 :
                                                        # the marker of left form 
                              key_counter += 1
                              break  # should not modify dictonary while looping
                     if key_counter > 0 :
                         possible_event_current.pop(key,None)
                         ktot                               -=  1/transitiontimes[key]
 

           elif SS_left == 0  :
               # need to add the beta to coil convert  
               key                                 =  str(str(strandpointleft-2) + ' 0 1 left SS')
               ktot                               +=  1/transitiontimes[key]
               possible_event_current[key]         =  'left_coil_to_beta'
               if strandpointleft  >= strandleftend + 2  :
                     key_counter = 0
                     for key  in  possible_event_current :
                         columns = key.split()
                         if possible_event_current[key] == 'left' and int(columns[4]) == 1 :
                                                        # the marker of left form 
                              key_counter += 1
                              break  # should not modify dictonary while looping
                     if key_counter > 0 :
                         possible_event_current.pop(key,None)
                         ktot                               -=  1/transitiontimes[key]


           else :
               print 'terminal SS could not found, marker left !  Waaaaaagggggghhhhhhh '
               print 'SS_left ', SS_left, 'strandpointleft ',strandpointleft
               TIME.sleep(42424242)

           # now, add event SS transiton  RIGHT, remove event if could not form Hbond
       
           # remeber in single bond ,   strandpointleft == 'strandpointright'
           if strandpointleft >=  strandrightend :
           #if strandpointleft  ==  strandrightend or strandpointleft  ==  strandrightend - 2  :
           #  previous we set  terminal always coil  , which not true
               # single bond,  pointleft is pointright
               pass

           elif SS_right == 1  :
                               
               #  need to add the beta to coil convert  
               key                                 =  str(str(strandpointleft+2) + ' 1 0 right SS')
               ktot                               +=  1/transitiontimes[key]
               possible_event_current[key]         =  'right_beta_to_coil'
               ## delete form from coil
               if strandpointleft  <=  strandrightend - 2  :
                     key_counter = 0
                     for key  in  possible_event_current :
                         columns = key.split()
                         if possible_event_current[key] == 'right'  and  int(columns[4]) == 0 :
                                                        # the marker of left form 
                             key_counter += 1
                             break  # should not modify dictonary while looping
                     if key_counter > 0 :
                         possible_event_current.pop(key,None)
                         ktot                               -=  1/transitiontimes[key]

           elif SS_right == 0  :
               # need to add the beta to coil convert  
               key                                 =  str(str(strandpointleft+2) + ' 0 1 right SS')
               ktot                               +=  1/transitiontimes[key]
               possible_event_current[key]         =  'right_coil_to_beta'
               ## delete form from beta
               if strandpointleft  <=  strandrightend - 2  :
                     key_counter = 0                   
                     for key  in  possible_event_current :
                         columns = key.split()
                         if possible_event_current[key] == 'right' and int(columns[4]) == 1 :
                                                        # the marker of left form 
                             key_counter += 1
                             break  # should not modify dictonary while looping
                     if key_counter > 0 :
                         possible_event_current.pop(key,None)
                         ktot                               -=  1/transitiontimes[key]

           else :
               print 'terminal SS could not found, marker right ! Better Good ! '
               TIME.sleep(42424242)
           print  'new ktot', ktot
           if verbose == 'true':
               print 'debuuuuuuuuuuuuuug , is ktot correct modified?'
               debug_ktot = float(0)
               for debug_key in possible_event_current :
                    print  'Event:\t', debug_key, ', transitiontime=', transitiontimes[debug_key], ', 1/transitiontime=', 1/transitiontimes[debug_key]
                    debug_ktot +=  1/transitiontimes[debug_key]
               print 'if ktot is correct it should be ', debug_ktot

           print   ' after SS modify '
           print   possible_event_current

           ##   now back to normal MSM state
           for key in possible_event_current:
               sumtransitiontimes = sumtransitiontimes + (1/transitiontimes[key])
               if   randomeventnumber <= (sumtransitiontimes/ktot) and keepchecking == 'true':
                    eventkey        = key
                    keepchecking    = 'false'
                    #  check which event happens here
                    if   possible_event_current[key] == 'left'  :   # left bond forms  # notice only use FCL could not confirm state as the bond can form in middle
                          leftevent = 'true'                     
                          # of couse one can set condition like  columns[0] != 'bound' and int(columns[0]) == FCLright  and int(columns[1]) < int(columns[0]) 
                          stateindexleft   =  str(strandpointleft) + '-' + str(corepointleft) + ' ' + 'left'     # but if you got a partly symetric peptide this will no longer work
                         ##  form a new H bond, now determine the  SS state of the new free left residue
                          if strandpointleft  ==  strandleftend  :
                              SS_left = 999
                              #  a null value for those reaching edges 
                         # elif strandpointleft  ==  strandleftend + 2  :
                         #     SS_left = 1
                         #  previous we set  terminal always coil  , which not true
                          else:
                              # determine the state, use random number (0-1) and probability from simulation   
                              SSrandom      = random.random()
                              if SSrandom   > SS_beta_proba[strandpointleft-2] :
                                  SS_left = 0
                              else:
                                  SS_left = 1
                           
                    elif possible_event_current[key] == 'right' :  # right bond form
                          rightevent = 'true'
                          stateindexleft   =  str(strandpointleft) + '-' + str(corepointleft) + ' ' + 'right'  
                         ##  form a new H bond, now determine the  SS state of the new free right residue
                          if strandpointright  ==  strandrightend  :
                              SS_right = 999
                              #  a null value for those reaching edges 
                         # elif strandpointright ==  strandrightend - 2  :
                         #     SS_right = 1
                         #  previous we set  terminal always coil  , which not true
                          else:
                              # determine the state, use random number (0-1) and probability from simulation   
                              SSrandom      = random.random()
                              if SSrandom   > SS_beta_proba[strandpointright+2] :
                                  SS_right = 0
                              else:
                                  SS_right = 1

                  #  elif possible_event_current[key] == 'nonspecific' :    #  change after v2018.2
                    elif possible_event_current[key] in allowablestate and specificstate != 'dissociate' :
                         specificstate    =  possible_event_current[key]   
                                #  different other event whih rely on outer loop to judge, the specific to non specific change need only a single line   
                         orientation      =  'none'   
                         strand = [0, 0, 0, 0, 0, 0, 0]                         #  clear the sequence

                         state_tracking_flag = 'end'
                         #  rise the flag of a trace process is finished, output code is below the time calculation part
                         if b2ubmarker  ==  'true'  and  b2ub == 'true' :
                              b2ubmarker  =  'false'   # only trace the first time
                              lasttime         =  time          
                              print 'first time enter nonspecific'
                         NC_last_time = time           # begin tracking how many real time spend on the nonspoecific state before dissociate/form a single hbond again
 
                    elif possible_event_current[key] == 'left_beta_to_coil' :
                         SS_left = 0 
                    elif possible_event_current[key] == 'left_coil_to_beta':
                         SS_left = 1
                    elif possible_event_current[key] == 'right_beta_to_coil':
                         SS_right = 0 
                    elif possible_event_current[key] == 'right_coil_to_beta':
                         SS_right = 1

                    else : 
                         wait = str(raw_input( 'uncommon situation happened for a single contact state, neither bond forming or break, results must be wrong'))
                         sys.exit()
                    if verbose == 'true':
                         print 'Event:\t', key, ',\tSum(Ki)/Kt=', (sumtransitiontimes/ktot)
                    break
               
        else:    
           print  ' > 1 hbond '
           stateindexleft   =  str(strandpointleft) + '-' + str(corepointleft) + ' ' + 'left' + ' ' + str(orientation)
           stateindexright  =  str(strandpointright) + '-' + str(corepointright) + ' ' + 'right' + ' ' + str(orientation)
           ktot        =  float(ktotsidestate[stateindexleft]) + float(ktotsidestate[stateindexright]) 


           print 'strandpointleft ', strandpointleft, 'strandleftend  ' ,  strandleftend
           print 'strandpointright ', strandpointright, 'strandrightend', strandrightend
           print ' SS left ', SS_left, ' SS_right ', SS_right

           #   now, append SS time, event to the pre-defined bound event
           print  'stateindexleft', stateindexleft, 'stateindexright', stateindexright
           print  'possiblesideevent[stateindexleft]',  possiblesideevent[stateindexleft]
           print   'possiblesideevent[stateindexright]', possiblesideevent[stateindexright]
           possible_event_current =   possiblesideevent[stateindexleft].copy()
           possible_event_current.update(possiblesideevent[stateindexright])             # python 2  way of merge andicitnary 

           print  'origin events'
           print  possiblesideevent[stateindexleft]
           print  possiblesideevent[stateindexright]
           print  'origin ktot' , ktot
           print  'ktotsidestate[stateindexleft]', ktotsidestate[stateindexleft], 'ktotsidestate[stateindexright]', ktotsidestate[stateindexright]


           # now, add event SS transiton, remove event if could not form Hbond
           if strandpointleft  ==  strandleftend :
           #if strandpointleft  ==  strandleftend or strandpointleft  ==  strandleftend + 2  :
           #  previous we set  terminal always coil  , which not true
               pass
           elif SS_left == 1  :
               #  need to add the beta to coil convert  
               key                                 =  str(str(strandpointleft-2) + ' 1 0 left SS')
               ktot                               +=  1/transitiontimes[key]
               possible_event_current[key]         =  'left_beta_to_coil'
               ## delete form from coil
               if strandpointleft  >= strandleftend + 2  :
                     key_counter = 0
                     for key  in  possible_event_current :
                         columns = key.split()
                         if possible_event_current[key] == 'left_form'  and  int(columns[4]) == 0 :
                                                        # the marker of left form 
                             key_counter += 1
                             break  # should not modify dictonary while looping
                     if key_counter > 0 :
                         possible_event_current.pop(key,None)
                         ktot                               -=  1/transitiontimes[key]


           elif SS_left == 0  :
               print  strandpointleft, strandleftend
               # need to add the beta to coil convert  
               key                                 =  str(str(strandpointleft-2) + ' 0 1 left SS')
               ktot                               +=  1/transitiontimes[key]
               possible_event_current[key]         =  'left_coil_to_beta'
               ## delete form from beta
               if strandpointleft  >= strandleftend + 2  :
                     key_counter = 0
                     for key  in  possible_event_current :
                         columns = key.split()
                         if possible_event_current[key] == 'left_form' and  int(columns[4]) == 1 :
                                                         # the marker of left form 
                             key_counter += 1
                             break  # should not modify dictonary while looping
                     if key_counter > 0 :
                         possible_event_current.pop(key,None)
                         ktot                               -=  1/transitiontimes[key]


           else :
               print 'terminal SS could not found, marker Waaaaaagggggghhhhhhh '
               TIME.sleep(42424242)

           # now, add event SS transiton  RIGHT, remove event if could not form Hbond
           print strandpointright
           print 'examing right right SS form event '
           if strandpointright  ==  strandrightend:
           #if strandpointright  ==  strandrightend or strandpointright  ==  strandrightend - 2  :
           #  previous we set  terminal always coil  , which not true
               pass
           elif SS_right == 1  :
               #  need to add the beta to coil convert  
               key                                 =  str(str(strandpointright+2) + ' 1 0 right SS')
               ktot                               +=  1/transitiontimes[key]
               possible_event_current[key]         =  'right_beta_to_coil'
               ## delete form from coil as  we are current in beta
               if strandpointright <= strandrightend - 2  :
                     key_counter = 0
                     for key  in  possible_event_current :
                         columns = key.split()
                         if possible_event_current[key] == 'right_form' and  int(columns[4]) == 0 :
                                                        # the marker of right form 
                             key_counter += 1 
                             break  # should not modify dictonary while looping
                     if key_counter > 0 :
                         possible_event_current.pop(key,None)
                         ktot                               -=  1/transitiontimes[key]


           elif SS_right == 0  :
               # need to add the beta to coil convert  
               key                                 =  str(str(strandpointright+2) + ' 0 1 right SS')
               ktot                               +=  1/transitiontimes[key]
               possible_event_current[key]         =  'right_coil_to_beta'
               if strandpointright <= strandrightend - 2    :
                     key_counter = 0
                     for key  in  possible_event_current :
                         columns = key.split()
                         if possible_event_current[key] == 'right_form' and  int(columns[4]) == 1 :
                                                        # the marker of right form 
                             key_counter += 1
                             break  # should not modify dictonary while looping
                     if key_counter > 0 :
                         possible_event_current.pop(key,None)
                         ktot                               -=  1/transitiontimes[key]


           else :
               print 'terminal SS could not found, marker right ! Better Good ! '
               TIME.sleep(42424242)

           print  'new    ktot', ktot
           if verbose == 'true':
               print 'debuuuuuuuuuuuuuug , is ktot correct modified?'
               debug_ktot = float(0)
               for debug_key in possible_event_current :
                    print  'Event:\t', debug_key, ', transitiontime=', transitiontimes[debug_key], ', 1/transitiontime=', 1/transitiontimes[debug_key]
                    debug_ktot +=  1/transitiontimes[debug_key]
               print 'if ktot is correct it should be ', debug_ktot
           print   ' after SS modify '
           print   possible_event_current


           #find where randomeventnumber lies in rate ladder, corresponding event will happen
           if verbose == 'true':
                   print '\tSi Sf Type'
                   print  'stateindexleft ', stateindexleft
           # previous first add the left side, then right side, so event occur a which side were directly known when the event chosen 
           # now after add SS we relay on the dictonry contents
           for key in  possible_event_current: 
                sumtransitiontimes = sumtransitiontimes + (1/transitiontimes[key])
                if verbose == 'true':
                    print 'Event:\t', key, ',\t\tSum(Ki)/Kt=', (sumtransitiontimes/ktot)
                #   now choose event  ,  also update SS state if there is no bond form/break
                #   note the bond step is updated at "update strand array"  part
                #   SS invovle bond break/form also update there , as it require boundary residues
                if randomeventnumber <= (sumtransitiontimes/ktot) and keepchecking == 'true':
                    eventkey = key
                    keepchecking = 'false'
                    if possible_event_current[key]  == 'left_break'  or  possible_event_current[key]  == 'left_form' :
                        leftevent = 'true'
                    elif possible_event_current[key]  == 'right_break'  or  possible_event_current[key]  == 'right_form' :
                        rightevent = 'true'
                    elif possible_event_current[key]  == 'left_beta_to_coil'  :
                        SS_left = 0
                    elif possible_event_current[key]  == 'right_beta_to_coil'  :
                        SS_left = 0
                    elif possible_event_current[key]  == 'left_coil_to_beta'  :
                        SS_left = 1
                    elif possible_event_current[key]  == 'right_coil_to_beta'  :
                        SS_right = 1
                    else :
                        print 'terminal SS could not found, check marker > 1 hbond! for the empor ! '
                        TIME.sleep(42424242)
                    break

        #if verbose == 'true':
        #    print
        #    print 'kt =', str(ktot)
        #    print 'Stateleft:', stateleft , 'restypeleft', restypeleft, 'restypeleftadj', restypeleftadj, 'leftevent', leftevent, 'rightevent', rightevent
        #    print 'Stateright:', stateright , 'restyperight', restyperight, 'restyperightadj', restyperightadj
        #    print 'randomeventnumber', randomeventnumber
        #    print 'randomtimenumber', randomtimenumber

        #  events over, if non choose, sth wrong
	if  keepchecking   == 'true' :
            wait = str(raw_input( 'no event choosed, results must be wrong'))
            sys.exit()

        #determine how much time elapsed and update total time
        timeelapsed = (-1/ktot) * math.log(1-randomtimenumber)
        time = time + timeelapsed

        # add in  simulation_2018_v2.6.py  , a nonspficic event is occur, we not change flag, we just record
        #                                    also, nonspcific event is always single step, so  " timeelapsed " will enough
        if state_tracking_flag == 'NA':
             state_residence_time[prev_specificstate] = state_residence_time[prev_specificstate]  + float(timeelapsed)
        if state_tracking_flag == 'end':
             state_residence_time[tempfullregister] =   state_residence_time[tempfullregister]  + float(time) - state_tracking_timestep
             print 'state tracer finish a record of ', tempfullregister, ' the time spend is ', str(float(time) - state_tracking_timestep)
             state_tracking_flag  = 'NA'
             # add in  simulation_2018_v2.6.py  , when a Hbond state end, 
             # !!!!!!! nothing need to done, because the state_occur[nonspcific_x] will be added when exit that state
        if state_tracking_flag == 'begin' :
             state_tracking_timestep = float(time)
             state_tracking_flag  = 'NA'
             # add in  simulation_2018_v2.6.py  , when a Hbond state begin, need add time to the spate it just expericed
             state_residence_time[prev_specificstate] = state_residence_time[prev_specificstate]  + float(timeelapsed)

        # add in  simulation_2018_v2.6.py  , sub - state record 


        # based on current orientitaion, determine how coreid are changed along strand id
        if orientation == 'antiparallel':
            coreshift = 2
        elif orientation == 'parallel':
            coreshift = -2

        #update strand array corresponding to left side of peptide, if applicable
        leftbondfound = 'false'
        #if leftevent == 'true'  and orientation == 'antiparallel':
        if leftevent == 'true' :
        # correct at Jun 3   2016
            columns = eventkey.split()
            for i in range(0,len(strand)):
                if leftbondfound == 'false' and strand[i] == 1:
                    if int(columns[1]) > int(columns[0]): #bond breaks
                        strand[i] = 0
                        strand[i+1] = 0
                        strandpointleft = strandpointleft + 2
                        corepointleft   = corepointleft   - coreshift   #  antiparallel   -2, parallel +2                 
                    if int(columns[1]) < int(columns[0]): #bond forms
                        strand[i-1] = 1
                        strand[i-2] = 1
                        strandpointleft = strandpointleft - 2
                        corepointleft   = corepointleft   + coreshift   #  antiparallel   +2, parallel -2 
                    leftbondfound = 'true'
            # update SS state use new terminal point  
            if possible_event_current[key]  == 'left_break' :
                #print key
                # in v2.5, merge break event # now last number is probability of beta
                SSrandom      = random.random()
                if SSrandom   <= columns[4] :
                    SS_left = 0
                else:
                    SS_left = 1
                # update SS state

            elif  possible_event_current[key]  == 'left_form' :
                # update SS state use new terminal point
                if strandpointleft  ==  strandleftend  :
                    SS_left = 999
                    #  a null value for those reaching edges 
                #elif strandpointleft  ==  strandleftend +2 :
                #    SS_left = 1
                #  previous we set  terminal always coil  , which not true
                else:
                    # determine the state, use random number (0-1) and probability from simulation   
                    SSrandom      = random.random()
                    if SSrandom   > SS_beta_proba[strandpointleft-2] :
                        SS_left = 0
                    else:
                        SS_left = 1

        #update strand array to right side of peptide, if applicable
        rightbondfound = 'false'
        if rightevent == 'true':
            columns = eventkey.split()
            for i in range(0,len(strand)):
                if rightbondfound == 'false' and strand[len(strand) - 1 - i] == 1:
                    if int(columns[1]) > int(columns[0]): #bond breaks
                        strand[len(strand) - 1 - i] = 0
                        strand[len(strand) - 2 - i] = 0
                        strandpointright            = strandpointright - 2
                        corepointright              = corepointright   + coreshift   #  antiparallel   +2, parallel -2
                    if int(columns[1]) < int(columns[0]): #bond forms
                        strand[len(strand) - i] = 1
                        strand[len(strand) - i + 1] = 1
                        strandpointright            = strandpointright + 2
                        corepointright              = corepointright   - coreshift   #  antiparallel   -2, parallel +2  
                        # update SS state use new terminal point
                    rightbondfound = 'true'
            # update SS state use new terminal point
            if possible_event_current[key]  == 'right_break' :
                # in v2.5, merge break event # now last number is probability of beta
                SSrandom      = random.random()
                if SSrandom   <= columns[4] :
                    SS_right = 0
                else:
                    SS_right = 1
                # update SS state
            elif  possible_event_current[key]  == 'right_form' :
                # update SS state use new terminal point
                if strandpointright  ==  strandrightend  :
                    SS_right = 999
                    #  a null value for those reaching edges 
                #elif strandpointright  ==  strandrightend -2 :
                #    SS_right = 1
                #  previous we set  terminal always coil  , which not true
                else:
                   # determine the state, use random number (0-1) and probability from simulation   
                    SSrandom      = random.random()
                    if SSrandom   > SS_beta_proba[strandpointright+2] :
                        SS_right = 0
                    else:
                        SS_right = 1


        #update strand array corresponding to unbinding event, if applicable
        # as it's un bind, no need update the strandpoint & corepoint
        if unbindevent == 'true':
            strand = [0, 0, 0, 0, 0, 0, 0]
            stateleft = 8
            stateright = 0
            timebound = time - timeofbinding     # timeofbinding  = 0, defined in beginning

        #update lastregisterstate
        lastregisterstate = 'misregister'
        if lastregister in inregisters and orientation == 'antiparallel' :
            lastregisterstate = 'inregister'

        if verbose != 'true':
            print 'current cycle', loopcounter
            print 'EVENT CHOSEN:', eventkey, 'sumtransitiontimes', sumtransitiontimes
            print 'Time elapsed(ps) = (-1/ktot) * math.log(1-randomtimenumber)',  '           -1/ktot :',  -1/ktot
            print 'Time elapsed(ps) =', str(-1/ktot), '*', str(math.log(1-randomtimenumber))
            print 'Time elapsed(ps) =', str(timeelapsed)
            print 'strand', strand
            print 'init register', register, 'now lastregister', lastregister, 'specificstate', specificstate
            print 'strandpointleft', strandpointleft, 'corepointleft', corepointleft, 'strandpointright', strandpointright, 'corepointright', corepointright
            print 'Time(ps):', time
            print 'cureent orientation', orientation
            print  'end cycle'
            print
            print

        #  add in  simulation_2018_v2.6.py 
        #  now comare the strand and prev_strand, update sub_state_residence_time[key] and sub_state_occur[key]
        prev_FCL    = 7 - prev_strand.count(1)
        current_FCL = 7 - strand.count(1)
        if strand == prev_strand :
            if strand == [0, 0, 0, 0, 0, 0, 0] :
                 # non specific transition , pass
                 pass
            else:
                 # SS transition, add time, NOT add occur
                 sub_state_residence_time[tempfullregister][current_FCL] += timeelapsed
        else:
            if strand == [0, 0, 0, 0, 0, 0, 0] :
                 # break event ... , add time, NOT add occur
                 sub_state_residence_time[tempfullregister][prev_FCL] += timeelapsed
            elif prev_strand == [0, 0, 0, 0, 0, 0, 0] :
                 # form from non specific event  ... , add occur  NOT time (also need use [tempfullregister] not prev_specificstate )
                 sub_state_occur[tempfullregister][6] += 1
            elif prev_strand.count(1) < strand.count(1) :
                 # form  event      
                 # add time AND event 
                 # but time to prev_FCL, occur to current_FCL
                 sub_state_occur[tempfullregister][current_FCL]         += 1
                 sub_state_residence_time[tempfullregister][prev_FCL]   += timeelapsed
            elif prev_strand.count(1) > strand.count(1) :
                 # break  event add time, NOT add occur
                 sub_state_residence_time[tempfullregister][prev_FCL]   += timeelapsed
                                                                            # typo fxed at oct 2 2019 , previous is "timeelaps"         
            
        #check if strand is fully bound or fully unbound and break simulation loop if appropriate
        # note, from bond to dissociate, the initial/last time may same as the nonspecific state may be only enterd at last step
        # also note b2ub can lead empty result--- back to bound again 
        if float(register[2:3])/2 - int(float(register[2:3])/2) == 0:
            tag1 = "e"
        else:
            tag1 = "o"    

         #   this  code for  MD compare  only
#        if strand == [0, 0, 0, 0, 0, 0, 0] and b2ub == 'true':
#            f.write( b2ubstate + ' generalnonspecifc ' + b2ubregisterstate  + ' ' + str(time) + ' ' + tag1 + ' for MD compare '  + '\n')
#            break
        
        if strand == [0, 0, 0, 0, 0, 0, 0] and specificstate == 'dissociate':
            if verbose == 'true':
                print
                print 'Peptide is unbound'
                print 'Time:', int(time), 'ps'

            #  add in  simulation_2018_v2.6.py  print non-specific state/sub-state summary
            for key in allowablestate: #  add in  simulation_2018_v2.6.py
                print  key, ' NA NA NA ', ' total time ', str(state_residence_time[key]), ' occur ', str(state_occur[key]), ' state tracer total time report  '
                if state_occur[key] > 0 :
                     print 'nonspecific NA NA NA ',   str(state_residence_time[key]/state_occur[key]),  ' state tracer average time report  '
                else :
                     print 'nonspecific NA NA NA ',   str(0) , ' state tracer average time report  (not occur event)'

            print 'total resdence time                                                        sub_state report '
            print '                                          FCL                              sub_state report '  
            print 'state                            6        5       4       3       2       1      0     sub_state report  ' 
            for key in allowableregister :
                print key, '        ',
                for FCL_unmber in [6,5,4,3,2,1,0]  :
                     print sub_state_residence_time[key][FCL_unmber],
                print '     sub_state report result  '
 
      

            if b2ub == 'true':
                #  start from fully bound, it come to unbind  in two different way (in reality may only has the first one)
                f.write('fullybinds ' +  'unbinds ' + ' inregister ' + str(timebound) + ' for k-off  ' + '\n')
                f.write('fullybinds ' +  'generalnonspecific ' + ' inregister ' + str(lasttime-0) + ' count from laststate' +  '\n')
                f.write('fullybinds ' +  'generalnonspecific ' + ' inregister ' + str(b2ub_fbcounter) + ' time re-experience fb state' +  '\n')
                f.write('fullybinds ' +  'generalnonspecific ' + ' inregister ' + str(NCcounter) + ' ' + str(NCtime) + ' time experience nonspcific state ratio ' + \
                         str(NCtime/time) +  '\n')
                if pdboutput == 'true' :
                    model = writepdb(strand, model, time, timeelapsed,  specificstate, orientation, strandpointleft, corepointleft)         
                for key in allowableregister :
                    print   key, ' total time ', str(state_residence_time[key]), ' occur ', str(state_occur[key]), ' state tracer total time report  '
                    if state_occur[key] > 0 :
                         print   key,  str(state_residence_time[key]/state_occur[key]),  ' state tracer average time report  '
                    else :
                         print   key,  str(0) , ' state tracer average time report  (not occur event)'
                break
            if onec2b == 'true'  :
                print register[2:3]
                if float(register[2:3])/2 - int(float(register[2:3])/2) == 0:
                    tag1 = "e"
                else:
                    tag1 = "o"
                f.write(tag1 + '.' + register + ' unbinds ' + registerstate + ' ' + str(timebound) + ' count from initialstate' + '\n')
                if pdboutput == 'true' :
                    model = writepdb(strand, model, time, timeelapsed,  specificstate, orientation, strandpointleft, corepointleft)
                # f.write(tag1 + '.' + lastregister + ' unbinds ' + lastregisterstate + ' ' + str(lasttimebound) + ' count from laststate' + '\n')
                for key in allowableregister :
                    print   key, ' total time ', str(state_residence_time[key]), ' occur ', str(state_occur[key]), ' state tracer total time report  '
                    if state_occur[key] > 0 :
                         print   key,  str(state_residence_time[key]/state_occur[key]),  ' state tracer average time report  '
                    else :
                         print   key,  str(0) , ' state tracer average time report  (not occur event)'
                break
            if ub2b == 'true'  :
                f.write(' nonspecific ' + ' unbinds ' + registerstate + ' ' + str(timebound) + ' for k-on , fail ' + '\n')
                f.write(' nonspecific ' + ' unbinds ' + registerstate + ' ' + str(NCcounter) + ' time experience nonspcific state ratio ' + \
                         str(NCtime/time) +  '\n')
                if pdboutput == 'true' :
                    model = writepdb(strand, model, time, timeelapsed,  specificstate, orientation, strandpointleft, corepointleft)
                for key in allowableregister :
                    print   key, ' total time ', str(state_residence_time[key]), ' occur ', str(state_occur[key]), ' state tracer total time report  '
                    if state_occur[key] > 0 :
                         print   key,  str(state_residence_time[key]/state_occur[key]),  ' state tracer average time report  '
                    else :
                         print   key,  str(0) , ' state tracer average time report  (not occur event)'
                break



        if strand == [1, 1, 1, 1, 1, 1, 1]  and orientation == 'antiparallel' :
            if verbose == 'true':
                print 'Peptide is bound'
                print 'Time:', int(time), 'ps'
            if onec2b == 'true' and orientation == 'antiparallel' :
                timebound = time - timeofbinding
                # timeofbinding  = 0, defined in beginning
                # lasttimebound = time - lasttime
                f.write(tag1 + '.' + register + ' fullybinds ' + registerstate + ' ' + str(timebound) + ' count from initialstate' + '\n')
                if pdboutput == 'true' :
                     model = writepdb(strand, model, time, timeelapsed,  specificstate, orientation, strandpointleft, corepointleft)

                # f.write(tag1 + '.' + lastregister + ' fullybinds ' + lastregisterstate + ' ' + str(lasttimebound) + ' count from laststate' + '\n')
                for key in allowableregister :
                    print   key, ' total time ', str(state_residence_time[key]), ' occur ', str(state_occur[key]), ' state tracer total time report  '
                    if state_occur[key] > 0 :
                         print   key,  str(state_residence_time[key]/state_occur[key]),  ' state tracer average time report  '
                    else :
                         print   key,  str(0) , ' state tracer average time report  (not occur event)'


                #  add in  simulation_2018_v2.6.py  print non-specific state/sub-state summary
                for key in allowablestate: #  add in  simulation_2018_v2.6.py
                    print  key, ' NA NA NA ', ' total time ', str(state_residence_time[key]), ' occur ', str(state_occur[key]), ' state tracer total time report  '
                    if state_occur[key] > 0 :
                         print 'nonspecific NA NA NA ',    str(state_residence_time[key]/state_occur[key]),  ' state tracer average time report  '
                    else :
                         print 'nonspecific NA NA NA ',    str(0) , ' state tracer average time report  (not occur event)'

                print 'total resdence time                                                        sub_state report '
                print '                                          FCL                              sub_state report '
                print 'state                            6        5       4       3       2       1      0     sub_state report  '
                for key in allowableregister :
                    print key, '        ',
                    for FCL_unmber in [6,5,4,3,2,1,0]  :
                         print sub_state_residence_time[key][FCL_unmber],
                    print '     sub_state report result  '

                break
            if ub2b == 'true'  and orientation == 'antiparallel' :
                timebound = time - timeofbinding
                lasttimebound = time - lasttime
                f.write(' nonspecific ' + ' fullybinds ' + registerstate + ' ' + str(timebound) + ' for k-on , success' + '\n')
                f.write(' nonspecific ' + ' fullybinds ' + registerstate + ' ' + str(time-lasttime) + ' success , count from laststate' + '\n')
                f.write(' nonspecific ' + ' fullybinds ' + registerstate + ' ' + str(NCcounter) + ' time experience nonspcific state ratio ' + \
                         str(NCtime/time) +  '\n')
                if pdboutput == 'true' :
                     model = writepdb(strand, model, time, timeelapsed,  specificstate, orientation, strandpointleft, corepointleft)
                #  as it's ub2b, it meands stat tracker now enter an in-register antiparallel state, which not end yet                
                # note the  occur need not to be treated as it has record this one when leave nb state
                # also note the  onec2b  not debug for this yet
                if tag1 == "e" :
                     state_residence_time[allowableregister[0]] = state_residence_time[allowableregister[0]] + lasttimebound
                     #  0 correspond to  antiparallel even even 0     1   correspond to antiparallel odd odd 0 
                elif tag1 == "o" :
                     state_residence_time[allowableregister[1]] = state_residence_time[allowableregister[1]] + lasttimebound
                for key in allowableregister :
                     print   key, ' total time ', str(state_residence_time[key]), ' occur ', str(state_occur[key]), ' state tracer total time report  '
                     if state_occur[key] > 0 :
                         print   key,  str(state_residence_time[key]/state_occur[key]),  ' state tracer average time report  '
                     else :
                         print   key,  str(0) , ' state tracer average time report  (not occur event)'
                #  add in  simulation_2018_v2.6.py  print non-specific state/sub-state summary
                for key in allowablestate: #  add in  simulation_2018_v2.6.py
                    print  key, ' NA NA NA ', ' total time ', str(state_residence_time[key]), ' occur ', str(state_occur[key]), ' state tracer total time report  '
                    if state_occur[key] > 0 :
                         print 'nonspecific NA NA NA ',  str(state_residence_time[key]/state_occur[key]),  ' state tracer average time report  '
                    else :
                         print 'nonspecific NA NA NA ',  str(0) , ' state tracer average time report  (not occur event)'

                print 'total resdence time                                                        sub_state report '
                print '                                          FCL                              sub_state report '
                print 'state                            6        5       4       3       2       1      0     sub_state report  '
                for key in allowableregister :
                    print key, '        ',
                    for FCL_unmber in [6,5,4,3,2,1,0]  :
                         print sub_state_residence_time[key][FCL_unmber],
                    print '     sub_state report result  '

                break

            if b2ub == 'true'  and orientation == 'antiparallel' :
                b2ub_fbcounter +=  1


        if strand == [0, 1, 1, 1, 1, 1, 0] and int(register[2:3]) % 2 == 1   and orientation == 'antiparallel' :
            if verbose == 'true':
                print 'Peptide is bound'
                print 'Time:', int(time), 'ps'
            if onec2b == 'true' and orientation == 'antiparallel':
                timebound = time - timeofbinding
                # lasttimebound = time - lasttime
                f.write(tag1 + '.' + register + ' fullybinds ' + registerstate + ' ' + str(timebound) + ' count from initialstate' + '\n')
                if pdboutput == 'true' :
                     model = writepdb(strand, model, time, timeelapsed,  specificstate, orientation, strandpointleft, corepointleft)

                # f.write(tag1 + '.' + lastregister + ' fullybinds ' + lastregisterstate + ' ' + str(lasttimebound) + ' count from laststate' + '\n')
                for key in allowableregister :
                    print   key, ' total time ', str(state_residence_time[key]), ' occur ', str(state_occur[key]), ' state tracer total time report  '
                    if state_occur[key] > 0 :
                         print   key,  str(state_residence_time[key]/state_occur[key]),  ' state tracer average time report  '
                    else :
                         print   key,  str(0) , ' state tracer average time report  (not occur event)'

                #  add in  simulation_2018_v2.6.py  print non-specific state/sub-state summary
                for key in allowablestate: #  add in  simulation_2018_v2.6.py
                    print  key, ' NA NA NA ', ' total time ', str(state_residence_time[key]), ' occur ', str(state_occur[key]), ' state tracer total time report  '
                    if state_occur[key] > 0 :
                         print 'nonspecific NA NA NA ',  str(state_residence_time[key]/state_occur[key]),  ' state tracer average time report  '
                    else :
                         print 'nonspecific NA NA NA ',  str(0) , ' state tracer average time report  (not occur event)'

                print 'total resdence time                                                        sub_state report '
                print '                                          FCL                              sub_state report '
                print 'state                            6        5       4       3       2       1      0     sub_state report  '
                for key in allowableregister :
                    print key, '        ',
                    for FCL_unmber in [6,5,4,3,2,1,0]  :
                         print sub_state_residence_time[key][FCL_unmber],
                    print '     sub_state report result  '
                break
            if ub2b == 'true'  and orientation == 'antiparallel' :
                timebound = time - timeofbinding
                lasttimebound = time - lasttime
                f.write(' nonspecific ' + ' fullybinds ' + registerstate + ' ' + str(timebound) + ' for k-on , success' + '\n')
                f.write(' nonspecific ' + ' fullybinds ' + registerstate + ' ' + str(time-lasttime) + ' success , count from laststate' + '\n')
                f.write(' nonspecific ' + ' fullybinds ' + registerstate + ' ' + str(NCcounter) + ' time experience nonspcific state' +  '\n')
                if pdboutput == 'true' :
                     model = writepdb(strand, model, time, timeelapsed,  specificstate, orientation, strandpointleft, corepointleft)
                #  as it's ub2b, it meands stat tracker now enter an in-register antiparallel state, which not end yet                
                # note the  occur need not to be treated as it has record this one when leave nb state
                # also note the  onec2b  not debug for this yet
                if tag1 == "e" :
                     state_residence_time[allowableregister[0]] = state_residence_time[allowableregister[0]] + lasttimebound
                elif tag1 == "o" :
                     state_residence_time[allowableregister[1]] = state_residence_time[allowableregister[1]] + lasttimebound

                for key in allowableregister :
                    print   key, ' total time ', str(state_residence_time[key]), ' occur ', str(state_occur[key]), ' state tracer total time report  '
                    if state_occur[key] > 0 :
                         print   key,  str(state_residence_time[key]/state_occur[key]),  ' state tracer average time report  '
                    else :
                         print   key,  str(0) , ' state tracer average time report  (not occur event)'

                #  add in  simulation_2018_v2.6.py  print non-specific state/sub-state summary
                for key in allowablestate: #  add in  simulation_2018_v2.6.py
                    print  key, ' NA NA NA ', ' total time ', str(state_residence_time[key]), ' occur ', str(state_occur[key]), ' state tracer total time report  '
                    if state_occur[key] > 0 :
                         print 'nonspecific NA NA NA ',  str(state_residence_time[key]/state_occur[key]),  ' state tracer average time report  '
                    else :
                         print 'nonspecific NA NA NA ',  str(0) , ' state tracer average time report  (not occur event)'

                print 'total resdence time                                                        sub_state report '
                print '                                          FCL                              sub_state report '
                print 'state                            6        5       4       3       2       1      0     sub_state report  '
                for key in allowableregister :
                    print key, '        ',
                    for FCL_unmber in [6,5,4,3,2,1,0]  :
                         print sub_state_residence_time[key][FCL_unmber],
                    print '     sub_state report result '

                break

            if b2ub == 'true'  and orientation == 'antiparallel' :
                b2ub_fbcounter +=  1


    #simulation loop has exited, check whether it's because peptide is fully bound or fully unbound, record result to working file

    print
    print 'Time:', int(time), 'ps'
    print 'Done'
    

    return


###################################################################################
#Process command line arguments

    
if __name__ == '__main__':

    #this array is filled with the names of required inputs missing from the command line, it will be printed after the entire command line is processed
    missingarg = []

    print
    print 'Input parameter comments:'

    if '-dir' in sys.argv:
        dir = str(sys.argv[(sys.argv).index('-dir') + 1])
    else:
        print '-dir < > is missing (options: wt10kcal, wt1kcal, cha19, cha20, cha1920, sasa_change_phe_*, etc.)'
        missingarg.append('dir')


    if '-transitiondatafile' in sys.argv:
        transitiondatafile = str(sys.argv[(sys.argv).index('-transitiondatafile') + 1])
    else:
        print '-transitiondatafile < > is missing, default of dir/averagetransitions.dat will be used'
        transitiondatafile = str(dir) + '/averagetransitions.dat'


    pdboutput = 'false'
    if '-pdboutput' in sys.argv:
        pdboutput = 'true'
        pdbfile = open(dir + '/simulation.pdb', 'w')
    else:
        print 'no pdb output (add -pdboutput to argument line if you want this)'

    verbose = 'false'
    if '-verbose' in sys.argv:
        verbose = 'true'
    else:
        print 'no verbose output (add -verbose to argument line if you want this)'

    onec2b = 'false'
    b2ub = 'false'
    ub2b = 'false'
    register = '2-8'
    if ('-onec2b' in sys.argv and '-b2ub' in sys.argv) or ('-ub2b' in sys.argv and '-b2ub' in sys.argv) or ('-onec2b' in sys.argv and '-randomconif' in sys.argv):
        print 'Error: you must only choose one of the following arguments: -onec2b, -b2ub, or -ub2b'
        sys.exit()
    elif '-onec2b' in sys.argv:
        onec2b = 'true'
        register = str(sys.argv[(sys.argv).index('-onec2b') + 1])
        if register not in allowableregisters:
            print 'Error: invalid register given for onec2b'
            sys.exit()
    elif '-b2ub' in sys.argv:
        register = str(sys.argv[(sys.argv).index('-b2ub') + 1])
        if register not in allowableregisters:
            print 'Error: invalid register'
            sys.exit()
        b2ub = 'true'
    elif '-ub2b' in sys.argv:
        ub2b = 'true'
        register = str(sys.argv[(sys.argv).index('-ub2b') + 1]) 
        # add above register input line in v4, previous only use 2-8
        if register not in allowableregisters:
            print 'Error: invalid register given for onec2b'
            sys.exit()
    else:
        print 'strand will start unbound (add -onec2b <register> if you would like it to initially start with one contact,'
        print '-b2ub <register> if it starts fully bound, '
        print 'or -ub2b <register> if you would like it to initially start with nonspecfic )'

    if len(missingarg) != 0:
        print
        print 'Error: the following command line arguments are required:'
        for arg in missingarg:
            print arg
        print
        sys.exit()
    
    main(model, register)
    
