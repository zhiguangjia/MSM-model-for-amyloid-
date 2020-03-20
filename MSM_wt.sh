#!/bin/csh

#set inreg_even_strand = ( 2 4 6 8 )
#set inreg_even_core = ( 8 6 4 2 )
#set inreg_odd_strand = ( 3 5 7 )
#set inreg_odd_core = ( 7 5 3 )
#set misreg_even_strand = ( 6 4 2 8 6 4 7 5 3 7 5 3 )
#set misreg_even_core = ( 2 4 6 4 6 8 2 4 6 4 6 8 )
#set misreg_odd_strand = ( 8 6 4 6 4 2 )
#set misreg_odd_core = ( 3 5 7 3 5 7 )
#set start = 1
#set stop = 10

set inreg_even_strand = ( 2  )
set inreg_even_core = ( 8  )
set inreg_odd_strand = ( 3  )
set inreg_odd_core = ( 7  )
set misreg_even_strand = ( 6 )
set misreg_even_core = ( 2 )
set misreg_odd_strand = ( 8 )
set misreg_odd_core = ( 3 )
#  change in v4, as the toff in misregister should not be included in toff


set register = in
set core = even

set mutant = wt
set trajnum = 1

set specialcase = wt

echo "mutant $mutant"
echo "specialcase $specialcase"


rm  $specialcase/toff.dat
rm  $specialcase/ton_success.dat
rm  $specialcase/ton_fail.dat
rm  $specialcase/toff_to_first.dat
rm  $specialcase/ton_success_from-last.dat
rm  $specialcase/toff_to_fb_counter.dat
rm   $specialcase/ton_cross_nonspecifc_counter.dat
rm   $specialcase/toff_cross_nonspecifc_counter.dat
rm   $specialcase/state_*
rm   $specialcase/2019_oct2_nb_state_time_debug_ub2b
rm   $specialcase/2019_oct2_nb_state_time_debug_b2ub


set finaldatapath = . 

cp  distribution_MSM.py    /scrach/zgjia/
cp  sub_state_report_header /scrach/zgjia/
cp  average_substate.py   /scrach/zgjia/
cp  sub_state_list          /scrach/zgjia/

repeat:

set i = 1

set sim = 1
set numjobs = 0

if ( $core == "even" && $register == "in" ) while ( $i <= $#inreg_even_strand )
if ( $core == "odd" && $register == "in" ) while ( $i <= $#inreg_odd_strand )
if ( $core == "even" && $register == "mis" ) while ( $i <= $#misreg_even_strand )
if ( $core == "odd" && $register == "mis" ) while ( $i <= $#misreg_odd_strand )

#if (`echo $protein[$i] | awk '{print substr( $0, 0, 1) }'` == "e") set core = even
#if (`echo $protein[$i] | awk '{print substr( $0, 0, 1) }'` == "o") set core = odd

    if ($core == "even") then
	set tag1 = e
	if ($register == "in") then
	    set tag2 = $inreg_even_strand[$i]
	    set tag3 = $inreg_even_core[$i]
	endif
	if ($register == "mis") then
	    set tag2 = $misreg_even_strand[$i]
	    set tag3 = $misreg_even_core[$i]
	endif
    endif
    if ($core == "odd") then
	set tag1 = o
	if ($register == "in") then
	    set tag2 = $inreg_odd_strand[$i]
	    set tag3 = $inreg_odd_core[$i]
	endif
	if ($register == "mis") then
	    set tag2 = $misreg_odd_strand[$i]
	    set tag3 = $misreg_odd_core[$i]
	endif
    endif

    
#echo $tag1.$tag2-$tag3


echo "$tag2-$tag3"



set k = 0
while ($k < 20)

echo $k

rm $specialcase/timedata
python MSM.py -b2ub $tag2-$tag3 -dir $specialcase -transitiondatafile $specialcase/averagetransitions.dat  >  temp
#cat $specialcase/timedata | grep 'for MD compare' >> $specialcase/MD_compare.dat

cat $specialcase/timedata  | grep 'for k-off'     >> $specialcase/toff.dat
cat $specialcase/timedata  | grep 'count from laststate'          >> $specialcase/toff_to_first.dat
cat  $specialcase/timedata  | grep 'time re-experience fb state'  >> $specialcase/toff_to_fb_counter.dat
cat  $specialcase/timedata  | grep ' time experience nonspcific state'  >> $specialcase/toff_cross_nonspecifc_counter.dat

rm $specialcase/timedata

./MSM_traj_to_time.py temp >> $specialcase/2019_oct2_nb_state_time_debug_b2ub 

echo b2ub done

    if ($k == 0 ) then
        cat temp | grep " state tracer total time report"  > tt
        awk '{ print " " $1 " " $2 " " $3 " " $4 " "  }' tt > State_list
        cp  State_list $specialcase/state_b2ub_totaltime_$core
        cp  State_list $specialcase/state_b2ub_averagetime_$core   
        cp  State_list $specialcase/state_b2ub_occur_$core
        cp  State_list $specialcase/state_ub2b_totaltime_$core
        cp  State_list $specialcase/state_ub2b_averagetime_$core
        cp  State_list $specialcase/state_ub2b_occur_$core

        # sub_state analysis does not need consider non-spcific   #  add in  MSM.py
        # so a previous MANUALLY generated list works fine  
        cp sub_state_report_header $specialcase/sub_state_b2ub_totaltime_$core
        cp sub_state_report_header $specialcase/sub_state_ub2b_totaltime_$core      
    endif

    # state  record 
    cat temp | grep " state tracer total time report"  > tt
    awk '{ print " " $7 }' tt >  ttt
    cp   $specialcase/state_b2ub_totaltime_$core  tt
    paste -d " " tt ttt >  $specialcase/state_b2ub_totaltime_$core

    cat temp | grep " state tracer total time report"  > tt
    awk '{ print " " $9 }' tt >  ttt
    cp   $specialcase/state_b2ub_occur_$core  tt
    paste -d " " tt ttt >  $specialcase/state_b2ub_occur_$core


    cat temp | grep " state tracer average time report"  > tt 
    awk '{ print " " $5 }' tt > ttt
    cp   $specialcase/state_b2ub_averagetime_$core tt
    paste -d " " tt ttt >      $specialcase/state_b2ub_averagetime_$core

    #  sub state record :
    #  note , treatment logic differenct from  state  record  , in state  record  , one state has one number (for one MSM sim)
    #                                                           in sub-state, one state  has 6 number (each for each FCL) (for one MSM sim)
    cat temp | grep "  sub_state report result"  >> $specialcase/sub_state_b2ub_totaltime_$core
 
   
python MSM.py -ub2b $tag2-$tag3 -dir $specialcase -transitiondatafile $specialcase/averagetransitions.dat   >  temp
#perform sim of peptide with 1 contact to find twait
#tail -n -1 $specialcase/timedata >>! $specialcase/twait.dat

cat $specialcase/timedata | grep 'for k-on , fail'     >> $specialcase/ton_fail.dat
cat $specialcase/timedata | grep 'for k-on , success'  >> $specialcase/ton_success.dat

cat $specialcase/timedata | grep 'success , count from laststate'     >> $specialcase/ton_success_from-last.dat
cat  $specialcase/timedata  | grep ' time experience nonspcific state'  >> $specialcase/ton_cross_nonspecifc_counter.dat

./MSM_traj_to_time.py temp >> $specialcase/2019_oct2_nb_state_time_debug_ub2b



    cat temp | grep " state tracer total time report"  > tt
    awk '{ print " " $7 }' tt >  ttt
    cp   $specialcase/state_ub2b_totaltime_$core  tt
    paste -d " " tt ttt >  $specialcase/state_ub2b_totaltime_$core


    cat temp | grep " state tracer total time report"  > tt
    awk '{ print " " $9 }' tt >  ttt
    cp   $specialcase/state_ub2b_occur_$core  tt
    paste -d " " tt ttt >  $specialcase/state_ub2b_occur_$core


    cat temp | grep " state tracer average time report"  > tt
    awk '{ print " " $5 }' tt > ttt
    cp   $specialcase/state_ub2b_averagetime_$core tt
    paste -d " " tt ttt >      $specialcase/state_ub2b_averagetime_$core

    #  sub state record :
    #  note , treatment logic differenct from  state  record  , in state  record  , one state has one number (for one MSM sim)
    #                                                           in sub-state, one state  has 6 number (each for each FCL) (for one MSM sim)
    cat temp | grep "  sub_state report result"  >> $specialcase/sub_state_ub2b_totaltime_$core


echo ub2b done

@ k++
#echo "Performed $k simulations for register $tag1.$tag2-$tag3"

end
# end k loop   



@ i++
end
# end i loop

if ( $core == "even") then
set core = odd
goto repeat
endif



#  v0 version read  $specialcase/state_ub2b_averagetime_$core  $specialcase/state_ub2b_totaltime_$core    in  script_history_version
./distribution_MSM.py  $specialcase/state_b2ub_occur_odd  $specialcase/state_b2ub_totaltime_odd   $specialcase/state_summary_b2ub_odd      $k   
./distribution_MSM.py  $specialcase/state_b2ub_occur_even $specialcase/state_b2ub_totaltime_even  $specialcase/state_summary_b2ub_even     $k
./distribution_MSM.py  $specialcase/state_ub2b_occur_odd  $specialcase/state_ub2b_totaltime_odd   $specialcase/state_summary_ub2b_odd      $k
./distribution_MSM.py  $specialcase/state_ub2b_occur_even $specialcase/state_ub2b_totaltime_even  $specialcase/state_summary_ub2b_even     $k


#  sub state record :
#  note , treatment logic differenct from  state  record  , in state  record  , one state has one number (for one MSM sim)
#                                                           in sub-state, one state  has 6 number (each for each FCL) (for one MSM sim)



##   get total time for sub state, classfied by FCL
./average_substate.py $specialcase/sub_state_b2ub_totaltime_odd ratio total
cat sub_state_report_header ratio > $specialcase/sub_summary_residence_time_ratio_b2ub_odd
cat sub_state_report_header total > $specialcase/sub_summary_residence_time_b2ub_odd
./average_substate.py $specialcase/sub_state_b2ub_totaltime_even ratio total
cat sub_state_report_header ratio > $specialcase/sub_summary_residence_time_ratio_b2ub_even
cat sub_state_report_header total > $specialcase/sub_summary_residence_time_b2ub_even
./average_substate.py $specialcase/sub_state_ub2b_totaltime_odd ratio total
cat sub_state_report_header ratio > $specialcase/sub_summary_residence_time_ratio_ub2b_odd
cat sub_state_report_header total > $specialcase/sub_summary_residence_time_ub2b_odd
./average_substate.py $specialcase/sub_state_ub2b_totaltime_even ratio total
cat sub_state_report_header ratio > $specialcase/sub_summary_residence_time_ratio_ub2b_even
cat sub_state_report_header total > $specialcase/sub_summary_residence_time_ub2b_even

# attach the time of nonspecific to it ..
cat  $specialcase/state_summary_b2ub_odd  | grep  nonspecific > tt
awk '{ print  $1 " " $2 " " $3 " " $4 " " $7  "    0.00    0.00    0.00    0.00    0.00    0.00    0.00" }'  tt  >> $specialcase/sub_summary_residence_time_b2ub_odd
#                #   state               residence time      FCL time (0 for nonspcific)                                        
cat  $specialcase/state_summary_b2ub_even | grep  nonspecific > tt
awk '{ print  $1 " " $2 " " $3 " " $4 " " $7  "    0.00    0.00    0.00    0.00    0.00    0.00    0.00" }'  tt  >> $specialcase/sub_summary_residence_time_b2ub_even
cat  $specialcase/state_summary_ub2b_odd  | grep  nonspecific > tt
awk '{ print  $1 " " $2 " " $3 " " $4 " " $7  "    0.00    0.00    0.00    0.00    0.00    0.00    0.00" }'  tt  >> $specialcase/sub_summary_residence_time_ub2b_odd
cat  $specialcase/state_summary_ub2b_even | grep  nonspecific > tt
awk '{ print  $1 " " $2 " " $3 " " $4 " " $7  "    0.00    0.00    0.00    0.00    0.00    0.00    0.00" }'  tt  >> $specialcase/sub_summary_residence_time_ub2b_even


#rm -r  /scrach/zgjia//$specialcase




exit
#  change in v4, as the toff in misregister should not be included in toff



if ( $register == "in") then
set register = mis
set core = even
goto repeat
endif

exit



if ($core == "odd" && $register == "mis") then
set register = in
set core = even
endif

#if ($mutant == wildtype) then
#set mutant = cha19
#goto repeat
#endif

#if ($mutant == cha19) then
#set mutant = cha20
#goto repeat
#endif

#if ($mutant == cha20) then
#set mutant = cha1920
#goto repeat
#endif
