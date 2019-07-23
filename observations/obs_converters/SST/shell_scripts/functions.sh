#!/bin/bash
#
# This code was donated to DART by Romain Escudier, who was at Rutgers at the time.
# It is not protected by the DART copyright agreement. Thanks Romain!
#
# DART $Id$
#
#---------------------------------------------------------------------------------------------#
#                                                                                             #
# Shell functions                                                                             #
#                                                                                             #
#---------------------------------------------------------------------------------------------#

# List of functions
# Time functions
#   - get_sum_from_array(array)                   : Compute sum of an array
#   - get_date_from_cycle(cycle,startdate,dt)     : Compute date from cycle and start date
#   - get_timediff_dates(date1,date2)             : Compute difference in days between two dates
#   - is_leap_year(year)                          : Determine if year is leap year
#   - print_time_dart(date)                       : Display date in format (YYYY-MM-DD HH:MM:SS)
#   - print_time_dart_list(date)                  : Display date in format (YYYY, MM, DD, HH, MM, SS)
# Netcdf functions
#   - get_param_from_nc(file,var)                 : Get parameter value from netCDF file
#   - get_ndim_from_nc(file,dim)                  : Get dimension size from netCDF file
# Miscellaneous
#   - compute_eq_integer_result(eq)               : Compute truncated result of equation
#   - comp_files_md5sum(file1,file2)              : Compare md5sum of two files
#   - str2num(str)                                : Convert string (with padded 0) to num
#   - prog_bar(i,N)                               : Display a progress bar for loop in codes


nb_days_months=(31 28 31 30 31 30 31 31 30 31 30 31)

#---------------------------------------------------------------------------------------------#
# Compute sum of an array

get_sum_from_array() {

declare -i TOTAL=0
for vari in $*; do
  TOTAL=${TOTAL}+${vari}
done
echo ${TOTAL}

}

#---------------------------------------------------------------------------------------------#
# Compute date from cycle and start date

get_date_from_cycle() { CYCLE=$1 ; STARTDATE=$2 ; DTCYCLE=$3

   if (( ${CYCLE}<0 )); then
      get_date_from_cycle_neg ${CYCLE} ${STARTDATE} ${DTCYCLE}
   elif (( ${CYCLE}==0 )); then
      echo ${STARTDATE}
   else

      # Local variable
      declare -i ISTP_TMP=$(( ${CYCLE}*$DTCYCLE )) # Remaining number of days
      declare -i YEAR_START=$(str2num ${STARTDATE:0:4})
      declare -i YEAR_TMP=$(str2num ${YEAR_START})
      declare -i MONT_TMP=$(str2num ${STARTDATE:4:2})
      declare -i DAYS_TMP=$(str2num ${STARTDATE:6:2})

      #######################################################
      #           FIND WHICH YEAR FOR CYCLE
      #######################################################

      while ((${ISTP_TMP}>=0))
      do
         declare -i NDAYS=365
         # check if it is a leap year
         IS_LEAP=$(is_leap_year ${YEAR_TMP})
         if [ "${IS_LEAP}" = true ] ; then
            NDAYS=366
         fi

         # Different if it is the first year (may not start on Jan 1st)
         if ((${YEAR_TMP} == ${YEAR_START})) ; then
            if ((${MONT_TMP}>3)) || [ "${IS_LEAP}" = false ]; then
               NDAYS=$(get_sum_from_array ${nb_days_months[*]:${MONT_TMP}-1})-${DAYS_TMP}+1
            else
               NDAYS=$(get_sum_from_array ${nb_days_months[*]:${MONT_TMP}-1})-${DAYS_TMP}+2
            fi
         fi
         ISTP_TMP=${ISTP_TMP}-${NDAYS}
         YEAR_TMP=YEAR_TMP+1
      done

      YEAR_TMP=YEAR_TMP-1
      if ((${YEAR_TMP} != ${YEAR_START})) ; then
         MONT_TMP=1
         DAYS_TMP=1
      fi

      ISTP_TMP=$((${ISTP_TMP}+${NDAYS}+$DAYS_TMP))

      #######################################################
      #           FIND WHICH MONTH FOR CYCLE
      #######################################################

      while ((${ISTP_TMP}>0))
      do
         NDAYS=${nb_days_months[${MONT_TMP}-1]}
         if ((${MONT_TMP} == 2)) ; then
            if [ "${IS_LEAP}" = true ] ; then
               NDAYS=29
            fi
         fi
         ISTP_TMP=${ISTP_TMP}-${NDAYS}
         MONT_TMP=MONT_TMP+1
      done

      MONT_TMP=MONT_TMP-1
      declare -i DAYS_TMP=$((${ISTP_TMP}+${NDAYS}))

      MONT_DISP=$( printf "%02d" ${MONT_TMP} )
      DAYS_DISP=$( printf "%02d" ${DAYS_TMP} )
      echo "${YEAR_TMP}${MONT_DISP}${DAYS_DISP}"
   fi

}

get_date_from_cycle_neg() { CYCLE=$1 ; STARTDATE=$2 ; DTCYCLE=$3

   # Local variable
   declare -i ISTP_TMP=$(( ${CYCLE}*$DTCYCLE )) # Remaining number of days
   declare -i YEAR_START=$(str2num ${STARTDATE:0:4})
   declare -i YEAR_TMP=$(str2num ${YEAR_START})
   declare -i MONT_START=$(str2num ${STARTDATE:4:2})
   declare -i MONT_TMP=$(str2num ${MONT_START})
   declare -i DAYS_TMP=$(str2num ${STARTDATE:6:2})

   #######################################################
   #           FIND WHICH YEAR FOR CYCLE
   #######################################################

   while ((${ISTP_TMP}<=0))
   do
      declare -i NDAYS=365
      # check if it is a leap year
      IS_LEAP=$(is_leap_year ${YEAR_TMP})
      if [ "${IS_LEAP}" = true ] ; then
         NDAYS=366
      fi

      # Different if it is the first year (may not start on Jan 1st)
      if ((${YEAR_TMP} == ${YEAR_START})) ; then
         if ((${MONT_TMP}>3)) || [ "${IS_LEAP}" = false ]; then
            NDAYS=$(get_sum_from_array ${nb_days_months[*]:0:${MONT_TMP}-1})+${DAYS_TMP}
         else
            NDAYS=$(get_sum_from_array ${nb_days_months[*]:0:${MONT_TMP}-1})+${DAYS_TMP}-1
         fi
      fi
      ISTP_TMP=${ISTP_TMP}+${NDAYS}
      YEAR_TMP=YEAR_TMP-1
   done

   YEAR_TMP=YEAR_TMP+1
   if ((${YEAR_TMP} != ${YEAR_START})) ; then
      MONT_TMP=12
      DAYS_TMP=31
   fi

   ISTP_TMP=$((${ISTP_TMP}-${NDAYS}))

   #######################################################
   #           FIND WHICH MONTH FOR CYCLE
   #######################################################

   while ((${ISTP_TMP}<=0))
   do
      NDAYS=${nb_days_months[${MONT_TMP}-1]}
      if ((${MONT_TMP} == ${MONT_START})) ; then
         NDAYS=${DAYS_TMP}
      else
         if ((${MONT_TMP} == 2)) ; then
            if [ "${IS_LEAP}" = true ] ; then
               NDAYS=29
            fi
         fi
      fi
      ISTP_TMP=${ISTP_TMP}+${NDAYS}
      MONT_TMP=MONT_TMP-1
   done

   MONT_TMP=MONT_TMP+1
   declare -i DAYS_TMP=${ISTP_TMP}

   MONT_DISP=$( printf "%02d" ${MONT_TMP} )
   DAYS_DISP=$( printf "%02d" ${DAYS_TMP} )
   echo "${YEAR_TMP}${MONT_DISP}${DAYS_DISP}"

}



#---------------------------------------------------------------------------------------------#
# Compute difference in days between two dates

get_timediff_dates() { DATE1=$1 ; DATE2=$2

   declare -i YEAR1=${DATE1:0:4} MONT1=$(str2num ${DATE1:4:2}) DAYS1=$(str2num ${DATE1:6:2})
   declare -i YEAR2=${DATE2:0:4} MONT2=$(str2num ${DATE2:4:2}) DAYS2=$(str2num ${DATE2:6:2})

   declare -i NTIME=0

   #######################################################
   #           Loop on years
   #######################################################
   for YEAR_TMP in $( seq ${YEAR1} $((${YEAR2}-1)) ) ; do
      declare -i NDAYS=365
      # check if it is a leap year
      IS_LEAP=$(is_leap_year ${YEAR_TMP})

      if [ "${IS_LEAP}" = true ] ; then
         NDAYS=366
      fi
      NTIME=${NTIME}+${NDAYS}

   done

   #######################################################
   #           Remove year1 days
   #######################################################
   LEAP1=$(is_leap_year ${YEAR1})
   if [[ "${LEAP1}" = true && ((${MONT1} > 2)) ]] ; then
      NDAYS=$(get_sum_from_array ${nb_days_months[*]:0:$((${MONT1}-1))})+1
   else
      NDAYS=$(get_sum_from_array ${nb_days_months[*]:0:$((${MONT1}-1))})
   fi
   NDAYS=${NDAYS}+${DAYS1}
   NTIME=${NTIME}-NDAYS

   #######################################################
   #           Add year2 days
   #######################################################
   LEAP2=$(is_leap_year ${YEAR2})
   if [[ "${LEAP2}" = true && ((${MONT2} > 2)) ]] ; then
      NDAYS=$(get_sum_from_array ${nb_days_months[*]:0:$((${MONT2}-1))})+1
   else
      NDAYS=$(get_sum_from_array ${nb_days_months[*]:0:$((${MONT2}-1))})
   fi
   NDAYS=${NDAYS}+${DAYS2}
   NTIME=${NTIME}+NDAYS

   echo ${NTIME}
}



#---------------------------------------------------------------------------------------------#
# Determine if year is leap year

is_leap_year() { YEAR_TMP=$1

   IS_LEAP=false
   # check if it is a leap year
   declare -i B4=0
   declare -i B100=0
   declare -i B400=0
   B4=$((${YEAR_TMP}/4))
   B4=$(($B4*4))
   B100=$((${YEAR_TMP}/100))
   B100=$(($B100*100))
   B400=$((${YEAR_TMP}/400))
   B400=$(($B400*400))
   if ((${YEAR_TMP} == $B4 )) ; then
      if ((${YEAR_TMP} == $B100)) ; then
         if ((${YEAR_TMP} == $B400)) ; then
            IS_LEAP=true
         fi
      else
         IS_LEAP=true
      fi
   fi
   echo $IS_LEAP
}


#---------------------------------------------------------------------------------------------#
# Compute difference in seconds between two standard dates

get_timediff_dates_std() { DATE1=$1 ; DATE2=$2

   YEAR1=${DATE1:0:4} MONT1=${DATE1:5:2} DAYS1=${DATE1:8:2}
   YEAR2=${DATE2:0:4} MONT2=${DATE2:5:2} DAYS2=${DATE2:8:2}
   declare -i HOUR1=$(str2num ${DATE1:11:2}) MIN1=$(str2num ${DATE1:14:2}) SEC1=$(str2num ${DATE1:17:2})
   declare -i HOUR2=$(str2num ${DATE2:11:2}) MIN2=$(str2num ${DATE2:14:2}) SEC2=$(str2num ${DATE2:17:2})

   NDAYS=$(get_timediff_dates $YEAR1$MONT1$DAYS1 $YEAR2$MONT2$DAYS2)

   NSEC1AFTER=$(( (24-$HOUR1-1) * 3600 + (60-$MIN1-1) * 60 + (60-$SEC1) ))
   NSEC2BEFOR=$(( $HOUR2 * 3600 + $MIN2 * 60 + $SEC2 ))

   NSEC_TOT=$(( $NSEC1AFTER+$NSEC2BEFOR+($NDAYS-1)*86400 ))

   echo $NSEC_TOT
}


#---------------------------------------------------------------------------------------------#
# Get parameter value from netCDF file

get_param_from_nc() { FILENAME=$1 ; VARNAME=$2

   VARVALUE=$(ncdump -v ${VARNAME} ${FILENAME}  | awk -F "data:" 'RS="ceci1est8une9valeur5impossible" {print $2}' | \
              awk -F "[=,;]" '{print $2}' | xargs);
   echo ${VARVALUE}

}

#---------------------------------------------------------------------------------------------#
# Get dimension size from netCDF file

get_ndim_from_nc() { FILENAME=$1 ; DIMNAME=$2

   DIMLINE=$(ncdump -h ${FILENAME} | awk -F "variables:" 'RS="ceci1est8une9valeur5impossible" {print $1}' | \
          grep ${DIMNAME})
   NDIM=$(awk -F "[=,;]" '{print $2}' <<< $DIMLINE | xargs)
   if [ "$NDIM" == "UNLIMITED" ]; then
      NDIM=$(awk -F "[(,)]" '{print $2}' <<< $DIMLINE | awk '{print $1}' | xargs)
   fi
   echo ${NDIM}

}

#---------------------------------------------------------------------------------------------#
# Display date in dart format for namelists (YYYY-MM-DD HH:MM:SS)

print_time_dart() { MYDATE=$1

   MYYEAR=${MYDATE:0:4}
   MYMONTH=${MYDATE:4:2}
   MYDAY=${MYDATE:6:2}

   MYHR=${MYDATE:8:2}
   MYMN=${MYDATE:10:2}
   MYSC=${MYDATE:12:2}
   if [ -z $MYHR ] ; then
      DATEOUT="${MYYEAR}-${MYMONTH}-${MYDAY} 00:00:00"
   else
      if [ -z $MYMN ] ; then
         DATEOUT="${MYYEAR}-${MYMONTH}-${MYDAY} ${MYHR}:00:00"
      else
         if [ -z $MYSC ] ; then
            DATEOUT="${MYYEAR}-${MYMONTH}-${MYDAY} ${MYHR}:${MYMN}:00"
         else
            DATEOUT="${MYYEAR}-${MYMONTH}-${MYDAY} ${MYHR}:${MYMN}:${MYSC}"
         fi
      fi
   fi
   echo ${DATEOUT}

}


#---------------------------------------------------------------------------------------------#
# Display date in dart list format for namelists (YYYY, MM, DD, HH, MM, SS)

print_time_dart_list() { MYDATE=$1

   MYYEAR=${MYDATE:0:4}
   MYMONTH=${MYDATE:4:2}
   MYDAY=${MYDATE:6:2}

   MYHR=${MYDATE:8:2}
   MYMN=${MYDATE:10:2}
   MYSC=${MYDATE:12:2}
   if [ -z $MYHR ] ; then
      DATEOUT="${MYYEAR}, ${MYMONTH}, ${MYDAY}, 0, 0, 0"
   else
      if [ -z $MYMN ] ; then
         DATEOUT="${MYYEAR}, ${MYMONTH}, ${MYDAY}, ${MYHR}, 0, 0"
      else
         if [ -z $MYSC ] ; then
            DATEOUT="${MYYEAR}, ${MYMONTH}, ${MYDAY}, ${MYHR}, ${MYMN}, 0"
         else
            DATEOUT="${MYYEAR}, ${MYMONTH}, ${MYDAY}, ${MYHR}, ${MYMN}, ${MYSC}"
         fi
      fi
   fi

   echo ${DATEOUT}

}

#---------------------------------------------------------------------------------------------#
# Compute truncated result of equation

compute_eq_integer_result() { EQUATION=$1

   res=$(echo "scale=0; $1" | bc -l)
   echo ${res%.*}

}

#---------------------------------------------------------------------------------------------#
# Compare md5sum of two files

comp_files_md5sum() { FILE1=$1 ; FILE2=$2

   sum1=$( md5sum ${FILE1} ); sum1=${sum1%${FILE1}}
   sum2=$( md5sum ${FILE2} ); sum2=${sum2%${FILE2}}

   if [ "${sum1}" == "${sum2}" ] ; then
      echo "true"
   else
      echo "false"
   fi

}

#---------------------------------------------------------------------------------------------#
# Convert string (with padded 0) to num

str2num() { STR=$1

echo $(( 10#$STR ))

}



#---------------------------------------------------------------------------------------------#
# Disp progress bar

plot_bar() { I=$1; N_TOT=$2

repl() { printf "$1"'%.s' $(eval "echo {1.."$(($2))"}"); }

if (( $I % $(($N_TOT/100)) == 0 )); then
   perc=$(($I*100/$N_TOT))
   size=$(($(tput cols)*6/10))
   if (( perc > 0 )); then
      size_bar=$(($size*perc/100))
      bar=$(repl "#" $size_bar )
      size_blank=$(($size-$size_bar))
      blank=$(repl "_" $size_blank)
      disp="["$bar$blank"]\t($perc%)\t                                    \r"
    else
      blank=$(repl "_" $size)
      disp="["$blank"]\t($perc%)\t                                     \r"
   fi
   echo -ne $disp
fi

}

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

