#!/bin/csh

unalias ls
set prog = $0
set prog = $prog:t


set tmp_flsp = tmp.$$
set pwd = `pwd`

if( $#argv > 0 ) then
#---------------------------------------
#  process arguments
#---------------------------------------
  set args1 = ($argv)
  set args = `echo $args1 | tr '[a-z]' '[A-Z]'`
#---------------------------------------
#  iterate over arguments
#---------------------------------------
  foreach argument ($args)
    foreach keyword (ENS MOZBC_INP TYPE JULIAN_DAYS TOD)
      set sarg = $keyword
      set sarg = "$sarg"=
      echo $argument | grep "$sarg" >& /dev/null
      if( ! $status ) then
        if( $keyword == MOZBC_INP ) then
          set mozbc_inp = `echo $args1[1] | cut -d= -f2`
        else if( $keyword == TYPE ) then
          set bnd_typ = `echo $args1[1] | cut -d= -f2`
          set bnd_typ = `echo $bnd_typ | tr '[a-z]' '[A-Z]'`
          if( $bnd_typ != IC && $bnd_typ != BC && $bnd_typ != BOTH ) then
            echo "type keyword $bnd_typ not in {ic,bc,both}; $prog terminating"
            exit -1
          endif
          echo " "
          if( $bnd_typ != BOTH ) then
            echo "Bndy type = $bnd_typ"
          else
            echo "Bndy type = IC and BC"
          endif
        else if( $keyword == ENS ) then
          set wrk = `echo $args1[1] | cut -d= -f2`
          if( "$wrk" != "" ) then
            set ens_id = `echo $wrk | sed 's/,/ /g'`
          endif
        else if( $keyword == TOD ) then
          set wrk = `echo $args1[1] | cut -d= -f2`
          if( "$wrk" != "" ) then
            set tod = `echo $wrk | sed 's/,/ /g'`
          endif
        else if( $keyword == JULIAN_DAYS ) then
          set wrk = `echo $args1[1] | cut -d= -f2`
          if( "$wrk" != "" ) then
            set jd = `echo $wrk | sed 's/,/ /g'`
          endif
        endif
        break
      else
        continue
      endif
    end
    shift args1
  end
endif

#---------------------------------------
#  check mozbc for existence
#---------------------------------------
if( ! -e mozbc ) then
  echo "Can not find file mozbc; $prog terminating"
  exit -1
endif
#---------------------------------------
#  check bnd type for assignment
#---------------------------------------
if( ! $?bnd_typ ) then
  echo "type keyword not set, must = {ic,bc,both}; $prog terminating"
  exit -1
endif
#---------------------------------------
#  if not assigned set mozbc_inp default
#---------------------------------------
if( ! $?mozbc_inp ) then
  set mozbc_inp = mozbc.inp
endif
#---------------------------------------
#  check mozbc input namelist file for existence
#---------------------------------------
if( ! -e $mozbc_inp ) then
  echo "Can not find file $mozbc_inp; $prog terminating"
  exit -1
endif
#---------------------------------------
#  check ens_id for valid number
#---------------------------------------
if( $?ens_id ) then
  foreach ndx ($ens_id)
    (echo $ndx | egrep "[^0-9]") >& /dev/null
    if( ! $status ) then
      echo "Ensemble entry $ndx is not a number"
      exit -1
    endif
  end
else
  set ens_id = all
endif
#---------------------------------------
#  check tod for valid number
#---------------------------------------
if( $?tod ) then
  foreach ndx ($tod)
    (echo $ndx | egrep "[^0-9]") >& /dev/null
    if( ! $status ) then
      echo "time of day entry $ndx is not a number"
      exit -1
    endif
    if( $ndx < 0 || $ndx > 86400 ) then
      echo "time of day entry $ndx is not in set [0,86400]"
      exit -1
    endif
  end
else
  set tod = all
endif
#---------------------------------------
#  check julian days for valid number
#---------------------------------------
if( $?jd ) then
  foreach ndx ($jd)
    (echo $ndx | egrep "[^0-9]") >& /dev/null
    if( ! $status ) then
      echo "julian day entry $ndx is not a number"
      exit -1
    endif
  end
else
  set jd = all
endif

if( $bnd_typ == IC || $bnd_typ == BOTH ) then
  set filestem = wrfinput_d01
else
  set filestem = wrfbdy
endif

#---------------------------------------
#  setup all files list
#---------------------------------------
set files = ""
foreach ensno ($ens_id)
  set file_list = ""
  if( $ensno == all ) then
    set ens_str = ".e*"
  else
    set ens_str = ".e${ensno}"
  endif
  foreach td ($tod)
    if( $td == all ) then
      set td_str = "_*"
    else
      set td_str = _$td
    endif
    foreach julday ($jd)
      if( $julday == all ) then
        set jd_str = "_*"
      else
        set jd_str = _$jd
      endif
      set add_list = "${filestem}${jd_str}${td_str}${ens_str}"
      ls -1 $add_list >& tmp.$$
      set line1 = `head -1 tmp.$$`
      rm -f tmp.$$
      if( $ensno == all || $td == all || $julday == all ) then
        if( "$line1" == "No match." ) then
          echo "Can not find any $add_list files; $prog terminating"
          exit -1
        endif
      else
        (echo $line1 | egrep "not exist") >& /dev/null
        if( ! $status ) then
          echo "Can not find $add_list file; $prog terminating"
          exit -1
        endif
      endif
      set file_list = "$file_list ${add_list}"
    end
  end
  set new_files = `ls -1 $file_list`
  set files = ($files $new_files)
end

set file_cnt = $#files
if( $bnd_typ == BOTH ) then
  @ file_cnt *= 2
endif

#echo "file count = $file_cnt"
#echo "file list = $files"

#---------------------------------------
#  check files
#---------------------------------------
foreach file ($files)
  if( $bnd_typ == BC ) then
    ls -1 wrfinput_d01_* >& tmp.$$
    set line1 = `head -1 tmp.$$`
    rm -f tmp.$$
    if( "$line1" == "No match." ) then
      echo "Can not find any wrfinput_d01_* files; $prog terminating"
      exit -1
    endif
  else if( $bnd_typ == BOTH ) then
    set date  = `echo $file | cut -d_ -f3`
    set secs  = `echo $file | cut -d_ -f4`
    set ensno = `echo $file | cut -d_ -f5`
    set bdy_file = "wrfbdy_d01.${date}_${secs}.e${ensno}"
    if( ! -e $bdy_file ) then
      echo "File $bdy_file not found; $prog terminating"
      exit -1
    endif
  endif
end

if( -e wrfinput_d01 ) then
  rm -f wrfinput_d01
endif

if( $bnd_typ == BC ) then
  set inp_file = `ls wrfinput_d01_*`
  echo $inp_file
  ln -s $inp_file[1] wrfinput_d01
  if( $status ) then
    echo "Failed to link $inp_file to wrfinput_d01"
    exit -1
  endif
endif

if( -e wrfbdy_d01 ) then
  rm -f wrfbdy_d01
endif
#---------------------------------------
#  loop over files
#---------------------------------------

foreach file ($files)
  if( $bnd_typ == IC || $bnd_typ == BOTH ) then
#---------------------------------------
#  IC file or BOTH
#---------------------------------------
    ln -s $file wrfinput_d01 >& /dev/null
    if( $status ) then
      echo "Failed to link $file to wrfinput_d01"
      exit -1
    endif
    if( $bnd_typ == BOTH ) then
      set date  = `echo $file | cut -d_ -f3`
      set secs  = `echo $file | cut -d_ -f4`
      set ensno = `echo $file | cut -d_ -f5`
      set bdy_file = "wrfbdy_d01.${date}_${secs}.e${ensno}"
      ln -s $bdy_file wrfbdy_d01 >& /dev/null
      if( $status ) then
        echo "Failed to link $bdy_file to wrfbdy_d01"
        exit -1
      endif
    endif
  else
#---------------------------------------
#  BC file
#---------------------------------------
    ln -s $file wrfbdy_d01
    if( $status ) then
      echo "Failed to link $file to wrfdby_d01"
      exit -1
    endif
    set bdy_file = $file
  endif

  set outfile = mozbc.out.$file
  if( -e $outfile ) then
    rm -f $outfile
  endif
  if( $bnd_typ == IC || $bnd_typ == BOTH ) then
    set msg = "Running mozbc for file $file"
    if( $bnd_typ == BOTH ) then
      set msg = "$msg and $bdy_file"
    endif
  else
    set msg = "Running mozbc for file $bdy_file"
  endif
  echo $msg

  if( -e run_mozbc_done ) then
    rm -f run_mozbc_done
  endif

#---------------------------------------
#  call mozbc utility
#---------------------------------------
    ./mozbc < $mozbc_inp > $outfile
#---------------------------------------
#  check signal file
#---------------------------------------
  if( ! -e run_mozbc_done ) then
    echo "There was an error running mozbc, check file $outfile; $prog terminating"
    rm -f wrfbdy_d01
    rm -f wrfinput_d01
    exit -1
  endif
#---------------------------------------
#  clean up
#---------------------------------------
  if( $bnd_typ == BC || $bnd_typ == BOTH ) then
    rm -f wrfbdy_d01
    if( $status ) then
      echo "Failed to rm wrfbdy_d01"
      exit -1
    endif
  endif
  if( $bnd_typ == IC ) then
    rm -f wrfinput_d01
    if( $status ) then
      echo "Failed to rm wrfinput_d01"
      exit -1
    endif
  endif
end

if( $bnd_typ == BC ) then
  rm -f wrfinput_d01 >& /dev/null
endif
if( -e run_mozbc_done ) then
  rm -f run_mozbc_done
endif

echo " "
echo "--------------------------------"
echo "run_mozbc successfully completed"
echo "--------------------------------"

exit 0
