; Find the mean state space increment, output the fields to a single mean file
; that can be used to make plots

; Incoming arguments are:
; 1) fname (analysis_increment.nc)
; 2) fout  (mean_increment.nc)

; This script is called from the diagnostic_obs.sh script as: 
; ncl fname=analysis_increment_{d??}.nc  fout=mean_increments_{d??}.nc mean_increment.ncl

begin

  if (.not. isvar("fname")) then
    print("ERROR in mean_increment.ncl: fname not provided.")
    exit
  end if
  if (.not. isvar("fout")) then
    print("ERROR in mean_increment.ncl: fout not provided.")
    exit
  end if

  flist = systemfunc("ls -1 ../*/" + fname + " 2>/dev/null")
  nfils = dimsizes(flist)

  if (nfils .le. 0) then
    print("ERROR in mean_increment.ncl: no files found matching ../*/" + fname)
    exit
  end if

; If  we only want say the last 7 days, then grab only the last 28
; here we practice with 3 days
 anl_days = 7
 ntimes = anl_days*4
 if (nfils .gt. ntimes) then
   tempf = flist(nfils-ntimes:nfils-1) 
   delete(flist)
   flist = tempf
   nfils = ntimes
   delete(tempf)
 end if 
 fil   = addfiles(flist, "r")
 ListSetType(fil, "join")

 pull_2D_field_names = (/"T2", "Q2", "U10", "V10", "PSFC"/)
 pull_3D_field_names = (/"U", "V", "THM", "QVAPOR"/)
 npulls = dimsizes(pull_2D_field_names)

; Below will dump out the data to a file for replotting later
    cnew = addfile(fout,"c")
; Work through 2D fields
  do i=0,npulls-1
   print("   Extracting 2d variable "+pull_2D_field_names(i))
   do fil_num=0,nfils-1
;   print(" reading file "+flist(fil_num))
; Dimensions are ncljoin, Time, south_north, west_east
; Copy zero is the ensemble mean 
    pull_var = fil[fil_num]->$pull_2D_field_names(i)$(:,:,:,:)
    dims = dimsizes(pull_var)
    if (fil_num .eq. 0) then  ; first iteration, make var
     alltimes_var = new ( (/nfils,dims(2),dims(3)/), typeof(pull_var) )
    end if
;    printVarSummary(pull_var)
    alltimes_var(fil_num,:,:) = pull_var(0,0,:,:)
;   printVarSummary(alltimes_var)
    delete(pull_var)
   end do
; average over time (first dimension)
  mean_alltimes_var = dim_avg_n(alltimes_var,0)
; standard deviation over time (first dimension)
  stdv_alltimes_var = dim_stddev_n(alltimes_var,0)
; write to new file
   varname ="mean_"+pull_2D_field_names(i)
   cnew->$varname$ = mean_alltimes_var
   delete(varname)
   varname ="stdv_"+pull_2D_field_names(i)
   cnew->$varname$ = stdv_alltimes_var
   delete(varname)
   delete(alltimes_var)
   delete(mean_alltimes_var)
   delete(stdv_alltimes_var)
   delete(dims)
  end do

; work through 3D fields
 npulls = dimsizes(pull_3D_field_names)
  do i=0,npulls-1
   print("   Extracting 3d variable "+pull_3D_field_names(i))
   do fil_num=0,nfils-1
;    print(" reading file "+flist(fil_num))
; dimensions are ncljoin, Time, level, south_north, west_east
; copy zero is the ensemble mean
    pull_var = fil[fil_num]->$pull_3D_field_names(i)$(:,:,:,:,:)
    dims = dimsizes(pull_var)
    if (fil_num .eq. 0) then  ; first iteration, make var
     alltimes_var = new ( (/nfils,dims(2),dims(3),dims(4)/), typeof(pull_var) )
    end if
;    printVarSummary(pull_var)
    alltimes_var(fil_num,:,:,:) = pull_var(0,0,:,:,:)
    delete(pull_var)
   end do
; average over time (first dimension)
  mean_alltimes_var = dim_avg_n(alltimes_var,0)
; standard deviation over time (first dimension)
  stdv_alltimes_var = dim_stddev_n(alltimes_var,0)
; write to new file
   varname ="mean_"+pull_3D_field_names(i)
   cnew->$varname$ = mean_alltimes_var
   delete(varname)
   varname ="stdv_"+pull_3D_field_names(i)
   cnew->$varname$ = stdv_alltimes_var
   delete(varname)
   delete(alltimes_var)
   delete(mean_alltimes_var)
   delete(stdv_alltimes_var)
   delete(dims)
  end do

end
