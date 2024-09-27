;=======================================================================
; subroutines for iasi_extract_svd_transform
; based on iasi_extract_intexb.pro
;=======================================================================
; Get the Vdata information from an HDF file
function get_vd, filename, varname
       file_id = hdf_open(filename, /read)
       vd_ref = hdf_vd_find(file_id, strtrim(varname,2))
       vdata=hdf_vd_attach(file_id,vd_ref)
       nread=hdf_vd_read(vdata,var)
       hdf_vd_detach, vdata
return,var
end

; Get the Scientific Data information from an HDF file
function read_iasi, filename, varname

  sd_id = hdf_sd_start(filename,/read)
  index = hdf_sd_nametoindex(sd_id, varname)
  sds_id = hdf_sd_select(sd_id,index)
  hdf_sd_getdata, sds_id, dat
  hdf_sd_endaccess, sds_id
  hdf_sd_end, sd_id
  return, dat
end

;=======================================================================
; main IDL routine
; needs read_iasi, get_vd, 
; needs calc_avgker_v3, iasi_v4_apriori.dat and read_aprior_dat_v3 for V3
;=======================================================================
pro iasi_co_cpsr_extract, inf, outf, bin_beg_sec, bin_end_sec, lon_min, lon_max, lat_min, lat_max
!EXCEPT=2
;=======================================================================
;  Code to read IASI data (similar to MOPITT V5)
;  Outputs data to an ascii file for DART input
;  Outputs station data to an ascii file for diagnostics
;  Notes: 
;      Need all the functions in the same directory.
;    to run:
;     IDL> .r iasi_extract_svd_transform.pro
;     IDL> iasi_extract_svd_transform
;       
;  But now, this is called by a shell script to process DART IASI obs
;
;  inf 		--> input file
;  outf 	--> output file
;  bin_beg      --> beginning hour of the bin  (follows DART 6-hourly bins)
;  bin_end      --> end hour of the bin (follows DART 6-hourly bins)
;  when saving IASI data
;=======================================================================
;
; floating underflow in la_svd routines (compared output with matlab)
; seems to be very similar --suppress exception for now
;
log10e = alog10(exp(1))
;
; debug level (debug=0 means no debug output)
debug = 0
;
; QUALITY CONTROLS
; set QC here -- based on IASI Data Quality Statement
; these settings are very ad hoc !!! see personal comments
; edit if necessary
;
; dofs (NEED TO ADJUST FOR IASI)
dofs_threshold_low = .5
dofs_threshold_hi  = 3.
;
; pick daytime or nighttime or all -- there appears to be some contention
; whether or not there is bias in nighttime retrievals --
sza_day = 90.0 ;all day profiles
sza_nit = 180.0 ;all profiles
;
; polar regions --there are potential for biases near the poles
day_lat_edge_1 = -70.0  ; 70S
day_lat_edge_2 =  70.0  ; 70N
nit_lat_edge_1 = -60.0  ; 60S
nit_lat_edge_2 =  60.0  ; 60N
;
; retrieval error as fraction of its  prior error
; this is very ad hoc based on percent apriori (post error/prior error)
max_error_reduction_qc = 0.05
;
; convert bins into seconds
;
; Define IASI file name
iasi_input_file   = inf
iasi_output_file  = outf
;
; Print input and output file names
print, 'APM: Input file: ', iasi_input_file
print, 'APM: Output file: ', iasi_output_file
print, 'APM: bin_beg_sec: ', bin_beg_sec
print, 'APM: bin_end_sec: ', bin_end_sec
;
;=======================================================================
; read IASI file
;=======================================================================
;
; Read Observation Time
name = 'Observation Time UTC HMS'
obs_time = get_vd(iasi_input_file, name)
;print, 'Observation Time UTC HMS ',obs_time[0]
;
; Read Seconds in Day
name = 'Seconds in Day'
obs_sec = get_vd(iasi_input_file, name)
nx = long(n_elements(obs_sec)-1)
;print, 'Seconds in Day ',obs_sec[0]
;
; Read Pressure grid;
name = 'Pressure Grid'
press = read_iasi(iasi_input_file, name)
print, 'APM: Read ',name
;
; Read Solar Zenith Angle
name = 'Solar Zenith Angle'
sza = get_vd(iasi_input_file, name)
print, 'APM: Read ',name
;
; Read Surface Pressure
name = 'Surface Pressure'
psurf = get_vd(iasi_input_file, name)
print, 'APM: Read ',name
;
; Read Surface Emissivity
name = 'Surface Emissivity'
semiss = get_vd(iasi_input_file, name)
print, 'APM: Read ',name
;
; Read Retrieved CO Mixing Ratio Profile
name = 'Retrieved CO Mixing Ratio Profile'
codata = read_iasi(iasi_input_file, name)
covmr = reform(codata[0,*,*])
covmr_err = reform(codata[1,*,*])
print, 'APM: Read ',name
;
; Read Retrieved CO Surface Mixing Ratio 
name = 'Retrieved CO Surface Mixing Ratio'
codata = read_iasi(iasi_input_file, name)
scovmr = reform(codata[0,*])
scovmr_err = reform(codata[1,*])
print, 'APM: Read ',name
;
; Read Retrieved CO Total Column
name  = 'Retrieved CO Total Column'
cocoldata = read_iasi(iasi_input_file, name)
cototcol = reform(cocoldata[0,*])
cototcol_err = reform(cocoldata[1,*])
print, 'APM: Read ',name
;
; Read Retrieval Averaging Kernel Matrix
name = 'Retrieval Averaging Kernel Matrix'
avgker = read_iasi(iasi_input_file, name)
print, 'APM: Read ',name
;
; Read Retrieval Error Covariance Matrix
name = 'Retrieval Error Covariance Matrix'
covmatrix = read_iasi(iasi_input_file, name)
print, 'APM: Read ',name
;
; Read A Priori CO Mixing Ratio Profile
name = 'A Priori CO Mixing Ratio Profile'
codata = read_iasi(iasi_input_file, name)
copriorvmr = reform(codata[0,*,*])
copriorvmr_err = reform(codata[1,*,*])
print, 'APM: Read ',name
;
; Read A Priori Surface CO Mixing Ratio 
name = 'A Priori CO Surface Mixing Ratio'
codata = read_iasi(iasi_input_file, name)
scoprior = reform(codata[0,*])
scoprior_err = reform(codata[1,*])
print, 'APM: Read ',name
;
; Read Degrees of Freedom for Signal
name = 'Degrees of Freedom for Signal'
dofs = get_vd(iasi_input_file, name)
print, 'APM: Read ',name
;
; Read DEM Altitude
name = 'DEM Altitude'
dem_alt = get_vd(iasi_input_file, name)
print, 'APM: Read ',name
;
; Read Time
name = 'Time'
time = get_vd(iasi_input_file, name)
print, 'Time ',time[0]
;
; Read Latitude;
name = 'Latitude'
lat  = get_vd(iasi_input_file, name)
print, 'APM: Read ',name
;
; Read Longitude;
name = 'Longitude'
lon  = get_vd(iasi_input_file, name)
print, 'APM: Read ',name
;
; Read Retrieved CO Column Profile
name = 'Retrieved CO Column Profile'
codata = read_iasi(iasi_input_file, name)
cocol = reform(codata[0,*,*])
cocol_err = reform(codata[1,*,*])
print, 'APM: Read ',name
;
; Read Apriori CO Column Profile
name = 'Apriori CO Column Profile'
codata = read_iasi(iasi_input_file, name)
copriorcol = reform(codata[0,*,*])
copriorcol_err = reform(codata[1,*,*])
print, 'APM: Read ',name
;   
print, 'APM: Complted data read '
;
; Open output file
unit=10
openw, unit, iasi_output_file, /get_lun
;
; Define/initialize other parameters
; log10 conversion for d ln(VMR) = d(VMR)/VMR = d log10(VMR) /log10(e)
log10e = alog10(exp(1)) 
;
; qc count (most of qc - dofs, time, sza, partial qc, high lat)
qc_count = 0.0
;
; qc for apriori contribution (subset of qc_count)
qc_count2 = 0.0
;
; qc for all (subset of qc_count2)
allqc_count = 0.0
;
; all pixel count
allpix_count = 0.0
;
; loop through each pixel
for k = 0L, nx do begin 
;
; Update counter
   allpix_count = allpix_count + 1.0
;  
; Set the status flag
   qstatus = 0
;
; get DEM altitude for this pixel
; this should be zero for sea level pixels
   alt=dem_alt(k)
   idx_str=fix(alt)
;
; iasi_dimm --> effective number of levels
   iasi_dim = 19   
   iasi_dimm = iasi_dim-idx_str
;
; AVERAGING KERNEL: change notation (for clarity)
   iasi_A = fltarr(iasi_dim,iasi_dim, /nozero)
   iasi_A = avgker[*,*,k]
;
; Truncate to effective number of levels (not needed for IASI)
   A = fltarr(iasi_dimm,iasi_dimm, /nozero)
   A = iasi_A[iasi_dim-iasi_dimm:iasi_dim-1,iasi_dim-iasi_dimm:iasi_dim-1]
   A = transpose(A)
;
; RETRIEVAL ERROR COVARIANCE: change notation (for clarity)
   iasi_Cx = covmatrix[*,*,k]
;
; Truncate to effective number of levels (not needed for IASI)
   Cx = fltarr(iasi_dimm,iasi_dimm, /nozero)
   Cx = iasi_Cx[iasi_dim-iasi_dimm:iasi_dim-1,iasi_dim-iasi_dimm:iasi_dim-1]
   Cx = transpose(Cx)
   Co_col = fltarr(iasi_dimm, /nozero)
   Co_col = cocol[iasi_dim-iasi_dimm:iasi_dim-1,k]
   Co_vmr = fltarr(iasi_dimm, /nozero)
   Co_vmr = covmr[iasi_dim-iasi_dimm:iasi_dim-1,k]
   Co_prior_col = fltarr(iasi_dimm, /nozero)
   Co_prior_col = copriorcol[iasi_dim-iasi_dimm:iasi_dim-1,k]
   Co_prior_vmr = fltarr(iasi_dimm, /nozero)
   Co_prior_vmr = copriorvmr[iasi_dim-iasi_dimm:iasi_dim-1,k]
   air_column = fltarr(iasi_dimm, /nozero)
   for i = 0,iasi_dimm-1 do begin
      air_column(i)=Co_col(i)/Co_vmr(i)
   endfor
;   print, 'APM air_column: ', air_column(0:iasi_dimm-1)
   for i = 0,iasi_dimm-1 do begin
      air_column(i)=Co_prior_col(i)/Co_prior_vmr(i)
   endfor
;   print, 'APM air_column: ',air_column(0:iasi_dimm-1)
;
; Convert averaging kernel to VMR (i,j error application corrected)
; Conversion needed but done in DART obs converter
;   for i = 0,iasi_dimm-1 do begin 
;      for j = 0,iasi_dimm-1 do begin
;         print, 'APM A, ac_i, ac_j ',A(i,j), air_column(i), air_column(j) 
;         A(i,j) = A(i,j) * air_column(i) / air_column(j) 
;      endfor
;   endfor
;
; Convert Cx from percent to VMR (not needed)
;   for i = 0,iasi_dimm-1 do begin 
;      for j = 0,iasi_dimm-1 do begin 
;         Cx(i,j) = Cx(i,j) * Co_prior_vmr(i) * Co_prior_vmr(j) 
;      endfor
;   endfor
;
; A PRIORI ERROR COVARIANCE: (this is a placeholder)
   Ca = fltarr(iasi_dimm,iasi_dimm, /nozero)
   Czero=(0.3*log10e)^2.
   Pref=100.0
   zmin=200.
   for i = 0,iasi_dimm-1 do begin
      for j = 0, iasi_dimm-1 do begin
         if( i eq j ) then begin
            Ca(i,j) = Czero
         endif else begin
            press_i=(press(i+idx_str,k)+press(i+idx_str+1,k))/2.
            press_j=(press(j+idx_str,k)+press(j+idx_str+1,k))/2.
;            if( press_i lt zmin) then begin
;               press_i = zmin
;            endif 
;            if( press_j lt zmin) then begin
;               press_j = zmin 
;            endif   
;            Ca(i,j) = Czero/exp(((press_i-press_j)/Pref)^2.)
            Ca(i,j) = Czero
         endelse
      endfor
   endfor
;
; DEGREES OF FREEDOM FOR SIGNAL:
   dfs = dofs[k]
;
; APM: at this point we have full averaging kernal
; APM: this is pre-DART QA/QC and may nedd to be revised to get more
; obs to DART
;      and covariance matrix
;  QC: most qc applied here
; 
   hr=float(fix(obs_time[k]/10000))
   min=float(fix(obs_time[k]/100)-hr*100)
   scc=float(obs_time[k]-hr*10000-min*100)
   tod_sec=hr*60*60+min*60+scc   
;
   if debug eq 1 then begin
      print, 'obs_time ',obs_time[k]
;      print, 'hr-mn-sc ',hr,min,scc
;      print, 'DFS: ',dfs, dofs_threshold_low, dofs_threshold_hi
;      print, 'ZAN: ',sza[k], sza_day
;      print, 'LAT: ',lat[k], day_lat_edge_1, day_lat_edge_2
;      print, 'LAT: ',lat[k], nit_lat_edge_1, nit_lat_edge_2
      print, 'SEC: ',obs_sec[k], bin_beg_sec, bin_end_sec
      print, 'SEC: ',tod_sec, bin_beg_sec, bin_end_sec
      print, 'LAT: ',lat[k], lat_min, lat_max
      print, 'LON: ',lon[k], lon_min, lon_max
      print, ' '
   endif
   if( $
      ( dfs ge dofs_threshold_low ) && ( dfs le dofs_threshold_hi ) && $
      ((( sza[k] lt sza_day ) && ( lat[k] gt day_lat_edge_1 ) && ( lat[k] lt day_lat_edge_2 )) || $
      (( sza[k] ge sza_day ) && ( lat[k] gt nit_lat_edge_1 ) && ( lat[k] lt nit_lat_edge_2 ))) && $
      ( obs_sec[k] ge bin_beg_sec ) and (obs_sec[k] le bin_end_sec ) and $
      ( lat[k] ge lat_min ) && ( lat[k] le lat_max ) && $
      ( lon[k] ge lon_min ) && lon[k] le ( lon_max ) && $
      ( qstatus eq 0 )) then begin 
;
; update qc count
      qc_count = qc_count + 1.0
;
; make iasi_dimm - element profiles
      co = fltarr(iasi_dimm,/nozero)
      coerr = fltarr(iasi_dimm,/nozero)
      prior = fltarr(iasi_dimm,/nozero)
      priorerr = fltarr(iasi_dimm,/nozero)
      priorerrb = fltarr(iasi_dimm,/nozero)
      pressure = fltarr(iasi_dimm+1,/nozero)
;
; save free atmospheric data in profiles
      for ik=idx_str,iasi_dim-1 do begin
         ik_lev = ik-idx_str
         co[ik_lev]=Co_vmr[ik_lev]
         coerr[ik_lev]=covmr_err[ik,k] 
         prior[ik_lev]=Co_prior_vmr[ik_lev]
         priorerrb[ik_lev]=copriorvmr_err[ik,k] 
         pressure[ik_lev]=press[ik,k]
      endfor
      if (iasi_dimm eq iasi_dim) then begin
         pressure[iasi_dim]=press[iasi_dim,k]
      endif else begin
         pressure[iasi_dimm]=press[iasi_dim,k]
      endelse 
;
; update qc_count2
      qc_count2 = qc_count2 + 1.0
;
; now that we have Ca, Cx, and A, we need to get xa, x and xe
; change notation and make sure youre in log space
; both v3 and v4 report VMR units for xa and x
      xa=transpose(prior)
      x=transpose(co)
      I = Identity(iasi_dimm)
;
; calculate prior term (I-A)xa  in x = A x_t + (I-A)xa + Gey
; needed for obs/forward operator (prior term of expression)
      ImA = (I-A)
      AmI = (A-I)
      ImAxa = ImA##xa
;
; calculate Cm --here Cx = Cm + Cs
; where Cx is the posterior error covariance
; Cm = < (Gey)(Gey)^T > = GSeG^T, Se is measurement noise covariance
; Cm -> error due to measurement noise
; Cs -> error due to smoothing (application of prior)
; Cm = (I-A) Ca (I + (A-I)^T) = Cx (K^T Se^-1 K ) Cx
; Cs = Cx Sa^-1 Cx 
; x = Cx K^T Se^-1 (y - Kxa)
; A = Cx K^T Se^-1 K = G K 
;
; THIS IS PLACE HOLDER
      Cm = ImA##Ca##( I + transpose(AmI) ) 
;
; A   - averaging kernal
; Cm  - measurement error covariance
; xa  - prior
; x   - retrieval
;
; output paramters 
      valid_nrows = iasi_dimm
      output_start_row = 0
      output_end_row   = iasi_dimm-1
;
; finally output the variables in ascii
; QC -> only output values with qstatus=0
      if (qstatus eq 0) then  begin
;
; update counter
         allqc_count = allqc_count + 1.0
;
; note that the format is for iasi_dim levels
         if ((k gt 0) && (obs_sec[k-1] lt obs_sec[k])) then begin
;            printf, unit, 'NO_SVD_TRANS', obs_sec[k], lat[k], lon[k], $
;            iasi_dimm, dfs, format='(a12,1x,3(g15.7,1x),i2,1x,g15.7)'
            print, 'APM SEC, TOD_SEC ',obs_sec[k],tod_sec
            printf, unit, 'NO_SVD_TRANS', obs_sec[k], lat[k], lon[k], $
            iasi_dimm, dfs
;
; effective height levels
;            printf, unit, pressure[0:iasi_dimm], format='(20(g15.7,1x))'
            printf, unit, pressure[0:iasi_dimm]
;
; retrieval
;            printf, unit, x[0:iasi_dimm-1], format='(19(g15.7,1x))'
            printf, unit, x[0:iasi_dimm-1]
;
; retrieval prior vmr
;            printf, unit, xa[0:iasi_dimm-1], format='(19(g15.7,1x))'
            printf, unit, xa[0:iasi_dimm-1]
;
; retrieval prior col
;            printf, unit, Co_prior_col[0:iasi_dimm-1], format='(19(g15.7,1x))'
            printf, unit, Co_prior_col[0:iasi_dimm-1]
;
; averaging kernal
;            printf, unit, A[0:iasi_dimm-1,0:iasi_dimm-1], format='(361(g15.7,1x))'
            printf, unit, A[0:iasi_dimm-1,0:iasi_dimm-1]
;   print, 'APM: ak ', A(0,0:iasi_dimm-1)
;   print, 'APM: prior_col ',Co_prior_col(0:iasi_dimm-1)
;   print, 'APM: prior_vmr ',Co_prior_vmr(0:iasi_dimm-1)
; prior error covariance (PLACE HOLDER)
;            printf, unit, Ca[0:iasi_dimm-1,0:iasi_dimm-1], format='(361(g15.7,1x))'
            printf, unit, Ca[0:iasi_dimm-1,0:iasi_dimm-1]
;
; retrieval error covariance
;            printf, unit, Cx[0:iasi_dimm-1,0:iasi_dimm-1], format='(361(g15.7,1x))'
            printf, unit, Cx[0:iasi_dimm-1,0:iasi_dimm-1]
;
; measurment error covariance (PLACE HOLDER)
;            printf, unit, Cm[0:iasi_dimm-1,0:iasi_dimm-1], format='(361(g15.7,1x))'
            printf, unit, Cm[0:iasi_dimm-1,0:iasi_dimm-1]
;
; total column data
;            printf, unit, cototcol[k],cototcol_err[k], format='(2(g15.7,1x))'
            printf, unit, cototcol[k],cototcol_err[k]
             endif
          endif   ; QC --qstatus for numerical issues       
;       endif      ; QC -- max error reduction
    endif       ; most qc testing   
endfor		; pixels (k data) 
;    
; close files
close, unit
;
; print counters
print, 'BIN TIME'
print, bin_beg_sec, bin_end_sec
print, 'QC Count'
print, allqc_count 
print, 'ALL Pixels Count'
print, allpix_count
print, '% Count'
print, allqc_count*100.0/allpix_count
print, 'ACCEPTED - first level qc/qc'
print, qc_count
print, 'ACCEPTED - max error reduction qc/qc'
print, qc_count2
;
print, '================================'
print, 'IDL SVD transformation DONE for ', iasi_input_file
print, '================================'
print, ' '
;
end    ; end of iasi_extract_no_svd_transform
