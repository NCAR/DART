;=======================================================================
; subroutines for mopitt_extract_svd_transform
; based on mopitt_extract_intexb.pro
;=======================================================================
; Get the Vdata information from an HDF file
function get_vd, filename, varname
   file_id = h5f_open(filename)
   data_id = h5d_open(file_id, varname)
   var     = h5d_read(data_id)
   h5d_close, data_id
   return,var
end
;
; Get the Scientific Data information from an HDF file
function read_mopitt, filename, varname
   sd_id = hdf_sd_start(filename,/read)
   index = hdf_sd_nametoindex(sd_id, varname)
   sds_id = hdf_sd_select(sd_id,index)
   hdf_sd_getdata, sds_id, dat
   hdf_sd_endaccess, sds_id
   hdf_sd_end, sd_id
   return, dat
end
;
;=======================================================================
; main IDL routine
; needs read_mopitt, get_vd, 
; needs calc_avgker_v3, mopitt_v4_apriori.dat and read_aprior_dat_v3 for V3
; afa change output_nrows_leading to output_rows_leading
;=======================================================================
pro mopitt_co_cpsr_extract, inf, outf, bin_beg_sec, bin_end_sec, lon_min, lon_max, lat_min, lat_max 
;=======================================================================
; Code to read MOPITT data V3 or V4
; Outputs data to an ascii file for DART input
; Outputs station data to an ascii file for diagnostics
; written by Ave Arellano (NCAR) 
;
; Notes: 
; Need all the functions in the same directory.
; to run:
;   IDL> .r mopitt_extract_svd_transform.pro
;   IDL> mopitt_extract_svd_transform
;       
; But now, this is called by a shell script to process DART MOPITT obs
;   inf 		--> input file
;   outf 	--> output file
;   bin_beg_sec --> beginning seconds of the bin  (follows DART 6-hourly bins)
;   bin_end_sec --> end seconds of the bin (follows DART 6-hourly bins)
;   num_version  --> integer for MOPITT version (3 or 4)
;   what_cov     --> covariance (3 or 5) MOPITT has two versions for v4 
;   output_rows_leading --> how many leading components to assimilate?
;   sband  --> spectral band (if tir, nir or tirnir)
;   apm_no_transform --> switch to ignore scaling/scd transformation
;   when saving MOPITT data
;=======================================================================
; floating underflow in la_svd routines (compared output with matlab)
; seems to be very similar --suppress exception for now
!EXCEPT = 0
;
   print, 'IDL file in  ', inf
   print, 'IDL file out ', outf
   print, 'IDL bin_str  ', bin_beg_sec
   print, 'IDL bin end  ', bin_end_sec
   print, 'IDL lon min  ', lon_min
   print, 'IDL lon max  ', lon_max
   print, 'IDL lat min  ', lat_min
   print, 'IDL lat max  ', lat_max
;
   num_version=5
   version = 'v5'
   what_cov=3
   output_rows_leading=2
   sband='tirnir'
   apm_no_transform='true'
;
; note that mopittlev [0] is psurf (hPa) the 1000 is a placeholder
   mopittpress=[1000., 900., 800., 700., 600., 500., 400., 300., 200., 100.]
;
; note that covariance in v4 are to calculated
; prior error covariance is fixed --see V4 User's Guide
; set prior covariance parameters
   if ( what_cov eq 3 ) then begin
      prior_error = 0.3 ; in log 
   endif else begin
      prior_error = 0.5 ; in log
   endelse
   delta_pressure_lev = 100. ; hPa
   log10e = alog10(exp(1))
   C0     = (prior_error*log10e)^2
   Pc2    = delta_pressure_lev^2
   mopitt_dim = n_elements(mopittpress)
;
;=======================================================================
; QUALITY CONTROLS
; set QC here -- based on MOPITT Data Quality Statement
; these settings are very ad hoc !!! see personal comments
; edit if necessary

; dofs (i dont think dof is higher than 2 for MOPITT)
; from the pdfs of dofs, the threshold below appears to be 'outliers'
;
   case sband of
     'tir':     begin
                   dofs_threshold_low = 0.5
                   dofs_threshold_hi  = 2.0
                end
     'nir':     begin
                   dofs_threshold_low = 0.5
                   dofs_threshold_hi  = 1.0
                end
     'tirnir': begin
                   dofs_threshold_low = 0.5
                   dofs_threshold_hi  = 3.0
                end
   endcase
;
; pick daytime or nighttime or all -- there appears to be some contention
; whether or not there is bias in nighttime retrievals --
   sza_day = 90.0          ; all day data
   sza_nit = 180.0         ; all day and night data
;
; polar regions --there are potential biases near the poles
   day_lat_edge_1 = -70.0  ; 70S
   day_lat_edge_2 =  70.0  ; 70N
   nit_lat_edge_1 = -60.0  ; 60S
   nit_lat_edge_2 =  60.0  ; 60N
;
; retrieval error as fraction of its prior error
; this is very ad hoc based on percent apriori (post error/prior error)
   max_error_reduction_qc = 1.00  ; it's difficult to set this in log 
                                  ; make this 95% and let dofs be its qc
;
;=======================================================================
;
; Define MOPITT file name
   mopitt_input_file   = inf
   mopitt_output_file  = outf
;
; echo what we are processing here 
   print, 'IDL accessing MOPITT file: ', mopitt_input_file
   print, 'IDL writing to ascii file: ', mopitt_output_file
;
;=======================================================================
; read MOPITT file
;=======================================================================
; Read Seconds in Day
   name = '/HDFEOS/SWATHS/MOP02/Geolocation Fields/SecondsinDay'
   sec1 = get_vd(mopitt_input_file, name)
   nx   = long(n_elements(sec1)-1)
   sec  = float(sec1)
;
; Read Latitude;
   name = '/HDFEOS/SWATHS/MOP02/Geolocation Fields/Latitude'
   lat  = get_vd(mopitt_input_file, name)
;
; Read Longitude;
   name = '/HDFEOS/SWATHS/MOP02/Geolocation Fields/Longitude'
   lon  = get_vd(mopitt_input_file, name)
;
; Read Cloud Description
   name  = '/HDFEOS/SWATHS/MOP02/Data Fields/CloudDescription'
   cloud = get_vd(mopitt_input_file, name)
;
; Read Surface Pressure
   name  = '/HDFEOS/SWATHS/MOP02/Data Fields/SurfacePressure'
   psurf = get_vd(mopitt_input_file, name)
;
; Read Solar Zenith Angle
   name = '/HDFEOS/SWATHS/MOP02/Data Fields/SolarZenithAngle'
   sza  = get_vd(mopitt_input_file, name)
;
; Read Surface Indicator
   name = '/HDFEOS/SWATHS/MOP02/Data Fields/SurfaceIndex'
   sind = get_vd(mopitt_input_file, name)
;
; Read CO Total Column
   name  = '/HDFEOS/SWATHS/MOP02/Data Fields/RetrievedCOTotalColumn'
   cocol = get_vd(mopitt_input_file, name)
   cocol0 = reform(cocol[0,*])
   cocol1 = reform(cocol[1,*])
;
; Read DOFS
   name = '/HDFEOS/SWATHS/MOP02/Data Fields/DegreesofFreedomforSignal'
   dofs = get_vd(mopitt_input_file, name)
;
; Read Retrieved Non-Surface CO Mixing Ratio
   name   = '/HDFEOS/SWATHS/MOP02/Data Fields/RetrievedCOMixingRatioProfile'
   codata = get_vd(mopitt_input_file, name)
   comix    = reform(codata[0,*,*])
   comixerr = reform(codata[1,*,*])
;
; Read Retrieved CO Surface Mixing Ratio 
   name   = '/HDFEOS/SWATHS/MOP02/Data Fields/RetrievedCOSurfaceMixingRatio'
   codata = get_vd(mopitt_input_file, name)
   scomix    = reform(codata[0,*])
   scomixerr = reform(codata[1,*])
;
; Read Retrieval Averaging Kernel Matrix
   name   = '/HDFEOS/SWATHS/MOP02/Data Fields/RetrievalAveragingKernelMatrix'
   avgker = get_vd(mopitt_input_file, name)
;
; Read A Priori Surface CO Mixing Ratio 
   name   = '/HDFEOS/SWATHS/MOP02/Data Fields/APrioriCOSurfaceMixingRatio'
   codata = get_vd(mopitt_input_file, name)
   sperc    = reform(codata[0,*])
   spercerr = reform(codata[1,*])
;
; Read A Priori CO Mixing Ratio Profile
   name   = '/HDFEOS/SWATHS/MOP02/Data Fields/APrioriCOMixingRatioProfile'
   codata = get_vd(mopitt_input_file, name)
   perc    = reform(codata[0,*,*])
   percerr = reform(codata[1,*,*])
;
; Read Retrieval Error Covariance Matrix
   name      = '/HDFEOS/SWATHS/MOP02/Data Fields/RetrievalErrorCovarianceMatrix'
   covmatrix = get_vd(mopitt_input_file, name)
   print, 'APM: Completed MOPITT data read'
;
;=======================================================================
; Open output file
;=======================================================================
;
   openw, unit, mopitt_output_file, /get_lun 
;
;=======================================================================
; define/initialize other variables here
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
;=======================================================================
;
; Now, loop through each pixels
   for k = 0L, nx do begin
      allpix_count = allpix_count + 1.0
;
;=====================================================
      qstatus = 0
      mopittlev = mopittpress
;
; check number of levels
; mop_dim --> effective number of levels
      mop_dim = -999
      if (psurf[k] gt 900.0) then begin
         mopittlev = mopittpress
         mopittlev[0]=psurf[k]
         mop_dim = 10
      endif
      if (psurf[k] le 900.0) &&  (psurf[k] gt 800.0) then begin
         mopittlev = mopittpress
         mopittlev[1] = psurf[k]
         mop_dim = 9
      endif
      if (psurf[k] le 800.0) &&  (psurf[k] gt 700.0) then begin
         mopittlev = mopittpress
         mopittlev[2] = psurf[k]
         mop_dim = 8
      endif
      if (psurf[k] le 700.0) &&  (psurf[k] gt 600.0) then begin
         mopittlev = mopittpress
         mopittlev[3] = psurf[k]
         mop_dim = 7
      endif
      if (psurf[k] le 600.0) &&  (psurf[k] gt 500.0) then begin
         mopittlev = mopittpress
         mopittlev[4] = psurf[k]
         mop_dim = 6
      endif
      if (psurf[k] le 500.0) &&  (psurf[k] gt 400.0) then begin
         mopittlev = mopittpress
         mopittlev [5] = psurf[k]
         mop_dim = 5
      endif
      if (mop_dim lt 0) then begin
         print, 'MOPITT Surface Level Too High'
         qstatus = 1
      endif
;
; Construct prior covariance matrix
; based on MOPITT V4 User Guide
      covmat_ap = fltarr(mop_dim,mop_dim, /nozero)
      d_mop_dim = mopitt_dim-mop_dim
      for ii=d_mop_dim,mopitt_dim-1 do begin
         for jj=d_mop_dim,mopitt_dim-1 do begin
            tmp=C0*( exp( -((mopittlev[ii]-mopittlev[jj])^2)/Pc2) )
            covmat_ap[ii-d_mop_dim,jj-d_mop_dim]=tmp
         endfor
      endfor
;
; change notation (for my case)
      Ca = covmat_ap
;
; Invert the a priori covariance matrix
      invCa = invert(Ca, status)
;
; status = 0: successful
; status = 1: singular array
; status = 2: warning that a small pivot element was used and that
;             significant accuracy was probably lost.
; Status is usually 2
;
      if status eq 1 then begin
         qstatus = 1
      endif else begin
         qstatus = 0
      endelse
;
; change notation ( for clarity )
      mop_A = avgker[*,*,k]
;
; need to do transpose (IDL is column-major)
      mop_A = transpose(mop_A)
;
; truncate to effective number of levels
      A = fltarr(mop_dim,mop_dim, /nozero)
      A = mop_A[mopitt_dim-mop_dim:mopitt_dim-1,mopitt_dim-mop_dim:mopitt_dim-1]
      I = Identity(mop_dim)
;
; change notation (for my case)
      mop_Cx = covmatrix[*,*,k]
;
; need to do transpose (IDL is column-major)
      mop_Cx = transpose(mop_Cx)
;
; truncate to effective number of levels
      Cx = fltarr(mop_dim,mop_dim, /nozero)
      Cx = mop_Cx[mopitt_dim-mop_dim:mopitt_dim-1,mopitt_dim-mop_dim:mopitt_dim-1]
;
; assign dfs
      dfs = dofs[k]
;
;=====================================================
; APM: at this point AVE has full averaging kernal
;      and covariance matrixes
; QC: most qc applied here 
;=====================================================
;      print, 'IDL lon, lat ', lon[k],lat[k]
      if (dfs gt dofs_threshold_low && dfs lt dofs_threshold_hi && $
      ((sza[k] le sza_day && lat[k] gt day_lat_edge_1 && lat[k] lt day_lat_edge_2) || $
      (sza[k] ge sza_day && lat[k] gt nit_lat_edge_1 && lat[k] lt nit_lat_edge_2)) && $
      sec[k] ge bin_beg_sec && sec[k] le bin_end_sec && $
      lat[k] ge lat_min && lat[k] le lat_max && $
      lon[k] ge lon_min && lon[k] le lon_max && $
      qstatus eq 0) then begin 
;
;         print, dfs, sza[k], lat[k], lon[k], sec[k], psurf[k], qstatus
         qc_count = qc_count + 1.0
         co = fltarr(mop_dim,/nozero)
         coerr = fltarr(mop_dim,/nozero)
         priorerr = fltarr(mop_dim,/nozero)
         priorerrb = fltarr(mop_dim,/nozero)
         prior = fltarr(mop_dim,/nozero)
;
; assign file variables to profile quantities
         co[0] = scomix[k]*1e-9		; in mixing ratio now
         coerr[0] = scomixerr[k]*1e-9 
         prior[0] = sperc[k]*1e-9
         priorerrb[0] = spercerr[k]*1e-9 ; this is from MOPITT product
;
         ik_lev = 0
         for ik = 1,mopitt_dim-1 do begin
            if (perc[ik-1,k] gt 0.0 ) then begin
               ik_lev = ik_lev + 1
               prior[ik_lev]=perc[ik-1,k]*1e-9
               priorerrb[ik_lev]=percerr[ik-1,k]*1e-9 ; this is from MOPITT product
               co[ik_lev]=comix[ik-1,k]*1e-9
               coerr[ik_lev]=comixerr[ik-1,k]*1e-9 ; this is from MOPITT product
            endif
         endfor
;
; here we change to co instead of prior for the normalization point
; see notes on calculation of retrieval error
; the reason for this is that the prior and posterior should have
; the same normalization point for consistency
         for ik = 0,mop_dim-1 do begin
            priorerr[ik] =sqrt(Ca[ik,ik])*co[ik]/log10e 
         endfor
;
; before we used 700hPa or 500 hPa level for apriori contrib
; i think it's better to use the maximum error reduction instead
         error_reduction = fltarr(mop_dim,/nozero) 
         for ik = 0, mop_dim-1 do begin
            error_reduction[ik] = (priorerr[ik]-coerr[ik])/priorerr[ik]
         endfor
         max_error_reduction = max(error_reduction, max_index)
         apriori_contrib = max_error_reduction
;
; QC: check error reduction qc
         if (apriori_contrib le max_error_reduction_qc)  then begin
            qc_count2 = qc_count2 + 1.0
;
; now that we have Ca, Cx, and A, we need to get xa, x and xe
; change notation and make sure youre in log space
            xa = transpose(alog10(prior))
            x = transpose(alog10(co))
            I = Identity(mop_dim)
;
;======================================================================
;
; convert Cx (VMR) to fractional
; convert Ca (VMR) to fractional          
; recalculate xe (same as coerr but in fractional form)
            xe = fltarr(mop_dim, /nozero)
            for ik=0,mop_dim-1 do begin
               xe[ik] = sqrt(Cx[ik,ik])
            endfor
;
; recalculate xe_a (same as priorerr but in fractional form)
            xe_a =  fltarr(mop_dim, /nozero)
            for ik=0,mop_dim-1 do begin
               xe_a[ik] = sqrt(Ca[ik,ik])
            endfor
;
; 
; calculate prior term (I-A)xa in x = A x_t + (I-A)xa + Gey
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
            Cm = ImA##Ca##( I + transpose(AmI) ) 
;
; Ok. Here's the SVD approach
; ========================================================================
; (1) We need to scale the retrieval 
;     x = xa + A(x_t-xa) + ex or x = A x_t + (I-A) xa + ex
;     Cm  = < e_x ex^T > , < (Gey)(Gey)^T > 
;     sCm = ex^-1 
;     sCm x = ( sCm A ) x_t + ( sCm (I-A) xa )
; ========================================================================

; to do the scaling, we get the square root of Cm
; either doing cholesky decomposition or eigenvalue/vector decomposition
; turns out Cm is sometimes singular or not positive-definite
; so for now, we'll use SVD

;=======================================================================
; but first, precondition Cm and make it symmetric
; i think this is valid since most of the differences i see are
; very small --hence errors associated with numerical linear algebra
;
; for some numerical reason, Cm is not symmetric
; have problems doing eigenvalue decomposition with non-symmetric matrix
            for ik = 0,mop_dim-1 do begin
               for ijk = ik,mop_dim-1 do begin
                  Cm[ik,ijk]=Cm[ijk,ik]
               endfor
            endfor
;
; save original Cm for debugging later on
;
; APM: Looks like I want to save Cm_dummy
; APM: Here are the unscaled, unrotated equivalents
;       A   - averaging kernal
;       Cm  - measurement error covariance
;       xa  - prior
;       x   - retrieval
;
            Cm_dummy = Cm
;
; get eigenvalues and eigenvectors
            eigenvalues = la_eigenql(Cm, EIGENVECTORS = eigenvectors, status=status)
            if status ne 0 then begin
               print, 'Cm la_eigengl did not converge'
               print, Cm, status
               qstatus = 1
            endif
;
; APM: not necessary because Cm is a covariance matrix which is 
; symmetric, positive semidefinite => all real eigenvalues >= 0
;
; check for complex values
; la_eigenql readme says it outputs real eigenvectors
; but in matlab, sometimes, eigenvalues can be complex 
; so might as well check for complex numbers
            for ik=0,mop_dim-1 do begin
               ftype = size(eigenvalues[ik], /type)
               if (ftype eq 6  ) then begin
                  qstatus = 1
               endif
            endfor
;
; APM: not necessary - see preceding APM comment
;
; precondition by removing all negative eigenvalues
; and replace by floating point precision
            eval = Identity(mop_dim) 
            for ik=0,mop_dim-1 do begin
               if (eigenvalues[ik] lt (Machar()).eps  ) then begin
                  eval[ik,ik]=(Machar()).eps
               endif else begin
                  eval[ik,ik]=eigenvalues[ik]
               endelse
            endfor
;
; reconstruct covariance matrix
            Cm = transpose(eigenvectors)##eval##(eigenvectors)
;
;=======================================================================
; whoops, that was tough.
; now take the inverse square-root of Cm
; from eigenvalue decomposition 
; S = L D L^T
; S^1/2 = L D^1/2 L^T
; S^-1/2 = L D^-1/2 L^T
; this is similar to SVD --we opt to use SVD here for convenience
; S = R D L^T 
; svd(S) = sqrt(eig(S#S))
; we can also do a cholesky decomposition which Rodgers suggested
; spectral decomposition
; S = T^T T
; S^-1/2 = T^-1
;=======================================================================
;
; APM: measurement covariance based scaling
; 
; status = 0: The computation was successful
; status =>0: The computation did not converge. The status value specifies how many superdiag did not
;             converge to zero
            la_svd, Cm, s_Cm, u_Cm, v_Cm, status=status
            if status ne 0 then begin
               print, 'Cm svd did not converge'
               print, Cm, status
               qstatus = 1
            endif
;
; get the square root of Cm
            sqrt_s_Cm = fltarr(mop_dim, /nozero)
            for ik = 0,mop_dim-1 do begin
               if ( s_Cm[ik] le 0.0 )  then begin 
                  s_Cm[ik] = (Machar()).eps
               endif
               sqrt_s_Cm[ik] = sqrt(s_Cm[ik])
            endfor
;
; get the inverse square root of Cm
            inv_s  = Identity(mop_dim)
            for ik = 0,mop_dim-1 do begin 
               inv_s[ik,ik] = 1.0/sqrt_s_Cm[ik]
            endfor
            sCm    = v_Cm##inv_s##transpose(u_Cm)
;
; APM: a priori covariance based scaling
; 
; status = 0: The computation was successful
; status =>0: The computation did not converge. The status value specifies how many superdiag did not
;             converge to zero
            la_svd, Ca, s_Ca, u_Ca, v_Ca, status=status
            if status ne 0 then begin
               print, 'Ca svd did not converge'
               print, Ca, status
               qstatus = 1
            endif
;
; get the square root of Ca
            sqrt_s_Ca = fltarr(mop_dim, /nozero)
            for ik = 0,mop_dim-1 do begin
               if (s_Ca[ik] lt 0)  then begin 
                  s_Ca[ik] = (Machar()).eps
               endif 
               sqrt_s_Ca[ik] = sqrt(s_Ca[ik])
            endfor
            sqrt_Ca = Identity(mop_dim)
            for ik = 0,mop_dim-1 do begin 
               sqrt_Ca[ik,ik] = sqrt_s_Ca[ik]
            endfor
            sCa = v_Ca##sqrt_Ca##transpose(u_Ca) 
;
; =================================================================
; scale the whole expression
; =================================================================
;
; APM this is the measurement covariance based scaling
;
; scale A using sCm
            sA = sCm##A
;
; scale Cm --> now should be identity 
; the line below is a check 
; sCms = sCm##Cm##transpose(sCm)
; scale Ca --> this wont be use but may be useful for debugging
            sCas = sCm##Ca##(sCm)
;
; force sCms to be identity --this has been checked
            sCms = Identity(mop_dim)
;
; scale x
            sx = sCm##x
;
; scale prior term
            sImAxa = sCm##ImAxa
;
; ========================================================================
; (2) Get SVD of scaled Ax and rotate retrieval to its maximum information
; sA##sCa = USV^T
; U^T ( sCm x ) = U^T ( sCm A ) x_t +  U^T (sCm (I-A) xa ) 
; U^T ( sx ) = U^T ( sA ) x_t +  U^T ( sImAxa ) 
; ========================================================================

; from Migliorini et al 2008, they rotated the scaled covariance
; <H'PH'^T> for EIG which is similar to getting the square root of Ca in SVD
; see Stefano's email regarding P (2/10/2010)
; i think it makes more sense to have the singular vectors of Ax rather than A
; it accounts for both the variability in x and sensitivity of the retreival
;
; APM: Ave calculates SVD for Ax (averaging kernal time retrieval) as
; opposed to A (averaging kernal).
;
; APM: NOTE TOO: sAx is the result of a RHS a priori covariance based scaling of
; sA (the scaled averaging kernal) where the first scaling is the
; measurement covariance based scaling
; 
            sAx = sA##sCa
;
; take the svd of sA
;
; status = 0: The computation was successful
; status =>0: The computation did not converge. The status value specifies how many superdiag did not
;             converge to zero
            la_svd, sAx, S, U, V, status=status
            if status ne 0 then begin
               print, 'sA svd did not converge'
               print, sA, status
               qstatus = 1
            endif
;
; calculate transpose(U)A
            transA = transpose(U)##sA 
;
; need to check U if it is reflected rather than translated
; note that transA is row major now
; that means you need to sum the columns --
; here we assumed that the sum of the averaging kernel rows
; should be at least positive, if it is negative then
; it is reflected and we should take the opposite sign
            sum_uA = total(transA,1)
;        
; change sign of rows of U (in matlab this is really the 
; singular vector in columns)
            for ik=0,mop_dim-1 do begin
               if (sum_uA[ik] le 0.0) then begin
                  for ikk=0,mop_dim-1 do begin
                     U[ik,ikk] = U[ik,ikk]*(-1.0)
                  endfor
               endif
            endfor
;
; recalculate transA
            transA = transpose(U)##sA
;        
; rotate scaled Cx and Ca
; again this is a check if it's really identity
; Cmn = transpose(U)##sCms##U
; rotated Ca is not really used but might be useful for debugging
            Can = transpose(U)##sCas##U
;
;force Cmn to be identity -- this has been checked 
            Cmn = Identity(mop_dim)
;
; calculate new y 
            yn = transpose(U)##sx
;           
; rotate prior term
            transImAxa = transpose(U)##sImAxa
; 
; get new errors
            e2 = diag_matrix(Cm)
            e2n = diag_matrix(Cmn)
            e2a = diag_matrix(Can)
            for ik=0, mop_dim-1 do begin
               if (e2n[ik] lt 0) then begin
                  qstatus = 1.0
               endif else begin
                  e2n[ik] = sqrt(e2n[ik])
                  e2a[ik] = sqrt(e2a[ik])
               endelse
            endfor
;
; ========================================================================
; (3) Truncate to dofs level
; now look at transA and pick n retrievals corresponding to dofs
; ========================================================================
            nrows_leading0 = ceil(dfs)
;
; alternatively find the point where the variance explain is >95%
; 95% is arbitrary
            varsA = fltarr(mop_dim, /nozero)
            sumsA = 0
            for ik=0,mop_dim-1  do begin
               sumsA = sumsA + transA[ik,ik]
            endfor
            for ik=0,mop_dim-1 do begin 
               varsA[ik] = transA[ik,ik]/sumsA
            endfor
            cumsumA = 0
            for ik=0,mop_dim-1 do begin 
               cumsumA = cumsumA + varsA[ik]
            endfor
            nrows_leading1 = where(cumsumA ge 0.95)
            if ( nrows_leading1 ne nrows_leading0 ) then begin
               nrows_leading = max([nrows_leading0,nrows_leading1])
            endif else begin
               nrows_leading = nrows_leading1
            endelse
;
; what needs to be assimilated
            case output_rows_leading of
               0: begin
                  valid_nrows = nrows_leading
                  output_start_row = 0
                  output_end_row = nrows_leading-1
               end
               1: begin
                  valid_nrows = 1
                  output_start_row = 0
                  output_end_row = 0
               end
               2: begin
                  valid_nrows = mop_dim
                  output_start_row = 0
                  output_end_row = mop_dim-1
               end
            endcase
;
; ========================================================================
; finally output the variables in ascii
; yes, in ascii -- this is really for my own convenience
; it's not computationally efficient but it makes easier debugging
; besides this is only a temp file for DART to use
; ========================================================================
; QC -> only output values with qstatus=0
            if (qstatus eq 0) then  begin
               allqc_count = allqc_count + 1.0
;
; Code to write no_transform data
               if (apm_no_transform eq 'true') then begin
                  printf, unit, 'NO_SVD_TRANS', sec[k], lat[k], lon[k], $
                  mop_dim, dfs, format='(a15,16(e14.6))'
; effective pressure levels
                  printf, unit, mopittlev(10-fix(mop_dim):9), format='(10(e14.6))'
; retrieval
                  printf, unit, x+9.0, format='(10(e14.6))'
; prior retrieval
                  printf, unit, xa+9.0, format='(10(e14.6))'
; averaging kernel
                  printf, unit, transpose(A), format='(100(e14.6))'
; prior error covariance
                  printf, unit, transpose(Ca), format='(100(e14.6))'
; retrieval error covariance
                  printf, unit, transpose(Cx), format='(100(e14.6))'
; measurment error covariance
                  printf, unit, transpose(Cm), format='(100(e14.6))'
; total column
                  printf, unit, cocol0[k],cocol1[k], format='(2(e14.6))'
               endif    ; apm_no_transform
            endif          ; QC - qstatus for numerical issues
         endif             ; apriori contribution
      endif                ; QC - most quality control
   endfor		   ; MOPITT pixels (k data) 
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
   print, 'Count (%)'
   print, allqc_count*100.0/allpix_count
;
   print, '================================'
   print, 'IDL SVD transformation DONE for ', mopitt_input_file
   print, '================================'
   print, ' '
end    ; end of mopitt_extract_svd_transform
