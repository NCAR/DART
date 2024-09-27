;=======================================================================
; main IDL routine

pro create_ascii_IASI_O3, inf_col, inf_err, inf_vmr, outf, bin_beg_sec, bin_end_sec, lon_min, lon_max, lat_min, lat_max, datestring
   !quiet=1
   close,/all
;
;=======================================================================
;  inf_col      --> input file
;  inf_err 	--> input file
;  inf_vmr 	--> input file
;  outf 	--> output file
;  bin_beg      --> beginning hour of the bin  (follows DART 6-hourly bins)
;  bin_end      --> end hour of the bin (follows DART 6-hourly bins)
;  when saving IASI data
;=======================================================================
;
;
; Print input and output file names
   print, 'APM: Input file Column: ', inf_col
   print, 'APM: Input file ERROR: ', inf_err
   print, 'APM: Input file VMR: ', inf_vmr
   print, 'APM: Output file: ', outf
   print, 'APM: bin_beg_sec: ', bin_beg_sec
   print, 'APM: bin_end_sec: ', bin_end_sec
   print, 'APM: lon_min: ', lon_min
   print, 'APM: lon_max: ', lon_max
   print, 'APM: lat_min: ', lat_min
   print, 'APM: lat_max: ', lat_max
   print, 'APM: datestring: ', datestring
   sza_day=90.
   sza_nit=180.
;
   day_lat_edge1=-70.
   day_lat_edge2=70.
   nit_lat_edge1=-60.
   nit_lat_edge2=60.
;
; Open and readd the observatino error covariance file
   nl=41
   nlp=nl+1
   apcov_tmp = fltarr(nl, nl)
   openr,1,'IASI_apcov.dat'
   readf,1,apcov_tmp
   close,1
   apcov=apcov_tmp[0:nl-1, 0:nl-1]   ; relative units
;
; Open the output file
   openw,11,outf
   print, 'IASI O3 Output file is ',outf
;
   n=100000.
   nsub=500
   lon=fltarr(n)
   lat=fltarr(n)
   hms=fltarr(n)
   z0=fltarr(n)
   sza=fltarr(n)
   rms=fltarr(n)
   dofs=fltarr(n)
   vecAKcol = fltarr(n,nl*nl)
   press=fltarr(n,nlp)
   vmr=fltarr(n,nl)
   apvmr=fltarr(n,nl)
   col=fltarr(n,nl)
   bias=fltarr(n)
   apmtmp=fltarr(n)
   apcol = fltarr(n,nl)
   error = fltarr(n,nl)
   surfT = fltarr(n)
   cloud = fltarr(n)
;
; Open and read the VMR input file
   openr,1,inf_vmr
;
   dummy=''
   count=0L
   isub=0
   while not eof(1) do begin
;   while isub ne nsub do begin
      readf,1, dummy
      parts  = str_sep(dummy,',')
      parts = strcompress(parts,/remove_all)
      k=where(parts ne '')
      parts=parts(k)
      lon[count] = parts[0]
      lat[count] = parts[1]
      hms[count] = parts[2]
      z0[count] = parts[3]
      rms[count] = parts[4]
      bias[count] = parts[5]
      dofs[count] = parts[6]
      sza[count] = parts[10]
      cloud[count] = parts[11]
      surfT[count] = parts[12]
      press[count,*] = parts[13:54]
      vmr[count,*] = parts[55:95]
      apvmr[count,*] = parts[96:136]
      count=count+1L
      isub=isub+1L
   endwhile
   close,1
;
; Open and read the COLUMN input file
   openr,1,inf_col
;
   dummy=''
   count2=0L
   isub=0
   while not eof(1) do begin
;   while isub ne nsub do begin
      readf,1, dummy
      parts  = str_sep(dummy,',')
      parts = strcompress(parts,/remove_all)
      k=where(parts ne '')
      parts=parts(k)
      vecAKcol[count2,*] = parts[55:1735]
      col[count2,*] = parts[1736:1776]
      apcol[count2,*] = parts[1777:1817]
      count2=count2+1L
      isub=isub+1L
   endwhile
   close,1
;
; Open and read the ERROR input file
   openr,1,inf_err
;
   dummy=''
   count3=0L
   isub=0
   while not eof(1) do begin
;   while isub ne nsub do begin
      readf,1, dummy
      parts  = str_sep(dummy,',')
      parts = strcompress(parts,/remove_all)
      k=where(parts ne '')
      parts=parts(k)
      error[count3,*] = parts[13:53]   ; [%]
      count3=count3+1L
      isub=isub+1L
   endwhile
   close,1
;
; COLUMN file can have double the entries, Catherine had a bug when
; extracting the data.
   if count ne count2 or count ne count3 then begin
      if count ne count2/2. or count ne count3 then stop
   endif
;
; Process the IASI O3 data
   nprof = count
;
;   iprof=nprof
;   print, 'bias ',bias[iprof]
;   print, 'rms ',rms[iprof]
;   print, 'cloud ',cloud[iprof]
;   print, 'hms ',hms[iprof]
;   print, 'sza ',sza[iprof]
;   print, 'dofs ',dofs[iprof]
;
   for iprof=0, nprof-1 do begin
       hh=fix(hms[iprof])/10000
       mm=(fix(hms[iprof])-hh*10000)/100
       ss=hms[iprof]-hh*10000-mm*100
       iasi_sec=hh*60.*60.+mm*60.+ss
;
; select good quality data over our domain
;      if bias[iprof] lt 3.0 and cloud[iprof] le 5 and rms[iprof] le 3.0e-8 and $
;      (sza[iprof] lt sza_day and lat[iprof] gt day_lat_edge1 and lat[iprof] lt day_lat_edge2) or $
;      (sza[iprof] ge sza_day and lat[iprof] gt nit_lat_edge1 and lat[iprof] lt nit_lat_edge2) and $
;      lat[iprof] ge 7 and lat[iprof] le 54 and lon[iprof] ge -176 and lon[iprof] le -50 then begin
;
      if lat[iprof] ge lat_min and lat[iprof] le lat_max $
      and lon[iprof] ge lon_min and lon[iprof] le lon_max $
      and iasi_sec ge bin_beg_sec and iasi_sec lt bin_end_sec  then begin
         nlev_use=nl
         nlev_usep=nlp
         lev_use=findgen(nlp)
         lev1=0
         not_a_lev= where(finite(press[iprof,*],/NAN))
         if not_a_lev[0] ne -1 then begin
            nlev_usep=nlp-n_elements(not_a_lev)
            nlev_use=nlev_usep-1
            lev1=max(not_a_lev)+1
         end
;
         retlev=reform(press[iprof,lev1:nlp-1])
         alt_retlev_all = findgen(nlp)
         alt_retlev_ag = alt_retlev_all[lev1:nlp-1]-alt_retlev_all[lev1]
         alt_retlev = alt_retlev_all[lev1:nlp-1]                 
         retprs_mid=fltarr(nlev_use)
         retalt_mid=fltarr(nlev_use)
         for li=0,nlev_usep-2 do begin
            retprs_mid(li)=(retlev(li)+retlev(li+1))/2.
            retalt_mid(li)=(alt_retlev_ag(li)+alt_retlev_ag(li+1))/2.
         end
;         print, ' '
;         print, 'iprof ',iprof,nprof-1
;         print, 'not_a_lev[0] ',not_a_lev[0]
;         print, 'n_elements ',n_elements(not_a_lev)
;         print, 'nlev_use ',nlev_use
;         print, 'lev1 ',lev1 
;         print, 'psf, prs_1 ',retlev(0),retprs_mid(0) 
;
         iasi_vmr=reform(vmr[iprof,lev1:nl-1])
         iasi_col=reform(col[iprof,lev1:nl-1])
         ap_vmr=reform(apvmr[iprof,lev1:nl-1])
         ap_col=reform(apcol[iprof,lev1:nl-1])
         AKcol=fltarr(nlev_use,nlev_use)
         AKvmr=fltarr(nlev_use,nlev_use)
         ind=0
         for li=0,nlev_use-1 do begin
            AKcol(0:nlev_use-1,li)=vecAKcol(iprof,ind:ind+nlev_use-1)
            ind=ind+nlev_use
         end
;
; Calculate air_column
         air_column=fltarr(nlev_use)
         for li=0,nlev_use-1 do begin
            air_column(li)=ap_col(li)/ap_vmr(li)
         end 
;
; Convert averaging kernel to VMR
         for li=0,nlev_use-1 do begin
            for lj=0,nlev_use-1 do begin
               AKvmr(li,lj)=AKcol(li,lj)/air_column(lj)*air_column(li)
            end
         end 
; 
         ineg = where(finite(AKcol,/NAN))
         if ineg[0] ne -1 then begin
            print,'lat,lon,hms,iprof ',lat[iprof],lon[iprof],hms[iprof],iprof
            print, 'nlev_use,lev1 ',nlev_use,lev1
            for lj=0,nlev_use-1 do begin
               print, 'column', lj
               print, 'AKcol ',AKcol(0:nlev_use-1,lj)
            endfor
            print, 'iasi_col ',iasi_col(0:nlev_use-1)
            print, 'ap_col ',ap_col(0:nlev_use-1)
            stop
         endif 
;
         prior=reform(apcol[iprof,lev1:nl-1])
         apcov_col=reform(apcov[lev1:nl-1, lev1:nl-1])
         apcov_vmr=reform(apcov[lev1:nl-1, lev1:nl-1])
;
; Convert the error covariance
         for i=0,nlev_use-1 do begin
            for j=0,nlev_use-1 do begin
               apcov_col(i,j)=apcov_col(i,j)*ap_col(i)*ap_col(j)
            endfor
         endfor
         for i=0,nlev_use-1 do begin
            for j=0,nlev_use-1 do begin
               apcov_vmr(i,j)=apcov_vmr(i,j)*ap_vmr(i)*ap_vmr(j)
            endfor
         endfor
;
         okay=1
         if okay eq 1 then begin                    
            nll=nl-lev1
            I = Identity(nll)
            ImAcol = I - AKcol
            ImAvmr = I - AKvmr
            AcolmI = AKcol - I
            AvmrmI = AKvmr - I
;            print, 'dim I',size(I)
;            print, 'dim ImAcol',size(ImAcol)
;            print, 'dim ImAvmr',size(ImAvmr)
;            print, 'dim AKcol',size(AKcol)
;            print, 'dim AKvmr',size(AKvmr)
;            print, 'dim ap_col',size(ap_col)
;            print, 'dim ap_vmr',size(ap_vmr)
;            print, 'dim apcov_col',size(apcov_use)
;            print, 'ImA_col ', size(ImAcol)
;            print, 'ap_col ', size(ap_col)
            ImAxa_col = ImAcol##ap_col
            ImAxa_vmr = ImAvmr##ap_vmr
            Cm_col = ImAcol##apcov_col##(I + transpose(AcolmI))
            Cx_col = ImAcol##apcov_col
            Cm_vmr = ImAvmr##apcov_vmr##(I + transpose(AvmrmI))
            Cx_vmr = ImAvmr##apcov_vmr
;            for i=0,nlev_use-1 do begin
;               print, 'Cm_col ',i,Cm_col(i,0:nlev_use-1)
;               print, 'Cx_col ',i,Cx_col(i,0:nlev_use-1)
;               print, 'Cm_vmr ',i,Cm_vmr(i,0:nlev_use-1)
;               print, 'Cx_vmr ',i,Cx_vmr(i,0:nlev_use-1)
;            endfor
;
; location and time info and also retrieval level and number of levels
; and dfs for surface -100 mb
            printf,11, 'NO_SVD_TRANS', float(iasi_sec), lat[iprof], lon[iprof], float(lev1), float(nlev_use)
;
; retrieval surface pressure       
            printf, 11, retlev(0)
; retrieval altitude mid-levels       
            printf, 11, retalt_mid(0:nll-1)
; retrieval pressure mid-levels
            printf, 11, retprs_mid(0:nll-1)
; AK matrix
            printf, 11, AKcol(0:nll-1,0:nll-1)
; a priori profile (column) 
            printf, 11, ap_col(0:nll-1)
; a priori profile (VMR) 
            printf, 11, ap_vmr(0:nll-1)
; IASI retrieval profile (column)
            printf, 11, iasi_col(0:nll-1)
; IASI retrieval profile (VMR)
            printf, 11, iasi_vmr(0:nll-1)
; prior error covariance (column)
            printf, 11, apcov_col(0:nll-1,0:nll-1)
; measurement error covariance (column)
            printf, 11, Cm_col(0:nll-1,0:nll-1)
; retrieval error covariance (column)
            printf, 11, Cx_col(0:nll-1,0:nll-1)
; prior error covariance (VMR)
            printf, 11, apcov_vmr(0:nll-1,0:nll-1)
; measurement error covariance (VMR)
            printf, 11, Cm_vmr(0:nll-1,0:nll-1)
; retrieval error covariance (VMR)
            printf, 11, Cx_vmr(0:nll-1,0:nll-1)
;            if iprof eq 0 then begin
;               for i=0,nlev_use-1 do begin
;                  print, 'AKcol ',i,AKcol(i,0:nlev_use-1)
;                  print, 'Cx_vmr ',i,Cx_vmr(i,0:nlev_use-1)
;               endfor
;            endif
         endif
      endif
   endfor
   close,11
end
