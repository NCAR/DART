; DART software - Copyright UCAR. This open source software is provided
; by UCAR, "as is", without charge, subject to all terms of use at
; http://www.image.ucar.edu/DAReS/DART/DART_download
;
; DART $Id$

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
; The following script is a substitute for the modified ideal.exe routine used
; to generate an ensemble initial condition set appropriate for idealized
; radar assimilation cases with WRF/DART. 
; The basic flow is as follows: 
; 1) run geogrid.exe (in WPS) to generate the grid information
; 2) run sounding perturbation routine to generate perturbed base state profiles
; 3) run the shell script using the 'stock' ideal.exe routine to generate template 
;    files for each perturbed sounding and then calling this script to swap the 
;    grid information into the perturb input files and add local perturbations
; 4) run wrf_to_dart to generate an initial condition set for filter
;
; Routines originally developed by David Dowell and Altug Aksoy. Conversion to ncl 
; by Glen Romine.
; Updated Dec. 2010 to account for the modified model_mod for wind locations
;
; Note that running the default ideal.exe inserts a thermal bubble (3K) in the 
; center of the domain
;
; Things to still improve:
;  1) Allow for a grid that doesn't start with location '1' (read from patch attributes)

; New global variables
;  member_num      = 1                                                     ; Ensemble member number - should be passed in at the command line
  wrfin_file      = "./wrfinput_d01"  ; Path to input file
  wrfgeo_file     = "../geo_em.d01.nc" ; file with the geography          ; location from prepare_ideal_temp directory
  bubble_dir      = "/Volumes/joshua/romine/assim/radar_ic"               ; directory with bubble location files (randominit=3)
                                                                          ; output from WPS geogrid.exe
                    ; later versions of this script might aim to combine the creation of the geo file 
                    ; with the script executing this routine
  morrison        = 1             ; set to '1' if you want to have REFL_10CM and FALL_SPD_Z_WEIGHTED
; geo_em.d01.nc
; original fields in ideal.input.filter
; dt is in config_flags%
  dt          = 2.0               ; dt     -> model time step (s)
  islat       = 0                 ; islat  -> specify lat/lon of grid? (1=yes - hardcoded, 0=no - read from files)
; xmin, ymin
  xymin       = (/0.0, 0.0/)      ; for islat=1 only, xmin, ymin   -> x, y coordinates (km) of southwest grid corner
; glat, glon
  glatlon     = (/0.0, 0.0/)      ; for islat=1 only, glat, glon   -> latitude (deg), longitude (deg) corresponding to 
                                  ; (x=0, y=0)
; clat, clon
  clatlon     = (/35.704, -97.321/) ; for islat=0 only, clat, clon -> latitude (deg), longitude (deg) corresponding to
                                  ; center of the domain (i.e. radar location)
  ishgt       = 1                 ; ishgt  -> specify grid height? (1=yes - specify below, 0=no - read from file)
; galt = 1000.*galtkm
  galtkm      = 0.350             ; for ishgt=1 only, galtkm       -> grid height (km MSL)
  ismapfct    = 0                 ; ismapfct -> specify map variables? (1=yes - hardcoded, 0=no - read from file)
  mapfactor   = 0.981             ; for ismapfct=1 only, mapfactor -> map factor to be hardcoded
  sinalpha    = 0.0               ; for ismapfct=1 only, sinalpha  -> local sine of map rotation to be hardcoded
  cosalpha    = 1.0               ; for ismapfct=1 only, cosalpha  -> local cosine of map rotation to be hardcoded
  iscorls     = 1                 ; iscorls -> specify Coriolis terms? (1=yes - specify below, 0=no - read from file)
  coriolis_e  = 0.0               ; for iscorls=1 only, coriolis_e -> Coriolis cosine latitude term to be hardcoded
  coriolis_f  = 0.0               ; for iscorls=1 only, coriolis_f -> Coriolis sine latitude term to be hardcoded
; the following 3 are in config_flags%
  julyr       = 1981              ; julyr  -> Julian year
  julday      = 137               ; julday -> Julian day
  gmt         = 16.3333           ; gmt    -> Time of day GMT (hours, including fraction?)
  randominit  = 1                 ; randominit  -> random bubble locations? (0=wrf default, 1=uniform, 2=gaussian,
                                  ; 3=read from file)
; rbubh, rbubv
  rbubhv      = (/7.5, 1.5/)      ; rbubh (km), rbubv (km)
  nb          = 3                 ; number of warm bubbles
  randomtemp  = 0                 ; randomtemp  -> random temperature perturbation? (0=no, 1=yes)
; below: delt, stdt, lowt, hight
  magtpert    = (/5.0, 0.0, 0.0, 0.0/)  ; For non-random T perturbations: magnitude of temperature perturbations (K)
                                  ; next two numbers ignored.  For random T perturbations: mean magnitude, 
                                  ; std deviation, lower threshold, upper threshold all in (K)
; code has a dummy variable here - no match
  numseed       = 0               ; For randominit=2 only, Number of most probable bubble locations (seeds)
                                  ; ignored for other options; up to 5 allowed currently
  stddist       = 0.0             ; For randominit=2 only, Standard deviation of Gaussian distribution of bubble 
                                  ; location in horizontal grid
; ibmin, ibmax
  ibval         = (/4, 14/)       ; For randominit=1: i grid indices of perturbation region (west-east direction)
                                  ; For randominit=2: i grid indices of locations of most probable bubbles (separate 
                                  ; with ","; as many as number of seeds)
; jbmin, jbmax
  jbval         = (/6, 16/)       ; For randominit=1: j grid indices of perturbation region (south-north direction)
                                  ; For randominit=2: j grid indices of locations of most probable bubbles (separate 
                                  ; with ","; as many as number of seeds)
; kbmin, kbmax
  kbval         = (/1, 6/)        ; k grid indices of perturbation region (height); Gaussian distribution not applied 
                                  ; in the vertical
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Goal here will be to take output from a standard (release version) ideal.exe from
; wrf v3.1.1 (wrfinput_d01) and create an ensemble initial condition set. The default
; ideal.exe sets a number of fields that need to be overwritten with the values from
; WPS or based on flags/fields set above. Also, a 3K thermal bubble is in the center
; of the domain.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Extract the individual variables with consistent names to original routine from 
; definitions provided above
;
  xmin  = xymin(0)
  ymin  = xymin(1)
  glat  = glatlon(0)
  glon  = glatlon(1)
  clat  = clatlon(0)
  clon  = clatlon(1)
  galt  = 1000.*galtkm
  rbubh = rbubhv(0)
  rbubv = rbubhv(1)
  delt  = magtpert(0)
  stdt  = magtpert(1)
  lowt  = magtpert(2)
  hight = magtpert(3)
  if (randominit .eq. 1) then
    ibmin = ibval(0)
    ibmax = ibval(1)
    jbmin = jbval(0)
    jbmax = jbval(1)
  end if
  if (randominit .eq. 2) then
; numseed specified guassian bubble locations
  end if
  kbmin = kbval(0)
  kbmax = kbval(1)
; 
; Other constants and definitions
  rearth       = 6367.0
  rad2deg      = 57.29578
  g            = 9.81
  r_d          = 287.04
  r_v          = 461.6
  rvovrd       = r_v/r_d
  p1000mb      = 100000.
  t0           = 300.
  cp           = 7.*r_d/2.
  cv           = cp - r_d 
  cvpm         = -1.*(cv/cp)
; Some of the field names in the geography file differ from the wrfinput
; field names - so these are paired up in the set below. IF YOU CHANGE
; THESE, YOU NEED TO CHECK ASSUMED ARRAY POSITIONS BELOW.
  geog_field_names = (/  "XLAT_M",   "XLONG_M",      "CLAT",     "CLONG",   \
                          "HGT_M", "MAPFAC_M",   "MAPFAC_U",  "MAPFAC_V",   \
                      "MAPFAC_MX", "MAPFAC_MY", "MAPFAC_UX", "MAPFAC_UY",   \
                      "MAPFAC_VX", "MAPFAC_VY",  "SINALPHA",  "COSALPHA",   \
                              "E",         "F", \
                        "XLONG_U",  "XLONG_V",     "XLAT_U",    "XLAT_V"/)
  tmpl_field_names = (/    "XLAT",     "XLONG",      "CLAT",     "CLONG",  \
                            "HGT",  "MAPFAC_M",  "MAPFAC_U",  "MAPFAC_V",  \
                      "MAPFAC_MX", "MAPFAC_MY", "MAPFAC_UX", "MAPFAC_UY",  \
                      "MAPFAC_VX", "MAPFAC_VY",  "SINALPHA",  "COSALPHA",  \
                              "E",         "F", \
                        "XLONG_U",  "XLONG_V",     "XLAT_U",    "XLAT_V"/)
; index values:                 0,           1,           2,           3, 
;                               4,           5,           6,           7,
;                               8,           9,          10,          11,
;                              12,          13,          14,          15,
;                              16,          17,          
;                              18,          19,          20,          21
; Open the template and geography files and replace the global fields
  in_template_file = addfile(wrfin_file+".nc","w")
  in_geog_file = addfile(wrfgeo_file,"r")
; Grab some need constants from the template file:
  global_atts=getvaratts(in_template_file)
  wrftitle = stringtochar(in_template_file@$global_atts(0)$)
  wrfver = chartostring(wrftitle(19:20))
  print ("wrf model version "+wrfver)
  nxm  = in_template_file@$global_atts(3)$-1
  nym  = in_template_file@$global_atts(4)$-1
  nzm  = in_template_file@$global_atts(5)$-1
  dx   = in_template_file@$global_atts(6)$
  dy   = in_template_file@$global_atts(7)$
; check PHB - staggered or unstaggered?
  ztop = in_template_file->PHB(0,nzm-1,0,0)/g  ; approximation - should be good enough
                                               ; since this is only used for thermal
                                               ; bubble size
  dxkm = 0.001*dx
  dykm = 0.001*dy
  dzkm = 0.001*ztop/(nzm-1)
  print("(dxkm, dykm, dzkm) : ("+dxkm+", "+dykm+", "+dzkm+")")
;
; Write in the date, time step, lat/lon information
      truelat1 = 30.
      truelat2 = 60.
      moad_cen_lat = clat
      stand_lon = clon
      map_proj = 1
      replace_att     = True
      replace_att@JULYR = julyr
      replace_att@JULDAY = julday
      replace_att@GMT = gmt
      replace_att@DT = dt
      replace_att@CEN_LAT = clat
      replace_att@CEN_LON = clon
      replace_att@TRUELAT1 = truelat1
      replace_att@TRUELAT2 = truelat2
      replace_att@MOAD_CEN_LAT = clat
      replace_att@STAND_LON = clon
      replace_att@MAP_PROJ = map_proj
      fileattdef(in_template_file,replace_att)
      delete(replace_att)

   if (islat .eq. 1) then
      print("Lat/lon values will be hardcoded...")
      print("xmin, ymin = "+xmin+" "+ymin)
      print("glat, glon = "+glat+" "+glon)
      x=fspan(xmin, (nxm-1)*dxkm,nxm)
      y=fspan(ymin, (nym-1)*dykm,nym)
; match up dimensions
; lat/lon approx, lat is a function of y only, but lon is a function
; of both x and y - so 2d.
      grid2d = new((/nxm,nym/),float,-999.)
      x2d    = conform(grid2d,x,0)
      y2d    = conform(grid2d,y,1)
      glat2d = conform(y2d,glat,-1)
      glon2d = conform(x2d,glon,-1)
      const2d = conform(y2d,rad2deg/rearth,-1)
      xlattmp = glat2d + const2d * y2d
      const2da = conform(x2d,rad2deg,-1)
      const2db = conform(x2d,rearth,-1)
      xlontmp = glon2d + const2da * (x2d / (const2db*cos(0.5*(glat2d+xlattmp)/const2da)) )
      delete(grid2d)
      delete(const2d)
      delete(const2da)
      delete(const2db)
      delete(x2d)
      delete(y2d)
      delete(glon2d)
      delete(glat2d)
; Have 2d lat/lon, but field is 3D (includes time dimension) but there is
; only one time in wrfinput file. Below, just grab whatever is in the file
; now and overwrite with new fields. Assuming XLAT = CLAT, XLON = CLON
      i=0
      swap_var = in_template_file->$tmpl_field_names(i)$
      swap_var(0,:,:) = (/xlattmp/)
      in_template_file->$tmpl_field_names(i)$=(/swap_var(:,:,:)/)
      if ( wrfver .eq. "V3" ) then
        in_template_file->$tmpl_field_names(i+2)$=(/swap_var(:,:,:)/)
      end if
      delete(swap_var)
      i=1
      swap_var = in_template_file->$tmpl_field_names(i)$
      swap_var(0,:,:) = (/xlontmp/)
      in_template_file->$tmpl_field_names(i)$=(/swap_var(:,:,:)/)
      if ( wrfver .eq. "V3" ) then
       in_template_file->$tmpl_field_names(i+2)$=(/swap_var(:,:,:)/)
      end if
      delete(swap_var)
; Below is another block of code setting more constants and attributes
; the lat/lon for the center of the grid box is found - or near the center
; in the event the grid index is an odd # 
      if (mod(nxm,2) .eq. 0) then
        i1 = (nxm) / 2
        i2 = i1
      else
        i1 = (nxm-1) / 2
        i2 = i1+1
      end if
      if (mod(nym,2) .eq. 0) then
        j1 = (nym) / 2
        j2 = j1
      else
        j1 = (nym-1) / 2
        j2 = j1+1
      end if
      cen_lat = 0.25*(xlattmp(i1,j1)+xlattmp(i2,j1)+xlattmp(i1,j2)+xlattmp(i2,j2))
      cen_lon = 0.25*(xlontmp(i1,j1)+xlontmp(i2,j1)+xlontmp(i1,j2)+xlontmp(i2,j2))
      truelat1 = 30.
      truelat2 = 60.
      moad_cen_lat = cen_lat
      stand_lon = cen_lon
      map_proj = 1
      replace_att     = True
      replace_att@CEN_LAT = cen_lat
      replace_att@CEN_LON = cen_lon
      replace_att@TRUELAT1 = truelat1
      replace_att@TRUELAT2 = truelat2
      replace_att@MOAD_CEN_LAT = moad_cen_lat
      replace_att@STAND_LON = stand_lon
      replace_att@MAP_PROJ = map_proj
      fileattdef(in_template_file,replace_att)
      delete(xlattmp)
      delete(xlontmp)
      delete(replace_att)
   else
      print("Lat/lon values will be read from geog file...")
      if (wrfver .eq. "V3") then
       nv = 3
      else
       nv = 1
      end if
      do i=0,nv
       swap_var = in_geog_file->$geog_field_names(i)$(:,:,:)
       in_template_file->$tmpl_field_names(i)$=(/swap_var(:,:,:)/)
       delete(swap_var)
      end do
   end if
   print(" ")
; done with islat
   if (ishgt .eq. 1) then
      print("Surface elevation values will be hardcoded...")
      print("Surface elevation set to "+galt+" m")
      i=4
      swap_var = in_template_file->$tmpl_field_names(i)$
      swap_var(:,:,:) = galt
      in_template_file->$tmpl_field_names(i)$=(/swap_var(:,:,:)/)
      delete(swap_var)
; Updates to grid hgt requires adjusting the lower boundary base state geopotential height
      swap_var = in_template_file->PHB
      swap_var(0,0,:,:) = g * galt
      in_template_file->PHB=(/swap_var(:,:,:,:)/)
      delete(swap_var)
   else
      print("Surface elevation values will be read from geog file...")
      i=4
      swap_var = in_geog_file->$geog_field_names(i)$(:,:,:)
      in_template_file->$tmpl_field_names(i)$=(/swap_var(:,:,:)/)
   end if
   print(" ")
; done with ishgt
   if (ismapfct .eq. 1) then
      print("Map-related variables will be hardcoded...")
      print("Map factor set to "+mapfactor+", sina set to "+sinalpha+  \
            ", cosa set to "+cosalpha)
      if (wrfver .eq. "V3") then
       nv = 13
      else
       nv = 7
      end if
      do i=5,nv
        swap_var = in_geog_file->$geog_field_names(i)$(:,:,:)
        swap_var(:,:,:) = mapfactor
        in_template_file->$tmpl_field_names(i)$=(/swap_var(:,:,:)/)
        delete(swap_var)
      end do
; MF_VX_INV may also need to be written - not sure what this field is used for, 
; but it's not in the geog file, and is supposed to just be the inverse of the 
; MAPFAC_VX field - so maybe just hard code to write this. 
      if (wrfver .eq. "V3") then
       swap_var(:,:,:) = 1./mapfactor
       in_template_file->MF_VX_INV=(/swap_var(:,:,:)/)
       delete(swap_var)
      end if
      i=14
      swap_var = in_geog_file->$geog_field_names(i)$(:,:,:)
      swap_var(:,:,:) = sinalpha
      in_template_file->$tmpl_field_names(i)$=(/swap_var(:,:,:)/)
      delete(swap_var)
      i=15
      swap_var = in_geog_file->$geog_field_names(i)$(:,:,:)
      swap_var(:,:,:) = cosalpha
      in_template_file->$tmpl_field_names(i)$=(/swap_var(:,:,:)/)
   else
      print("Map-related factors will be read from geog file...")
      if (wrfver .eq. "V3") then
       nv = 13
      else
       nv = 7
      end if
      do i=5,nv
        swap_var = in_geog_file->$geog_field_names(i)$(:,:,:)
        in_template_file->$tmpl_field_names(i)$=(/swap_var(:,:,:)/)
        delete(swap_var)
      end do
; MF_VX_INV may also need to be written - not sure what this field is used for, 
; but it's not in the geog file, and is supposed to just be the inverse of the 
; MAPFAC_VX field - so maybe just hard code to write this. 
      if (wrfver .eq. "V3") then
       swap_var = in_template_file->MAPFAC_VX(:,:,:)
       swap_var(:,:,:) = 1.0/(swap_var(:,:,:))
       in_template_file->MF_VX_INV=(/swap_var(:,:,:)/)
       delete(swap_var)
      end if
      do i=14,15
        swap_var = in_geog_file->$geog_field_names(i)$(:,:,:)
        in_template_file->$tmpl_field_names(i)$=(/swap_var(:,:,:)/)
        delete(swap_var)
      end do
   end if
   print(" ")
; done with ismapfct
   if (iscorls .eq. 1) then
      print("Coriolis values will be hardcoded...")
      i=16
      swap_var = in_template_file->$tmpl_field_names(i)$
      swap_var(:,:,:) = coriolis_e
      in_template_file->$tmpl_field_names(i)$=(/swap_var(:,:,:)/)
      delete(swap_var)
      i=17
      swap_var = in_template_file->$tmpl_field_names(i)$
      swap_var(:,:,:) = coriolis_f
      in_template_file->$tmpl_field_names(i)$=(/swap_var(:,:,:)/)
      delete(swap_var)
   else
      print("Coriolis values will be read from geog file...")
      do i=16,17
        swap_var = in_geog_file->$geog_field_names(i)$(:,:,:)
        in_template_file->$tmpl_field_names(i)$=(/swap_var(:,:,:)/)
        delete(swap_var)
      end do
   end if
; copy over the U and V point lat/lon
   do i=18,21
      swap_var = in_geog_file->$geog_field_names(i)$(:,:,:)
      in_template_file->$tmpl_field_names(i)$=(/swap_var(:,:,:)/)
      delete(swap_var)
   end do
   print(" ")
; 
; Time to add in the bubbles/noise:
; randominit = 0, do nothing, standard wrf bubble which is already in place
; any other randominit value and need to clear out the bubble and return
; to a horizontally homogeneous base state 
; randominit = 1, throw in a specified number of random bubbles to a region
; randominit = 2, specified gaussian bubbles
; randominit = 3, read file 'bubble_locations' 
;
; After any bubble changes, recalculate the density and balance hydrostatically
; 
; 
; Initializes seeds for the random generators "gaussian" and "uniform".

  quote = inttochar(34)
  stringval = "date "+quote+"+%s"+quote
  iseed1 = stringtointeger(systemfunc("date "+quote+"+%s"+quote))  ; sec since 010170 00:00
  iseed2 = stringtointeger(systemfunc("date "+quote+"+%S"+quote))+member_num  ; 0-60 +ens member num
  random_setallseed(iseed1,iseed2)
  if (randominit .eq. 0) then
    print("=====================================")
    print("Standard WRF bubble.")
  else
; remove the thermal bubble in fields T, ph (note ideal also updates alt, al but these
; fields are not placed in the wrfinput file by default.)
; Will just grab the profile from the southwest corner and overwrite the 3D field
    swap_var = in_template_file->T
    vprofile = swap_var(0,:,0,0)
    vprofile3d = conform(swap_var(0,:,:,:),vprofile,0)
    swap_var(0,:,:,:) = (/vprofile3d(:,:,:)/)
    in_template_file->T = swap_var
    delete(swap_var)
    delete(vprofile)
    delete(vprofile3d)
    swap_var = in_template_file->PH
    vprofile = swap_var(0,:,0,0)
    vprofile3d = conform(swap_var(0,:,:,:),vprofile,0)
    swap_var(0,:,:,:) = (/vprofile3d(:,:,:)/)
    in_template_file->PH = swap_var
    delete(swap_var)
    delete(vprofile)
    delete(vprofile3d)
  end if
  if (randominit .eq. 1) then
   do n=0, nb-1
     print("Bubble number "+n)
     if (randomtemp .eq. 1) then
      waiting = True
      do while (waiting)
        v1 = 2.*random_uniform(0,1,1)-1.
        v2 = 2.*random_uniform(0,1,1)-1.
        rsq = v1^2.+v2^2.
        if (rsq .ge. 1. .or. rsq .eq. 0.) then
           waiting = True
        else
           waiting = False
        end if
      end do
       fac = sqrt( -2.*log(rsq)/rsq )
       gaussdev = v2*fac
       temppert = stdt*gaussdev
       vals = (/delt+temppert,lowt/)
       deltemp = max(vals)
       vals = (/deltemp,hight/)
       deltemp = min(vals)
       delete(vals)
     else
       deltemp = delt
     end if
     print("adding temperature perturbation "+deltemp)
     unidev = random_uniform(0,1,1)
     ib = floattointeger(ibmin+unidev*(ibmax-ibmin))
     unidev = random_uniform(0,1,1)
     jb = floattointeger(jbmin+unidev*(jbmax-jbmin))
     unidev = random_uniform(0,1,1)
     kb = floattointeger(kbmin+unidev*(kbmax-kbmin))
     print("ib, jb, kb for perturbation "+ib+", "+jb+", "+kb)
     swap_var = in_template_file->T
;
     do i = 0, nxm-1
       do j = 0, nym-1
         do k = 0, nzm-1
           dh = dxkm^2. * (i-ib)^2. + dykm^2. * (j-jb)^2.
           dv = dzkm^2. * (k-kb)^2.
           wgt = (1.0 - dh/rbubh^2. - dh/rbubh^2. - dv/rbubv^2.) \
               / (1.0 + dh/rbubh^2. + dh/rbubh^2. + dv/rbubv^2.)
           if (wgt .gt. 0.0) then
              swap_var(0,k,j,i)=swap_var(0,k,j,i)+deltemp*wgt
           end if
         end do
       end do
     end do
     in_template_file->T = (/swap_var(:,:,:,:)/)
     delete(swap_var)
   end do
  end if ; randominit == 1
  if (randominit .eq. 2) then
; Gaussian locations for seeds
    hxboundarylim = floattointeger(ceil(rbubh/dxkm) + 1)
    hyboundarylim = floattointeger(ceil(rbubh/dykm) + 1)
    print("=====================================")
    do ns = 0, numseed-1
      print("Seed number "+ns)
      print("Gaussian locations for seed i, j: "+ibval(ns)+" "+jbval(ns))
       do n=0, nb-1
         print("Bubble number "+n)
         waiting = True
         do while (waiting)
           v1 = 2.*random_uniform(0,1,1)-1.
           v2 = 2.*random_uniform(0,1,1)-1.
           rsq = v1^2.+v2^2.
           if (rsq .ge. 1. .or. rsq .eq. 0.) then
             waiting = True
           else
             waiting = False
           end if
         end do
         fac = sqrt( -2.*log(rsq)/rsq )
         gaussdev = v2*fac
         distperti = floattointeger(stddist*gaussdev)
         gaussdev = v1*fac
         distpertj = floattointeger(stddist*gaussdev)
         if (randomtemp.eq.1) then
           waiting = True
           do while (waiting)
             v1 = 2.*random_uniform(0,1,1)-1.
             v2 = 2.*random_uniform(0,1,1)-1.
             rsq = v1^2.+v2^2.
             if (rsq .ge. 1. .or. rsq .eq. 0.) then
               waiting = True
             else
               waiting = False
             end if
           end do
           fac = sqrt( -2.*log(rsq)/rsq )
           gaussdev = v2*fac
           temppert = stdt*gaussdev
           vals = (/delt+temppert,lowt/)
           deltemp = max(vals)
           vals = (/deltemp,hight/)
           deltemp = min(vals)
           delete(vals)
         else
            deltemp = delt
         end if
         print("adding temperature perturbation "+deltemp)
         ib = ibval(ns) + distperti
         vals = (/ib,hxboundarylim/)
         printVarSummary(vals)
         printVarSummary(ib)
         ib = max(vals)
         vals = (/ib,nxm-hxboundarylim/)
         ib = min(vals)
         jb = jbval(ns) + distpertj
         vals = (/jb,hyboundarylim/)
         jb = max(vals)
         vals = (/jb,nym-hyboundarylim/)
         jb = min(vals)
         delete(vals)
         unidev = random_uniform(0,1,1)
         kb = floattointeger(kbmin + unidev*(kbmax-kbmin))
         print("ib, jb, kb for perturbation "+ib+" "+jb+" "+kb)
         swap_var = in_template_file->T
         do i = 0, nxm-1
           do j = 0, nym-1
             do k = 0, nzm-1
               dh = dxkm^2. * (i-ib)^2. + dykm^2. * (j-jb)^2.
               dv = dzkm^2. * (k-kb)^2.
               wgt = (1.0 - dh/rbubh^2. - dh/rbubh^2. - dv/rbubv^2.) \
                   / (1.0 + dh/rbubh^2. + dh/rbubh^2. + dv/rbubv^2.)
               if (wgt .gt. 0.0) then
                 swap_var(0,k,j,i)=swap_var(0,k,j,i)+deltemp*wgt
               end if
             end do
           end do
         end do
         in_template_file->T = (/swap_var(:,:,:,:)/)
         delete(swap_var)
       end do
    end do
    print("=====================================")
  end if
  if (randominit .eq. 3) then
; File specified locations for seeds
; Test read of a specified bubble location file:
; the format is as follows: bubble mag (deltemp) [cr] bubble loc i  j  k
  bubfil = asciiread(bubble_dir+"/bubble_locations"+member_num,-1,"float")
; figure out the number of values in the array, and sort
  nv = dimsizes(bubfil) / 4      ; number of bubbles in file
  deltemp = new(nv,float,-999.)
  ib      = new(nv,integer,-999)
  jb      = new(nv,integer,-999)
  kb      = new(nv,integer,-999)
  do i=0,nv-1
    deltemp(i) = bubfil(i*4)
    ib(i)      = floattointeger(bubfil(i*4+1))
    jb(i)      = floattointeger(bubfil(i*4+2))
    kb(i)      = floattointeger(bubfil(i*4+3))
  end do
  delete(bubfil)
; Have bubble info - step through numseeds and 
; add in bubbles
;      do ns = 0, numseed-1       ; I don't see why this loop is here - not using
;       print("Seed number "+ns)
    if (nb .gt. nv) then
       print("too many bubbles requested: nb = "+nb+" nv = "+nv)
       exit
    end if
    do n=0, nb-1
      print("Bubble number "+n)
      print("adding temperature perturbation "+deltemp)
      print("ib, jb, kb for perturbation "+ib+", "+jb+", "+kb)
      swap_var = in_template_file->T
      do i = 0, nxm-1
        do j = 0, nym-1
          do k = 0, nzm-1
            dh = dxkm^2. * (i-ib(n))^2. + dykm^2. * (j-jb(n))^2.
            dv = dzkm^2. * (k-kb(n))^2.
            wgt = (1.0 - dh/rbubh^2. - dh/rbubh^2. - dv/rbubv^2.) \
                / (1.0 + dh/rbubh^2. + dh/rbubh^2. + dv/rbubv^2.)
            if (wgt .gt. 0.0) then
              swap_var(0,k,j,i)=swap_var(0,k,j,i)+deltemp*wgt
            end if 
          end do
        end do
      end do
      in_template_file->T = (/swap_var(:,:,:,:)/)
      delete(swap_var)
    end do ; end step thorugh bubbles
;     end do
  end if
; if we added bubbles here (randominit other than 0) need to update the hydrostatic
; balance for the thermal bubbles changes. To do this, need modified temp field, water 
; vapor mixing ratio, base state dry air mass in column, perturbation dry air mass
; in column, perturbation pressure, base state pressure
;
  if (randominit .ne. 0) then   
   qvap = in_template_file->QVAPOR
   thet = in_template_file->T
; If v3 - get from wrfinput, if v2 calculate
   if (wrfver .eq. "V3") then
     pres = in_template_file->P
     pbas = in_template_file->PB
   end if
;
   altt = new((/nzm,nym,nxm/),float,-999.)
   alll = new((/nzm,nym,nxm/),float,-999.)
   alba = new((/nzm,nym,nxm/),float,-999.)
   rdnw = in_template_file->RDNW
   muba = in_template_file->MUB
   muuu = in_template_file->MU
   tini = in_template_file->T_INIT
   swap_var = in_template_file->PH
;
   if (wrfver .eq. "V2") then
    pres = new((/1,nzm,nym,nxm/),float,-999.)
    pbas = new((/1,nzm,nym,nxm/),float,-999.)
    rdnn = in_template_file->RDN
    do i = 0, nxm-1
      do j = 0, nym-1
        k = nzm-1
        qvf1 = 0.5*(qvap(0,k,j,i)+qvap(0,k,j,i))
        qvf2 = 1./(1.+qvf1)
        qvf1 = qvf1*qvf2
        pres(0,k,j,i) = - 0.5*(muuu(0,j,i)+qvf1*muba(0,j,i))/rdnw(0,k)/qvf2
        do k = nzm-2,0   ; step from top down
          qvf1 = 0.5*(qvap(0,k,j,i)+qvap(0,k+1,j,i))
          qvf2 = 1./(1.+qvf1)
          qvf1 = qvf1*qvf2
          pres(0,k,j,i) = pres(0,k+1,j,i) - (muuu(0,j,i) + qvf1*muba(0,j,i))  \
                          /qvf2/rdnn(k+1)
        end do
      end do
    end do
    delete(rdnn)
   end if   ; V2 diff
;
   do i = 0, nxm-1
     do j = 0, nym-1
       do k = 0, nzm-1
        qvf = 1. + rvovrd*qvap(0,k,j,i)
        alba(k,j,i) = (r_d/p1000mb)*(tini(0,k,j,i)+t0)*             \
                      (pbas(0,k,j,i)/p1000mb)^cvpm
        altt(k,j,i) = (r_d/p1000mb)*(thet(0,k,j,i)+t0)*qvf*         \ 
                    (((pres(0,k,j,i)+pbas(0,k,j,i))/p1000mb)^cvpm)
        alll(k,j,i) = altt(k,j,i) - alba(k,j,i)
       end do
       do k = 1, nzm  ; skip lowest boundary
;        prev = swap_var(0,k,j,i)
        swap_var(0,k,j,i) = swap_var(0,k-1,j,i) - (1./rdnw(0,k-1))*(  \ 
                     (muba(0,j,i)+muuu(0,j,i))*alll(k-1,j,i)+         \ 
                      muuu(0,j,i)*alba(k-1,j,i)  ) 
;       part1 = (1./rdnw(0,k-1))*( (muba(0,j,i)+muuu(0,j,i))*alll(k-1,j,i) )
;       part2 = (1./rdnw(0,k-1))*( muuu(0,j,i)*alba(k-1,j,i) )
;       print("k  "+k+" part1 "+part1+"   part2 "+part2+" sum "+  \
;             (part1+part2)+" prev "+prev+"  new "+ swap_var(0,k,j,i))
       end do
     end do
   end do
   in_template_file->PH = (/swap_var(:,:,:,:)/)
   delete(swap_var)
   delete(muuu)
   delete(muba)
   delete(rdnw)
   delete(alba)
   delete(alll)
   delete(altt)
   delete(pbas)
   delete(pres)
   delete(thet)
   delete(qvap)
  end if
; if you need to add REFL_10CM and FALL_SPD_Z_WEIGHTED to the wrfinput file, flag morrison to '1'
  if (morrison .eq. 1) then
; grab h_diabatic as a template, which should be 'clean' of values, but overwrite with zeroes just in case
   swap_var = in_template_file->H_DIABATIC
   swap_var(:,:,:,:) = 0.0
   swap_var@description = "Reflectivity-weighted fallspeed (lamda = 10 cm)"
   swap_var@units       = "m s-1"
   swap_var@stagger     = " "
   in_template_file->FALL_SPD_Z_WEIGHTED = swap_var
   swap_var@description = "Radar reflectivity (lamda = 10 cm)"
   swap_var@units       = "dBZ"
   swap_var@stagger     = " "
   in_template_file->REFL_10CM           = swap_var 
   delete(swap_var)
  end if
end

 
; <next few lines under version control, do not edit>
; $URL$
; $Revision$
; $Date$
