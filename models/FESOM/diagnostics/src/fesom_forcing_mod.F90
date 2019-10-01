!
! Includes tools to process forcing files
! Tools are used to produce the analysis in:
! Aydoğdu, A., Pinardi, N., Özsoy, E., Danabasoglu, G., Gürses, Ö., and Karspeck, A.: Circulation of the Turkish
! Straits System under interannual atmospheric forcing, Ocean Sci., 14, 999-1019, https://doi.org/10.5194/os-14-999-2018, 2018b.
!
! Provided by ali.aydogdu@cmcc.it
!
module fesom_forcing_mod

  use g_config,   only : runyear, day2ext, resultpath, runid, save_count, level_number, thalweg_directory, endday
  use o_param,    only : km3yr2m3sec, rad, vcpw, g, marmax_lon, marmin_lon, marmax_lat, marmin_lat, rho0r
  use o_elements, only : marm_area, voltriangle, bafux_2d, bafuy_2d, elem2D_nodes
  use o_mesh,     only : coord_nod3d, coord_nod2d, nod3d_below_nod2d
  use o_array,    only : tracer, stress_x, stress_y, uf
  use g_parfe,    only : myDim_nod2d, myDim_nod3d, mype, myDim_elem2D
  use utilities,  only : r4, r8, i4, i8, inside

  implicit none

  public  :: compute_wind_stress_curl, & !! computes the wind stress curl from forcing.nc
             compute_surface_buoyancy, & !! computes the buoyancy from forcig 
             compute_wind_work, &        !! computes wind work
             compute_forcing_monthly_timeseries, & !! computes monthly timeseries for a given variable
             forcing_array_setup, &                !! sets arrays for forcing variables
             read_forcing_input, &                 !! reads forcing file
             cal_nodal_alpha_beta                  !! computes thermal expansion and saline contraction coefficients

  real(r8)    :: Ce_atm_oce=1.75e-3 !! exchange coeff. of latent heat over open water
  real(r8)    :: Ch_atm_oce=1.75e-3 !! exchange coeff. of sensible heat over open water
  real(r8)    :: Cd_atm_oce=1.0e-3  !! drag coefficient between atmosphere and water
  real(r8)    :: rho_air   =1.2     !! air density

  real(r8)    :: Ce_atm_ice=1.75e-3 !! exchange coeff. of latent heat over ice
  real(r8)    :: Ch_atm_ice=1.75e-3 !! exchange coeff. of sensible heat over ice
  real(r8)    :: Cd_atm_ice=1.32e-3 !! drag coefficient between atmosphere and ice

  namelist /forcing_exchange_coeff/ Ce_atm_oce, Ch_atm_oce, Cd_atm_oce, &
       Ce_atm_ice, Ch_atm_ice, Cd_atm_ice

  ! forcing arrays
  real(r8), allocatable, dimension(:)         :: wind_u, wind_v
  real(r8), allocatable, dimension(:)         :: srf_Tair, srf_shum
  real(r8), allocatable, dimension(:)         :: srf_Tdew
  real(r8), allocatable, dimension(:)         :: srf_shortwave, srf_longwave
  real(r8), allocatable, dimension(:)         :: srf_long_w
  real(r8), allocatable, dimension(:)         :: srf_prec_rain, srf_prec_snow
  real(r8), allocatable, dimension(:)         :: srf_runoff, srf_evaporation
  real(r8), allocatable, dimension(:)         :: srf_water_flux, srf_heat_flux
  real(r8), allocatable, dimension(:)         :: srf_salt_virtual, srf_salt_relax
  real(r8), allocatable, dimension(:)         :: wind_stress_x,wind_stress_y
  real(r8), allocatable, dimension(:)         :: ccvr, Pair

  ! shortwave penetration
  real(r8), allocatable, dimension(:)         :: chl, sw_3d

  real(r8), allocatable, dimension(:)         :: thdgr, thdgrsn, flice
  real(r8), allocatable, dimension(:)         :: srf_olat_heat, srf_osen_heat, srf_olwout

  ! drag coefficient Cd_atm_oce and transfer coefficients for evaporation
  ! Ce_atm_oce and sensible heat Ch_atm_oce between atmosphere and ocean
  real(r8), allocatable, dimension(:)	  :: Cd_atm_oce_arr
  real(r8), allocatable, dimension(:)	  :: Ch_atm_oce_arr
  real(r8), allocatable, dimension(:)	  :: Ce_atm_oce_arr

  ! drag coefficient Cd_atm_oce between atmosphere and ice
  ! real(kind=8), allocatable, dimension(:)	  :: Cd_atm_ice_arr
  real, dimension(:), allocatable      :: ustar    ! surface friction velocity       (m/s)
  real, dimension(:), allocatable      :: Bo       ! surface turb buoy. forcing  (m^2/s^3)
  real, dimension(:), allocatable      :: Talpha   ! -d(rho)/ d(pot.temperature)  (kg/m^3/C)
  real, dimension(:), allocatable      :: Sbeta    ! d(rho)/ d(salinity)       (kg/m^3/PSU)


  contains
  !----------------------------------------------------------------------------

  subroutine compute_wind_stress_curl

    character(6)      :: cday,clevel
    character(4)      :: cyear
    character(9)      :: out_folder
    character(120)    :: out_filename
    integer(i8)       :: i,j,layer,elem,row,q
    integer           :: elnodes(3)
    integer           :: totday
    real(r8)      :: udx, udy, u_el(3), v_el(3), usum, vsum
    real(r8)      :: dx(3), dy(3), vol, inv12, inv3
    real(r8),allocatable       :: ws_curl(:,:)
    real(r8),allocatable       :: ws_curl_tot(:)
    real(r8),allocatable       :: stress_xy(:,:), wind_xy(:,:)

    inv12=1.0_8/12.0_8
    inv3=1.0_8/3.0_8
    totday=endday
    out_folder='.'
    layer=level_number
    allocate(stress_xy(myDim_nod2D,2))
    allocate(wind_xy(myDim_nod2D,2))
    allocate(ws_curl(2,myDim_nod2D))
    allocate(ws_curl_tot(myDim_nod2D))
    stress_xy=0
    wind_xy=0
    ws_curl=0

    write (clevel,'(a,i3.3)')'LEV',layer
    do day2ext=1,totday

      write (cday,'(a,i3.3)')'DAY',day2ext
      if(mype==0) write(*,*) 'DAY: ', cday
      if(mype==0) write(*,*)
      call read_forcing_input

      if(mype==0) write(*,*) "forcing data is read"
      if(mype==0) write(*,*)

      do elem=1,myDim_elem2D
          elnodes=elem2D_nodes(:,elem)
          vol=voltriangle(elem)/3
          dx=bafux_2d(:,elem)
          dy=bafuy_2d(:,elem)
          u_el=wind_stress_x(elnodes)
          v_el=wind_stress_y(elnodes)
          usum=sum(dy*u_el)
          vsum=sum(dx*v_el)
          do q=1,3
            row=elnodes(q)
            ws_curl(1,row) = usum
            ws_curl(2,row) = vsum
          enddo
      end do

      do i=1,myDim_nod2D
      stress_xy(i,1)=stress_xy(i,1)+wind_stress_x(i)
      stress_xy(i,2)=stress_xy(i,2)+wind_stress_y(i)
      wind_xy(i,1)=wind_xy(i,1)+wind_u(i)
      wind_xy(i,2)=wind_xy(i,2)+wind_v(i)
      ws_curl_tot(i)=( ws_curl(2,i)-ws_curl(1,i) ) * inv3
      end do

    end do

    out_filename=trim(out_folder)//'/WSTC_'//runid//'_'//runyear//'_ANNUAL_'//clevel//'.asc'
      if(mype==0) write(*,*) "Output Filename: ", out_filename
      if(mype==0) write(*,*)

    open(unit=111,file=out_filename,status='replace',access='append',form='formatted')
    do i=1,myDim_nod2D
      write(111,'(2F10.4, 7F16.9)') &
      coord_nod2D(1,i)/rad, coord_nod2D(2,i)/rad, &
      ws_curl_tot(i)/totday, &
      stress_xy(i,1)/totday, stress_xy(i,2)/totday, &
      sqrt(stress_xy(i,1)**2+stress_xy(i,2)**2)/totday, &
      wind_xy(i,1)/totday, wind_xy(i,2)/totday, &
      sqrt(wind_xy(i,1)**2+wind_xy(i,2)**2)/totday
    end do
    close(111)


  end subroutine compute_wind_stress_curl

  subroutine compute_surface_buoyancy


    character(6)      :: cday,clevel
    character(4)      :: cyear
    character(9)      :: out_folder
    character(120)    :: out_filename
    integer           :: i,j,layer,elem
    integer           :: elnodes(3)
    integer           :: totday, row
    real, dimension(:), allocatable      :: Bo_sum(:)    ! surface buoyancy flux m2/s-3
    real, dimension(:), allocatable      :: HF_sum(:)
    real, dimension(:), allocatable      :: WF_sum(:)

    allocate(Bo_sum(myDim_nod2d))
    allocate(HF_sum(myDim_nod2d))
    allocate(WF_sum(myDim_nod2d))

    Bo_sum  = 0
    HF_sum  = 0
    WF_sum  = 0
    out_folder  = '.'
    totday  = endday
    layer   = 1

    out_filename=trim(out_folder)//'/BUOY_'//runid//'_'//runyear//'_ANNUAL_LEV001.asc'
      if(mype==0) write(*,*) "Output Filename: ", out_filename
      if(mype==0) write(*,*)

    open(unit=131,file=out_filename,status='replace',access='append',form='formatted')

      allocate(ustar(myDim_nod3d))
      allocate(Bo(myDim_nod3d))
      allocate(Talpha(myDim_nod3d))
      allocate(Sbeta(myDim_nod3d))

    do day2ext =1,totday

      call oce_input

      call read_forcing_input


      !-----------------------------------------------------------------------
      !     compute thermal and haline expansion coefficients (without factor of rho).
      !     thermal expansion coefficient without 1/rho factor       (kg/m3/C)
      !           talpha= -d(rho{k,k})/d(T(k))
      !     salt expansion coefficient without 1/rho factor        (kg/m3/PSU)
      !           sbeta = d(rho{k,k})/d(S(k))
      !-----------------------------------------------------------------------

      if(mype==0) write(*,*) "calling cal_nodal_alpha_beta"
      if(mype==0) write(*,*)

      call cal_nodal_alpha_beta

      if(mype==0) write(*,*) "alpha beta are calculated"
      if(mype==0) write(*,*)

      !-----------------------------------------------------------------------
      ! friction velocity, turbulent sfc buoyancy forcing
      ! ustar = sqrt( sqrt( stress_x^2 + stress_y^2 ) / rho ) (m/s)
      ! bo =  -g * ( Talpha*heat_flux/vcpw + Sbeta * salinity*water_flux ) (m^2/s^3)
      !-----------------------------------------------------------------------

      do i=1,myDim_nod2d
         row=nod3d_below_nod2d(1,i)
         ustar(i) = sqrt( sqrt(stress_x(i)**2 + stress_y(i)**2)*rho0r)
         Bo(i)    = -g * (Talpha(row) * srf_heat_flux(i)/vcpw  &   !heat_flux & water_flux: positive up
              + Sbeta(row) * srf_water_flux(i) * tracer(row,2))
         Bo_sum(i) = Bo_sum(i) + Bo(i)
         HF_sum(i) = HF_sum(i) + srf_heat_flux(i)
         WF_sum(i) = WF_sum(i) + srf_water_flux(i)
      enddo
    end do
         deallocate(ustar, Bo, Talpha, Sbeta)
      do i=1,myDim_nod2d
      write(131,'(2F10.4, 3F16.9)')  &
      coord_nod2D(1,i)/rad, coord_nod2D(2,i)/rad, Bo_sum(i)/totday, &
      HF_sum(i)/totday, WF_sum(i)/totday
      end do
    close(131)
    deallocate(Bo_sum, HF_sum, WF_sum)


  end subroutine compute_surface_buoyancy


  subroutine compute_wind_work

    character(9)   :: out_folder
    character(120) :: out_filename

    integer        :: totday, row, layer, i, k, elnodes2(3)

    real(kind=8)   :: inv3, lat_o, lon_o
    real(kind=8)   :: tau_vel(2), wind_work
    real(kind=8)   :: t_stress, wind_stress
    real(kind=8)   :: totvol, area, volmar
    real(kind=8)   :: rho_0 = 1028
    real(kind=8),allocatable       :: stress_xy(:,:)
    real(kind=8),allocatable       :: wind_xy(:,:)
    real(kind=8),allocatable       :: curr_xy(:,:)

    inv3 = 1.0_8/3.0_8

    totday     = endday
    out_folder = '.'
    layer      = level_number

    allocate(stress_xy(myDim_nod2D,2))
    allocate(wind_xy(myDim_nod2D,2))
    allocate(curr_xy(myDim_nod2D,2))

    stress_xy  = 0
    wind_xy    = 0
    curr_xy    = 0

       out_filename=trim(out_folder)//'/WWRK_'//runid//'_'//runyear//'_DAILY.asc'
       open(unit=101,file=out_filename,status='replace',access='append',form='formatted')

  dloop:  do day2ext = 1,totday

    totvol     = 0
    area       = 0
    wind_work  = 0
    wind_stress= 0

      call oce_input
      call read_forcing_input


      do i=1,myDim_nod2d

        row            = nod3d_below_nod2d(1,i)

        wind_xy(i,1)   = wind_u(i)
        wind_xy(i,2)   = wind_v(i)
        curr_xy(i,1)   = uf(row)
        curr_xy(i,2)   = uf(row+myDim_nod3d)

        stress_xy(i,1) = rho_air*Cd_atm_oce * &
                         sqrt(wind_xy(i,1)**2+wind_xy(i,2)**2)*wind_xy(i,1)

        stress_xy(i,2) = rho_air*Cd_atm_oce * &
                         sqrt(wind_xy(i,1)**2+wind_xy(i,2)**2)*wind_xy(i,2)
      end do

  eloop:    do i=1,myDim_elem2D
      elnodes2  = elem2D_nodes(:,i)
      tau_vel   = 0
      t_stress  = 0
      area      = 0

      do k=1,3
        lat_o = coord_nod2D(2,elnodes2(k))/rad
        lon_o = coord_nod2D(1,elnodes2(k))/rad

        if ( lat_o.gt.marmax_lat.or.lon_o.gt.marmax_lon &
           .or.lat_o.lt.marmin_lat.or.lon_o.lt.marmin_lon) then
         go to 1234
        endif
      end do

      do k=1,3
         lat_o = coord_nod2D(2,elnodes2(k))/rad
         lon_o = coord_nod2D(1,elnodes2(k))/rad

         if ( lat_o.gt.41.and.lon_o.gt.28.7 ) go to 1234
      end do

      do k=1,3
          tau_vel(1) = tau_vel(1) + stress_xy(elnodes2(k),1) * curr_xy(elnodes2(k),1)
          tau_vel(2) = tau_vel(2) + stress_xy(elnodes2(k),2) * curr_xy(elnodes2(k),2)
          t_stress   = t_stress + sqrt(stress_xy(elnodes2(k),1)**2+stress_xy(elnodes2(k),2)**2)
      enddo
         totvol    = totvol + voltriangle(i)
         tau_vel   = tau_vel * inv3
         wind_stress = wind_stress + t_stress * voltriangle(i) * inv3
         wind_work   = wind_work + (tau_vel(1) + tau_vel(2)) * voltriangle(i)

        1234 continue
     enddo eloop
     write(101,'(I5,4F30.8)') day2ext, wind_work / ( totvol * rho_0 ), wind_stress / totvol, totvol

    end do dloop


  end subroutine compute_wind_work

  subroutine compute_forcing_monthly_timeseries

    logical           :: check
    character(6)      :: cday
    character(4)      :: cyear
    character(9)      :: out_folder
    character(100)    :: POLYFILE, BONDFILE
    character(126)    :: dout_filename, mout_filename
    integer           :: i, j, k, layer,totalday, start_date, final_date
    integer           :: mmonth(13), myear, nline, inc

    integer           :: elnodes2(3)

    real(r8)          :: inv3=1.0/3.0_8
    real(r8)          :: inv4=1.0/4.0_8
    real(r8)          :: inv12=1.0/12.0_8
    real(r8)          :: nodarea
    real(r8)          :: lat(3), lon(3)
    real(r8), allocatable, dimension(:,:) :: polygon

    real(r8)          :: el_shortwave, el_longwave, el_long_w
    real(r8)          :: el_prec_rain, el_prec_snow, el_runoff
    real(r8)          :: el_Tair, el_Tdew, el_shum, el_evaporation
    real(r8)          :: el_heat_flux, el_water_flux
    real(r8)          :: el_olat, el_osen, el_olwout
    real(r8)          :: el_salt_virtual, el_salt_relax
    real(r8)          :: el_wind_u, el_wind_v
    real(r8)          :: el_wind_stress_x, el_wind_stress_y

    real(r8)          :: n2_shortwave, n2_longwave, n2_long_w
    real(r8)          :: n2_prec_rain, n2_prec_snow, n2_runoff
    real(r8)          :: n2_Tair, n2_Tdew, n2_shum, n2_evaporation
    real(r8)          :: n2_heat_flux, n2_water_flux
    real(r8)          :: n2_olat, n2_osen, n2_olwout
    real(r8)          :: n2_salt_virtual, n2_salt_relax
    real(r8)          :: n2_wind_u, n2_wind_v
    real(r8)          :: n2_wind_stress_x, n2_wind_stress_y

    real(r8)          :: m2_shortwave, m2_longwave, m2_long_w
    real(r8)          :: m2_prec_rain, m2_prec_snow, m2_runoff
    real(r8)          :: m2_Tair, m2_Tdew, m2_shum, m2_evaporation
    real(r8)          :: m2_heat_flux, m2_water_flux
    real(r8)          :: m2_olat, m2_osen, m2_olwout
    real(r8)          :: m2_salt_virtual, m2_salt_relax
    real(r8)          :: m2_wind_u, m2_wind_v
    real(r8)          :: m2_wind_stress_x, m2_wind_stress_y


    BONDFILE='BOUND_MARM'
    POLYFILE=trim(thalweg_directory)//'/'//trim(BONDFILE)//'.lst'

    open(unit=111,file=POLYFILE)
    read(111,*) nline
    print*, 'nline:', nline

    allocate(polygon(2,nline))
    do inc = 1,nline
    read(111,*) polygon(1,inc), polygon(2,inc)
    end do
    close(111)

    if ( runyear == '2008' ) then
            mmonth=(/ 1, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
    else
            mmonth=(/ 1, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
    end if

    out_folder='.'

    layer=level_number

    mout_filename=trim(out_folder)//'/FORC_'//runid//'_'//runyear//'_MNTLY_MARM.asc'
    dout_filename=trim(out_folder)//'/FORC_'//runid//'_'//runyear//'_DAILY_MARM.asc'
    open(unit=101,file=dout_filename,status='replace',access='append',form='formatted')
    open(unit=102,file=mout_filename,status='replace',access='append',form='formatted')


  lpmnt: do j=1,12

       m2_shortwave     =  0
       m2_longwave      =  0
       m2_long_w        =  0
       m2_runoff        =  0
       m2_prec_rain     =  0
       m2_prec_snow     =  0
       m2_wind_u        =  0
       m2_wind_v        =  0
       m2_Tair          =  0
       m2_Tdew          =  0
       m2_shum          =  0
       m2_evaporation   =  0
       m2_heat_flux     =  0
       m2_olat          =  0
       m2_osen          =  0
       m2_olwout        =  0
       m2_water_flux    =  0
       m2_salt_virtual  =  0
       m2_salt_relax    =  0
       m2_wind_stress_x =  0
       m2_wind_stress_y =  0


       start_date = sum(mmonth(1:j))
       final_date = sum(mmonth(2:j+1))
       totalday   = final_date-start_date+1


  lpday: do day2ext = start_date,final_date


       n2_shortwave     =  0
       n2_longwave      =  0
       n2_long_w        =  0
       n2_runoff        =  0
       n2_prec_rain     =  0
       n2_prec_snow     =  0
       n2_wind_u        =  0
       n2_wind_v        =  0
       n2_Tair          =  0
       n2_Tdew          =  0
       n2_shum          =  0
       n2_evaporation   =  0
       n2_heat_flux     =  0
       n2_olat          =  0
       n2_osen          =  0
       n2_olwout        =  0
       n2_water_flux    =  0
       n2_salt_virtual  =  0
       n2_salt_relax    =  0
       n2_wind_stress_x =  0
       n2_wind_stress_y =  0

         write (cday,'(a,i3.3)')'DAY',day2ext

         call read_forcing_input


  lpele: do i=1,myDim_elem2D

         elnodes2=elem2D_nodes(:,i)
         lon = ( coord_nod2D(1,elnodes2(:)) )/rad
         lat = ( coord_nod2D(2,elnodes2(:)) )/rad

          do k=1,3
            check = inside(lon(k), lat(k), polygon(1,:), polygon(2,:), nline)
            if (.not.check) go to 1335
          end do

          el_shortwave     =  0
          el_longwave      =  0
          el_long_w        =  0
          el_runoff        =  0
          el_prec_rain     =  0
          el_prec_snow     =  0
          el_wind_u        =  0
          el_wind_v        =  0
          el_Tair          =  0
          el_Tdew          =  0
          el_shum          =  0
          el_evaporation   =  0
          el_heat_flux     =  0
          el_olat          =  0
          el_osen          =  0
          el_olwout        =  0
          el_water_flux    =  0
          el_salt_virtual  =  0
          el_salt_relax    =  0
          el_wind_stress_x =  0
          el_wind_stress_y =  0

  lpnod:  do k=1,3
           el_shortwave     = el_shortwave     + (srf_shortwave(elnodes2(k)))   * voltriangle(i)
  !         el_longwave      = el_longwave      + (srf_longwave(elnodes2(k)))    * voltriangle(i)
           el_longwave      = el_longwave      + (srf_heat_flux(elnodes2(k)) - &
             ( srf_shortwave(elnodes2(k)) + srf_olat_heat(elnodes2(k)) + srf_osen_heat(elnodes2(k)) ) )  * voltriangle(i)
           el_long_w        = el_long_w        + (srf_long_w(elnodes2(k)))      * voltriangle(i)
           el_runoff        = el_runoff        + (srf_runoff(elnodes2(k)))      * voltriangle(i)
           el_prec_rain     = el_prec_rain     + (srf_prec_rain(elnodes2(k)))   * voltriangle(i)
           el_prec_snow     = el_prec_snow     + (srf_prec_snow(elnodes2(k)))   * voltriangle(i)
           el_wind_u        = el_wind_u        + (wind_u(elnodes2(k)))          * voltriangle(i)
           el_wind_v        = el_wind_v        + (wind_v(elnodes2(k)))          * voltriangle(i)
           el_Tair          = el_Tair          + (srf_Tair(elnodes2(k)))        * voltriangle(i)
           el_Tdew          = el_Tdew          + (srf_Tdew(elnodes2(k)))        * voltriangle(i)
           el_shum          = el_shum          + (srf_shum(elnodes2(k)))        * voltriangle(i)
           el_evaporation   = el_evaporation   + (srf_evaporation(elnodes2(k))) * voltriangle(i)
           el_heat_flux     = el_heat_flux     + (srf_heat_flux(elnodes2(k)))   * voltriangle(i)
           el_olat          = el_olat          + (srf_olat_heat(elnodes2(k)))   * voltriangle(i)
           el_osen          = el_osen          + (srf_osen_heat(elnodes2(k)))   * voltriangle(i)
           el_osen          = el_olwout        + (srf_olwout(elnodes2(k)))      * voltriangle(i)
           el_water_flux    = el_water_flux    + (srf_water_flux(elnodes2(k)))  * voltriangle(i)
           el_salt_virtual  = el_salt_virtual  + (srf_salt_virtual(elnodes2(k)))* voltriangle(i)
           el_salt_relax    = el_salt_relax    + (srf_salt_relax(elnodes2(k)))  * voltriangle(i)
           el_wind_stress_x = el_wind_stress_x + (wind_stress_x(elnodes2(k)))   * voltriangle(i)
           el_wind_stress_y = el_wind_stress_y + (wind_stress_y(elnodes2(k)))   * voltriangle(i)
          end do lpnod
  !         el_water_flux    = el_water_flux    + (srf_runoff(elnodes2(k))+srf_evaporation(elnodes2(k))+srf_prec_rain(elnodes2(k))+srf_prec_snow(elnodes2(k)))  * voltriangle(i)
          n2_shortwave     = n2_shortwave     + el_shortwave     * inv3
          n2_longwave      = n2_longwave      + el_longwave      * inv3
          n2_long_w        = n2_long_w        + el_long_w        * inv3
          n2_runoff        = n2_runoff        + el_runoff        * inv3
          n2_prec_rain     = n2_prec_rain     + el_prec_rain     * inv3
          n2_prec_snow     = n2_prec_snow     + el_prec_snow     * inv3
          n2_wind_u        = n2_wind_u        + el_wind_u        * inv3
          n2_wind_v        = n2_wind_v        + el_wind_v        * inv3
          n2_Tair          = n2_Tair          + el_Tair          * inv3
          n2_Tdew          = n2_Tdew          + el_Tdew          * inv3
          n2_shum          = n2_shum          + el_shum          * inv3
          n2_evaporation   = n2_evaporation   + el_evaporation   * inv3
          n2_heat_flux     = n2_heat_flux     + el_heat_flux     * inv3
          n2_olat          = n2_olat          + el_olat          * inv3
          n2_osen          = n2_osen          + el_osen          * inv3
          n2_olwout        = n2_olwout        + el_olwout        * inv3
          n2_water_flux    = n2_water_flux    + el_water_flux    * inv3
          n2_salt_virtual  = n2_salt_virtual  + el_salt_virtual  * inv3
          n2_salt_relax    = n2_salt_relax    + el_salt_relax    * inv3
          n2_wind_stress_x = n2_wind_stress_x + el_wind_stress_x * inv3
          n2_wind_stress_y = n2_wind_stress_y + el_wind_stress_y * inv3

         1335 continue
        end do lpele

             write(101,'(22F25.9)')  &
                   n2_shortwave     / marm_area, &
                   n2_longwave      / marm_area, &
                   n2_long_w        / marm_area, &
                   n2_runoff      / km3yr2m3sec, &
                   n2_prec_rain   / km3yr2m3sec, &
                   n2_prec_snow   / km3yr2m3sec, &
                   n2_wind_u        / marm_area, &
                   n2_wind_v        / marm_area, &
                   n2_Tair          / marm_area, &
                   n2_Tdew          / marm_area, &
                   n2_shum          / marm_area, &
                   n2_evaporation / km3yr2m3sec, &
                   n2_heat_flux     / marm_area, &
                   n2_olat          / marm_area, &
                   n2_osen          / marm_area, &
                   n2_olwout        / marm_area, &
                   n2_water_flux  / km3yr2m3sec, &
                   n2_salt_virtual  / marm_area, &
                   n2_salt_relax    / marm_area, &
                   n2_wind_stress_x / marm_area, &
                   n2_wind_stress_y / marm_area, &
                   marm_area

          m2_shortwave     = m2_shortwave     + n2_shortwave     / marm_area
          m2_longwave      = m2_longwave      + n2_longwave      / marm_area
          m2_long_w        = m2_long_w        + n2_long_w        / marm_area
          m2_runoff        = m2_runoff        + n2_runoff      / km3yr2m3sec
          m2_prec_rain     = m2_prec_rain     + n2_prec_rain   / km3yr2m3sec
          m2_prec_snow     = m2_prec_snow     + n2_prec_snow   / km3yr2m3sec
          m2_wind_u        = m2_wind_u        + n2_wind_u        / marm_area
          m2_wind_v        = m2_wind_v        + n2_wind_v        / marm_area
          m2_Tair          = m2_Tair          + n2_Tair          / marm_area
          m2_Tdew          = m2_Tdew          + n2_Tdew          / marm_area
          m2_shum          = m2_shum          + n2_shum          / marm_area
          m2_evaporation   = m2_evaporation   + n2_evaporation / km3yr2m3sec
          m2_heat_flux     = m2_heat_flux     + n2_heat_flux     / marm_area
          m2_olat          = m2_olat          + n2_olat          / marm_area
          m2_osen          = m2_osen          + n2_osen          / marm_area
          m2_olwout        = m2_olwout        + n2_olwout        / marm_area
          m2_water_flux    = m2_water_flux    + n2_water_flux  / km3yr2m3sec
          m2_salt_virtual  = m2_salt_virtual  + n2_salt_virtual  / marm_area
          m2_salt_relax    = m2_salt_relax    + n2_salt_relax    / marm_area
          m2_wind_stress_x = m2_wind_stress_x + n2_wind_stress_x / marm_area
          m2_wind_stress_y = m2_wind_stress_y + n2_wind_stress_y / marm_area


      end do lpday

             write(102,'(22F25.9)')  &
                   m2_shortwave     / totalday, &
                   m2_longwave      / totalday, &
                   m2_long_w        / totalday, &
                   m2_runoff        / totalday, &
                   m2_prec_rain     / totalday, &
                   m2_prec_snow     / totalday, &
                   m2_wind_u        / totalday, &
                   m2_wind_v        / totalday, &
                   m2_Tair          / totalday, &
                   m2_Tdew          / totalday, &
                   m2_shum          / totalday, &
                   m2_evaporation   / totalday, &
                   m2_heat_flux     / totalday, &
                   m2_olat          / totalday, &
                   m2_osen          / totalday, &
                   m2_olwout        / totalday, &
                   m2_water_flux    / totalday, &
                   m2_salt_virtual  / totalday, &
                   m2_salt_relax    / totalday, &
                   m2_wind_stress_x / totalday, &
                   m2_wind_stress_y / totalday, &
                   marm_area

     end do lpmnt
     close(101)
     close(102)
  end subroutine compute_forcing_monthly_timeseries

  subroutine forcing_array_setup

    integer    :: n2

    n2=myDim_nod2D
    ! Allocate memory for atmospheric forcing
    allocate(srf_shortwave(n2), srf_longwave(n2))
    allocate(srf_long_w(n2))
    allocate(srf_prec_rain(n2), srf_prec_snow(n2))
    allocate(wind_u(n2), wind_v(n2))
    allocate(srf_Tair(n2), srf_shum(n2))
    allocate(srf_Tdew(n2), srf_heat_flux(n2))
    allocate(srf_runoff(n2), srf_evaporation(n2), srf_water_flux(n2))
    allocate(srf_salt_virtual(n2), srf_salt_relax(n2))
    allocate( wind_stress_x(n2), wind_stress_y(n2) )


    srf_shortwave=0.
    srf_longwave=0.
    srf_long_w=0.
    srf_prec_rain=0.
    srf_prec_snow=0.
    wind_u=0.
    wind_v=0.
    srf_Tair=0.
    srf_Tdew=0.
    srf_shum=0.
    srf_runoff=0.
    srf_evaporation=0
    srf_heat_flux=0
    srf_water_flux=0
    srf_salt_virtual=0
    srf_salt_relax=0
    wind_stress_x=0
    wind_stress_y=0

   ! if(rad_data_source=='MFS') then
    allocate(ccvr(n2), Pair(n2))
    ccvr=0.
    Pair=0.
  !  end if


    ! shortwave penetration
  #ifdef use_sw_pene
    allocate(chl(n2))
    allocate(sw_3d(myDim_nod3d))
    chl=0.0
  #endif

    !for ice diagnose
  !#ifdef use_ice
    allocate(thdgr(n2), thdgrsn(n2), flice(n2))
    allocate(srf_olat_heat(n2), srf_osen_heat(n2), srf_olwout(n2))
    thdgr=0.
    thdgrsn=0.
    flice=0.
    srf_olat_heat=0.
    srf_osen_heat=0.
    srf_olwout=0.
  !#endif

    ! drag coefficient and transfer coefficients for latent and sensible heat
    allocate(Cd_atm_oce_arr(n2))
    allocate(Ce_atm_oce_arr(n2))
    allocate(Ch_atm_oce_arr(n2))
    Cd_atm_oce_arr=Cd_atm_oce
    Ce_atm_oce_arr=Ce_atm_oce
    Ch_atm_oce_arr=Ch_atm_oce

    if(mype==0) write(*,*) 'forcing arrays have been set up'
    if(mype==0) write(*,*)

  end subroutine forcing_array_setup

  subroutine read_forcing_input

  #include "netcdf.inc"

    real(kind=8), allocatable :: aux2(:)
    integer                   :: istart(2), icount(2), n3
    integer                   :: status, ncid, dimid_rec, j, nrec
    integer                   :: dimid_2d, dimid_3d, dimids(2)
    integer                   :: time_varid, iter_varid
    integer                   :: thdgr_varid, thdgrsn_varid
    integer                   :: tair_varid, shum_varid, uwind_varid, vwind_varid
    integer                   :: rain_varid, snow_varid, runoff_varid
    integer                   :: evap_varid, lwrd_varid, swrd_varid
    integer                   :: qnet_varid, wnet_varid
    integer                   :: olat_varid, osen_varid, olwout_varid
    integer                   :: virtual_salt_varid, relax_salt_varid
    integer                   :: stress_x_varid, stress_y_varid
    character(100)            :: longname
    character(100)            :: filename
    character(1)              :: trind


    allocate(aux2(myDim_nod2D))
    nrec=day2ext
    if(mype==0) write(*,*) 'nrec: ', nrec
    if(mype==0) write(*,*)

    ! open files
    filename=trim(ResultPath)//runid//'.'//runyear//'.forcing.diag.nc'

    if(mype==0) write(*,*) 'Input filename: ',filename
    if(mype==0) write(*,*)

    status = nf_open(filename, nf_nowrite, ncid)
    if (status .ne. nf_noerr) call handle_err(status)
    ! inquire variable id

    status = nf_inq_varid(ncid, 'tair', tair_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_inq_varid(ncid, 'shum', shum_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_inq_varid(ncid, 'uwind', uwind_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_inq_varid(ncid, 'vwind', vwind_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_inq_varid(ncid, 'rain', rain_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_inq_varid(ncid, 'snow', snow_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_inq_varid(ncid, 'runoff', runoff_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_inq_varid(ncid, 'evap', evap_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_inq_varid(ncid, 'lwrd', lwrd_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_inq_varid(ncid, 'swrd', swrd_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_inq_varid(ncid, 'qnet', qnet_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_inq_varid(ncid, 'olat', olat_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_inq_varid(ncid, 'osen', osen_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_inq_varid(ncid, 'olwout', olwout_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_inq_varid(ncid, 'wnet', wnet_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_inq_varid(ncid, 'virtual_salt', virtual_salt_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_inq_varid(ncid, 'relax_salt', relax_salt_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_inq_varid(ncid, 'stress_x', stress_x_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_inq_varid(ncid, 'stress_y', stress_y_varid)
    if (status .ne. nf_noerr) call handle_err(status)

       ! Define the netCDF variables for 2D fields.
       ! In Fortran, the unlimited dimension must come
       ! last on the list of dimids.

    ! read variables

    istart=(/1,nrec/)
    icount=(/myDim_nod2D, 1/)
    status = nf_get_vara_double(ncid, tair_varid, istart, icount, aux2)
    if (status .ne. nf_noerr) call handle_err(status)
    srf_Tair=aux2(1:myDim_nod2D)
    status = nf_get_vara_double(ncid, shum_varid, istart, icount, aux2)
    if (status .ne. nf_noerr) call handle_err(status)
    srf_shum=aux2(1:myDim_nod2D)
    status = nf_get_vara_double(ncid, uwind_varid, istart, icount, aux2)
    if (status .ne. nf_noerr) call handle_err(status)
    wind_u=aux2(1:myDim_nod2D)
    status = nf_get_vara_double(ncid, vwind_varid, istart, icount, aux2)
    if (status .ne. nf_noerr) call handle_err(status)
    wind_v=aux2(1:myDim_nod2D)
    status = nf_get_vara_double(ncid, rain_varid, istart, icount, aux2)
    if (status .ne. nf_noerr) call handle_err(status)
    srf_prec_rain=aux2(1:myDim_nod2D)
    status = nf_get_vara_double(ncid, snow_varid, istart, icount, aux2)
    if (status .ne. nf_noerr) call handle_err(status)
    srf_prec_snow=aux2(1:myDim_nod2D)
    status = nf_get_vara_double(ncid, runoff_varid, istart, icount, aux2)
    if (status .ne. nf_noerr) call handle_err(status)
    srf_runoff=aux2(1:myDim_nod2D)
    status = nf_get_vara_double(ncid, evap_varid, istart, icount, aux2)
    if (status .ne. nf_noerr) call handle_err(status)
    srf_evaporation=aux2(1:myDim_nod2D)
    status = nf_get_vara_double(ncid, lwrd_varid, istart, icount, aux2)
    if (status .ne. nf_noerr) call handle_err(status)
    srf_longwave=aux2(1:myDim_nod2D)
    status = nf_get_vara_double(ncid, swrd_varid, istart, icount, aux2)
    if (status .ne. nf_noerr) call handle_err(status)
    srf_shortwave=aux2(1:myDim_nod2D)
    status = nf_get_vara_double(ncid, qnet_varid, istart, icount, aux2)
    if (status .ne. nf_noerr) call handle_err(status)
    srf_heat_flux=aux2(1:myDim_nod2D)
    status = nf_get_vara_double(ncid, olat_varid, istart, icount, aux2)
    if (status .ne. nf_noerr) call handle_err(status)
    srf_olat_heat=aux2(1:myDim_nod2D)
    status = nf_get_vara_double(ncid, osen_varid, istart, icount, aux2)
    if (status .ne. nf_noerr) call handle_err(status)
    srf_osen_heat=aux2(1:myDim_nod2D)
    status = nf_get_vara_double(ncid, olwout_varid, istart, icount, aux2)
    if (status .ne. nf_noerr) call handle_err(status)
    srf_olwout=aux2(1:myDim_nod2D)
    status = nf_get_vara_double(ncid, wnet_varid, istart, icount, aux2)
    if (status .ne. nf_noerr) call handle_err(status)
    srf_water_flux=aux2(1:myDim_nod2D)
    status = nf_get_vara_double(ncid, virtual_salt_varid, istart, icount, aux2)
    if (status .ne. nf_noerr) call handle_err(status)
    srf_salt_virtual=aux2(1:myDim_nod2D)
    status = nf_get_vara_double(ncid, relax_salt_varid, istart, icount, aux2)
    if (status .ne. nf_noerr) call handle_err(status)
    srf_salt_relax=aux2(1:myDim_nod2D)
    status = nf_get_vara_double(ncid, stress_x_varid, istart, icount, aux2)
    if (status .ne. nf_noerr) call handle_err(status)
    wind_stress_x=aux2(1:myDim_nod2D)
    status = nf_get_vara_double(ncid, stress_y_varid, istart, icount, aux2)
    if (status .ne. nf_noerr) call handle_err(status)
    wind_stress_y=aux2(1:myDim_nod2D)

    status=nf_close(ncid)
    if (status .ne. nf_noerr) call handle_err(status)

    if(mype==0) write(*,*) "forcing variables are read"
    if(mype==0) write(*,*)

    ! the next record to be saved
    save_count=nrec+1

    deallocate(aux2)

  end subroutine read_forcing_input


  subroutine cal_nodal_alpha_beta



      !   A function to calculate the thermal expansion coefficient
      !   and saline contraction coefficient.
      !
      ! REFERENCE:
      !    McDougall, T.J. 1987.  Neutral Surfaces
      !    Journal of Physical Oceanography, vol 17, 1950-1964,
      !
      ! INPUT:
      !   tracer(:,1) = potential temperature [degree C (ITS-90)]
      !   tracer(:,2) = salinity              [psu      (PSS-78)]
      !   z           = pressure (or -depth)  [db]
      !
      ! OUTPUT:
      !   Talpha = Thermal expansion coeff (alpha) [degree_C.^-1]
      !   Sbeta  = Saline contraction coeff (beta) [psu.^-1]
      !
      ! Qiang Wang, 25,11,2004
      ! Modified to compute nodal values for KPP, Feb. 2011, Qiang
      !-----------------------------------------------------------------
      ! CHECK VALUE:
      !    sw_beta=0.72088e-3 psu^-1 @ S=40.0psu, ptmp=10.0C (ITS-90), p=4000db
      !    a_over_b=0.34765 psu*C^-1 @ S=40.0psu, ptmp=10.0C, p=4000db
      !-----------------------------------------------------------------

      integer    :: i, ni, row
      real       :: t1,t1_2,t1_3,t1_4,p1,p1_2,p1_3,s1,s35,s35_2
      real       :: a_over_b

      ! we should ensure that we are evoluting potential temperature in the model.

      ! cycle over node
      do i = 1, myDim_nod3d
         ! prepare values that will be used below
         t1 = tracer(i,1)*1.00024_8
         s1 = tracer(i,2)
         p1 = abs(coord_nod3D(3,i))

         t1_2 = t1*t1
         t1_3 = t1_2*t1
         t1_4 = t1_3*t1
         p1_2 = p1*p1
         p1_3 = p1_2*p1
         s35 = s1-35.0_8
         s35_2 = s35*s35

         ! calculate beta
         Sbeta(i) = 0.785567e-3 - 0.301985e-5*t1 &
              + 0.555579e-7*t1_2 - 0.415613e-9*t1_3 &
              + s35*(-0.356603e-6 + 0.788212e-8*t1 &
              + 0.408195e-10*p1 - 0.602281e-15*p1_2) &
              + s35_2*(0.515032e-8) &
              + p1*(-0.121555e-7 + 0.192867e-9*t1 - 0.213127e-11*t1_2) &
              + p1_2*(0.176621e-12 - 0.175379e-14*t1) &
              + p1_3*(0.121551e-17)

         ! calculate the thermal expansion / saline contraction ratio
         a_over_b = 0.665157e-1 + 0.170907e-1*t1 &
              - 0.203814e-3*t1_2 + 0.298357e-5*t1_3 &
              - 0.255019e-7*t1_4 &
              + s35*(0.378110e-2 - 0.846960e-4*t1 &
              - 0.164759e-6*p1 - 0.251520e-11*p1_2) &
              + s35_2*(-0.678662e-5) &
              + p1*(0.380374e-4 - 0.933746e-6*t1 + 0.791325e-8*t1_2) &
              + p1_2*t1_2*(0.512857e-12) &
              - p1_3*(0.302285e-13)

         ! calculate alpha
         Talpha(i) = a_over_b*Sbeta(i)
      end do
  end subroutine cal_nodal_alpha_beta

end module fesom_forcing_mod
