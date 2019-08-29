
! FESOM 2 (Finite-volumE Sea ice-Ocean Model)
! multi-resolution ocean general circulation model
! FESOM/fesom2 is licensed under the GNU General Public License v2.0
! Copyright (C) 2018  FESOM team
!
! This source file was taken from  the FESOM V1.4 modules


subroutine read_namelist
  ! Routine reads namelist files to overwrite default parameters.

  use g_config
  use o_param
  use g_parfe
  use g_clock, only: timenew, daynew, yearnew
  implicit none

  character(len=100)   :: nmlfile
  namelist /clockinit/ timenew, daynew, yearnew

  nmlfile ='namelist.config'    ! name of general configuration namelist file
  open (20,file=nmlfile)
  read (20,NML=modelname)
  read (20,NML=dart_nml)
  read (20,NML=timestep)
  read (20,NML=timeseries)
  read (20,NML=clockinit) 
  read (20,NML=paths)
  read (20,NML=initialization)  
  read (20,NML=inout)
  read (20,NML=mesh_def)
  read (20,NML=geometry)
  read (20,NML=calendar)
  read (20,NML=postproc)
  close (20)
  print*, "namelist is read"
end subroutine read_namelist


subroutine ocean_array_setup
  ! Sets up the arrays needed by the ocean part
  use o_param
  use o_array
  use o_mesh
  use o_read
  use o_elements
  use g_config
  use g_parfe
  integer         :: size3D, size2D         

  size3d=myDim_nod3D
  size2d=myDim_nod2D
  num_tracer=2

  ! density and pressure arrays
  allocate(density_insitu(size3D), density_ref(size3D)) 
  density_insitu=0.
  density_ref=0.
  if(grid_type/=2) then
     allocate(hpressure(size3D))
     hpressure=0.
  end if
  if(grid_type/=1) then
     allocate(PGF(2,max_num_layers-1,myDim_elem2D))
     PGF=0.
!     call init_pressure_force
  end if

  ! Arrays used for ocean boundary forcing
  allocate(stress_x(size2d))
  allocate(stress_y(size2d)) 
  allocate(heat_flux(size2d)) 
  allocate(water_flux(size2d))
  allocate(Tsurf(size2d))
  allocate(Ssurf(size2d)) 
  allocate(ts_sfc_force(size2d, 2))
  allocate(uv_sfc_force(size2d, 2))
  allocate(uv_bott_force(size2d, 2))
  stress_x=0.
  stress_y=0.
  heat_flux=0.
  water_flux=0.
  ts_sfc_force=0.0
  uv_sfc_force=0.0
  uv_bott_force=0.0

  ! T, S fields, their increments and rhs
  allocate(tracer(size3d,num_tracer))     
  allocate(dtracer(size3d,num_tracer))   
  allocate(tracer_rhs(size3d,num_tracer)) 
  allocate(tracer0(size3d,num_tracer))
  tracer=0.0
  tracer0=0.0
  dtracer=0.0
  tracer_rhs=0.0

  ! u, v fields and their rhs 
#ifndef use_non_hydrostatic
  allocate(uf(2*size3D))             
  allocate(uf0(2*size3D))
  allocate(duf(2*size3D))           
  allocate(uv_rhs(2*size3D))
#else
  allocate(nhp(size3D), nhp0(size3D), nhp_rhs(size3D))   
  allocate(uf(3*size3D))          
  allocate(uf0(3*size3D))
  allocate(duf(3*size3D))           
  allocate(uv_rhs(3*size3D)) 
  nhp=0.
  nhp0=0.
#endif
  uf=0.
  uf0=0.
  duf=0.
  uv_rhs=0.

  !ssh
  allocate(ssh(size2D), ssh0(size2D), dssh(size2d), ssh_rhs(size2D)) 
  print*,"sizes: ",size2D, size3D
  ssh=0.
  ssh0=0.
  dssh=0.
  uf=0.       ! ::OG::  not needed
  uf0=0.      ! ::OG::  not needed

  ! arrays for the AB2 coriolis case
  if(.not.use_cori_semi) then
     allocate(ucori(size3d), vcori(size3d))
     allocate(ucori_back(size3D), vcori_back(size3D))
     ucori=0.
     vcori=0.
     ucori_back=0.
     vcori_back=0.
  endif

  ! rhs of w equation and w-potential field
#ifndef use_non_hydrostatic
  allocate(wrhs(size3D),w(size3D))   
  wrhs=0.
  w=0.
#endif

  ! arrays for salt fluxes
  allocate(virtual_salt(size2d), relax_salt(size2d))
  virtual_salt=0.
  relax_salt=0.
#ifdef use_fullfreesurf
  allocate(real_salt_flux(size2d))
  real_salt_flux=0.
#endif  

  ! Redi/GM
  if (Redi_GM) then    
     if(nslope_version==1) then 
        allocate(neutral_slope(3,max_num_layers-1,myDim_elem2d))
        neutral_slope=0.0
     else
        allocate(neutral_slope_elem(3,myDim_elem3d))
        neutral_slope_elem=0.0
     end if
  end if

!  ! vertical mixing 
!  allocate(Av(ToDim_nod3d))
!  Av=0.0
!  if(trim(mix_scheme)=='KPP') then
!     allocate(Kv(ToDim_nod3d,2))
!     Kv=0.0
!     call oce_mixing_kpp_init
!  else
!     allocate(Kv(ToDim_nod3d,1))
!     Kv=0.0
!  end if
!  if(trim(mix_scheme)=='MY2p5') then
!     call oce_mixing_MY2p5_init
!  end if
!  if(tidal_mixing) call oce_mixing_tidal_init 
! 
!  ! initialize the fct scheme
!#ifdef use_tracer_fct    
!  call fct_init   
!#endif
!
!  if(mype==0) write(*,*) 'Ocean arrays have been set up'
  print*, 'Arrays are set: '

end subroutine ocean_array_setup

subroutine oce_input
  ! read restart fields for ocean dynamics and active tracer variables
  use o_param
  use o_mesh
  use o_array
  use g_clock
  use g_config
  use g_PARFE
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, j, dimid_rec, nrec
  integer                   :: ssh_varid, tra_varid(2)
  integer                   :: u_varid, v_varid, w_varid, wpot_varid
  integer                   :: istart(2), icount(2), n3
  character(100)            :: filename
  character(1)              :: trind
  real(kind=8), allocatable :: aux2(:), aux3(:) 

  allocate(aux2(myDim_nod2D), aux3(myDim_nod3D)) 
  n3=myDim_nod3D           
  nrec=day2ext
  print*, 'nrec: ', nrec

  ! open files
  filename=trim(ResultPath)//runid//'.'//runyear//'.oce.mean.nc'
  print*,'Input filename: ',filename
  status = nf_open(filename, nf_nowrite, ncid)
  if (status .ne. nf_noerr) call handle_err(status)
  ! inquire variable id
  status=nf_inq_varid(ncid, 'ssh', ssh_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'u', u_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'v', v_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'w', w_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'temp', tra_varid(1))
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'salt', tra_varid(2))
  if (status .ne. nf_noerr) call handle_err(status)

  ! read variables

  ! 2d fields
  istart=(/1,nrec/)
  icount=(/myDim_nod2D, 1/)
  status=nf_get_vara_double(ncid, ssh_varid, istart, icount, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
  ssh=aux2(1:myDim_nod2D)         

  ! 3d fields
  istart=(/1,nrec/)
  icount=(/myDim_nod3D, 1/)

  do j=1,2
     status=nf_get_vara_double(ncid, tra_varid(j), istart, icount, aux3) 
     if (status .ne. nf_noerr) call handle_err(status)
     tracer(:,j)=aux3(1:myDim_nod3D)  
  end do

  status=nf_get_vara_double(ncid, u_varid, istart, icount, aux3)
  if (status .ne. nf_noerr) call handle_err(status)
  uf(1:n3)=aux3(1:myDim_nod3D)     

  status=nf_get_vara_double(ncid, v_varid, istart, icount, aux3) 
  if (status .ne. nf_noerr) call handle_err(status)
  uf(1+n3:2*n3)=aux3(1:myDim_nod3D)  

  status=nf_get_vara_double(ncid, w_varid, istart, icount, aux3) 
  uf(1+2*n3:3*n3)=aux3(myDim_nod3D)
  if (status .ne. nf_noerr) call handle_err(status)

  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! the next record to be saved
  save_count=nrec+1

  deallocate(aux3, aux2)   

end subroutine oce_input

subroutine handle_err(errcode) 
use g_parfe 
implicit none 
       
#include "netcdf.inc"  

integer errcode 
         
write(*,*) 'Error: ', nf_strerror(errcode) 
!  call par_ex 
stop
end subroutine handle_err
