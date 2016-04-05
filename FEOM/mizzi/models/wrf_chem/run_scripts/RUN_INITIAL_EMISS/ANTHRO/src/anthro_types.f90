
   module anthro_types

   implicit none

   integer, parameter :: linemax = 500
   integer, parameter :: linsize = 132
   integer, parameter :: namsize = 32
   integer, parameter :: maxsrc = 50

   type glb_att
     integer :: len
     integer :: type
     character(len=132)  :: name
     integer(1), pointer :: attr_byte(:)
     integer(2), pointer :: attr_short(:)
     integer, pointer    :: attr_int(:)
     real, pointer       :: attr_real(:)
     real(8), pointer    :: attr_dbl(:)
     character(len=256)  :: attr_char
   end type

   type anthro_map_type
     integer           :: src_cnt                         ! count of src species
     character(len=namsize) :: emis_name                  ! emission species name
     character(len=namsize) :: src_var(maxsrc)            ! src species names
     real              :: src_wght(maxsrc)                ! multiplier for each src species
     real, allocatable :: cat_wght(:,:)                   ! multiplier for sub cats for each src species
     real, allocatable :: emission(:,:)                   ! emission (kg/m^2/s)
     logical           :: is_gas                          ! .t. => gas phase, .f. => aerosol
   end type anthro_map_type

   type dates
     integer :: date
     integer :: secs
   end type dates

   type data_file_type
     integer :: ntimes
     integer :: ncid_lo, ncid_hi
     integer :: lo_tndx, hi_tndx
     integer :: lo_buf_ndx, hi_buf_ndx
     integer :: gap_date, gap_secs
     integer :: grid_ndx
     integer, allocatable :: date(:)
     integer, allocatable :: secs(:)
     real    :: missing_value = 1.e36
     real    :: dels
     real    :: molecw
     real    :: con_fac(2)
     real, allocatable :: emis(:,:,:,:)
     real, allocatable :: src_data(:,:)
     character(len=256) :: filespec
     character(len=128) :: filename
     logical :: read_lo_tndx
     logical :: read_hi_tndx
     logical :: in_gap
     logical :: t_interp
     logical :: active
     logical, allocatable :: cat_active(:)
   end type data_file_type

   end module anthro_types
