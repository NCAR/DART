
   module attr_types

   implicit none

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

   type dom_glb_att
     integer :: ngatts
     type(glb_att), pointer :: attrs(:)
   end type

   end module attr_types
