!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef __sgi
#ifdef _COMPILER_VERSION
!the MIPSPro compiler defines _COMPILER_VERSION
#define sgi_mipspro
#else
#define sgi_generic
#endif
#endif

#if defined(_CRAY) || defined(sgi_mipspro)
#define SGICRAY
#endif

!parallel machine types
#if defined(_CRAY) && !defined(_CRAYT3E) && !defined(_CRAYT3D)
#define CRAYPVP
#endif

#if defined(_CRAYT3E) || defined(_CRAYT3D) || defined(sgi_mipspro)
#define SGICRAY_MPP
#endif

!most compilers support Cray pointers
!if you find a compiler that doesn't, #undef this inside a suitable #ifdef
#define use_CRI_pointers

#if defined __SX || defined __SXdbl4 || __IFC
! NEC-SX and Intel fortran compiler have limited Cray pointer support
!for NEC, user must include the cpp flag -D__SX by hand, not automatic.
!The intel compiler automatically sets __IFC.
#undef use_CRI_pointers
#endif

!values of kind: double and long are 8-byte, float and int are 4-byte
!pointer_kind is used for storing addresses as integers
#if defined(SGICRAY)
#define DOUBLE_KIND 8
#define FLOAT_KIND 4
#define LONG_KIND 8
#define INT_KIND 4
#define SHORT_KIND 2
#define POINTER_KIND 8
#else
!these might be different on non-SGICRAY, I believe
#define DOUBLE_KIND 8
#define FLOAT_KIND 4
#define LONG_KIND 8
#define INT_KIND 4
#define SHORT_KIND 2
#define POINTER_KIND 8
#endif

#ifdef sgi_generic
!this is for the Edinburgh n32/o32 compiler, which won't accept 8-byte ints at any price
#define no_8byte_integers
#define LONG_KIND 4
#endif

#ifdef __SXdbl4
!When -A dbl4 is used on NEC-SX both 4-byte reals become 8-byte reals.
!(and 8-byte reals stay 8-byte reals, so they are both the same)
!by forbidding 4-byte reals, 4-byte cmplx is also forbidden
#define no_4byte_reals
!I think by redefining FLOAT_KIND to 8, I no longer need to redefine NF_*
!but I'll leave these in for now.
#define FLOAT_KIND 8
#define NF_GET_VAR_REAL nf_get_var_double
#define NF_GET_VARA_REAL nf_get_vara_double
#define NF_GET_ATT_REAL nf_get_att_double
#endif
