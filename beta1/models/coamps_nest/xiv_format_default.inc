! This code is not protected by the DART copyright agreement.
! DART $Id$

!     innovation vector flat file format
! 
! rcs keywords: &RCSfile& 
!               &Revision: 2479 & &Date: 2006-08-21 00:41:09 +0000 (Mon, 21 Aug 2006) & 
! SCCS IDENTIFICATION:  %W% %G% 
!                       %U% %P% 
!
!     innovation file headers:
!
!        header 1: grid parameters in lines 1-13 
!            'igrid=', igrid
!            'iref=', iref
!            'jref=', jref
!            'im=', im
!            'jm=', jm
!            'lm=', lm
!            'reflat=', reflat
!            'reflon=', reflon
!            'stdlt1=', stdlt1
!            'stdlt2=', stdlt2
!            'stdlon=', stdlon
!            'delx=', delx
!            'dely=', dely
!
!        header 2: pressure levels in lines 1-13 
!            'pranal', pranal
!
!        header 3: box number and location in lines 1-13 
!        number of boxes and the ending index of each box, ne_ob
!        nb_ob is easily computed from ne_ob, so it is not written
!
!            'n_boxm=', n_boxm,'ne_ob', ne_ob
!
!        header 4: number of observations,
!                background dtg, hour   in lines 1-13 
!            'number of obs = ',nob
!            ' cdtg_bk = ',cdtg_bk,' tau_bk = ',int(tau_bk)
!
!        header 5: variable list
!
!
!     innovation file variables:
!           n                      !  ob number
!           ob(n)                  !  ob value
!           bk_ob(n)               !  background value
!           t_bk_ob(n)             !  temperature of background
!           xiv_ob(n)              !  innovation (ob - background)
!           err_ob(n)              !  observation error
!           etc_ob(n)              !  specialized parameter, e.g. v for ssmi ff
!           rlat_ob(n)             !  latitude
!           rlon_ob(n)             !  longitude
!           p_ob(n)                !  pressure
!           jvarty_ob(n)           !  variable type
!           insty_ob(n)            !  instrument type
!           nvprof_ob(n)           !  number of values in profile
!           ichk_ob(n)             !  quality control values
!           idt_ob(n)              !  time in seconds from analysis time
!           c_pf_ob(n)             !  platform identifier
!           c_db_ob(n)             !  data base address of original val
!           idp_ob(n)              !  sfc pressure change at ob location
!
! ------------------------------------------------------------------------------
!  sample header
!12345678901234567890123456789012345678901234567890
! ------------------------------------------------------------------------------
!
!    igrid=              2
!     iref=              1
!     jref=             46
!       im=          16761
!       jm=              1
!       lm=             21
!   reflat=        33.0000
!   reflon=       230.0000
!   stdlt1=        60.0000
!   stdlt2=        30.0000
!   stdlon=       270.0000
!     delx=     12000.0000
!     dely=     12000.0000
!
!    pranal
!1013.20
!1000.00
! 975.00
! 950.00
! 925.00
! 900.00
! 850.00
! 800.00
! 700.00
! 500.00
! 400.00
! 300.00
! 250.00
! 200.00
! 150.00
! 100.00
!  70.00
!  50.00
!  30.00
!  20.00
!  10.00
!
!   n_boxm=      1
!     ne_ob
!     0
!    number of obs =   117973  cdtg_bk = 2003051206 tau_bk =      6
!      n       ob       bk     t_bk       iv      err       etc       lat       lon        p  vty  ity   nvp  chk      dt                pf         org  idp
!      1    52.00   170.74   295.91  -118.74     9.00     52.00      4.35    305.58  1000.00    1    1     3  -66  -10800  81209  msfc_lnd   d_surface     0
! ------------------------------------------------------------------------------
! end sample header
!1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
! ------------------------------------------------------------------------------
      character (len=29), parameter :: fmt_hdr1 = & 
       '(6(/,a10,i15),7(/,a10,f15.4))'
!
      character (len=16), parameter :: fmt_hdr2 =  &
       '(/,a10,/,(f7.2))'
!
      character (len=25), parameter :: fmt_hdr3 = & 
       '(/,a10,i7,/,a10,/,(12i6))'
!
      character (len=23), parameter :: fmt_hdr4 =  &
       '(a20,i8,a12,a10,a10,i6)'
!
      character (len=88), parameter :: fmt_hdr5 = &
       '(a7,5(1x,a8),3(1x,a9),1x,a8,2(1x,a4)'// &
       ',1x,a5,1x,a4,1x,a7,1x,a17,1x,a11,1x,a4,1x,a8)'
!
      character (len=91), parameter :: fmt_xiv =  & 
       '(i7,5(1x,f8.2),3(1x,f9.2),1x,f8.2,2(1x,i4)'// &
       ',1x,i5,1x,i4,1x,i7,1x,a17,1x,a11,1x,i4,1x,f8.2)'

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
