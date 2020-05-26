! This code may (or may not) be part of the COAMPS distribution,
! So it is not protected by the DART copyright agreement.
!
! DART $Id: $

! get_name_info
! -------------
! Given the value of module-specific constant values and a name of 
! a variable in the COAMPS restart file, returns an integer 
! representing the dimension type and a two-dimensional integer 
! array containing the position of that variable in that particular 
! dimension set in the COAMPS restart file written by either 
! single-processor or multi-processor I/O. 
!  PARAMETERS
!   IN  DIM_TYPE_2D         Constant for 2-D dimension
!   IN  DIM_TYPE_3D         Constant for 3-D dimension
!   IN  DIM_TYPE_3DW        Constant for 3-D (w level) dimension
!   IN  SINGLEIO            Constant for single-processor I/O
!   IN  MULTIIO             Constant for multi-processor I/O
!   IN  var_name            the name of the variable to look up
!   OUT var_dim_type        the dimension type (2d/3d/3dw)
!                           1: 2D  2: 3D  3: 3DW
!   OUT var_record_num      the position of the variable in its
!                           particular dimension
! -------------

subroutine get_name_info(DIM_TYPE_2D, DIM_TYPE_3D, DIM_TYPE_3DW,   &
                         SINGLEIO, MULTIIO, var_name, var_dim_type,&
                         var_record_num)
  integer, intent(in)                :: DIM_TYPE_2D
  integer, intent(in)                :: DIM_TYPE_3D
  integer, intent(in)                :: DIM_TYPE_3DW
  integer, intent(in)                :: SINGLEIO
  integer, intent(in)                :: MULTIIO
  character(len=*), intent(in)       :: var_name
  integer, intent(out)               :: var_dim_type
  integer, dimension(2), intent(out) :: var_record_num

  select case (trim(var_name))
  case('aalhs')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 1
    var_record_num(MULTIIO)  = 1
  case('acx')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 104
  case('acy')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 105
  case('akhm')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 2
    var_record_num(MULTIIO)  = 2
  case('akhs')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 118
  case('akhu')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 3
    var_record_num(MULTIIO)  = 3
  case('akhv')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 4
    var_record_num(MULTIIO)  = 4
  case('akms')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 119
  case('alb')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 114
  case('alb_usgs')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 115
  case('albclim')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 116
  case('albed')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 5
    var_record_num(MULTIIO)  = 5
  case('albveg')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 117
  case('aln')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 6
    var_record_num(MULTIIO)  = 6
  case('am_urban')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 63
  case('amoob')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 7
    var_record_num(MULTIIO)  = 7
  case('anca')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 1
    var_record_num(MULTIIO)  = 1
  case('aoz')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 8
    var_record_num(MULTIIO)  = 8
  case('asl')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 61
  case('asmax')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 120
  case('atl')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 62
  case('auz')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 9
    var_record_num(MULTIIO)  = 9
  case('avlhf')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 11
    var_record_num(MULTIIO)  = 11
  case('avshf')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 12
    var_record_num(MULTIIO)  = 12
  case('avstx')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 13
    var_record_num(MULTIIO)  = 13
  case('avsty')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 14
    var_record_num(MULTIIO)  = 14
  case('avz')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 10
    var_record_num(MULTIIO)  = 10
  case('bblhs')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 2
    var_record_num(MULTIIO)  = 2
  case('bcr')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 106
  case('bcx')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 107
  case('bcy')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 108
  case('beta_s')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 121
  case('bexp')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 122
  case('bflux')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 15
    var_record_num(MULTIIO)  = 15
  case('blht')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 16
    var_record_num(MULTIIO)  = 16
  case('br')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 123
  case('cclhs')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 3
    var_record_num(MULTIIO)  = 3
  case('ccx')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 109
  case('ccy')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 110
  case('cd')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 124
  case('charnock')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 125
  case('clds')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 64
  case('cmc')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 126
  case('coht')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 65
  case('cond')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 1
    var_record_num(MULTIIO)  = 1
  case('cosz_urb2d')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 127
  case('cqq')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 128
  case('cvclds')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 66
  case('cwave')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 129
  case('deltax')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 130
  case('deltay')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 131
  case('dew')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 132
  case('diab')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 2
    var_record_num(MULTIIO)  = 2
  case('dimpl')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 4
    var_record_num(MULTIIO)  = 4
  case('dish')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 67
  case('distx')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 17
    var_record_num(MULTIIO)  = 17
  case('disty')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 18
    var_record_num(MULTIIO)  = 18
  case('div3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 3
    var_record_num(MULTIIO)  = 3
  case('dksat')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 133
  case('dqcdt')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 29
    var_record_num(MULTIIO)  = 29
  case('dqidt')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 69
  case('dqrdt')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 30
    var_record_num(MULTIIO)  = 30
  case('dqsdt')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 68
  case('dqvdt')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 31
    var_record_num(MULTIIO)  = 31
  case('drag_urban')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 134
  case('drip')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 135
  case('dswtop')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 136
  case('dtdt')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 32
    var_record_num(MULTIIO)  = 32
  case('dtdts')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 19
    var_record_num(MULTIIO)  = 19
  case('dwave')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 137
  case('dwsat')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 138
  case('dxav1')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 20
    var_record_num(MULTIIO)  = 20
  case('dxav2')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 21
    var_record_num(MULTIIO)  = 21
  case('dxm')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 22
    var_record_num(MULTIIO)  = 22
  case('dxu')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 23
    var_record_num(MULTIIO)  = 23
  case('dxv')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 24
    var_record_num(MULTIIO)  = 24
  case('dyav1')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 25
    var_record_num(MULTIIO)  = 25
  case('dyav2')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 26
    var_record_num(MULTIIO)  = 26
  case('dym')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 27
    var_record_num(MULTIIO)  = 27
  case('dyu')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 28
    var_record_num(MULTIIO)  = 28
  case('dyv')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 29
    var_record_num(MULTIIO)  = 29
  case('dzsdx')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 30
    var_record_num(MULTIIO)  = 30
  case('dzsdxm')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 31
  case('dzsdy')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 31
    var_record_num(MULTIIO)  = 32
  case('dzsdym')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 33
  case('e1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 8
    var_record_num(MULTIIO)  = 8
  case('e2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 9
    var_record_num(MULTIIO)  = 9
  case('e3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 10
    var_record_num(MULTIIO)  = 10
  case('ebi')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 48
    var_record_num(MULTIIO)  = 89
  case('ebs')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 49
    var_record_num(MULTIIO)  = 90
  case('ec')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 139
  case('edir')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 140
  case('edis')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 32
    var_record_num(MULTIIO)  = 34
  case('eimpl')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 5
    var_record_num(MULTIIO)  = 5
  case('ekh')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 4
    var_record_num(MULTIIO)  = 4
  case('ekm')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 5
    var_record_num(MULTIIO)  = 5
  case('eks')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 6
    var_record_num(MULTIIO)  = 6
  case('embck')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 145
  case('emiss')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 144
  case('esnow')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 143
  case('etp')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 141
  case('ett')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 142
  case('exbm')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 7
    var_record_num(MULTIIO)  = 7
  case('exbw')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 6
    var_record_num(MULTIIO)  = 6
  case('f1')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 152
  case('f_roof')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 147
  case('f_urban')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 148
  case('fcd')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 111
  case('fcx')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 112
  case('fcy')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 113
  case('fimpl')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 7
    var_record_num(MULTIIO)  = 7
  case('flx1')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 149
  case('flx2')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 150
  case('flx3')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 151
  case('fm')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 33
    var_record_num(MULTIIO)  = 35
  case('frzx')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 146
  case('fu')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 34
    var_record_num(MULTIIO)  = 36
  case('fv')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 35
    var_record_num(MULTIIO)  = 37
  case('g_urb2d')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 155
  case('gflux')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 153
  case('gimpl')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 8
    var_record_num(MULTIIO)  = 8
  case('grdi')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 36
    var_record_num(MULTIIO)  = 38
  case('grdj')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 37
    var_record_num(MULTIIO)  = 39
  case('grdrot')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 38
    var_record_num(MULTIIO)  = 40
  case('grntyp')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 154
  case('gwet')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 39
    var_record_num(MULTIIO)  = 41
  case('gwet_usgs')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 158
  case('gwetcl')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 42
  case('h_urban')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 159
  case('hflxl')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 40
    var_record_num(MULTIIO)  = 43
  case('hflxs')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 41
    var_record_num(MULTIIO)  = 44
  case('himpl')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 9
    var_record_num(MULTIIO)  = 9
  case('hlong')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 42
    var_record_num(MULTIIO)  = 45
  case('hs')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 156
  case('hwave')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 157
  case('hxm')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 43
    var_record_num(MULTIIO)  = 46
  case('hxu')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 46
    var_record_num(MULTIIO)  = 49
  case('hxv')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 44
    var_record_num(MULTIIO)  = 47
  case('hym')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 45
    var_record_num(MULTIIO)  = 48
  case('hyu')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 47
    var_record_num(MULTIIO)  = 50
  case('hyv')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 48
    var_record_num(MULTIIO)  = 51
  case('imaskm')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 229
  case('imasku')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 230
  case('imaskv')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 231
  case('k_urban')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 161
  case('kdt')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 160
  case('ktop')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 228
  case('lh_urb2d')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 162
  case('nc1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 48
  case('nc2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 49
  case('nc3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 50
  case('ncbi')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 75
  case('ncbs')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 102
  case('ncn1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 54
  case('ncn2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 55
  case('ncn3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 56
  case('ni1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 70
  case('ni2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 71
  case('ni3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 72
  case('nibi')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 73
  case('nibs')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 103
  case('nnbi')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 74
  case('nnbs')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 105
  case('nr1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 51
  case('nr2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 52
  case('nr3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 53
  case('nrbi')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 76
  case('nrbs')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 104
  case('nroot')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 227
  case('omg_urb2d')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 163
  case('p1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 11
    var_record_num(MULTIIO)  = 11
  case('p1init')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 78
  case('p2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 12
    var_record_num(MULTIIO)  = 12
  case('p3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 13
    var_record_num(MULTIIO)  = 13
  case('pc')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 164
  case('phi')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 49
    var_record_num(MULTIIO)  = 52
  case('phim')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 50
    var_record_num(MULTIIO)  = 53
  case('phis')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 51
    var_record_num(MULTIIO)  = 54
  case('ppbi')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 50
    var_record_num(MULTIIO)  = 91
  case('ppbs')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 51
    var_record_num(MULTIIO)  = 92
  case('precp')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 52
    var_record_num(MULTIIO)  = 55
  case('prlcl')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 53
    var_record_num(MULTIIO)  = 56
  case('prtop')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 54
    var_record_num(MULTIIO)  = 57
  case('psfc')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 55
    var_record_num(MULTIIO)  = 58
  case('psih')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 56
    var_record_num(MULTIIO)  = 59
  case('psim')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 57
    var_record_num(MULTIIO)  = 60
  case('psisat')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 165
  case('psr')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 58
    var_record_num(MULTIIO)  = 61
  case('psrbm')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 166
  case('psrpb')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 167
  case('q2_urb2d')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 170
  case('q_urban')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 169
  case('qbm')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 79
  case('qc1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 33
    var_record_num(MULTIIO)  = 33
  case('qc2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 34
    var_record_num(MULTIIO)  = 34
  case('qc3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 35
    var_record_num(MULTIIO)  = 35
  case('qc_urb2d')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 168
  case('qcbi')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 52
    var_record_num(MULTIIO)  = 93
  case('qcbs')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 53
    var_record_num(MULTIIO)  = 94
  case('qg1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 45
    var_record_num(MULTIIO)  = 45
  case('qg2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 46
    var_record_num(MULTIIO)  = 46
  case('qg3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 47
    var_record_num(MULTIIO)  = 47
  case('qgbi')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 60
    var_record_num(MULTIIO)  = 77
  case('qgbs')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 61
    var_record_num(MULTIIO)  = 101
  case('qi1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 36
    var_record_num(MULTIIO)  = 36
  case('qi2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 37
    var_record_num(MULTIIO)  = 37
  case('qi3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 38
    var_record_num(MULTIIO)  = 38
  case('qibi')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 54
    var_record_num(MULTIIO)  = 95
  case('qibs')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 55
    var_record_num(MULTIIO)  = 96
  case('qr1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 39
    var_record_num(MULTIIO)  = 39
  case('qr2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 40
    var_record_num(MULTIIO)  = 40
  case('qr3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 41
    var_record_num(MULTIIO)  = 41
  case('qrbi')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 56
    var_record_num(MULTIIO)  = 97
  case('qrbs')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 57
    var_record_num(MULTIIO)  = 98
  case('qs1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 42
    var_record_num(MULTIIO)  = 42
  case('qs2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 43
    var_record_num(MULTIIO)  = 43
  case('qs3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 44
    var_record_num(MULTIIO)  = 44
  case('qsbi')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 58
    var_record_num(MULTIIO)  = 99
  case('qsbs')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 59
    var_record_num(MULTIIO)  = 100
  case('qstar')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 59
    var_record_num(MULTIIO)  = 62
  case('quartz')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 171
  case('qv1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 14
    var_record_num(MULTIIO)  = 14
  case('qv2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 15
    var_record_num(MULTIIO)  = 15
  case('qv3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 16
    var_record_num(MULTIIO)  = 16
  case('qvbi')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 62
    var_record_num(MULTIIO)  = 106
  case('qvbs')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 63
    var_record_num(MULTIIO)  = 107
  case('qvsea')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 60
    var_record_num(MULTIIO)  = 63
  case('rainc')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 61
    var_record_num(MULTIIO)  = 64
  case('rainc_ts')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 173
  case('raink')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 62
    var_record_num(MULTIIO)  = 65
  case('rainks')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 172
  case('rainp')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 64
    var_record_num(MULTIIO)  = 67
  case('rains')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 63
    var_record_num(MULTIIO)  = 66
  case('rains_ts')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 174
  case('rainsz')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 80
  case('raint')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 65
    var_record_num(MULTIIO)  = 68
  case('raintend')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 175
  case('rbm')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 17
    var_record_num(MULTIIO)  = 17
  case('rbw')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 10
    var_record_num(MULTIIO)  = 10
  case('rc')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 176
  case('rcq')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 177
  case('rcs')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 178
  case('rcsoil')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 179
  case('rct')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 180
  case('rgl')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 181
  case('rib')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 66
    var_record_num(MULTIIO)  = 69
  case('rn_urb2d')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 182
  case('rsmin')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 183
  case('runoff1')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 184
  case('runoff2')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 185
  case('runoff3')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 186
  case('sh_urb2d')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 191
  case('shdfac')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 187
  case('shdmax')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 188
  case('shdmin')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 189
  case('sigdt')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 11
    var_record_num(MULTIIO)  = 11
  case('slpeta')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 192
  case('slpr')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 67
    var_record_num(MULTIIO)  = 70
  case('slw')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 68
    var_record_num(MULTIIO)  = 71
  case('smcdry')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 193
  case('smcmax')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 194
  case('smcref')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 195
  case('smcwlt')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 196
  case('sncovr')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 197
  case('sneqv')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 198
  case('snoalb')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 190
  case('snomlt')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 199
  case('snotime')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 200
  case('snow')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 69
    var_record_num(MULTIIO)  = 72
  case('snow_ts')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 201
  case('snowd')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 70
    var_record_num(MULTIIO)  = 73
  case('snowx')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 71
    var_record_num(MULTIIO)  = 74
  case('snup')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 202
  case('soilm')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 203
  case('soiltyp')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 205
  case('soilw')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 204
  case('spdmax')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 101
    var_record_num(MULTIIO)  = 226
  case('ssw')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 72
    var_record_num(MULTIIO)  = 75
  case('stres')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 73
    var_record_num(MULTIIO)  = 76
  case('strx')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 74
    var_record_num(MULTIIO)  = 77
  case('stry')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 75
    var_record_num(MULTIIO)  = 78
  case('t1000')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 80
    var_record_num(MULTIIO)  = 83
  case('tb_urb2d')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 206
  case('tc_urb2d')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 207
  case('tg_urb2d')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 208
  case('th1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 19
    var_record_num(MULTIIO)  = 19
  case('th1init')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 84
  case('th1p')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 81
  case('th2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 20
    var_record_num(MULTIIO)  = 20
  case('th2_urb2d')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 211
  case('th2p')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 82
  case('th3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 21
    var_record_num(MULTIIO)  = 21
  case('th3p')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 83
  case('thbi')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 64
    var_record_num(MULTIIO)  = 108
  case('thbm')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 18
    var_record_num(MULTIIO)  = 18
  case('thbs')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 65
    var_record_num(MULTIIO)  = 109
  case('thbw')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 12
    var_record_num(MULTIIO)  = 12
  case('thdiff')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 85
  case('tr_urb2d')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 209
  case('trad')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 22
    var_record_num(MULTIIO)  = 22
  case('troof')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 86
  case('ts_urb2d')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 210
  case('tsea')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 76
    var_record_num(MULTIIO)  = 79
  case('tsoil')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 77
    var_record_num(MULTIIO)  = 80
  case('tstar')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 78
    var_record_num(MULTIIO)  = 81
  case('ttau0')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 79
    var_record_num(MULTIIO)  = 82
  case('u1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 23
    var_record_num(MULTIIO)  = 23
  case('u10m')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 97
    var_record_num(MULTIIO)  = 100
  case('u1init')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 87
  case('u2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 24
    var_record_num(MULTIIO)  = 24
  case('u3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 25
    var_record_num(MULTIIO)  = 25
  case('ubi')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 66
    var_record_num(MULTIIO)  = 110
  case('ubs')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 67
    var_record_num(MULTIIO)  = 111
  case('uc_urb2d')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 212
  case('ulwtop')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 213
  case('ustar')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 81
    var_record_num(MULTIIO)  = 84
  case('uzgeu')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 57
  case('uzgev')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 58
  case('v1')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 26
    var_record_num(MULTIIO)  = 26
  case('v10m')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 98
    var_record_num(MULTIIO)  = 101
  case('v1init')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 88
  case('v2')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 27
    var_record_num(MULTIIO)  = 27
  case('v3')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 28
    var_record_num(MULTIIO)  = 28
  case('vbi')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 68
    var_record_num(MULTIIO)  = 112
  case('vbs')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = 69
    var_record_num(MULTIIO)  = 113
  case('vegtyp')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 214
  case('vzgeu')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 59
  case('vzgev')
    var_dim_type = DIM_TYPE_3D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 60
  case('w0avg')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 16
    var_record_num(MULTIIO)  = 16
  case('w1')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 13
    var_record_num(MULTIIO)  = 13
  case('w2')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 14
    var_record_num(MULTIIO)  = 14
  case('w3')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 15
    var_record_num(MULTIIO)  = 15
  case('wbi')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 17
    var_record_num(MULTIIO)  = -1
  case('wbs')
    var_dim_type = DIM_TYPE_3DW
    var_record_num(SINGLEIO) = 18
    var_record_num(MULTIIO)  = -1
  case('wtm')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 82
    var_record_num(MULTIIO)  = 85
  case('wtu')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 83
    var_record_num(MULTIIO)  = 86
  case('wtv')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 84
    var_record_num(MULTIIO)  = 87
  case('wvst')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 217
  case('wvsu')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 215
  case('wvsv')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 216
  case('xlai')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 218
  case('xland')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 85
    var_record_num(MULTIIO)  = 88
  case('xlulc')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 219
  case('xnr')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 99
    var_record_num(MULTIIO)  = 102
  case('xpos')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 86
    var_record_num(MULTIIO)  = 89
  case('xxxb_urb2d')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 221
  case('xxxc_urb2d')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 223
  case('xxxg_urb2d')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 222
  case('xxxr_urb2d')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 220
  case('ynr')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 100
    var_record_num(MULTIIO)  = 103
  case('ypos')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 87
    var_record_num(MULTIIO)  = 90
  case('z0')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 96
    var_record_num(MULTIIO)  = 99
  case('z0_usgs')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 225
  case('z0c')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = -1
    var_record_num(MULTIIO)  = 224
  case('zaol')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 88
    var_record_num(MULTIIO)  = 91
  case('zerom')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 89
    var_record_num(MULTIIO)  = 92
  case('zerou')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 90
    var_record_num(MULTIIO)  = 93
  case('zerov')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 91
    var_record_num(MULTIIO)  = 94
  case('zsfc')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 92
    var_record_num(MULTIIO)  = 95
  case('zsfcu')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 93
    var_record_num(MULTIIO)  = 96
  case('zsfcv')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 94
    var_record_num(MULTIIO)  = 97
  case('zvkrm')
    var_dim_type = DIM_TYPE_2D
    var_record_num(SINGLEIO) = 95
    var_record_num(MULTIIO)  = 98
  case default
    write (*,*) "Can't match name " // var_name // "in restart"
  end select
end subroutine get_name_info

! <next few lines under version control, do not edit>
! $URL: $
! $Id: $
! $Revision: $
! $Date: $
