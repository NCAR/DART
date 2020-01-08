
! FESOM 2 (Finite-volumE Sea ice-Ocean Model)
! multi-resolution ocean general circulation model
! FESOM/fesom2 is licensed under the GNU General Public License v2.0
! Copyright (C) 2018  FESOM team
!
! This source file was taken from  the FESOM V1.4 modules

subroutine ocean_mesh_setup
  use o_param
  use o_elements
  use o_mesh
  use g_config
  use g_parfe
  implicit none
  


  call mesh_scaling                         ! long., lat. are transf. into rad

  call standard_element_definition_2D       ! Basis functions and scalar
  call standard_element_definition_3D        
  call basisfunctions_3D                    ! Local basis f. deriv. are computed
  call basisfunctions_2D                     

  print*, 'element volume and basis function derivatives prepared' 

  call check_mesh_quality_resolution

  call build_nghbr_arrays                   ! Builds arrays nod_in_elem2D

  call find_bottom_nodes

!=$0  if(grid_type/=1 .or. Redi_GM) call find_layer_elem3d

!=$0
!#if defined(use_opbnd_tide) || defined(use_opbnd_restoring)
!  call check_mesh_quality_opbd
!  call build_open_boundary_arrays 	    ! build OB arrays
!#endif

!#ifdef use_fullfreesurf
  call find_updating_element    
!#endif
!=$0

end subroutine ocean_mesh_setup
!
!---------------------------------------------------------------------------
!
subroutine check_mesh_quality_resolution
  use o_PARAM
  use o_MESH
  use o_ELEMENTS
  use g_parfe
  implicit none

  real(kind=8)   :: max_area, min_area

  ! display minimal and maximal triangle area
  max_area=maxval(voltriangle)
  min_area=minval(voltriangle)
     print*, 'The minimal triangle area is: ', min_area, ' m**2'
     print*, 'The maximal triangle area is: ', max_area, ' m**2'    

end subroutine check_mesh_quality_resolution
!
!--------------------------------------------------------------
!
subroutine check_mesh_quality_opbd
  ! check open boundary nodes
  use o_PARAM
  use o_MESH
  use o_ELEMENTS
  use g_parfe
  implicit none

  integer        :: elem, elnodes(3), ind(3), cnt

  cnt=0
  do elem=1,myDim_elem2D
     elnodes=elem2D_nodes(:,elem)
     ind=index_nod3d(nod3d_below_nod2d(1,elnodes))
     if(count(ind==10)==0 .and. any(ind==12)) then
        cnt=cnt+1
     end if
  end do
  if(cnt>0) then
     print*, '-----------------------------------------------------------------------'
     print*, 'Warning:'
     print*, 'There are ', cnt, 'surface triangles without interior nodes at the'
     print*, 'open boundary. Problems can occur when defining/searching open boundary'
     print*, 'nodes. You are suggested to modify your mesh before continuing.'
     print*, '-----------------------------------------------------------------------'
  end if

  ! note: more check could be updated here, e.g., check the quality of triangles (angle size).
end subroutine check_mesh_quality_opbd
!
!--------------------------------------------------------------
!
subroutine mesh_scaling  
  !
  ! Transforms degrees in rad in coord_nod2D(2,myDim_nod2D)     
  ! Constructs the arrays cos_elem2D(myDim_elem2D)
  ! Constructs num_layers_below_nod2D, 
  ! and does transform to rad in coord_nod3D

  use g_config
  use o_param
  use o_MESH
  use o_ELEMENTS
  use g_PARFE

  implicit none

  integer         :: i,j,ind2,ind3
  integer         :: n,  node, nodeup
  real(kind=8)    :: lon,lat,rlon,rlat

  ! =======================
  !  Lon and lat to radians
  ! =======================
  coord_nod2D(1,:)=coord_nod2D(1,:)*rad
  coord_nod2D(2,:)=coord_nod2D(2,:)*rad

  ! =======================
  ! Mean cos on 2D elements
  ! This sets spherical geometry!
  ! =======================
  allocate(cos_elem2D(myDim_elem2D))
  do i=1, myDim_elem2D  
     cos_elem2D(i)=sum(cos(coord_nod2D(2,elem2D_nodes(:,i))))/3.0
  end do
  if(cartesian) cos_elem2D=1.0  

  ! =======================
  ! number of layers 
  ! =======================
  allocate(num_layers_below_nod2D(myDim_nod2D))
  num_layers_below_nod2D=-1
  do n=1,myDim_nod2D
     do j=1,max_num_layers
        node=nod3D_below_nod2D(j,n)
        if (node > 0) then
           num_layers_below_nod2D(n)=num_layers_below_nod2D(n) + 1
        else
           exit
        end if
     end do
  end do
  ! ========================
  ! Lon and lat to radians for 3D nodes
  ! ========================

  do n=1,myDim_nod2D
     !   coord_nod3D correction:
     do j=1,max_num_layers
        node=nod3D_below_nod2D(j,n)
        if (node < 1) exit
        coord_nod3D(1,node)=coord_nod2D(1,n)
        coord_nod3D(2,node)=coord_nod2D(2,n)
     end do
  end do

  ! setup geolat, which contains geographic latitude
  allocate(geolat(myDim_nod3d))
     geolat=coord_nod3d(2,:)

  print*, 'converting coord. from degree to rad. & defining num. layers done'
end subroutine mesh_scaling
!
!=====================================================================
!
subroutine standard_element_definition_2D
  !
  !    - BASISFUNCTIONS ON 2D STANDARD ELEMENTS 
  !         stdbafunc(1)=1-x-y  stdbafunc(2)=x  stdbafunc(3)=y
  !
  use o_ELEMENTS
  use o_mesh
  !
  implicit none
  !
  integer :: i,j
  !
  Vol2D =  1.0_8/2.0_8               ! Vol2D = < 1,1>
  !        
  sProd_2Di = 1.0_8/6.0_8            ! <stdbafunc(i),1.>
  !
  allocate(sProd_2Dij(3,3))
  sProd_2Dij=1.0_8/24.0_8            ! <stdbafunc(i),stdbafunc(j)>
  do j=1,3
     sProd_2Dij(j,j)=1.0_8/12.0_8 
  end do

  ! Scalar products are only required as sProd_2D/Vol2D:
  sProd_2Dij=sProd_2Dij/Vol2D
  sProd_2Di=sProd_2Di/Vol2D

  !   derivative_stdbafu_x(i,j) = d(Fi(j))/dx(i) on the standard element 
  allocate(derivative_stdbafu_x_2D(2,3))
  !
  derivative_stdbafu_x_2D= 0.
  derivative_stdbafu_x_2D(:,1)= -1.
  derivative_stdbafu_x_2D(1,2)= 1.
  derivative_stdbafu_x_2D(2,3)= 1.
  !
end subroutine standard_element_definition_2D
!
!--------------------------------------------------------------
!
subroutine standard_element_definition_3D
  !
  !    - BASISFUNCTIONS ON 3D STANDARD ELEMENTS (stdbafu)
  !      stdbafunc(1)=1-x-y-z, stdbafunc(2)=x, stdbafunc(3)=y, stdbafunc(4)=z     
  use o_ELEMENTS
  use o_MESH 
  !
  implicit none 
  !
  integer :: i,j
  !
  Vol3D = 1.0_8/6.0_8                  ! Vol3D = < 1.,1.>  
  !
  sProd_3Di = 1.0_8/24.0_8             ! sProd_3Di=<stdbafunc(i),1.>
  !
  allocate(sProd_3Dij(4,4))
  sProd_3Dij=1./120.                   ! sProd_3Dij=<stbafunc(i),stdbafunc(j)>
  do j=1,4
     sProd_3Dij(j,j)=1./60.
  end do


  ! Scalar products are only required as sProd_3D/Vol3D:
  sProd_3Dij=sProd_3Dij/Vol3D
  sProd_3Di=sProd_3Di/Vol3D

  !   derivative_stdbafu_x(j,i) = d(Fi(j))/dx(i) on the standard element
  allocate(derivative_stdbafu_x_3D(3,4))
  !
  do j=1,4
     do i=1,3
        derivative_stdbafu_x_3D(i,j)= 0.0
        if(j==1)   derivative_stdbafu_x_3D(i,j)=-1.
        if(i==j-1) derivative_stdbafu_x_3D(i,j)= 1.
     end do
  end do
  !
end subroutine standard_element_definition_3D
!
!=====================================================================
! 
subroutine basisfunctions_2D
  use o_ELEMENTS
  use o_MESH
  use g_PARFE
  implicit none
  !
  real(kind=8)                           :: DET2D
  real(kind=8), dimension(2,3)           :: derivative_locbafu_x_2D
  real(kind=8), dimension(2,2)           :: jacobian2D, jacobian2D_inv
  integer                                :: elem,i


  allocate(bafux_2d(3,myDim_elem2d), bafuy_2d(3,myDim_elem2d))
  allocate(voltriangle(myDim_elem2d))

  bafux_2d = 0.0
  bafuy_2d = 0.0
  voltriangle = 0.0

  do elem=1,myDim_elem2d
     call local_element_def_2D(elem, DET2D, derivative_locbafu_x_2D)
     do i=1,3
        bafux_2d(i,elem) = derivative_locbafu_x_2D(1,i)
        bafuy_2d(i,elem) = derivative_locbafu_x_2D(2,i)
     enddo
     voltriangle(elem) = abs(DET2D) * Vol2D
  enddo

end subroutine basisfunctions_2D
!
!=====================================================================
! 
subroutine basisfunctions_3D
  !  Derivatives of basis functions and volumes of tetrahedra are computed
  !  through transform from a standard element to a real one.
  !  Local Cartesian metrics is used
  use o_ELEMENTS
  use o_MESH
  use g_PARFE
  implicit none
  !
  real(kind=8), dimension(3,3)           :: jacobian3D
  real(kind=8), dimension(3,3)           :: jacobian3D_inv
  real(kind=8)                           :: DET3D
  real(kind=8), dimension(3,4)           :: derivative_locbafu_x_3D
  integer                                :: elem,i

  allocate(bafux_3d(4,myDim_elem3d), bafuy_3d(4,myDim_elem3d))
  allocate(bafuz_3d(4,myDim_elem3d),voltetra(myDim_elem3d))
  bafux_3d = 0.0
  bafuy_3d = 0.0
  bafuz_3d = 0.0
  voltetra = 0.0

  do elem=1,myDim_elem3d
     call local_element_def_3D(elem, DET3D, derivative_locbafu_x_3D)
     do i=1,4
        bafux_3d(i,elem) = derivative_locbafu_x_3D(1,i)
        bafuy_3d(i,elem) = derivative_locbafu_x_3D(2,i)
        bafuz_3d(i,elem) = derivative_locbafu_x_3D(3,i)
     enddo
     voltetra(elem) = abs(DET3D) * Vol3D
  enddo

end subroutine basisfunctions_3D
!
! =======================================================================
!
subroutine local_element_def_2D(element, DET, derivative_locbafu_x_2D)

  use o_ELEMENTS
  use o_MESH
  use o_param
  use g_config
  implicit none
  !
  integer, intent(IN)                        :: element
  real(kind=8), dimension(2,2)               :: jacobian2D
  real(kind=8), dimension(2,2)               :: jacobian2D_inv
  real(kind=8), intent(OUT)                  :: DET
  real(kind=8), dimension(2,3), intent(OUT)  :: derivative_locbafu_x_2D
  !
  real(kind=8), dimension(2,3)               :: local_cart
  real(kind=8), dimension(3,2)               :: der_transp
  integer                                    :: i, node
  real(kind=8)                               :: meancos
  !
  meancos=cos_elem2D(element)
  do i=1,3
     node=elem2D_nodes(i,element)
     !
     !  scaled cartesian coordinates
     !
     local_cart(1,i)=coord_nod2D(1,node) 
     local_cart(2,i)=r_earth * coord_nod2D(2,node) 
  end do
  !
  !  jacobian
  !
  do i=1,2
     jacobian2D(:,i)= local_cart(:,i+1)-local_cart(:,1)
     if (jacobian2D(1,i)> domain_length/2.0) jacobian2D(1,i)=jacobian2D(1,i)-domain_length
     if (jacobian2D(1,i)<-domain_length/2.0) jacobian2D(1,i)=jacobian2D(1,i)+domain_length
  end do
  jacobian2D(1,:)=jacobian2D(1,:)*meancos *r_earth
  !
  !  inverse of jacobian
  !
  call matrix_inverse_2x2(jacobian2D, jacobian2D_inv, DET)
  !
  der_transp=matmul(transpose(derivative_stdbafu_x_2D), jacobian2D_inv)
  derivative_locbafu_x_2D=transpose(der_transp)
  !
  !
end subroutine local_element_def_2D
!
!============================================================================
!
subroutine local_element_def_3D(element, DET, derivative_locbafu_x_3D)
  !  Auxilliary routine to basisfunctions_3D
  use o_ELEMENTS
  use o_MESH
  use o_param
  use g_config
  !
  implicit none
  !
  integer, intent(IN)                                 :: element
  real(kind=8), dimension(3,3)                        :: jacobian3D
  real(kind=8), dimension(3,3)                        :: jacobian3D_inv
  real(kind=8), intent(OUT)                           :: DET
  real(kind=8), dimension(3,4), intent(OUT)           :: derivative_locbafu_x_3D
  !
  real(kind=8), dimension(3,4)                        :: local_cart
  real(kind=8), dimension(4,3)                        :: der_transp
  integer                                             :: i, node
  real(kind=8)                                        :: meancos 

  !
  !  scaled cartesian coordinates on elements
  !

  do i=1,4
     node = elem3D_nodes(i,element)
     local_cart(1,i) = coord_nod3D(1,node)
     local_cart(2,i) = r_earth * coord_nod3D(2,node)
     local_cart(3,i) = coord_nod3D(3,node)
  end do
  meancos=cos_elem2D(elem2D_corresp_to_elem3D(element))
  !
  !  transformation matrix Xl(k)=jacobian(i,k)*Xs(i)
  !
  do i=1,3
     jacobian3D(:,i) = local_cart(:,i+1)-local_cart(:,1)
     if (jacobian3D(1,i)> domain_length/2.0) jacobian3D(1,i)=jacobian3D(1,i)-domain_length
     if (jacobian3D(1,i)<-domain_length/2.0) jacobian3D(1,i)=jacobian3D(1,i)+domain_length
  end do

  jacobian3D(1,:)=jacobian3D(1,:)*meancos*r_earth
  !
  !  inverse jacobian
  !
  call matrix_inverse_3x3(jacobian3D, jacobian3D_inv, DET)
  !
  !  derivatives of local basis functions
  !
  der_transp=matmul(transpose(derivative_stdbafu_x_3D), jacobian3D_inv)
  derivative_locbafu_x_3D=transpose(der_transp)

end subroutine local_element_def_3D
!
!=======================================================================
!
subroutine  matrix_inverse_2x2 (A, AINV, DET)
  !
  !
  implicit none
  !
  real(kind=8), dimension(2,2), intent(IN)  :: A
  real(kind=8), dimension(2,2), intent(OUT) :: AINV
  real(kind=8), intent(OUT)                 :: DET
  !
  integer                                   :: i,j
  !
  DET  = A(1,1)*A(2,2) - A(1,2)*A(2,1)
  if ( DET .eq. 0.0 )  then
     do j=1,2
        print*, (A(i,j),i=1,2)
     end do
     stop 'SINGULAR 2X2 MATRIX'
  else
     AINV(1,1) =  A(2,2)/DET
     AINV(1,2) = -A(1,2)/DET
     AINV(2,1) = -A(2,1)/DET
     AINV(2,2) =  A(1,1)/DET
  endif
end subroutine matrix_inverse_2x2
!
!=======================================================================
!
subroutine matrix_inverse_3x3(A, AINV, DET)
  !
  implicit none
  !
  real(kind=8), dimension(3,3), intent(IN)  :: A
  real(kind=8), dimension(3,3), intent(OUT) :: AINV
  real(kind=8), intent(OUT)                 :: DET
  !
  integer                                   :: i,j
  !
  AINV(1,1) =  A(2,2)*A(3,3) - A(3,2)*A(2,3)
  AINV(2,1) = -A(2,1)*A(3,3) + A(3,1)*A(2,3)
  AINV(3,1) =  A(2,1)*A(3,2) - A(3,1)*A(2,2)
  AINV(1,2) = -A(1,2)*A(3,3) + A(3,2)*A(1,3)
  AINV(2,2) =  A(1,1)*A(3,3) - A(3,1)*A(1,3)
  AINV(3,2) = -A(1,1)*A(3,2) + A(3,1)*A(1,2)
  AINV(1,3) =  A(1,2)*A(2,3) - A(2,2)*A(1,3)
  AINV(2,3) = -A(1,1)*A(2,3) + A(2,1)*A(1,3)
  AINV(3,3) =  A(1,1)*A(2,2) - A(2,1)*A(1,2)
  DET = A(1,1)*AINV(1,1) + A(1,2)*AINV(2,1) + A(1,3)*AINV(3,1)
  !
  if ( DET .eq. 0.0 )  then
     do j=1,3
        print*, (A(i,j),i=1,3)
     end do
     stop 'SINGULAR 3X3 MATRIX'
  else
     AINV = AINV/DET
  endif
  !
end subroutine matrix_inverse_3x3
!
! =======================================================================
!
subroutine build_nghbr_arrays
  !
  ! Assembles additional arrays which list for each node the elements 
  ! containing the node and node neighbours

  use o_MESH
  use o_ELEMENTS
  use g_PARFE
  implicit none

  integer                            :: j,k,m,a,tr(3),tet(4), counter, el,ml(1)
  integer, allocatable, dimension(:) :: ind
  integer, dimension(100)            :: AUX=0

  !--------------- 2D mesh:
  ! Builds nod_in_elem2D
     
  allocate(ind(myDim_nod3D))

  ind=0
  do j=1,myDim_elem2D
     tr=elem2D_nodes(:,j)
     ind(tr)=ind(tr)+1
  end do
  allocate(nod_in_elem2D(myDim_nod2D))
  nod_in_elem2D%nmb=ind(1:myDim_nod2D)    
  do j=1,myDim_nod2D   
     allocate(nod_in_elem2D(j)%addresses(ind(j)))
  end do
  ind=0
  do j=1,myDim_elem2D   
     tr=elem2D_nodes(:,j)
     ind(tr)=ind(tr)+1
     do k=1,3
        if(tr(k)<=myDim_nod2D) then 
           nod_in_elem2D(tr(k))%addresses(ind(tr(k)))=j
        end if
     end do
  end do
  ! the list of elements is ordered, and no sorting is needed

  ! Builds nghbr_nod2D
  allocate(nghbr_nod2D(myDim_nod2D))
  ind=0
  do j=1, myDim_nod2D
     counter=0
     do m=1,nod_in_elem2D(j)%nmb
        el=nod_in_elem2D(j)%addresses(m)
        do k=1, 3
           a=elem2D_nodes(k,el)       
           if (ind(a)==0) then  
              ind(a)=1 
              counter=counter+1         
              aux(counter)=a
           end if
        end do
     end do
     nghbr_nod2D(j)%nmb=counter
     allocate(nghbr_nod2D(j)%addresses(counter))

     ! we need to sort array aux(1:counter)
     do m=counter,1,-1
        ml=maxloc(aux(1:counter))
        a=ml(1)
        nghbr_nod2D(j)%addresses(m)=aux(a)
        ind(aux(a))=0
        aux(a)=-999
     end do
  end do

  !---------------3D mesh---------------------
  ! Construction of nod_in_elem3D
  ind=0
  do j=1,myDim_elem3D
     tet=elem3D_nodes(:,j)
     ind(tet)=ind(tet)+1
  end do
  allocate(nod_in_elem3D(myDim_nod3D))
  nod_in_elem3D%nmb=ind(1:myDim_nod3D)
  do j=1,myDim_nod3D   
     allocate(nod_in_elem3D(j)%addresses(ind(j)))
  end do
  ind=0
  do j=1,myDim_elem3D   
     tet=elem3D_nodes(:,j)
     ind(tet)=ind(tet)+1
     do k=1,4
        if(tet(k)<=myDim_nod3D) then                         
           nod_in_elem3D(tet(k))%addresses(ind(tet(k)))=j     
        end if
     end do
  end do
  ! the list of elements is ordered, and no sorting is needed

  ! Builds nghbr_nod3D

  allocate(nghbr_nod3D(myDim_nod3D))

  ind=0
  do j=1, myDim_nod3D
     counter=0
     do m=1,nod_in_elem3D(j)%nmb
        el=nod_in_elem3D(j)%addresses(m)
        do k=1, 4
           a=elem3D_nodes(k,el)       
           if (ind(a)==0) then  
              ind(a)=1 
              counter=counter+1         
              aux(counter)=a
           end if
        end do
     end do
     nghbr_nod3D(j)%nmb=counter
     allocate(nghbr_nod3D(j)%addresses(counter))

     ! we need to sort array(aux(1:counter))
     do m=counter,1,-1
        ml=maxloc(aux(1:counter))
        a=ml(1)
        nghbr_nod3D(j)%addresses(m)=aux(a)
        ind(aux(a))=0
        aux(a)=-999
     end do
  end do

  deallocate(ind)

end subroutine build_nghbr_arrays
!
!==========================================================================
!
subroutine find_bottom_nodes
  use o_MESH
  use o_ELEMENTS
  use g_PARFE
  implicit none
  integer         :: j
  !
  allocate(bt_nds(myDim_nod2D))
  do j=1, myDim_nod2D
     bt_nds(j)=nod3D_below_nod2D(num_layers_below_nod2D(j)+1,j)
  end do
end subroutine find_bottom_nodes
!
!==========================================================================
!
subroutine find_cluster_area
  use o_param
  use o_MESH
  use o_ELEMENTS
  use g_PARFE
  implicit none
  
  integer         :: i, elem, elnodes(3)
  real(kind=8)    :: inv3, vol, fluxlat

  inv3=1.0/3.0_8

  allocate(cluster_area_2d(myDim_nod2D))
  cluster_area_2d=0.0

  do elem=1,myDim_elem2d
     elnodes=elem2d_nodes(:,elem)
     vol=voltriangle(elem)*inv3
     cluster_area_2d(elnodes)=cluster_area_2d(elnodes)+vol
  end do


  vol=0.0
  do i=1,myDim_nod2d
     vol=vol+cluster_area_2d(i)
  end do
  ocean_area=vol 
  domain_area=ocean_area
  print*, 'Model whole surface AREA: ', domain_area

  vol=0.0
  do i=1,myDim_nod2d
     if ( (coord_nod2D(2,i)/rad).lt.aegflx_lat ) then
     vol=vol+cluster_area_2d(i)
     endif
  end do
  aegflx_area=vol
  print*, 'AegSea Flux AREA: ', aegflx_area

  vol=0.0
  do i=1,myDim_nod2d
     if ( (coord_nod2D(2,i)/rad).gt.blkflx_lat ) then
     vol=vol+cluster_area_2d(i)
     endif
  end do
  blkflx_area=vol
  print*, 'BlkSea Flux AREA: ', blkflx_area

  
end subroutine find_cluster_area
!
!==========================================================================
!
subroutine find_layer_elem3d
  !find the layer number of 3d elements
  use o_MESH
  use o_ELEMENTS
  use g_parfe

  implicit none

  integer                   :: i, j, k, elem3, tetra_nodes(4)
  integer, allocatable      :: auxind(:)

  allocate(elem3d_layer(myDim_elem3D))
  allocate(auxind(myDim_nod3d))	

  do i=1,myDim_nod2d
     do k=1,max_num_layers
        j=nod3d_below_nod2d(k,i)  
        if(j<0) exit
        auxind(j)=k
     end do
  end do
  do elem3=1, myDim_elem3D
     tetra_nodes=elem3d_nodes(:,elem3)
     elem3d_layer(elem3)=minval(auxind(tetra_nodes))
  end do

  deallocate(auxind)
end subroutine find_layer_elem3d
!
!==========================================================================
!
subroutine build_open_boundary_arrays
  use o_MESH
  use o_ELEMENTS
  use o_param
  use o_array
  use g_config
  use g_PARFE
  implicit none
  integer         	:: j, k, m, cnt, col, nodes(4), nodes2(3), nodes23(4)
  integer         	:: tri(3,4), fourth_node(4), ind(4)
  integer               :: ind2(3), edg(2,3), third_node(3)
  real(kind=8)       	:: vec1(3), vec2(3), nvec(3), vec_dir(3)
  real(kind=8)          :: meancos, node_cart(3,4)
  real(kind=8)          :: vec2d(2), nvec2d(2), vec_dir2d(2), node_cart2(2,3)
  integer, allocatable, dimension(:,:) :: auxtri, auxedg
  integer, allocatable, dimension(:)   :: aux, aux2

  ! =============
  ! First, find the list of triangles on the open boundary
  ! =============

  cnt=0
  do j=1, myDim_elem3D        
     nodes=elem3D_nodes(:,j)
     nodes23=nod2D_corresp_to_nod3D(nodes)        
     ind=index_nod3D(nod3D_below_nod2D(1,nodes23))
     if (count(ind==12)==0) cycle
     ! form all possible triangles
     tri(:,1)=nodes(1:3)
     tri(:,2)=nodes(2:4)
     tri(:,3)=nodes((/1, 3, 4/))
     tri(:,4)=nodes((/1, 2, 4/))
     fourth_node(1)=nodes(4)
     fourth_node(2)=nodes(1)
     fourth_node(3)=nodes(2)
     fourth_node(4)=nodes(3)
     do k=1, 4
        !(1) if tri(:,k) belongs to open boundary it should project into
        ! a line at the surface (the triangle lies in the vertical plane).

        nodes2= nod2D_corresp_to_nod3D(tri(:,k))
	if((nodes2(1)==nodes2(2)).or.(nodes2(1)==nodes2(3)).or.(nodes2(2)==nodes2(3))) then

           !(2) there are no wet nodes and at least one node belongs
           !to the open boundary

           ind(1:3)=index_nod3D(nod3d_below_nod2D(1,nodes2))
           if((count(ind(1:3)/=10)==3).and.(count(ind(1:3)==12)>0)) then
              ! (3) the segment given by two different nodes from nodes2
              !    can be shared by 2 triangles (in the corners on bad
              !    meshes; this however implies that there are triangles
              !    formed by boundary nodes). We assume that such triangles are
              !    absent. 
              cnt=cnt+1
              exit
           end if
	end if
     end do
  end do

  allocate(auxtri(4,cnt))
  cnt=0	
  do j=1, myDim_elem3D          
     nodes=elem3D_nodes(:,j)
     nodes23=nod2D_corresp_to_nod3D(nodes)        
     ind=index_nod3D(nod3D_below_nod2D(1,nodes23))
     if (count(ind==12)==0) cycle
     ! form all possible triangles
     tri(:,1)=nodes(1:3)
     tri(:,2)=nodes(2:4)
     tri(:,3)=nodes((/1, 3, 4/))
     tri(:,4)=nodes((/1, 2, 4/))
     fourth_node(1)=nodes(4)
     fourth_node(2)=nodes(1)
     fourth_node(3)=nodes(2)
     fourth_node(4)=nodes(3)
     do k=1, 4
        !(1) if tri(:,k) belongs to open boundary it should project into
        !a line at the surface (the triangle lies in the vertical plane).

        nodes2= nod2D_corresp_to_nod3D(tri(:,k))
	if((nodes2(1)==nodes2(2)).or.(nodes2(1)==nodes2(3)).or.(nodes2(2)==nodes2(3))) then
           !(2) there are no wet nodes and at least one node belongs
           !to the open boundary

           ind(1:3)=index_nod3D(nod3d_below_nod2D(1,nodes2))
           if((count(ind(1:3)/=10)==3).and.(count(ind(1:3)==12)>0)) then
              ! (3) the segment given by two different nodes from nodes2
              !    can be shared by 2 triangles (in the corners on bad
              !    meshes; this however implies that there are triangles
              !    formed by boundary nodes). We assume that such triangles are
              !    absent. 
              cnt=cnt+1
              auxtri(1:3,cnt)=tri(:,k)
              auxtri(4,cnt)=fourth_node(k)
              exit
           end if
	end if
     end do
  end do

  nmbr_opbnd_tri=cnt
  !if(nmbr_opbnd_tri==0) then
  !   deallocate(auxtri)     
  !   return                           
  !end if
  allocate(opbnd_tri(nmbr_opbnd_tri,4))
  do j=1,nmbr_opbnd_tri
     opbnd_tri(j,:)=auxtri(:,j)
  end do
  deallocate(auxtri)  

  ! =============
  ! Second, compute outer normal to open boundary triangles and triangle areas
  ! =============
  allocate(opbnd_nv(nmbr_opbnd_tri,4))
  do j=1, nmbr_opbnd_tri
     ! (a) compute cartesian coordinates
     do m=1,4
        node_cart(:,m)=coord_nod3D(:,opbnd_tri(j,m))
     end do
     meancos=0.
     do m=1,3
        meancos=meancos+cos(node_cart(2,m))
     end do
     meancos=meancos/3.0

     ! (b) two vectors in the triangle plane and the vector to the fourth node
     vec1=node_cart(:,2)-node_cart(:,1)
     vec2=node_cart(:,3)-node_cart(:,1)
     vec_dir=node_cart(:,4)-node_cart(:,1)

     ! (c) check for cyclicity:
     if(vec1(1)>domain_length/2.0_8) vec1(1)=vec1(1)-domain_length
     if(vec1(1)<-domain_length/2.0_8) vec1(1)=vec1(1)+domain_length
     if(vec2(1)>domain_length/2.0_8) vec2(1)=vec2(1)-domain_length
     if(vec2(1)<-domain_length/2.0_8) vec2(1)=vec2(1)+domain_length
     if(vec_dir(1)>domain_length/2.0_8) vec_dir(1)=vec_dir(1)-domain_length
     if(vec_dir(1)<-domain_length/2.0_8) vec_dir(1)=vec_dir(1)+domain_length

     vec1(1:2)=vec1(1:2)*r_earth
     vec2(1:2)=vec2(1:2)*r_earth
     vec_dir(1:2)=vec_dir(1:2)*r_earth
     ! Local Cartesian coordinates:
     vec1(1)=vec1(1)*meancos
     vec2(1)=vec2(1)*meancos
     vec_dir(1)=vec_dir(1)*meancos

     ! (d) compute vector product of vec1 and vec2
     nvec(1)=vec1(2)*vec2(3)-vec1(3)*vec2(2)
     nvec(2)=vec1(3)*vec2(1)-vec1(1)*vec2(3)
     nvec(3)=vec1(1)*vec2(2)-vec1(2)*vec2(1)
     opbnd_nv(j,4)=0.5*sqrt(sum(nvec*nvec))  ! area of surface tri in sca coord.
     ! Return to non-scaled nvec:
     nvec=nvec/sqrt(sum(nvec*nvec))

     ! (e) outer normal
     if (sum(nvec*vec_dir) > 0) nvec = -nvec
     opbnd_nv(j,1:3)=nvec
  end do

  ! ==================
  ! Third, find the list of 3d and 2d nodes on the open boundary
  ! ==================
  ! list of 3d nodes on the open boundary
  allocate(aux2(myDim_nod3D))  
  aux2=0
  do j=1, nmbr_opbnd_tri
     do m=1,3
        if(aux2(opbnd_tri(j,m))==0) aux2(opbnd_tri(j,m))=1
     end do
  end do
  nmbr_opbnd_n3D=sum(aux2)
  allocate(aux(nmbr_opbnd_n3D))
  cnt=0
  do j=1, myDim_nod3D           
     if(aux2(j)/=0) then
        cnt=cnt+1
        aux(cnt)=j
     end if
  end do
  allocate(opbnd_n3D(nmbr_opbnd_n3D))
  opbnd_n3D=aux

  ! number of 2D nodes that have index_nod3D=12 on opbnd
  cnt=0
  do j=1, myDim_nod2D        
     m=nod3D_below_nod2D(1,j)           
     if (index_nod3D(m)==12) cnt=cnt+1
  end do
  nmbr_opbnd_n2D=cnt
  ! number of 2D nodes that have index_nod3D=11 on opbnd
  cnt=0
  do j=1,nmbr_opbnd_n3D
     col=aux(j)
     if (index_nod3D(col)==11) cnt=cnt+1
  end do
  nmbr_opbnd_t2D=cnt
  ! total 2D nodes on opbnd
  nmbr_opbnd_t2D=nmbr_opbnd_n2D+nmbr_opbnd_t2D

  allocate(opbnd_n2D(nmbr_opbnd_t2D))
  allocate(mapping_opbnd_n2d(myDim_nod2d)) 
  cnt=0
  do j=1, myDim_nod2D    
     m=nod3D_below_nod2D(1,j)       
     if(index_nod3D(m)==12) then
        cnt=cnt+1
        opbnd_n2D(cnt)=j
        mapping_opbnd_n2d(j)=cnt
     end if
  end do
  cnt=nmbr_opbnd_n2D
  do j=1,nmbr_opbnd_n3D
     col=aux(j)
     if (index_nod3D(col)==11) then
        col=nod2d_corresp_to_nod3d(col)
        cnt=cnt+1
        opbnd_n2D(cnt)=col
        mapping_opbnd_n2d(col)=cnt
     end if
  end do
  deallocate(aux, aux2)

  ! ================
  ! 4th, depth at the open boundary nodes
  ! ================
  allocate(opbnd_dep(nmbr_opbnd_t2d))
  do k=1, nmbr_opbnd_t2d
     m=opbnd_n2d(k)
     cnt=num_layers_below_nod2D(m)+1
     j=nod3d_below_nod2d(cnt,m)
     opbnd_dep(k)=abs(coord_nod3d(3,j))
  end do


#ifdef use_ice
  ! ================
  ! 5th, find the list of edges on the open boundary
  ! ================
  allocate(auxedg(3,myDim_elem2D))
  cnt=0
  do j=1, myDim_elem2D
     nodes2=elem2D_nodes(:,j)
     nodes(1:3)=nod3D_below_nod2D(1,nodes2)   
     ind2=index_nod3D(nodes(1:3))             
     if (count(ind2==12)==0) cycle
     ! form all possible edges
     edg(:,1)=nodes2(1:2)
     edg(:,2)=nodes2(2:3)
     edg(:,3)=nodes2((/3, 1/))

     third_node(1)=nodes2(3)
     third_node(2)=nodes2(1)
     third_node(3)=nodes2(2)

     do k=1, 3
        ind2(1:2)=index_nod3d(nod3D_below_nod2D(1,edg(:,k)))  
        if(count(ind2(1:2)/=10)==2) then
	   if(count(ind(1:2)==12)>0) then
              cnt=cnt+1
              auxedg(1:2,cnt)=edg(:,k)
              auxedg(3,cnt)=third_node(k) 
              exit
           end if
        end if
     end do
  end do
  nmbr_opbnd_edg=cnt
  allocate(opbnd_edg(nmbr_opbnd_edg,3))
  do j=1,nmbr_opbnd_edg
     opbnd_edg(j,:)=auxedg(:,j)
  end do
  deallocate(auxedg)  

  ! ==================
  ! 6th, compute outer normal to open boundary edges and edge length
  ! ==================
  allocate(opbnd_edg_nv(nmbr_opbnd_edg, 3))
  do j=1, nmbr_opbnd_edg
     ! (a) compute cartesian coordinates
     do m=1,3
        node_cart2(:,m)=coord_nod2D(:,opbnd_edg(j,m))
     end do
     meancos=0.
     do m=1,2
        meancos=meancos+cos(node_cart2(2,m))
     end do
     meancos=meancos/2.0

     ! (b) vectors of the edge and the vector to the third node
     vec2d=node_cart2(:,2)-node_cart2(:,1)
     vec_dir2d=node_cart2(:,3)-node_cart2(:,1)

     ! (c) check for cyclicity:
     if(vec2d(1)>domain_length/2.0_8) vec2d(1)=vec2d(1)-domain_length
     if(vec2d(1)<-domain_length/2.0_8) vec2d(1)=vec2d(1)+domain_length
     if(vec_dir2d(1)>domain_length/2.0_8) vec_dir2d(1)=vec_dir2d(1)-domain_length
     if(vec_dir2d(1)<-domain_length/2.0_8) vec_dir2d(1)=vec_dir2d(1)+domain_length

     vec2d(1:2)=vec2d(1:2)*r_earth
     vec_dir2d(1:2)=vec_dir2d(1:2)*r_earth
     ! Local Cartesian coordinates:
     vec2d(1)=vec2d(1)*meancos
     vec_dir2d(1)=vec_dir2d(1)*meancos

     ! (d) normal vector
     nvec2d(1)=vec2d(2)
     nvec2d(2)=vec2d(1)
     opbnd_edg_nv(j,3)=sqrt(nvec2d(1)**2+nvec2d(2)**2)
     nvec2d=nvec2d/opbnd_edg_nv(j,3)

     ! (e) outer normal
     if (sum(nvec2d*vec_dir2d) > 0.0) nvec2d = -nvec2d
     opbnd_edg_nv(j,1:2)=nvec2d
  end do
#endif

  print*, 'open boundary grid prepared'

end subroutine build_open_boundary_arrays
!
!=========================================================================
!
#ifdef use_fullfreesurf
subroutine find_updating_element
  use o_param
  use o_ELEMENTS
  use o_MESH
  use g_parfe
  implicit none
  !
  integer                                :: el, n, n_el, elnodes(4)
  !
  allocate(map_elem(myDim_elem3d))   
  map_elem=0
  n=0
  do el=1,myDim_elem3d              
     elnodes=elem3d_nodes(:,el)
     if(any(index_nod3D(elnodes)<=12)) then  
        n=n+1
        map_elem(el)=n
     end if
  end do
  !
  allocate(voltetra_new(n), bafux_3d_new(4,n))
  allocate(bafuy_3d_new(4,n), bafuz_3d_new(4,n))
  !
  do el=1,myDim_elem3D
     n_el=map_elem(el)
     if(n_el==0) cycle
     voltetra_new(n_el)=voltetra(el)
     bafux_3d_new(:,n_el)=bafux_3d(:,el)
     bafuy_3d_new(:,n_el)=bafuy_3d(:,el)
     bafuz_3d_new(:,n_el)=bafuz_3d(:,el)
  end do
  
  print*, 'index of first layer of elements recorded' 
end subroutine find_updating_element
#endif
!
!--------------------------------------------------------------------------
!
#ifdef use_fullfreesurf
subroutine update_mesh
  use o_param
  use o_array
  use o_mesh
  use o_ELEMENTS
  use g_parfe
  implicit none
  !
  integer                              :: m, el, n_el, i, row
  real(kind=8), dimension(3,3)         :: jacobian3D
  real(kind=8), dimension(3,3)         :: jacobian3D_inv
  real(kind=8), dimension(3,4)         :: derivative_locbafu_x_3D
  real(kind=8)                         :: DET3D
  
  ! update z coordiante
  do m=1,myDim_nod2d     
     row=nod3d_below_nod2d(1,m)
     coord_nod3d(3,row)=ssh(m) 
  end do
  
  ! update new vol and derivatives
  do el=1, myDim_elem3d                
     n_el=map_elem(el)
     if(n_el==0) cycle
     bafux_3d(:,el)=bafux_3d_new(:,n_el)
     bafuy_3d(:,el)=bafuy_3d_new(:,n_el)
     bafuz_3d(:,el)=bafuz_3d_new(:,n_el)
     voltetra(el)=voltetra_new(n_el)
     call local_element_def_3D(el, DET3D, derivative_locbafu_x_3D)
     do i=1,4
        bafux_3d_new(i,n_el) = derivative_locbafu_x_3D(1,i)
        bafuy_3d_new(i,n_el) = derivative_locbafu_x_3D(2,i)
        bafuz_3d_new(i,n_el) = derivative_locbafu_x_3D(3,i)
     enddo
     voltetra_new(n_el) = abs(DET3D)*Vol3D
     voltetra(el)=voltetra_new(n_el)
  enddo
end subroutine update_mesh
#endif
