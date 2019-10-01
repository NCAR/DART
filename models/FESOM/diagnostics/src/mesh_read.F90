
! FESOM 2 (Finite-volumE Sea ice-Ocean Model)
! multi-resolution ocean general circulation model
! FESOM/fesom2 is licensed under the GNU General Public License v2.0
! Copyright (C) 2018  FESOM team
!
! This source file was taken from  the FESOM V1.4 modules


subroutine read_elem
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!! READ elem3d.out !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use o_param
  use o_elements
  use o_mesh
  use o_read
  use g_parfe

  mype=0
  fileID=mype+10
  open(fileID, file=(trim(MeshPath)//'elem3d.out'))
    read(fileID,*) myDim_elem3D   
    allocate(elem3D_nodes(4, myDim_elem3D))
      read(fileID,*) elem3d_nodes
  close(fileID)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!! READ elem2d.out !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  fileID=mype+10
  open(fileID, file=(trim(MeshPath)//'elem2d.out'))
    read(fileID,*) myDim_elem2D   
    allocate(elem2D_nodes(3, myDim_elem2D))
      read(fileID,*) elem2d_nodes
  close(fileID)
  print*, 'Element tables are read: '

end subroutine read_elem
subroutine read_node
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!! READ nod3d.out !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use o_param
  use o_elements
  use o_mesh
  use o_read
  use g_parfe

  mype=0
  fileID=mype+10
  open(fileID, file=(trim(MeshPath)//'nod3d.out'))
    read(fileID,*) myDim_nod3D
    allocate(coord_nod3D(3,myDim_nod3D))
    allocate(index_nod3D(myDim_nod3D))
    allocate(myList_nod3D(myDim_nod3D))
      do n=1,myDim_nod3D
        read(fileID,*) m, x, y, z, ind
        coord_nod3D(1,n)=x
        coord_nod3D(2,n)=y
        coord_nod3D(3,n)=z
        index_nod3D(n)=ind
        myList_nod3D(n)=m
     end do
  close(fileID)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!! READ nod2d.out !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  fileID=mype+10
  open(fileID, file=(trim(MeshPath)//'nod2d.out'))
    read(fileID,*) myDim_nod2D
    allocate(coord_nod2D(2,myDim_nod2D))
    allocate(index_nod2D(myDim_nod2D))
    allocate(myList_nod2D(myDim_nod3D))
      do n=1,myDim_nod2D
        read(fileID,*) m, x, y, ind
        coord_nod2D(1,n)=x
        coord_nod2D(2,n)=y
        index_nod2D(n)=ind
        myList_nod2D(n)=m
      end do
  close(fileID)
  print*, 'Node tables are read: '

end subroutine read_node

subroutine read_aux3

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!! READ aux3d.out !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use o_param
  use o_elements
  use o_mesh
  use o_read
  use g_parfe

  mype=0
  fileID=mype+10
  open(fileID, file=trim(meshpath)//'aux3d.out')
  read(fileID, *) max_num_layers

!!!!===========================================================!!!!
!!!!=========nod3D_below_nod2D=================================!!!!
!!!!===========================================================!!!!

  allocate(nod3D_below_nod2D(max_num_layers,myDim_nod2D))       
  do n=1, myDim_nod2D
     read(fileID, *) vert_nodes(1:max_num_layers)
        nod3D_below_nod2D(:,n)=vert_nodes(1:max_num_layers)
  end do

  do n=1, myDim_nod2D
     do m=1, max_num_layers
        n1=nod3D_below_nod2D(m,n)
        if(n1>0)  then
           nod3D_below_nod2D(m,n)=n1
        end if
     end do
  end do

!!!!===========================================================!!!!
!!!!=========nod2D_corresp_to_nod3D============================!!!!
!!!!===========================================================!!!!

  allocate(nod2D_corresp_to_nod3D(myDim_nod3D)) 

  do n=1, myDim_nod3D
     read(fileID, *) m
        nod2D_corresp_to_nod3D(n)=m
  end do

  do n=1,myDim_nod3D
     m=nod2D_corresp_to_nod3D(n)
     nod2D_corresp_to_nod3D(n)=m
  end do

!!!!===========================================================!!!!
!!!!========elem2D_corresp_to_elem3D===========================!!!!
!!!!===========================================================!!!!

  allocate(elem2D_corresp_to_elem3D(myDim_elem3D)) 
  do n=1, myDim_elem3D
     read(fileID,*) m
        elem2D_corresp_to_elem3D(n)=m
  end do

  close(fileID)
  print*, 'Aux3D tables are read: '

end subroutine read_aux3

subroutine read_depth
  use o_param
  use o_elements
  use o_mesh
  use o_read
  use g_parfe

  real(kind=8) m
  fileID=mype+10
  allocate(layerdepth(max_num_layers)) 
  open(fileID, file=(trim(MeshPath)//'m3d.ini'))
  do n=1,max_num_layers-1
     read(fileID,*) m
       layerdepth(n)=m
  end do 
  close(fileID)
end subroutine read_depth
