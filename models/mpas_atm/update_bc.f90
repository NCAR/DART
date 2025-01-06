! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program update_bc

!----------------------------------------------------------------------
! purpose: interface between DART and the model model
!
! method: Read DART analysis vector in netcdf and replace the corresponding
!         field in the mpas file to advance model after running this program.
!         Updated to process all ensemble members.
!
!         The update_bc_nml namelist defines the input and output file
!         name lists for all ensemble members.
!         The input list should be matched with output_state_file_list in &filter_nml.
!         
! author: Soyoung Ha 23 Aug 16
!         Updated in  4 May 2017 for the Manhatten release
!         Updated in 28 Jul 2020 for checking dimension sizes in analysis and lbc files.
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, open_file, close_file, &
                             get_next_filename, E_ERR, E_MSG, error_handler
use time_manager_mod, only : time_type, print_time, print_date, operator(-), &
                             get_time, get_date, operator(/=)
use direct_netcdf_mod,only : read_variables  !HK read_variables?!
use        model_mod, only : static_init_model, &
                             get_model_size,    &
                             get_analysis_time, &
                             anl_domid, lbc_domid, &
                             uv_increments_cell_to_edges, &
                             uv_field_cell_to_edges, &
                             get_analysis_weight, &
                             on_boundary_cell, &
                             on_boundary_edge, &
                             cell_next_to_boundary_edge
                             

use state_structure_mod, only : get_num_variables, get_domain_size, &
                                get_variable_name


use netcdf_utilities_mod, only : nc_open_file_readonly, &
                                 nc_open_file_readwrite, &
                                 nc_get_dimension_size,   &
                                 nc_close_file, &
                                 NF90_MAX_NAME

implicit none

character(len=*), parameter :: source   = 'models/mpas_atm/update_bc.f90'

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------
character(len=256)  :: update_analysis_file_list = 'filter_in.txt'
character(len=256)  :: update_boundary_file_list = 'boundary_inout.txt'
logical             :: lbc_update_from_reconstructed_winds = .true.
logical             :: lbc_update_winds_from_increments    = .true.


namelist /update_bc_nml/ update_analysis_file_list, update_boundary_file_list, &
                         lbc_update_from_reconstructed_winds, lbc_update_winds_from_increments

!----------------------------------------------------------------------
character (len=256)   :: next_infile, next_outfile
character (len=256)   :: bdy_template_filename
character (len=512)   :: string1, string2
integer               :: iunit, io
integer               :: ncAnlID, ncBdyID
integer               :: filenum
type(time_type)       :: model_time
type(time_type)       :: state_time
integer :: nCellsA      = -1  ! Total number of cells  in ncAnlID
integer :: nCellsB      = -1  ! Total number of cells  in ncBdyID
integer :: nVertLevelsA = -1  ! Total number of levels in ncAnlID
integer :: nVertLevelsB = -1  ! Total number of levels in ncBdyID

integer :: nCells        = -1  ! Total number of cells making up the grid
integer :: nVertices     = -1  ! Unique points in grid that are corners of cells
integer :: nEdges        = -1  ! Straight lines between vertices making up cells
integer :: nVertLevels   = -1  ! Vertical levels; count of vert cell centers
integer :: vertexDegree  = -1  ! Max number of cells/edges that touch any vertex
integer :: nSoilLevels   = -1  ! Number of soil layers

integer :: cellid, edgeid, varid

real(r8), allocatable, dimension(:,:) :: lbc_u, lbc_ucell, lbc_vcell, inc_lbc_ucell, inc_lbc_vcell
real(r8), allocatable, dimension(:,:) :: old_lbc_ucell, old_lbc_vcell, delta_u
real(r8), allocatable, dimension(:,:) :: a_var_data, b_var_data, var_data
real(r8) :: weight

character(len=NF90_MAX_NAME) :: bvarname, avarname

!----------------------------------------------------------------------

call initialize_utilities(progname=source)

! Read the namelist to get the input filename. 
call find_namelist_in_file("input.nml", "update_bc_nml", iunit)
read(iunit, nml = update_bc_nml, iostat = io)
call check_namelist_read(iunit, io, "update_bc_nml")

call static_init_model()
call get_grid_dims(nCells, nVertices, nEdges, nVertLevels, vertexDegree, nSoilLevels)
!----------------------------------------------------------------------
! Reads lists of input mpas (prior) and filter (analysis) files 
!HK @todo loop around files, why not run this code in parallel?
!----------------------------------------------------------------------
filenum = 1
fileloop: do        ! until out of files  

   ! get a file name from the list, one at a time.
   next_infile  = get_next_filename(update_analysis_file_list, filenum)
   next_outfile = get_next_filename(update_boundary_file_list, filenum)
   if (next_infile == '' .or. next_outfile == '') exit fileloop

   !----------------------------------------------------------------------
   ! Reads input lbc (prior) and filter (analysis) files 
   !----------------------------------------------------------------------
   ncAnlID = nc_open_file_readonly(next_infile,  'update_bc - open readonly')   ! analysis file from DART (ouput from filter for anl_dom)
   ncBdyID = nc_open_file_readwrite(next_outfile, 'update_bc - open readwrite') ! prior boundary, original mpas file

   !----------------------------------------------------------------------
   ! Read the model time
   !----------------------------------------------------------------------
   state_time = get_analysis_time(ncAnlID, next_infile)
   model_time = get_analysis_time(ncBdyID, next_outfile)
   call print_time(state_time,'DART current time')
   call print_time(model_time,'mpas current time')

   if ( model_time /= state_time ) then
      call print_time(state_time,'DART current time',logfileunit)
      call print_time(model_time,'mpas current time',logfileunit)
      write(string1,*) trim(next_infile),' current time must equal model time'
      call error_handler(E_ERR,source,string1,source)
   endif

   !----------------------------------------------------------------------
   ! Check dimension size in both files
   !----------------------------------------------------------------------
   nCellsA      = nc_get_dimension_size(ncAnlID, 'nCells',      source)
   nCellsB      = nc_get_dimension_size(ncBdyID, 'nCells',      source)
   nVertLevelsA = nc_get_dimension_size(ncAnlID, 'nVertLevels', source)
   nVertLevelsB = nc_get_dimension_size(ncBdyID, 'nVertLevels', source)
  
   if((nCellsA /= nCellsB) .or. (nVertLevelsA /= nVertLevelsB)) then  
      ! HK @todo also check against static_init_model values
      write(string1,*) 'Domain size mismatches'
      call error_handler(E_ERR,'update_bc',string1,source)
   endif

   if (lbc_update_from_reconstructed_winds) then ! save a copy of the reconstruced cell winds 
      
      if (.not. lbc_file_has_reconstructed_winds) then
         write(string1, *) 'Cannot update edge winds from increments because the boundary file does not contain the reconstructed winds (lbc_ur, lbc_vr)'
         write(string2, *) 'lbc_update_winds_from_increments should be .false.'
         call error_handler(E_MSG,'statevector_to_boundary_file',string1,&
                            source, text2=string2)
      endif

      allocate(old_lbc_ucell(nVertLevels, nCells))
      allocate(old_lbc_vcell(nVertLevels, nCells))
      allocate(      delta_u(nVertLevels, nEdges))

      call nc_get_variable(ncBdyID, 'uReconstructMeridional', old_lbc_ucell)
      call nc_get_variable(ncBdyID, 'uReconstructZonal', old_lbc_vcell)

   endif

 
   ! Update variables except 'u' from the analysis
   allocate(a_var_data(nVertLevels, nCells))
   allocate(b_var_data(nVertLevels, nCells))
   allocate(var_data(nVertLevels, nCells))

   VARLOOP: do varid = 1, get_num_variables(lbc_domid)
      bvarname = get_variable_name(lbc_domid, varid)
      avarname = trim(bvarname(5:)) !corresponding field in analysis domain

      ! skip edge normal winds
      if (bvarname == 'u') cycle VARLOOP

      ! reconstructed cell-center winds have different names in the lbc file.
      if (bvarname == 'lbc_ur') avarname = 'uReconstructZonal'
      if (bvarname == 'lbc_vr') avarname = 'uReconstructMeridional'    

      if (bvarname(1:4) /= 'lbc_') then
         write(string1, *) 'skipping update of boundary variable ', trim(bvarname)  !HK @todo why do you need to tell people this?
         write(string2, *) 'because the name does not start with "lbc"'
         call error_handler(E_MSG,'statevector_to_boundary_file',string1,&
                              source, text2=string2)
         cycle VARLOOP
      endif

      call nc_get_variable(ncAnlID, avarname, a_var_data)
      call nc_get_variable(ncBdyID, bvarname, b_var_data)
      
      ! for each cell in the grid, find the analysis in the
      ! boundary region and blend them with prior lbc values.
      CELLS: do cellid = 1, nCells

         if (.not. on_boundary_cell(cellid)) cycle CELLS

         weight = get_analysis_weight(cellid) ! 1.0 is interior, 0.0 is exterior boundary
         var_data(:, cellid) = (1.0_r8 - weight) * b_var_data(:, cellid) + weight * a_var_data(:, cellid)
 
      enddo CELLS
 
      call nc_put_variable(ncBdyID, bvarname, var_data)
  
   enddo VARLOOP

   deallocate(a_var_data, b_var_data, var_data)

   ! deal with u edge winds
   if (.not. lbc_update_from_reconstructed_winds) then  ! update u edge winds directly

      allocate(a_var_data(nVertLevels, nEdges))
      allocate(b_var_data(nVertLevels, nEdges))
      allocate(var_data(nVertLevels, nEdges))

      call nc_get_variable(ncAnlID, 'u', a_var_data)      ! analysis edge winds
      call nc_get_variable(ncBdyID, 'lbc_u', b_var_data)  ! prior edge winds in the lbc file
       
     ! for each edge in the grid, find the ones which are in the
     ! boundary region and blend their values.
      EDGES: do edgeid = 1, nEdges
  
        if (.not. on_boundary_edge(edgeid)) cycle EDGES
  
        weight = get_analysis_weight(edgeid, .false.) ! 1.0 is interior, 0.0 is exterior boundary
        var_data(:, edgeid) = (1.0_r8 - weight) * b_var_data(:,edgeid) + weight * a_var_data(:,edgeid)
  
      enddo EDGES
  
      call nc_put_variable(ncBdyID, 'lbc_u', var_data)

   else ! update u edge winds from reconstructed winds

      allocate(        lbc_u(nVertLevels, nEdges))
      allocate(    lbc_ucell(nVertLevels, nCells))
      allocate(    lbc_vcell(nVertLevels, nCells))
  
      call nc_get_variable(ncBdyID, 'lbc_ur', lbc_vcell)  ! already blended in VARLOOP
      call nc_get_variable(ncBdyID, 'lbc_vr', lbc_ucell)  ! already blended in VALOOP
  
      if (lbc_update_winds_from_increments) then
  
         call nc_get_variable(ncBdyID, 'lbc_u', lbc_u)  ! not blended
         allocate(inc_lbc_ucell(nVertLevels, nCells))
         allocate(inc_lbc_vcell(nVertLevels, nCells))
         allocate(delta_u(nVertLevels, nEdges))

         inc_lbc_ucell = lbc_ucell - old_lbc_ucell 
         inc_lbc_vcell = lbc_vcell - old_lbc_vcell 
   
         delta_u =  uv_increments_cell_to_edges(inc_lbc_ucell, inc_lbc_vcell)
  
         IEDGE: do edgeid = 1, nEdges

            ! Soyoung: Add the blended u increments back to lbc_u in the boundary zone.
            !          We should not change the analysis u in the interior domain, but
            !          We should also check bdyMaskCell for the two adjacent cells as
            !          bdyMaskEdge is assigned with the lower mask value between the two
            !          cells.
            ! Ex) An edge between cell1 (w/ bdyMaskCell = 0) and cell2 (w/ bdyMaskCell = 1)
            ! has bdyMaskEdge = 0. In this case, even if bdyMaskEdge of the edge is zero,
            ! cell2 has been updated in the CELLS loop above, so the edge has to be updated.
   
            if (.not. on_boundary_edge(edgeid) .and. .not. cell_next_to_boundary_edge(edgid)) cycle IEDGE
  
            lbc_u(:,edgeid) = lbc_u(:,edgeid) + delta_u(:,edgeid)
  
         enddo IEDGE
  
         call nc_put_variable(ncBdyID, 'lbc_u', lbc_u)
         deallocate(old_lbc_ucell, old_lbc_vcell, inc_lbc_ucell,inc_lbc_vcell, delta_u)    

      else ! just replace, no increments

        call uv_field_cell_to_edges(lbc_ucell, lbc_vcell, lbc_u)
        call nc_put_variable(ncBdyID, 'lbc_u', var_data) 

      endif

   endif    
    
  call print_date( model_time,'update_bc:model date')
  call print_time( model_time,'update_bc:model time')
  call print_date( model_time,'update_bc:model date',logfileunit)
  call print_time( model_time,'update_bc:model time',logfileunit)

  call nc_close_file(ncAnlID, source)
  call nc_close_file(ncBdyID, source)

  filenum = filenum + 1

end do fileloop

call finalize_utilities()

end program update_bc
