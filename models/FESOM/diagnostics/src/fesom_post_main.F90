
!> FESOM diagnostic post-processor

!> Program to extract the data for the different graphics from either FESOM
!> or the DART diagnostic output of an experiment with FESOM.
!> fesom_post_main is intended to be called by the 'master' scripts
!> in diagnostics/scripts which modify the namelists appropriately.
!>
!> fesom_post_main reads the namelist.config:&postproc namelist to select 
!> which data to produce. The following table relates the value of
!> postproc:tool to the data being produced.
!>
!>  tool =  1   basin_mean_evolution
!>  tool =  2   read_thalweg_from_nc
!>  tool =  3   read_section_from_netcdf
!>  tool =  4   marmara_mean_evolution
!>  tool =  5   calc_section_monthly_mean
!>  tool =  6   calc_thalweg_monthly_mean
!>  tool =  7   read_ensemble_from_netcdf
!>  tool =  8   synthetic_ferrybox_from_nr
!>  tool =  9   read_section_from_NR_diff
!>  tool = 10   read_ctd_data
!>  tool = 11   profile_from_netcdf
!>  tool = 12   velocity_at_the_exit
!>  tool = 13   read_section_from_inc
!>  tool = 14   dardanelles_for_MFS
!>  tool = 15   total_kinetic_energy
!>  tool = 16   surface_kinetic_energy
!>  tool = 17   calc_section_annual_mean
!>  tool = 18   calc_thalweg_annual_mean
!>  tool = 19   compute_vorticity
!>  tool = 20   compute_wind_stress_curl
!>  tool = 21   compute_net_flux
!>  tool = 22   compute_surface_buoyancy
!>  tool = 23   compute_forcing_monthly_timeseries
!>  tool = 24   compute_wind_work
!>  tool = 25   read_ship_track
!>  tool = 26   bosphorus_for_blk_mfs
!>  tool = 27   compute_volume_transport

program fesom_post_main

  use fesom_forcing_mod,         only : compute_forcing_monthly_timeseries, &
                                        compute_wind_stress_curl, &
                                        compute_surface_buoyancy, &
                                        compute_wind_work, &
                                        forcing_array_setup

  use fesom_ocean_mod,           only : find_surface_area, &
                                        basin_mean_evolution, &
                                        read_thalweg_from_nc, &
                                        read_section_from_netcdf, &
                                        marmara_mean_evolution, &
                                        calc_section_monthly_mean, &
                                        calc_thalweg_monthly_mean, &
                                        velocity_at_the_exit, &
                                        dardanelles_for_MFS, &
                                        total_kinetic_energy, &
                                        surface_kinetic_energy, &
                                        calc_section_annual_mean, &
                                        calc_thalweg_annual_mean, &
                                        compute_vorticity, &
                                        compute_net_flux, &
                                        bosphorus_for_blk_mfs, &
                                        compute_volume_transport

  use fesom_dart_mod,            only : read_ensemble_from_netcdf, &
                                        read_section_from_NR_diff, &
                                        read_section_from_inc

  use fesom_observation_mod,     only : synthetic_ferrybox_from_nr, &
                                        read_ctd_data, &
                                        profile_from_netcdf, &
                                        read_ship_track

  use g_config,                  only : tool

  real              :: t0, t1, t2, t3, t4, t5, t6, t7, t8
  call read_namelist
  call cpu_time(t0); print*, 'time elapsed:',t0
  call read_elem
  call cpu_time(t1); print*, 'time elapsed:',t1
  call read_node
  call cpu_time(t2); print*, 'time elapsed:',t2
  call read_aux3
  call cpu_time(t3); print*, 'time elapsed:',t3
  call read_depth
  call cpu_time(t4); print*, 'time elapsed:',t4
  call ocean_array_setup
  call cpu_time(t5); print*, 'time elapsed:',t5
  call ocean_mesh_setup
  call cpu_time(t6); print*, 'time elapsed:',t6
  call find_cluster_area
  call cpu_time(t7); print*, 'time elapsed:',t7
  call find_surface_area
  if ( tool.eq.1 ) then
     call basin_mean_evolution
  else if ( tool.eq.2 ) then
     call read_thalweg_from_nc
  else if ( tool.eq.3 ) then
     call read_section_from_netcdf
  else if ( tool.eq.4 ) then
     call marmara_mean_evolution
  else if ( tool.eq.5 ) then
     call calc_section_monthly_mean
  else if ( tool.eq.6 ) then
     call calc_thalweg_monthly_mean
  else if ( tool.eq.7 ) then
     call read_ensemble_from_netcdf
  else if ( tool.eq.8 ) then
     call synthetic_ferrybox_from_nr
  else if ( tool.eq.9 ) then
     call read_section_from_NR_diff
  else if ( tool.eq.10 ) then
     call read_ctd_data
  else if ( tool.eq.11 ) then
     call profile_from_netcdf
  else if ( tool.eq.12 ) then
     call velocity_at_the_exit
  else if ( tool.eq.13) then
     call read_section_from_inc
  else if ( tool.eq.14 ) then
     call dardanelles_for_MFS
  else if ( tool.eq.15 ) then
     call total_kinetic_energy
  else if ( tool.eq.16 ) then
     call surface_kinetic_energy
  else if ( tool.eq.17 ) then
     call calc_section_annual_mean
  else if ( tool.eq.18 ) then
     call calc_thalweg_annual_mean
  else if ( tool.eq.19 ) then
     call compute_vorticity
  else if ( tool.eq.20 ) then
     call forcing_array_setup
     call compute_wind_stress_curl
  else if ( tool.eq.21 ) then
     call compute_net_flux
  else if ( tool.eq.22 ) then
     call forcing_array_setup
     call compute_surface_buoyancy
  else if ( tool.eq.23 ) then
     call forcing_array_setup
     call compute_forcing_monthly_timeseries
  else if ( tool.eq.24 ) then
     call forcing_array_setup
     call compute_wind_work
  else if ( tool.eq.25 ) then
     call read_ship_track
  else if ( tool.eq.26 ) then
     call bosphorus_for_blk_mfs
  else if ( tool.eq.27 ) then
     call compute_volume_transport
  end if
  call cpu_time(t8)
end program fesom_post_main
