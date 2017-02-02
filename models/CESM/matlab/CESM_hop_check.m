function CESM_hop_check(testname)
%% CESM_hop_check quantifies the difference between two runs of a model.
%  The complete filenames are contained within this script. YOU MUST EDIT THEM.
%

%% DART software - Copyright UCAR. This open source software is provided
% by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

switch ( testname )

   case 'rof_test' % RESULT: 

      file1 = '/glade/scratch/thoar/cesm_startup/run/cesm_startup.cam_0002.i.2004-01-04-00000.nc';
      file2 = '/glade/scratch/thoar/cesm_hybrid0/run/cesm_hybrid0.cam_0002.i.2004-01-04-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/cesm_startup/run/cesm_startup.pop_0002.r.2004-01-04-00000.nc';
      file2 = '/glade/scratch/thoar/cesm_hybrid0/run/cesm_hybrid0.pop_0002.r.2004-01-04-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/cesm_startup/run/cesm_startup.clm2_0002.r.2004-01-04-00000.nc';
      file2 = '/glade/scratch/thoar/cesm_hybrid0/run/cesm_hybrid0.clm2_0002.r.2004-01-04-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/cesm_startup/run/cesm_startup.cice_0002.r.2004-01-04-00000.nc';
      file2 = '/glade/scratch/thoar/cesm_hybrid0/run/cesm_hybrid0.cice_0002.r.2004-01-04-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/cesm_startup/run/cesm_startup.rtm_0002.r.2004-01-04-00000.nc';
      file2 = '/glade/scratch/thoar/cesm_hybrid0/run/cesm_hybrid0.rtm_0002.r.2004-01-04-00000.nc';
      Compare_netCDF_files(file1,file2)

   case 'clm_vanilla' % RESULT: mlaidiff differs

      file1 = '/glade/scratch/thoar/clm_hop_test/2day1hop/clm_hop_test.clm2_0001.r.2000-02-02-00000.nc';
      file2 = '/glade/scratch/thoar/clm_hop_test/2day2hop/clm_hop_test.clm2_0001.r.2000-02-02-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/clm_hop_test/2day1hop/clm_hop_test.clm2_0002.r.2000-02-02-00000.nc';
      file2 = '/glade/scratch/thoar/clm_hop_test/2day2hop/clm_hop_test.clm2_0002.r.2000-02-02-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/clm_hop_test/2day1hop/clm_hop_test.clm2_0003.r.2000-02-02-00000.nc';
      file2 = '/glade/scratch/thoar/clm_hop_test/2day2hop/clm_hop_test.clm2_0003.r.2000-02-02-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/clm_hop_test/2day1hop/clm_hop_test.clm2_0004.r.2000-02-02-00000.nc';
      file2 = '/glade/scratch/thoar/clm_hop_test/2day2hop/clm_hop_test.clm2_0004.r.2000-02-02-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/clm_hop_test/2day1hop/clm_hop_test.cpl.r.2000-02-02-00000.nc';
      file2 = '/glade/scratch/thoar/clm_hop_test/2day2hop/clm_hop_test.cpl.r.2000-02-02-00000.nc';
      Compare_netCDF_files(file1,file2)
      
   case 'clm_sourcemods' % RESULT: mlaidiff differs

      file1 = '/glade/scratch/thoar/clm_test/2day1hop_sourcemods/clm_test.clm2_0001.r.2000-02-02-00000.nc';
      file2 = '/glade/scratch/thoar/clm_test/2day2hop_sourcemods/clm_test.clm2_0001.r.2000-02-02-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/clm_test/2day1hop_sourcemods/clm_test.clm2_0002.r.2000-02-02-00000.nc';
      file2 = '/glade/scratch/thoar/clm_test/2day2hop_sourcemods/clm_test.clm2_0002.r.2000-02-02-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/clm_test/2day1hop_sourcemods/clm_test.clm2_0003.r.2000-02-02-00000.nc';
      file2 = '/glade/scratch/thoar/clm_test/2day2hop_sourcemods/clm_test.clm2_0003.r.2000-02-02-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/clm_test/2day1hop_sourcemods/clm_test.clm2_0004.r.2000-02-02-00000.nc';
      file2 = '/glade/scratch/thoar/clm_test/2day2hop_sourcemods/clm_test.clm2_0004.r.2000-02-02-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/clm_test/2day1hop_sourcemods/clm_test.cpl.r.2000-02-02-00000.nc';
      file2 = '/glade/scratch/thoar/clm_test/2day2hop_sourcemods/clm_test.cpl.r.2000-02-02-00000.nc';
      Compare_netCDF_files(file1,file2)
      
   case 'clm_overkill'  % RESULT: mlaidiff and pfts1d_ci differ

      file1 = '/glade/scratch/thoar/clm_hop_test/2day1hop/clm_hop_test.clm2_0001.r.2000-02-02-00000.nc';
      file2 = '/glade/scratch/thoar/clm_test/2day1hop_sourcemods/clm_test.clm2_0001.r.2000-02-02-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/clm_hop_test/2day1hop/clm_hop_test.clm2_0002.r.2000-02-02-00000.nc';
      file2 = '/glade/scratch/thoar/clm_test/2day1hop_sourcemods/clm_test.clm2_0002.r.2000-02-02-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/clm_hop_test/2day1hop/clm_hop_test.clm2_0003.r.2000-02-02-00000.nc';
      file2 = '/glade/scratch/thoar/clm_test/2day1hop_sourcemods/clm_test.clm2_0003.r.2000-02-02-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/clm_hop_test/2day1hop/clm_hop_test.clm2_0004.r.2000-02-02-00000.nc';
      file2 = '/glade/scratch/thoar/clm_test/2day1hop_sourcemods/clm_test.clm2_0004.r.2000-02-02-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/clm_hop_test/2day1hop/clm_hop_test.cpl.r.2000-02-02-00000.nc';
      file2 = '/glade/scratch/thoar/clm_test/2day1hop_sourcemods/clm_test.cpl.r.2000-02-02-00000.nc';
      Compare_netCDF_files(file1,file2)
      
   case 'clm_nullDART' % RESULT: mlaidiff differs

      file1 = '/glade/scratch/thoar/clm_test/2day2hop_sourcemods/clm_test.clm2_0001.r.2000-02-02-00000.nc';
      file2 = '/glade/scratch/thoar/clm_test/2day2hop_nullDART/clm_test.clm2_0001.r.2000-02-02-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/clm_test/2day2hop_sourcemods/clm_test.clm2_0002.r.2000-02-02-00000.nc';
      file2 = '/glade/scratch/thoar/clm_test/2day2hop_nullDART/clm_test.clm2_0002.r.2000-02-02-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/clm_test/2day2hop_sourcemods/clm_test.clm2_0003.r.2000-02-02-00000.nc';
      file2 = '/glade/scratch/thoar/clm_test/2day2hop_nullDART/clm_test.clm2_0003.r.2000-02-02-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/clm_test/2day2hop_sourcemods/clm_test.clm2_0004.r.2000-02-02-00000.nc';
      file2 = '/glade/scratch/thoar/clm_test/2day2hop_nullDART/clm_test.clm2_0004.r.2000-02-02-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/clm_test/2day2hop_sourcemods/clm_test.cpl.r.2000-02-02-00000.nc';
      file2 = '/glade/scratch/thoar/clm_test/2day2hop_nullDART/clm_test.cpl.r.2000-02-02-00000.nc';
      Compare_netCDF_files(file1,file2)
      
   case 'popnull'

      file1 = '/glade/scratch/thoar/pop_test/2day1hop_noFilter/pop_test.pop_0001.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/pop_test/2day2hop_noFilter/pop_test.pop_0001.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/pop_test/2day1hop_noFilter/pop_test.pop_0002.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/pop_test/2day2hop_noFilter/pop_test.pop_0002.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/pop_test/2day1hop_noFilter/pop_test.pop_0003.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/pop_test/2day2hop_noFilter/pop_test.pop_0003.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/pop_test/2day1hop_noFilter/pop_test.pop_0004.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/pop_test/2day2hop_noFilter/pop_test.pop_0004.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/pop_test/2day1hop_noFilter/pop_test.cice_0001.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/pop_test/2day2hop_noFilter/pop_test.cice_0001.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/pop_test/2day1hop_noFilter/pop_test.cice_0002.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/pop_test/2day2hop_noFilter/pop_test.cice_0002.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/pop_test/2day1hop_noFilter/pop_test.cice_0003.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/pop_test/2day2hop_noFilter/pop_test.cice_0003.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/pop_test/2day1hop_noFilter/pop_test.cice_0004.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/pop_test/2day2hop_noFilter/pop_test.cice_0004.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/pop_test/2day1hop_noFilter/pop_test.cpl.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/pop_test/2day2hop_noFilter/pop_test.cpl.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)

   case 'pop'

      file1 = '/glade/scratch/thoar/pop_test/2day1hop/pop_test.pop_0001.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/pop_test/2day2hop/pop_test.pop_0001.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/pop_test/2day1hop/pop_test.pop_0002.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/pop_test/2day2hop/pop_test.pop_0002.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/pop_test/2day1hop/pop_test.pop_0003.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/pop_test/2day2hop/pop_test.pop_0003.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/pop_test/2day1hop/pop_test.pop_0004.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/pop_test/2day2hop/pop_test.pop_0004.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/pop_test/2day1hop/pop_test.cice_0001.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/pop_test/2day2hop/pop_test.cice_0001.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/pop_test/2day1hop/pop_test.cice_0002.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/pop_test/2day2hop/pop_test.cice_0002.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/pop_test/2day1hop/pop_test.cice_0003.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/pop_test/2day2hop/pop_test.cice_0003.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/pop_test/2day1hop/pop_test.cice_0004.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/pop_test/2day2hop/pop_test.cice_0004.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/pop_test/2day1hop/pop_test.cpl.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/pop_test/2day2hop/pop_test.cpl.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)

   case 'cam'

      file1 = '/glade/scratch/thoar/cam_test2/2day1hop_initial/cam_test2.cam_0001.i.2008-11-06-00000.nc';
      file2 = '/glade/scratch/thoar/cam_test2/2day2hop_initial/cam_test2.cam_0001.i.2008-11-06-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/cam_test2/2day1hop_initial/cam_test2.cam_0002.i.2008-11-06-00000.nc';
      file2 = '/glade/scratch/thoar/cam_test2/2day2hop_initial/cam_test2.cam_0002.i.2008-11-06-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/cam_test2/2day1hop_initial/cam_test2.clm2_0001.r.2008-11-06-00000.nc';
      file2 = '/glade/scratch/thoar/cam_test2/2day2hop_initial/cam_test2.clm2_0001.r.2008-11-06-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/cam_test2/2day1hop_initial/cam_test2.clm2_0002.r.2008-11-06-00000.nc';
      file2 = '/glade/scratch/thoar/cam_test2/2day2hop_initial/cam_test2.clm2_0002.r.2008-11-06-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/cam_test2/2day1hop_initial/cam_test2.cice_0001.r.2008-11-06-00000.nc';
      file2 = '/glade/scratch/thoar/cam_test2/2day2hop_initial/cam_test2.cice_0001.r.2008-11-06-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/cam_test2/2day1hop_initial/cam_test2.cice_0002.r.2008-11-06-00000.nc';
      file2 = '/glade/scratch/thoar/cam_test2/2day2hop_initial/cam_test2.cice_0002.r.2008-11-06-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/cam_test2/2day1hop_initial/cam_test2.cpl.r.2008-11-06-00000.nc';
      file2 = '/glade/scratch/thoar/cam_test2/2day2hop_initial/cam_test2.cpl.r.2008-11-06-00000.nc';
      Compare_netCDF_files(file1,file2)

   case 'Bstartdate'

      file1 = '/glade/scratch/thoar/cesm_test/old_startdate/cesm_test.cam_0001.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/cesm_test/new_startdate/cesm_test.cam_0001.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/cesm_test/old_startdate/cesm_test.cam_0002.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/cesm_test/new_startdate/cesm_test.cam_0002.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/cesm_test/old_startdate/cesm_test.clm2_0001.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/cesm_test/new_startdate/cesm_test.clm2_0001.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/cesm_test/old_startdate/cesm_test.clm2_0002.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/cesm_test/new_startdate/cesm_test.clm2_0002.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/cesm_test/old_startdate/cesm_test.pop_0001.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/cesm_test/new_startdate/cesm_test.pop_0001.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/cesm_test/old_startdate/cesm_test.pop_0002.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/cesm_test/new_startdate/cesm_test.pop_0002.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/cesm_test/old_startdate/cesm_test.rtm_0001.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/cesm_test/new_startdate/cesm_test.rtm_0001.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/cesm_test/old_startdate/cesm_test.rtm_0002.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/cesm_test/new_startdate/cesm_test.rtm_0002.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/cesm_test/old_startdate/cesm_test.cice_0001.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/cesm_test/new_startdate/cesm_test.cice_0001.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/cesm_test/old_startdate/cesm_test.cice_0002.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/cesm_test/new_startdate/cesm_test.cice_0002.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      file1 = '/glade/scratch/thoar/cesm_test/old_startdate/cesm_test.cpl.r.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/cesm_test/new_startdate/cesm_test.cpl.r.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)

   case 'Bhoptest'

      disp('Testing CAM initial files.')
      file1 = '/glade/scratch/thoar/cesm_test/3day1hop/cesm_test.cam_0001.i.2004-01-07-00000.nc';
      file2 = '/glade/scratch/thoar/cesm_test/3day3hop/cesm_test.cam_0001.i.2004-01-07-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      disp('Testing CAM restart files.')
      file1 = '/glade/scratch/thoar/cesm_test/3day1hop/cesm_test.cam_0001.r.2004-01-07-00000.nc';
      file2 = '/glade/scratch/thoar/cesm_test/3day3hop/cesm_test.cam_0001.r.2004-01-07-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      disp('Testing CLM restart files on day 0 ... identical?')
      file1 = '/glade/scratch/thoar/cesm_test/3day1hop/cesm_test.clm2_0001.r.2004-01-04-00000.nc';
      file2 = '/glade/scratch/thoar/cesm_test/3day3hop/cesm_test.clm2_0001.r.2004-01-04-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause
      disp('Testing CLM restart files on day +3 ... changed by CAM state?')
      file1 = '/glade/scratch/thoar/cesm_test/3day1hop/cesm_test.clm2_0001.r.2004-01-07-00000.nc';
      file2 = '/glade/scratch/thoar/cesm_test/3day3hop/cesm_test.clm2_0001.r.2004-01-07-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      disp('Testing CLM history files on day +1 ... identical?')
      file1 = '/glade/scratch/thoar/cesm_test/3day1hop/cesm_test.clm2_0001.h0.2004-01-05-00000.nc';
      file2 = '/glade/scratch/thoar/cesm_test/3day3hop/cesm_test.clm2_0001.h0.2004-01-05-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause
      disp('Testing CLM history files on day +2 ... changed by CAM state?')
      file1 = '/glade/scratch/thoar/cesm_test/3day1hop/cesm_test.clm2_0001.h0.2004-01-06-00000.nc';
      file2 = '/glade/scratch/thoar/cesm_test/3day3hop/cesm_test.clm2_0001.h0.2004-01-06-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      disp('Testing POP restart files.')
      file1 = '/glade/scratch/thoar/cesm_test/3day1hop/cesm_test.pop_0001.r.2004-01-07-00000.nc';
      file2 = '/glade/scratch/thoar/cesm_test/3day3hop/cesm_test.pop_0001.r.2004-01-07-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      disp('Testing RTM restart files.')
      file1 = '/glade/scratch/thoar/cesm_test/3day1hop/cesm_test.rtm_0001.r.2004-01-07-00000.nc';
      file2 = '/glade/scratch/thoar/cesm_test/3day3hop/cesm_test.rtm_0001.r.2004-01-07-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      disp('Testing CICE restart files.')
      file1 = '/glade/scratch/thoar/cesm_test/3day1hop/cesm_test.cice_0001.r.2004-01-07-00000.nc';
      file2 = '/glade/scratch/thoar/cesm_test/3day3hop/cesm_test.cice_0001.r.2004-01-07-00000.nc';
      Compare_netCDF_files(file1,file2)
      disp('Pausing, hit any key to continue ...')
      pause

      disp('Testing CPL restart file.')
      file1 = '/glade/scratch/thoar/cesm_test/3day1hop/cesm_test.cpl.r.2004-01-07-00000.nc';
      file2 = '/glade/scratch/thoar/cesm_test/3day3hop/cesm_test.cpl.r.2004-01-07-00000.nc';
      Compare_netCDF_files(file1,file2)

   case 'B3day_h'

      Compare_netCDF_files('/glade/scratch/thoar/cesm_test/cesm_test.pop_0001.hv.diff.nc')
      Compare_netCDF_files('/glade/scratch/thoar/cesm_test/cesm_test.pop_0001.h.once.diff.nc')
      Compare_netCDF_files('/glade/scratch/thoar/cesm_test/cam_0001.h0.2004-01-04-00000.diff.nc')
      Compare_netCDF_files('/glade/scratch/thoar/cesm_test/cam_0001.h0.2004-01-06-00000.diff.nc')
      Compare_netCDF_files('/glade/scratch/thoar/cesm_test/cam_0001.h0.2004-01-07-00000.diff.nc')
      Compare_netCDF_files('/glade/scratch/thoar/cesm_test/cam_0001.h0.2004-01-05-00000.diff.nc')
      Compare_netCDF_files('/glade/scratch/thoar/cesm_test/clm2_0001.h1.2004-01-04-00000.diff.nc')
      Compare_netCDF_files('/glade/scratch/thoar/cesm_test/clm2_0001.h0.2004-01-05-00000.diff.nc')
      Compare_netCDF_files('/glade/scratch/thoar/cesm_test/clm2_0001.h1.2004-01-05-00000.diff.nc')
      Compare_netCDF_files('/glade/scratch/thoar/cesm_test/clm2_0001.h1.2004-01-06-00000.diff.nc')
      Compare_netCDF_files('/glade/scratch/thoar/cesm_test/clm2_0001.h1.2004-01-07-00000.diff.nc')

   otherwise

      error('no known case of name %s',testname)

end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$
